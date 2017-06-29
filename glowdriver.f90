program glowdriver

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Version 0.981, 6/2017

! Stan Solomon and Ben Foster, 1/2015
! Stan Solomon, 12/2015, 1/2016
! Stan Solomon, 3/2016, MPI parallel version

! Main multi-processor driver for the GLOW model.
! Uses TIE-GCM history files or MSIS/IRI for input.
! Requires MPI and netCDF libraries

! For definitions of use-associated variables, see subroutine GLOW and module CGLOW.
! For definitions of TGCM input variables see module READTGCM
! For definitions of output arrays see module OUTPUT

! Other definitions:
! f107p   Solar 10.7 cm flux for previous day
! ap      Ap index of geomagnetic activity
! z       altitude array, km

! Array dimensions:
! jmax    number of altitude levels
! nbins   number of energetic electron energy bins
! lmax    number of wavelength intervals for solar flux
! nmaj    number of major species
! nst     number of states produced by photoionization/dissociation
! nei     number of states produced by electron impact
! nex     number of ionized/excited species
! nw      number of airglow emission wavelengths
! nc      number of component production terms for each emission

  use mpi

  use cglow,only: cglow_init      ! subroutine to allocate use-associated variables
  use cglow,only: jmax,nbins,lmax,nmaj,nei,nex,nw,nc,nst
  use cglow,only: idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec
  use cglow,only: iscale,jlocal,kchem,xuvfac
  use cglow,only: sza,dip,efrac,ierr
  use cglow,only: zz,zo,zn2,zo2,zns,znd,zno,ztn,ze,zti,zte
  use cglow,only: ener,del,phitop,wave1,wave2,sflux,pespec,sespec,uflx,dflx,sion
  use cglow,only: photoi,photod,phono,aglw,ecalc,zxden,zeta,zceta,zlbh
  use cglow,only: data_dir

  use readtgcm,only: read_tgcm           ! subroutine to read tgcm history file
  use readtgcm,only: read_tgcm_coords    ! subroutine to read tgcm history coordinates
  use readtgcm,only: find_mtimes         ! subroutine to find model times on history file
  use readtgcm,only: nlon_tgcm=>nlon, glon_tgcm=>glon
  use readtgcm,only: nlat_tgcm=>nlat, glat_tgcm=>glat
  use readtgcm,only: nlev_tgcm=>nlev
  use readtgcm,only: iyear_tgcm=>iyear,iday_tgcm=>iday,ut_tgcm=>ut
  use readtgcm,only: f107_tgcm=>f107d,f107a_tgcm=>f107a,hpower
  use readtgcm,only: alfacusp=>alfac,ecusp=>ec,alfadrizzle=>alfad,edrizzle=>ed
  use readtgcm,only: eflux,nflux,alfa,drizzle,cusp
  use readtgcm,only: zg,tn,un,vn,o2,o1,n2,n4s,n2d,no,ti,te,ne
  use readtgcm,only: mxtimes

  use output,only: output_init    ! subroutine to allocate output arrays
  use output,only: create_ncfile  ! subroutine to create netcdf output file
  use output,only: write_ncfile   ! subroutine to write data to netcdf output file
  use output,only: nlon,nlat,nlev ! grid dimensions (set here in glow_drv)
  use output,only: lon,lat,lev    ! grid coordinates (set here in glow_drv)
  use output,only: writelbh,writered        ! switches
  use output,only: zzz,ao,ao2,an2,ano,an4s,an2d,atn,ati,ate,ane,aun,avn, &
                   nec,ped,hall,xden,eta,lbh,redline ! output arrays

  implicit none

  character(len=1024) :: &
    tgcm_ncfile,         &    ! path to tgcm history file (tiegcm or timegcm)
    glow_ncfile,         &    ! name for of netCDF glow output files
    glow_ncfileit,       &    ! path to individual netCDF glow output file (integer appended)
    iri90_dir                 ! directory containing iri data files
  character(len=7) :: ifile

  real,allocatable :: z(:)            ! glow height coordinate in km (jmax)
  real,allocatable :: zun(:),zvn(:)   ! winds on glow grid
  real,allocatable :: pedcond(:), hallcond(:)  ! Pederson and Hall conductivities in S/m (mho)
  real,allocatable :: outf(:,:)       ! iri output (11,jmax)
  real,allocatable :: recbuf(:,:)     ! receive buffer for MPI gather
  real :: utstart,utstep,utstop
  real :: rz12,stl,fmono,emono,kp
  real :: d(8), t(2), sw(25), oarr(30)
  integer nproc,itask,mpierr,lat0,lat1,size2d,size3d,sizeb
  integer :: l,j,jj,ijf,jmag,iday,mmdd,i,ii,n,k,ix,itail,m
  integer :: nlat_msis, nlon_msis  ! number of lats and lons in grid for MSIS/IRI runs
  integer :: start_mtime(3)        ! model start time (day,hour,minute)
  integer :: stop_mtime(3)         ! model stop time (day,hour,minute)
  integer :: indate,itime,ntimes,mtimes(3,mxtimes),itimes(mxtimes)
  logical :: tgcm, first, jf(12)
  data sw/25*1./, first/.true./

  namelist /glow_input/ &
    indate,utstart,utstep,utstop,nlat_msis,nlon_msis,f107a,f107,f107p,ap, &
    iscale,jlocal,kchem,xuvfac,ef,ec,itail,fmono,emono, &
    tgcm_ncfile,iri90_dir,jmax,glow_ncfile, &
    start_mtime,stop_mtime,data_dir,writelbh,writered

! Execute:

!
! Set up MPI:
!
  call mpi_init(mpierr)
  call mpi_comm_rank(MPI_COMM_WORLD,itask,mpierr)
  call mpi_comm_size(MPI_COMM_WORLD,nproc,mpierr)
!
! Initialize tgcm_ncfile and start/stop times:
! If start_mtime and/or stop_mtime are not read from namelist, model days of -999 
! will flag find_mtimes to find the first and/or last histories on the file.
!
  tgcm_ncfile = ' '
  start_mtime = (/-999,0,0/) 
  stop_mtime  = (/-999,0,0/) 
!
! Root task only:  Read namelist inputs from input file.
! Read times, coordinates, 1D vars, from tgcm history file (tiegcm or timegcm), if provided.
! If tgcm history file is not provided, will use namelist inputs for MSIS/IRI/NOEM:
!
  if (itask == 0) then
    read (5,nml=glow_input)
    if (len_trim(tgcm_ncfile) > 0) then
      tgcm = .true.
      call read_tgcm_coords(tgcm_ncfile)
      ntimes = find_mtimes(tgcm_ncfile,start_mtime,stop_mtime,mtimes,itimes)
    else
      tgcm = .false.
      ntimes = ifix((utstop-utstart)/utstep) + 1
      write(6,"('glow_drv: tgcm_ncfile not provided, will use MSIS/IRI')")
    endif
  endif
!
! Broadcast namelist inputs, tgcm coordinates, times, etc., to all processors:
!
  call mpi_bcast(tgcm,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(data_dir,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(iscale,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(jlocal,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(kchem,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(xuvfac,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(itail,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(ntimes,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  if (tgcm) then
    call mpi_bcast(itimes,mxtimes,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call mpi_bcast(nlon_tgcm,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call mpi_bcast(nlat_tgcm,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    if (itask /= 0) then
      allocate(glon_tgcm(nlon_tgcm))
      allocate(glat_tgcm(nlat_tgcm))
      allocate(eflux(nlon_tgcm,nlat_tgcm))
      allocate(alfa(nlon_tgcm,nlat_tgcm))
    endif
    call mpi_bcast(glon_tgcm,nlon_tgcm,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
    call mpi_bcast(glat_tgcm,nlat_tgcm,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
  else
    call mpi_bcast(nlon_msis,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call mpi_bcast(nlat_msis,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call mpi_bcast(ef,1,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
    call mpi_bcast(ec,1,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
    call mpi_bcast(fmono,1,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
    call mpi_bcast(emono,1,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
  endif
!
! Set grid dimensions:
!
  if (tgcm) then
    jmax=68
    nlev = jmax
    nlat = nlat_tgcm
    nlon = nlon_tgcm
  else
    jmax = 102
    nlev = jmax
    nlat = nlat_msis
    nlon = nlon_msis
  endif
!
! Allocate arrays in other modules (formerly in common blocks):
!
  call cglow_init
  call output_init
!
! Set up grid:
!
  if (tgcm) then
    lon(:) = glon_tgcm(:)
    lat(:) = glat_tgcm(:)
    do k=1,jmax
      lev(k)=float(k)/4.-10.25
    enddo
  else
    do i=1,nlon
      lon(i) = float(i-1)/float(nlon_msis)*360. - 180.
    enddo
    do j=1,nlat
      lat(j) = float(j-1)/float(nlat_msis)*180. - 90. + 90./float(nlat_msis)
    enddo
    lev(:)=0.
  endif
!
! Assign latitude bands to processors:
!
  if (itask == 0 .and. (nlat/nproc)*nproc /= nlat) then 
    write(6,"('glow_drv: number of latitudes must be an integer multiple of number of processors')")
    write(6,"('NLAT =',i3,'   NPROC =',i3)") nlat,nproc
    stop
  endif
  size2d=nlon*nlat
  size3d=nlon*nlat*nlev
  sizeb=nlon*nlat/nproc
  lat0 = 1 + itask * nlat/nproc
  lat1 = (1+itask) * nlat/nproc
!
! Allocate local arrays:
!
  allocate(z(jmax))
  allocate(zun(jmax))
  allocate(zvn(jmax))
  allocate(pedcond(jmax))
  allocate(hallcond(jmax))
  allocate(outf(11,jmax))
  allocate(recbuf(nlon,nlat))
!
! Set electron energy grid:
!
  call egrid (ener, del, nbins)
!
! Time loop:
!
  do itime=1,ntimes
!
! Get input fields if this is a tgcm run, otherwise use namelist inputs to MSIS/IRI/NOEM.
! Only do this section if this is the root task:
!
  if (itask == 0) then 
    if (tgcm) then
      call read_tgcm(tgcm_ncfile,itimes(itime))
      idate=iyear_tgcm*1000+iday_tgcm
      ut   =ut_tgcm*3600.
      f107 =f107_tgcm(itime)
      f107p=f107_tgcm(itime)
      f107a=f107a_tgcm(itime)
      kp=alog((hpower(itime)+4.86)/16.82)/.32
      ap=2.2*exp(.6*kp)
    else
      idate=indate
      ut   =utstart+(itime-1)*utstep
    endif
    write(6,"('glow_drv: idate=',i7,' ut=',f7.1)") idate,ut
    write(6,"('glow_drv: F107=',f5.1,' F107a=',f5.1,' F107p=',f5.1,' Ap=',f5.1)")f107,f107a,f107p,ap
!
! Loop over latitude and longitude to interpolate input fields:
!
    do l=1,nlat
      glat = lat(l)
      do i=1,nlon
        glong = lon(i)
!
! If this is a tgcm run, use altitude grid, fields, and auroral inputs from tgcm history.
! Otherwise, use default altitude grid, MSIS/NOEM/IRI fields, and namelist auroral inputs:
!
        if (tgcm) then
          call tzgrid(i,l,jmax,z,zo,zo2,zn2,zns,znd,zno,ztn,zun,zvn,ze,zti,zte)
        else
          stl = ut/3600. + glong/15.
          if (stl < 0.) stl = stl + 24.
          if (stl >= 24.) stl = stl - 24.
          call mzgrid (jmax,nex,idate,ut,glat,glong,stl,f107a,f107,f107p,ap,iri90_dir, &
                     z,zo,zo2,zn2,zns,znd,zno,ztn,zun,zvn,ze,zti,zte,zxden)
        endif 
!
! Fill global arrays:
!
        zzz(i,l,:)  = z(:) * 1.e5        !km to cm at all altitude levels
        ao(i,l,:)   = zo(:)
        ao2(i,l,:)  = zo2(:)
        an2(i,l,:)  = zn2(:)
        ano(i,l,:)  = zno(:)
        an4s(i,l,:) = zns(:)
        an2d(i,l,:) = znd(:)
        atn(i,l,:)  = ztn(:)
        ati(i,l,:)  = zti(:)
        ate(i,l,:)  = zte(:)
        aun(i,l,:)  = zun(:)
        avn(i,l,:)  = zvn(:)
        ane(i,l,:)  = ze(:)

      enddo    ! longitude loop
    enddo   ! latitude loop
  endif   ! bottom of root-task-only conditional for input/interpolation section
!
! Broadcast idate, ut, f107, f107a, 3D fields, (and 2D auroral fields), to all processors:
!
  call mpi_bcast(idate,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(ut,1,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(f107,1,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(f107a,1,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(zzz,size3d,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(ao,size3d,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(ao2,size3d,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(an2,size3d,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(ano,size3d,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(an4s,size3d,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(an2d,size3d,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(atn,size3d,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(ati,size3d,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(ate,size3d,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(aun,size3d,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(avn,size3d,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(ane,size3d,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
  if (tgcm) then
    call mpi_bcast(eflux,size2d,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
    call mpi_bcast(alfa,size2d,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
  endif
!
! Loop over all subdomain latitudes and longitudes for calls to GLOW:
!
  do l=lat0,lat1
    glat = lat(l)
    do i=1,nlon
      glong = lon(i)
      phitop(:) = 0.
!
! Get auroral inputs from TGCM history or from namelist inputs:
!
      if (tgcm) then
        ef=eflux(i,l)
        ec=alfa(i,l)*1000.       ! keV to eV
        if (ef>.001 .and. ec>1.) call maxt(ef,ec,ener,del,nbins,itail,fmono,emono,phitop)
      else
        if (ef>.001 .and. ec>1.) call maxt (ef,ec,ener,del,nbins,itail,fmono,emono,phitop)
      endif 
!
! Transfer global fields back to altitude arrays at specific lat/lon:
!
      zz(:)  = zzz(i,l,:)
      zo(:)  = ao(i,l,:)
      zo2(:) = ao2(i,l,:)
      zn2(:) = an2(i,l,:)
      zno(:) = ano(i,l,:)
      zns(:) = an4s(i,l,:)
      znd(:) = an2d(i,l,:)
      ztn(:) = atn(i,l,:)
      zti(:) = ati(i,l,:)
      zte(:) = ate(i,l,:)
      zun(:) = aun(i,l,:)
      zvn(:) = avn(i,l,:)
      ze(:)  = ane(i,l,:)
!
! Call GLOW to calculate ionized and excited species, and airglow emission rates:
!
      call glow
!
! Write error code and energy conservation to standard output if it is out of range:
!
      if (ierr > 0) write(6,"('glow_drv: IERR = ',i1)") ierr
      if (abs(efrac) > 0.2) &
        write(6,"('glow_drv:  EFRAC =',f5.2,'  SZA =',f5.2,'  DIP=',f5.2)") efrac,sza,dip
!
! Call CONDUCT to calculate Pederson and Hall conductivities:
!
      do j=1,jmax
        call conduct (glat, glong, z(j), zo(j), zo2(j), zn2(j), &
                      zxden(3,j), zxden(6,j), zxden(7,j), ztn(j), zti(j), zte(j), &
                      pedcond(j), hallcond(j))
      enddo
!
! Collect arrays for output:
!
      do j=1,jmax
        nec(i,l,j)  = ecalc(j)
        ped(i,l,j)  = pedcond(j)
        hall(i,l,j)  = hallcond(j)
        do ix=1,nex
          xden(i,l,j,ix) = zxden(ix,j)
        enddo
        do ix=1,nw
          eta(i,l,j,ix) = zeta(ix,j)
        enddo
        do ix=1,nc
          lbh(i,l,j,ix) = zlbh(ix,j)
          redline(i,l,j,ix)=zceta(ix,5,j)
        enddo
      enddo
    enddo    ! longitude loop
  enddo   ! latitude loop
!
! Gather fields to root task for output:
!
  do j=1,jmax
    call mpi_gather &
         (nec(:,lat0:lat1,j),sizeb,MPI_REAL,recbuf,sizeb,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
    if (itask == 0) nec(:,:,j)=recbuf(:,:)
    call mpi_gather &
         (ped(:,lat0:lat1,j),sizeb,MPI_REAL,recbuf,sizeb,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
    if (itask == 0) ped(:,:,j)=recbuf(:,:)
    call mpi_gather &
         (hall(:,lat0:lat1,j),sizeb,MPI_REAL,recbuf,sizeb,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
    if (itask == 0) hall(:,:,j)=recbuf(:,:)
    do ix=1,nex
      call mpi_gather &
           (xden(:,lat0:lat1,j,ix),sizeb,MPI_REAL,recbuf,sizeb,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
      if (itask == 0) xden(:,:,j,ix) = recbuf(:,:)
    enddo
    do ix=1,nw
      call mpi_gather &
           (eta(:,lat0:lat1,j,ix),sizeb,MPI_REAL,recbuf,sizeb,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
      if (itask == 0) eta(:,:,j,ix) = recbuf(:,:)
    enddo
    do ix=1,nc
      call mpi_gather &
           (lbh(:,lat0:lat1,j,ix),sizeb,MPI_REAL,recbuf,sizeb,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
      if (itask == 0) lbh(:,:,j,ix) = recbuf(:,:)
      call mpi_gather &
           (redline(:,lat0:lat1,j,ix),sizeb,MPI_REAL,recbuf,sizeb,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
      if (itask == 0) redline(:,:,j,ix) = recbuf(:,:)
    enddo
  enddo
!
! Output section:
! Create and define a new netCDF output file for each time (root task only):
!
  if (itask == 0) then 
    write (ifile,"('.',i3.3,'.nc')"),itime
    glow_ncfileit = trim(glow_ncfile) // ifile
    call create_ncfile(glow_ncfileit,tgcm_ncfile)
    call write_ncfile(glow_ncfileit)
  endif

  enddo ! bottom of time loop

  call mpi_finalize(mpierr)

end program glowdriver
