module output

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Ben Foster and Stan Solomon, 1/15
! Stan Solomon, 1/16

! Handles output from the GLOW model to netCDF files.
! For usage, see program glowdriver.f90.
 
! The GLOW grid dimensions and coordinates must be set before
! any routines in this module are called. They are currently
! set by the main driver (glowdriver.f90) and module cglow.f90.


  use netcdf

  use cglow,only: nbins,lmax,nmaj,nei,nex,nw,nc,nst
  use cglow,only: idate,ut,f107a,f107,f107p,ap
  use cglow,only: iscale,jlocal,kchem,xuvfac

  implicit none

  integer :: nlon,nlat,nlev                       ! dimensions
  real,allocatable,save :: lon(:),lat(:),lev(:)   ! coordinates 
  logical :: writelbh,writered                    ! switches
!
! Output arrays:
!
  real,allocatable :: zzz(:,:,:)    ! geometric altitude (nlon,nlat,nlev)
  real,allocatable :: ao(:,:,:)     ! atomic oxygen
  real,allocatable :: ao2(:,:,:)    ! molecular oxygen
  real,allocatable :: an2(:,:,:)    ! molecular nitrogen
  real,allocatable :: ano(:,:,:)    ! nitric oxide
  real,allocatable :: an4s(:,:,:)   ! N(4S)
  real,allocatable :: an2d(:,:,:)   ! N(2D)
  real,allocatable :: atn(:,:,:)    ! neutral temperature
  real,allocatable :: ati(:,:,:)    ! ion temperature
  real,allocatable :: ate(:,:,:)    ! electron temperature
  real,allocatable :: aun(:,:,:)    ! zonal neutral wind
  real,allocatable :: avn(:,:,:)    ! meridional neutral wind
  real,allocatable :: ane(:,:,:)    ! electron density from input model
  real,allocatable :: nec(:,:,:)    ! electron density calculated by glow
  real,allocatable :: ped(:,:,:)    ! Pederson conductivity calculated by glow (not output)
  real,allocatable :: hall(:,:,:)   ! Hall conductivity calculated by glow (not output)
  real,allocatable :: xden(:,:,:,:) ! ionized and excited states (nlon,nlat,nlev,nex)
  real,allocatable :: eta(:,:,:,:)  ! airglow emission rates (nlon,nlat,nlev,nw)
  real,allocatable :: lbh(:,:,:,:)  ! lbh (v'=0-6) excitation rates (nlon,nlat,nlev,nc)
  real,allocatable :: redline(:,:,:,:)  ! OI 6300 emission components (nlon,nlat,nlev,nc)
!
! Netcdf file unit:
!
  integer :: ncid

  contains

!-----------------------------------------------------------------------

  subroutine output_init
!
! Allocate output arrays (called from glowdriver.F90).
!
  allocate(lon(nlon))
  allocate(lat(nlat))
  allocate(lev(nlev))
  allocate(zzz(nlon,nlat,nlev))
  allocate(ao(nlon,nlat,nlev))
  allocate(ao2(nlon,nlat,nlev))
  allocate(an2(nlon,nlat,nlev))
  allocate(ano(nlon,nlat,nlev))
  allocate(an4s(nlon,nlat,nlev))
  allocate(an2d(nlon,nlat,nlev))
  allocate(atn(nlon,nlat,nlev))
  allocate(ati(nlon,nlat,nlev))
  allocate(ate(nlon,nlat,nlev))
  allocate(aun(nlon,nlat,nlev))
  allocate(avn(nlon,nlat,nlev))
  allocate(ane(nlon,nlat,nlev))
  allocate(nec(nlon,nlat,nlev))
  allocate(ped(nlon,nlat,nlev))
  allocate(hall(nlon,nlat,nlev))
  allocate(xden(nlon,nlat,nlev,nex))
  allocate(eta(nlon,nlat,nlev,nw))
  allocate(lbh(nlon,nlat,nlev,nc))
  allocate(redline(nlon,nlat,nlev,nc))

  end subroutine output_init

!-----------------------------------------------------------------------

  subroutine create_ncfile(ncfile,tgcm_ncfile)
!
! Create and define a netcdf glow output file, with dimensions, attributes, 
! variables, etc., but do not write data yet. 
!
! Args:
  character(len=*),intent(in) :: ncfile,tgcm_ncfile
!
! Local:
  integer :: istat,id,idv
  integer :: id_lon, id_lat, id_lev, id_nex, id_nw, id_nc
  integer :: idv_lon,idv_lat,idv_lev
  integer :: idv_idate,idv_ut,idv_f107a,idv_f107,idv_f107p,idv_ap
  integer :: idv_iscale,idv_jlocal,idv_kchem,idv_xuvfac
  character(len=16) :: logname
  character(len=24) :: create_date,create_time
  character(len=64) :: create_datetime

  istat = nf90_create(ncfile,NF90_CLOBBER,ncid) 
  if (istat /= NF90_NOERR) then
    write(6,"('>>> Error creating netcdf output file ',a)") trim(ncfile)
    call handle_ncerr(istat,'Error from nf90_create',1)
  endif
  write(6,"('Created netcdf dataset ',a)") trim(ncfile)
!
! Define dimensions:
!
  istat = nf90_def_dim(ncid,'lon',nlon,id_lon)
  istat = nf90_def_dim(ncid,'lat',nlat,id_lat)
  istat = nf90_def_dim(ncid,'lev',nlev,id_lev)
  istat = nf90_def_dim(ncid,'nex',nex,id_nex)
  istat = nf90_def_dim(ncid,'nw',nw,id_nw)
  istat = nf90_def_dim(ncid,'nc',nc,id_nc)
!
! Define scalars:
!
  istat = nf90_def_var(ncid,'IDATE',NF90_INT,idv_idate)
  istat = nf90_put_att(ncid,idv_idate,'long_name','date in yyyyddd format')
  istat = nf90_put_att(ncid,idv_idate,'units','days')

  istat = nf90_def_var(ncid,'UT',NF90_FLOAT,idv_ut)
  istat = nf90_put_att(ncid,idv_ut,'long_name','universal time')
  istat = nf90_put_att(ncid,idv_ut,'units','hours')

  istat = nf90_def_var(ncid,'F107A',NF90_FLOAT,idv_f107a)
  istat = nf90_put_att(ncid,idv_f107a,'long_name','F10.7 index 81-day centered average')
  istat = nf90_put_att(ncid,idv_f107a,'units','sfu')

  istat = nf90_def_var(ncid,'F107',NF90_FLOAT,idv_f107)
  istat = nf90_put_att(ncid,idv_f107,'long_name','F10.7 index')
  istat = nf90_put_att(ncid,idv_f107,'units','sfu')

  istat = nf90_def_var(ncid,'F107P',NF90_FLOAT,idv_f107p)
  istat = nf90_put_att(ncid,idv_f107p,'long_name','F10.7 index for previous day')
  istat = nf90_put_att(ncid,idv_f107p,'units','sfu')

  istat = nf90_def_var(ncid,'AP',NF90_FLOAT,idv_ap)
  istat = nf90_put_att(ncid,idv_ap,'long_name','Ap index (daily)')
  istat = nf90_put_att(ncid,idv_ap,'units','none')

  istat = nf90_def_var(ncid,'ISCALE',NF90_INT,idv_iscale)
  istat = nf90_put_att(ncid,idv_iscale,'long_name','Solar flux scaling method')
  istat = nf90_put_att(ncid,idv_iscale,'units','none')

  istat = nf90_def_var(ncid,'JLOCAL',NF90_INT,idv_jlocal)
  istat = nf90_put_att(ncid,idv_jlocal,'long_name','Electron transport local switch')
  istat = nf90_put_att(ncid,idv_jlocal,'units','none')

  istat = nf90_def_var(ncid,'KCHEM',NF90_INT,idv_kchem)
  istat = nf90_put_att(ncid,idv_kchem,'long_name','Ion chemistry selction switch')
  istat = nf90_put_att(ncid,idv_kchem,'units','none')

  istat = nf90_def_var(ncid,'XUVFAC',NF90_FLOAT,idv_xuvfac)
  istat = nf90_put_att(ncid,idv_xuvfac,'long_name','XUV scaling factor')
  istat = nf90_put_att(ncid,idv_xuvfac,'units','none')
!
! Define coordinate variables):
!
  istat = nf90_def_var(ncid,'lon',NF90_FLOAT,(/id_lon/),idv_lon)
  istat = nf90_put_att(ncid,idv_lon,'long_name',&
          'geographic longitude (-west, +east)')
  istat = nf90_put_att(ncid,idv_lon,'units','degrees_east')

  istat = nf90_def_var(ncid,'lat',NF90_FLOAT,(/id_lat/),idv_lat)
  istat = nf90_put_att(ncid,idv_lat,'long_name',&
          'geographic latitude (-south, +north)')
  istat = nf90_put_att(ncid,idv_lat,'units','degrees_north')

  istat = nf90_def_var(ncid,'lev',NF90_FLOAT,(/id_lev/),idv_lev)
  istat = nf90_put_att(ncid,idv_lev,'long_name',&
          'pressure/altitude level number')
  istat = nf90_put_att(ncid,idv_lat,'units','dimensionless')
!
! Define data variable arrays:
!
  istat = nf90_def_var(ncid,'ZZZ',NF90_FLOAT,(/id_lon,id_lat,id_lev/),idv)
  istat = nf90_put_att(ncid,idv,'long_name','geometric altitude')
  istat = nf90_put_att(ncid,idv,'units','cm')

  istat = nf90_def_var(ncid,'AO',NF90_FLOAT,(/id_lon,id_lat,id_lev/),idv)
  istat = nf90_put_att(ncid,idv,'long_name','atomic oxygen')
  istat = nf90_put_att(ncid,idv,'units','cm-3')

  istat = nf90_def_var(ncid,'AO2',NF90_FLOAT,(/id_lon,id_lat,id_lev/),idv)
  istat = nf90_put_att(ncid,idv,'long_name','molecular oxygen')
  istat = nf90_put_att(ncid,idv,'units','cm-3')

  istat = nf90_def_var(ncid,'AN2',NF90_FLOAT,(/id_lon,id_lat,id_lev/),idv)
  istat = nf90_put_att(ncid,idv,'long_name','molecular nitrogen')
  istat = nf90_put_att(ncid,idv,'units','cm-3')

  istat = nf90_def_var(ncid,'ANO',NF90_FLOAT,(/id_lon,id_lat,id_lev/),idv)
  istat = nf90_put_att(ncid,idv,'long_name','nitric oxide')
  istat = nf90_put_att(ncid,idv,'units','cm-3')

  istat = nf90_def_var(ncid,'AN4S',NF90_FLOAT,(/id_lon,id_lat,id_lev/),idv)
  istat = nf90_put_att(ncid,idv,'long_name','atomic nitrogen')
  istat = nf90_put_att(ncid,idv,'units','cm-3')

  istat = nf90_def_var(ncid,'AN2D',NF90_FLOAT,(/id_lon,id_lat,id_lev/),idv)
  istat = nf90_put_att(ncid,idv,'long_name','atomic nitrogen doublet-D')
  istat = nf90_put_att(ncid,idv,'units','cm-3')

  istat = nf90_def_var(ncid,'ATN',NF90_FLOAT,(/id_lon,id_lat,id_lev/),idv)
  istat = nf90_put_att(ncid,idv,'long_name','netural temperature')
  istat = nf90_put_att(ncid,idv,'units','K')

  istat = nf90_def_var(ncid,'ATI',NF90_FLOAT,(/id_lon,id_lat,id_lev/),idv)
  istat = nf90_put_att(ncid,idv,'long_name','ion temperature')
  istat = nf90_put_att(ncid,idv,'units','K')

  istat = nf90_def_var(ncid,'ATE',NF90_FLOAT,(/id_lon,id_lat,id_lev/),idv)
  istat = nf90_put_att(ncid,idv,'long_name','electron temperature')
  istat = nf90_put_att(ncid,idv,'units','K')

  istat = nf90_def_var(ncid,'AUN',NF90_FLOAT,(/id_lon,id_lat,id_lev/),idv)
  istat = nf90_put_att(ncid,idv,'long_name','zonal neutral wind')
  istat = nf90_put_att(ncid,idv,'units','cm/s')

  istat = nf90_def_var(ncid,'AVN',NF90_FLOAT,(/id_lon,id_lat,id_lev/),idv)
  istat = nf90_put_att(ncid,idv,'long_name','meridional neutral wind')
  istat = nf90_put_att(ncid,idv,'units','cm/s')

  istat = nf90_def_var(ncid,'ANE',NF90_FLOAT,(/id_lon,id_lat,id_lev/),idv)
  istat = nf90_put_att(ncid,idv,'long_name','electron density input')
  istat = nf90_put_att(ncid,idv,'units','cm-3')

  istat = nf90_def_var(ncid,'NEC',NF90_FLOAT,(/id_lon,id_lat,id_lev/),idv)
  istat = nf90_put_att(ncid,idv,'long_name','electron density calculated')
  istat = nf90_put_att(ncid,idv,'units','cm-3')

  istat = nf90_def_var(ncid,'XDEN',NF90_FLOAT,(/id_lon,id_lat,id_lev,id_nex/),idv)
  istat = nf90_put_att(ncid,idv,'long_name','excited and ionized state densities')
  istat = nf90_put_att(ncid,idv,'units','cm-3')

  istat = nf90_def_var(ncid,'ETA',NF90_FLOAT,(/id_lon,id_lat,id_lev,id_nw/),idv)
  istat = nf90_put_att(ncid,idv,'long_name','volume emission rates')
  istat = nf90_put_att(ncid,idv,'units','cm-3 s-1')

  if (writelbh) then
    istat = nf90_def_var(ncid,'LBH',NF90_FLOAT,(/id_lon,id_lat,id_lev,id_nc/),idv)
    istat = nf90_put_att(ncid,idv,'long_name','LBH (v=0-6) excitation rates')
    istat = nf90_put_att(ncid,idv,'units','cm-3 s-1')
  endif

  if (writered) then
    istat = nf90_def_var(ncid,'REDLINE',NF90_FLOAT,(/id_lon,id_lat,id_lev,id_nc/),idv)
    istat = nf90_put_att(ncid,idv,'long_name','OI 6300 emission components')
    istat = nf90_put_att(ncid,idv,'units','cm-3 s-1')
  endif
!
! Global file attributes:
!
! Create date and time:
!
  call datetime(create_date,create_time)
  create_datetime = trim(create_date)//' '//trim(create_time)
  istat = nf90_put_att(ncid,NF90_GLOBAL,"create_date",trim(create_datetime))
!
! Log name (author):
!
  call getenv('LOGNAME',logname)
  istat = nf90_put_att(ncid,NF90_GLOBAL,'logname',trim(logname))
!
! tgcm history file, if one has been read:
!
  if (len_trim(tgcm_ncfile) > 0) then
    istat = nf90_put_att(ncid,NF90_GLOBAL,'tgcm_ncfile',trim(tgcm_ncfile))
  else
    istat = nf90_put_att(ncid,NF90_GLOBAL,'tgcm_ncfile','[none]')
  endif
!
! Take out of define mode and go into data mode:
!
  istat = nf90_enddef(ncid)
!
! Write scalars:
!
  istat = nf90_put_var(ncid,idv_idate,idate)
  istat = nf90_put_var(ncid,idv_ut,ut)
  istat = nf90_put_var(ncid,idv_f107a,f107a)
  istat = nf90_put_var(ncid,idv_f107,f107)
  istat = nf90_put_var(ncid,idv_f107p,f107p)
  istat = nf90_put_var(ncid,idv_ap,ap)
  istat = nf90_put_var(ncid,idv_iscale,iscale)
  istat = nf90_put_var(ncid,idv_jlocal,jlocal)
  istat = nf90_put_var(ncid,idv_kchem,kchem)
  istat = nf90_put_var(ncid,idv_xuvfac,xuvfac)
!
! Write coordinate values:
!
  istat = nf90_put_var(ncid,idv_lon,lon)
  istat = nf90_put_var(ncid,idv_lat,lat)
  istat = nf90_put_var(ncid,idv_lev,lev)
!
! Close the file:
!
  istat = nf90_close(ncid)

  end subroutine create_ncfile

!-----------------------------------------------------------------------

  subroutine write_ncfile(ncfile)
!
! Write data to netcdf glow output file at current time. The file has 
! already been created w/ dimensions by sub create_ncfile.
!
! Args:
  character(len=*),intent(in) :: ncfile
!
! Local:
  integer :: istat,id

  istat = nf90_open(ncfile,NF90_WRITE,ncid)
  if (istat /= NF90_NOERR) then
    write(6,"('>>> Error opening file for writing: ncfile=',a)") trim(ncfile)
    call handle_ncerr(istat,'Error from nf90_open',1)
  endif

  istat = nf90_inq_varid(ncid,'ZZZ',id)
  istat = nf90_put_var(ncid,id,zzz)
  istat = nf90_inq_varid(ncid,'AO',id)
  istat = nf90_put_var(ncid,id,ao)
  istat = nf90_inq_varid(ncid,'AO2',id)
  istat = nf90_put_var(ncid,id,ao2)
  istat = nf90_inq_varid(ncid,'AN2',id)
  istat = nf90_put_var(ncid,id,an2)
  istat = nf90_inq_varid(ncid,'ANO',id)
  istat = nf90_put_var(ncid,id,ano)
  istat = nf90_inq_varid(ncid,'AN4S',id)
  istat = nf90_put_var(ncid,id,an4s)
  istat = nf90_inq_varid(ncid,'AN2D',id)
  istat = nf90_put_var(ncid,id,an2d)
  istat = nf90_inq_varid(ncid,'ATN',id)
  istat = nf90_put_var(ncid,id,atn)
  istat = nf90_inq_varid(ncid,'ATI',id)
  istat = nf90_put_var(ncid,id,ati)
  istat = nf90_inq_varid(ncid,'ATE',id)
  istat = nf90_put_var(ncid,id,ate)
  istat = nf90_inq_varid(ncid,'AUN',id)
  istat = nf90_put_var(ncid,id,aun)
  istat = nf90_inq_varid(ncid,'AVN',id)
  istat = nf90_put_var(ncid,id,avn)
  istat = nf90_inq_varid(ncid,'ANE',id)
  istat = nf90_put_var(ncid,id,ane)
  istat = nf90_inq_varid(ncid,'NEC',id)
  istat = nf90_put_var(ncid,id,nec)
  istat = nf90_inq_varid(ncid,'XDEN',id)
  istat = nf90_put_var(ncid,id,xden)
  istat = nf90_inq_varid(ncid,'ETA',id)
  istat = nf90_put_var(ncid,id,eta)
  if (writelbh) then
    istat = nf90_inq_varid(ncid,'LBH',id)
    istat = nf90_put_var(ncid,id,lbh)
  endif
  if (writered) then
    istat = nf90_inq_varid(ncid,'REDLINE',id)
    istat = nf90_put_var(ncid,id,redline)
  endif

  istat = nf90_close(ncid)

  end subroutine write_ncfile

!-----------------------------------------------------------------------

      subroutine datetime(curdate,curtime)
!
! Return character*8 values for current date and time.
! (sub date_and_time is an f90 intrinsic)
!
      implicit none
!
! Args:
      character(len=*),intent(out) :: curdate,curtime
!
! Local:
      character(len=8) :: date
      character(len=10) :: time
      character(len=5) :: zone
      integer :: values(8)
!
      curdate = ' '
      curtime = ' '
      call date_and_time(date,time,zone,values)
!
!     write(6,"('datetime: date=',a,' time=',a,' zone=',a)")
!    |  date,time,zone
!     write(6,"('datetime: values=',8i8)") values
!
      curdate(1:2) = date(5:6)
      curdate(3:3) = '/'
      curdate(4:5) = date(7:8)
      curdate(6:6) = '/'
      curdate(7:8) = date(3:4)
!
      curtime(1:2) = time(1:2)
      curtime(3:3) = ':'
      curtime(4:5) = time(3:4)
      curtime(6:6) = ':'
      curtime(7:8) = time(5:6)
!
      end subroutine datetime

!-----------------------------------------------------------------------

  subroutine handle_ncerr(istat,msg,ifatal)
!
! Handle a netcdf lib error:
!
    integer,intent(in) :: istat,ifatal
    character(len=*),intent(in) :: msg
!
    write(6,"(/72('-'))")
    write(6,"('>>> Error from netcdf library:')")
    write(6,"(a)") trim(msg)
    write(6,"('istat=',i5)") istat
    write(6,"(a)") nf90_strerror(istat)
    write(6,"(72('-')/)")
    if (ifatal > 0) stop('Fatal netcdf error')
  end subroutine handle_ncerr

!-----------------------------------------------------------------------

end module output
