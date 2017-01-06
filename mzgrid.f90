! Subroutine MZGRID gets fields from eMpirical models on default GLOW altitude grid.

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Stan Solomon, 12/15, 1/16
! Extracted from glowdriver.f90 into separate file mzgrid.f90, SCS, 12/16

! Neutral densities from NRLMSISE-00 model (a.k.a. MSIS2K).
! Nitric oxide densities from NOEM model (via snoem.f and snoemint.f).
! Electron densities from IRI-90.

! Inputs:
!         jmax   Number of altitude levels used by GLOW (should = 86 for MSIS/IRI/NOEM runs)
!         nex    Number of ionized/excited species (for array zxden)
!         idate  Nate in yyyyddd or yyddd format
!         ut     Universal time, seconds)
!         glat   Latitude, degrees
!         glong  Longitude, degrees
!         stl    Local solar time, hours
!         f107a  F10.7 index 81-day centered average
!         f107   F10.7 index of day
!         f107p  F107. index of previous day
!         ap     Ap index daily value
!         iri90_dir  Directory containing IRI input files (set in namelist inputs)
          
! Outputs:
!         z      Geographic altitude (km)
!         zo     O number density, cm-3
!         zo2    O2    "
!         zn2    N2    "
!         zns    N(4S) "
!         n2d    N(2D) "      (set to zero since this is calculated by GLOW)
!         no     NO    "
!         ztn    Tn, K
!         zun    Zonal wind velocity, cm/s, currently = 0 (not used by GLOW)
!         zvn    Meridional wind velocity, cm/s currently = 0 (not used by GLOW)
!         zti    Ti, K
!         zte    Te, K
!         zxden  Array of ionized/excited species density, cm-3, must be dimensioned (nex,jmax)

subroutine mzgrid (jmax,nex,idate,ut,glat,glong,stl,f107a,f107,f107p,ap,iri90_dir, &
                   z,zo,zo2,zn2,zns,znd,zno,ztn,zun,zvn,ze,zti,zte,zxden)

  implicit none

  integer,intent(in) :: jmax,nex,idate
  real,intent(in) :: ut,glat,glong,stl,f107a,f107,f107p,ap
  character(len=1024),intent(in) :: iri90_dir
  real,intent(out) :: z(jmax),zo(jmax),zo2(jmax),zn2(jmax),zns(jmax),znd(jmax), &
       zno(jmax),ztn(jmax),zti(jmax),zte(jmax),zun(jmax),zvn(jmax),ze(jmax),zxden(nex,jmax)

  integer :: j,ijf,jmag,iday,mmdd
  real :: rz12, d(8), t(2), sw(25), oarr(30)
  logical :: jf(12)
  real,allocatable :: outf(:,:)              ! iri output (11,jmax)
  data sw/25*1./

  if (jmax /= 86) then
    write(6,"('mzgrid: unknown JMAX = ',i5)") jmax
    stop 'mzgrid'
  endif

  allocate(outf(11,jmax))

!
! Set default altitudes:
!
      z = (/ 80.0,81.5,83.0,84.5,86.0,87.5,89.0,90.5,92.0,93.5, &
             95.0,97.5,99.0,100.5,102.0,103.5,105.0,107.5,109.,111., &
             113.,115.,117.5,120.,122.5,125.,128.,131.,135.,139., &
             143.,148.,153.,158.,164.,170.,176.,183.,190.,197., &
             205.,213.,221.,229.,237.,245.,254.,263.,272.,281., &
             290.,300.,310.,320.,330.,340.,350.,360.,370.,380., &
             390.,400.,410.,420.,430.,440.,450.,460.,470.,480., &
             490.,500.,510.,520.,530.,540.,550.,560.,570.,580., &
             590.,600.,610.,620.,630.,640. /)
!
! Call MSIS-2K to get neutral densities and temperature:
!
        call tselec(sw)

        do j=1,jmax ! levels
          call gtd7(idate,ut,z(j),glat,glong,stl,f107a,f107p,ap,48,d,t)
          zo(j) = d(2)
          zn2(j) = d(3)
          zo2(j) = d(4)
          zns(j) = d(8)
          ztn(j) = t(2)
          znd(j)  = 0.
        enddo
!
! Call SNOEMINT to obtain NO profile from the Nitric Oxide Empirical Model (NOEM):
!
        call snoemint(idate,glat,glong,f107,ap,jmax,z,ztn,zno)
!
! Call International Reference Ionosphere-1990 subroutine to get
! electron density, electron temperature, and ion temperature:
! The directory iri90_dir is the location of the ccirnn.asc and ursinn.asc files.
!
        jf(:) = .true.
!       jf(12) = .false.
        jf(5) = .false.
        jmag = 0
        rz12 = -f107a
        iday = idate - idate/1000*1000
        mmdd = -iday
        outf = 0.

        call iri90(jf,jmag,glat,glong,rz12,mmdd,stl,z,jmax,trim(iri90_dir),outf,oarr)

        do j=1,jmax
          ze(j) = outf(1,j) / 1.E6
          if (ze(j) .lt. 100.) ze(j) = 100.
          zti(j) = outf(3,j)
          if (zti(j) .lt. ztn(j)) zti(j) = ztn(j)
          zte(j) = outf(4,j)
          if (zte(j) .lt. ztn(j)) zte(j) = ztn(j)
          zxden(3,j) = ze(j) * outf(5,j)/100.
          zxden(6,j) = ze(j) * outf(8,j)/100.
          zxden(7,j) = ze(j) * outf(9,j)/100.
        enddo
!
! Until implementation of an empirical neutral wind model, winds are set to zero:
!
        zun(:) = 0.
        zvn(:) = 0.

end subroutine mzgrid
