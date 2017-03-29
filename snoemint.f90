! Subroutine SNOEMINT gets NO estimate from the NOEM emperical model and
! INTerpolates it onto an altitude grid.  Extrapolation is done above 150
! km assuming a scale height approximation, and below 100 km
! assuming a constant number density profile.

! Stan Solomon, 12/2014
! Refactored to f90, scs, 6/2016
! Fixed bug in Kp (estimated from Ap) so that xkp >= 0, scs, 1/2017

! Input:
!   IDATE  Date in yyddd or yyyyddd format
!   GLAT   Geographic latitude in degrees
!   GLONG  Geographic longitude in degrees
!   F107   10.7 cm radio flux index
!   AP     Ap index
!   JMAX   Number of points in altitude grid
!   Z      Altitude grid in km
!   ZTN    Temperature at Z in K
! Output:
!   ZNO    Nitric oxide density at Z in cm-3


    subroutine snoemint(idate,glat,glong,f107,ap,jmax,z,ztn,zno)

      implicit none

      integer,intent(in) :: idate, jmax
      real,intent(in) :: glat, glong, f107, ap, z(jmax), ztn(jmax)
      real,intent(out) :: zno(jmax)

      real,parameter :: pi=3.1415926536

      integer :: iday, klat1, klat2, kz1, kz2, j
      real :: zg(16), xmlatno(33), zmno(33,16), zmnoi(16)
      real :: xmlong, xmlat, xkp, rat, h

! Find magnetic latitude:
 
      call geomag(0,glong,glat,xmlong,xmlat)
 
! Get zonal mean no profiles:
 
      iday=idate-idate/1000*1000
      xkp=1.75*alog(0.4*ap)
      if (xkp < 0.0) xkp=0.0
      call snoem(iday,xkp,f107,zg,xmlatno,zmno)
 
! Interpolate altitude profile at magnetic latitude:
 
      klat1=ifix(xmlat+80.)/5+1
      klat2=klat1+1
      if (klat1 < 1) klat1=1
      if (klat1 > 33) klat1=33
      if (klat2 < 1) klat1=1
      if (klat2 > 33) klat2=33
      rat=xmlat/5.-ifix(xmlat)/5

      do j=1,16
        zmnoi(j) = alog(zmno(klat1,j)*(1.-rat)+zmno(klat2,j)*rat)
      end do

! Interpolate onto altitude grid:
! Use constant value below 100 km and scale height assumption above 150 km:

      h=0.03*ztn(jmax)
      do j=1,jmax
        if (z(j) <= 100.) zno(j)=exp(zmnoi(16))
        if (z(j) > 100. .and. z(j) <= 150.) then
          kz2=ifix((150.-z(j))*.3)+1
          kz1=kz2+1
          zno(j)=exp(zmnoi(kz1) + (zmnoi(kz2)-zmnoi(kz1)) * (z(j)-zg(kz1)) / (zg(kz2)-zg(kz1)))
        endif
        if (z(j) > 150.) zno(j)=exp(zmnoi(1)+(150.-z(j))/h)
      end do

      return

    end
