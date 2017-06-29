! Subroutine SNOEM calculates nitric oxide zonal mean altitude profile
! as function of magnetic latitude for specified day of year, Kp, and F10.7.

! The NOEM empirical model is based on data from the SNOE ultraviolet
! spectrometer during 1998-2000, using empirical orthogonal function analysis.
! Altitude range is from 100 to 150 km.

! Marsh et al., JGR, 109, A07301, doi:10.1029/2003JA010199, 2004.

! Adapted by Stan Solomon, 5/2014, from IDL and F90 code supplied by Dan Marsh. 
! Minor revisions to make compatible with gfortran, SCS, 4/2017

    subroutine snoem(doy, kp, f107, z, mlat, nozm)

      use cglow,only: data_dir

      implicit none

      integer,intent(in) :: doy
      real,intent(in) :: kp, f107
      real,intent(out) :: z(16), mlat(33), nozm(33,16)

      real, save :: zin(16)          ! altitude grid
      real, save :: mlatin(33)       ! magnetic latitude grid
      real, save :: no_mean(33,16)   ! mean nitric oxide distribution
      real, save :: eofs(33,16,3)    ! empirical orthogonal functions
      real :: theta0                 ! day number in degrees
      real :: dec                    ! solar declination angle
      real :: m1, m2, m3             ! coefficients for first 3 eofs
      real, parameter :: pi=3.1415926536
      integer, save :: ifirst=1
      integer :: j, k, n
      character(len=1024) :: filepath 

!... read eof file on first call

      if (ifirst == 1) then
        ifirst = 0
        filepath = trim(data_dir)//'snoem_eof.dat'
        open(unit=1,file=filepath,status='old',action='read')
        read(1,*) (zin(k),k=1,16)
        read(1,*) (mlatin(j),j=1,33)
        read(1,*) ((no_mean(j,k),j=1,33),k=1,16)
        read(1,*) (((eofs(j,k,n),j=1,33),k=1,16),n=1,3)
        close(unit=1)
      endif

!... calculate coefficients (m1 to m3) for eofs based on geophysical parameters

!... eof1 - kp 

      m1 =  kp * 0.689254 - 1.53366

!... eof2 - declination

      theta0 = 2.*pi * float(doy - 1) / 365.

      dec = 0.006918 &
          - 0.399912 * cos(theta0)   + 0.070257 * sin(theta0) &
          - 0.006758 * cos(2*theta0) + 0.000907 * sin(2*theta0) & 
          - 0.002697 * cos(3*theta0) + 0.001480 * sin(3*theta0)

      dec = dec * 180./pi

      m2 = -0.31978 + dec*0.097309 + dec**2*0.00048979 - dec**3*0.00010360
      
!... eof3 - f107 

      m3 =  alog10(f107) * 6.35777 - 13.8163 

!... zonal mean distrib. is sum of mean and eofs

      do k=1,16
        do j=1,33
          nozm(j,k) = no_mean(j,k)-m1*eofs(j,k,1)+m2*eofs(j,k,2)-m3*eofs(j,k,3) 
        end do
      end do

      z(:) = zin(:)
      mlat(:) = mlatin(:)

      return

    end subroutine snoem
