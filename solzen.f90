! Subroutine SOLZEN

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Stan Solomon, 1988.
! Temporary Y2K fix-up, SCS, 2005.
! Refactored to f90, SCC, 2016.

! Returns Solar Zenith Angle SZA in degrees for specified date in form yyddd or yyyyddd,
! universal time in seconds, geographic latitude and longitude in degrees.

subroutine solzen (idate, ut, glat, glong, sza)

  implicit none

  integer,intent(in) :: idate
  real,intent(in) :: ut, glat, glong
  real,intent(out) :: sza

  real,parameter :: pi=3.1415926536
  real :: rlat, rlong, sdec, srasn, gst, rh, cossza

  rlat = glat * pi/180.
  rlong = glong * pi/180.
  call suncor (idate, ut, sdec, srasn, gst)
  rh = srasn - (gst+rlong)
  cossza = sin(sdec)*sin(rlat) + cos(sdec)*cos(rlat)*cos(rh)
  sza = acos(cossza) * 180./pi
  return

end subroutine solzen


! Subroutine SUNCOR returns the declination SDEC and right ascension
! SRASN of the sun in GEI coordinates, radians, for a given date IDATE
! in yyddd or yyyyddd format, universal time UT in seconds, and Greenwich sidereal
! time GST in radians.  Reference:  C.T. Russell, Geophysical Coordinate Transforms.

subroutine suncor (idate, ut, sdec, srasn, gst)

  implicit none

  integer,intent(in) :: idate
  real,intent(in) :: ut
  real,intent(out) :: sdec, srasn, gst

  real,parameter :: pi=3.1415926536
  real :: fday, dj, t, vl, g, slong, obliq, slp, sind, cosd
  integer :: iyr, iday

  fday=ut/86400.
  iyr=idate/1000
  iday=idate-iyr*1000

! Temporary Y2K fix-up:
! Should work with either yyddd or yyyyddd format from 1950 to 2050.
! Note deteriorating accuracy after ~2050 anyway.
! Won't work after 2100 due to lack of a leap year.

  if (iyr >= 1900) iyr=iyr-1900
  if (iyr < 50) iyr=iyr+100

  dj=365*iyr+(iyr-1)/4+iday+fday-0.5
  t=dj/36525.
  vl=amod(279.696678+.9856473354*dj,360.)
  gst=amod(279.696678+.9856473354*dj+360.*fday+180.,360.) * pi/180.
  g=amod(358.475845+.985600267*dj,360.) * pi/180.
  slong=vl+(1.91946-.004789*t)*sin(g)+.020094*sin(2.*g)
  obliq=(23.45229-0.0130125*t) *pi/180.
  slp=(slong-.005686) * pi/180.
  sind=sin(obliq)*sin(slp)
  cosd=sqrt(1.-sind**2)
  sdec=atan(sind/cosd)
  srasn=pi-atan2(1./tan(obliq)*sind/cosd,-cos(slp)/cosd)

  return

end subroutine suncor
