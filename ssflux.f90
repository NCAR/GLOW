! Subroutine SSFLUX
!
! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.
!
! Subroutine SSFLUX calculates the solar EUV and FUV flux in the range
! 0.5 to 1750 Angstroms for a specified level of solar activity.
!
! The calling routine supplies a scaling switch ISCALE, the daily 10.7
! cm flux F107, its 81-day centered average F107A, and an XUV
! enhancement factor XUVFAC.  XUVFAC is ignored if set to zero.
!
! XUVFAC is applied from 18-250 A for the Hinteregger model (ISCALE=0),
! from 18-50 A for the EUVAC model (ISCALE=1), and not at all for
! user-supplied data (ISCALE=2)
!
! The subroutine returns the longwave boundary WAVE1 and shortwave
! boundary WAVE2 of the wavelenth bins, and the solar flux in each bin
! SFLUX.  Bins are arranged in energy multiples of 2 from 0.5 to 8 A,
! aligned with k-shell boundaries at 23, 32, and 44 A from 18 to 44 nm,
! 10 A in width from 60 to 1050 A, and 50 A in width from 1050 to
! 1750 A with the exception of Lyman-alpha which has its own bin from
! 1210 to 1220 A. 
!
! Methods used:
!   If ISCALE=0 the flux is scaled using parameterization methods based
! on F107 and F107A.  For ionizing EUV, Hinteregger's contrast ratio
! method (Hinteregger et al., GRL, 8, 1147, 1981) is used, based on the
! reference spectrum SC#21REFW re-binned at 1 nm bin resolution.
! Enhancement ratios for H Ly-B and Fe XVI are calculated
! from F107 and F107A using Hinteregger's formula, employing
! coefficients which reduce to the reference values at F107=67.6,
! F107A=71.5.  The 'best fit' coefficients are not used as they produce
! some negative values at low solar activity, but remain in a
! 'commented out' data statement for reference.  The EUV spectrum is
! then scaled from these modeled ratios, using contrast ratios from SC#21REFW.
!   If ISCALE=1, the EUV flux (50-1050A) is scaled using the EUVAC model
! (Richards et al., JGR 99, 8981, 1994) re-binned onto ~1 nm intervals.
! The Hinteregger spectrum, scaled using the EUVAC algorithm, is used
! from 18 to 50A.  Note that Richards et al. specified that the flux
! would not go below 0.8 of the reference spectrum; here that value is
! changed to 0.1.
!   Neither of these models extends shortward of 18A, so from 1-18 A
! an amalgam of sources are used to derive an estimated flux, e.g.,
! DeJager, in Astronomical Observations from Space Vehicles, Steinberg,
! ed., 1964; Smith & Gottlieb, SSR 16, 771, 1974; Manson, in The Solar
! Output and its Variation, White, ed., 1977; Kreplin et al, ibid;
! Horan & Kreplin, Solar Phys. 74, 265, 1981; Wagner, ASR 8, (7)67, 1988.
!   For FUV from 1050A-1750A, 50A interval bins from the Woods and
! Rottman [2002] reference spectrum and scale factors based on
! UARS SOLSTICE data are used.  The scaling method follows the
! Hinteregger or EUVAC algorithm, whichever is selected, so as to
! linearly scale the spectrum between the reference value and maximum
! value calculated with F10.7=F10.7A=200.
!   If ISCALE=2, the solar flux (0-1750A) is read from a file named
! ssflux_user.dat in the current working directory.  The file must
! contain three columns:  WAVES, WAVEL, SFLUX (Angstroms and cm-2 s-1)
! in order of increasing wavelength.  The number of lines in the file
! must match the value of LMAX in glow.h.  (Note: still need to implement a
! method for time-dependent user-supplied inputs, subroutine will only read
! the input file on the first call, unless ISCALE changes.)
!
! Modification history:
!   Stan Solomon, 12/1988  Basic Hinteregger EUV, approx. SME FUV
!   Chris Gaskill, 7/1989  Added early Tobiska model
!   Stan Solomon,  8/1989  Corrections to above
!   Stan Solomon,  1/1990  Tobiska SERF2; added W & R spectra
!   Stan Solomon,  6/1991  Tobiska EUV 91; Hntggr Ly-B, Fe XVI scaling
!   Stan Solomon,  2/1992  Updated Tobiska EUV91; corrected SME FUV
!   Scott Bailey, 12/1993  Initial one-nm bins version
!   Stan Solomon,  6/2004  Added EUVAC option, cleaned up artifacts
!   Stan Solomon,  9/2004  Added ability to specify input data file
!   Stan Solomon,  3/2005  Changed all to photon units
!   Stan Solomon,  1/2015  Updated for f90; only read file on first call
!   Stan Solomon,  4/2016  Changed EUVAC floor to 0.1
!   Stan Solomon,  4/2016  Removed obsolete parameters
!   Stan Solomon,  6/2016  Completed change to lower case and f90
!   Stan Solomon,  7/2016  Reads file if ISCALE has changed
!
! Calling parameters:
! ISCALE   =0 for Hinteregger contrast ratio method
!          =1 for EUVAC
!          =2 for user-supplied data
! F107     daily 10.7 cm flux (1.E-22 W m-2 Hz-1)
! F107A    81-day centered average 10.7 cm flux
! XUVFAC   factor for scaling flux 18-250A or 18-50A (optional)

! Returned parameters:
! WAVE1    longwave bound of spectral intervals (Angstroms)
! WAVE2    shortwave bound of intervals
! SFLUX    scaled solar flux returned by subroutine (photons cm-2 s-1)
!
! Other definitions:
! LMAX     dimension of flux and scaling arrays, currently = 123
! WAVEL    = WAVE1
! WAVES    = WAVE2
! RFLUX    low solar activity flux
! SCALE1   scaling factors for H LyB-keyed chromospheric emissions
! SCALE2   scaling factors for FeXVI-keyed coronal emissions
! B1       fit coefficients for H LyB
! B2       fit coefficients for FeXVI
! R1       enhancement ratio for H LyB
! R2       enhancement ratio for FeXVI
! P107     average of F107 and F107A
! A        scaling factor for EUVAC model


    subroutine ssflux (iscale,f107,f107a,xuvfac,wave1,wave2,sflux)

      use cglow,only: lmax,data_dir
      save

      dimension wave1(lmax), wave2(lmax), sflux(lmax), &
                wavel(lmax), waves(lmax), rflux(lmax), uflux(lmax), &
                scale1(lmax), scale2(lmax), a(lmax), b1(3), b2(3)
      data epsil/1.0E-6/
      data islast/-1/
      character(len=1024) :: filepath

! regression coefficients which reduce to solar min. spectrum:
      data b1/1.0, 0.0138, 0.005/, b2/1.0, 0.59425, 0.3811/

! 'best fit' regression coefficients, commented out, for reference:
!     DATA B1/1.31, 0.01106, 0.00492/, B2/-6.618, 0.66159, 0.38319/


! Hinteregger contrast ratio method:

      if (iscale .eq. 0) then
        if (islast .ne. iscale) then
          filepath = trim(data_dir)//'ssflux_hint.dat'
          open(unit=1,file=filepath,status='old',readonly)
          read(1,*)
          do l=lmax,1,-1
            read(1,*) waves(l),wavel(l),rflux(l),scale1(l),scale2(l)
          enddo
          close(unit=1)
        endif
!
        r1 =  b1(1) + b1(2)*(f107a-71.5) + b1(3)*(f107-f107a+3.9)
        r2 =  b2(1) + b2(2)*(f107a-71.5) + b2(3)*(f107-f107a+3.9)
!
        do l=1,lmax
          sflux(l) = rflux(l) + (r1-1.)*scale1(l) + (r2-1.)*scale2(l)
          if (sflux(l) .lt. 0.0) sflux(l) = 0.0
          if (xuvfac .gt. epsil .and. wavel(l).lt.251.0 .and. waves(l).gt.17.0) &
            sflux(l)=sflux(l)*xuvfac
        enddo
      endif

! EUVAC Method:

      if (iscale .eq. 1) then
        if (islast .ne. iscale) then
          filepath = trim(data_dir)//'ssflux_euvac.dat'
          open(unit=1,file=filepath,status='old',readonly)
          read(1,*)
          do l=lmax,1,-1
            read(1,*) waves(l),wavel(l),rflux(l),a(l)
          enddo
          close(unit=1)
        endif

      p107 = (f107+f107a)/2.

        do l=1,lmax
          sflux(l) = rflux(l) * (1. + a(l)*(p107-80.))
          if (sflux(l) .lt. 0.1*rflux(l)) sflux(l) = 0.1*rflux(l)
          if (xuvfac .gt. epsil .and. wavel(l).lt.51.0 .and. waves(l).gt.17.0) &
            sflux(l)=sflux(l)*xuvfac
        enddo
      endif

! User-supplied data:

      if (iscale .eq. 2) then
        if (islast .ne. iscale) then
          filepath = trim(data_dir)//'ssflux_user.dat'
          open(unit=1,file=filepath,status='old',readonly)
          read(1,*)
          do l=lmax,1,-1
            read(1,*) waves(l),wavel(l),uflux(l)
          enddo
          close(unit=1)
        endif
        do l=1,lmax
          sflux(l)=uflux(l)
        enddo
      endif

! Fill wavelength arrays, substitute in H Lyman-alpha if provided:

      do l=1,lmax
        wave1(l) = wavel(l)
        wave2(L) = waves(l)
      enddo

      islast=iscale

      return

    end
