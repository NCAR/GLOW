! Subroutine GLOW

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Version 0.98, 1/2017
! Version 0.981, 6/2017

! Stan Solomon, 1988, 1989, 1990, 1991, 1992, 1994, 2000, 2002, 2005, 2015, 2016
!
! Subroutine GLOW is the master routine of the /glow package.
! It receives input parameters from the calling program using use-associated variables
! defined in module CGLOW, calls the other subroutines, and returns results.
! CGLOW also defines the array dimensions for the altitude grid, the electron
! energy grid, solar spectrum wavelengh bins, etc.

! Call subroutine CGLOW_INIT first to set array dimensions and allocate variables
! Call subroutine EGRID before call to GLOW to set up electron energy grid
! Call subroutine MAXT before call to GLOW to specify auroral electron flux (if any)

! Subroutines called by GLOW are:
!   FIELDM  calculates magnetic dip angle
!   SOLZEN  calculates solar zenith angle
!   SSFLUX  scales solar flux for activity level
!   RCOLUM  calculates slant column density of major species
!   EPHOTO  calculates photoionization and photoelectron production
!   QBACK   estimates background ('night time') ionization
!   ETRANS  computes electron transport, ionization, excitation
!           calls EXSECT for cross-sections, first call only
!   GCHEM   finds electron/ion/metastable densities, airglow emissions
!   BANDS   calculates vibrational distributions for selected band systems (currently only LBH)

! Supplied to subroutine using use-associated data defined in module CGLOW:
! IDATE   Date, in form yyddd
! UT      Universal Time; seconds
! GLAT    Geographic latitude; degrees
! GLONG   Geographic longitude; degrees
! ISCALE  Solar flux scaling switch, see subroutine SSFLUX
! JLOCAL  =0 for electron transport calculation, =1 for local calc only
! KCHEM   Ion/electron chemistry switch, see subroutine GCHEM
! F107    Solar 10.7 cm flux for day being modeled, 1.E-22 W m-2 Hz-1
! F107A   Solar 10.7 cm flux 81-day centered average
! XUVFAC  Factor by which to multiply to solar flux 16-250 A or 16-50 A.
! ZZ      altitude array; cm
! ZO      O number density at each altitude; cm-3
! ZN2     N2  "      "      "   "     "       "
! ZO2     O2         "
! ZNO     NO         "
! ZNS     N(4S)      "
! ZND     N(2D)      "
! ZRHO    mass density at each altitude; gm cm-3 (not currently in use)
! ZE      electron density at each alt; cm-3
! ZTN     neutral temperature at each alt; K
! ZTI     ion temperature at each alt; K
! ZTE     electron temp at each alt; K
! PHITOP  energetic electron flux into top of atmosphere; cm-2 s-1 eV-1

! Calculated by subroutine and returned using use-associated data defined in module CGLOW:
! SZA     solar zenith angle; radians
! DIP     magnetic field dip angle; radians
! EFRAC   energy conservation check from ETRANS, (out-in)/in
! IERR    error code returned from ETRANS:
!           0=normal, 1=local problem, 2=transport problem
! ZMAJ    major species density array, O, O2, N2; cm-3
! ZCOL    major species slant column density array, O, O2, N2; cm-2
! WAVE1   longwave edge of solar flux wavelength range; A
! WAVE2   shortwave edge of solar flux wavelength range; A
! SFLUX   scaled solar flux in each wavelength range; photons cm-2 s-1
! ENER    electron energy grid; eV
! DEL     width of each bin in electron energy grid; eV
! PESPEC  photoelectron production rate at energy, altitude; cm-3 s-1
! PIA     proton aurora ionization rate (not currently in use); cm-3 s-1.
! SESPEC  proton aurora secondary electron production rate (not currently in use); cm-3 s-1
! PHOTOI  photoionization rates for state, species, altitude; cm-3 s-1
!           O+ states: 4S, 2Do, 2Po, 4Pe, 2Pe
!           O2+ states: X, a+A, b, dissoc.
!           N2+ states: X, A, B, C, F, dissoc.
! PHOTOD  photodissoc. & exc. rates for state, species, alt.; cm-3 s-1
!           (1,2,J) = O2 -> O(3P) + O(1D))
!           (2,2,J) = O2 -> O(3P) + O(1S)
!           (1,3,J) = N2 -> N + N
! PHONO   photoionization/dissociation/excitation rates for NO, cm-3 s-1
!         (1,J) = NO+ from H Ly-a ionization
! SION    electron impact ionization rates calculated by ETRANS; cm-3 s-1
! UFLX    upward hemispherical electron flux; cm-2 s-1 eV-1
! DFLX    downward hemispherical electron flux; cm-2 s-1 eV-1
! AGLW    Electron impact exc. rates; state, species, alt.; cm-3 s-1
!           O states: 1D, 1S, 5S, 3S, 3p5P, 3p3P, 3d3D, 3s'3D
!           O2 states: a, b, (A+A'+c), B(SRC), 9.9eV, Ryds., vib.
!           N2 states: (A+B+W), B', C, (a+a'+w), 1Pu, b', Ryds., vib.
! EHEAT   ambient electron heating rate, eV cm-3 s-1
! TEZ     total energetic electron energy deposition, eV cm-3 s-1
! TPI     total photoionization rate at each altitude, cm-3 s-1
! TEI     total electron impact ionization rate at each altitude, cm-3 s-1
! TIR     total ionization rate at each altitude (TPI+TEI), cm-3 s-1
! ECALC   electron density, calculated below 200 km, cm-3
! ZXDEN   array of excited and and/or ionized state densities at each altitude:
!           O+(2P), O+(2D), O+(4S), N+, N2+, O2+, NO+, N2(A), N(2P),
!           N(2D), O(1S), O(1D); cm-3
! ZETA    array of volume emission rates at each altitude:
!           3371A, 4278A, 5200A, 5577A, 6300A, 7320A, 10400A, 3466A,
!           7774A, 8446A, 3726A, LBH, 1356, 1493, 1304; cm-3 s-1
! ZCETA   array of contributions to each v.e.r at each alt; cm-3 s-1
! VCB     array of vertical column brightnesses (as above); Rayleighs

! Array dimensions:
! JMAX    number of altitude levels
! NBINS   number of energetic electron energy bins
! LMAX    number of wavelength intervals for solar flux
! NMAJ    number of major species
! NST     number of states produced by photoionization/dissociation
! NEI     number of states produced by electron impact
! NEX     number of ionized/excited species
! NW      number of airglow emission wavelengths
! NC      number of component production terms for each emission


    subroutine glow

      use cglow,only: jmax,lmax,nw,nst,nmaj,nbins,iscale,nei,ierr, &
                      glat,glong,idate,ut,f107,f107a, &
                      ener,del,dip,sza,xuvfac, &
                      wave1,wave2,sflux,zmaj,zo,zo2,zn2,zz,ztn,zcol, &
                      photoi,photod, phono,pespec,pia,sespec,phitop, &
                      uflx,dflx,sion,aglw,eheat,tez, efrac,zno

      implicit none

      real :: zvcd(nmaj,jmax),xf,yf,zf,ff,dec,sdip,teflux

      real,parameter :: pi=3.1415926536
      integer,save :: ifirst=1
      integer :: j,i,ist,n,iei

! First call only: set up energy grid:

      if (ifirst == 1) then
        ifirst = 0
        call egrid (ener, del, nbins)
      endif

! Find magnetic dip angle and solar zenith angle (radians):

      call fieldm (glat, glong, 300., xf, yf, zf, ff, dip, dec, sdip)
      dip = abs(dip) * pi/180.
      if (dip < 0.01) dip=0.01

      call solzen (idate, ut, glat, glong, sza)
      sza = sza * pi/180.

! Scale solar flux:

      call ssflux (iscale, f107, f107a, xuvfac, wave1, wave2, sflux)

! Pack major species density array:

      do j=1,jmax
        zmaj(1,j) = zo(j)
        zmaj(2,j) = zo2(j)
        zmaj(3,j) = zn2(j)
      enddo

! Calculate slant path column densities of major species in the
! direction of the sun:

      call rcolum (sza, zz, zmaj, ztn, zcol, zvcd, jmax, nmaj)

! Call subroutine EPHOTO to calculate the photoelectron production
! spectrum and photoionization rates as a function of altitude,
! unless all altitudes are dark, in which case zero arrays:

      if (sza < 1.85) then
        call ephoto
      else
        photoi(:,:,:) = 0.0
        photod(:,:,:) = 0.0
        phono(:,:) = 0.0
        pespec(:,:) = 0.0
      endif

! Zero proton aurora ionization rates and secondary electron production (not currently in use):

      pia(:,:) = 0.0
      sespec(:,:) = 0.0

! Add background ionization to photoionization:

      call qback (zmaj, zno, zvcd, photoi, phono, f107, jmax, nmaj, nst)

! Call subroutine ETRANS to calculate photoelectron and auroral
! electron transport and electron impact excitation rates, unless
! there are no energetic electrons, in which case zero arrays:

      teflux = 0.
      do n=1,nbins
        teflux = teflux + phitop(n)
      enddo

      if (teflux > 0.001 .or. sza < 1.85) then
        call etrans
      else
        uflx(:,:) = 0.0
        dflx(:,:) = 0.0
        sion(:,:) = 0.0
        aglw(:,:,:) = 0.0
        eheat(:) = 0.0
        tez(:) = 0.0
        efrac = 0.0
        ierr = 0
      endif

! Call subroutine GCHEM to calculate the densities of excited and
! ionized consituents, airglow emission rates, and vertical column
! brightnesses:

      call gchem

! Call subroutine BANDS to calculate band-specific airglow emission
! rates (currently only LBH upper state distribution):

      call bands

      return

    end subroutine glow
