! Subroutine CONDUCT

! Calculates Pedersen and Hall conductivity from ion/neutral densities and temperatures.
 
! Ryan McGranaghan, 2014
! Modifications by Stan Solomon, 6/2016:
!   corrected Tr=(Ti+Tn)/2 instead of (Ti+Te)/2
!   changed to Ne = O+ + NO+ + O2+ for charge neutrality and correct Hall conductance
!   added Tn to call, removed Ne, and re-ordered calling parameters
!   changed densities to cm^-3 (and multiplied results by 1.e6)
!   cleaned up comments
!   refactored for f90

! Input parameters:
!  lat, lon, alt      : latitude [deg], longitude [deg], altitude [km]
!  nO, nO2, nN2       : neutral constituent densities [cm^-3]
!  n0p, nO2p, nNOp    : ion consituent densities [cm^-3]
!  Tn, Ti, Te         : neutral, electron and ion temperatures [K]

! Output parameters:
!  PedCond            : Pedersen conductivity [S/m]
!  HallCond           : Hall conductivity [S/m] 

! Variable definitions:
!   Ion-neutral momentum transfer collision frequencies (from TIE-GCM lamdas.F routine):
!     nu_o2po2 ! O2+ ~ O2 collision freq (resonant, temperature dependent)
!     nu_opo2  ! O+  ~ O2 collision freq (non-resonant)
!     nu_nopo2 ! NO+ ~ O2 collision freq (non-resonant)
!     nu_o2po  ! O2+ ~ O  collision freq (non-resonant)
!     nu_opo   ! O+  ~ O  collision freq (resonant, temperature dependent)
!     nu_nopo  ! NO+ ~ O  collision freq (non-resonant)
!     nu_o2pn2 ! O2+ ~ N2 collision freq (non-resonant)
!     nu_opn2  ! O+  ~ N2 collision freq (non-resonant)
!     nu_nopn2 ! NO+ ~ N2 collision freq (non-resonant)
!     nu_o2p   ! [[o2p~o2]n(o2)+[o2p~o]n(o)+[o2p~n2]n(n2)]
!     nu_op    ! [[op ~o2]n(o2)+[op ~o]n(o)+[op ~n2]n(n2)]
!     nu_nop   ! [[nop~o2]n(o2)+[nop~o]n(o)+[nop~n2]n(n2)]
!     nu_ne    ! electron~neutral

    SUBROUTINE CONDUCT(lat, long, alt, nO, nO2, nN2, nOp, nO2p, nNOp, &
                         Tn, Ti, Te, PedCond, HallCond)

      implicit none

      real, intent(in) :: lat,long,alt,nOp,nO2p,nNOp
      real, intent(in) :: Tn, Te, Ti, nO, nO2, nN2
      real, intent(out) :: PedCond, HallCond

      real :: XFmag, YFmag, ZFmag, Bmagn, DIPmag, DECmag, SDIPmag
      real :: om_Op, om_O2p, om_NOp, om_e
      real :: nu_O2pO2, nu_OpO2, nu_NOpO2, nu_OpO, nu_NOpO, nu_O2pO
      real :: nu_O2pN2, nu_NOpN2, nu_OpN2, nu_O2p, nu_Op, nu_NOp, nu_en
      real :: r_O2p, r_Op, r_NOp, r_e
      real :: nE

      real, parameter :: Na = 6.0221413e23   ! Avagadros number
      real, parameter :: me = 9.10938291e-31 ! electron mass in kg
      real, parameter :: qe = 1.60217657e-19 ! electron charge in C
      real, parameter :: Mbar_O2 = 0.031999  ! molecular mass of O2 in kg/mol
      real, parameter :: Mbar_O  = 0.0159995 ! molecular mass of O in kg/mol 
      real, parameter :: Mbar_NO = 0.030     ! molecular mass of NO in kg/mol
      real, parameter :: Burnfac = 1.5       ! Burnside Factor

! Calculate magnetic field strength and convert to from Gauss to Tesla:

      call FIELDM (lat,long,alt,XFmag,YFmag,ZFmag,Bmagn,DIPmag,DECmag, SDIPmag)
      Bmagn = Bmagn*1.e-4
      
! Calculate gyro frequencies Using Tesla, Coulomb, and m^-3:
          
      om_Op = (qe*Bmagn*Na)/(Mbar_O)
      om_O2p = (qe*Bmagn*Na)/(Mbar_O2)
      om_NOp = (qe*Bmagn*Na)/(Mbar_NO)
      om_e = ( (qe*Bmagn)/me )
      
! Calculate collision frequencies using method adapted from the TIE-GCM
! which includes the Burnside factor correction for O+-O collision frequency.
! NOTE: densities are in cm^-3 and coefficients are in cm^3 s^-1.

      nu_O2pO2 = nO2*(2.59e-11)*sqrt((Ti+Tn)/2.)*((1.-0.073*log10((Ti+Tn)/2.))**2)      
      nu_OpO2 = nO2*(6.64e-10)
      nu_NOpO2 = nO2*(4.27e-10)
      nu_OpO = nO*((3.67e-11)*sqrt((Ti+Tn)/2.)*((1.-0.064*log10((Ti+Tn)/2.))**2)*BurnFac)
      nu_NOpO = nO*(2.44e-10)
      nu_O2pO = nO*(2.31e-10)
      nu_O2pN2 = nN2*(4.13e-10)
      nu_NOpN2 = nN2*(4.34e-10)
      nu_OpN2 = nN2*(6.82e-10)

! Sum the specific collision frequencies and calculate electron
! collision frequency with neutrals (ignoring the contribution from
! collisions with ions):

      nu_O2p = nu_O2pO2 + nu_O2pO + nu_O2pN2
      nu_Op = nu_OpO2 + nu_OpO + nu_OpN2
      nu_NOp = nu_NOpO2 + nu_NOpO + nu_NOpN2
      nu_en = (2.33e-11) * nN2*Te*(1.-(1.21e-4)*Te) &
             + (1.82e-10)*nO2*sqrt(Te)*(1.+(3.6e-2)*sqrt(Te)) &
             + (8.9e-11)*nO*sqrt(Te)*(1.+(5.7e-4)*Te) 

! Comment from TIE-GCM source code in file lamdas.F:
!  6/2/06 btf: Multiply rnu_ne by 4, as per Art Richmond:
!  The effective electron-neutral collision frequency is increased in 
!  an ad-hoc manner by a factor of 4 in order for the model to produce
!  electric fields and currents below 105 km that agree better with
!  observations, as recommended by Gagnepain et al. (J. Atmos. Terr. 
!  Phys., 39, 1119-1124, 1977).

      nu_en = nu_en * 4.
         
! Calculate the ratios of collision frequency to gyro frequency:

      r_O2p = nu_O2p/om_O2p
      r_Op = nu_Op/om_Op
      r_NOp = nu_NOp/om_NOp
      r_e = nu_en/om_e

! Define electron density for purpose of conductivity calculations to be
! the sum of the major ions:

      nE = nOp+nO2p+nNOp

! Calculate Pederson and Hall conductivity using the TIE-GCM formulation:
! NOTE: multiplied by 1.e6 to convert densities to m^-3, so that
! conductivities are in Siemens/meter, a.k.a. mhos.


      PedCond = (1.e6*qe/Bmagn) *                                 &
                                  ( nOp*(r_Op/(1+r_Op**2)) +      &
                                    nO2p*(r_O2p/(1+r_O2p**2)) +   &
                                    nNOp*(r_NOp/(1+r_NOp**2)) +   &
                                    nE*(r_e/(1+r_e**2)) )

      HallCond = (1.e6*qe/Bmagn) *                                &
                                   ( nE/(1+r_e**2) -              &
                                     nOp/(1+r_Op**2) -            &
                                     nO2p/(1+r_O2p**2) -          &
                                     nNOp/(1+r_NOp**2) )

      if (HallCond < 0.) HallCond=0.

      return
    end subroutine conduct
