! Subroutine QBACK

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Stan Solomon, 11/1988, 11/1992, 3/2005
! New version uses updated TIE-GCM and TIME-GCM qinite.F formulation, SCS, 1/2013
! Refactored to f90, SCS, 6/2016

! Estimates background ("nighttime") ionization rates.
! Four components are used:
!   Geocoronal Lyman-beta 102.6 nm (ionizes O2 only)
!   Geocoronal He I 58.4 nm
!   Geocoronal He II 30.4 nm
!   Geocoronal Lyman-alpha 121.6 nm (ionizes NO only)


! Inputs:
!   zmaj    major species O, O2, N2 at each altitude
!   zno     nitric oxide at each altitude
!   zvcd    vertical column density for each major species above each altitude
!   photoi  array of photoionization rates for each state, species, altitude
!   f107    solar 10.7 cm radio flux activity index
!   jm      number of altitude levels
!   nmaj    number of major species (3)
!   nst     number of states
! Output:
!   Array photoi is incremented by the background rates

!Other definitions:
!   al      photon flux at 102.6 nm, 58.4 nm, 30.4 nm
!   flyan   photon flux at 121.6 nm
!   sa      absorption cross sections for O, O2, N2 at each wavelength
!   si      ionization cross sections for O, O2, N2 at each wavelength
!   salyao2 absorption cross section for O2 at 121.6 nm
!   silyano ionization cross section for NO at 121.6 nm
!   bn2p    branching ratio for N2+ from ionization of N2
!   bn1p    branching ratio for N+ from ionization of N2
!   tau     optical depth
!   qbo1    production rate of O+
!   qbo2    production rate of O2+
!   qbn2    production rate of N2+
!   qbn1    production rate of N+
!   qbno    production rate of NO+

! All units cgs.


  subroutine qback (zmaj,zno,zvcd,photoi,phono,f107,jm,nmaj,nst)

    implicit none

    integer,intent(in) :: jm, nmaj, nst
    real,intent(in) :: zmaj(nmaj,jm), zno(jm), zvcd(nmaj,jm), f107
    real,intent(inout) :: photoi(nst,nmaj,jm), phono(nst,jm)

    integer :: j, l
    real :: al(3), sa(3,3), si(3,3), salyao2, silyano, bn2p, bn1p
    real :: flyan, tau, qbo1, qbo2, qbn2, qbn1, qbno

    data al /1.5e7, 1.5e6, 1.5e6/
    data sa /      0.,  1.6e-18,       0., &
             10.2e-18, 22.0e-18, 23.1e-18, &
              8.4e-18, 16.0e-18, 11.6e-18/
    data si /      0.,  1.0e-18,       0., &
             10.2e-18, 22.0e-18, 23.1e-18, &
              8.4e-18, 16.0e-18, 11.6e-18/ 
    data salyao2/8.0e-21/
    data silyano/2.0e-18/
    data bn2p/0.86/
    data bn1p/0.14/

! Calculate Lyman-alpha 121.6 nm geocoronal flux as a function of F10.7:

    flyan = 5.E9*(1.+0.002*(f107-65.))

! Loop over altitudes:

    do j=1,jm
      qbo1 = 0.
      qbo2 = 0.
      qbn2 = 0.
      qbn1 = 0.

! Calculate optical depth and ionization rates for major species:

      do l=1,3
        tau=(sa(1,l)*zvcd(1,j)+sa(2,l)*zvcd(2,j)+sa(3,l)*zvcd(3,j))
        qbo1 = qbo1+al(l)*si(1,l)*zmaj(1,j)*exp(-tau)
        qbo2 = qbo2+al(l)*si(2,l)*zmaj(2,j)*exp(-tau)
        qbn2 = qbn2+bn2p*(al(l)*si(3,l)*zmaj(3,j)*exp(-tau))
        qbn1 = qbn1+bn1p*(al(l)*si(3,l)*zmaj(3,j)*exp(-tau))
      enddo

! Calculate optical depth of Ly-alpha, and ionization rate of NO:

      tau = salyao2*zvcd(2,j)
      qbno = flyan*silyano*zno(j)*exp(-tau)

! Increment ionization rate array with background ionization:

      photoi(1,1,j) = photoi(1,1,j) + qbo1
      photoi(1,2,j) = photoi(1,2,j) + qbo2
      photoi(1,3,j) = photoi(1,3,j) + qbn2
      photoi(6,3,j) = photoi(6,3,j) + qbn1
      phono(1,j) = phono(1,j) + qbno

    enddo

    return
  end subroutine qback
