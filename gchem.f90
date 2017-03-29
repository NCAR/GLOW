! Subroutine GCHEM

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Stan Solomon, 1988, 1989, 1992, 1999, 2005, 2016

! 3/2016: Replaced quartic solution to electron density equation
! with iterative method.  Also added constraint O+ > 0 for KCHEM=3.
! 1/2017: Removed 7774 -> 1356 fudge factor (and reduced 7774 cross section in exsect.f)

! Electron density must be supplied above 200 km in array ZE, a priori
! values are expected but not necessarily required at and below 200 km,
! depending on the value of KCHEM.  An initial guess of the N(2D)
! density in array ZND is also expected but not required.
! Other neutral species (O, N2, O2, NO, N(4S)) must be supplied,
! see subroutine GLOW.

! Chemical calculations are controlled by switch KCHEM:
! 0 = no calculations at all are performed.
! 1 = electron density, O+(4S), N+, N2+, O2+, NO+ supplied at all
!     altitudes; O+(2P), O+(2D), excited neutrals, and emission rates
!     are calculated.
! 2 = electron density, O+(4S), O2+, NO+ supplied at all altitudes;
!     O+(2P), O+(2D), N+, N2+, excited neutrals, emissions calculated.
! 3 = electron density supplied at all altitudes; everything else
!     calculated.  Note that this may violate charge neutrality and/or
!     lead to other unrealistic results below 200 km, if the electron
!     density supplied is significantly different from what the model
!     thinks it should be.  If it is desired to use a specified ionosphere,
!     KCHEM=2 is probably a better option.
! 4 = electron density supplied above 200 km; electron density below
!     200 km is calculated, everything else calculated at all altitudes.
!     Electron density for the next two levels above J200 is log interpolated
!     between E(J200) and E(J200+3).

! For definitions of use-associated variables see subroutine GLOW and module CGLOW.

! Other definitions:
! A        Einstein coefficients; s-1
! B        Branching ratios
! BZ       Altitude-dependent branching ratios
! G        Resonant scattering g-factors at each altitude; s-1
! KZ       Temperature dependent rate coeffs at each altitude; cm3s-1
! OEI      O electron impact ionization rates; cm-3s-1
! O2EI     O2   "       "        "       "      "
! RN2EI    N2   "       "        "       "      "
! TEI      Total "      "        "       "      "
! OPI      O photoionization rate, cm-3s-1
! O2PI     O2       "          "      "
! RN2PI    N2       "          "      "
! TPI      Total    "          "      "
! TIR      Total ionization rate; cm-3 s-1
! RN2ED    N2 electron impact dissociation rate; cm-3s-1
! SRCED    O2     "      "         "         " (SR continuum); cm-3s-1
! P        Volume production rate for each species, altitude; cm-3s-1
! L        Loss rate for each species, altitude; s-1
! OMINUS   O- density for mutual neutralization contribution to O*
! T1       Effective temperature divided by 300 for O+ + N2; K
! T2          "          "        "                 O+ + O2; K
! T3          "          "        "                 N2+ + O; K
! T4          "          "        "                 N2+ + O2; K
! T5          "          "        "                 O+ + NO; K
! QQ, RR, SS, TT, UU, VV, WW, XX:  Combined terms for calculation of
!                              O+(4S) given e

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
! NR      number of rate coefficients, branching ratios, A and G factors

! References for rate coefficients, transition coefficients, branching ratios, and g-factors:
!      k1  O+(4S) + N2             St. Maurice & Torr, 1978 (from Albritton et al., 1977)
!      k2  O+(4S) + O2             Combination of (Chen et al, 1978 and St. Maurice & Torr, 1978)
!      k3  N2+ + O -> NO+ + O      New fit to McFarland et al, 1974
!      k4  N2+ + O2                McFarland et al, 1973
!      k5  N(2D) + O2              Lin & Kaufman, 1971, cf. Piper et al, 1987
!      k6  N(2D) + O               Fell et al., 1990
!      k7  N(2D) + e               Frederick & Rusch, 1977; Queffelec et al, 1985
!      k8  O(1D) + N2              Streit et al, 1976
!      k9  O(1D) + O2              Streit et al, 1976
!      k10 O(1D) + e               Link, 1982
!      k11 O(1S) + O               Slanger & Black, 1981
!      k12 O+(2D) + N2             Johnsen & Biondi, 1980
!      k13 O+(2D) + O2             Johnsen & Biondi, 1980
!      k14 O+(2D) + e              Henry et al., 1969
!      k15 O+(2D) + O              Torr & Torr, 1980
!      k16 O+(2P) + N2             Rusch et al, 1977
!      k17 O+(2P) + O2             Link, 1982
!      k18 O+(2P) + e              Henry et al., 1969
!      k19 O+(2P) + O              Rusch et al, 1977
!      k20 O2(c) + O               Solheim & Llewellyn, 1979
!      k21 O2(c) + N2              Solheim & Llewellyn, 1979
!      k22 NO+ + e                 Walls & Dunn, 1974; Torr et al, 1977; Alge et al, 1983;
!                                  Dulaney et al,1987; Davidson & Hobson, 1987
!      k23 N2+ + e                 Mehr & Biondi, 1969
!      k24 O2+ + e                 Mehr & Biondi; Walls & Dunn; Torr et al; Alge et al, 1983
!      k25 N+ + O2                 Langford et al, 1985
!      k26 N2(A) + O               Piper et al, 1981b (av. v=1,2)
!      k27 O(1D) + O               Abreu et al, 1986; Yee, pc, 1991
!      k28 O + et                  Link, 1982
!      k29 N2(A) + O2              Piper et al, 1981a
!      k30 O2+ + NO                Lindeger & Ferguson, 1974; G&R, 1979
!      k31 N(2D) + NO              Black et al, 1969; fr. Roble, 1986
!      k32 N+ + O                  Torr, 1985 (Levine, Photochemistry)
!      k33 N(2P) + O               Zipf et al, 1980; Young & Dunn, 1975
!      k34 N(2P) + O2              Zipf et al, 1980; Rawlins, 1988
!                                  (cf Ianuzzi & Kaufman, 1980, 3.5E-12)
!      k35 N(2P) + NO              Rees & Jones, 1973, Zipf et al, 1980
!      k36 O(1S) + O2              Slanger et al, 1972; fr. Bates, 1978
!      k37 O2+ + N                 Fehsenfeld (1977)
!      k38 O+ + N(2D)              Bates, 1989 (PSS 37, 363)
!      k39 N2+ + O -> N2 + O+      Torr, 1985; Torr et al, 1988; Knutsen et al, 1988
!      k40 O+ + NO                 St. Maurice & Torr, 1978
!      k41 O+ + e -> 7774, 1356    Melendez, 1999; Qin, 2015 (cf. Tinsley 1973; Julienne, 1974)
!      k42 O + e -> O-             Melendez, 1999; Qin, 2015 (cf. Tinsley 1973)
!      k43 O- + O+ -> O* + O       Melendez, 1999; Qin, 2015 (cf. Tinsley 1973)
!      k44 O- + O+ -> O2 + e       Melendez, 1999; Qin, 2015 (cf. Tinsley 1973)
!      k45 O+ + e -> 8446, 1304    Estimate from Tinsley, 1973; Julienne, 1974)
!      A1  5200       N(4S-2D)     Wiese et al, 1966
!      A2  6300       O(3P-1D)     Baluja and Zeippen, 1988
!      A3  6364       O(3P-1D)     Baluja and Zeippen, 1988
!      A4  2972       O(3P-1S)     Kernahan & Pang, 1975
!      A5  5577       O(1D-1S)     Kernahan & Pang, 1975
!      A6  3726       O+(4S-2D)    Kernahan & Pang, 1975
!      A7  2470       O+(4S-2P)    Weise et al, 1966
!      A8  7319-30    O+(2D-2P)    Weise et al, 1966
!      A9  (Hertz II) O2(X-c)      Solheim & Llewellyn, 1978
!      A10 (Veg-Kap)  N2(X-A)      Shemansky, 1969
!      A11 3466       N(4S-2P)     Chamberlain, 1961
!      A12 10400      N(2D-2P)     Chamberlain, 1961
!      B1  O(1S) from O2+ + e      Yee et al, 1988
!      B2  O(1D) from O2+ + e      Abreu et al, 1986
!      B3  N(2D) from NO+ + e      Kley et al, 1976
!      B4  N(2D) from N2+ + e      Queffelec et al, 1985
!      B5  N(2D) from N2+ + O      Frederick & Rusch, 1977
!      B6  O(1D) from N(2D) + O2   Link, 1983; Langford et al, 1985
!      B7  O(1D) from O+(2D) + O   ?
!      B8  O+(2D) from O+(2P) + e  Link, 1982
!      B9  O+(2P) from O + e*      Gerard & Rusch, 1979; Jones, 1975
!      B10 O+(2D) from O + e*      Gerard & Rusch, 1979; Jones, 1975
!      B11 O+(4S) from O + e*      Gerard & Rusch, 1979; Jones, 1975
!      B12 O+(2P) from O2 + e*     Link, 1982; guess of .3/3
!      B13 O+(2D) from O2 + e*     Link, 1982; guess of .3/3
!      B14 O+(4S) from O2 + e*     Link, 1982; guess of .3/3
!      B15 N+ from N2 + e*         Richards & Torr, 1985
!      B16 N(2D) from above        Zipf et al, 1980
!      B17 O(1D) from N+ + O2      Langford et al, 1985
!      B18 O(1S) from N2(A) + O    Sharp & Torr, 1979
!      B19 O(1S) from O2(*) + O    ? (= 0 at present)
!      B20 O(1D) from N(2D) + O    ?
!      B21 NO+ from N+ + O2        Langford et al, 1985
!      B22 O2+ from N+ + O2        Langford et al, 1985
!      B23 N(2P) from N2+ + e      Queffelec et al, 1985
!      B24 N2 + protons -> N + N   ?
!      B25 N(2D) from N2 + e* dis  Zipf et al, 1980
!      B26 N(2P) from N2 + e* dis  Zipf et al, 1980
!      B27 N(2D) from N2 + hv      Richards et al, 1981 (add to B28)
!      B28 N(2P) from N2 + hv      ? (cf Zipf & McGlaughlin, 1978)
!      B29 N(2D) from N(2P) + O    ?
!      B30 O+(2P) from O2 + hv     ?
!      B31 O+(2D)  "               ?
!      B32 O+(4S)  "               ?
!      B33 O(1S) from O2+ + N      Frederick et al, 1976; Kopp ea, 1977
!      B34 O(1S) from N(2D) + NO   Frederick et al, 1976; Kopp ea, 1977
!      B35 O2 + protons -> (O1D)   ?
!      B36 N2+(B) from N2 + e*     Borst & Zipf, 1970; Shemansky & Broadfoot, 1971
!      B37 (0,0) (3914) fr. N2+(B) Shemansky & Broadfoot, 1971
!      B38 (0,1) (4278) fr. N2+(B) Shemansky & Broadfoot, 1971
!      B39 (0,0) (3371) fr. N2(C)  Conway, 1983; Benesch et al, 1966
!      B40 (0,9) (3352) fr. N2(A)  Cartwright, 1978; Shemansky, 1969
!      B41 O+(2Po) fr. O+(2Pe)     Kirby et al, 1979
!      B42 O+(2Do) fr. O+(2Pe)     Kirby et al, 1979
!      B43 N2(C) bound fraction    ?
!      B44 7990 fr. O(3s'3D)       appx. fr. Hecht, p.c.
!      B45                         not currently in use 
!      B46 N 1493 fr. N2+hv DI     guess 
!      B47 N 1493 fr. N2+e* DI     guess, cf. Mumma and Zipf (1973), Meier (1991)
!      B48 N2(a) from (a,a',w)     estimate from comparison with Ajello & Shemansky, GUVI data, etc.
!      B49 7774, 1356 fr. O-+O+    Melendez, 1999; Qin, 2015 (cf. Tinsley 1973; Julienne, 1974)
!      G1  N2+B(0,0) (3914)        Broadfoot, 1967
!      G2  N2+B(0,1) (4278)        Broadfoot, 1967


  SUBROUTINE GCHEM
!
    use cglow,only: jmax, nmaj, nex, nw, nc, kchem, sza, &
                    zz, zo, zn2, zo2, zno, zns, znd, ze, ztn, zti, zte, &
                    photoi, photod, phono, pia, sion, aglw, &
                    tei, tpi, tir, e=>ecalc, den=>zxden, zeta, zceta, vcb
!
    implicit none
    integer,parameter :: nr=50
    real,parameter :: re=6.37E8
!
    real ::   A(NR), B(NR), BZ(NR,JMAX), G(NR,JMAX), KZ(NR,JMAX), &
              OEI(JMAX), O2EI(JMAX), RN2EI(JMAX), &
              OPI(JMAX), O2PI(JMAX), RN2PI(JMAX), &
              RN2ED(JMAX), SRCED(JMAX), P(NEX,JMAX), L(NEX,JMAX), OMINUS(JMAX), &
              T1(JMAX), T2(JMAX), T3(JMAX), T4(JMAX), T5(JMAX), &
              QQ(JMAX), RR(JMAX), SS(JMAX), TT(JMAX), UU(JMAX), &
              VV(JMAX), WW(JMAX), XX(JMAX)
    real ::   gh,dz,tatomi,alphaef,toti
    integer :: i,iw,ic,ix,n,j200,iter
!
    DATA A/ 1.07E-5, 0.00585, 0.00185, 0.0450, 1.0600, 9.7E-5, 0.0479, 0.1712, 0.0010, 0.7700, &
            0.00540, 0.07900,     0.0,    0.0,    0.0,     0.0,   0.0,    0.0,    0.0,    0.0, &
            30*0.0 /
    DATA B/ 0.07, 1.20, 0.76, 1.85, 1.00, 0.10, 0.50, 0.81, 0.20, 0.32, &
            0.48, 0.10, 0.10, 0.10, 0.16, 0.50, 0.30, 0.19, 0.00, 0.10, &
            0.43, 0.51, 0.10, 0.60, 0.54, 0.44, 0.80, 0.20, 1.00, 0.33, &
            0.33, 0.34, 0.21, 0.20, 0.10, 0.11, 0.65, 0.20, 0.24, 0.02, &
            0.18, 0.72, 0.75, 0.10, 0.00, 0.05, 0.02, 0.70, 0.54, 0.00  /
!
!
    IF (KCHEM .EQ. 0) RETURN
!
!
! Zero airglow and density arrays:
!
    zeta(:,:) = 0.
    zceta(:,:,:) = 0.
    vcb(:) = 0.
    if (kchem .ge. 3) den(:,:) = 0.
    g(:,:) = 0.
!
!
! Assign g-factors at altitudes which are sunlit:
!
    DO I=1,JMAX
      GH = (RE+ZZ(I)) * SIN(SZA)
      IF (SZA .LT. 1.6 .OR. GH .GT. RE) THEN
        G(1,I) = 0.041
        G(2,I) = 0.013
      ENDIF
    ENDDO
!
!
! Calculate rate coefficients as a function of altitude:
!
    DO I=1,JMAX
      T1(I) = (16.*ZTN(I)+28.*ZTI(I)) / (16.+28.) / 300.
      T2(I) = (16.*ZTN(I)+32.*ZTI(I)) / (16.+32.) / 300.
      T3(I) = (28.*ZTN(I)+16.*ZTI(I)) / (28.+16.) / 300.
      T4(I) = (28.*ZTN(I)+32.*ZTI(I)) / (28.+32.) / 300.
      T5(I) = (16.*ZTN(I)+30.*ZTI(I)) / (16.+30.) / 300.
      IF (T1(I) .LT. 5.6667) THEN
        KZ(1,I) = 1.533E-12 - 5.92E-13*T1(I) + 8.6E-14*T1(I)**2
      ELSE
        KZ(1,I) = 2.73E-12 - 1.155E-12*T1(I) + 1.483E-13*T1(I)**2
      ENDIF
      IF (T2(I) .LT. 6.6667) THEN
        KZ(2,I) = 3.53E-11 - 1.84E-11*T2(I) + 4.62E-12*T2(I)**2 &
                           - 4.95E-13*T2(I)**3 + 2.00E-14*T2(I)**4
      ELSE
        KZ(2,I) = 2.82E-11 - 7.74E-12*T2(I) + 1.073E-12*T2(I)**2 &
                           - 5.17E-14*T2(I)**3 + 9.65E-16*T2(I)**4
      ENDIF
      KZ(3,I) = 1.793e-10 - 6.242e-11*T3(I) + 1.225e-11*T3(I)**2 &
                          - 1.016e-12*T3(I)**3 + 3.087e-14*T3(I)**4
      KZ(4,I) = 5.0E-11 * (1./T4(I)) ** 0.8
      KZ(5,I) = 6.0E-12
      KZ(6,I) = 6.9E-13
      KZ(7,I) = 5.5E-10 * (ZTE(I)/300.) ** 0.5
      KZ(8,I) = 2.0E-11 * EXP(107.8/ZTN(I))
      KZ(9,I) = 2.9E-11 * EXP(67.5 /ZTN(I))
      KZ(10,I) = 8.1E-10 * (ZTE(I)/300.) ** 0.5
      KZ(11,I) = 2.0E-14
      KZ(12,I) = 8.0E-10
      KZ(13,I) = 7.0E-10
      KZ(14,I) = 6.6E-08 * (300./ZTE(I)) ** 0.5
      KZ(15,I) = 1.0E-11
      KZ(16,I) = 4.8E-10
      KZ(17,I) = 4.8E-10
      KZ(18,I) = 1.7E-07 * (300./ZTE(I)) ** 0.5
      KZ(19,I) = 5.2E-11
      KZ(20,I) = 2.1E-11 * EXP(-1136./ZTN(I))
      KZ(21,I) = 1.0E-13
      KZ(22,I) = 4.2E-07 * (300./ZTE(I)) ** 0.85
      KZ(23,I) = 1.8E-07 * (300./ZTE(I)) ** 0.39
      KZ(24,I) = 1.95E-07 * (300./ZTE(I)) ** 0.70
      IF (ZTE(I) .GE. 1200.) KZ(24,I) = 1.6E-07 * (300./ZTE(I)) ** 0.55
      KZ(25,I) = 6.0E-10
      KZ(26,I) = 3.1E-11
      KZ(27,I) = 3.0E-12
      IF (ZTE(I) .LT. 500.) THEN
        KZ(28,I) = 1.0E-29
      ELSE
        KZ(28,I) = 2.6E-11 * ZTE(I)**0.5 * EXP(-22740./ZTE(I))
      ENDIF
      KZ(29,I) = 4.1E-12
      KZ(30,I) = 4.4E-10
      KZ(31,I) = 7.0E-11
      KZ(32,I) = 1.0E-12
      KZ(33,I) = 1.2E-11
      KZ(34,I) = 2.0E-12
      KZ(35,I) = 1.8E-10
      KZ(36,I) = 4.0E-12 * EXP(-865./ZTN(I))
      KZ(37,I) = 1.2E-10
      KZ(38,I) = 1.3E-10
      KZ(39,I) = 2.0E-11
      IF (T5(I) .LT. 5) THEN
        KZ(40,I) = 8.36E-13 - 2.02E-13*T5(I) + 6.95E-14*T5(I)**2
      ELSE
        KZ(40,I) = 5.33E-13 - 1.64E-14*T5(I) + 4.72E-14*T5(I)**2 &
                            - 7.05E-16*T5(I)**3
      ENDIF
      KZ(41,I) = 7.3E-13
      KZ(42,I) = 1.3E-15
      KZ(43,I) = 1.0E-7
      KZ(44,I) = 1.4E-10
      KZ(45,I) = 4.0E-13
    ENDDO
!
!
! Calculate Electron impact ionization, photoionization, and electron
! impact dissociation rates at each altitude; put a priori electron
! density in calculated electron density array: put a priori N(2D) in DEN array:
!
    DO I=1,JMAX
      OEI(I)   = SION(1,I)+PIA(1,I)
      O2EI(I)  = SION(2,I)+PIA(2,I)
      RN2EI(I) = SION(3,I)+PIA(3,I)
      TEI(I)   = OEI(I)+O2EI(I)+RN2EI(I)
      OPI(I)   = PHOTOI(1,1,I)+PHOTOI(2,1,I)+PHOTOI(3,1,I)+PHOTOI(4,1,I)+PHOTOI(5,1,I)
      O2PI(I)  = PHOTOI(1,2,I)+PHOTOI(2,2,I)+PHOTOI(3,2,I)
      RN2PI(I) = PHOTOI(1,3,I)+PHOTOI(2,3,I)+PHOTOI(3,3,I)+PHOTOI(4,3,I)+PHOTOI(5,3,I)
      TPI(I)   = OPI(I)+O2PI(I)+RN2PI(I)+PHOTOI(4,2,I)+PHOTOI(6,3,I)+PHONO(1,I)
      TIR(I)   = TEI(I)+TPI(I)
      RN2ED(I) = AGLW(5,3,I)+AGLW(6,3,I)+AGLW(7,3,I)+B(24)*PIA(3,I)
      SRCED(I) = AGLW(4,2,I) + B(35)*PIA(2,I)
      E(I)     = ZE(I)
      DEN(10,I)= ZND(I)
    ENDDO
!
!
! Find level below which electron density will be calculated:
!
    IF (KCHEM .GE. 4) THEN
      DO I=JMAX,1,-1
        IF (ZZ(I) .GT. 2.0001E7) J200=I-1
      ENDDO
    ELSE
      J200=0
    ENDIF
!
!
! Iterative loop assures that electron density and feedback reactions
! (O+(2P,2D)+e, O+(4S)+N(2D), N2++O) are correctly computed:
!
    DO ITER=1,5
!
!
! Calculate atomic ion densities at each altitude:
!
    DO I=1,JMAX
!
!
! O+(2P):
!
      P(1,I)= PHOTOI(3,1,I) &
            + B(41) * PHOTOI(5,1,I) &
            + B(30) * PHOTOI(4,2,I) &
            + B(9)  * OEI(I) &
            + B(12) * O2EI(I)
      L(1,I)= KZ(16,I) * ZN2(I) &
            + KZ(17,I) * ZO2(I) &
            + KZ(19,I) * ZO(I) &
            + KZ(18,I) * E(I) &
            + A(8) &
            + A(7)
      DEN(1,I) = P(1,I) / L(1,I)
!
!
! O+(2D):
!
      P(2,I)= PHOTOI(2,1,I) &
            + B(42) * PHOTOI(5,1,I) &
            + B(31) * PHOTOI(4,2,I) &
            + B(10) * OEI(I) &
            + B(13) * O2EI(I) &
            + B(8)  * KZ(18,I) * DEN(1,I) * E(I) &
            + A(8)  * DEN(1,I)
      L(2,I)= KZ(12,I) * ZN2(I) &
            + KZ(13,I) * ZO2(I) &
            + KZ(15,I) * ZO(I) &
            + KZ(14,I) * E(I) &
            + A(6)
      DEN(2,I) = P(2,I) / L(2,I)
!
!
! N+:
!
      IF (KCHEM .GE. 2) THEN
        P(4,I) = PHOTOI(6,3,I) &
               + B(15) * RN2EI(I) &
               + KZ(38,I) * DEN(3,I) * DEN(10,I)
        L(4,I) = KZ(25,I) * ZO2(I) &
               + KZ(32,I) * ZO(I)
        DEN(4,I) = P(4,I) / L(4,I)
      ENDIF
!
!
! O+(4S):
!
      IF (KCHEM .GE. 3) THEN
        P(3,I)= PHOTOI(1,1,I) + PHOTOI(4,1,I) &
              + B(32) * PHOTOI(4,2,I) &
              + B(11) * OEI(I) &
              + B(14) * O2EI(I)  &
              + KZ(14,I) * DEN(2,I) * E(I) & 
              + KZ(15,I) * DEN(2,I) * ZO(I) & 
              + A(6) * DEN(2,I) &
              + (1.-B(8)) * KZ(18,I) * DEN(1,I) * E(I) &
              + KZ(19,I) * DEN(1,I) * ZO(I) & 
              + A(7) * DEN(1,I) &
              + KZ(32,I) * DEN(4,I) * ZO(I) &
              + KZ(39,I) * DEN(5,I) * ZO(I)
        L(3,I)= KZ(1,I) * ZN2(I) &
              + KZ(2,I) * ZO2(I) &
              + KZ(38,I) * DEN(10,I)
        DEN(3,I) = P(3,I) / L(3,I)
      ENDIF
!
    ENDDO    ! bottom of atomic ion loop
!
!
! Above 200 km, (or at all altitudes if KCHEM=3) use a priori
! electron density to calculate O+(4S):
!
    IF (KCHEM .GE. 3) THEN
!
      DO I=J200+1,JMAX
        P(5,I)= RN2PI(I) &
                + (1.-B(15)) * RN2EI(I) &
                + KZ(12,I) * DEN(2,I) * ZN2(I) &
                + KZ(16,I) * DEN(1,I) * ZN2(I)
        L(5,I)= KZ(3,I)  * ZO(I) &
                + KZ(4,I)  * ZO2(I) &
                + KZ(23,I) * E(I) &
                + KZ(39,I) * ZO(I)
        DEN(5,I) = P(5,I) / L(5,I)
        QQ(I) = PHONO(1,I) &
                + KZ(3,I)  * DEN(5,I) * ZO(I) &
                + B(21) * KZ(25,I) * DEN(4,I) * ZO2(I)
        RR(I) = KZ(30,I) * ZNO(I) &
                + KZ(37,I) * ZNS(I)
        SS(I) = KZ(1,I) * ZN2(I) &
                + KZ(40,I) * ZNO(I)
        TT(I) = KZ(22,I) * E(I)
        UU(I) = O2PI(I) &
                + (1.-B(12)-B(13)-B(14)) * O2EI(I) &
                + KZ(13,I) * DEN(2,I) * ZO2(I) &
                + KZ(17,I) * DEN(1,I) * ZO2(I) &
                + KZ(4,I)  * DEN(5,I) * ZO2(I) &
                + B(22) * KZ(25,I) * DEN(4,I) * ZO2(I)
        VV(I) = KZ(2,I) * ZO2(I)
        WW(I) = KZ(24,I) * E(I) &
                + KZ(30,I) * ZNO(I) &
                + KZ(37,I) * ZNS(I)
        XX(I) = DEN(1,I) + DEN(2,I) + DEN(4,I) + DEN(5,I)
        DEN(3,I) = (TT(I)*WW(I)*E(I) - TT(I)*WW(I)*XX(I) - TT(I)*UU(I) &
                   - QQ(I)*WW(I) - RR(I)*UU(I) ) / &
                      (TT(I)*WW(I) + TT(I)*VV(I) + RR(I)*VV(I) + SS(I)*WW(I))
        if (den(3,i) .lt. 0.) den(3,i)=0.
      ENDDO
!
    ENDIF
!
!
! If KCHEM=4, calculate electron density below 200 km using iterative method:
! First time: approximate electron density using effective recombination rate.
! Subsequent iterations: update electron density from sum of ions.
!
    if (kchem .ge. 4) then
!
      if (iter .eq. 1) then
        do i=1,j200
          tatomi=den(1,i)+den(2,i)+den(3,i)+den(4,i)
          alphaef=(kz(22,i)+kz(24,i))/2.
          e(i)=(tatomi+sqrt(tatomi**2+4.*tir(i)/alphaef))/2.
        enddo
      else
        do i=1,j200
          toti = den(1,i)+den(2,i)+den(3,i)+den(4,i) &
                +den(5,i)+den(6,i)+den(7,i)
          e(i) = (toti + e(i)) / 2.
        enddo
      endif
!
!
! Smoothly transition to electron density above 200 km:
!
      E(J200+1) = E(J200) * ( E(J200+3) / E(J200) ) &
                  ** ( (ZZ(J200+1)-ZZ(J200)) / (ZZ(J200+3)-ZZ(J200)) )
      E(J200+2) = E(J200) * (E(J200+3)/E(J200)) &
                  ** ( (ZZ(J200+2)-ZZ(J200)) / (ZZ(J200+3)-ZZ(J200)) )
!
    endif
!
!
! Calculate molecular ion densities and excited species densites:
!
    DO I=1,JMAX
!
!
! N2+:
!
      IF (KCHEM .GE. 2) THEN
        P(5,I)= RN2PI(I) &
              + (1.-B(15)) * RN2EI(I) &
              + KZ(12,I) * DEN(2,I) * ZN2(I) &
              + KZ(16,I) * DEN(1,I) * ZN2(I)
        L(5,I)= KZ(3,I)  * ZO(I) &
              + KZ(4,I)  * ZO2(I) &
              + KZ(23,I) * E(I) &
              + KZ(39,I) * ZO(I)
          DEN(5,I) = P(5,I) / L(5,I)
      ENDIF
!
!
! O2+:
!
      IF (KCHEM .GE. 3) THEN
        P(6,I)= O2PI(I) &
              + (1.-B(12)-B(13)-B(14)) * O2EI(I) &
              + KZ(2,I)  * DEN(3,I) * ZO2(I) &
              + KZ(13,I) * DEN(2,I) * ZO2(I) &
              + KZ(17,I) * DEN(1,I) * ZO2(I) &
              + KZ(4,I)  * DEN(5,I) * ZO2(I) &
              + B(22) * KZ(25,I) * DEN(4,I) * ZO2(I)
        L(6,I)= KZ(24,I) * E(I) &
              + KZ(30,I) * ZNO(I) &
              + KZ(37,I) * ZNS(I)
        DEN(6,I) = P(6,I)/ L(6,I)
      ENDIF
!
!
! NO+:
!
      IF (KCHEM .GE. 3) THEN
        P(7,I)= PHONO(1,I) &
              + KZ(1,I)  * DEN(3,I) * ZN2(I) &
              + KZ(40,I) * DEN(3,I) * ZNO(I) &
              + KZ(3,I)  * DEN(5,I) * ZO(I) &
              + B(21) * KZ(25,I) * DEN(4,I) * ZO2(I) &
              + KZ(30,I) * DEN(6,I) * ZNO(I) &
              + KZ(37,I) * DEN(6,I) * ZNS(I)
        L(7,I)= KZ(22,I) * E(I)
        DEN(7,I) = P(7,I) / L(7,I)
      ENDIF
!
!
! N2(A):
!
      P(8,I)= AGLW(1,3,I) + AGLW(2,3,I) + B(43)*AGLW(3,3,I)
      L(8,I)= KZ(26,I) * ZO(I) &
            + KZ(29,I) * ZO2(I) &
            + A(10)
      DEN(8,I) = P(8,I) / L(8,I)
!
!
! N(2P):
!
      P(9,I)= B(28) * PHOTOD(1,3,I) &
            + B(28) * PHOTOI(6,3,I) &
            + B(26) * RN2ED(I) &
            + B(23) * KZ(23,I) * DEN(5,I) * E(I)
      L(9,I)= KZ(33,I) * ZO(I) &
            + KZ(34,I) * ZO2(I) &
            + KZ(35,I) * ZNO(I) &
            + A(11) &
            + A(12)
      DEN(9,I) = P(9,I) / L(9,I)
!
!
! N(2D):
!
      P(10,I)= B(27) * PHOTOD(1,3,I) &
             + B(27) * PHOTOI(6,3,I) &
             + B(25) * RN2ED(I) &
             + B(16) * B(15) * RN2EI(I) &
             + B(3)  * KZ(22,I) * DEN(7,I) * E(I) &
             + B(4)  * KZ(23,I) * DEN(5,I) * E(I) &
             + B(5)  * KZ(3,I)  * DEN(5,I) * ZO(I) &
             + B(29) * KZ(33,I) * DEN(9,I) * ZO(I) &
             + A(12) * DEN(9,I)
      L(10,I)= KZ(5,I)  * ZO2(I) &
             + KZ(6,I)  * ZO(I) &
             + KZ(7,I)  * E(I) &
             + KZ(31,I) * ZNO(I) &
             + KZ(38,I) * DEN(3,I) &
             + A(1)
      DEN(10,I) = P(10,I) / L(10,I)
!
!
! O(1S):
!
      BZ(1,I) = 0.12 + 0.02 * ALOG10 (E(I)/ZO(I)*(300./ZTE(I))**0.7)
      IF (BZ(1,I) .LT. 0.03) BZ(1,I)=0.03
      P(11,I)= AGLW(2,1,I) &
             + BZ(1,I) * KZ(24,I) * DEN(6,I)  * E(I) &
             + B(18) * KZ(26,I) * DEN(8,I) * ZO(I) &
             + B(33) * KZ(37,I) * DEN(6,I) * ZNS(I) &
             + B(34) * KZ(31,I) * DEN(10,I) * ZNO(I) &
             + PHOTOD(2,2,I)
      L(11,I)= KZ(11,I) * ZO(I) &
             + KZ(36,I) * ZO2(I) &
             + A(5) &
             + A(4)
      DEN(11,I) = P(11,I) / L(11,I)
!
!
! O(1D):
!
      P(12,I)= AGLW(1,1,I) &
             + KZ(28,I) * E(I)  * ZO(I) &
             + B(2)  * KZ(24,I) * DEN(6,I)  * E(I) &
             + B(6)  * KZ(5,I)  * DEN(10,I) * ZO2(I) &
             + B(20) * KZ(6,I)  * DEN(10,I) * ZO(I) &
             + B(17) * KZ(25,I) * DEN(4,I)  * ZO2(I) &
             + B(7)  * KZ(15,I) * DEN(2,I)  * ZO(I) &
             + SRCED(I) &
             + PHOTOD(1,2,I) &
             + A(5)  * DEN(11,I)
      L(12,I)= KZ(8,I)  * ZN2(I)  &
             + KZ(9,I)  * ZO2(I) &
             + KZ(10,I) * E(I) &
             + KZ(27,I) * ZO(I) &
             + A(2) &
             + A(3)
      DEN(12,I) = P(12,I) / L(12,I)
!
    ENDDO   ! bottome of molecular ion / excited species loop
!
    ENDDO   ! bottom of iterative looop
!
!
! Impose charge neutrality:
!
    do i=1,j200
      e(i)=den(1,i)+den(2,i)+den(3,i)+den(4,i)+den(5,i)+den(6,i)+den(7,i)
    enddo
!
!
! Calculate O- for mutual neutralization source of O*
!
    do i=1,jmax
      ominus(i) = (kz(42,i)*zo(i)*e(i)) / (kz(43,i)*den(3,i)+kz(44,i)*zo(i))
    enddo
!
!
! Calculate airglow emission rates; fill ZCETA array with partial rates
! from each source; fill ZETA array with total rate for each emission:
!
    DO I=1,JMAX
!
      ZCETA(1,1,I) = B(39) * AGLW(3,3,I)
      ZCETA(2,1,I) = B(40) * A(10) * P(8,I) / L(8,I)
!
      ZCETA(1,2,I) = B(38) * B(36) * RN2EI(I)
      ZCETA(2,2,I) = B(38) * PHOTOI(3,3,I)
      ZCETA(3,2,I) = G(2,I) * DEN(5,I)
!
      ZCETA(1,3,I) = A(1) * B(27) * PHOTOD(1,3,I) / L(10,I)
      ZCETA(2,3,I) = A(1) * B(27) * PHOTOI(6,3,I) / L(10,I)
      ZCETA(3,3,I) = A(1) * B(25) * RN2ED(I) / L(10,I)
      ZCETA(4,3,I) = A(1) * B(16) * B(15) * RN2EI(I) / L(10,I)
      ZCETA(5,3,I) = A(1) * B(3)  * KZ(22,I) * DEN(7,I) * E(I) /L(10,I)
      ZCETA(6,3,I) = A(1) * B(4)  * KZ(23,I) * DEN(5,I) * E(I) /L(10,I)
      ZCETA(7,3,I) = A(1) * B(5)  * KZ(3,I)  * DEN(5,I) * ZO(I) /L(10,I)
      ZCETA(8,3,I) = A(1) * B(29) * KZ(33,I) * DEN(9,I) * ZO(I) /L(10,I)
      ZCETA(9,3,I) = A(1) * A(12) * DEN(9,I) / L(10,I)
!
      ZCETA(1,4,I) = A(5) * AGLW(2,1,I) / L(11,I)
      ZCETA(2,4,I) = A(5) * BZ(1,I)*KZ(24,I) * DEN(6,I) * E(I)  /L(11,I)
      ZCETA(3,4,I) = A(5) * B(18) * KZ(26,I) * DEN(8,I) * ZO(I) /L(11,I)
      ZCETA(4,4,I) = A(5) * B(33) * KZ(37,I) * DEN(6,I) * ZNS(I)/L(11,I)
      ZCETA(5,4,I) = A(5) * B(34) * KZ(31,I) * DEN(10,I)* ZNO(I)/L(11,I)
      ZCETA(6,4,I) = PHOTOD(2,2,I) / L(11,I)
!
      ZCETA(1,5,I) = A(2) * AGLW(1,1,I) / L(12,I)
      ZCETA(2,5,I) = A(2) * KZ(28,I) * E(I)  * ZO(I) / L(12,I)
      ZCETA(3,5,I) = A(2) * B(2)  * KZ(24,I) * DEN(6,I)  * E(I)/L(12,I)
      ZCETA(4,5,I) = A(2) * B(6)  * KZ(5,I)  * DEN(10,I) *ZO2(I)/L(12,I)
      ZCETA(5,5,I) = A(2) * B(20) * KZ(6,I)  * DEN(10,I) * ZO(I)/L(12,I)
      ZCETA(6,5,I) = A(2) * B(17) * KZ(25,I) * DEN(4,I) * ZO2(I)/L(12,I)
      ZCETA(7,5,I) = A(2) * B(7)  * KZ(15,I) * DEN(2,I)  * ZO(I)/L(12,I)
      ZCETA(8,5,I) = A(2) * SRCED(I) / L(12,I)
      ZCETA(9,5,I) = A(2) * PHOTOD(1,2,I) / L(12,I)
      ZCETA(10,5,I)= A(2) * A(5)  * DEN(11,I) / L(12,I)
!
      ZCETA(1,6,I) = A(8) * (PHOTOI(3,1,I)+B(41)*PHOTOI(5,1,I)) / L(1,I)
      ZCETA(2,6,I) = A(8) * B(30) * PHOTOI(4,2,I) / L(1,I)
      ZCETA(3,6,I) = A(8) * B(9) * OEI(I) / L(1,I)
      ZCETA(4,6,I) = A(8) * B(12) * O2EI(I) / L(1,I)
!
      ZCETA(1,7,I) = A(12) * B(28) * PHOTOD(1,3,I) / L(9,I)
      ZCETA(2,7,I) = A(12) * B(28) * PHOTOI(6,3,I) / L(9,I)
      ZCETA(3,7,I) = A(12) * B(26) * RN2ED(I) / L(9,I)
      ZCETA(4,7,I) = A(12) * B(23) * KZ(23,I) * DEN(5,I) * E(I) / L(9,I)
!
      ZCETA(1,8,I) = A(11) * B(28) * PHOTOD(1,3,I) / L(9,I)
      ZCETA(2,8,I) = A(11) * B(28) * PHOTOI(6,3,I) / L(9,I)
      ZCETA(3,8,I) = A(11) * B(26) * RN2ED(I) / L(9,I)
      ZCETA(4,8,I) = A(11) * B(23) * KZ(23,I) * DEN(5,I) * E(I) / L(9,I)
!
      ZCETA(1,9,I) = AGLW(5,1,I)
      ZCETA(2,9,I) = KZ(41,I) * DEN(3,I) * E(I)
      ZCETA(3,9,I) = B(49) * KZ(43,I) * OMINUS(I) * DEN(3,I)
!
      ZCETA(1,10,I) = AGLW(6,1,I)
      ZCETA(2,10,I) = AGLW(7,1,I)
      ZCETA(3,10,I) = B(44) * AGLW(8,1,I)
!
      ZCETA(1,11,I) = A(6) * (PHOTOI(2,1,I)+B(42)*PHOTOI(5,1,I))/ L(2,I)
      ZCETA(2,11,I) = A(6) * B(31) * PHOTOI(4,2,I) / L(2,I)
      ZCETA(3,11,I) = A(6) * B(10) * OEI(I) / L(2,I)
      ZCETA(4,11,I) = A(6) * B(13) * O2EI(I) / L(2,I)
      ZCETA(5,11,I) = A(6) * B(8)  * KZ(18,I) * DEN(1,I) * E(I) / L(2,I)
      ZCETA(6,11,I) = A(6) * A(8)  * DEN(1,I) / L(2,I)
!
      ZCETA(1,12,I) = AGLW(4,3,I) * B(48)
!
      ZCETA(1,13,I) = AGLW(3,1,I)
      ZCETA(2,13,I) = AGLW(5,1,I)
      ZCETA(3,13,I) = KZ(41,I) * DEN(3,I) * E(I)
      ZCETA(4,13,I) = B(49) * KZ(43,I) * OMINUS(I) * DEN(3,I)
!
      ZCETA(1,14,I) = B(46)*PHOTOI(6,3,I)
      ZCETA(2,14,I) = B(47)*B(15)*RN2EI(I)
!
      ZCETA(1,15,I) = AGLW(4,1,I)
      ZCETA(2,15,I) = AGLW(6,1,I)
      ZCETA(3,15,I) = AGLW(7,1,I)
      ZCETA(4,15,I) = B(44) * AGLW(8,1,I)
      ZCETA(5,15,I) = KZ(45,I) * DEN(3,I) * E(I)
!
      ZETA(1,I)  = ZCETA(1,1,I)+ZCETA(2,1,I)
      ZETA(2,I)  = ZCETA(1,2,I)+ZCETA(2,2,I)+ZCETA(3,2,I)
      ZETA(3,I)  = ZCETA(1,3,I)+ZCETA(2,3,I)+ZCETA(3,3,I) &
                  +ZCETA(4,3,I)+ZCETA(5,3,I)+ZCETA(6,3,I) &
                  +ZCETA(7,3,I)+ZCETA(8,3,I)+ZCETA(9,3,I)
      ZETA(4,I)  = ZCETA(1,4,I)+ZCETA(2,4,I)+ZCETA(3,4,I) &
                  +ZCETA(4,4,I)+ZCETA(5,4,I)+ZCETA(6,4,I)
      ZETA(5,I)  = ZCETA(1,5,I)+ZCETA(2,5,I)+ZCETA(3,5,I) &
                  +ZCETA(4,5,I)+ZCETA(5,5,I)+ZCETA(6,5,I) &
                  +ZCETA(7,5,I)+ZCETA(8,5,I)+ZCETA(9,5,I) &
                  +ZCETA(10,5,I)
      ZETA(6,I)  = ZCETA(1,6,I)+ZCETA(2,6,I)+ZCETA(3,6,I)+ZCETA(4,6,I)
      ZETA(7,I)  = ZCETA(1,7,I)+ZCETA(2,7,I)+ZCETA(3,7,I)+ZCETA(4,7,I)
      ZETA(8,I)  = ZCETA(1,8,I)+ZCETA(2,8,I)+ZCETA(3,8,I)+ZCETA(4,8,I)
      ZETA(9,I)  = ZCETA(1,9,I)+ZCETA(2,9,I)+ZCETA(3,9,I)
      ZETA(10,I) = ZCETA(1,10,I)+ZCETA(2,10,I)+ZCETA(3,10,I)
      ZETA(11,I) = ZCETA(1,11,I)+ZCETA(2,11,I)+ZCETA(3,11,I) &
                  +ZCETA(4,11,I)+ZCETA(5,11,I)+ZCETA(6,11,I)
      ZETA(12,I) = ZCETA(1,12,I)
      ZETA(13,I) = ZCETA(1,13,I)+ZCETA(2,13,I)+ZCETA(3,13,I)+ZCETA(4,13,I)
      ZETA(14,I) = ZCETA(1,14,I)+ZCETA(2,14,I)
      ZETA(15,I) = ZCETA(1,15,I)+ZCETA(2,15,I)+ZCETA(3,15,I) &
                  +ZCETA(4,15,I)+ZCETA(5,15,I)
!
    ENDDO ! bottom of airglow loop
!
!
! Calculate vertical column brightnesses:
!
    DO I=1,JMAX
      IF (I .EQ. JMAX) THEN
        DZ = (ZZ(I) - ZZ(I-1))
      ELSE
        IF (I .EQ. 1) THEN
          DZ = (ZZ(I+1) - ZZ(I))
        ELSE
          DZ = (ZZ(I+1) - ZZ(I-1)) / 2.0
        ENDIF
      ENDIF
      DO IW=1,NW
        VCB(IW) = VCB(IW) + ZETA(IW,I) * DZ
      ENDDO
    ENDDO
!
!
! Convert brightnesses to Rayleighs:
!
    DO IW=1,NW
      VCB(IW) = VCB(IW) / 1.E6
    ENDDO
!
!
    RETURN
!
  END SUBROUTINE GCHEM

