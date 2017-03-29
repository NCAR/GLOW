! Subroutine MAXT
!
! This software is part of the glow model.  Use is governed by the open source
! academic research license agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.
!
! Stan Solomon, 11/1989, 9/1991, 1/1994, 3/2005
! Refactored to f90, 6/6/6/6/6/6/2016
!
! Generates Maxwellian electron spectra with, optionally, a low energy tail
! of the form used by Meier et al., JGR 94, 13541, 1989.
! Also can generate a monoenergetic flux in a single bin.
!
! Supplied by calling routine:
!     eflux  total energy flux in erg cm-2 s-1
!     ezer   characteristic energy in ev
!     ener   energy grid in ev
!     del    energy bin width in ev
!     nbins  number of energy bins (dimension of ener, del, and phi)
!     itail  1 = maxwellian with low-energy tail, 0 = regular maxwellian
!     fmono  additional monoenergetic energy flux in erg cm-2 s-1
!     emono  characteristic enerngy of fmono in ev
!
! Returned by subroutine:
!     phi    hemispherical flux in cm-2 s-1 ev-1


  subroutine maxt (eflux,ezer,ener,del,nbins,itail,fmono,emono,phi)

    implicit none

    integer,intent(in) :: nbins, itail
    real,intent(in) :: eflux, ezer, ener(nbins), del(nbins), fmono, emono
    real,intent(out) :: phi(nbins)

    integer :: k
    real :: te, b, phimax, erat

    te = 0.

    if (ezer < 500.) then
      b = 0.8*ezer
    else
      b = 0.1*ezer + 350.
    endif

    phimax = exp(-1.)

    do k=1,nbins
      erat = ener(k) / ezer
      if (erat > 60.) erat = 60.
      phi(k) = erat * exp(-erat)
      if (itail > 0) phi(k) = phi(k) + 0.4*phimax*(ezer/ener(k))*exp(-ener(k)/b)
      te = te + phi(k) * del(k) * ener(k) * 1.6022e-12
    enddo

    do k=1,nbins
      phi(k) = phi(k) * eflux / te
    enddo

    if (fmono > 0.) then
      do k=1,nbins
        if (emono > ener(k)-del(k)/2. .and. emono < ener(k)+del(k)/2.) &
          phi(k)=phi(k)+fmono/(1.6022e-12*del(k)*ener(k))
      enddo
    endif

    return

  end subroutine maxt
