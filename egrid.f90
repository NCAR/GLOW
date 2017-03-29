! Subroutine EGRID sets up electron energy grid

! This software is part of the GLOW model.  Use is governed by the open source
! academic research license agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Stan Solomon, 1/1992
! Refactored to f90, SCS, 6/2016

! Inputs:
!   nbins   number of bins in the electron energy grid
! Outputs:
!   ener    energy at center of each bin, eV
!   del     width of each bin, eV

 
subroutine egrid (ener, del, nbins)

  implicit none

  integer,intent(in) :: nbins
  real,intent(out) :: ener(nbins), del(nbins)

  integer :: n

  do n=1,nbins
    if (n <= 21) then
      ener(n) = 0.5 * float(n)
    else
      ener(n) = exp (0.05 * float(n+26))
    endif
  enddo

  del(1) = 0.5

  do n=2,nbins
    del(n) = ener(n)-ener(n-1)
  enddo

  do n=1,nbins
    ener(n) = ener(n) - del(n)/2.0
  enddo

  return

end subroutine egrid
