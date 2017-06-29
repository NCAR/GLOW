! Subroutine RCOLUM

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Stan Solomon, 1988, 1991
! Stan Solomon, 2016: removed problematic extrapolation below lower
! boundary.  If grazing height is below lower boundary of atmosphere
! supplied, column density is set to 1.0e30.
! Stan Solomon, 2016: refactored for f90.
! Stan Solomon, 2017: corrected chi to pi-chi at line 79.  Note this has no effect.

! Calculates the column density ZCOL for each species ZMAJ above height
! ZZ at zenith angle CHI.  Calls subroutine VCD to calculate the
! vertical column density, and then uses a fit to the Chapman Grazing
! Incidence Integral [Smith and Smith, JGR 77, 3592, 1972] to calculate 
! the slant column density.  If CHI is less than 90 degrees, column
! densities are calculated directly; if CHI is greater than 90 degrees
! the column density at grazing height for 90 degrees is calculated and
! doubled, and the column density above ZZ(J) is subtracted.  If the
! grazing height is lower than the bottom of the atmosphere supplied, 
! column densities are set to 'infinity', i.e., 1.0e30.


    subroutine rcolum (chi, zz, zmaj, tn, zcol, zvcd, jmax, nmaj)

      implicit none

      integer,intent(in) :: jmax, nmaj
      real,intent(in) :: chi, zz(jmax), zmaj(nmaj,jmax), tn(jmax)
      real,intent(out) :: zcol(nmaj,jmax), zvcd(nmaj,jmax)

      integer,parameter :: nm=3
      real,parameter :: pi=3.1415926536
      real,parameter :: re=6.37e8

      integer :: i, j, k
      real :: zcg(nm), ghrg, ghz, tng
      real,external :: chap

      call vcd (zz, zmaj, zvcd, jmax, nmaj)

      if (chi >= 2.) then 
        do i=1,nmaj
          do j=1,jmax
            zcol(i,j) = 1.0e30
          enddo
        enddo
        return
      endif

      if (chi <= pi/2.) then
        do i=1,nmaj
          do j=1,jmax
            zcol(i,j) = zvcd(i,j) * chap(chi,zz(j),tn(j),i)
          enddo
        enddo
      else
        do j=1,jmax
          ghrg=(re+zz(j))*sin(chi) 
          ghz=ghrg-re 
          if (ghz <= zz(1)) then
            do i=1,nmaj
              zcol(i,j) = 1.0e30
            enddo
          else
            do k=1,j-1
              if (zz(k) <= ghz .and. zz(k+1) > ghz) then
                tng = tn(k)+(tn(k+1)-tn(k))*(ghz-zz(k))/(zz(k+1)-zz(k))
                do i=1,nmaj
                  zcg(i) = zvcd(i,k) * (zvcd(i,k+1) / zvcd(i,k)) ** &
                           ((ghz-zz(k)) / (zz(k+1)-zz(k)))
                enddo
              endif
            enddo
            do i=1,nmaj
              zcol(i,j) = 2. * zcg(i) * chap(pi/2.,ghz,tng,i) &
                        - zvcd(i,j) * chap(pi-chi,zz(j),tn(j),i)
            enddo
          endif
        enddo
      endif

      return 

    end subroutine rcolum

!----------------------------------------------------------------------

    real function chap (chi, z, t, i)

      implicit none

      real,intent(in) :: chi, z, t
      integer,intent(in) :: i

      integer,parameter ::  nmaj=3
      real,parameter :: pi=3.1415926536
      real,parameter :: re=6.37e8
      real,parameter :: g=978.1

      real :: am(nmaj), gr, hn, hg, hf, sqhf
      real,external :: sperfc

      data am/16., 32., 28./

      gr=g*(re/(re+z))**2 
      hn=1.38e-16*t/(am(i)*1.662e-24*gr)
      hg=(re+z)/hn 
      hf=0.5*hg*(cos(chi)**2) 
      sqhf=sqrt(hf) 
      chap=sqrt(0.5*pi*hg)*sperfc(sqhf) 

      return

    end function chap

!----------------------------------------------------------------------

    real function sperfc(dummy) 

      implicit none

      real,intent(in) :: dummy

      if (dummy <= 8.) then
        sperfc = (1.0606963+0.55643831*dummy) / &
                 (1.0619896+1.7245609*dummy+dummy*dummy)
      else
        sperfc=0.56498823/(0.06651874+dummy) 
      endif 

      return 

    end function sperfc

!----------------------------------------------------------------------

    subroutine vcd(zz,zmaj,zvcd,jmax,nmaj)

      implicit none

      integer,intent(in) :: jmax,nmaj
      real,intent(in) :: zz(jmax), zmaj(nmaj,jmax)

      real,intent(out) :: zvcd(nmaj,jmax)

      integer :: i, j
      real :: rat

      do i=1,nmaj
        zvcd(i,jmax) =   zmaj(i,jmax) &
                       * (zz(jmax)-zz(jmax-1)) &
                       / alog(zmaj(i,jmax-1)/zmaj(i,jmax))
        do j=jmax-1,1,-1
          rat = zmaj(i,j+1) / zmaj(i,j)
          zvcd(i,j) = zvcd(i,j+1)+zmaj(i,j)*(zz(j)-zz(j+1))/alog(rat)*(1.-rat)
        enddo
      enddo

      return

    end subroutine vcd
