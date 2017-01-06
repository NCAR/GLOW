! Subroutine TZGRID maps fields from a TIE-GCM or TIME-GCM history onto the
! GLOW altitude grid, interpolating and/or extrapolating as necessary.

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Stan Solomon, 12/15, 1/16
! Extracted from glowdriver.f90 into separate file tzgrid.f90, SCS, 12/16

! Inputs:
!         i      Longitude index
!         l      Latitude index
!         jmax   Number of altitude levels used by GLOW (should = 68 for TIE/TIME runs)
! Outputs:
!         z      Geographic altitude of pressure surface (km)
!         zo     O number density, cm-3
!         zo2    O2    "
!         zn2    N2    "
!         zns    N(4S) "
!         n2d    N(2D) "      (optional - could be zero since this is calculated by GLOW)
!         no     NO    "
!         ztn    Tn, K
!         zun    Zonal wind velocity, cm/s (not used by GLOW but passed into output))
!         zvn    Meridional wind velocity, cm/s (not used by GLOW but passed into output))
!         zti    Ti, K
!         zte    Te, K

! GLOW altitude grid extends from lev=ln(P0/P) -10 to +6.75, in intervals of H/4.
! Lower boundary extrapolation for TIE-GCM from lev =-10 to -7 assumes H = 6 km
! low-res TIE-GCM: 29 levels; extrapolate lower boundary, use both midpoints and interfaces
! high-res TIE-GCM: 57 levels; extrapolate lower boundary, use interfaces
! low-res TIME-GCM: 49 levels; use both midpoints and interfaces
! high-res TIME-GCM: 97 levels; use interfaces

subroutine tzgrid(i,l,jmax,z,zo,zo2,zn2,zns,znd,zno,ztn,zun,zvn,ze,zti,zte)

  use readtgcm,only: nlev,zg,tn,un,vn,o2,o1,n2,n4s,n2d,no,ti,te,ne

  implicit none

  integer,intent(in) :: i,l,jmax
  real,intent(out) :: z(jmax),zo(jmax),zo2(jmax),zn2(jmax),zns(jmax),znd(jmax), &
                      zno(jmax),ztn(jmax),zti(jmax),zte(jmax),zun(jmax),zvn(jmax),ze(jmax)
  integer :: j,jj

  if (nlev/=29 .and. nlev/=57 .and. nlev/=49 .and. nlev/=97) then
    write(6,"('zgrid: unknown NLEV = ',i5)") nlev
    stop 'zgrid'
  endif

  if (nlev == 29) then                              ! low-res TIE-GCM
    do j=1,jmax
      jj=(j-11)/2
      if (j <= 13) then                             ! at or below lower boundary
        z(j)   = zg(i,l,1)-float(13-j)*6./4.        ! extrapolate
        zo(j)  = o1(i,l,1)*exp(-float(14-j)/4.)
        zo2(j) = o2(i,l,1)*exp(float(14-j)/4.)
        zn2(j) = n2(i,l,1)*exp(float(14-j)/4.)
        zns(j) = n4s(i,l,1)
        znd(j) = n2d(i,l,1)
        zno(j) = no(i,l,1)
        ztn(j) = tn(i,l,1)
        zti(j) = ti(i,l,1)
        zte(j) = te(i,l,1)
        zun(j) = un(i,l,1)
        zvn(j) = vn(i,l,1)
        ze(j)  = ne(i,l,1)
      else
        if (j/2*2 == j) then                        ! at midpoint
          z(j)   = (zg(i,l,jj)+zg(i,l,jj+1))/2      ! interpolate zg to midpoint
          zo(j)  = o1(i,l,jj)                       ! use fields at midpoint
          zo2(j) = o2(i,l,jj)
          zn2(j) = n2(i,l,jj)
          zns(j) = n4s(i,l,jj)
          znd(j) = n2d(i,l,jj)
          zno(j) = no(i,l,jj)
          ztn(j) = tn(i,l,jj)
          zti(j) = ti(i,l,jj)
          zte(j) = te(i,l,jj)
          zun(j) = un(i,l,jj)
          zvn(j) = vn(i,l,jj)
          ze(j)  = sqrt(ne(i,l,jj)*ne(i,l,jj+1))
        else
          z(j)   = zg(i,l,jj)                       ! at interface
          zo(j)  = sqrt(o1(i,l,jj-1)*o1(i,l,jj))    ! interpolate fields to interface
          zo2(j) = sqrt(o2(i,l,jj-1)*o2(i,l,jj))
          zn2(j) = sqrt(n2(i,l,jj-1)*n2(i,l,jj))
          zns(j) = sqrt(n4s(i,l,jj-1)*n4s(i,l,jj))
          znd(j) = sqrt(n2d(i,l,jj-1)*n2d(i,l,jj))
          zno(j) = sqrt(no(i,l,jj-1)*no(i,l,jj))
          ztn(j) = (tn(i,l,jj-1)+tn(i,l,jj))/2.
          zti(j) = (ti(i,l,jj-1)+ti(i,l,jj))/2.
          zte(j) = (te(i,l,jj-1)+te(i,l,jj))/2.
          zun(j) = (un(i,l,jj-1)+un(i,l,jj))/2.
          zvn(j) = (vn(i,l,jj-1)+vn(i,l,jj))/2.
          ze(j)  = ne(i,l,jj)
        endif
      endif
    enddo 
    z(jmax) = z(jmax-1) + (z(jmax-1)-z(jmax-2))     !patch top level
    ze(jmax) = ze(jmax-1)**2 / ze(jmax-2)
  endif

  if (nlev == 57) then                              ! high-res TIE-GCM
    do j=1,jmax
      jj=j-12
      if (j <= 13) then                             ! at or below lower boundary
        z(j)   = zg(i,l,1)-float(13-j)*6./4.        ! extrapolate
        zo(j)  = o1(i,l,1)*exp( 0.125-float(14-j)/4.)
        zo2(j) = o2(i,l,1)*exp(-0.125+float(14-j)/4.)
        zn2(j) = n2(i,l,1)*exp(-0.125+float(14-j)/4.)
        zns(j) = n4s(i,l,1)
        znd(j) = n2d(i,l,1)
        zno(j) = no(i,l,1)
        ztn(j) = tn(i,l,1)
        zti(j) = ti(i,l,1)
        zte(j) = te(i,l,1)
        zun(j) = un(i,l,1)
        zvn(j) = vn(i,l,1)
        ze(j)  = ne(i,l,1)
      else
          z(j)   = zg(i,l,jj)                       ! at interface
          zo(j)  = sqrt(o1(i,l,jj-1)*o1(i,l,jj))    ! interpolate fields to interface
          zo2(j) = sqrt(o2(i,l,jj-1)*o2(i,l,jj))
          zn2(j) = sqrt(n2(i,l,jj-1)*n2(i,l,jj))
          zns(j) = sqrt(n4s(i,l,jj-1)*n4s(i,l,jj))
          znd(j) = sqrt(n2d(i,l,jj-1)*n2d(i,l,jj))
          zno(j) = sqrt(no(i,l,jj-1)*no(i,l,jj))
          ztn(j) = (tn(i,l,jj-1)+tn(i,l,jj))/2.
          zti(j) = (ti(i,l,jj-1)+ti(i,l,jj))/2.
          zte(j) = (te(i,l,jj-1)+te(i,l,jj))/2.
          zun(j) = (un(i,l,jj-1)+un(i,l,jj))/2.
          zvn(j) = (vn(i,l,jj-1)+vn(i,l,jj))/2.
          ze(j)  = ne(i,l,jj)
      endif
    enddo 
  endif

  if (nlev == 49) then                              ! low-res TIME-GCM
    do j=1,jmax
        jj=(j+29)/2
        if (j/2*2 == j) then                        ! at midpoint
          z(j)   = (zg(i,l,jj)+zg(i,l,jj+1))/2      ! interpolate zg to midpoint
          zo(j)  = o1(i,l,jj)                       ! use fields at midpoint
          zo2(j) = o2(i,l,jj)
          zn2(j) = n2(i,l,jj)
          zns(j) = n4s(i,l,jj)
          znd(j) = n2d(i,l,jj)
          zno(j) = no(i,l,jj)
          ztn(j) = tn(i,l,jj)
          zti(j) = ti(i,l,jj)
          zte(j) = te(i,l,jj)
          zun(j) = un(i,l,jj)
          zvn(j) = vn(i,l,jj)
          ze(j)  = sqrt(ne(i,l,jj)*ne(i,l,jj+1))
        else
          z(j)   = zg(i,l,jj)                       ! at interface
          zo(j)  = sqrt(o1(i,l,jj-1)*o1(i,l,jj))    ! interpolate fields to interface
          zo2(j) = sqrt(o2(i,l,jj-1)*o2(i,l,jj))
          zn2(j) = sqrt(n2(i,l,jj-1)*n2(i,l,jj))
          zns(j) = sqrt(n4s(i,l,jj-1)*n4s(i,l,jj))
          znd(j) = sqrt(n2d(i,l,jj-1)*n2d(i,l,jj))
          zno(j) = sqrt(no(i,l,jj-1)*no(i,l,jj))
          ztn(j) = (tn(i,l,jj-1)+tn(i,l,jj))/2.
          zti(j) = (ti(i,l,jj-1)+ti(i,l,jj))/2.
          zte(j) = (te(i,l,jj-1)+te(i,l,jj))/2.
          zun(j) = (un(i,l,jj-1)+un(i,l,jj))/2.
          zvn(j) = (vn(i,l,jj-1)+vn(i,l,jj))/2.
          ze(j)  = ne(i,l,jj)
        endif
    enddo 
    z(jmax) = z(jmax-1) + (z(jmax-1)-z(jmax-2))     !patch top level
    ze(jmax) = ze(jmax-1)**2 / + ze(jmax-2)
  endif

  if (nlev == 97) then                              ! high-res TIME-GCM
    do j=1,jmax
          jj=j+28
          z(j)   = zg(i,l,jj)                       ! at interface
          zo(j)  = sqrt(o1(i,l,jj-1)*o1(i,l,jj))    ! interpolate fields to interface
          zo2(j) = sqrt(o2(i,l,jj-1)*o2(i,l,jj))
          zn2(j) = sqrt(n2(i,l,jj-1)*n2(i,l,jj))
          zns(j) = sqrt(n4s(i,l,jj-1)*n4s(i,l,jj))
          znd(j) = sqrt(n2d(i,l,jj-1)*n2d(i,l,jj))
          zno(j) = sqrt(no(i,l,jj-1)*no(i,l,jj))
          ztn(j) = (tn(i,l,jj-1)+tn(i,l,jj))/2.
          zti(j) = (ti(i,l,jj-1)+ti(i,l,jj))/2.
          zte(j) = (te(i,l,jj-1)+te(i,l,jj))/2.
          zun(j) = (un(i,l,jj-1)+un(i,l,jj))/2.
          zvn(j) = (vn(i,l,jj-1)+vn(i,l,jj))/2.
          ze(j)  = ne(i,l,jj)
    enddo 
  endif

end subroutine tzgrid
