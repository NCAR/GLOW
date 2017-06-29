! Subroutine ETRANS

! This software is part of the glow model.  Use is governed by the open source
! academic research license agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Banks & Nagy 2-stream electron transport code
! Adapted by Stan Solomon, 1986, 1988
! Uses variable altitude and energy grids
! Updated comments and removed artifacts, scs, 2005
! Moved common blocks into use-associated variables defined in cglow.f90, btf, 2015
! Changed input to subroutine impit to argument list, scs, 2015
! Modernized to remove upper case and numbered line statements, scs, 2016
! Refactored for f90, scs, 2016

! Subroutine EXSECT called first time only to calculate electron impact cross sections.

! Definitions:
! Use-associated variables (formerly in common blocks):  See glow.f90, cglow.f90, and exsect.f
! psi    first term of parabolic d.e., = 1
! alpha  second term "; cm-1
! beta   third term  "; cm-2
! gamma  forth term  "; cm-4 s-1 ev-1
! delz   altitude increments; cm
! del2   sum of altitude increment and next higher increment; cm
! dela   average of "
! delp   product of dela and next higher delz
! delm   product of dela and delz
! dels   product of delz and next higer delz
! den    dummy array for transfer of calculated downward flux
! fac    factor for extrapolating production rate, = 0
! prod   sum of photoelectron production and secondary electrons from protons; cm-3 s-1 ev-1
! eprod  energy of "; ev cm-3
! t1     elastic collision term; cm-1
! t2     elastic + inelastic collision term; cm-1
! tsa    total energy loss cross section for each species; cm2
! produp upward cascade + secondary production; cm-3 s-1 ev-1
! prodwn downward "
! phiup  upward flux; cm-2 s-1 ev-1
! phidwn downward "
! tsigne thermal electron collision term; cm-1
! secion total ionization rate; cm-3 s-1
! secp   secondary electron production; cm-3 s-1 ev-1
! r1     ratio term for calculating upward flux; cm-2 s-1 ev-1
! expt2  exponential term for calculating upward flux
! produa collection array for calculating produp; cm-3 s-1 ev-1
! prodda  "                               prodwn
! phiinf downward flux at top of atmos., divided by avmu; cm-2 s-1 ev-1
! potion ionizaition potential for each species; ev
! avmu   cosine of the average pitch angle

! Array dimensions:
! jmax    number of altitude levels
! nbins   number of energetic electron energy bins
! nmaj    number of major species
! nei     number of states produced by electron impact


  subroutine etrans

    use cglow,only: nmaj,nbins,jmax,nei,ierr,jlocal, &
                    dip,ener,del,aglw,eheat,sion,phitop,zz,pespec, &
                    sespec,zte,ze,zmaj,uflx,dflx,tez,efrac           ! formerly /cglow/
    use cglow,only: siga,sigs,pe,sigex,sec,iimaxx,pin                ! formerly /cxsect/
    use cglow,only: ww                                               ! formerly /cxpars/

    implicit none

    integer,save :: ifirst=1
    integer :: ii,ib,ibb,i,n,jj,j,k,jjj4,iv,ll,kk,im,iq
    real :: prod(jmax), eprod(jmax), t1(jmax), t2(jmax), tsa(nmaj), &
            produp(jmax,nbins), prodwn(jmax,nbins), &
            phiup(jmax), phidwn(jmax), tsigne(jmax), taue(jmax), &
            secion(jmax), secp(nmaj,jmax), r1(jmax), expt2(jmax), &
            produa(jmax), prodda(jmax), phiinf(nbins), potion(nmaj), &
            alpha(jmax),beta(jmax),gamma(jmax),psi(jmax),del2(jmax), &
            delp(jmax),delm(jmax),dels(jmax),den(jmax), &
            delz(jmax),dela(jmax)
    real :: sindip,rmusin,phiout,dag,et,eet,fluxj,edep,epe,ephi,aprod,ein,eout,fac
    real,parameter :: avmu=0.5


    potion = (/16.,16.,18./)
    ierr = 0
    fac = 0.
    sindip = sin(dip)
    rmusin = 1. / sindip / avmu
    psi(1)   = 1.
!
! First call only:  calculate cross-sectons:
!
    if (ifirst == 1) then
      call exsect (ener, del)
      ifirst = 0
    endif
!
! Zero variables:
!
    alpha(1) = 0.
    beta(1) = 0.
    gamma(1) = 0.
    phiout = 0.0
    eheat(:) = 0.0
    eprod(:) = 0.0
    secion(:) = 0.0
    sion(:,:) = 0.0
    aglw(:,:,:) = 0.0
    produp(:,:) = 1.0e-20
    prodwn(:,:) = 1.0e-20
!
! Divide downward flux at top of atmos. by average pitch angle cosine:
!
    phiinf(:) = phitop(:) / avmu
!
! Calcualte delta z's:
!
    delz(1) = zz(2)-zz(1)
    do i=2,jmax
      delz(i) = zz(i)-zz(i-1)
    enddo
    do i=1,jmax-1
      del2(i) = delz(i)+delz(i+1)
      dela(i) = del2(i)/2.
      delp(i) = dela(i)*delz(i+1)
      delm(i) = dela(i)*delz(i)
      dels(i) = delz(i)*delz(i+1)
    enddo
    del2(jmax) = del2(jmax-1)
    dela(jmax) = dela(jmax-1)
    delp(jmax) = delp(jmax-1)
    delm(jmax) = delp(jmax-1)
    dels(jmax) = dels(jmax-1)
!
! Top of energy loop:
!
    do j=nbins,1,-1
!
! Calculate production:
!
      do i = 1, jmax
        prod(i) = (pespec(j,i)+sespec(j,i)) * rmusin / del(j)
        eprod(i) = eprod(i) + prod(i) * ener(j) * del(j) / rmusin
      enddo
!
! Total energy loss cross section for each species:
!
      tsa(:) = 0.0
      if (j > 1) then
        do k = 1, j-1
          do i = 1, nmaj
            tsa(i) = tsa(i) + siga(i,k,j) * (del(j-k)/del(j))
          enddo
        enddo
      else
        do i=1,nmaj
          tsa(i) = tsa(i) + siga(i,1,j) + 1.e-18
        enddo
      endif
!
! Thermal electron energy loss:
!
      jjj4 = j - 1
      if (j == 1) jjj4 = 1
      dag = ener(j) - ener(jjj4)
      if (dag <= 0.0) dag = del(1)
!
      do i = 1, jmax
        et = 8.618e-5 * zte(i)
        eet = ener(j) - et
        if (eet <= 0.0) then
          tsigne(i) = 0.0
        else  
          tsigne(i) = ((3.37e-12*ze(i)**0.97)/(ener(j)**0.94)) &
                    * ((eet)/(ener(j) - (0.53*et))) ** 2.36
        endif
        tsigne(i) = tsigne(i) * rmusin / dag
      enddo
!
! Collision terms:
!
      do i = 1, jmax
        t1(i) = 0.0
        t2(i) = 0.0
        do iv = 1, nmaj
          t1(i) = t1(i) + zmaj(iv,i) * sigs(iv,j) * pe(iv,j)
          t2(i) = t2(i) + zmaj(iv,i) * (sigs(iv,j)*pe(iv,j) + tsa(iv))
        enddo
        t1(i) = t1(i) * rmusin
        t2(i) = t2(i) * rmusin + tsigne(i)
      enddo
!
! Bypass next section if local calculation was specified:
!
      if (jlocal /= 1) then
!
! Solve parabolic d.e. by Crank-Nicholson method to find downward flux:
!
        do i = 2, jmax-1
          psi(i) = 1.
          alpha(i) = (t1(i-1) - t1(i+1)) / (del2(i) * t1(i))
          beta(i) = t2(i) * (t1(i+1) - t1(i-1)) / (t1(i) * del2(i)) &
                  - (t2(i+1) - t2(i-1)) / del2(i) - t2(i)**2 + t1(i)**2
          if (prod(i) < 1.e-30) prod(i) = 1.e-30
          if (prodwn(i,j) < 1.e-30) prodwn(i,j) = 1.e-30
          gamma(i) = (prod(i)/2.0) * (-t1(i) - t2(i) - alpha(i) &
                   - (prod(i+1) - prod(i-1))/prod(i)/del2(i)) &
                   + prodwn(i,j) * (-alpha(i) - t2(i) &
                   - (prodwn(i+1,j)-prodwn(i-1,j))/prodwn(i,j)/del2(i)) &
                   - produp(i,j) * t1(i)
        enddo
        if (abs(beta(2)) < 1.e-20) then
          beta(2) = 1.e-20
          ierr = 2
        endif
        phidwn(2) = gamma(2) / beta(2)
        den(1) = phidwn(2)
        fluxj = phiinf(j)
        call impit(jmax,fluxj,fac,alpha,beta,gamma,psi,del2,delp,delm,dels,den)
        phidwn(:) = den(:)
!
! Apply lower boundary condition: phiup=phidwn.  Should be nearly zero.
!
        phiup(1) = phidwn(1)
!
! Integrate back upward to calculate upward flux:
!
        do i = 2, jmax
          r1(i) = (t1(i)*phidwn(i) + (prod(i)+2.*produp(i,j))/2.) / t2(i)
          taue(i) = t2(i)*delz(i)
          if (taue(i) > 60.) taue(i)=60.
          expt2(i) = exp(-taue(i))
        enddo
        do i=2,jmax
          phiup(i) = r1(i) + (phiup(i-1)-r1(i)) * expt2(i)
        enddo

      else
!
! Local calculation only:
!
        do i = 1, jmax
          if (t2(i) <= t1(i)) then
            ierr = 1
            t2(i) = t1(i) * 1.0001
          endif
          phiup(i) = (prod(i)/2.0 + produp(i,j)) / (t2(i) - t1(i))
          phidwn(i) = (prod(i)/2.0 + prodwn(i,j)) / (t2(i) - t1(i))
        enddo

      endif
!
! Multiply fluxes by average pitch angle cosine and put in arrays:
!
      do i=1,jmax
        uflx(j,i) = phiup(i) * avmu
        dflx(j,i) = phidwn(i) * avmu
      enddo
!
!  Calculate outgoing electron energy flux for conservation check:
!
      phiout = phiout + phiup(jmax) * del(j) * ener(j)
!
! Cascade production:
!
      if (j > 1) then
        do k = 1, j-1
          ll = j - k
          produa(:)=0.
          prodda(:)=0.
          do n = 1, nmaj
            do i=1,jmax
              produa(i) = produa(i) &
                         + zmaj(n,i) * (siga(n,k,j)*pin(n,j)*phidwn(i) &
                         + (1. - pin(n,j))*siga(n,k,j)*phiup(i))
              prodda(i) = prodda(i) &
                         + zmaj(n,i) * (siga(n,k,j)*pin(n,j)*phiup(i) &
                         + (1. - pin(n,j))*siga(n,k,j)*phidwn(i))
            enddo
          enddo
          do i=1,jmax
            produp(i,ll) = produp(i,ll) + produa(i) * rmusin
            prodwn(i,ll) = prodwn(i,ll) + prodda(i) * rmusin
          enddo
        enddo
      endif
      kk = j - 1
      if (kk > 0) then
        do i = 1, jmax
          produp(i,kk) = produp(i,kk)+tsigne(i)*phiup(i)*(del(j)/del(kk))
          prodwn(i,kk) = prodwn(i,kk)+tsigne(i)*phidwn(i)*(del(j)/del(kk))
        enddo
      endif
!
! Electron heating rate:
!
      dag = del(j)
      do i = 1, jmax
        eheat(i) = eheat(i) + tsigne(i) * (phiup(i)+phidwn(i)) * dag**2
      enddo
!
! Electron impact excitation rates:
!
      do ii = 1, jmax
        do i = 1, nmaj
          do ibb = 1, nei
            aglw(ibb,i,ii) = aglw(ibb,i,ii) + (phiup(ii) + phidwn(ii)) &
                             * sigex(ibb,i,j) * del(j) * zmaj(i,ii)
          enddo
        enddo
      enddo
!
! Calculate production of secondaries into k bin for energy j bin and add to production:
!
      do k = 1, iimaxx(j)
        do n = 1, nmaj
          do i = 1, jmax
            secp(n,i) = sec(n,k,j) * zmaj(n,i) * (phiup(i) + phidwn(i))
            sion(n,i) = sion(n,i) + secp(n,i) * del(k)
            secion(i) = secion(i) + secp(n,i) * del(k)
            produp(i,k) = produp(i,k) + (secp(n,i)*.5*rmusin)
            prodwn(i,k) = prodwn(i,k) + (secp(n,i)*.5*rmusin)
          enddo
        enddo
      enddo

    enddo       ! bottom of energy loop

      eheat(:) = eheat(:) / rmusin
!
! Calculate energy deposited as a function of altitude and total energy deposition:
!
    edep = 0.
    do im=1,jmax
      tez(im) = eheat(im)
      do ii=1,nmaj
        tez(im) = tez(im) + sion(ii,im)*potion(ii)
        do iq=1,nei
          tez(im) = tez(im) + aglw(iq,ii,im)*ww(iq,ii)
        enddo
      enddo
      edep = edep + tez(im) * dela(im)
    enddo
!
! Calculate energy input, output, and fractional conservation:
!
    epe = 0.0
    ephi = 0.0
    do i = 2, jmax
      aprod = sqrt(eprod(i)*eprod(i-1))
      epe = epe + aprod * delz(i)
    enddo
    do jj = 1, nbins
      ephi = ephi + phiinf(jj) * ener(jj) * del(jj) / rmusin
    enddo
    ein = ephi + epe
    phiout = phiout / rmusin
    eout = edep + phiout
    efrac = (eout - ein) / ein

    return

  end subroutine etrans

!-----------------------------------------------------------------------

! Subroutine impit solves parabolic differential equation by implicit Crank-Nicholson method

  subroutine impit(jmax,fluxj,fac,alpha,beta,gamma,psi,del2,delp,delm,dels,den)

    implicit none

    integer,intent(in) :: jmax
    real,intent(in) :: fluxj,fac,alpha(jmax),beta(jmax),gamma(jmax), &
               psi(jmax),del2(jmax),delp(jmax),delm(jmax),dels(jmax)
    real,intent(out) :: den(jmax)

    integer :: i1,i,kk,jk
    real :: dem,k(jmax),l(jmax),a(jmax),b(jmax),c(jmax),d(jmax)

    i1 = jmax - 1
    do i = 1, i1
      a(i) = psi(i) / delp(i) + alpha(i) / del2(i)
      b(i) = -2. * psi(i) / dels(i) + beta(i)
      c(i) = psi(i) / delm(i) - alpha(i) / del2(i)
      d(i) = gamma(i)
    enddo
    k(2) = (d(2) - c(2)*den(1)) / b(2)
    l(2) = a(2) / b(2)
    do i = 3, i1
      dem = b(i) - c(i) * l(i-1)
      k(i) = (d(i) - c(i)*k(i-1)) / dem
      l(i) = a(i) / dem
    enddo
    den(i1) = (k(i1) - l(i1)*fluxj) / (1. + l(i1)*fac)
    den(jmax) = den(i1)
    do kk = 1, jmax-3
      jk = i1 - kk
      den(jk) = k(jk) - l(jk) * den(jk + 1)
    enddo

    return

  end subroutine impit
