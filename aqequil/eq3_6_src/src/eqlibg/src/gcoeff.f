      subroutine gcoeff(abar,acflgc,actwlc,adh,adhh,adhv,al10,
     $ aphi,azero,a3bar,a3bars,bdh,bdhh,bdhv,bdot,bdoth,bdotv,
     $ cco2,conc,delam,dgpit,dpelm,dpslm,dselm,elam,fje,fxi,gpit,
     $ ielam,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,
     $ ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta,insgf,
     $ iopg,ipbtmx,izmax,jcsort,nalpha,napmax,napt,narn1,narn2,
     $ natmax,nazmmx,nazpmx,nchlor,nhydr,nmut,nmutmx,nmux,nmxi,
     $ nmxmax,nmxx,nopgmx,noutpt,nslt,nsltmx,nslx,nstmax,nsxi,
     $ nsxmax,nsxx,nttyo,omega,palpha,pelm,pmu,press,pslamn,
     $ pslm,qhawep,qpit75,selm,sigmam,tempk,uspec,xbarwc,xbrwlc,
     $ zchar,zchsq2,zchcu6)
c
c     This subroutine computes activity coefficients of aqueous species
c     using various models. The model used is determined by the option
c     flag iopg(1):
c
c        -1  Davies equation
c         0  B-dot equation plus others
c         1  Pitzer's equations. The use or non-use of higher order
c            electrostatic terms is determined by the E-theta flag
c            on the supporting data file.
c         2  HC + DHC
c
c     The computed activity coefficients may be normalized to make
c     them consistent with a given pH convention. This is determined
c     by the option flag iopg(2):
c
c        -1  No normalization
c         0  Normalization to the NBS pH scale
c         1  Normalization to the Mesmer pH scale (pH is numerically
c            equal to -log m(H+))
c
c     This subroutine is called by:
c
c       EQLIB/ngcadv.f
c       EQ3NR/arrset.f
c       EQ6/exivar.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       adh    = Debye-Huckel A(gamma) parameter
c       azero  = array of hard core diameters of aqueous species
c       bdh    = Debye-Huckel B parameter
c       bdot   = B-dot parameter
c       cco2   = coefficients of the Drummond (1981) equation
c       conc   = array of species concentrations
c       fje    = the 3rd-order electrostatic moment function J
c       fxi    = the ionic strength (the 2nd-order electrostatic
c                  moment function I)
c       ielam  = flag to not use (-1) or use (0) "higher order"
c                  (higher than 2nd order) electrostatic terms in
c                  those activity coefficient models which contain
c                  provision for such terms
c       insgf  = array of activity coefficient flags for aqueous
c                  species, used in the case of neutral solute
c                  species when using the B-dot model
c       iopg   = array of activity coefficient option switches
c       izmax  = max norm of the electrical charges of the aqueous
c                  species
c       jcsort = array of species indices, in order of increasing
c                  concentration, but with sorting restricted to within
c                  phase ranges
c       narn1  = index of the first aqueous species; this is also
c                  the index of the solvent, water
c       narn2  = index of the last aqueous species
c       nchlor = index of the aqueous chloride ion
c       nhydr  = index of the aqueous hydrogen ion
c       omega  = water constant; ~55.51
c       press  = pressure, bars
c       sigmam = sum of solute molalities
c       tempk  = temperature, K
c       xbarwc = mole fraction of water (calculated)
c       xbrwlc = log mole fraction of water (calculated)
c       zchar  = array of electrical charge numbers for the various
c                  species
c       zchsq2 = array of z**2/2 values
c       zchcu6 = array of z**3/6 values
c
c     Principal Output:
c
c       acflgc = array of log activity coefficients (calculated)
c                  of the various species
c       actwlc = log activity of water (calculated)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations, for all variables except
c     those intimately related to Pitzer's equations.
c
      integer ipbtmx,natmax,nopgmx,nstmax
c
      integer insgf(natmax),iopg(nopgmx),jcsort(nstmax)
c
      integer ielam,izmax,narn1,narn2,nchlor,nhydr,noutpt,nttyo
c
      integer ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,
     $ ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta
c
      real(8) acflgc(nstmax),azero(natmax),a3bars(natmax),conc(nstmax),
     $ cco2(5),zchar(nstmax),zchsq2(nstmax),zchcu6(nstmax)
c
      character(len=48) uspec(nstmax)
c
      real(8) adh,adhh,adhv,aphi,bdh,bdhh,bdhv,bdot,bdoth,bdotv
c
      real(8) abar,actwlc,al10,a3bar,fje,fxi,omega,press,sigmam,tempk,
     $ xbarwc,xbrwlc
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations, for all variables
c     intimately related to Pitzer's equations.
c
c     Integer variables and arrays.
c
      integer napmax,nazmmx,nazpmx,nmutmx,nmxmax,nsltmx,nsxmax
c
      integer napt,nmut,nslt
c
      logical qhawep,qpit75
c
      integer nalpha(nsltmx),nmux(3,nmutmx),nmxi(2,natmax),
     $ nmxx(3,nmxmax),nslx(2,nsltmx),nsxi(2,natmax),nsxx(2,nsxmax)
c
c     Empirical interaction coefficients and related functions at
c     temperature and pressure.
c
      real(8) dgpit(2,ipbtmx,napmax),dpslm(2,nsltmx),
     $ gpit(ipbtmx,napmax),palpha(ipbtmx,napmax),pmu(nmutmx),
     $ pslamn(0:ipbtmx,nsltmx),pslm(nsltmx)
c
c     Theoretical higher-order electrostatic interaction coefficients
c     and related functions at temperature and pressure.
c
      real(8) delam(2,nazpmx,nazpmx),dpelm(2,nazpmx,nazpmx),
     $ dselm(2,nazmmx:nazpmx),elam(nazpmx,nazpmx),pelm(nazpmx,nazpmx),
     $ selm(nazmmx:nazpmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer iz,j2,n,na,nref,ns,nss,ns1,ns2,ns3
c
      integer ilnobl
c
      character(len=8) ux8
c
      real(8) acfnbs,afac,alm,bt,btp,chc,delacf,elsum,elsump,elsums,
     $ elsumw,f,fa,flsum,fp,fpp,lnaw,lngam,muterm,musum,musumw,slsum,
     $ slsump,spsum,spsump,ssumw,xx,xy,zchtot
c
      real(8) csum,csumca,csum1,csum2,psum1,psum2,zesum
c
      real(8) tlg
c
c-----------------------------------------------------------------------
c
c     Here btp is the effective product of Debye-Huckel B and the hard
c     core diameter used in Pitzer's equations.
c
      data btp /1.2/
c
c     Here chc is the hard core repulsion constant.
c
      data chc /1.2577e-3/
c
c-----------------------------------------------------------------------
c
c     Note: the following statements don't actually do anything except
c     cause the compiler not to complain that fje, press, and zchcu6
c     are not used.
c
      xx = fje
      xy = xx
      fje = xy
c
      xx = press
      xy = xx
      press = xy
c
      xx = zchcu6(1)
      xy = xx
      zchcu6(1) = xy
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the mole fraction of water from sigmam. Note that
c     this value is subject to under-relaxation truncation constraints
c     placed on sigmam.
c
      xbarwc = omega/(omega + sigmam)
      xbrwlc = tlg(xbarwc)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the total charge per kg H2O (Z in Pitzer's equations).
c     Note that this is a sorted summation. It is not necessary to
c     exclude water, because it has zero charge.
c
      zchtot = 0.
      do nss = narn1,narn2
        ns = jcsort(nss)
        zchtot = zchtot + conc(ns)*abs(zchar(ns))
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopg(1) .eq. -1) then
c
c       Use the Davies equation.
c
        call gdavie(acflgc,actwlc,adh,al10,fxi,narn1,narn2,nstmax,
     $  omega,sigmam,xbrwlc,zchsq2)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopg(1) .eq. 0) then
c
c       Use the B-dot equation and associated approximations.
c
        call gbdot(acflgc,actwlc,adh,al10,azero,bdh,bdot,cco2,
     $  fxi,insgf,narn1,narn2,natmax,nstmax,omega,sigmam,tempk,
     $  xbrwlc,zchar,zchsq2)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopg(1) .eq. 1) then
c
c       Use Pitzer's equations.
c
c       Compute the Debye-Huckel function f and its ionic strength
c       derivatives.
c
        bt = btp
        call gfdho(aphi,bt,f,fp,fpp,fxi)
c
        elsumw = 0.
        elsums = 0.
        elsump = 0.
        if (ielam .ge. 0) then
c
c         Get E-lambda (elam) and its derivatives (delam) for the
c         various charge combinations.
c
          call gelam(aphi,delam,dpelm,elam,fxi,izmax,nazpmx,
     $    pelm,qpit75)
c
c         Compute the following second order sums in E-lambda and its
c         ionic strength derivatives:
c
c           elsumw: SUM(ij) [E-lambda(ij) + I*E-lambda'(ij)]*m(i)*m(j)
c
c           elsums: SUM(ij) E-lambda'(ij)*m(i)*m(j)
c
c           elsump: SUM(ij) E-lambda''(ij)*m(i)*m(j) [Not used]
c
c         Here ij refers to any species pair, but non-zero
c         contributions only arise from cc' and aa' pairs. Note that:
c
c           elsumw = 2*
c                 [SUM(c'>c) {E-theta(cc') + I*E-theta'(cc')}m(c)m(c')
c                + SUM(a'>a) {E-theta(aa') + I*E-theta'(aa')}m(a)m(a')]
c
c           elsums = 2*[SUM(c'>c) E-theta'(cc')m(c)m(c')
c                     + SUM(a'>a) E-theta'(aa')m(a)m(a')]
c
c           elsump = 2*[SUM(c'>c) E-theta''(cc')m(c)m(c')
c                     + SUM(a'>a) E-theta''(aa')m(a)m(a')] [Not used]
c
c         Note also that "flim" in EQ3/6 is 2 * F, as F is commonly
c         defined in the Pitzer literature. This is because
c         flim will be multipled by z(i)^2/2 instead of z(i)^2.
c
          call gesum(conc,delam,elam,elsump,elsums,elsumw,
     $    fxi,narn1,narn2,nazpmx,nstmax,zchar)
c
c         Compute the following first order sum arrays in E-lambda
c         and its ionic strength derivative:
c
c           selm(i):  SUM(j) E-lambda(ij)*m(j)
c           dselm(1,i): SUM(j) E-lambda'(ij)*m(j) [Not used]
c
          call gselm(conc,delam,dselm,elam,izmax,narn1,narn2,nazmmx,
     $    nazpmx,nstmax,selm,zchar)
        endif
c
c       Compute the g(x) function (gpit) and its derivatives (dgpit)
c       for pairs of alpha coefficients.
c
        call gdd(dgpit,fxi,gpit,ipbtmx,napmax,napt,palpha)
c
c       Compute the S-lambda functions and their ionic strength
c       derivatives.
c
        call gslam(dgpit,dpslm,gpit,ipbtmx,nalpha,napmax,nslt,
     $  nsltmx,pslamn,pslm)
c
c       Compute the following second order sum in S-lambda and its
c       ionic strength derivatives:
c
c         ssumw: SUM(ij) [ S-lambda(ij) + I*S-lambda'(ij) ]*m(i)*m(j)
c
        call gssum(conc,dpslm,fxi,nslt,nsltmx,nslx,nstmax,pslm,
     $  ssumw,uspec)
c
        if (qhawep) then
c
c         Use the form corresponding to C (from Cphi), psi, zeta,
c         and any mu not derived from Cphi, psi, or zeta.
c         Calculate the following sum that appears in the activity
c         coefficients of both water and solutes:
c
c           csumca:  SUM(ca) m(c)m(a)C(ca)
c
c         This sum is technically a double sum, but it is evaluated
c         here as an equivalent single sum over all ca combinations.
c
          csumca = 0.
          do n = ifcphi1,ilcphi1
            ns1 = nmux(1,n)
            ns3 = nmux(3,n)
            csumca = csumca + pmu(n)*conc(ns1)*conc(ns3)
          enddo
        endif
c
c       Compute the contributions to the activity coefficient of
c       water from third-order terms.
c
        muterm = 0.
c
        if (qhawep) then
c
c         Use the form corresponding to C (from Cphi), psi, zeta,
c         and any mu not derived from Cphi, psi, or zeta.
c
          csum = zchtot*csumca
c
          psum1 = 0.
c
c         This sum is technically a triple sum, but is evaluated
c         as a single sum over all cc'a combinations.
c
          do n = ifpsi1,ilpsi1
            ns1 = nmux(1,n)
            ns2 = nmux(2,n)
            ns3 = nmux(3,n)
            psum1 = psum1 + pmu(n)*conc(ns1)*conc(ns2)*conc(ns3)
          enddo
c
          psum2 = 0.
c
c         This sum is technically a triple sum, but is evaluated
c         as a single sum over all aa'c combinations.
c
          do n = ifpsi2,ilpsi2
            ns1 = nmux(1,n)
            ns2 = nmux(2,n)
            ns3 = nmux(3,n)
            psum2 = psum2 + pmu(n)*conc(ns1)*conc(ns2)*conc(ns3)
          enddo
c
          zesum = 0.
c
c         This sum is technically a triple sum, but is evaluated
c         as a single sum over all nca combinations.
c
          do n = ifzeta,ilzeta
            ns1 = nmux(1,n)
            ns2 = nmux(2,n)
            ns3 = nmux(3,n)
            zesum = zesum + pmu(n)*conc(ns1)*conc(ns2)*conc(ns3)
          enddo
c
          musum = 0.
c
c         The following involves sums in mu(nnn) and mu(nnn').
c
          do n = ifnnn,ilnnn
            musum = musum + 2.*pmu(n)*conc(ns1)*conc(ns1)*conc(ns1)
          enddo
c
          do n = ifn2n,iln2n
            ns1 = nmux(1,n)
            ns2 = nmux(2,n)
            ns3 = nmux(3,n)
c
c           Note that it doesn't matter here whether ns1 = ns2 or
c           ns2 = ns3.
c
            musum = musum + 6.*pmu(n)*conc(ns1)*conc(ns2)*conc(ns3)
          enddo
c
          muterm = csum + psum1 + psum2 + zesum + musum
c
c         Multiply muterm by 2 for use in calculating the activity
c         of water directly.
c
          muterm = 2.*muterm
c
        else
c
c         Use the original pure mu form.
c
c         Compute the following third order sum in mu:
c
c           musumw: SUM(ijk) mu(ijk)*m(i)*m(j)*m(k)
c
          call gmsum(conc,musumw,nmut,nmutmx,nmux,nstmax,pmu,uspec)
          muterm = 2.*musumw
c
        endif
c
c       Compute the activity of water.
c
        actwlc = -(1./(al10*omega)) *
     $  (sigmam + fxi*fp - f + elsumw + ssumw + muterm)
c
c       Compute the activity coefficient of water. Note that this
c       is mole-fraction based.
c
        acflgc(narn1) = actwlc - xbrwlc
c
c       Compute the following second order sums in the first and
c       second ionic strength derivatives of S-lambda:
c
c         spsum:  SUM(jk) S-lambda'(jk)*m(j)*m(k)
c         spsump: SUM(jk) S-lambda''(jk)*m(j)*m(k) [Not used]
c
c       Here jk always refers to a cation-ion pair. Note that:
c
c         spsum = 2 * SUM(ca) B'(ca)m(a)m(c)
c         spsump = 2 * SUM(ca) B''(ca)m(a)m(c)
c
c       Note also that "flim" in EQ3/6 is 2 * F, as F is commonly
c       defined in the Pitzer literature. This is because
c       flim will be multipled by z(i)^2/2 instead of z(i)^2.
c
        call gsdsm(conc,dpslm,nslt,nsltmx,nslx,nstmax,
     $  spsum,spsump,uspec)
c
        flsum = fp + elsums + spsum
c
c       Compute the activity coefficients of the aqueous
c       solute species.
c
        do ns = narn1 + 1,narn2
          na = ns - narn1 + 1
c
c         Logically, na (or ns) denotes aqueous species i.
c
          acflgc(ns) = 0.
          iz = nint(zchar(ns))
c
          elsum = 0.
          if (ielam .ge. 0) then
c
c           Get the E-lambda sum (this only depends on the
c           magnitude of the charge number):
c
c             selm(i): SUM(j) E-lambda(ij)*m(j)
c
            elsum = selm(iz)
          endif
c
c         Compute the following first order sums in S-lambda
c         and its ionic strength derivative:
c
c           slsum(i):  SUM(j) S-lambda(ij)*m(j)
c           slsump(i): SUM(j) S-lambda'(ij)*m(j)
c
          call gsgsm(conc,dpslm,na,natmax,nsltmx,nstmax,nsxi,nsxmax,
     $    nsxx,pslm,slsum,slsump,uspec)
c
c         Compute the contributions from third-order terms.
c
          muterm = 0.
c
          if (qhawep) then
c
c           Use the form corresponding to C (from Cphi), psi, zeta,
c           and any mu not derived from Cphi, psi, or zeta.
c
            csum1 = 0.
            if (zchar(ns) .gt. 0.) then
c
c             Have a cation. Loop over C(ca).
c
              do n = ifcphi1,ilcphi1
                ns1 = nmux(1,n)
                if (ns .eq. ns1) then
                  ns3 = nmux(3,n)
                  csum1 = csum1 + pmu(n)*conc(ns3)
                endif
              enddo
            elseif (zchar(ns) .lt. 0.) then
c
c             Have an anion. Loop over C(ac).
c
              do n = ifcphi2,ilcphi2
                ns1 = nmux(1,n)
                if (ns .eq. ns1) then
                  ns3 = nmux(3,n)
                  csum1 = csum1 + pmu(n)*conc(ns3)
                endif
              enddo
            endif
            csum1 = zchtot*csum1
c
            psum1 = 0.
            if (zchar(ns) .gt. 0.) then
c
c             Have a cation. Loop over psi(cc'a).
c
              do n = ifpsi1,ilpsi1
                ns1 = nmux(1,n)
                ns2 = nmux(2,n)
                ns3 = nmux(3,n)
                if (ns .eq. ns1) then
                  psum1 = psum1 + pmu(n)*conc(ns2)*conc(ns3)
                elseif (ns .eq. ns2) then
                  psum1 = psum1 + pmu(n)*conc(ns1)*conc(ns3)
                endif
              enddo
            elseif (zchar(ns) .lt. 0.) then
c
c             Have an anion. Loop over psi(aa'c).
c
              do n = ifpsi2,ilpsi2
                ns1 = nmux(1,n)
                ns2 = nmux(2,n)
                ns3 = nmux(3,n)
                if (ns .eq. ns1) then
                  psum1 = psum1 + pmu(n)*conc(ns2)*conc(ns3)
                elseif (ns .eq. ns2) then
                  psum1 = psum1 + pmu(n)*conc(ns1)*conc(ns3)
                endif
              enddo
            endif
c
            psum2 = 0.
c
c           This sum is technically a double sum, but is evaluated
c           as a single sum over all distinct values of psi(aa'c)
c           or psi(cc'a).
c
            if (zchar(ns) .gt. 0.) then
c
c             Have a cation. Loop over psi(aa'c).
c
              do n = ifpsi2,ilpsi2
                ns1 = nmux(1,n)
                ns2 = nmux(2,n)
                ns3 = nmux(3,n)
                if (ns .eq. ns3) then
                  psum2 = psum2 + pmu(n)*conc(ns1)*conc(ns2)
                endif
              enddo
            elseif (zchar(ns) .lt. 0.) then
c
c             Have an anion. Loop over psi(cc'a).
c
              do n = ifpsi1,ilpsi1
                ns1 = nmux(1,n)
                ns2 = nmux(2,n)
                ns3 = nmux(3,n)
                if (ns .eq. ns3) then
                  psum2 = psum2 + pmu(n)*conc(ns1)*conc(ns2)
                endif
              enddo
            endif
c
            csum2 = abs(zchar(ns))*csumca
c
            zesum = 0.
c
c           This sum is technically a double sum, but is evaluated
c           as a single sum over all distinct values of zeta(nca).
c
c           Have a cation, an anion, or a neutral. Loop over zeta(nca).
c
            do n = ifzeta,ilzeta
              ns1 = nmux(1,n)
              ns2 = nmux(2,n)
              ns3 = nmux(3,n)
              if (ns .eq. ns1) then
                zesum = zesum + pmu(n)*conc(ns2)*conc(ns3)
              elseif (ns .eq. ns2) then
                zesum = zesum + pmu(n)*conc(ns1)*conc(ns3)
              elseif (ns .eq. ns3) then
                zesum = zesum + pmu(n)*conc(ns1)*conc(ns2)
              endif
            enddo
c
            musum = 0.
c
            if (iz .eq. 0) then
c
c             The following involves sums in mu(nnn) and mu(nnn').
c
              do n = ifnnn,ilnnn
                ns1 = nmux(1,n)
                if (ns .eq. ns1) then
                  musum = musum + 3.*pmu(n)*conc(ns1)*conc(ns1)
                endif
              enddo
c
              do n = ifn2n,iln2n
                ns1 = nmux(1,n)
                ns2 = nmux(2,n)
                ns3 = nmux(3,n)
                if (ns .eq. ns1) then
                  if (ns1 .eq. ns2) then
                    musum = musum + 6.*pmu(n)*conc(ns1)*conc(ns3)
                  else
                    musum = musum + 3.*pmu(n)*conc(ns3)*conc(ns3)
                  endif
                elseif (ns .eq. ns3) then
                  if (ns2 .eq. ns3) then
                    musum = musum + 6.*pmu(n)*conc(ns1)*conc(ns3)
                  else
                    musum = musum + 3.*pmu(n)*conc(ns1)*conc(ns1)
                  endif
                endif
              enddo
            endif
c
            muterm = csum1 + psum1 + psum2 + csum2 + zesum + musum
c
          else
c
c           Use the original pure mu form.
c
c           Compute the following second order sum in mu:
c
c             musum(i): SUM(jk) mu(ijk)*m(j)*m(k)
c
            call gmdsm(conc,musum,na,natmax,nmutmx,nmxi,nmxmax,
     $      nmxx,ns,nstmax,pmu,uspec)
c
            muterm = 3.*musum
          endif
c
c         Note that zchsq2(ns) is z(i)^2/2.
c
          acflgc(ns) = (zchsq2(ns)*flsum + 2.0*elsum + 2.0*slsum
     $    + muterm) / al10
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopg(1) .eq. 2) then
c
c       Use the HC + DHC equations
c
c       Get the average ion size (abar), the characteristic average
c       cube of the distance of closest approach (a3bars), and the
c       average cube of the distance of closest approach (a3bar).
c
        call cabar(abar,azero,conc,jcsort,fxi,narn1,narn2,natmax,
     $  nstmax,zchsq2)
c
        call ca3bar(azero,a3bar,a3bars,conc,jcsort,narn1,
     $  narn2,natmax,nstmax,sigmam)
c
c       Get the Debye-Huckel function f and its ionic strength
c       derivatives.
c
        bt = bdh*abar
        call gfdhc(adh,bt,f,fp,fpp,fxi)
c
c       Compute the activity of water.
c
        xx = 1. + (sigmam/omega)
        alm = log(xx)
        lnaw = - alm - (fxi*fp - f)/omega
     $  - (chc*a3bar*(sigmam**2))/(omega)
        actwlc = lnaw/al10
        acflgc(narn1) = actwlc - xbrwlc
c
c       Compute the activity coefficients of the aqueous
c       solute species.
c
        fa = 3.*fxi*fp - 2.*f
        do ns = narn1 + 1,narn2
          acflgc(ns) = 0.
          na = ns - narn1 + 1
          afac = (azero(na) - abar)/abar
          lngam = -alm + zchsq2(ns)*fp + (zchsq2(ns)*fa*afac)/fxi
     $    + 2.*chc*sigmam**a3bars(na)
          acflgc(ns) = lngam/al10
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopg(1).le.-2 .or. iopg(1).ge.3) then
        write (ux8,'(i5)') iopg(1)
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1000) ux8(1:j2)
        write (nttyo,1000) ux8(1:j2)
 1000   format(/' * Error - (EQLIBG/gcoeff) Programming error trap:',
     $  ' Have an',/7x,'illegal value of ',a,' for the option switch',
     $  ' iopg(1), which determines',/7x,'the model to use for the',
     $  ' activity coefficients of aqueous species.')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make activity coefficients consistent with a specified pH
c     scale.
c
c       iopg(2) = -1   Use activity coefficients as is
c                = 0   Make consistent with NBS pH scale
c                = 1   Make consistent with the Mesmer pH scale
c                      (log gamma H+ = 0)
c
      if (iopg(2) .eq. 0)  then
        call nbsgam(acfnbs,adh,fxi,nchlor,noutpt,nttyo)
        delacf = acfnbs - acflgc(nchlor)
        nref = nchlor
        call gcscal(acflgc,delacf,narn1,narn2,nref,nstmax,zchar)
      endif
c
      if (iopg(2) .eq. 1) then
        delacf = -acflgc(nhydr)
        nref = nhydr
        call gcscal(acflgc,delacf,narn1,narn2,nref,nstmax,zchar)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
