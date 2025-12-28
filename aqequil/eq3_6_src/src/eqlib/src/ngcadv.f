      subroutine ngcadv(abar,acflg,acflgo,actwlc,adh,adhh,adhv,
     $ afcnst,al10,aphi,azero,a3bar,a3bars,bacfmx,bdh,bdhh,bdhv,
     $ bdot,bdoth,bdotv,bgamx,bpx,bsigmm,bfje,bfxi,cco2,cgexj,
     $ chfacf,chfsgm,conc,delam,dgpit,dpelm,dpslm,dselm,elam,
     $ eps100,fje,fjeo,fxi,fxio,gpit,ibpxt,ielam,ifcphi1,ifcphi2,
     $ ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,ilcphi1,ilcphi2,ilnnn,
     $ iln2n,ilpsi1,ilpsi2,ilzeta,insgf,iopg,iter,ipndx1,ixrn1,
     $ ixrn2,izmax,jcsort,jern1,jern2,jgext,jsol,kx1,kxt,nalpha,
     $ napt,narn1,narn2,nchlor,ncmpr,net,nhydr,nmut,nmux,nmxi,
     $ nmxx,noutpt,nslt,nslx,nst,nsxi,nsxx,nttyo,omega,palpha,
     $ pelm,pmu,press,pslamn,pslm,qhawep,qpit75,qpracf,q6mode,
     $ rlxgam,selm,sigmam,sigmmo,tempk,ubacmx,ubgamx,uphase,
     $ uspec,wfac,xbar,xbarlg,xbarwc,xbrwlc,zchar,zchcu6,zchsq2)
c
c     This subroutine recalculates the ionic strength, etc., and the
c     activity of water and the molal activity coefficients of aqueous
c     species. If necessary, it also recalculates the mole fraction
c     activity coefficients of solid solution components. The
c     values for the above parameters returned by this subroutine are
c     subject to under-relaxation controls intended to aid convergence,
c     particularly in the case of concentrated aqueous solutions.
c     Associated residual functions such as bsigmm, bfxi, and bgamx
c     are not affected by this under-relaxation.
c
c     This subroutine no longer recalculates the concentrations and
c     other properties of dependent species.
c
c     This subroutine is called by:
c
c       EQLIB/newton.f
c       EQ3NR/arrset.f
c       EQ6/optmzr.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       adh    = Debye-Huckel A(gamma) parameter
c       acflgo = array of old values of the activity coefficients of
c                  the various species
c       aphi   = Debye-Huckel A(phi) parameter
c       azero  = array of hard core diameters of aqueous species
c       bdh    = Debye-Huckel B(gamma) parameter
c       bdot   = B-dot parameter
c       cco2   = coefficients of the Drummond (1981) equation
c       conc   = array of species concentrations
c       fjeo   = old value of the ionic asymmetry (J)
c       fxio   = old value of the ionic strength (I)
c       ibpxt  = array of the numbers of non-zero site-mixing
c                  paramters for computing activity coefficients in
c                  solid solutions
c       ielam  = flag to not use (-1) or use (0) "higher order"
c                  (higher than 2nd order) electrostatic terms in
c                  those activity coefficient models which contain
c                  provision for such terms
c       insgf  = array of activity coefficient flags for aqueous
c                  species, used in the case of neutral solute
c                  species when using the B-dot model
c       iopg   = array of activity coefficient option switches
c       iter = Newton-Raphson iteration number
c       izmax  = max norm of the electrical charges of the aqueous
c                  species
c       jcsort = array of species indices, in order of increasing
c                  concentration, but with sorting restricted to within
c                  phase ranges
c       jsol   = array of identifiers for solid solution activity
c                  coefficient models
c       narn1  = index of the first aqueous species; this is also
c                  the index of the solvent, water
c       narn2  = index of the last aqueous species
c       nchlor = index of the aqueous chloride ion
c       ncmpr  = array giving the range of species belonging to a given
c                  phase: ncmpr(1,np) is the first such species for the
c                  np-th phase, and ncmpr(2,np) is the last such species
c       nhydr  = index of the aqueous hydrogen ion
c       omega  = water constant; ~55.51
c       press  = pressure, bars
c       qpracf = debugging print flag, causes print of activity
c                  coefficient update
c       q6mode = flag denoting usage for EQ3NR or EQ6:
c                  .false. = EQ3NR
c                  .true.  = EQ6NR
c       rlxgam = reported under-relaxation factor applied to corrections
c                  in the activity coefficients of aqueous species
c       sigmmo = old value of sigmam
c       tempk  = temperature, K
c       uphase = array of phase names
c       uspec  = array of species names
c       wfac   = array of non-ideality parameters for solid solutions
c       xbarwc = mole fraction of water (calculated)
c       xbrwlc = log mole fraction of water (calculated)
c       zchar  = array of electrical charge numbers for the various
c                  species
c       zchsq2 = array of (z**2)/2 values
c       zchcu6 = array of (z**3/)6 values
c
c     Principal output:
c
c       abar     average ion size
c       acflg  = array of log activity coefficients of the various
c                  species
c       actwlc = log activity of water (calculated)
c       a3bar    average cube of distance of closest apporach
c       a3bars   characteristic average cube of distance of closest
c                  apporach for each solute species
c       fje    = the ionic asymmetry (the 3rd-order electrostatic
c                  moment function J)
c       fxi    = the ionic strength (the 2nd-order electrostatic
c                  moment function I)
c       sigmam = the sum of solute molalities
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      include 'eqlib/eqldv.h'
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer noutpt,nttyo
c
      integer ibpxt(nxtmax),insgf(natmax),iopg(nopgmx),ipndx1(kmax),
     $ jcsort(nstmax),jern1(jetmax,netmax),jern2(jetmax,netmax),
     $ jgext(netmax),jsol(nxtmax),ncmpr(2,nptmax)
c
      integer nalpha(nsltmx),nmux(3,nmutmx),nmxi(2,natmax),
     $ nmxx(3,nmxmax),nslx(2,nsltmx),nsxi(2,natmax),nsxx(2,nsxmax)
c
      integer napt,nmut,nslt
c
      integer ielam,iter,ixrn1,ixrn2,izmax,kx1,kxt,narn1,narn2,nchlor,
     $ net,nhydr,nst
c
      integer ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,
     $ ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta
c
      logical qhawep,qpit75,qpracf,q6mode
c
      character*48 uspec(nstmax)
      character*48 ubacmx,ubgamx
      character*24 uphase(nptmax)
c
      real*8 dgpit(2,ipbtmx,napmax),dpslm(2,nsltmx),
     $ gpit(ipbtmx,napmax),palpha(ipbtmx,napmax),pmu(nmutmx),
     $ pslamn(0:ipbtmx,nsltmx),pslm(nsltmx)
c
      real*8 delam(2,nazpmx,nazpmx),dpelm(2,nazpmx,nazpmx),
     $ dselm(2,nazmmx:nazpmx),elam(nazpmx,nazpmx),pelm(nazpmx,nazpmx),
     $ selm(nazmmx:nazpmx)
c
      real*8 acflg(nstmax),acflgo(nstmax),azero(natmax),a3bars(natmax),
     $ bpx(ibpxmx,nxtmax),cco2(5),cgexj(jetmax,netmax),conc(nstmax),
     $ wfac(iktmax,nxtmax),xbar(nstmax),xbarlg(nstmax),zchar(nstmax),
     $ zchsq2(nstmax),zchcu6(nstmax)
c
      real*8 adh,adhh,adhv,aphi,bdh,bdhh,bdhv,bdot,bdoth,bdotv
c
      real*8 abar,actwlc,afcnst,al10,a3bar,bacfmx,bgamx,bsigmm,bfje,
     $ bfxi,chfacf,chfsgm,eps100,fje,fjeo,fxi,fxio,omega,press,rlxgam,
     $ sigmam,sigmmo,tempk,xbarwc,xbrwlc
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer itgfix,j2,kcol,np,ns
c
      integer ilnobl
c
      logical qconc1,qconc2,qconc3
c
      real*8 adfx,adfxmx,afjea,ax,bgamxo,chfsmi,dfje,dfx,fjec,fxic,
     $ rlxgcf,rlxgmo,rlxgml,sigmmc,sxl,sxu
c
c-----------------------------------------------------------------------
c
c     Compute the inverse change limit on "Sigma m" (the sum of solute
c     molalities), the ionic strength, and the like.
c
      chfsmi = 1./chfsgm
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Recompute sigma m and compute the associated residual.
c
      call csigm(conc,jcsort,narn1,narn2,nstmax,sigmmc)
      bsigmm = 0.
      if (sigmmo .gt. 0.) bsigmm = (sigmmc - sigmmo)/sigmmo
      sxu = chfsgm*sigmmo
      if (sigmmo .le. 0.) sxu = sigmmc
      sigmam = min(sxu,sigmmc)
      sxl = chfsmi*sigmmo
      sigmam = max(sxl,sigmam)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Recompute the ionic strength (the 2nd-order electrostatic
c     moment function I) and compute the associated residual.
c
      call cfxi(conc,fxic,jcsort,narn1,narn2,nstmax,zchsq2)
      bfxi = 0.
      if (fxio .gt. 0.) bfxi = (fxic - fxio)/fxio
      sxu = chfsgm*fxio
      if (fxio .le. 0.) sxu = fxic
      fxi = min(sxu,fxic)
      sxl = chfsmi*fxio
      fxi = max(sxl,fxi)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Recompute the 3rd-order electrostatic moment function J and
c     the associated residual. Unlike sigma m and I, J can be zero
c     or negative, so the treatment is slightly different.
c
      call cfje(conc,fjec,jcsort,narn1,narn2,nstmax,zchcu6)
      afjea = 0.5*(abs(fjeo) + abs(fjec))
      dfje = fjec - fjeo
      bfje = 0.
      if (afjea .gt. 0.) bfje = dfje/afjea
      if (fjec.gt.0. .and. fjeo.gt.0.) then
        sxu = chfsgm*fjeo
        if (fjeo .le. 0.) sxu = fjec
        fje = min(sxu,fjec)
        sxl = 0.
        fje = max(sxl,fje)
      elseif (fjec.lt.0. .and. fjeo.lt.0.) then
        sxl = chfsgm*fjeo
        if (fjeo .le. 0.) sxl = fjec
        fje = max(sxl,fjec)
        sxu = 0.
        fje = min(sxu,fje)
      elseif (fjec.lt.0. .and. fjeo.gt.0.) then
        fje = 0.
      elseif (fjec.gt.0. .and. fjeo.lt.0.) then
        fje = 0.
      else
        fje = fjec
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the activity coefficients of aqueous species.
c
c     Calling sequence substitutions:
c       acflg for acflgc
c
      call gcoeff(abar,acflg,actwlc,adh,adhh,adhv,al10,
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
c     Calculate the activity coefficients of exchanger species.
c
c     Calling sequence substitutions:
c       acflg for acflgc
c
      call lamgex(acflg,cgexj,jern1,jern2,jetmax,jgext,net,
     $ netmax,nstmax,xbarlg)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute activity coefficient residual norms.
c
      bgamxo = bgamx
c
cXX   The call to EQLIB/betacf.f is below. Reconcile the usage of
cXX   EQLIB/betgam.f and EQLIB/betacf.f when the activity coefficients
cXX   of species in non-aqueous phases are updated numerically the same
cXX   as those of aqueous species.
      call betgam(acflg,acflgo,bgamx,narn1,narn2,nstmax,
     $ ubgamx,uspec)
c
      ubacmx = ubgamx
      bacfmx = bgamx
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Algorithm for under-relaxing changes to the activity coefficients
c     of aqueous species in concentrated solutions. This is designed to
c     damp out oscillation in activity coefficient corrections. Here
c     rlxgam is the under-relaxation factor. Note that activity
c     coefficients may kept fixed for the first few iterations.
c
c     Save the previous value of the under-relaxation factor.
c
      rlxgmo = rlxgam
      rlxgam = 1.0
c
c     Compute a limit (rlxgml)on the under-relaxation factor which
c     is intended to slow its growth from small values, especially
c     in the case where the activity coefficients have been fixed
c     (rlxgam = 0.)
c
      rlxgml = rlxgmo
      if (rlxgml .lt. 0.10) then
        rlxgml = 0.10
      else
        rlxgml = rlxgml + 0.10
      endif
      sxu = 1.0
      rlxgml = min(sxu,rlxgml)
c
c     Determine the course of action according to how concentrated
c     the solution is, based on the sum of solute molalities
c     parameter (Sigma m).
c
      qconc1 = sigmam.ge.2.0 .or. sigmmo.ge.2.0
      qconc2 = sigmam.ge.8.0 .or. sigmmo.ge.8.0
      qconc3 = sigmam.ge.12.0 .or. sigmmo.ge.12.0
c
c     Here itgfix is the number of iterations to hold the activity
c     coefficients constant.
c
      itgfix = 0
      if (qconc1) itgfix = 6
      if (qconc2) itgfix = 8
      if (qconc3) itgfix = 10
c
      if (iter .le. itgfix) then
c
c       Fix the activity coefficients.
c
        rlxgam = 0.
      else
c
c       Apply other under-relaxation controls if called for.
c
        if (bgamxo .gt. 0.) then
c
c         Apply the limit based on the analysis of residual
c         function behavior.
c
          ax = -bgamx/bgamxo
          if (ax .ge. -0.5) then
            rlxgam = rlxgmo*(1./(1. + ax))
            sxl = 0.05
            rlxgam = max(sxl,rlxgam)
          endif
        endif
c
c       Apply the growth from small values limit.
c
        rlxgam = min(rlxgml,rlxgam)
c
c       If the under-relaxation factor is close to one, make it one.
c
        if (abs(rlxgam - 1.0) .le. 0.05) rlxgam = 1.0
c
        if (qconc2) then
c
c         Force some minimum under-relaxation for more concentrated
c         solutions.
c
          if (iter .le. 30) then
            sxu = 0.50
            rlxgam = min(sxu,rlxgam)
          endif
        endif
c
        if (qconc3) then
c
c         Force some minimum under-relaxation for highly concentrated
c         solutions.
c
          if (iter .le. 60) then
            sxu = 0.25
            rlxgam = min(sxu,rlxgam)
          endif
        endif
      endif
c
      if (abs(rlxgam) .gt. eps100) then
c
c       Apply a limit to the change in the activity coefficients.
c       This is yet another form of under-relaxation.
c
        adfx = abs(acflg(narn1) - acflgo(narn1))
        adfxmx = adfx
        do ns = narn1 + 1,narn2
          adfx = abs(acflg(ns) - acflgo(ns))
          adfxmx = max(adfxmx,adfx)
        enddo
        if (adfxmx .gt. chfacf) then
          rlxgcf = chfacf/adfxmx
          rlxgam = min(rlxgcf,rlxgam)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Apply under-relaxation, if any, to the activity coefficients.
c
      if (abs(rlxgam) .le. eps100) then
c
c       Keep activity coefficients fixed.
c
        do ns = narn1,narn2
          acflg(ns) = acflgo(ns)
        enddo
      elseif (abs(rlxgam - 1.0) .gt. eps100) then
c
c       Apply under-relaxation to activity coefficients.
c
        do ns = narn1,narn2
          dfx = acflg(ns) - acflgo(ns)
          acflg(ns) = acflgo(ns) + rlxgam*dfx
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qpracf) then
        write (noutpt,2300) sigmam,fxi,fje,xbrwlc,xbarwc
 2300   format(/3x,'sigmam= ',1pe12.5,/6x,'fxi= ',e12.5,
     $  /6x,'fje= ',e12.5,
     $  /3x,'xbrwlc= ',0pf10.5,/3x,'xbarwc= ',1pe12.5)
c
        write (noutpt,1005)
 1005   format(//11x,'Activity Coefficients of Aqueous Species',
     $  //7x,'Species',23x,'New',11x,'Old',/)
        do ns = narn1,narn2
          write (noutpt,1010) ns,uspec(ns),acflg(ns),acflgo(ns)
 1010     format(1x,i3,2x,a24,2x,1pe12.5,2x,e12.5)
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (q6mode) then
c
c       Recalculate the activity coefficients of any solid solution
c       components present in the ES.
c
        do kcol = kx1,kxt
          np = ipndx1(kcol)
c
c         Calling sequence substitutions:
c           acflg for acflgc
c
          call lambda(acflg,afcnst,bpx,ibpxmx,ibpxt,iktmax,ixrn1,
     $    ixrn2,jsol,ncmpr,noutpt,np,nptmax,nstmax,nttyo,nxtmax,wfac,
     $    xbar,xbarlg,uphase,uspec)
c
          if (qpracf) then
            j2 = ilnobl(uphase(np))
            write (noutpt,1020) uphase(np)(1:j2)
 1020       format(/6x,'Activity Coefficients of Species in ',a,
     $      //7x,'Species',23x,'New',11x,'Old',/)
            do ns = narn1,narn2
              write (noutpt,1010) ns,uspec(ns),acflg(ns),acflgo(ns)
            enddo
          endif
        enddo
c
        call betacf(acflg,acflgo,bacfmx,nst,nstmax,ubacmx,uspec)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
