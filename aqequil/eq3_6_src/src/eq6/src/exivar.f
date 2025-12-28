      subroutine exivar(abar,acflg,acflgo,act,actlg,actwlc,adh,
     $ adhh,adhv,al10,aphi,azero,a3bar,a3bars,bdh,bdhh,bdhv,bdot,
     $ bdoth,bdotv,cco2,cdrs,cegexs,cgexj,conc,conclg,cpgexs,
     $ egexjc,egexjf,egexs,eps100,fje,fjeo,fo2,fo2lg,fsort,fugac,
     $ fugalg,fxi,fxio,ielam,iern1,iern2,ietmax,ifcphi1,ifcphi2,
     $ ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,
     $ iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,
     $ ilrn2,ilzeta,imrn1,imrn2,insgf,iodb,iopg,ipbtmx,istack,
     $ ixrn1,ixrn2,izmax,jcsort,jern1,jern2,jetmax,jflag,jgext,
     $ jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jssort,jstack,
     $ kbt,kdim,kelect,kmax,km1,ko2gaq,kwater,kxt,loph,losp,
     $ lsort,mgext,moph,mosp,mrgexs,mtb,napmax,narn1,narn2,natmax,
     $ nazmmx,nazpmx,nbasp,nbaspd,nbt,nbtmax,nchlor,ncmpr,ndrs,
     $ ndrsmx,ndrsr,nelect,nern1,nern2,net,netmax,ngexsa,ngext,
     $ ngrn1,ngrn2,ngt,ngtmax,nhydr,nmutmx,nmxmax,nodbmx,nopgmx,
     $ noutpt,no2gaq,nphasx,npt,nptmax,nsltmx,nst,nstmax,nsxmax,
     $ nttyo,omega,omeglg,press,qhawep,qpit75,q6mode,sigmam,sigmmo,
     $ tempk,ugexj,ugexmo,uphase,uspec,xbar,xbarlg,xbarw,xbarwc,
     $ xbrwlc,xbrwlg,xlks,zchar,zchcu6,zchsq2,zgexj,zvclg1,zvec1)
c
c     This subroutine expands the system description from the data
c     read from the input file. This includes estimating the numbers
c     of moles of all phases and species present, the concentrations,
c     activity coefficients, and activities of all the species, the
c     ionic strength, and the sum of the molalities of all aqueous
c     aqueous solute species.
c
c     Only a minimum expansion is carried out by the present
c     subroutine. Further refinements, if any, are made by EQ6/optmzr.
c     The iodb(3) debug print option switch is used in this subroutine
c     to control certain prints that are analogous to those in
c     EQ6/optmzr.
c
c     This subroutine is called by:
c
c       EQ6/path.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       adh    = Debye-Huckel A(gamma) parameter
c       aphi   = Debye-Huckel A(phi) parameter
c       azero  = array of hard core diameters of aqueous species
c       bdh    = Debye-Huckel B(gamma) parameter
c       bdot   = B-dot parameter
c       cco2   = coefficients of the Drummond (1981) equation
c       cdrs   = array of reaction coefficients
c       ielam  = flag to not use (-1) or use (0) higher order
c                  electrostatic terms in relevant activity coefficient
c                  models
c       insgf  = array of activity coefficient flags for aqueous
c                  species, used in the case of neutral solute
c                  species when using the B-dot model
c       iopg   = array of activity coefficient option switches
c       izmax  = max norm of the electrical charges of the aqueous
c                  species
c       nbt    = the number of basis species
c       nchlor = index of the aqueous chloride ion
c       ncmpr  = array giving the range in arrays corresponding to
c                  species of those species which belong to a given
c                  phase
c       ndrs   = array parallel to cdrs giving the index of the
c                  corresponding species
c       ndrsr  = array giving the range in the cdrs/ndrs arrays
c                  corresonding to the reaction associated with a
c                  given species
c       nhydr  = index of the aqueous hydrogen ion
c       press  = pressure, bars
c       tempk  = temperature, K
c       uphase = array of phase names
c       uspec  = array of species names
c       xlks   = array of equilibrium constants
c       zchar  = array of electrical charge numbers for the various
c                  species
c       zchsq2 = array of (z**2)/2 values
c       zchcu6 = array of (z**3/)6 values
c       zvec1  = array of master variables
c       zvclg1 = array of logarithms of master variables
c
c     Principal output:
c
c       abar   = average ion size
c       acflg  = array of logarithms of activity coefficients
c       acflgo = array of old values of the activity coefficients of
c                  the various species
c       act    = array of species activities
c       actlg  = array of logarithms of species activities
c       actwlc = log activity of water (calculated)
c       a3bar    average cube of distance of closest apporach
c       a3bars   characteristic average cube of distance of closest
c                  apporach for each solute species
c       conc   = array of species concentrations
c       conclg = array of logarithms of species concentrations
c       fo2    = the oxygen fugacity
c       fo2lg  = logarithm of the oxygen fugacity
c       fugac  = array of gas fugacities
c       fugalg = array of logarithms of gas fugacities
c       fxi    = the ionic strength
c       fxio   = old value of the ionic strength
c       moph   = array of numbers of moles of phases
c       loph   = array of logarithms of numbers of moles of phases
c       mosp   = array of numbers of moles of species
c       losp   = array of logarithms of numbers of moles of species
c       sigmam = the sum of solute molalities
c       sigmmo = old value of sigmam
c       xbar   = array of mole fractions of species
c       xbarlg = array of logarithms of mole fractions of species
c       xbarw  = mole fraction of water
c       xbarwc = mole fraction of water (calculated)
c       xbrwlc = log mole fraction of water (calculated)
c       xbrwlg = log mole fraction of water
c
c-----------------------------------------------------------------------
c
c     Modules.
c
c     The module mod6pt contains data required to evaluate Pitzer's
c     equations. Only a subset of these data is required here.
c
      use mod6pt
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ietmax,ipbtmx,jetmax,kmax,napmax,natmax,nazmmx,nazpmx,
     $ nbtmax,ndrsmx,netmax,ngtmax,nmutmx,nmxmax,nodbmx,nopgmx,nptmax,
     $ nsltmx,nstmax,nsxmax
c
      integer noutpt,nttyo
c
      integer igstak(ngtmax),iindx1(kmax),insgf(natmax),iodb(nodbmx),
     $ iopg(nopgmx),istack(nstmax),jcsort(nstmax),jern1(jetmax,netmax),
     $ jern2(jetmax,netmax),jflag(nstmax),jgext(netmax),
     $ jgsort(ngtmax),jgstak(ngtmax),jjsort(nstmax),jpflag(nptmax),
     $ jsflag(nstmax),jsitex(nstmax),jssort(nstmax),jstack(nstmax),
     $ nbasp(nbtmax),nbaspd(nbtmax),ncmpr(2,nptmax),ndrs(ndrsmx),
     $ ndrsr(2,nstmax),ngexsa(ietmax,jetmax,netmax),
     $ ngext(jetmax,netmax),nphasx(nstmax)
c
      integer ielam,iern1,iern2,ifrn1,ifrn2,igas,ilrn1,ilrn2,imrn1,
     $ imrn2,ixrn1,ixrn2,izmax,kbt,kdim,kelect,km1,ko2gaq,kwater,kxt,
     $ narn1,narn2,nbt,nchlor,nelect,nern1,nern2,net,ngrn1,ngrn2,ngt,
     $ nhydr,no2gaq,npt,nst
c
      integer ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,
     $ ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta
c
      logical qhawep,qpit75,q6mode
c
      character*48 uspec(nstmax)
      character*24 ugexmo(netmax),uphase(nptmax)
      character*8 ugexj(jetmax,netmax)
c
      real*8 acflg(nstmax),acflgo(nstmax),act(nstmax),actlg(nstmax),
     $ azero(natmax),a3bars(natmax),cco2(5),cdrs(ndrsmx),
     $ cegexs(ietmax,jetmax,netmax),cgexj(jetmax,netmax),conc(nstmax),
     $ conclg(nstmax),cpgexs(ietmax,jetmax,netmax),
     $ egexjc(jetmax,netmax),egexjf(jetmax,netmax),
     $ egexs(ietmax,jetmax,netmax),fsort(ngtmax),fugac(ngtmax),
     $ fugalg(ngtmax),loph(nptmax),losp(nstmax),lsort(nstmax),
     $ mgext(jetmax,netmax),moph(nptmax),mosp(nstmax),
     $ mrgexs(ietmax,jetmax,netmax),mtb(nbtmax),xbar(nstmax),
     $ xbarlg(nstmax),xlks(nstmax),zchar(nstmax),zchcu6(nstmax),
     $ zchsq2(nstmax),zgexj(jetmax,netmax),zvclg1(kmax),zvec1(kmax)
c
      real*8 adh,adhh,adhv,aphi,bdh,bdhh,bdhv,bdot,bdoth,bdotv
c
      real*8 abar,actwlc,al10,a3bar,eps100,fje,fjeo,fo2,fo2lg,fxi,
     $ fxio,omega,omeglg,press,sigmam,sigmmo,tempk,xbarw,xbarwc,
     $ xbrwlc,xbrwlg
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jlen,kcol,nb,ns
c
      logical qxbarw
c
      character*56 uspn56
c
      real*8 zx1,zx2
c
      real*8 texp,tlg
c
c-----------------------------------------------------------------------
c
      data qxbarw/.false./
c
c-----------------------------------------------------------------------
c
      if (iodb(3) .ge. 1) then
        write (noutpt,1000)
 1000   format(/11x,' --- Starting initializing calculations ---',/)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize log mole fraction of water in aqueous solution,
c     using a plausible value. This will later be corrected as needed.
c
      xbarwc = 1.0
      xbrwlc = 0.
      xbar(narn1) = xbarwc
      xbarlg(narn1) = xbrwlc
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the concentrations, etc., of basis and non-basis
c     species. Here the activity coefficients are all zero.
c
      call ncmpex(acflg,act,actlg,cdrs,cegexs,cgexj,conc,
     $ conclg,cpgexs,egexjc,egexjf,egexs,eps100,fo2,fo2lg,fsort,
     $ fugac,fugalg,iern1,iern2,ietmax,ifrn1,ifrn2,igas,igstak,
     $ iindx1,ilrn1,ilrn2,imrn1,imrn2,istack,ixrn1,ixrn2,jcsort,
     $ jern1,jern2,jetmax,jflag,jgext,jgsort,jgstak,jjsort,jpflag,
     $ jsflag,jsitex,jssort,jstack,kbt,kdim,kelect,kmax,km1,ko2gaq,
     $ kwater,kxt,loph,losp,lsort,mgext,mrgexs,mtb,moph,mosp,narn1,
     $ narn2,nbasp,nbt,nbtmax,ncmpr,ndrs,ndrsmx,ndrsr,nelect,nern1,
     $ nern2,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,noutpt,
     $ no2gaq,nphasx,npt,nptmax,nst,nstmax,nttyo,omega,omeglg,
     $ press,qxbarw,q6mode,ugexj,ugexmo,uphase,uspec,xbar,xbarlg,
     $ xbarw,xbarwc,xbrwlc,xbrwlg,xlks,zchar,zgexj,zvclg1,zvec1)
c
c     Clamp estimates of concentrations and numbers of moles of
c     dependent aqueous species.
c
      do ns = narn1,narn2
        if (jflag(ns) .ge. 30) then
          conc(ns) = 0.
          mosp(ns) = 0.
          conclg(ns) = -99999.
          losp(ns) = -99999.
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make the first estimates of SIGMA m, ionic strength, and J taking
c     into account the concentrations of only data file basis species.
c     This generally provides a set of safe values, as the computed
c     concentrations of non-basis species at this point could be
c     extremely large.
c
      sigmam = 0.
      fxi = 0.
      fje = 0.
      do nb = 1,nbt
        ns = nbaspd(nb)
        if (ns.gt.narn1 .and. ns.le.narn2) then
          sigmam = sigmam + conc(ns)
          fxi = fxi + conc(ns)*zchsq2(ns)
          fje = fje + conc(ns)*zchcu6(ns)
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make the first real estimate of the mole fraction of water in
c     aqueous solution.
c
      xbarwc = omega/(omega + sigmam)
      xbrwlc = tlg(xbarwc)
      xbarw = xbarwc
      xbrwlg = xbrwlc
      xbar(narn1) = xbarwc
      xbarlg(narn1) = xbrwlc
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make the first estimates of the aqueous species activity
c     coefficients.
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
c     Make the first estimates of the exchanger species activity
c     coefficients.
c
c     Calling sequence substitutions:
c       acflg for acflgc
c
      call lamgex(acflg,cgexj,jern1,jern2,jetmax,jgext,net,
     $ netmax,nstmax,xbarlg)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(3) .ge. 1) then
c
c       Print the attempted phase assemblage.
c
        write (noutpt,1010)
 1010   format(/' Initial phase assemblage:',/)
        do kcol = 1,kbt
          nb = iindx1(kcol)
          ns = nbaspd(nb)
c
c         Calling sequence substitutions:
c           uspec(ns) for unam48
c
          call fmspnx(jlen,uspec(ns),uspn56)
          write (noutpt,1020) kcol,uspn56(1:jlen)
 1020     format(2x,i3,2x,a)
        enddo
        do kcol = km1,kxt
          ns = iindx1(kcol)
c
c         Calling sequence substitutions:
c           uspec(ns) for unam48
c
          call fmspnm(jlen,uspec(ns),uspn56)
          write (noutpt,1020) kcol,uspn56(1:jlen)
        enddo
        write (noutpt,1030)
 1030   format(1x)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(3) .ge. 2) then
        write (noutpt,1050)
 1050   format(/16x,'--- Initialization Summary ---',
     $  //2x,'kcol   Name',32x,'zvclg1      zvec1',/)
        do kcol = 1,kbt
          nb = iindx1(kcol)
          ns = nbasp(nb)
          zx1 = zvclg1(kcol)
          zx2 = texp(zx1)
c
c         Calling sequence substitutions:
c           uspec(ns) for unam48
c
          call fmspnx(jlen,uspec(ns),uspn56)
          write (noutpt,1060) kcol,uspn56,zx1,zx2
 1060     format(1x,i4,2x,a32,2x,f10.4,2x,1pe12.5)
        enddo
        do kcol = km1,kxt
          ns = iindx1(kcol)
          zx1 = zvclg1(kcol)
          zx2 = texp(zx1)
c
c         Calling sequence substitutions:
c           uspec(ns) for unam48
c
          call fmspnm(jlen,uspec(ns),uspn56)
          write (noutpt,1060) kcol,uspn56,zx1,zx2
        enddo
        write (noutpt,1030)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(3) .ge. 1) then
        write (noutpt,1120) sigmam,fxi,fje
 1120   format(/8x,'sigmam= ',1pe12.5,/11x,'fxi= ',e12.5,
     $  /11x,'fje= ',e12.5,/)
      endif
c
      if (iodb(3) .ge. 4) then
        write (noutpt,1130)
 1130   format(/7x,'Species',20x,'gamma',/)
        do ns = narn1,narn2
          write (noutpt,1140) uspec(ns),acflg(ns)
 1140     format(5x,a24,3x,1pe12.5)
        enddo
        write (noutpt,1030)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      sigmmo = sigmam
      fxio = fxi
      fjeo = fje
      do ns = narn1,narn2
        acflgo(ns) = acflg(ns)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(3) .ge. 1) then
        write (noutpt,1200)
 1200   format(/11x,' --- Finished initializing calculations ---',/)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
