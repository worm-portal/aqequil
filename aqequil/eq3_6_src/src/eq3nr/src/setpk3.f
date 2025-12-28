      subroutine setpk3(electr,iindx1,jflag,jflgi,kbt,kdim,kmax,
     $ kmt,kprs,kwater,kxt,mtb,mtbi,mtbaq,mtbaqi,narn1,narn2,nbasp,
     $ nbaspd,nbaspi,nbti,nbtmax,ndrsrd,nern1,nern2,nobswt,nstmax,
     $ ntitl,ntitl2,ntitmx,omeglg,press,pressi,scamas,sigzi,tempc,
     $ tempci,ubmtbi,uobsw,uspec,utitl,utitl2,uzveci,uzvec1,
     $ zvclgi,zvclg1)
c
c     This subroutine sets up certain variables and arrays for writing
c     on the pickup file. This subroutine must be called prior to
c     writing a pickup file.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c     Principal output:
c
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer kmax,nbtmax,nstmax,ntitmx
c
      integer iindx1(kmax),jflag(nstmax),jflgi(nbtmax),nbasp(nbtmax),
     $ nbaspd(nbtmax),nbaspi(nbtmax),ndrsrd(2,nstmax)
c
      integer kbt,kdim,kmt,kprs,kwater,kxt,narn1,narn2,nern1,nern2,
     $ nbti,nobswt,ntitl,ntitl2
c
      character(len=80) utitl(ntitmx),utitl2(ntitmx)
      character(len=48) ubmtbi(nbtmax),uobsw(2,nbtmax),uspec(nstmax),
     $ uzveci(kmax),uzvec1(kmax)
c
      real(8) mtb(nbtmax),mtbi(nbtmax),mtbaq(nbtmax),mtbaqi(nbtmax),
     $ zvclgi(kmax),zvclg1(kmax)
c
      real(8) electr,omeglg,press,pressi,scamas,sigzi,tempc,tempci
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer kcol,krow,n,nb,nbi,nerr,nsi,ns1,ns2,nt1
c
      real(8) lscama
c
      real(8) tlg
c
c-----------------------------------------------------------------------
c
      nerr = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Old title.
c
      ntitl2 = ntitl
      do n = 1,ntitl
        utitl2(n) = utitl(n)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Original temperature.
c
      tempci = tempc
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Original pressure.
c
      pressi = press
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Index limits.
c
      kprs = 0
c
      kmt = kdim
      kxt = kdim
      nbti = kbt
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Species for which mass balances are defined.
c
c       ubmtbi = names of the data file basis species for which mass
c                  balances are defined
c       jflgi  = jflag input for basis species
c          0 = Retain as an active basis species
c         30 = Convert to a dependent species; fold the mass balance
c                total for this species into the mass balance totals
c                of basis species which remain active
c
      do krow = 1,kbt
        nb = iindx1(krow)
        ns1 = nbaspd(nb)
        ns2 = nbasp(nb)
        nbi = krow
        nsi = nbaspi(nb)
        ubmtbi(nbi) = uspec(nsi)
c
        if (ns1.ge.narn1 .and. ns1.le.narn2) then
          mtbi(nbi) = scamas*mtb(nb)
          mtbaqi(nbi) = scamas*mtbaq(nb)
          nt1 = ndrsrd(2,ns1) - ndrsrd(1,ns1) + 1
          if (nt1 .lt. 2) then
c
c           Have a strict basis species associated with the current
c           mass balance.
c
            jflgi(nbi) = jflag(ns2)
          else
c
c           Have an auxiliary basis species associated with the current
c           mass balance.
c
            jflgi(nbi) = 30
          endif
        elseif (ns1.ge.nern1 .and. ns1.le.nern2) then
          mtbi(nbi) = mtb(nb)
          mtbaqi(nbi) = mtbaq(nb)
          jflgi(nbi) = jflag(ns2)
        else
          mtbi(nbi) = mtb(nb)
          mtbaqi(nbi) = mtbaq(nb)
          jflgi(nbi) = jflag(ns2)
        endif
c
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Electrical balance.
c
      electr = scamas*sigzi
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Ordinary basis switching directives.
c
      n = 0
      do kcol = 1,kbt
        nb = iindx1(kcol)
        ns1 = nbaspd(nb)
        ns2 = nbasp(nb)
        if (ns1 .ne. ns2) then
          n = n + 1
          uobsw(1,n) = uspec(ns1)
          uobsw(2,n) = uspec(ns2)
        endif
      enddo
      nobswt = n
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Matrix column variables and corresponding values.
c
      lscama = tlg(scamas)
      zvclg1(kwater) = omeglg
      do kcol = 1,kdim
        nb = iindx1(kcol)
        ns2 = nbasp(nb)
        uzveci(kcol) = uzvec1(kcol)
        if (ns2.ge.narn1 .and. ns2.le.narn2) then
          zvclgi(kcol) = zvclg1(kcol) + lscama
        elseif (ns2.ge.nern1 .and. ns2.le.nern2) then
          zvclgi(kcol) = zvclg1(kcol)
        else
          zvclgi(kcol) = zvclg1(kcol)
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nerr .gt. 0) stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
