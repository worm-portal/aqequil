      subroutine absswb(adhfs,adhfsx,advfs,advfsx,avcnst,axhfs,
     $ axhfsx,axlks,axlksx,axvfs,axvfsx,cdrs,cdrtw,cdrsx,cdrw,
     $ csts,dhfs,dvfs,eps100,ibswx,iindx1,iodb,ipch,ipchmx,ipcv,
     $ ipcvmx,jcsort,jflag,jsflag,kbt,kmax,mosp,mtb,narn1,narn2,
     $ narxmx,narxt,nbasp,nbaspd,nbaspx,nbt,nbtmax,nbw,ndrs,ndrsmx,
     $ ndrsr,ndrsrx,ndrsx,nelect,nhydr,nodbmx,no2gaq,noutpt,nst,
     $ nstmax,nsts,nstsmx,nstsr,nswtch,ntpr,ntprmx,nttyo,presg,
     $ press,qbassw,qbswx,tempc,uspec,uzvec1,weight,xhfs,xvfs,xlks)
c
c     This subroutine carries out automatic basis switching to optimize
c     optimize the basis set after Newton-Raphson iteration has
c     converged at the current point of reaction progress. Such
c     optimization tends to improve the code numerics (finite
c     difference stability, matrix conditioning, etc.) This subroutine
c     is similar to EQLIB/absswa.f, which carries out automatic basis
c     switching as a pre-Newton-Raphson optimization technique.
c
c     This subroutine is called by:
c
c       EQ6/path.f
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
      integer ipchmx,ipcvmx,kmax,narxmx,nbtmax,ndrsmx,nodbmx,nstmax,
     $ nstsmx,ntprmx
c
      integer noutpt,nttyo
c
      integer ibswx(nbtmax),iindx1(kmax),iodb(nodbmx),jcsort(nstmax),
     $ jflag(nstmax),jsflag(nstmax),narxt(ntprmx),nbasp(nbtmax),
     $ nbaspd(nbtmax),nbaspx(nbtmax),ndrs(ndrsmx),ndrsr(2,nstmax),
     $ ndrsrx(2,nstmax),ndrsx(ndrsmx),nsts(nstsmx),nstsr(2,nstmax)
c
      integer ipch,ipcv,kbt,narn1,narn2,nbt,nbw,nelect,nhydr,no2gaq,
     $ nst,nswtch,ntpr
c
      logical qbassw,qbswx
c
      character*48 uspec(nstmax),uzvec1(kmax)
c
      real*8 adhfs(narxmx,ntprmx,ipchmx,nstmax),
     $ adhfsx(narxmx,ntprmx,ipchmx,nstmax),
     $ advfs(narxmx,ntprmx,ipcvmx,nstmax),
     $ advfsx(narxmx,ntprmx,ipcvmx,nstmax),
     $ axhfs(narxmx,ntprmx,nstmax),axhfsx(narxmx,ntprmx,nstmax),
     $ axlks(narxmx,ntprmx,nstmax),axlksx(narxmx,ntprmx,nstmax),
     $ axvfs(narxmx,ntprmx,nstmax),axvfsx(narxmx,ntprmx,nstmax),
     $ cdrs(ndrsmx),cdrsx(ndrsmx),cdrtw(nstmax),cdrw(nstmax),
     $ csts(nstsmx),dhfs(ipchmx,nstmax),dvfs(ipcvmx,nstmax),
     $ mosp(nstmax),mtb(nbtmax),weight(nstmax),xhfs(nstmax),
     $ xlks(nstmax),xvfs(nstmax)
c
      real*8 avcnst,presg,press,tempc
c
      real*8  eps100
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jlen1,jlen2,jlen3,krow,nb,nb1,nb2,ns,nse,nsi,nsj,ns1,ns2
c
      character*56 usp156,usp256,usp356
c
      real*8 api,apj,cx,fx1,fx2,rx,wb1,wb2,wsi
c
      real*8 coefdr,coefst
c
c-----------------------------------------------------------------------
c
c     Initialize the number of basis switches made.
c
      nswtch = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Pick candidates for automatic basis switching. The indices of
c     candidates are stored in the array ibswx. Note that the method
c     employed here is similar to that in EQLIB/abswpk.f. However,
c     there are some differences. Here, basis switching is not being
c     used as a method of pre-Newton-Raphson optimization. Hence, there
c     are no residual functions used in defining criteria for making
c     a switch. the relative contributions to mass balances are used
c     instead.
c
      qbswx = .false.
      do nb = 1,nbt
        ibswx(nb) = 0
        nse = nbaspd(nb)
        nsj = nbasp(nb)
        if (nse.ne.narn1 .and. nse.ne.nhydr .and. nse.ne.no2gaq
     $    .and.nse.ne.nelect) then
c
c         If the current species is not H2O, H+, O2(g,aq), or e-,
c         consider a switch.
c
c         Get weights for the current mass balance.
c
          do ns = narn1,narn2
            weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
          enddo
c
c         Screen out species in the mass balance that are not linked
c         by the current reaction with the current basis species.
c         Do this by seting the corresponding weights to zero.
c
          do ns = narn1,narn2
            if (jflag(ns) .eq. 30) then
c
c             Calling sequence substitutions:
c               nsj for nse
c
              cx = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsj,ns,nstmax)
              if (cx .eq. 0.) weight(ns) = 0.
            endif
          enddo
c
c         Find a candidate.
c
          call fbassw(jcsort,jflag,mosp,narn1,narn2,nse,nsi,nsj,
     $    nstmax,weight,wsi)
c
c         If a candidate was found, apply a filter.
c
          if (nsi .gt. 0) then
c
c           Use the current candidate species if its contribution
c           to the mass balance exceeds that of the current basis
c           species by a reasonable factor, here 10.
c
            api = weight(nsi)*mosp(nsi)
            apj = weight(nsj)*mosp(nsj)
            rx = api/apj
            if (rx .gt. 10.) then
              ibswx(nb) = nsi
              qbswx = .true.
            endif
          endif
        endif
      enddo
c
c     If there are no candidates for basis switching, exit the
c     present subroutine.
c
      if (.not.qbswx) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(1) .ge. 2) then
        write (noutpt,1000)
 1000   format(/5x,'--- Candidate Basis Switches ---',/)
        do nb = 1,nbt
          nse = nbaspd(nb)
          nsj = nbasp(nb)
          nsi = ibswx(nb)
          if (nsi .gt. 0) then
c
c           Calling sequence substitutions:
c             jlen1 for jlen
c             uspec(nsj) for unam48
c             usp156 for uspn56
c
            call fmspnx(jlen1,uspec(nsj),usp156)
c
c           Calling sequence substitutions:
c             jlen2 for jlen
c             uspec(nsi) for unam48
c             usp256 for uspn56
c
            call fmspnx(jlen2,uspec(nsi),usp256)
c
c           Calling sequence substitutions:
c             jlen3 for jlen
c             uspec(nse) for unam48
c             usp356 for uspn56
c
            call fmspnx(jlen3,uspec(nse),usp356)
c
            write (noutpt,1010) usp156(1:jlen1),usp256(1:jlen2),
     $      usp356(1:jlen3)
 1010       format(/3x,'Could replace ',a,' in the active basis set',
     $      ' with',/3x,a,' as the species associated with the mass',
     $      ' balance',/3x,'of ',a,'.')
          endif
        enddo
        write (noutpt,1020)
 1020   format(1x)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Resolve any conflicts in candidate basis switches. The algorithm
c     used here is similar to but not identical to that employed by
c     EQLIB/gabswx.f, which resolves conflicts in the context of using
c     automatic basis switching as a form of pre-Newton-Raphson
c     optimization.
c
      do nb1 = 1,nbt - 1
        ns1 = ibswx(nb1)
        if (ns1 .gt. 0) then
          do nb2 = nb1 + 1,nbt
            ns2 = ibswx(nb2)
            if (ns2 .eq. ns1) then
c
c             Calling sequence substitutions:
c               nb1 for nb
c
              wb1 = coefst(csts,nsts,nstsmx,nstsr,nb1,ns,nstmax)
              fx1 = wb1*mosp(ns1)/mtb(nb1)
c
c             Calling sequence substitutions:
c               nb2 for nb
c
              wb2 = coefst(csts,nsts,nstsmx,nstsr,nb2,ns,nstmax)
              fx2 = wb2*mosp(ns2)/mtb(nb2)
              if (fx1 .gt. fx2) then
                ibswx(nb2) = 0
              else
                ibswx(nb1) = 0
                go to 100
              endif
            endif
          enddo
        endif
  100   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      nswtch = 0
      do nb = 1,nbt
        nsi = ibswx(nb)
        if (nsi .gt. 0) nswtch = nswtch + 1
      enddo
c
      if (nswtch .le. 0) then
        if (iodb(1) .ge. 2) write (noutpt,1070)
 1070   format(/10x,'No switches will be made.',/)
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(1) .ge. 2) then
        write (noutpt,1100)
 1100   format(/5x,'--- Automatic Basis Switches ---',/)
        do nb = 1,nbt
          nse = nbaspd(nb)
          nsj = nbasp(nb)
          nsi = ibswx(nb)
          if (nsi .gt. 0) then
c
c           Calling sequence substitutions:
c             jlen1 for jlen
c             uspec(nsj) for unam48
c             usp156 for uspn56
c
            call fmspnx(jlen1,uspec(nsj),usp156)
c
c           Calling sequence substitutions:
c             jlen2 for jlen
c             uspec(nsi) for unam48
c             usp256 for uspn56
c
            call fmspnx(jlen2,uspec(nsi),usp256)
c
c           Calling sequence substitutions:
c             jlen3 for jlen
c             uspec(nse) for unam48
c             usp356 for uspn56
c
            call fmspnx(jlen3,uspec(nse),usp356)
c
            write (noutpt,1110) usp156(1:jlen1),usp256(1:jlen2),
     $      usp356(1:jlen3)
 1110       format(/3x,'Will replace ',a,' in the active basis set',
     $      ' with',/3x,a,' as the species associated with the mass',
     $      ' balance',/3x,'of ',a,'.')
          endif
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      call autosw(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,
     $ axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,eps100,ibswx,iindx1,
     $ ipch,ipchmx,ipcv,ipcvmx,jflag,jsflag,kbt,kmax,narn1,narxmx,
     $ nbasp,nbaspd,nbaspx,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsx,ndrsr,
     $ ndrsrx,noutpt,nst,nstmax,ntprmx,nttyo,qbassw,uspec,uzvec1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Recompute the cdrw array.
c
      call gcdrw(cdrs,cdrw,narn1,ndrs,ndrsmx,ndrsr,nst,nstmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Recompute the cdrtw array.
c
      call gcdrtw(cdrs,cdrtw,narn1,narn2,ndrs,ndrsmx,ndrsr,
     $ nelect,no2gaq,nst,nstmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Update the thermodynamic data to correspond to the new
c     active basis set. First, recompute the log K, etc., data for
c     the various reactions.
c
      call evdatr(adhfs,advfs,axhfs,axlks,axvfs,dhfs,dvfs,
     $ ipch,ipchmx,ipcv,ipcvmx,narxmx,narxt,nst,nstmax,ntpr,ntprmx,
     $ tempc,xhfs,xlks,xvfs)
c
c     Then make pressure corrections to these thermodynamic data.
c
      if (ipcv .ge. 0) then
        call pcorrx(avcnst,dhfs,dvfs,ipch,ipchmx,ipcv,ipcvmx,
     $  nbasp,nbt,nbtmax,ndrsr,nst,nstmax,presg,press,xhfs,
     $  xlks,xvfs)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
c
      end
