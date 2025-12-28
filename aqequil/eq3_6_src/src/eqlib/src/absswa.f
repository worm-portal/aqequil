      subroutine absswa(adhfs,adhfsx,advfs,advfsx,avcnst,axhfs,
     $ axhfsx,axlks,axlksx,axvfs,axvfsx,beta,cdrs,cdrsx,cdrtw,cdrw,
     $ csts,dhfs,dvfs,efac,eps100,ibswx,iebal,iindx1,iodb,ipch,ipchmx,
     $ ipcv,ipcvmx,jcsort,jflag,jsflag,jssort,kbt,kmax,mosp,narn1,
     $ narn2,narxmx,narxt,nbasp,nbaspd,nbaspx,nbt,nbtmax,nbw,ncosp,
     $ ndrs,ndrsmx,ndrsr,ndrsrx,ndrsx,nelect,nhydr,nodbmx,no2gaq,
     $ noutpt,nst,nstmax,nsts,nstsmx,nstsr,nswtch,ntpr,ntprmx,nttyo,
     $ presg,press,qbassw,qbswx,q6mode,tempc,uspec,uzvec1,weight,
     $ xvfs,xlks,xhfs)
c
c     This subroutine carries out automatic basis switching as a
c     pre-Newton-Raphson optimization technique (to reduce the
c     magnitude of very large mass balance residuals). It is similar
c     to EQ6/absswb.f, which carries out automatic basis switching to
c     optimize the basis set after the system has been solved at the
c     latest point of reaction progress.
c
c     This subroutine is called by:
c
c       EQ3NR/arrset.f
c       EQ6/optmzr.f
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
     $ jflag(nstmax),jsflag(nstmax),jssort(nstmax),narxt(ntprmx),
     $ nbasp(nbtmax),nbaspd(nbtmax),nbaspx(nbtmax),ncosp(nbtmax),
     $ ndrs(ndrsmx),ndrsr(2,nstmax),ndrsrx(2,nstmax),ndrsx(ndrsmx),
     $ nsts(nstsmx),nstsr(2,nstmax)
c
      integer iebal,ipch,ipcv,kbt,narn1,narn2,nbt,nbw,nelect,nhydr,
     $ no2gaq,nst,nswtch,ntpr
c
      logical qbassw,qbswx,q6mode
c
      character(len=48) uspec(nstmax),uzvec1(kmax)
c
      real(8) adhfs(narxmx,ntprmx,ipchmx,nstmax),
     $ adhfsx(narxmx,ntprmx,ipchmx,nstmax),
     $ advfs(narxmx,ntprmx,ipcvmx,nstmax),
     $ advfsx(narxmx,ntprmx,ipcvmx,nstmax),
     $ axhfs(narxmx,ntprmx,nstmax),axhfsx(narxmx,ntprmx,nstmax),
     $ axlks(narxmx,ntprmx,nstmax),axlksx(narxmx,ntprmx,nstmax),
     $ axvfs(narxmx,ntprmx,nstmax),axvfsx(narxmx,ntprmx,nstmax),
     $ cdrs(ndrsmx),cdrsx(ndrsmx),cdrtw(nstmax),cdrw(nstmax),
     $ csts(nstsmx),beta(kmax),dhfs(ipchmx,nstmax),dvfs(ipcvmx,nstmax),
     $ efac(nbtmax),mosp(nstmax),weight(nstmax),xhfs(nstmax),
     $ xlks(nstmax),xvfs(nstmax)
c
      real(8) avcnst,presg,press,tempc
c
      real(8) eps100
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jfl,jlen1,jlen2,jlen3,j2,nb,nb1,nb2,ncount,
     $ nse,nsi,nsj,nss,ns2
c
      integer ilnobl,nbasis
c
      character(len=56) usp156,usp256,usp356
      character(len=8) ux8
c
c-----------------------------------------------------------------------
c
c     Initialize the number of basis switches made.
c
      nswtch = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Pick candidates for automatic basis switching. This involves
c     a search similar to the one above, but with a variety of special
c     constraints. The indices of candidates are stored in the array
c     ibswx.
c
      call abswpk(beta,cdrs,csts,efac,ibswx,iebal,iindx1,
     $ jcsort,jflag,jssort,kbt,kmax,mosp,narn1,narn2,nbasp,nbaspd,
     $ nbt,nbtmax,ndrs,ndrsmx,ndrsr,nelect,nhydr,no2gaq,nstmax,
     $ nsts,nstsmx,nstsr,qbswx,q6mode,weight)
c
      if (.not.qbswx) then
        nswtch = 0
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(3) .ge. 1) then
        write (noutpt,1000)
 1000   format(/5x,'--- Candidate Basis Switches ---')
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
      if (.not.q6mode) then
c
c       Running EQ3NR. Wipe out any proposed basis switches in which
c       the species to be switched out is the other species involved
c       in a jflag = 17 18, or 21 option.
c
        ncount = 0
        do nb1 = 1,nbt
          nse = nbaspd(nb1)
          nsj = nbasp(nb1)
          jfl = jflag(nse)
          if (jfl.eq.17 .or. jfl.eq.18 .or. jfl.eq.21) then
            ns2 = ncosp(nb1)
            do nb2 = 1,nbt
              nss = nbasp(nb2)
              if (ns2 .eq. nss) then
                if (ibswx(nb2) .ne. 0) then
                  if (iodb(3) .ge. 2) then
c
c                   Calling sequence substitutions:
c                     jlen2 for jlen
c                     uspec(ns2) for unam48
c                     usp256 for uspn56
c
                    call fmspnx(jlen2,uspec(ns2),usp256)
c
c                   Calling sequence substitutions:
c                     jlen3 for jlen
c                     uspec(nse) for unam48
c                     usp356 for uspn56
c
                    call fmspnx(jlen3,uspec(nse),usp356)
c
                    write (ux8,'(i5)') jfl
                    call lejust(ux8)
                    j2 = ilnobl(ux8)
c
                    write (noutpt,1050) usp256(1:jlen2),ux8(1:j2),
     $              usp356(1:jlen3)
                    write (nttyo,1050) usp256(1:jlen2),ux8(1:j2),
     $              usp356(1:jlen3)
 1050               format(/" * Note - (EQLIB/absswa) Can't switch ",
     $              ' the species ',a,/7x,'out of the active basis set',
     $              ' because it is tied up in a jflag = ',a,' option',
     $              /7x,'for ',a,'.')
                    ncount = ncount + 1
                  endif
                  ibswx(nb2) = 0
                endif
                go to 100
              endif
            enddo
  100       continue
          endif
        enddo
        if (ncount .gt. 0) write (noutpt,1020)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Resolve any conflicts in candidate basis switches.
c
      call gabswx(beta,ibswx,iindx1,kbt,kmax,nbt,nbtmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Count the number of switches to make.
c
      nswtch = 0
      do nb = 1,nbt
        nsi = ibswx(nb)
        if (nsi .gt. 0) nswtch = nswtch + 1
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nswtch .le. 0) then
        if (iodb(3) .ge. 1) write (noutpt,1070)
 1070   format(/10x,'No switches will be made.',/)
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(3) .ge. 1) then
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
        write (noutpt,1020)
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
