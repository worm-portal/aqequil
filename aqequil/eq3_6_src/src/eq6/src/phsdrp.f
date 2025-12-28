      subroutine phsdrp(d1zvc1,iindx0,iindx1,iodb,ipndx1,iter,kmax,
     $ km1,km1s,kxt,kxts,nodbmx,nord,noutpt,npadd,npdel,nptmax,ntry,
     $ uphase,zvclgs,zvclg1)
c
c     This subroutine picks a phase to drop from the equilibrium
c     system. Two independent algorithms are used to find candidates,
c     and a final choice is made from these.
c
c     This subroutine is called by:
c
c       EQ6/eqcalc.f
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
      integer kmax,nodbmx,nptmax
c
      integer noutpt
c
      integer iindx0(kmax),iindx1(kmax),ipndx1(kmax),iodb(nodbmx)
c
      integer iter,km1,km1s,kxt,kxts,nord,npadd,npdel,ntry
c
      character*24 uphase(nptmax)
c
      real*8 d1zvc1(kmax),zvclgs(kmax),zvclg1(kmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,j3,kcol,kj,krow,npdc1,npdc1a,npdc2,np,ns
c
      integer ilnobl
c
      character*24 updc1,updc1a,updc2
c
      real*8 crit1,crit1a,crit2,criter,zx
c
c-----------------------------------------------------------------------
c
      npdc1 = 0
      npdc2 = 0
      crit1 = 0.
      crit2 = 0.
      updc1 = 'None'
      updc2 = 'None'
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Search for a candidate solid according to the first criterion
c     (rapid decrease in mass).
c
      npdc1a = 0
      updc1a = 'None'
      crit1a = 0.
c
      do kcol = km1,kxt
        zx = zvclg1(kcol) - zvclgs(kcol)
        if (zx .ge. -1.0) zx = 0.
        if (iter .le. 1) zx = zvclgs(kcol) + 10.
        if (zvclgs(kcol) .le. -100.) zx = zvclgs(kcol)
        if (zx .lt. crit1) then
          ns = iindx1(kcol)
          np = ipndx1(kcol)
          crit1a = zx
          npdc1a = np
          updc1a = uphase(np)
          if (np .ne. npadd) then
            crit1 = zx
            npdc1 = np
            updc1 = uphase(np)
          endif
        endif
      enddo
c
      if (npdc1.le.0) crit1 = 0.
      if (npdc1a.gt.0 .and. ntry.eq.2) then
        if (crit1a. lt. (crit1 -1.)) then
          npdc1 = npdc1a
          updc1 = updc1a
          crit1 = crit1a
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Search for a candidate solid according to the second criterion
c     (negative derivative of mass).
c
      if (nord.gt.0 .and. ntry.le.0) then
        crit2 = -100.
        if (kxt.ge.km1 .and. kxts.ge.km1) then
          do kcol = km1,kxt
            kj = 0
            ns = iindx1(kcol)
            np = ipndx1(kcol)
            do krow = km1s,kxts
              if (iindx0(krow) .eq. ns) then
                kj = krow
                go to 100
              endif
            enddo
  100       if (kj .gt. 0) then
              zx = d1zvc1(kj)
              if (zx .lt. crit2) then
                if (np .ne. npadd) then
                  crit2 = zx
                  npdc2 = np
                  updc2 = uphase(np)
                endif
              endif
            endif
          enddo
        endif
      endif
      if (npdc2.le.0) crit2 = 0.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(1) .ge. 1) then
        j2 = ilnobl(updc1)
        j3 = ilnobl(updc2)
        write (noutpt,1000) crit1,updc1(1:j2),crit2,updc2(1:j3)
 1000   format(/'   --- Phase Drop Search ---',
     $  /7x,'Criterion 1. Negative divergence of log mass variable',
     $  /11x,'Value= ',e12.5,4x,/11x,'Name= ',a,
     $  /7x,'Criterion 2. Negative derivative of log mass variable',
     $  /11x,'Value= ',e12.5,4x,/11x,'Name= ',a,/)
      endif
c
      criter = 0.
      npdel = 0
      if (crit1 .lt. criter) then
        criter = crit1
        npdel = npdc1
      endif
      if (crit2 .lt. criter) then
        criter = crit2
        npdel = npdc2
      endif
c
      end
