      subroutine combmb(cdrs,iindx1,ipndx1,jflag,kbt,kdim,km1,kmax,
     $ kmt,kx1,kxt,mtb,mtbaq,ndrsmx,nbasp,nbt,nbtmax,ndrs,ndrsr,
     $ noutpt,nstmax,nttyo,uspec,uzvec1,zvclg1,zvec1)
c
c     This subroutine combines mass balances totals so that the total
c     mass of an active auxiliary basis species whose jflag value
c     is 30 is combined into the total masses of the other basis
c     species which appear in its reaction. This species is then
c     removed from the active basis set.
c
c     This subroutine is called by:
c
c       EQ6/eq6.f
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
      integer kmax,ndrsmx,nbtmax,nstmax
c
      integer noutpt,nttyo
c
      integer iindx1(kmax),ipndx1(kmax),jflag(nstmax),nbasp(nbtmax),
     $ ndrs(ndrsmx),ndrsr(2,nstmax)
c
      integer kbt,kdim,km1,kmt,kx1,kxt,nbt
c
      character*48 uspec(nstmax),uzvec1(kmax)
c
      real*8 cdrs(ndrsmx),mtb(nbtmax),mtbaq(nbtmax),zvclg1(kmax),
     $ zvec1(kmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations with global dimensions.
c
      integer isv_kmax
c
      SAVE isv_kmax
c
      integer, dimension(:), allocatable :: iindxs,ipndxs
c
      SAVE iindxs,ipndxs
c
      real(8), dimension(:), allocatable :: zvclgs
c
      SAVE zvclgs
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jlen,k,kbts,km1s,kmts,krow,kx1s,kxts,n,nb,nbb,nerr,np,
     $ nr1,nr2,ns,nss,nt
c
      integer nbasis
c
      character*56 uspn56
c
      real*8 cx,mx,lx
c
      real*8 texp
c
c-----------------------------------------------------------------------
c
c     Allocate or reallocate local work arrays as needed.
c
      if (.not.ALLOCATED(iindxs)) then
c
c       Local work arrays are not allocated. Zero the saved
c       array size variables. Note that only one array is tested
c       to see if it is allocated. It is assumed that all local
c       work arrays are either allocated or not.
c
        isv_kmax = 0
      else
c
c       Local work arrays are allocated. Check to see if any of the
c       array size variables have changed. If so, deallocate
c       the corresponding local work arrays and zero the corresponding
c       saved size variables.
c
        if (kmax .ne. isv_kmax) then
          DEALLOCATE(iindxs,ipndxs,zvclgs)
          isv_kmax = 0
        endif
      endif
c
c     At this point, the saved array size values are zero if the
c     corresponding arrays need to be allocated.
c
      if (isv_kmax .eq. 0) then
        ALLOCATE(iindxs(kmax),ipndxs(kmax))
        ALLOCATE(zvclgs(kmax))
        isv_kmax = kmax
      endif
c
c     Zero the contents of the local work arrays.
c
      do k = 1,kmax
        iindxs(k) = 0
        ipndxs(k) = 0
      enddo
c
      do k = 1,kmax
        zvclgs(k) = 0.
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      nerr = 0
      kbts = kbt
      km1s = km1
      kmts = kmt
      kx1s = kx1
      kxts = kxt
c
      do krow = 1,kdim
        iindxs(krow) = iindx1(krow)
        ipndxs(krow) = ipndx1(krow)
        zvclgs(krow) = zvclg1(krow)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      kbt = 0
      do krow = 1,kbts
        nb = iindxs(krow)
        np = ipndxs(krow)
        ns = nbasp(nb)
        if (jflag(ns) .eq. 0) then
          kbt = kbt + 1
          uzvec1(kbt) = uspec(ns)
          iindx1(kbt) = nb
          ipndx1(kbt) = np
          zvclg1(kbt) = zvclgs(krow)
        else
          nr1 = ndrsr(1,ns)
          nr2 = ndrsr(2,ns)
          nt = nr2 - nr1 + 1
          if (nt .lt. 2) then
c
c           Calling sequence substitutions:
c             uspec(ns) for unam48
c
            call fmspnm(jlen,uspec(ns),uspn56)
            write (noutpt,1000) uspn56(1:jlen),jflag(ns)
            write (nttyo,1000) uspn56(1:jlen),jflag(ns)
 1000       format(/' * Error- (combmb) The species ',a,
     $      /7x,'has jflag value of ',i2,", but it's a strict basis",
     $      ' species.',/7x,"Its mass balance can't be combined into",
     $      ' that of another basis species.')
            nerr = nerr + 1
            go to 100
          endif
          do n = nr1 + 1,nr2
            cx = -cdrs(n)/cdrs(nr1)
            nss = ndrs(n)
c
c           Calling sequence substitutions:
c             nss for ns
c
            nbb = nbasis(nbasp,nbt,nbtmax,nss)
            mtb(nbb) = mtb(nbb) + cx*mtb(nb)
            mtbaq(nbb) = mtbaq(nbb) + cx*mtbaq(nb)
          enddo
          mtb(nb) = 0.
          mtbaq(nb) = 0.
        endif
  100   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      kdim = kbt
      km1 = kdim + 1
      kmt = kdim
      if (kmts .ge. km1s) then
        do krow = km1s,kmts
          kmt = kmt + 1
          ns = iindxs(krow)
          uzvec1(kmt) = uspec(ns)
          iindx1(kmt) = ns
          ipndx1(kmt) = ipndxs(krow)
          zvclg1(kmt) = zvclgs(krow)
        enddo
        kdim = kmt
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      kx1 = kdim + 1
      kxt = kdim
      if (kxts .ge. kx1s) then
        do krow = kx1s,kxts
          kxt = kxt + 1
          ns = iindxs(krow)
          uzvec1(kxt) = uspec(ns)
          iindx1(kxt) = ns
          ipndx1(kxt) = ipndxs(krow)
          zvclg1(kxt) = zvclgs(krow)
        enddo
        kdim = kxt
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do krow = 1,kdim
        lx = zvclg1(krow)
        mx = texp(lx)
        zvec1(krow) = mx
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nerr .gt. 0) stop
c
      end
