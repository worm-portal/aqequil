      subroutine chdxgz(delxi,dlxmx0,fdlim,fdzv0,iodb,iopt,jordlm,
     $ kdim,kmax,km1,kord,kxt,nodbmx,noptmx,nordmx,nordz,noutpt,
     $ nrd1mx,nsscmx,scalim,scfcz,sscrew,qmin,smp100,uzvec1,zklogu,
     $ zvec0,zvclg0)
c
c     This subroutine chooses a step size and order according to the
c     Gear accuracy criterion, examining the z vector and its
c     associated finite differences. Subroutine chdxgr.f performs
c     the same function for the r vector and its associated finite
c     differences.
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
      integer kmax,nodbmx,noptmx,nordmx,nrd1mx,nsscmx
c
      integer noutpt
c
      integer iodb(nodbmx),iopt(noptmx)
c
      integer jordlm,kdim,km1,kord,kxt,nordz
c
      logical qmin
c
      character(len=48) uzvec1(kmax)
c
      real(8) fdzv0(nrd1mx,kmax),sscrew(nsscmx),zvclg0(kmax),zvec0(kmax)
c
      real(8) delxi,dlxmx0,fdlim,scalim,scfcz,smp100,zklogu
c
c-----------------------------------------------------------------------
c
c     Local variable declarations with global dimensions.
c
      integer isv_nordmx
c
      SAVE isv_nordmx
c
      integer, dimension(:), allocatable :: kkk
c
      SAVE kkk
c
      character(len=24), dimension(:), allocatable :: uscal
c
      SAVE uscal
c
      real(8), dimension(:), allocatable :: scmax,scraw
c
      SAVE scmax,scraw
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ip1,j2,j3,k,kcol,kcoli
c
      integer ilnobl
c
      character(len=24) ubas,uzt,uzn
c
      real*8 adx,delchg,dlxii,fctr,fctrmn,lxx0,mxx0,relchg,xip1
c
c-----------------------------------------------------------------------
c
      data ubas /'Basis variable          '/
c
c-----------------------------------------------------------------------
c
c     Allocate or reallocate local work arrays as needed.
c
      if (.not.ALLOCATED(kkk)) then
c
c       Local work arrays are not allocated. Zero the saved
c       array size variables. Note that only one array is tested
c       to see if it is allocated. It is assumed that all local
c       work arrays are either allocated or not.
c
        isv_nordmx = 0
      else
c
c       Local work arrays are allocated. Check to see if any of the
c       array size variables have changed. If so, deallocate
c       the corresponding local work arrays and zero the corresponding
c       saved size variables.
c
        if (nordmx .ne. isv_nordmx) then
          DEALLOCATE(kkk,uscal,scmax,scraw)
          isv_nordmx = 0
        endif
      endif
c
c     At this point, the saved array size values are zero if the
c     corresponding arrays need to be allocated.
c
      if (isv_nordmx .eq. 0) then
        ALLOCATE(kkk(0:nordmx))
        ALLOCATE(uscal(0:nordmx))
        ALLOCATE(scmax(0:nordmx),scraw(0:nordmx))
        isv_nordmx = nordmx
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Zero the contents of the local work arrays.
c
      do k = 0,nordmx
        kkk(k) = 0
      enddo
c
      do k = 0,nordmx
        uscal(k) = ' '
      enddo
c
      do k = 0,nordmx
        scmax(k) = 0.
        scraw(k) = 0.
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (kord .le. 0) then
        nordz = 0
        scfcz = dlxmx0/delxi
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set up the regular choice of orders. Here i is the order, which
c     technically runs from 0 to jordlm. Here jordlm + 1 is the hidden
c     extra order, which is used to constrain the error in using order
c     jordlm.
c
      do i = 0,jordlm
        scmax(i) = 0.
        scraw(i) = 0.
      enddo
c
c     For i = 0, take delxi = dlxmx0. This keeps delxi greater
c     than or equal to dlxmx0, except within searches, as for phase
c     boundaries.
c
      scmax(0) = dlxmx0/delxi
      scraw(0) = scmax(0)
      kkk(0) = 0
      uscal(0) = 'None'
c
c     Process the higher orders.
c
      do i = 1,jordlm
        ip1 = i + 1
        dlxii = delxi**(ip1)
        fctrmn = 1.e+38
        kkk(i) = 0
        uscal(i) = 'None'
c
c       Constrain the step size to satisfy the tolerance on the
c       fractional error of matrix variables.
c
        do kcol = 1,kdim
          qmin = kcol.ge.km1 .and. kcol.le.kxt
          lxx0 = zvclg0(kcol)
c
c         Skip the constraint for a mineral with a very small mass.
c
          if (.not.qmin .or. lxx0.gt.zklogu) then
c
            adx = abs(fdzv0(ip1,kcol))
            delchg = adx*dlxii
            if (delchg .ne. 0.) then
c
c             Taylor's series in a linear quantity (e.g., mass).
c
              mxx0 = zvec0(kcol)
              if (mxx0 .ge. smp100) then
                relchg = delchg/mxx0
              else
                relchg = fdlim
              endif
c
              fctr = sscrew(1)/relchg
              if (fctr .lt. fctrmn) then
                fctrmn = fctr
                uscal(i)(1:24) = ubas(1:24)
                kkk(i) = kcol
              endif
            endif
          endif
        enddo
c
        xip1 = ip1
        scmax(i) = fctrmn**(1./xip1)
        scraw(i) = scmax(i)
        if (scmax(i) .gt. scalim) scmax(i) = scalim
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(5) .ge. 2) then
        write (noutpt,1000)
 1000   format(/' --- Order and Scale Factor Data (chdxgz) ---',/)
        do i = 0,jordlm
          uzt = uscal(i)
          kcoli = kkk(i)
          if (uzt(1:24) .eq. ubas(1:24)) then
            uzn = uzvec1(kcoli)(1:24)
          else
            uzn = ' '
          endif
c
          write (noutpt,1010) i,scmax(i),scraw(i)
 1010     format('   Order= ',i2,', scale= ',f10.5,', raw scale= ',
     $    f10.5)
          j2 = ilnobl(uzt)
          j3 = ilnobl(uzn)
          if (j3 .gt. 0) then
            write (noutpt,1020) uzt(1:j2),uzn(1:j3)
 1020       format(5x,'Constrained by ',a,1x,a)
          else
            write (noutpt,1030) uzt(1:j2)
 1030       format(5x,'Constrained by ',a)
          endif
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Choose the order that gives the highest scale factor.
c
      scfcz = scmax(0)
      nordz = 0
      do i = 1,jordlm
        if (scmax(i) .gt. scfcz) then
          scfcz = scmax(i)
          nordz = i
        endif
      enddo
c
  999 continue
      end
