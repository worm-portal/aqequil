      subroutine chdxgr(delxi,dlxmx0,fdri0,fdrr0,iodb,jordlm,jreac,
     $ kord,nodbmx,nordmx,nordr,noutpt,nrct,nrctmx,nrd1mx,nsscmx,
     $ scalim,scfcr,sscrew,qriinf,rirec0,rrelr0,ureac)
c
c     This subroutine chooses a step size and order according to the
c     Gear accuracy criterion, examining the r vector and its
c     associated finite differences. Subroutine chdxgz.f performs
c     the same function for the z vector and its associated finite
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
      integer nodbmx,nrctmx,nordmx,nrd1mx,nsscmx
c
      integer noutpt
c
      integer iodb(nodbmx),jreac(nrctmx)
c
      integer jordlm,kord,nordr,nrct
c
      logical qriinf
c
      character(len=24) ureac(nrctmx)
c
      real(8) fdri0(nrd1mx),fdrr0(nrd1mx,nrctmx),rrelr0(nrctmx),
     $ sscrew(nsscmx)
c
      real(8) delxi,dlxmx0,rirec0,scalim,scfcr
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
      integer i,ip1,j2,j3,k,nrc
c
      integer ilnobl
c
      character*24 uinv,urel,uzt,uzn
c
      real(8) adx,delchg,dlxii,fctr,fctrmn,relchg,rxx0,xip1
c
c-----------------------------------------------------------------------
c
      data uinv /'Inverse rate            '/,
     $     urel /'Relative rate           '/
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
      do k = 1,nordmx
        kkk(k) = 0
      enddo
c
      do k = 1,nordmx
        uscal(k) = ' '
      enddo
c
      do k = 1,nordmx
        scmax(k) = 0.
        scraw(k) = 0.
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (kord .le. 0) then
        nordr = 0
        scfcr= dlxmx0/delxi
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
c       Get results based on the inverse rate.
c
        if (.not.qriinf) then
          rxx0 = rirec0
          if (rxx0 .ne. 0.) then
            adx = abs(fdri0(ip1))
            delchg = adx*dlxii
            relchg = delchg/rxx0
            if (relchg .ne. 0.) then
              fctr = sscrew(3)/relchg
              if (fctr .lt. fctrmn) then
                fctrmn = fctr
                uscal(i)(1:24) = uinv(1:24)
                kkk(i) = 0
              endif
            endif
          endif
        endif
c
c       Get results based on relative rates.
c
        do nrc = 1,nrct
          if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1)  then
            rxx0 = abs(rrelr0(nrc))
            if (rxx0 .ne. 0.) then
              adx = abs(fdrr0(ip1,nrc))
              delchg = adx*dlxii
              if (delchg .ne. 0.) then
                fctr = sscrew(3)/delchg
                if (fctr .lt. fctrmn) then
                  fctrmn = fctr
                  uscal(i)(1:24) = urel(1:24)
                  kkk(i) = nrc
                endif
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
 1000   format(/' --- Order and Scale Factor Data (chdxgr) ---',/)
        do i = 0,jordlm
          uzt = uscal(i)
          nrc = kkk(i)
          if (uzt(1:24) .eq. urel(1:24)) then
            uzn = ureac(nrc)
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
      scfcr = scmax(0)
      nordr = 0
      do i = 1,jordlm
        if (scmax(i) .gt. scfcr) then
          scfcr = scmax(i)
          nordr = i
        endif
      enddo
c
  999 continue
      end
