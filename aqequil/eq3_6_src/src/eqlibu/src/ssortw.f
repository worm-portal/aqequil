      subroutine ssortw(asort,aval,jsort,nmax,noutpt,nttyo,nval)
c
c     This subroutine sorts the first "nval" elements of the array
c     "aval," using the shell sort method as taken from the "Numerical
c     Recipes" book.
c
c     The sorted elements of "aval" are placed in the array "asort,"
c     and the corresponding indices are placed in the array "jsort."
c
c     If jsort(1) is 0, the sort starts from scratch. Otherwise,
c     the "jsort" array is presumed to contain a prior result which
c     should be a good starting point for the present sorting.
c     All the arrays are dimensioned 1:nmax.
c
c     EQLIBU/qsortw.f is a possible alternative to this subroutine. It
c     uses the "quicksort" method.
c
c     Caution: to be safe, set the first element of the actual array
c     referenced here as jsort to 0 before calling this routine. If
c     this is not done, be sure that nothing has changed the actual
c     array since the previous call for this array. In particular,
c     the length of the array (here represented as nmax should be the
c     same. If it has changed, there is a very good possibility that
c     one or more of the elements of this array are invalid. The present
c     subroutine does only minimal validity checking. Note that with
c     regard to previous calls, include calls to either qsortw.f or
c     ssortw.f, as they are interchangeable in function. Note that
c     even if only one of these subroutines is employed, it is generally
c     used to sort various arrays. Thus, it is not a useful check to
c     save the previous value of nmax for comparison.
c
c     This subroutine is called by:
c
c       Any
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       aval   = array of values to be sorted
c       nmax   = dimension of the asort, aval, istack, jsort, and
c                jstack arrays
c       nval   = the number of elements in the aval array to be sorted,
c                beginning with the first element
c
c     Output:
c
c       asort  = array of sorted values, corresponding to the range
c                composed of the first nval elements fo the aval array
c       jsort  = array of sorted indices, corresponding to the values
c                in the asort array
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nmax
c
      integer noutpt,nttyo
c
      integer jsort(nmax)
      integer nval
c
      real*8 asort(nmax),aval(nmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ii,ileft,ilognp,j,jj,jpt,j2,j3,k,nn,n1
c
      integer ilnobl,jiarmn,jiarmx
c
      logical qx
c
      character*8 ux8,ux8a
c
      real*8 aln2i
c
c-----------------------------------------------------------------------
c
      parameter (aln2i = 1.0/0.63147185)
c
c-----------------------------------------------------------------------
c
      if (jsort(1) .ne. 0) then
c
c       Perform minimal checking of the incoming jsort array.
c
        qx = .false.
        i = jiarmn(jsort,nval)
        if (i .ne. 1) then
          write (ux8,'(i5)') i
          call lejust(ux8)
          j2 = ilnobl(ux8)
          write (noutpt,1000) ux8(1:j2)
          write (nttyo,1000) ux8(1:j2)
 1000     format(/' * Warning - (EQLIBU/ssortw) Programming error',
     $    ' trap: The minimum',/7x,'element of the incoming jsort',
     $    ' array is ',a,', not 1 as is required.',/7x,'Will ignore',
     $    ' this jsort array and create a new one from scratch.')
          qx = .true.
        endif
c
        i = jiarmx(jsort,nval)
        if (i .ne. nval) then
          write (ux8,'(i5)') i
          call lejust(ux8)
          j2 = ilnobl(ux8)
          write (ux8a,'(i5)') nval
          call lejust(ux8a)
          j3 = ilnobl(ux8a)
          write (noutpt,1010) ux8(1:j2),ux8a(1:j3)
          write (nttyo,1010) ux8(1:j2),ux8a(1:j3)
 1010     format(/' * Warning - (EQLIBU/ssortw) Programming error',
     $    ' trap: The maximum',/7x,'element of the incoming jsort',
     $    ' array is ',a,', not nval (',a,'),',/7x,'as is required.',
     $    ' Will ignore this jsort array and create a new one',
     $    /7x,'from scratch.')
          qx = .true.
        endif
c
        if (qx) jsort(1) = 0
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      ileft = (nval/8)*8
      if (jsort(1) .eq. 0) then
c
c       Start from scratch. Note that the loop is unrolled.
c
        do jpt = 1,ileft,8
          jsort(jpt) = jpt
          asort(jpt) = aval(jpt)
          jsort(jpt + 1) = jpt + 1
          asort(jpt + 1) = aval(jpt + 1)
          jsort(jpt + 2) = jpt + 2
          asort(jpt + 2) = aval(jpt + 2)
          jsort(jpt + 3) = jpt + 3
          asort(jpt + 3) = aval(jpt + 3)
          jsort(jpt + 4) = jpt + 4
          asort(jpt + 4) = aval(jpt + 4)
          jsort(jpt + 5) = jpt + 5
          asort(jpt + 5) = aval(jpt + 5)
          jsort(jpt + 6) = jpt + 6
          asort(jpt + 6) = aval(jpt + 6)
          jsort(jpt + 7) = jpt + 7
          asort(jpt + 7) = aval(jpt + 7)
        enddo
c
        do jpt = ileft + 1,nval
          jsort(jpt) = jpt
          asort(jpt) = aval(jpt)
        enddo
c
      else
c
c       Use the pre-existing jsort array as a starting point.
c       Note that the loop is unrolled.
c
        do jpt = 1,ileft,8
          asort(jpt) = aval(jsort(jpt))
          asort(jpt + 1) = aval(jsort(jpt + 1))
          asort(jpt + 2) = aval(jsort(jpt + 2))
          asort(jpt + 3) = aval(jsort(jpt + 3))
          asort(jpt + 4) = aval(jsort(jpt + 4))
          asort(jpt + 5) = aval(jsort(jpt + 5))
          asort(jpt + 6) = aval(jsort(jpt + 6))
          asort(jpt + 7) = aval(jsort(jpt + 7))
        enddo
c
        do jpt = ileft + 1,nval
          asort(jpt) = aval(jsort(jpt))
        enddo
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      nn = nval
      ilognp = int(log(real(nval))*aln2i + 1.e-5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do n1 = 1,ilognp
        nn = nn/2
        k = nval - nn
        do j = 1,k
          i = j
c
c         The label below is the return point for a loop which
c         is run through until no out-of-order elements are found.
c
  120     continue
          jj = i + nn
          if (aval(jsort(jj)) .lt. aval(jsort(i))) then
c
c           Have found an out of order element. Exchange elements.
c
            ii = jsort(i)
            jsort(i) = jsort(jj)
            jsort(jj) = ii
            i = i - nn
            if (i .eq. 1) go to 120
          endif
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Load the asort array. Note that the loop is unrolled.
c
      do jpt = 1,ileft,8
        asort(jpt) = aval(jsort(jpt))
        asort(jpt + 1) = aval(jsort(jpt + 1))
        asort(jpt + 2) = aval(jsort(jpt + 2))
        asort(jpt + 3) = aval(jsort(jpt + 3))
        asort(jpt + 4) = aval(jsort(jpt + 4))
        asort(jpt + 5) = aval(jsort(jpt + 5))
        asort(jpt + 6) = aval(jsort(jpt + 6))
        asort(jpt + 7) = aval(jsort(jpt + 7))
      enddo
c
      do jpt = ileft + 1,nval
        asort(jpt) = aval(jsort(jpt))
      enddo
c
      end
