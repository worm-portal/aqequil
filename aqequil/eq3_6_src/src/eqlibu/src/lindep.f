      subroutine lindep(aamatr,eps100,irow1,irow2,jcol1,jcol2,kmax,
     $ qldep)
c
c     This subroutine determines whether the rows of a submatrix of
c     matrix aamatr are linearly dependent or not. The submatrix runs
c     from rows irow1 to irow2 and columns jcol1 to jcol2. If there is
c     linear dependence, qldep is returned with a value of ".true.".
c     This subroutine factors all rows if possible. It is assumed that
c     the calling subroutine may use these rows if qldep is returned
c     with a value of ".true.".
c
c     This subroutine is called by:
c
c       EQ3NR/dawfix.f
c       EQ6/eqcalc.f
c       EQ6/jgibbs.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       aamatr = the matrix to be tested
c       eps100 = 100 times the real*8 machine epsilon
c       irow1  = start of range of rows to examine
c       irow2  = end of range of rows to examine
c       jcol1  = start of range of columns to examine
c       jcol2  = end of range of columns to examine
c       kmax   = dimension of aamatr
c
c     Output:
c
c       qldep  = flag variable
c                  .true.  = linear dependence
c                  .false. = no linear dependence
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer kmax
c
      integer irow1,irow2,jcol1,jcol2
c
      logical qldep
c
      real*8 aamatr(kmax,kmax)
      real*8 eps100
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,idim,ipiv,irow,j,jcol,jdim,nrowi,nrowp
c
      real*8 aax,avmax,av,factor
c
c-----------------------------------------------------------------------
c
      idim = irow2 - irow1 + 1
      jdim = jcol2 - jcol1 + 1
c
      do i = irow1,irow2
        avmax = 0.
        do j = jcol1,jcol2
          av = abs(aamatr(i,j))
          if (av .gt. avmax) avmax = av
        enddo
        if (avmax .eq. 0.) then
c
c         Have a row containing all zeros.
c
          qldep = .true.
          go to 999
        endif
      enddo
c
      if (idim .gt. jdim) then
c
c       Have linear dependence because the number of rows is greater
c       than the number of columns.
c
        qldep = .true.
        go to 999
      endif
c
c     Note: nrowp is the number of rows processed; irow is the
c     index of the row currently being tested.
c
      nrowp = 0
c
      do jcol = jcol1,jcol2
        irow = irow1 + nrowp
        avmax = 0.
        ipiv = irow
c
        do i = irow,irow2
          av = abs(aamatr(i,jcol))
          if (av .gt. avmax) then
            avmax = av
            ipiv = i
          endif
        enddo
c
        if (avmax .gt. 0.) then
c
c         The present column is non-empty for the present and following
c         rows. This suffices to establish the existence of another
c         independent row.
c
          nrowp = nrowp + 1
c
          if (ipiv .ne. irow) then
c
c           Pivot- exchange rows irow and ipiv.
c
            do j = jcol,jcol2
              aax = aamatr(irow,j)
              aamatr(irow,j) = aamatr(ipiv,j)
              aamatr(ipiv,j) = aax
            enddo
          endif
c
c         Move down the present column from the current row to the last
c         row. If a non-zero entry is found, multiply the row by the
c         factor necessary to set this entry to unity. Note that
c         entries in any preceding columns of the submatrix are zero
c         at this point.
c
          do i = irow,irow2
            if (aamatr(i,jcol) .ne. 0.) then
              factor = 1./aamatr(i,jcol)
              do j = jcol,jcol2
                aamatr(i,j) = factor*aamatr(i,j)
              enddo
            endif
          enddo
c
          do i = irow + 1,irow2
            if (aamatr(i,jcol) .ne. 0.) then
c
c             Have a following row whose entry in the current column is
c             not zero (hence it has a value of unity). Make it zero by
c             subtracting the irow-th row (whose entry in the current
c             column is also unity) from the following row now under
c             consideration.
c
              do j = jcol,jcol2
                aamatr(i,j) = aamatr(i,j) - aamatr(irow,j)
                if (abs(aamatr(i,j)) .le. eps100) aamatr(i,j) = 0.
              enddo
            endif
          enddo
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Count the number of independent rows.
c
      nrowi = 0
      do i = irow1,irow2
        do j = jcol1,jcol2
          if (aamatr(i,j) .ne. 0.) then
            nrowi = nrowi + 1
            go to 100
          endif
        enddo
  100   continue
      enddo
c
      qldep = nrowi .lt. idim
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
