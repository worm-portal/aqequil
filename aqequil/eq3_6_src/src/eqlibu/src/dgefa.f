      subroutine dgefa(gmmatr,kmax,kdim,ipivot,info)
c
c     This subroutine factors the real*8 array gmmatr, using the method
c     of L-U decomposition. It is an adaptation of the 1979 Linpack
c     subroutine of the same name. The size and order of the original
c     calling sequence has been preserved:
c     (a,lda,n,ipvt,info) = (gmmatr,kmax,kdim,ipivot,info)
c
c     Note: detection of singularity using the test corresponding to
c     the condition info = -1 was not implemented in the original
c     Linpack version of the subroutine and may be missing from other
c     modern equivalents.
c
c     This subroutine uses the following pseudo-Linpack BLAS (Basic
c     Linear Algebra Subsystem) subroutines:
c
c       EQLIBU/idamax.f
c       EQLIBU/daxpy.f
c       EQLIBU/dscal.f
c
c     This subroutine is called by:
c
c       EQLIBU/msolvr.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       gmmatr = the matrix (unfactored)
c       kmax   = the declared (maximum) dimension of the matrix
c       kdim   = the used dimension of the matrix
c       ipivot = the pivot vector (used by EQLIBU/dgesl.f)
c
c     Output:
c
c       gmmatr = the matrix, factored
c       info   = error flag:
c                  = -1      singularity
c                  =  0      no singularity
c                  =  k > 0  singularity
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
      integer info,kdim
c
      integer ipivot(kmax)
c
      real*8 gmmatr(kmax,kmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j,k,kp1,il
      integer idamax
c
      real*8 adiv,div,tx
c
c-----------------------------------------------------------------------
c
c     Gaussian elimination with partial pivoting.
c
      info = 0
      if (kdim .ge. 2) then
        do k = 1,kdim - 1
        kp1 = k + 1
c
c       Find il = pivot index.
c
        il = idamax(kdim - k + 1,gmmatr(k,k),1) + k - 1
        ipivot(k) = il
c
c       Zero pivot implies this column already triangularized.
c
        if (gmmatr(il,k) .ne. 0.) then
c
c         Interchange if necessary.
c
          if (il .ne. k) then
            tx = gmmatr(il,k)
            gmmatr(il,k) = gmmatr(k,k)
            gmmatr(k,k) = tx
          endif
c
c         Compute multipliers.
c
          div = gmmatr(k,k)
c
c         Test for potential zero divide.
c
          adiv=abs(div)
          if (adiv.le.0.) then
            info = -1
            go to 999
          endif
          tx = -1./div
          call dscal(kdim - k,tx,gmmatr(k + 1,k),1)
c
c         Row elimination with column indexing.
c
          do j = kp1,kdim
            tx = gmmatr(il,j)
            if (il .ne. k) then
              gmmatr(il,j) = gmmatr(k,j)
              gmmatr(k,j) = tx
            endif
            call daxpy(kdim - k,tx,gmmatr(k + 1,k),1,gmmatr(k + 1,j),1)
          enddo
          go to 100
        endif
        info = k
  100   continue
      enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      ipivot(kdim) = kdim
      if (gmmatr(kdim,kdim) .eq. 0.) info = kdim
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
