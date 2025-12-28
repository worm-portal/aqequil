      subroutine dgesl(gmmatr,kmax,kdim,ipivot,delvec)
c
c     This subroutine solves the real*8 system (gmmatr)*(x) = (b),
c     where the matrix gmmatr has previously been factored by
c     EQLIBU/dgefa.f. It is expected that singularity of the matrix
c     has been trapped by that subroutine. Initially, the array delvec
c     contains the right-hand-side vector (b); this is replaced by
c     the solution vector (x). This subroutine is an adaptation of the
c     1979 Linpack subroutine of the same name. The size and order in
c     the original calling sequence has been preserved:
c     (a,lda,n,ipvt,b) = (gmmatr,kmax,kdim,ipivot,delvec)
c
c     This subroutine uses the following pseudo-Linpack BLAS (Basic
c     Linear Algebra Subsystem) subroutines:
c
c       EQLIBU/daxpy.f
c
c     This subroutine is called by:
c
c       EQLIBU/msolvr.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       gmmatr = the matrix
c       kmax   = the declared (maximum) dimension of the matrix
c       kdim   = the used dimension of the matrix
c       ipivot = the pivot vector (from EQLIBU/dgefa.f)
c       delvec = initially, the right hand side vector (b)
c
c     Output:
c
c       delvec = the solution vector (x).
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
      integer kdim
c
      integer ipivot(kmax)
c
      real*8 gmmatr(kmax,kmax),delvec(kmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer il,k,kb
c
      real*8 tx
c
c----------------------------------------------------------------------
c
c     First solve (L)*(y) = (b).
c
      if (kdim .ge. 2) then
        do k = 1, kdim - 1
          il = ipivot(k)
          tx = delvec(il)
          if (il .ne. k) then
            delvec(il) = delvec(k)
            delvec(k) = tx
          endif
          call daxpy(kdim - k,tx,gmmatr(k + 1,k),1,delvec(k + 1),1)
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Now solve  (U)*(x) = (y).
c
      do kb = 1, kdim
        k = kdim + 1 - kb
        delvec(k) = delvec(k)/gmmatr(k,k)
        tx = -delvec(k)
        call daxpy(k - 1,tx,gmmatr(1,k),1,delvec(1),1)
      enddo
c
      end
