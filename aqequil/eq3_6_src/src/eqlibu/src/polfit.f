      subroutine polfit(aamatr,cof,gmmatr,ier,ipivot,npfmax,npft,
     $ noutpt,nttyo,xvec,yvec)
c
c     This subroutine fits an exact polynomial through npft points.
c     Each point is an x,y pair, where the values of x are in the
c     xvec array, those of y in the yvec array. The method is to solve
c     the matrix equation (A)(c) = (y) where:
c
c               | 1  x(1)  x(1)**2  ... |
c               | 1  x(2)  x(2)**2  ... |
c       A =     | 1  x(3)  x(3)**2  ... |
c               | ...                   |
c               | 1  x(n)  x(n)**2  ... |
c
c      and solve for the coefficients (c), which give the
c      representation:
c
c        y = c(1) + c(2)*x + c(3)*x**2 + ... + c(n)*x**n-1
c
c     It is assumed that the xvec array contains distinct values
c     (i.e., there are no duplicates). The matrix A is the
c     array aamatr, the matrix dimension n is npft, and the
c     vector c is the array cof.
c
c     If the xvec array was scaled using EQLIBU/scalx1.f, the resulting
c     yvec array should be rescaled using EQLIBU/rscaly.f.
c
c     This subroutine is called by:
c
c       EQPT/intrp.f
c       Any
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       xvec   = the x vector
c       yvec   = the y vector
c       npft   = number of x,y points
c
c     Output:
c
c       cof    = the array of coefficients fo the fitted polynomial
c       ier    = error flag (returned by EQLIBU/msolvr.f):
c                  =  0  Okay
c                  =  1  Encountered a zero matrix
c                  =  2  Encountered a non-zero, computationally
c                          singular matrix
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer npfmax
c
      integer ipivot(npfmax)
      integer ier,noutpt,npft,nttyo
c
      real*8 aamatr(npfmax,npfmax),gmmatr(npfmax,npfmax)
      real*8 cof(npfmax),xvec(npfmax),yvec(npfmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,j
c
      logical qpr
c
      real*8 aax
c
c-----------------------------------------------------------------------
c
      data qpr    /.false./
c
c-----------------------------------------------------------------------
c
c     Initialize error flag.
c
      ier = 0
c
c     Set up a matrix equation.
c
      do j = 1,npft
        cof(j) = 0.
        aax = 1.
        aamatr(j,1) = 1.
        do i = 2,npft
          aax = aax*xvec(j)
          aamatr(j,i) = aax
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Solve the matrix equation.
c
c     Calling sequence substitutions:
c       cof for delvec
c       npft for kdim
c       npfmax for kmax
c       yvec for rhsvec
c
      call msolvr(aamatr,cof,gmmatr,ier,ipivot,npft,npfmax,
     $ noutpt,nttyo,qpr,yvec)
c
      end
