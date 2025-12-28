      subroutine ginvrt(aimatr,delvec,gmmatr,ipivot,kdim,kmax)
c
c     This subroutine inverts the matrix aamatr. The inverted matrix
c     is returned in aimatr. Thus
c
c       [aamatr][aimatr] = I
c
c     where I is the unit matrix (1's on the diagonal, 0's elsewhere).
c
c     This subroutine assumes that aamatr has already been factored
c     into its L-U decomposition, gmmatr. The inverse matrix aimatr
c     is obtained by solving a series of matrix equations, each using
c     the same decomposed matrix but with a different unit vector as
c     the right hand side. The original aamatr matrix itself is not
c     used here.
c
c     A companion subroutine minvrt.f has the same function as the
c     present routine. However, it does its own factoring of aamatr
c     to obtain gmmatr.
c
c     This subroutine uses the following subroutines, which are
c     adaptations of 1979 Linpack subroutines of the same names:
c
c       EQLIBU/dgesl.f
c
c     The calling sequences have been preserved to facilitate
c     substituting equivalents optimized for specific platforms,
c     for anyone desiring to do so.
c
c     This subroutine is called by:
c
c       EQ6/garmat.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       gmmatr = the L-U decomposition of the matrix to be inverted
c       kdim   = the actual dimention of the matrix
c       kmax   = the maximum dimension of the matrix
c
c     Workspace:
c
c       delvec = at one point, the right-hand-side vector; at another,
c                  the solution vector
c       ipivot = the pivot vector
c
c     Output:
c
c       aimatr = the inverted matrix
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
      integer ipivot(kmax)
      integer kdim
c
      real*8 aimatr(kmax,kmax),gmmatr(kmax,kmax),delvec(kmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j,k
c
c-----------------------------------------------------------------------
c
c     Invert the matrix aamatr by looping over columns for the inverted
c     matrix (aimatr).
c
      do k = 1,kdim
c
c       Put the corresponding unit vector in delvec. This serves to
c       set the right-hand-side vector.
c
        do j = 1,kdim
          delvec(j) = 0.
        enddo
        delvec(k) = 1.
c
c       Solve for the unknown vector (delvec). Note that delvec
c       initially contains the right-hand-side vector. Normally
c       separate arrays would be used for the two vectors. This
c       example of bad programming is made here to retain consistency
c       with the original Linpack routine.
c
        call dgesl(gmmatr,kmax,kdim,ipivot,delvec)
c
c       Load the result in the inverted matrix.
c
        do j = 1,kdim
          aimatr(j,k) = delvec(j)
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
