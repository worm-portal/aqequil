      subroutine minvrt(aamatr,aimatr,delvec,gmmatr,ier,ipivot,
     $ kdim,kmax,noutpt,nttyo,qpr)
c
c     This subroutine inverts the matrix aamatr. The inverted matrix
c     is returned in aimatr. Thus
c
c       [aamatr][aimatr] = I
c
c     where I is the unit matrix (1's on the diagonal, 0's elsewhere).
c
c     The method is L-U decomposition of the matrix, in combination
c     with the subsequent solution of a series of matrix equations,
c     each using the same decomposed matrix but with a different unit
c     vector as the right hand side. The original aamatr matrix is
c     preserved. The workspace vector gmmatr holds the L-U
c     decomposition.
c
c     A companion subroutine ginvrt.f has the same function as the
c     present routine. However, it assumes that aamatr has already
c     been factored into gmmatr.
c
c     This subroutine uses the following subroutines, which are
c     adaptations of 1979 Linpack subroutines of the same names:
c
c       EQLIBU/dgefa.f
c       EQLIBU/dgesl.f
c
c     The calling sequences have been preserved to facilitate
c     substituting equivalents optimized for specific platforms,
c     for anyone desiring to do so.
c
c     This subroutine is called by:
c
c       EQ6/path.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       aamatr = the matrix
c       kdim   = the actual dimention of the matrix
c       kmax   = the maximum dimension of the matrix
c       noutpt = the unit number of the output file
c       nttyo  = the unit number of the screen file
c       qpr    = logical switch to print messages from this subroutine
c
c     Workspace:
c
c       delvec = at one point, the right-hand-side vector; at another,
c                  the solution vector
c       gmmatr = work space for a copy of the matrix; holds the
c                  L-U decomposition
c       ipivot = the pivot vector
c
c     Output:
c
c       aimatr = the inverted matrix
c       ier    = error flag:
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
      integer kmax
c
      integer noutpt,nttyo
c
      integer ipivot(kmax)
      integer ier,kdim
c
      logical qpr
c
      real*8 aamatr(kmax,kmax),aimatr(kmax,kmax),gmmatr(kmax,kmax),
     $ delvec(kmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer info,j,k
c
c-----------------------------------------------------------------------
c
c     Copy aamatr to gmmatr.
c
      do k = 1,kdim
        do j = 1,kdim
          gmmatr(j,k) = aamatr(j,k)
        enddo
      enddo
c
c     OR
c
c     do k = 1,kdim
c       call copyaa(aamatr(1,k),gmmatr(1,k),kdim)
c     enddo
c
c     OR
c
c     integer nmax
c     nmax = kmax*kdim
c     call copyaa(aamatr,gmmatr,nmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Avoid trying to decompose a zero matrix.
c
      do k = 1,kdim
        do j = 1,kdim
          if (gmmatr(j,k) .ne. 0.) go to 100
        enddo
      enddo
      ier = 1
      if (qpr) then
        write (noutpt,1000)
        write (nttyo,1000)
 1000   format(/' * Note - (EQLIBU/minvrt) Have a zero matrix to',
     $  ' invert.')
      endif
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Factor the matrix (make the L-U decomposition).
c
  100 call dgefa(gmmatr,kmax,kdim,ipivot,info)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Test info.
c
      if (info .ne. 0) then
        ier = 2
        if (qpr) then
          write (noutpt,1010)
          write (nttyo,1010)
 1010     format(/' * Note - (EQLIBU/minvrt) Have a computationally',
     $    ' singular matrix.')
        endif
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Complete the inversion. Loop over columns for the inverted
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
