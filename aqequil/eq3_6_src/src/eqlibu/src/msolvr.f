      subroutine msolvr(aamatr,delvec,gmmatr,ier,ipivot,kdim,kmax,
     $ noutpt,nttyo,qpr,rhsvec)
c
c     This subroutine solves the matrix equation:
c
c       (aamatr)(delvec) = (rhsvec)
c
c     The method employed is L-U decomposition. The original aamatr
c     matrix is preserved. The workspace vector gmmatr holds the
c     L-U decomposition.
c
c     This subroutine uses
c     the following subroutines, which are adaptations of 1979 Linpack
c     subroutines of the same names:
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
c       EQLIBU/polfit.f
c       EQLIBU/minvrt.f
c       EQLIB/nrstep.f
c       EQ3NR/arrsim.f
c       EQ6/jgibbs.f
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
c       rhsvec = the right hand side vector
c       qpr    = logical switch to print messages from this subroutine
c
c     Workspace:
c
c       gmmatr = work space for a copy of the matrix; holds the
c                  L-U decomposition
c       ipivot = the pivot vector
c
c     Output:
c
c       delvec = the solution vector
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
      real*8 aamatr(kmax,kmax),gmmatr(kmax,kmax),rhsvec(kmax),
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
c     Initialize error flag.
c
      ier = 0
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
 1000   format(/' * Note - (EQLIBU/msolvr) Have a zero matrix to',
     $  ' factor.')
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
 1010     format(/' * Note - (EQLIBU/msolvr) Have a computationally',
     $    ' singular matrix.')
        endif
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c      Copy rhsvec to delvec.
c
      do j = 1,kdim
        delvec(j) = rhsvec(j)
      enddo
c
c     OR
c
c     call copyaa(rhsvec,delvec,kdim)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Solve for the unknown vector (delvec). Note that delvec initially
c     contains a copy of the right-hand-side vector.
c
      call dgesl(gmmatr,kmax,kdim,ipivot,delvec)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
