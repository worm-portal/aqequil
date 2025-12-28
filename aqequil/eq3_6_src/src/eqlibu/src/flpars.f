      subroutine flpars(eps100,irang,noutpt,nttyo,smp100)
c
c     This subroutine obtains the following real*8 parameters: the
c     real*8 machine epsilon (eps), the smallest positive real*8
c     number (smpos), and the exponent range (irang). It returns
c     eps100 = 100*eps and smp100 = 100*smpos, not of eps and smpos.
c     This subroutine calls EQLIBU/chump.f to test the adequacy of the
c     machine epsilon and the exponent range. If inadequacies are
c     noted, warnings will be written to the output and screen files.
c
c     This subroutine is called by:
c
c       EQPT/eqpt.f
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       noutpt = unit number of the output file
c       nttyo  = unit number of the screen file
c
c     Output:
c
c       eps100 = 100*eps
c       smp100 = 100*smpos
c       irang  = the real*8 exponent range, int( -log10(smpos ))
c
c     Local:
c
c       eps    = the real*8 machine epsilon
c       smpos  = the smallest positive real*8 number
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer noutpt,nttyo
c
      integer irang
c
      real(8) eps100,smp100
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      real(8) eps,smpos,x
c
c-----------------------------------------------------------------------
c
      data x / 1.0 /
c
c-----------------------------------------------------------------------
c
c     Get eps (real*8 machine epsilon).
c
      eps = epsilon(x)
c
c     Get smpos(real*8 smallest positive number).
c
      smpos = tiny(x)
c
      eps100 = 100*eps
      smp100 = 100*smpos
c
c     Compute the exponent range.
c
      irang = int(-log10(smpos))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Clear the IEEE flag for floating-point underflow, if such a
c     flag is present, to avoid getting an unnecessary system
c     warning message. Underflow is a normal condition in EQ3/6.
c     Make any porting changes in the EQLIBU subroutine cliefu.
c     Do not make the porting changes here.
c
c     Note: the calculation of smpos above involves forcing an
c     underflow.
c
      call cliefu()
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Test the adequacy of the machine epsilon and the exponent range.
c
      call chump(eps,irang,noutpt,nttyo)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
