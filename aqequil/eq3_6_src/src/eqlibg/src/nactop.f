      subroutine nactop(iopg,nopgmx,noutpt,nttyo,uactop)
c
c     This subroutine sets the name of the option ("uactop") for
c     computing the activity coefficients of aqueous species. It
c     also sets associated logical flags concerning the generic type
c     of activity coefficient option. The variable "iopg(1)" determines
c     the exact activity coefficient model.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       iopg   = array of activity coefficient option switches
c
c     Principal output:
c
c       uactop = string describing the activity coefficient model
c                  corresponding to iopg(1)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nopgmx
c
      integer iopg(nopgmx)
      integer noutpt,nttyo
c
      character*32 uactop
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
c       None
c
c-----------------------------------------------------------------------
c
      if (iopg(1) .eq. -1) then
c
c       The Davies equation.
c
        uactop = 'the Davies equation'
      elseif (iopg(1) .eq. 0) then
c
c       The B-dot equation.
c
        uactop = 'the B-dot equation'
      elseif (iopg(1) .eq. 1) then
c
c       Pitzer's equations.
c
        uactop = "Pitzer's equations"
      elseif (iopg(1) .eq. 2) then
c
c       The HC + DHC equations.
c
        uactop = 'the HC + DH equations'
      else
c
c       Have an error.
c
        write (noutpt,1000) iopg(1)
        write (nttyo,1000) iopg(1)
 1000   format(/' * Error - (EQLIBG/nactop) Have iopg(1) = ',i3,'.',
     $  /7x,'This does not correspond to a valid option for computing',
     $  /7x,'the activity coefficients of aqueous species.')
          stop
      endif
c
      end
