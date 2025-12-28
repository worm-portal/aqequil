      subroutine nbsgam(acfnbs,adh,fxi,nchlor,noutpt,nttyo)
c
c     This subroutine calculates the log activity coefficient of the
c     chloride ion according to the NBS pH convention (e.g., Covington,
c     Bates, and Durst, 1985). The convention itself may be extended
c     outside the specified limit of 0.1 molal on the ionic strength.
c
c     Reference:
c
c       Covington, A.K., Bates, R.G., and Durst, R.A., 1985,
c       Definition of pH scales, standard reference values, measure-
c       ment of pH and related terminology (recommendations, 1984):
c       Pure and Applied Chemistry, v. 57, p. 533-542.
c
c     This subroutine is called by:
c
c       EQLIB/gpheh.f
c       EQLIBG/gcoeff.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       adh    = Debye-Huckel A(gamma) parameter
c       fxi    = the ionic strength (the 2nd-order electrostatic
c                  moment function I)
c       nchlor = the species index of the chloride ion
c
c     Principal output:
c
c       acfnbs = log gamma(Cl-), according to the Bates-Guggenheim
c                convention
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nchlor,noutpt,nttyo
c
      real*8 acfnbs,adh,fxi
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      real*8 fxisqt
c
c-----------------------------------------------------------------------
c
c     Test the species index of the chloride ion.
c
      if (nchlor .le. 0) then
        write (noutpt,1000)
        write (nttyo ,1000)
 1000   format(/' * Error - (EQLIBG/nbsgam) Have no index for the',
     $  /7x,"chloride ion. Can't use the extended NBS pH scale.")
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Evaluate the NBS expression for the molal activity coefficient
c     of the chloride ion.
c
      fxisqt = sqrt(fxi)
      acfnbs = - ( adh * fxisqt ) / ( 1.0 + ( 1.5 * fxisqt ) )
c
      end
