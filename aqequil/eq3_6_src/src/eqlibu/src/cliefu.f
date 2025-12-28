      subroutine cliefu()
c
c     This subroutine clears the IEEE flag for floating-point underflow,
c     if such a flag is present, to avoid getting an unnecessary system
c     warning message. Underflow is a normal condition in EQ3/6.
c
c     The presence of an IEEE flag for underflow is platform and
c     compiler dependent. If the flag is present, the need to clear
c     it may also depend on the platform and compiler. Currently,
c     the flag is present and needs to be cleared on Sun SPARCstations
c     in the case of Sun's Fortran 77 compiler. This is not the case
c     with any Fortran 90 compiler that has been seen to date,
c     including Sun's.
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
c       None
c
c     Output:
c
c       None
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
c       None
c
c-----------------------------------------------------------------------
c
com   BEGIN_SPARC_DEPENDENT_CODE
com     BEGIN_F77_DEPENDENT_CODE
com
com       Here cieeef() is a C subroutine that is distributed with
com       EQ3/6. The object file is usually attached to the library
com       (.a) file for the EQLIBU library, for convenience in linking.
com
c         call cieeef()
c         go to 999
com
com     END_F77_DEPENDENT_CODE
com   END_SPARC_DEPENDENT_CODE
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
