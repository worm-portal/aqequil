      subroutine undflw
c
c     This subroutine disables underflow trapping. This is not
c     necessary on most modern platform/OS/compiler combinations, as
c     they do not engage in underflow trapping by default at the
c     system level. More commonly, underflow trapping is now a compiler
c     option (which should be turned off in the case of EQ3/6).
c
c     On some (mainly older) systems, underflow trapping is built into
c     the system in a way that can only be overridden for a given
c     program by putting a platform-dependent system subroutine call
c     in that program. If such a call is necessary, it should be made
c     in the present subroutine. If it appears that such a call is
c     necessary, consult your local system documentation, compiler
c     documentation, compiler vendor, or local computer support person
c     to find out what this call should be. Do not call LLNL to find
c     out.
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
com   BEGIN_MACHINE_DEPENDENT_CODE
com
com     If a system subroutine call is necessary to disable underflow
com     trapping, add it here.
com
com
com   END_MACHINE_DEPENDENT_CODE
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
