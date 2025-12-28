      subroutine setrcp(aftarg,dlxmax,dlxmin,dlxmx0,npslmx,nsscmx,
     $ nsslmx,prcinf,sscrew,tolaft,tolsat,tolsst,zkfac,zklgmn,
     $ zklogl,zklogu)
c
c     This subroutine sets the values of various run control parameters
c     which are not read from the input file.
c
c     This subroutine is called by:
c
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c     Principal output:
c
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nsscmx
c
      integer npslmx,nsslmx
c
      real*8 sscrew(nsscmx)
c
      real*8 aftarg,dlxmax,dlxmin,dlxmx0,prcinf,tolaft,tolsat,tolsst,
     $ zkfac,zklgmn,zklogl,zklogu
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n
c
      real*8 tlg
c
c-----------------------------------------------------------------------
c
c     Supersaturation tolerance (if phase boundary searches are being
c     made and this is exceeded, the step size is cut). Note that tolsst
c     must be greater than tolsat (see aftarg and tolaft below).
c
      tolsst = 2.*tolsat
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Affinity target (the target for a search for a phase appearance
c     boundary).
c
      aftarg = 0.5*(tolsst + tolsat)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Affinity tolerance (the tolerance about the affinity target).
c     This requires the affinity at the located phase appearnce
c     boundary to lie between tolsat and tolsst.
c
      tolaft = tolsst - aftarg
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Setscrew parameters.
c
c       sscrew(1) = Step size control parameter, bounds the absolute
c                     value of the estimated error in any Taylor's
c                     series for a basis species variable. This is used
c                     to choose the order and the step size in normal
c                     reaction path computational mode. It is used to
c                     choose only the order in economy mode.
c
c       sscrew(2) = Not used.
c
c       sscrew(3) = Step size control parameter, bounds the absolute
c                     value of the estimated error in Taylor's
c                     series for rate functions (kinetic mode only).
c                     It also serves a function similar to that of
c                     sscrew(4) in testing the estimated error in the
c                     absolute time or the reaction progress of any
c                     individual irreversible reaction. If either of
c                     the sscrew(3) or sscrew(4) tests for corrector
c                     iteration or step size cuts is satisfied, no
c                     corrective action is taken.
c
c       sscrew(4) = Tolerance for corrector iteration or step size
c                     cut for differences between actual and predicted
c                     rate function values (kinetic mode only).
c
c       sscrew(5) = Maximum magnitude of Newton-Raphson correction
c                     term per iteration.
c
c       sscrew(6) = Step size parameter for economy mode; this bounds
c                     the change in a basis variable.
c
      do n = 1,nsscmx
        sscrew(n) = 0.
      enddo
c
      sscrew(1) = 1.e-4
c     if (sscrew(1) .lt. 1.e-6) sscrew(1) = 1.e-6
c     if (sscrew(1) .gt. 1.e-4) sscrew(1) = 1.e-4
c
      sscrew(3) = 1.e-4
c     if (sscrew(3) .lt. 1.e-6) sscrew(3) = 1.e-6
c     if (sscrew(3) .gt. 1.e-4) sscrew(3) = 1.e-4
c
      sscrew(4) = 1.e-6
c     if (sscrew(4) .lt. 1.e-8) sscrew(4) = 1.e-8
c     if (sscrew(4) .gt. 1.e-4) sscrew(4) = 1.e-4
c
      sscrew(5) = 4.0
c     if (sscrew(5) .lt. 1.0) sscrew(5) = 1.0
c     if (sscrew(5) .gt. 6.0) sscrew(5) = 6.0
c
      sscrew(6) = 4.0
c     if (sscrew(6) .lt. 0.5) sscrew(6) = 0.5
c     if (sscrew(6) .gt. 1.e+38) sscrew(6) = 1.e+38
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Z parameters.
c
c       zklogu = Threshhold and/or target value for the log number
c                  of moles of a phase; it is used in the following
c                  ways:
c
c                  (1) When the log number of moles of a phase is less
c                  than this value, EQ6 no longer tries to limit the
c                  step size to force the corresponding Taylor's
c                  series to meet the sscrew(1) accuracy criterion.
c                  try to keep the corresponding taylor's series
c                  accurate according to the sscrew(1) criterion.
c
c                  (2) This is the target value employed for the
c                  log number of moles when using Taylor's series
c                  to predict a phase disappearance boundary.
c
c                  (3) In the fluid-centered, flow-through open
c                  system model (iopt(1) = 2), this also determines
c                  (in conjunction with the zkfac parameter, see
c                  below) the maximum number of moles of a phase which
c                  is permitted to re-dissolve. It limits the number
c                  of moles which can be moved from the ES to the PRS
c                  (see zkfac below).
c
c                  Recommended values lie in the range -4.0 to -10.0.
c
c       zklogl = The amount by which the log number of moles of a phase
c                  is decremented by a transfer from the ES to the PRS
c                  in the fluid-centered, flow-through open system
c                  model (iopt(1) = 2). A value of 2.0 means that 99%
c                  of the number of moles is transferred, a value of
c                  3.0, 99.9%. Recommended that values lie in the range
c                  of 2.0 to 4.0.
c
c       zkfac  = Parameter which determines the minimum number of moles
c                  of a phase left in the ES after a shift to the PRS
c                  in the fluid-centered, flow-through open system
c                  model (iopt(1) = 2). This minimum value is calculated
c                  as zkfac*10**zklogu (equivalent to 10**zklgmn, see
c                  below). Recommended values lie in the range 0.9
c                  to 0.98.
c
c       zklgmn = Secondary parameter equivalent to the minimum log
c                  number of moles of a phase left in the ES after a
c                  shift to the PRS in the fluid-centered, flow-through
c                  open system.
c
      zklogu = -7.0
      zklogl = 2.0
      zkfac = 0.98
      zklgmn = tlg(zkfac) + zklogu
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The minimum step size. No mechanism for cutting the step size will
c     reduce it to less than this value. Typically, the minimum step
c     size is less than the zero-order step size. It is undesirable to
c     let the minimum step size get too big, as this determines the
c     resolution of certain events of interest (such as hitting a
c     desired pH value at the end of the run). The resolution of events
c     in reaction progress space and time space is not limited by this
c     parameter.
c
      dlxmin = 0.01*dlxmx0
      dlxmin = min(dlxmin,1.e-14)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Upper limit on the step size.
c
      dlxmax = prcinf
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum number of attempts to slide over a phase boundary.
c
      npslmx = 8
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum number of attempts to slide over a redox jump.
c
      nsslmx = 8
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
