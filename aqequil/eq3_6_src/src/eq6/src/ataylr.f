      subroutine ataylr(delxi,daffp0,nord,nordmx,npt,nptmax,affp0,affp)
c
c     This subroutine evaluates Taylor's series expansions for phase
c     affinities. These expansions are used to find phase boundaries
c     at which new phases appear in the ES.
c
c     Compare with:
c
c       EQ6/ptaylr.f
c       EQ6/rtaylr.f
c       EQ6/ztaylr.f
c
c     This subroutine is called by:
c
c       EQ6/fpbnpp.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       delxi  = the step size (in reaction progress)
c       daffp0 = the phase affinity derivative vector at the base point
c       nord   = the order of the truncated Taylor's series
c       nordmx = the maximum order of the truncated Taylor's series
c       npt    = the number of phases
c       nptmax = the maximum number of phases
c       affp0  = the phase affinity vector at the base point
c
c     Principal output:
c
c       affp   = the phase affinity vector at the new point
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nordmx,nptmax
c
      integer nord,npt
c
      real*8 affp(nptmax),affp0(nptmax),daffp0(nordmx,nptmax)
c
      real*8 delxi
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n,np
c
      real*8 ax,dxp
c
      real*8 fctrl
c
c-----------------------------------------------------------------------
c
      do np = 1,npt
        ax = affp0(np)
        dxp = 1.
        do n = 1,nord
          dxp = dxp*delxi
          ax = ax + ( daffp0(n,np)/fctrl(n) )*dxp
        enddo
        affp(np) = ax
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
