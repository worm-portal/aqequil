      subroutine rtaylr(delxi,drer0,drir0,jreac,nord,nrct,nrctmx,
     $ nrd1mx,rirec0,rirecp,rrelr0,rrelrp)
c
c     This subroutine evaluates Taylor's series expansions for the
c     inverse rate and the relative rates of all irreversible reactions.
c
c     Compare with:
c
c       EQ6/ataylr.f
c       EQ6/ptaylr.f
c       EQ6/ztaylr.f
c
c     This subroutine is called by:
c
c       EQ6/path.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       delxi  = the step size (in reaction progress)
c       drer0  = the relative rate derivative vector at the base point
c       drir0  = the inverse rate derivative vector at the base point
c       nord   = the order of the truncated Taylor's series
c       nrd1mx = the maximum order of the truncated Taylor's series + 1
c       nrct   = the number of reactants
c       nrctmx = the maximum number of reactants
c       rirec0 = the inverse rate at the base point
c       rrelr0 = the relative rate vector at the base point
c
c     Principal output:
c
c       rirecp = the inverse rate at the new point
c       rrelrp = the relative rate vector at the new point
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nrctmx,nrd1mx
c
      integer jreac(nrctmx)
c
      integer nord,nrct
c
      real*8 drer0(nrd1mx,nrctmx),drir0(nrd1mx),rrelr0(nrctmx),
     $ rrelrp(nrctmx)
c
      real*8 delxi,rirecp,rirec0
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n,nrc
c
      real*8 dxp,rx
c
      real*8 fctrl
c
c-----------------------------------------------------------------------
c
      rirecp = rirec0
      dxp = 1.
      do n = 1,nord
        dxp = dxp*delxi
        rirecp = rirecp + ( drir0(n)/fctrl(n) )*dxp
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do nrc = 1,nrct
        rx = 0.
        if (jreac(nrc) .eq. 0) then
          rx = rrelr0(nrc)
          dxp = 1.
          do n = 1,nord
            dxp = dxp*delxi
            rx = rx + ( drer0(n,nrc)/fctrl(n) )*dxp
          enddo
        endif
        rrelrp(nrc) = rx
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
