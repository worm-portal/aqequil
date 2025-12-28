      subroutine ptaylr(delxi,demop0,emop0,emop,nord,nordmx,npet,npetmx)
c
c     This subroutine evaluates Taylor's series expansions for the
c     number of moles of phases in the ES. These expansions are used
c     to find phase boundaries at which phases disappear from the ES.
c
c     Compare with:
c
c       EQ6/ataylr.f
c       EQ6/rtaylr.f
c       EQ6/ztaylr.f
c
c     See also:
c
c       EQ6/d1ptay.f
c       EQ6/d2ptay.f
c
c     This subroutine is called by:
c
c       EQ6/fpbdpp.f
c       EQ6/fpbflo.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       delxi  = the step size (in reaction progress)
c       demop0 = the number of moles derivative vector for phases in
c                  the ES at the base point
c       nord   = the order of the truncated Taylor's series
c       nordmx = the maximum order of the truncated Taylor's series
c       npet   = the number of phases in the ES
c       npetmx = the maximum number of phases in the ES
c       emop0  = the number of moles vector for phases in the ES at
c                  at the base point
c
c     Principal output:
c
c       emop   = the number of moles vector for phases in the ES at
c                  at the new point
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nordmx,npetmx
c
      integer nord,npet
c
      real*8 emop(npetmx),emop0(npetmx),demop0(nordmx,npetmx)
c
      real*8 delxi
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n,npe
c
      real*8 dxp,ex
c
      real*8 fctrl
c
c-----------------------------------------------------------------------
c
      do npe = 1,npet
        ex = emop0(npe)
        dxp = 1.
        do n = 1,nord
          dxp = dxp*delxi
          ex = ex + ( demop0(n,npe)/fctrl(n) )*dxp
        enddo
        emop(npe) = ex
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
