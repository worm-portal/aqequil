      subroutine d2ptay(delxi,demop0,d2emp1,nord,nordmx,npet,npetmx)
c
c     This subroutine computes the Taylor's series expansion for the
c     second derivative of the number of moles of the npe-th phase
c     in the ES. This second derivative is used to test whether or not
c     a critical point corresponds to a maximum.
c
c     See also:
c
c       EQ6/d1ptay.f
c       EQ6/ptaylr.f
c
c     This subroutine is called by:
c
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
c       npe    = the ES index of the phase for which the second
c                  derivative of the number of moles is to be computed
c       npetmx = the maximum number of phases in the ES
c
c     Principal output:
c
c       d2emp1 = the second derivative of the number of moles vector
c                  for phases in the ES at the new point
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
      real*8 demop0(nordmx,npetmx),d2emp1(npetmx)
c
      real*8 delxi
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j,jk,npe
c
      real*8 dxp,d2p
c
      real*8 fctrl
c
c-----------------------------------------------------------------------
c
      do npe = 1,npet
        d2emp1(npe) = 0.
        if (nord .gt. 1) then
          d2p = demop0(2,npe)
          if (nord .gt. 2) then
            jk = nord - 2
            dxp = 1.
            do j = 1,jk
              dxp = dxp*delxi
              d2p = d2p + ( demop0(j + 2,npe )/fctrl(j) )*dxp
            enddo
          endif
          d2emp1(npe) = d2p
        endif
      enddo
c
      end
