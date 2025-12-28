      subroutine d1ptay(delxi,demop0,d1emp1,nord,nordmx,npet,npetmx)
c
c     This subroutine evaluates the first derivative of Taylor's series
c     expansions for the number of moles of phases in the ES. These
c     expansions are used to find points of reaction progress at which
c     such variables maximize. At such points, for example in the
c     fluid-centered flow-through open system model,it may be necessary
c     to transfer mass of a corresponding phase from the ES to the PRS.
c
c     See also:
c
c       EQ6/d2ptay.f
c       EQ6/ztaylr.f
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
c       npet   = the number of phases in the ES
c       npetmx = the maximum number of phases in the ES
c
c     Principal output:
c
c       d1emp1 = the first derivative of the number of moles vector
c                  for phases in the ES at at the new point
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
      real*8 demop0(nordmx,npetmx),d1emp1(npetmx)
c
      real*8 delxi
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j,jk,npe
c
      real*8 dxp,d1p
c
      real*8 fctrl
c
c-----------------------------------------------------------------------
c
      if (nord .gt. 0) then
        do npe = 1,npet
          d1p = demop0(1,npe)
          if (nord .gt. 1) then
            jk = nord - 1
            dxp = 1.
            do j = 1,jk
              dxp = dxp*delxi
              d1p = d1p + ( demop0(j + 1,npe)/fctrl(j) )*dxp
            enddo
          endif
          d1emp1(npe) = d1p
        enddo
      endif
c
      end
