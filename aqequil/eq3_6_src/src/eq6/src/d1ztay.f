      subroutine d1ztay(delxi,dzvc0,d1zvc1,kdim,kmax,nord,nrd1mx)
c
c     This subroutine computes the Taylor's series expansions for the
c     first derivatives of the master variables. These expansions are
c     used to find points of reaction progress at which such variables
c     maximize. At such points, for example in the fluid-centered
c     flow-through open system model,it may be necessary to transfer
c     mass of a corresponding phase from the ES to the PRS.
c
c     See also:
c
c       EQ6/d2ztay.f
c       EQ6/ztaylr.f
c
c     This subroutine is called by:
c
c       EQ6/eqcalc.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       delxi  = the step size (in reaction progress)
c       dzvc0  = the dz/d(xi) vector at the base point
c       kdim   = the number of elements in the z vector
c       kmax   = the maximum number of elements in the z vector
c       nord   = the order of the truncated Taylor's series
c       nrd1mx = the maximum order of the truncated Taylor's series + 1
c
c     Principal output:
c
c       d1zvc1 = the dz/d(xi) vector at the new point
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer kmax,nrd1mx
c
      integer kdim,nord
c
      real*8 dzvc0(nrd1mx,kmax),d1zvc1(kmax)
c
      real*8 delxi
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j,jk,kcol
c
      real*8 dxp,d1z
c
      real*8 fctrl
c
c-----------------------------------------------------------------------
c
      if (nord .gt. 0) then
        do kcol = 1,kdim
          d1z = dzvc0(1,kcol)
          if (nord .gt. 1) then
            jk = nord - 1
            dxp = 1.
            do j = 1,jk
              dxp = dxp*delxi
              d1z = d1z + ( dzvc0(j + 1,kcol)/fctrl(j) )*dxp
            enddo
          endif
          d1zvc1(kcol) = d1z
        enddo
      endif
c
      end
