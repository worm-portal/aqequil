      subroutine d2ztay(delxi,dzvc0,d2zvc1,kcol,kmax,nord,nrd1mx)
c
c     This subroutine computes the Taylor's series expansion for the
c     second derivative of the kcol-th master variables. This second
c     derivative is used to test whether or not a critical point
c     corresponds to a maximum.
c
c     See also:
c
c       EQ6/d1ztay.f
c       EQ6/ztaylr.f
c
c     This subroutine is called by:
c
c       None
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       delxi  = the step size (in reaction progress)
c       dzvc0  = the dz/d(xi) vector at the base point
c       kcol   = the index of the master variable whose second
c                  derivative is to be calculated
c       kmax   = the maximum number of elements in the z vector
c       nord   = the order of the truncated Taylor's series
c       nrd1mx = the maximum order of the truncated Taylor's series + 1
c
c     Principal output:
c
c       d2zvc1 = the d2z/d(xi)2 element for the kcol-th master
c                  variable
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
      integer kcol,nord
c
      real*8 dzvc0(nrd1mx,kmax)
c
      real*8 delxi,d2zvc1
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j,jk
c
      real*8 dxp,d2z
c
      real*8 fctrl
c
c-----------------------------------------------------------------------
c
      d2zvc1 = 0.
      if (nord .gt. 1) then
        d2z = dzvc0(2,kcol)
        if (nord .gt. 2) then
          jk = nord - 2
          dxp = 1.
          do j = 1,jk
            dxp = dxp*delxi
            d2z = d2z + ( dzvc0(j + 2,kcol )/fctrl(j) )*dxp
          enddo
        endif
        d2zvc1 = d2z
      endif
c
      end
