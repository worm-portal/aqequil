      subroutine sfncge(delxi,xval0,xtargv,dxval0,nord,nrd1mx,resx)
c
c     This subroutine computes a general search function for
c     EQ6/search.f. This search function is a residual defined as the
c     difference between the calculated value of a function (xvalc)
c     described by a truncated Taylor's series and a target value
c     (xtargv).
c
c     This subroutine is called by:
c
c       EQ6/search.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       delxi  = step size
c       xtargv = target value
c       xval0  = the function value at delxi = 0
c       dxval0 = the derivatives of this function at delxi = 0
c       nord   = the order of the truncated Taylor's series
c
c     Principal output:
c
c       resx   = value of residual function
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nrd1mx
c
      integer nord
c
      real*8 dxval0(nrd1mx)
c
      real*8 delxi,resx,xtargv,xval0
c
      real*8 fctrl
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n
c
      real*8 dxp,xvalc
c
c-----------------------------------------------------------------------
c
      xvalc = xval0
      if (nord .gt. 0) then
        dxp = 1.
        do n = 1,nord
          dxp = dxp*delxi
          xvalc = xvalc + ( dxval0(n)/fctrl(n) )*dxp
        enddo
      endif
      resx = xvalc - xtargv
c
      end
