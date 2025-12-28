      subroutine gwhcfa(akmat1,delxi,dxsm11,hhcvec,nord,nordmx,
     $ nrd1mx,whcfac,xhcvec)
c
c     This subroutine computes the w factor (whcfac), which is
c     required for higher-order (stiff) ODE corrections. Logically,
c     w = g dot FCh where g and h are vectors, and (FC) is the matrix
c     required to convert finite differences utilizing the new point
c     of reaction progress to derivatives at the base point. This
c     factor depends only on the recent step size history, including
c     the current step size (delxi).
c
c     This subroutine is called by:
c
c       EQ6/path.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       akmat1 = the (FC) matrix
c       delxi  = the current step size
c       dxsm11 = the array of cumulative step sizes going back from
c                  new point, centered on the new point
c       nord   = the order of the finite difference method
c       nordmx = the maximum order of the finite difference method
c       nrd1mx = nordmx + 1
c       hhcvec = the h vector, here effectively work space
c       xhcvec = the FCh vector, here effectively work space
c
c     Principal output:
c
c       whcfac = the w factor
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nordmx,nrd1mx
c
      integer nord
c
      real*8 akmat1(nrd1mx,nrd1mx),dxsm11(nrd1mx),hhcvec(nordmx),
     $ xhcvec(nordmx)
c
      real*8 delxi,whcfac
c
c----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,j
c
      real(8) gx
c
c----------------------------------------------------------------------
c
c     First calculate h (hhcvec).
c
      hhcvec(1) = 1./delxi
      do j = 2,nord
        hhcvec(j) = hhcvec(j - 1)/dxsm11(j)
      enddo
c
c     Now calculate the (FCh) vector (xhcvec).
c
      do i = 1,nord
        xhcvec(i) = 0.
        do j = 1,nord
          xhcvec(i) = xhcvec(i) + akmat1(i,j)*hhcvec(j)
        enddo
      enddo
c
c     Finally, emulate w = g dot (FCh).
c
      whcfac = 0.
      gx = delxi
      do j = 1,nord
        gx = gx*delxi/(j + 1)
        whcfac = whcfac + gx*xhcvec(j)
      enddo
c
      end
