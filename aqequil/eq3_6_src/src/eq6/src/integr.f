      subroutine integr(delxi,dlxrct,drer0,nord,nrc,nrctmx,
     $ nrd1mx,rrelr0)
c
c     This subroutine integrates a Taylor's series for the relative
c     rate of a reaction to calculate the advancement in the
c     corresponding reaction progress variable (dlxrct).
c
c     This subroutine is called by:
c
c       EQ6/reacts.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       delxi  = overall reaction progress variable
c       drer0  = array of derivatives of relative rates for individual
c                  irreversible reactions
c       nord   = order of the finite difference approximation
c
c     Principal output:
c
c       dlxrct = array of increments of reaction progress for
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
      integer nrctmx,nrd1mx
c
      integer nord,nrc
c
      real*8 drer0(nrd1mx,nrctmx),rrelr0(nrctmx)
      real*8 delxi,dlxrct
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n
c
      real*8 dxp
c
      real*8 fctrl
c
c-----------------------------------------------------------------------
c
      dlxrct = rrelr0(nrc)*delxi
      if (nord .gt. 0) then
        dxp = delxi
        do n = 1,nord
          dxp = dxp*delxi
          dlxrct = dlxrct + ( drer0(n,nrc)/fctrl(n + 1) )*dxp
        enddo
      endif
c
      end
