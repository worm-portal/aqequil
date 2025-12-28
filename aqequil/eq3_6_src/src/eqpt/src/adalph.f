      subroutine adalph(alpha,ipbtmx,iz1,iz2,npx2,npx2mx)
c
c     This subroutine assigns standard values of the Pitzer
c     alpha parameters for the npx2-th species pair. The standard
c     values depend on the charge combination.
c
c     This subroutine is called by:
c
c       EQPT/rdpca.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ndat0s = unit number of the stripped DATA0 file
c
c     Principal output:
c
c       nazt   = the number of specified hard core diameters
c       uazp   = array of aqueous species names used to specify
c                  hard core diamters on the data file
c       azero  = array of corresponding hard core diameters
c       insgf  = array of corresponding neutral species
c                  activity coefficient flags
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ipbtmx,npx2mx
c
      integer iz1,iz2,npx2
c
      real(8) alpha(ipbtmx,npx2mx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,iaz1,iaz2
c
c-----------------------------------------------------------------------
c
c     Assign standard alpha values for the npx2-th species pair,
c     according the charge combination.
c
      iaz1 = abs(iz1)
      iaz2 = abs(iz2)
c
      do i = 1,ipbtmx
        alpha(i,npx2) = 0.
      enddo
c
      if (iaz1.gt.0 .or. iaz2.gt.0) then
        if (iaz1.eq.1 .or. iaz2.eq.1) then
          alpha(1,npx2) = 2.
        else
          alpha(1,npx2) = 1.4
          alpha(2,npx2) = 12.0
        endif
      endif
c
      end
