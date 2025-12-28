      subroutine gakmat(akmat0,dxsm00,nord,nrd1mx)
c
c     This subroutine computes the akmat0 matrix, which allows
c     calculation of estimates of derivatives from corresponding finite
c     differences. Note: this subroutine presumes that the diagonal
c     elements have been previously set to factorials and that the
c     lower triangle elements have previously been set to zero. These
c     parts of the matrix are invariant.
c
c     This routine can also be used to calculate the akmat1 matrix from
c     the dlxsm1 array, provided that in the calling sequence akmat1
c     is substituted for akmat0 and dlxsm1 for dxsm00. The akmat0 matrix
c     and the dxsm00 array are used with finite differences employed
c     to generate predictor functions, while akmat1 and dlxsm1 are used
c     are used with finite differences employed to generate corrector
c     functions.
c
c     This subroutine is called by:
c
c       EQ6/path.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c     Principal output:
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
      integer nrd1mx
c
      integer nord
c
      real*8 akmat0(nrd1mx,nrd1mx),dxsm00(nrd1mx)
c
c----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,im1,ip1,j,jm1
c
c----------------------------------------------------------------------
c
      if (nord .ge. 2) then
c
        i = 1
        ip1 = i + 1
        do j = ip1,nord
          jm1 = j - 1
          akmat0(i,j) = akmat0(i,jm1)*dxsm00(jm1)
        enddo
c
        do i = 2,nord - 1
          ip1 = i + 1
          im1 = i - 1
          do j = ip1,nord
            jm1 = j - 1
            akmat0(i,j) = akmat0(i,jm1)*dxsm00(jm1) + akmat0(im1,jm1)*i
          enddo
        enddo
      endif
c
      end
