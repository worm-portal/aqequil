      subroutine ztaylr(delxi,dzvc0,kdim,kmax,km1,kxt,nord,nrd1mx,
     $ qztayl,zklogu,zvclg0,zvclg1,zvec0,zvec1)
c
c     This subroutine calculates new values for algebraic master
c     variables (the z vector elements) from the finite-difference-
c     based truncated Taylor's series. If qztayl = .true., change
c     limits are applied to the results. The purpose of these limits
c     is to assist the hybrid Newton-Raphson iteration either to
c     converge or to generate useful divergence diagnostics.
c
c     Compare with:
c
c       EQ6/ataylr.f
c       EQ6/ptaylr.f
c       EQ6/rtaylr.f
c
c      See also:
c
c       EQ6/d1ztay.f
c       EQ6/d2ztay.f
c
c     This subroutine is called by:
c
c       EQ6/eqshel.f
c       EQ6/ldlxrc.f
c       EQ6/path.f
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
c       qztayl = flag to apply change limits
c       zvclg0 = the log z vector at the base point
c       zvec0  = the z vector at the base point
c
c     Principal output:
c
c       zvec1  = the z vector at the new point
c       zvclg1 = the log z vector at the new point
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
      integer kdim,km1,kxt,nord
c
      logical qztayl
c
      real*8 dzvc0(nrd1mx,kmax),zvclg0(kmax),zvclg1(kmax),zvec0(kmax),
     $ zvec1(kmax)
c
      real*8 delxi,zklogu
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer kcol,n
c
      real*8 dxp,zx,zx0,zxl,zxu
c
      real*8 fctrl,tlg
c
c-----------------------------------------------------------------------
c
c
c     Compute the expansions from the Taylor's series.
c
      do kcol = 1,kdim
        zx = zvec0(kcol)
        dxp = 1.
        do n = 1,nord
          dxp = dxp*delxi
          zx = zx + ( dzvc0(n,kcol)/fctrl(n) )*dxp
        enddo
        zvec1(kcol) = zx
c
c       Compute the corresponding logarithmic variable. Provide
c       protection if zx is negative. If this is the case, the
c       logarithmic variable won't be used.
c
        if (zx .lt. 0) then
          zx = 0.
        endif
        zvclg1(kcol) = tlg(zx)
      enddo
c
      if (qztayl) then
c
c       Apply change limits.
c
        do kcol = 1,kdim
          zx0 = zvec0(kcol)
          zxl = 1.e-20*zx0
          zxu = 1.e+20*zx0
          zx = zvec1(kcol)
          zx = max(zx,zxl)
          zx = min(zx,zxu)
          if (zx .lt. 0.) then
            zx = 0.
          endif
          zvec1(kcol) = zx
          zvclg1(kcol) = tlg(zx)
        enddo
c
        do kcol = km1,kxt
          if (zvclg0(kcol) .lt. zklogu) then
            zvclg1(kcol) = zvclg0(kcol)
            zvec1(kcol) = zvec0(kcol)
          endif
        enddo
      endif
c
      end
