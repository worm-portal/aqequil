      subroutine gsmdez(delxia,dzvc0,dzvc0s,kdim,kmax,nord,nrd1mx)
c
c     This subroutine computes smoothed or averaged derivatives of the
c     elements of the z vector. The actual derivatives (of various
c     orders) are averaged over the interval (-delxia,+delxia). The
c     base point (point 0) is at the center of this interval.
c
c     Subroutine gsmder.f performs the same function for elements of
c     the r vector.
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
      integer kmax,nrd1mx
c
      integer nord,kdim
c
      real(8) dzvc0(nrd1mx,kmax),dzvc0s(nrd1mx,kmax)
c
      real(8) delxia
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer k,kcol,n,nmax
c
      real(8) dx,dzx
c
c-----------------------------------------------------------------------
c
c     Zero the averaged derivative arrays.
c
      nmax = nrd1mx*kmax
      call initaz(dzvc0s,nmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The averaged derivatives are calculated over the interval
c     (-delxia,+delxia), about the base point (point 0).
c
      dx = 2.*delxia
c
c     Calculate average derivatives for elements of the z vector.
c
      do kcol = 1,kdim
        do n = 1,nord
          dzx = dzvc0(nord,kcol)
          do k = nord - 1,n,-1
            dzx = dzvc0(k,kcol) + (dzx*dx/(nord - n + 1))
          enddo
          dzvc0s(n,kcol) = dzx
        enddo
      enddo
c
  999 continue
      end
