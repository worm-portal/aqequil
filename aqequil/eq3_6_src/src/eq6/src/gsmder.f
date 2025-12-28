      subroutine gsmder(delxia,drer0,drer0s,drir0,drir0s,jreac,nord,
     $ nrct,nrctmx,nrd1mx)
c
c     This subroutine computes smoothed or averaged derivatives of the
c     elements of the r vector. The actual derivatives (of various
c     orders) are averaged over the interval (-delxia,+delxia). The
c     base point (point 0) is at the center of this interval.
c
c     Subroutine gsmdez.f performs the same function for elements of
c     the z vector.
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
      integer nrctmx,nrd1mx
c
      integer jreac(nrctmx)
c
      integer nord,nrct
c
      real(8) drer0(nrd1mx,nrctmx),drer0s(nrd1mx,nrctmx),drir0(nrd1mx),
     $ drir0s(nrd1mx)
c
      real(8) delxia
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer k,n,nmax,nrc
c
      real(8) drx,dx
c
c-----------------------------------------------------------------------
c
c     Zero the averaged derivative arrays.
c
      do n = 1,nrd1mx
        drir0s(n) = 0.
      enddo
c
      nmax = nrd1mx*nrctmx
      call initaz(drer0s,nmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The averaged derivatives are calculated over the interval
c     (-delxia,+delxia), about the base point (point 0).
c
      dx = 2.*delxia
c
c     Calculate average derivatives for the inverse rate.
c
      do n = 1,nord
        drx = drir0(nord)
        do k = nord - 1,n,-1
          drx = drir0(k) + (drx*dx/(nord - n + 1))
        enddo
        drir0s(n) = drx
      enddo
c
c     Calculate average derivatives for relative rates.
c
      do nrc = 1,nrct
        if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1)  then
          do n = 1,nord
            drx = drer0(nord,nrc)
            do k = nord - 1,n,-1
              drx = drer0(k,nrc) + (drx*dx/(nord - n + 1))
            enddo
            drer0s(n,nrc) = drx
          enddo
        endif
      enddo
c
  999 continue
      end
