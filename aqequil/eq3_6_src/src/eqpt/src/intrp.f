      subroutine intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,
     $ narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $ xvec,yvec)
c
c     This subroutine fits interpolating polynomials to data (avgrid)
c     on a temperature grid (tempc). The grid is divided into ranges.
c     A separate polynomial is fitted to the data in each range.
c     Its coefficients are obtained in the cof array. The coefficients
c     for all ranges are returned in the apr array.
c
c     This subroutine is called by:
c
c       EQPT/pcraq.f
c       EQPT/pcrsg.f
c       EQPT/wrpar.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       avgrid = array containing the data on the temperature grid
c
c     Principal output:
c
c       apr    = array of polynomial coefficients (for all ranges)
c       avgrid = array containing the data on the temperature grid
c       narxmx = the maximum number of points or coefficients per
c                  temperature range
c       narxt  = the actual number of points or coefficients per
c                  temperature range
c       nptrmx = the maximum number of temperature ranges
c       nptrt  = the actual number of temperature ranges
c
c     Workspace:
c
c       aamatr = matrix used to calculate the polynomial coefficients
c       cof    = array of fitted polynomial coefficients (for a
c                  single temperature range)
c       xvec   = array of scaled temperatures corresponding to the
c                  data in the yvec array
c       yvec   = array of data to be fitted (for a single temperature
c                  range)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer narxmx,ntprmx
c
      integer noutpt,nttyo
c
      integer ntprt
c
      integer ipivot(narxmx),narxt(ntprmx)
c
      real*8 apr(narxmx,ntprmx),avgrid(narxmx,ntprmx)
      real*8 tempc(narxmx,ntprmx),tempcs(narxmx,ntprmx),tmpcmx(ntprmx)
      real*8 aamatr(narxmx,narxmx),gmmatr(narxmx,narxmx)
      real*8 cof(narxmx),xvec(narxmx),yvec(narxmx)
c
      real*8 eps100
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ier,n,nmax,npft,ntpr
c
      real*8 tvecmx
c
c-----------------------------------------------------------------------
c
c     Initialize apr to 0.
c
      nmax = narxmx*ntprmx
      call initaz(apr,nmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Loop on temperature ranges.
c
      do ntpr = 1,ntprt
c
c       Put all real values (9999999. indicates no data) in the yvec
c       array. Put the corresponding scaled temperature values in
c       the xvec array.
c
        npft = 0
c
        do n = 1,narxt(ntpr)
          if (avgrid(n,ntpr) .lt. 9999999.) then
            npft = npft + 1
            xvec(npft) = tempcs(n,ntpr)
            yvec(npft) = avgrid(n,ntpr)
          endif
        enddo
c
c       Here npft is the number of usable values.
c
        if (npft .le. 0) then
c
c         There are no usable values.
c
          apr(1,ntpr) = 9999999.
          go to 110
        endif
c
        if (npft .gt. 1) then
c
c         Check for constant yvec.
c
          do n = 2,npft
            if(abs(yvec(n) -yvec(1)) .gt. eps100) go to 100
          enddo
          npft = 1
  100     continue
        endif
c
        if (npft .le. 1) then
          apr(1,ntpr) = yvec(1)
          go to 110
        endif
c
c       Fit the polynomial.
c
        call polfit(aamatr,cof,gmmatr,ier,ipivot,narxmx,npft,
     $  noutpt,nttyo,xvec,yvec)
c
        if (ier .gt. 0) then
          write (noutpt,1010) ntpr,tempc(1,ntpr),tempc(narxt(ntpr),ntpr)
          write (nttyo,1010) ntpr,tempc(1,ntpr),tempc(narxt(ntpr),ntpr)
 1010     format(/' * Error - (EQPT/intrp) Could not compute the',
     $    /7x,'coefficients of an interpolating polynomial in',
     $    /7x,'temperature range ',i2,' (',f6.2,' - ',f6.2,' C.')
          stop
        endif
c
c       Rescale the coefficients.
c
        tvecmx = tmpcmx(ntpr)
c
c       Calling sequence substitutions:
c         tvecmx for avxmax
c         cof for avy
c         cof for avys
c         npft for nmax
c
        call rscaly(tvecmx,cof,cof,eps100,npft)
c
c       Store the fitted coefficients for this range in the apr array.
c
        do n = 1,npft
          apr(n,ntpr) = cof(n)
        enddo
c
  110   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
