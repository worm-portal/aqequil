      subroutine fpbnpp(affp,affp0,aftarg,daffp0,delxi,dlxmin,dxval0,
     $ eps100,iodb,iopt,jpflag,nodbmx,noptmx,nord,nordmx,noutpt,npchk,
     $ npt,nptmax,nrd1mx,nttyo,tolaft,tolsat,uphase,xi0,xi1,xval0)
c
c     This subroutine finds the phase boundary at which a phase appears
c     in the equilibrium system (ES). Note that the phase affinities
c     are tracked using finite differences.
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
      integer nodbmx,nordmx,noptmx,nptmax,nrd1mx
c
      integer noutpt,nttyo
c
      integer iodb(nodbmx),iopt(noptmx),jpflag(nptmax),npchk(nptmax)
c
      integer nord,npt
c
      character*24 uphase(nptmax)
c
      real*8 affp(nptmax),affp0(nptmax),daffp0(nordmx,nptmax),
     $ dxval0(nrd1mx)
c
      real*8 aftarg,delxi,dlxmin,eps100,tolaft,tolsat,xi0,xi1,xval0
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      character*48 usearch
      character*24 unam24
c
      integer icount,ier,ilsign,j2,kpsst,n,np,npsst
c
      integer ilnobl
c
      real*8 afx,dxsv,tolsx,xtargv,xval
c
c-----------------------------------------------------------------------
c
      data usearch /'a phase just supersaturates                     '/
c
c-----------------------------------------------------------------------
c
      if (nord .le. 0) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The search target is not where the affinity is zero, but aftarg.
c     The interval for convergence is (tolsat,tolsst).  The target
c     (aftarg) is the midpoint of this interval.
c
      ilsign = -1
      xtargv = aftarg
      tolsx = tolaft
      dxsv = delxi
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      icount = -1
  100 xi1 = xi0 + delxi
      icount = icount + 1
      if (icount .gt. 50) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
          write (noutpt,1000)
 1000     format(/3x,'A scan for product phases just supersaturating',
     $    ' has failed.',/3x,'Dropping to the minimum step size.')
        endif
        delxi = dlxmin
        go to 999
      endif
c
c     Estimate the affinities of formation of the various phases from
c     Taylor's series expansions.
c
      call ataylr(delxi,daffp0,nord,nordmx,npt,nptmax,affp0,affp)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find any predicted supersaturations. Note that these count only
c     if they exceed the tolerance tolsat.
c
      xval = 0.
      npsst = 0
      kpsst = 0
      unam24 = ' '
      do np = 1,npt
        afx = affp(np)
        if (afx .ge. tolsat) then
          if (jpflag(np).le.0 .and. jpflag(np).ne.-1) then
            if (npchk(np) .ne. 1) then
              kpsst = kpsst + 1
              if (afx .gt. xval) then
                xval = affp(np)
                npsst = np
                unam24 = uphase(np)
              endif
            endif
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (kpsst .gt. 0) then
        if (delxi .gt. dlxmin) then
          if (abs(xval - xtargv) .gt. tolsx) then
c
            xval0 = affp0(npsst)
            do n = 1,nord
              dxval0(n) = daffp0(n,npsst)
            enddo
c
            call search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,
     $      nodbmx,nord,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,
     $      xtargv,xval0)
c
            if (ier .le. 0) go to 100
            if (ier .ge. 2) then
c
c             Note: if ier = 1, the returned "safe" value of delxi
c             is used.
c
              delxi = dlxmin
            endif
            go to 999
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(1).gt.0 .or. iodb(5).gt.0 .or. iopt(3).ge.1) then
        if (kpsst .gt. 0) then
          write (noutpt,1010) xi1,delxi
 1010     format(/3x,"Taylor's series predict the appearance of the",
     $    ' following',/3x,'phases in the ES at Xi= ',1pe11.4,
     $    ', delxi= ',e11.4,':',/)
c
          do np = 1,npt
            afx = affp(np)
            if (abs(afx - xtargv) .le. tolsx) then
              if (jpflag(np).le.0 .and. jpflag(np).ne.-1) then
                if (npchk(np) .ne. 1) then
                  j2 = ilnobl(uphase(np))
                  write (noutpt,1020) uphase(np)(1:j2)
 1020             format(5x,a)
                endif
              endif
            endif
          enddo
c
        endif
      endif
c
      if (iopt(3) .ge. 1) delxi = dxsv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
