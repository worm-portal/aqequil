      subroutine chksir(delxi,dlxmin,drir0,dxval0,eps100,iodb,
     $ nodbmx,nord,noutpt,nrd1mx,nttyo,rirec0,rirecp,xi0,xi1,xval0)
c
c     This subroutine checks the sign of the inverse rate. It finds the
c     point of reaction progress at which the inverse rate becomes
c     zero. Technically, this should never happen. The inverse rate
c     should start with a positive value and increase toward infinity.
c     The inverse rate is tracked using finite differences.
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
      integer nodbmx,nrd1mx
c
      integer noutpt,nttyo
c
      integer iodb(nodbmx)
c
      integer nord
c
      real*8 drir0(nrd1mx),dxval0(nrd1mx)
c
      real*8 delxi,dlxmin,eps100,rirec0,rirecp,xi0,xi1,xval0
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer icount,ier,ilsign,n
c
      character*48 usearch
      character*24 unam24
c
      real*8 dxp,tolsx,xtargv
c
      real*8 fctrl
c
c-----------------------------------------------------------------------
c
      data usearch /'the inverse rate becomes zero                   '/
      data unam24 /'                        '/
c
c-----------------------------------------------------------------------
c
      if (nord .le. 0) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Note that the search target is not where the inverse rate is zero,
c     but eps100. The interval for convergence is (0,2*eps100).
c
      ilsign = 1
      xtargv = eps100
      tolsx = eps100
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      icount = -1
  100 xi1 = xi0 + delxi
      icount = icount + 1
      if (icount .gt. 10) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
          write (noutpt,1000)
 1000     format(/3x,'A scan for the inverse rate becoming zero',
     $    ' has failed.',/3x,'Dropping to the minimum step size.')
        endif
        delxi = dlxmin
        go to 999
      endif
c
c     Estimate the inverse rate from a Taylor's series expansion.
c
      rirecp = rirec0
      dxp = 1.
      do n = 1,nord
        dxp = dxp*delxi
        rirecp = rirecp + ( drir0(n)/fctrl(n) )*dxp
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        if (rirecp .ge. 0.) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (delxi .gt. dlxmin) then
        if (rirec0 .le. (2.*eps100)) then
c
c         The inverse rate is too small at the base point.
c
          delxi = dlxmin
          if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1010) xi0,delxi
 1010       format(/' --- The inverse rate is already very close to',
     $      /7x,'zero at Xi= ',1pe11.4,'. Have set delxi equal to the',
     $      /7x,'minimum value of ',e11.4,'.')
            go to 100
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (delxi .gt. dlxmin) then
c
        xval0 = rirec0
        do n = 1,nord
          dxval0(n) = drir0(n)
        enddo
c
        call search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,
     $  nodbmx,nord,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,
     $  xtargv,xval0)
c
        if (ier .le. 0) go to 100
        if (ier .ge. 2) then
c
c         Note: if ier = 1, the returned "safe" value of delxi
c         is used.
c
          delxi = dlxmin
        endif
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(1).gt.0 .or. iodb(5).gt.0) then
        if (abs(rirecp - xtargv) .le. tolsx) then
          write (noutpt,1020) xi1,delxi
 1020     format(/3x,"Taylor's series predict that the inverse rate",
     $    /3x,'is zero at Xi= ',1pe11.4,', delxi= ',1pe11.4,'.')
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
