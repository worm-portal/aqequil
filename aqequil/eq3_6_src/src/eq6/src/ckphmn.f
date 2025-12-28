      subroutine ckphmn(delxi,dlxmin,dph0,dxh0mx,dxval0,eps100,iodb,
     $ nodbmx,nord,noutpt,nrd1mx,nttyo,ph0,ph1,phmin,prcinf,qdump,
     $ tolxsu,xi0,xi1,xval0)
c
c     This subroutine checks to see that the requested minimum value of
c     pH is not exceeded.
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
      logical qdump
c
      real*8 dph0(nrd1mx),dxval0(nrd1mx)
c
      real*8 delxi,dlxmin,dxh0mx,eps100,prcinf,phmin,ph0,ph1,tolxsu,
     $ xi0,xi1,xval0,xtargv
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
      real*8 dxp,dxsv,tolsx
c
      real*8 fctrl
c
c-----------------------------------------------------------------------
c
      data usearch /'pH reaches the requested minimum value          '/
      data unam24  /'                       '/
c
c-----------------------------------------------------------------------
c
      dxsv = delxi
c
      dxh0mx = prcinf
      if (nord .le. 0) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The search target is the requested limit on the pH.
c
      ilsign = +1
      xtargv = phmin
      tolsx = tolxsu
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      icount = -1
  100 xi1 = xi0 + delxi
      icount = icount + 1
      if (icount .gt. 10) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
          write (noutpt,1000)
 1000     format(/3x,'A scan for where the pH reaches the requested',
     $    ' minimum value',/3x,'has failed. Dropping to the minimum',
     $    ' step size.')
        endif
        delxi = dlxmin
        go to 990
      endif
c
c     Estimate the pH from a Taylor's series expansion.
c
      ph1 = ph0
      dxp = 1.
      do n = 1,nord
        dxp = dxp*delxi
        ph1 = ph1 + ( dph0(n)/fctrl(n) )*dxp
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (delxi .gt. dlxmin) then
        if (ph0 .le. (1. + tolxsu)*phmin) then
c
c         The pH at the base point is too close to the requested
c         minimum value.
c
          delxi = dlxmin
          dxh0mx = delxi
          if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1010) ph0,phmin,xi0,delxi
 1010       format(/3x,'The base point pH of ',1pe11.4,' is already',
     $      ' very close',/3x,'to the requested value of ',
     $      e11.4,/3x,'at Xi= ',e11.4,'. Have set delxi',
     $      ' equal to the minimum',/7x,'value of ',e11.4,'.')
            go to 100
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (delxi .gt. dlxmin) then
        if ((ph1 - phmin) .lt. -tolsx) then
c
          xval0 = ph0
          do n = 1,nord
            dxval0(n) = dph0(n)
          enddo
c
          call search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,
     $    nodbmx,nord,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,
     $    xtargv,xval0)
c
          dxh0mx = delxi
          if (ier .le. 0) go to 100
          if (ier .ge. 2) then
c
c           Note: if ier = 1, the returned "safe" value of delxi
c           is used.
c
            delxi = dlxmin
          endif
          dxh0mx = delxi
          go to 990
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(1).gt.0 .or. iodb(5).gt.0) then
        if (abs(ph1 - phmin) .le. tolsx) then
          write (noutpt,1020) phmin,xi1,delxi
 1020     format(/3x,"Taylor's series predict that the pH is",
     $    ' at the requested',/3x,'minimum value of ',1pe11.4,
     $    ' at Xi= ',e11.4,',',/3x,'delxi= ',e11.4,'.')
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  990 if (dxsv .gt. delxi) qdump = .false.
c
  999 continue
      end
