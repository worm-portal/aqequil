      subroutine ckehmn(delxi,dlxmin,deh0,dxe0mx,dxval0,eps100,iodb,
     $ nodbmx,nord,noutpt,nrd1mx,nttyo,eh0,eh1,ehmin,prcinf,qdump,
     $ tolxsu,xi0,xi1,xval0)
c
c     This subroutine checks to see that the requested minimum value of
c     Eh is not exceeded.
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
      real*8 deh0(nrd1mx),dxval0(nrd1mx)
c
      real*8 delxi,dlxmin,dxe0mx,eps100,prcinf,ehmin,eh0,eh1,tolxsu,
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
      data usearch /'Eh reaches the requested minimum value          '/
      data unam24  /'                       '/
c
c-----------------------------------------------------------------------
c
      dxsv = delxi
c
      dxe0mx = prcinf
      if (nord .le. 0) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The search target is the requested limit on the Eh.
c
      ilsign = +1
      xtargv = ehmin
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
 1000     format(/3x,'A scan for where the Eh reaches the requested',
     $    ' minimum value',/3x,'has failed. Dropping to the minimum',
     $    ' step size.')
        endif
        delxi = dlxmin
        go to 990
      endif
c
c     Estimate the Eh from a Taylor's series expansion.
c
      eh1 = eh0
      dxp = 1.
      do n = 1,nord
        dxp = dxp*delxi
        eh1 = eh1 + ( deh0(n)/fctrl(n) )*dxp
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (delxi .gt. dlxmin) then
        if (eh0 .le. (1. + tolxsu)*ehmin) then
c
c         The Eh at the base point is too close to the requested
c         minimum value.
c
          delxi = dlxmin
          dxe0mx = delxi
          if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1010) eh0,ehmin,xi0,delxi
 1010       format(/3x,'The base point Eh of ',1pe11.4,' v is already',
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
        if ((eh1 - ehmin) .lt. -tolsx) then
c
          xval0 = eh0
          do n = 1,nord
            dxval0(n) = deh0(n)
          enddo
c
          call search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,
     $    nodbmx,nord,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,
     $    xtargv,xval0)
c
          dxe0mx = delxi
          if (ier .le. 0) go to 100
          if (ier .ge. 2) then
c
c           Note: if ier = 1, the returned "safe" value of delxi
c           is used.
c
            delxi = dlxmin
          endif
          dxe0mx = delxi
          go to 990
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(1).gt.0 .or. iodb(5).gt.0) then
        if (abs(eh1 - ehmin) .le. tolsx) then
          write (noutpt,1020) ehmin,xi1,delxi
 1020     format(/3x,"Taylor's series predict that the Eh is",
     $    ' at the requested',/3x,'minimum value of ',1pe11.4,
     $    ' v at Xi= ',e11.4,',',/3x,'delxi= ',e11.4,'.')
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  990 if (dxsv .gt. delxi) qdump = .false.
c
  999 continue
      end
