      subroutine cko2mx(delxi,dlxmin,do20,dxo1mx,dxval0,eps100,iodb,
     $ nodbmx,nord,noutpt,nrd1mx,nttyo,fo2lg0,fo2lg1,o2max,prcinf,
     $ qdump,tolxsu,xi0,xi1,xval0)
c
c     This subroutine checks to see that the requested maximum value of
c     log fO2 is not exceeded.
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
      real*8 do20(nrd1mx),dxval0(nrd1mx)
c
      real*8 delxi,dlxmin,dxo1mx,eps100,prcinf,o2max,fo2lg0,fo2lg1,
     $ tolxsu,xi0,xi1,xval0,xtargv
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
      data usearch /'log fO2 reaches the requested maximum value     '/
      data unam24  /'                       '/
c
c-----------------------------------------------------------------------
c
      dxsv = delxi
c
      dxo1mx = prcinf
      if (nord .le. 0) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The search target is the requested limit on the log fO2.
c
      ilsign = -1
      xtargv = o2max
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
 1000     format(/3x,'A scan for where the log fO2 reaches the',
     $    ' requested maximum value',/3x,'has failed. Dropping to the',
     4    ' minimum step size.')
        endif
        delxi = dlxmin
        go to 990
      endif
c
c     Estimate the log fO2 from a Taylor's series expansion.
c
      fo2lg1 = fo2lg0
      dxp = 1.
      do n = 1,nord
        dxp = dxp*delxi
        fo2lg1 = fo2lg1 + ( do20(n)/fctrl(n) )*dxp
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (delxi .gt. dlxmin) then
        if (fo2lg0 .ge. (1. - tolxsu)*o2max) then
c
c         The log fO2 at the base point is too close to the requested
c         maximum value.
c
          delxi = dlxmin
          dxo1mx = delxi
          if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1010) fo2lg0,o2max,xi0,delxi
 1010       format(/3x,'The base point log fO2 of ',1pe11.4,' is',
     $      ' already very close',/3x,'to the requested value of ',
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
        if ((fo2lg1 - o2max) .gt. tolsx) then
c
          xval0 = fo2lg0
          do n = 1,nord
            dxval0(n) = do20(n)
          enddo
c
          call search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,
     $    nodbmx,nord,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,
     $    xtargv,xval0)
c
          dxo1mx = delxi
          if (ier .le. 0) go to 100
          if (ier .ge. 2) then
c
c           Note: if ier = 1, the returned "safe" value of delxi
c           is used.
c
            delxi = dlxmin
          endif
          dxo1mx = delxi
          go to 990
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(1).gt.0 .or. iodb(5).gt.0) then
        if (abs(fo2lg1 - o2max) .le. tolsx) then
          write (noutpt,1020) o2max,xi1,delxi
 1020     format(/3x,"Taylor's series predict that the log fO2 is",
     $    ' at the requested',/3x,'maximum value of ',1pe11.4,
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
