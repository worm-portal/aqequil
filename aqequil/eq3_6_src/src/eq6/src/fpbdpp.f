      subroutine fpbdpp(delxi,demop0,dlxmin,dxval0,emop,emop0,
     $ eps100,iemop,iodb,iopt,nodbmx,noptmx,nord,nordmx,noutpt,
     $ npet,nrd1mx,npetmx,nptmax,nttyo,uphase,xi0,xi1,xval0)
c
c     This subroutine finds the boundary at which a phase or phases
c     disappear from the equilibrium system (ES). The numbers of moles
c     of phases in the ES are tracked using finite differences.
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
      integer nodbmx,noptmx,nordmx,npetmx,nptmax,nrd1mx
c
      integer noutpt,nttyo
c
      integer iemop(npetmx),iodb(nodbmx),iopt(noptmx)
c
      integer nord,npet
c
      character*24 uphase(nptmax)
c
      real*8 demop0(nordmx,npetmx),dxval0(nrd1mx),emop(npetmx),
     $ emop0(npetmx)
c
      real*8 delxi,dlxmin,eps100,xi0,xi1,xval0
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer icount,ier,ilsign,j2,kpdis,n,np,npe,npedis
c
      integer ilnobl
c
      character*48 usearch
      character*24 unam24
c
      real*8 dxsv,mxx,mxx0,tolsx,xtargv,xval
c
c-----------------------------------------------------------------------
c
      data usearch /'a product phase vanishes                        '/
c
c-----------------------------------------------------------------------
c
      if (nord .le. 0) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The search target is not where the number of moles is
c     zero, but -eps100. The interval for convergence is (-2*eps100,0.).
c
      ilsign = 1
      xtargv = -eps100
      tolsx = eps100
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
 1000     format(/3x,'A scan for product phases vanishing',
     $    ' has failed.',/3x,'Dropping to the minimum step size.')
        endif
        delxi = dlxmin
        go to 999
      endif
c
c     Make a Taylor's series expansion of the number of moles of the
c     phases in the ES.
c
      call ptaylr(delxi,demop0,emop0,emop,nord,nordmx,npet,npetmx)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find any disappearing phases. Note that these count only if the
c     predicted number of moles is less than or equal to zero.
c
      xval = 1.0
      npedis = 0
      kpdis = 0
      unam24 = ' '
      do npe = 1,npet
        mxx = emop(npe)
        if (mxx .le. 0) then
          mxx0 = emop0(npe)
          if (mxx.lt.mxx0 .and. mxx0.gt.eps100) then
            kpdis = kpdis + 1
            if (mxx .lt. xval) then
              xval = mxx
              npedis = npe
              np = iemop(npe)
              unam24 = uphase(np)
            endif
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (kpdis .gt. 0) then
        if (delxi .gt. dlxmin) then
          if (abs(xval - xtargv) .gt. tolsx) then
c
            xval0 = emop0(npedis)
            do n = 1,nord
              dxval0(n) = demop0(n,npedis)
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
        if (kpdis .gt. 0) then
          write (noutpt,1010) xi1,delxi
 1010     format(/3x,"Taylor's series predict disappearance of the",
     $    ' following',/3x,'ES phases at Xi= ',1pe11.4,', delxi= ',
     $    e11.4,':',/)
c
          do npe = 1,npet
            mxx = emop(npe)
            if (abs(mxx - xtargv) .le. tolsx) then
              mxx0 = emop0(npe)
              if (mxx.lt.mxx0 .and. mxx0.gt.eps100) then
                np = iemop(npe)
                j2 = ilnobl(uphase(np))
                write (noutpt,1020) uphase(np)(1:j2)
 1020           format(5x,a)
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
