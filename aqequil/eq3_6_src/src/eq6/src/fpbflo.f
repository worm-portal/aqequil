      subroutine fpbflo(al10,delxi,demop0,dlxmin,dxval0,d1emp1,
     $ d2emp1,emop,emop0,eps100,fdpe0,iemop,ier,iodb,nodbmx,nord,
     $ nordmx,noutpt,npet,npetmx,nptmax,nrd1mx,nttyo,qdump,toldl,
     $ uaqsln,ufixf,uphase,xim1,xi0,xi1,xval0,zklogu)
c
c     This subroutine limits delxi by approximate position of
c     significant maxima in the masses of non-aqeuous species that are
c     in partial equilibrium with the aqueous solution. The numbers of
c     moles of such phases are tracked using finite differences.
c
c     This subroutine is called by:
c
c       EQ6/path.f.
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
      integer nodbmx,nordmx,npetmx,nptmax,nrd1mx
c
      integer noutpt,nttyo
c
      integer iemop(npetmx),iodb(nodbmx)
c
      integer ier,nord,npet
c
      logical qdump
c
      character*24 uphase(nptmax)
      character*24 uaqsln
      character*8 ufixf
c
      real*8 demop0(nordmx,npetmx),dxval0(nrd1mx),d1emp1(npetmx),
     $ d2emp1(npetmx),emop(npetmx),emop0(npetmx),fdpe0(nordmx,npetmx)
c
      real*8 al10,delxi,dlxmin,eps100,toldl,xim1,xi0,xi1,xval0,zklogu
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer icount,ilsign,j2,n,np,npe,npej,npj
c
      integer ilnobl
c
      character*48 usearch
      character*24 unam24
c
      real*8 dfmrxx,dfmxx,dfmxxj,dpx,dpx0,mxx,mxxu,mxx0,tolsx,xtargv
c
      real*8 texp
c
c-----------------------------------------------------------------------
c
      data usearch /'a product phase maximizes                       '/
c
c-----------------------------------------------------------------------
c
      if (nord .le. 1) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The search target is not where the first derivative is zero, but
c     -eps100. The interval for convergence is (-2*eps100,0.).
c
      ilsign = 1
      xtargv = -eps100
      tolsx = eps100
      mxxu = texp(zklogu)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      icount = -1
  100 xi1 = xi0 + delxi
      icount = icount + 1
      if (icount .gt. 50) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
          write (noutpt,1000)
 1000     format(/3x,'A scan for product phase masses maximizing',
     $    ' has failed.',/3x,'Dropping to the minimum step size.')
        endif
        delxi = dlxmin
        go to 999
      endif
c
c     Make a Taylor's series expansion of the first derivatives of the
c     numbers of moles of the phases in the ES.
c
      call d1ptay(delxi,demop0,d1emp1,nord,nordmx,npet,npetmx)
c
c     Make a Taylor's series expansion of the numbers of moles of the
c     phases in the ES.
c
      call ptaylr(delxi,demop0,emop0,emop,nord,nordmx,npet,npetmx)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find any phases whose mole numbers are decreasing. Pick the
c     one whose mole number would be decreased the most at the
c     current step size. Some phase types including the aqueous
c     solution phase and any fictive fugacity-fixing phases are
c     exempt. Phases whose mole numbers are already sufficiently
c     small (.le. 10^zklogu at the base point) point are ignored.
c
      dfmxxj = 0.
      npej = 0
      npj = 0
      unam24 = ' '
      do npe = 1,npet
        np = iemop(npe)
        if (uphase(np)(1:24) .ne. uaqsln(1:24)) then
          if (uphase(np)(1:5) .ne. ufixf(1:5)) then
c
c           Not an excluded phase type.
c
            mxx0 = emop0(npe)
            if (mxx0 .gt. mxxu) then
c
c             The base point mole number is sufficiently large
c             to consider.
c
                mxx = emop(npe)
                dfmxx = mxx - mxx0
                if (dfmxx .lt. 0.) then
c
c                 A decrease in the mole number is predicted.
c
                  dfmrxx = dfmxx/mxx0
                  if (dfmrxx .lt. -0.0001) then
c
c                   The predicted decrease in mole number exceeds
c                   a small relative amount.
c
                    dpx0 = demop0(1,npe)
                    dpx = d1emp1(npe)
c
                    if (dpx0.ge.eps100 .and. dpx.le.-eps100) then
c
c                   The derivative of the mole number changes
c                   from positive at the based point to negative
c                   at the new point.
c
                    if (fdpe0(1,npe) .lt. -eps100) then
c
c                     The derivative is confirmed by the first
c                     finite difference.
c
                      if (dfmxx .lt. dfmxxj) then
c
c                       Choose the phase with the largest predicted
c                       decrease in mole number.
c
                        dfmxxj = dfmxx
                        npej = npe
                        npj = np
                      endif
                    endif
                  endif
                endif
              endif
            endif
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (npej .gt. 0) then
        qdump = .true.
        if (delxi .gt. dlxmin) then
          dpx = d1emp1(npej)
          if (abs(dpx + eps100) .gt. eps100) then
c
            unam24 = uphase(npj)
            xval0 = demop0(1,npej)
            do n = 1,nord - 1
              dxval0(n) = demop0(n + 1,npej)
            enddo
            if (nord .gt. 0) dxval0(nord) = 0.
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
      if (qdump) then
        if (iodb(1).gt.0 .or.iodb(5).gt.0) then
          write (noutpt,1010) xi1,delxi
 1010     format(/3x,"Taylor's series predict maxima for the number of",
     $    ' moles of',/3x,'the following phases at Xi= ',1pe11.4,
     $    ', delxi= ',e11.4,':',/)
c
          do npe = 1,npet
            mxx0 = emop0(npe)
            dpx = d1emp1(npe)
            np = iemop(npe)
c
            if (mxx0 .gt. mxxu) then
              if (abs(dpx + eps100) .gt. eps100) then
                if (uphase(np)(1:24) .ne. uaqsln(1:24)) then
                  if (uphase(np)(1:5) .ne. ufixf(1:5)) then
                    j2 = ilnobl(uphase(np))
                    write (noutpt,1020) uphase(np)(1:j2)
 1020               format(5x,a)
                  endif
                endif
              endif
            endif
          enddo
c
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
