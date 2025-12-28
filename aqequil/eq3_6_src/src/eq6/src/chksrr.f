      subroutine chksrr(delxi,dlxmin,drer0,dxval0,eps100,iodb,jreac,
     $ nodbmx,noutpt,nord,nrct,nrctmx,nrd1mx,nttyo,rrelr0,rrelrp,
     $ tolsrr,ureac,xi0,xi1,xval0)
c
c     This subroutine checks the signs of the relative rates. It finds
c     the point of reaction progress at which any relative rate of an
c     irreversible reaction becomes zero. The relative rates are
c     tracked using finite differences.
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
      integer nodbmx,nrctmx,nrd1mx
c
      integer noutpt,nttyo
c
      integer iodb(nodbmx),jreac(nrctmx)
c
      integer nord,nrct
c
      character*24 ureac(nrctmx)
c
      real*8 drer0(nrd1mx,nrctmx),dxval0(nrd1mx),rrelr0(nrctmx),
     $ rrelrp(nrctmx)
c
      real*8 delxi,dlxmin,eps100,tolsrr,xi0,xi1,xval0
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer icount,ier,ilsign,j2,krzero,n,nrc,nrzero
c
      integer ilnobl
c
      character*48 usearch
      character*24 unam24
c
      real*8 arrxp,atgsrr,dxp,rrx,rrx0,rrxp,tolsx,xtargv,xval
c
      real*8 fctrl
c
c-----------------------------------------------------------------------
c
      data usearch /'a relative rate changes sign                    '/
c
c-----------------------------------------------------------------------
c
      if (nord .le. 0) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The search target is not where a relative rate is zero. The
c     magnitude of the search target is atgsrr = 0.5*tolsrr. The search
c     tolerance has the same value. Thus, in the case of a relative rate
c     going from positive to negative values, the search target is
c     -atgsrr and the interval for convergence is (-tolsrr,0.). In the
c     case of a relative rate going from negative to positive values,
c     the search target is atgsrr, and the interval for convergence is
c     (0.,tolsrr).
c
      atgsrr = 0.5*tolsrr
      tolsx = atgsrr
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      icount = -1
  100 xi1 = xi0 + delxi
      icount = icount + 1
      if (icount .gt. 50) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
          write (noutpt,1000)
 1000     format(/3x,'A scan for reactant relative rates changing sign',
     $    ' has failed.',/3x,'Dropping to the minimum step size.')
        endif
        delxi = dlxmin
        go to 999
      endif
c
c     Estimate the relative rates from Taylor's series expansions.
c
      do nrc = 1,nrct
        rrx = 0.
        if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
          rrx = rrelr0(nrc)
          dxp = 1.
          do n = 1,nord
            dxp = dxp*delxi
            rrx = rrx + ( drer0(n,nrc)/fctrl(n) )*dxp
          enddo
        endif
        rrelrp(nrc) = rrx
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find any cross-overs. Note that there are two kinds, positive to
c     negative, and negative to positive.
c
      xval = 0.
      nrzero = 0
      krzero = 0
      unam24 = ' '
      do nrc = 1,nrct
        if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
          rrx0 = rrelr0(nrc)
          rrxp = rrelrp(nrc)
          if (rrx0 .gt. 0.) then
            if (rrxp .lt. -tolsrr) then
c
c             A case has been found in which the cross over tolerance
c             is exceeded for the positive to negative case.
c
              krzero = krzero + 1
              arrxp = -rrxp
              if (arrxp .gt. xval) then
                xval = arrxp
                nrzero = nrc
                unam24 = ureac(nrc)
                ilsign = 1
                xtargv = -atgsrr
              endif
            endif
          elseif (rrx0 .lt. 0.) then
            if (rrxp .gt. tolsrr) then
c
c             A case has been found in which the cross over tolerance
c             is exceeded for the negative to positive case.
c
              krzero = krzero + 1
              arrxp = rrxp
              if (arrxp .gt. xval) then
                xval = arrxp
                nrzero = nrc
                unam24 = ureac(nrc)
                ilsign = -1
                xtargv = atgsrr
              endif
            endif
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (krzero .gt. 0) then
        if (delxi .gt. dlxmin) then
          if (abs(xval - xtargv) .gt. tolsx) then
c
            xval0 = rrelr0(nrzero)
            do n = 1,nord
              dxval0(n) = drer0(n,nrzero)
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
      if (iodb(1).gt.0 .or. iodb(5).gt.0) then
        if (krzero .gt. 0) then
          if (abs(xval - xtargv) .le. tolsx) then
            j2 = ilnobl(unam24)
            write (noutpt,1010) unam24(1:j2),xi1,delxi
 1010       format(/" --- Taylor's series predict a zero relative",
     $      /7x,'rate for ',a,' at Xi= ',1pe11.4,',',
     $      /7x,'delxi= ',1pe11.4,' ---')
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
