      subroutine fpexrc(delxi,dlxmin,drer0,dxval0,eps100,iodb,
     $ jreac,morr,morr0,nodbmx,nord,noutpt,nrct,nrctmx,nrd1mx,
     $ nttyo,qdump,rrelr0,ureac,xi0,xi1,xval0)
c
c     This subroutine finds the point of reaction progress at which a
c     reactant becomes exhausted.
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
      logical qdump
c
      character*24 ureac(nrctmx)
c
      real*8 drer0(nrd1mx,nrctmx),dxval0(nrd1mx),morr(nrctmx),
     $ morr0(nrctmx),rrelr0(nrctmx)
c
      real*8 delxi,dlxmin,dlxrct,dxp,eps100,mxx,xi0,xi1,xval0
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer icount,ier,ilsign,j2,kexrc,n,nexrc,nordp1,nrc
c
      integer ilnobl
c
      character*48 usearch
      character*24 unam24
c
      real*8 dxsv,tolsx,xtargv,xval
c
      real*8 fctrl
c
c-----------------------------------------------------------------------
c
      data usearch /'a reactant is exhausted                         '/
c
c-----------------------------------------------------------------------
c
      dxsv = delxi
c
c     The search target is not where the number of moles remaining is
c     zero, but -eps100. The interval for convergence is (-2*eps100,0).
c
      nordp1 = nord + 1
      ilsign = 1
      xtargv = -eps100
      tolsx = eps100
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      icount = -1
  100 xi1 = xi0 + delxi
      icount = icount + 1
      if (icount .gt. 50) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
          write (noutpt,1000)
 1000     format(/3x,'A scan for reactant masses becoming exhausted',
     $    ' has failed.',/3x,'Dropping to the minimum step size.')
        endif
        delxi = dlxmin
        go to 990
      endif
c
c     Estimate the number of moles remaining of each reactant from
c     Taylor's series expansions.
c
      do nrc = 1,nrct
        dlxrct = 0.
        if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
          dlxrct = rrelr0(nrc)*delxi
          dxp = delxi
          do n = 1,nord
            dxp = dxp*delxi
            dlxrct = dlxrct + ( drer0(n,nrc)/fctrl(n + 1) )*dxp
          enddo
        endif
c
cXX     The following assumes that each reactant has a reaction
cXX     coefficient of 1.
c
        morr(nrc) = morr0(nrc) - dlxrct
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find any cases of exhausted reactants.
c
      xval = 1.0
      nexrc = 0
      kexrc = 0
      unam24 = ' '
      do nrc = 1,nrct
        mxx = morr(nrc)
        if (mxx .le. 0.) then
          if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
            kexrc = kexrc + 1
            if (mxx .lt. xval) then
              xval = morr(nrc)
              nexrc = nrc
              unam24 = ureac(nrc)
            endif
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (kexrc .gt. 0) then
        if (delxi .gt. dlxmin) then
          if (abs(xval - xtargv) .gt. tolsx) then
c
            xval0 = morr0(nexrc)
c
cXX         The following assumes that each reactant has a reaction
cXX         coefficient of 1.
c
            dxval0(1) = -rrelr0(nexrc)
            do n = 2,nordp1
              dxval0(n) = -drer0(n - 1,nexrc)
            enddo
c
c           Calling sequence substitutions:
c             nordp1 for nord
c
            call search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,
     $      nodbmx,nordp1,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,
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
            go to 990
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(1).gt.0 .or. iodb(5).gt.0) then
        if (kexrc .gt. 0) then
          write (noutpt,1010) xi1,delxi
 1010     format(/3x,"Taylor's series predict exhaustion of the",
     $    ' following',/3x,'reactants at Xi= ',1pe11.4,
     $    ', delxi= ',e11.4,':',/)
c
          do nrc = 1,nrct
            mxx = morr(nrc)
            if (abs(mxx - xtargv) .le. tolsx) then
              if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
                j2 = ilnobl(ureac(nexrc))
                write (noutpt,1020) ureac(nexrc)(1:j2)
 1020           format(5x,a)
              endif
            endif
          enddo
c
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  990 if (dxsv .gt. delxi) qdump = .false.
c
  999 continue
      end
