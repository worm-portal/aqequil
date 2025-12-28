      subroutine chksar(afrc0,afrcp,dafrc0,delxi,dlxmin,dxval0,
     $ eps100,iodb,jreac,nodbmx,noutpt,nord,nordmx,nrct,nrctmx,
     $ nrd1mx,nrk,nttyo,tolsar,ureac,xi0,xi1,xval0)
c
c     This subroutine checks the signs of the affinities of the
c     reactants. It finds the point of reaction progress at which the
c     affinity of any irreversible reaction changes sign. The reactant
c     affinities are tracked using finite differences.
c
c     This mechanism is used for two purposes. One is to insure that
c     a rate law for say the forward direction is not extrapolated
c     to the region of unfavorable affinity when a distinct rate
c     law is specified for the opposite direction. The second purpose
c     of this mechanism is to find the point at which a formerly
c     exhausted reactant must be reactivated so that it can be
c     precipitated according to a specified rate law (not a partial
c     equilibrium condition).
c
c     Note the following:
c
c       jreac = reactant status flag:
c
c                  0 = set to react
c                 -1 = saturated, but the remaining reactant mass
c                        continues to react irreversibly
c                  1 = exhausted
c                  2 = saturated; the status of any remaining reactant
c                        mass is changed to that of a product phase
c
c     When jreac = -1 or 2, the affinity is zero by definition. In this
c     case, finite-difference expressions of the affinity should be
c     avoided because the affinity values used to build them may
c     reflect somewhat random non-zero values consistent with the
c     applied convergence tolerances. The corresponding finite
c     differences may therefore behave wildly.
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
      integer nodbmx,nordmx,nrctmx,nrd1mx
c
      integer noutpt,nttyo
c
      integer iodb(nodbmx),jreac(nrctmx),nrk(2,nrctmx)
c
      integer nord,nrct
c
      character*24 ureac(nrctmx)
c
      real*8 afrc0(nrctmx),afrcp(nrctmx),dafrc0(nordmx,nrctmx),
     $ dxval0(nrd1mx)
c
      real*8 delxi,dlxmin,eps100,tolsar,xi0,xi1,xval0
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
      real*8 aaxxp,atgsar,axx,axx0,axxp,dxp,tolsx,xtargv,xval
c
      real*8 fctrl
c
c-----------------------------------------------------------------------
c
      data usearch /'a reactant affinity changes sign                '/
c
c-----------------------------------------------------------------------
c
      if (nord .le. 0) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The search target is not where a reactant affinity is zero. The
c     magnitude of the search target is atgsar = 0.5*tolsar. The search
c     tolerance has the same value. Thus, in the case of a reactant
c     affinity going from positive to negative values, the search target
c     is -atgsar and the interval for convergence is (-tolsar,0.). In
c     the case of a reactant affinity going from negative to positive
c     values, the search target is atgsar, and the interval for
c     convergence is (0.,tolsar).
c
      atgsar = 0.5*tolsar
      tolsx = atgsar
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      icount = -1
  100 xi1 = xi0 + delxi
      icount = icount + 1
      if (icount .gt. 50) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
          write (noutpt,1000)
 1000     format(/3x,'A scan for reactant affinities changing sign',
     $    ' has failed.',/3x,'Dropping to the minimum step size.')
        endif
        delxi = dlxmin
        go to 999
      endif
c
c     Estimate the reactant affinities from Taylor's series expansions.
c     Avoid using Taylor's series for reactants whose affinities must
c     be zero.
c
      do nrc = 1,nrct
        afrcp(nrc) = 0.
        if (jreac(nrc).eq.0 .or. jreac(nrc).eq. 1) then
          axx = afrc0(nrc)
          dxp = 1.
          do n = 1,nord
            dxp = dxp*delxi
            axx = axx + ( dafrc0(n,nrc)/fctrl(n) )*dxp
          enddo
          afrcp(nrc) = axx
        endif
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
        if (jreac(nrc).eq.0 .or. jreac(nrc).eq. 1) then
          axx0 = afrc0(nrc)
          axxp = afrcp(nrc)
          if (axx0 .gt. 0.) then
            if (axxp .lt. -tolsar) then
c
c             A case has been found in which the cross over tolerance
c             is exceeded for the positive to negative case.
c
              if (nrk(2,nrc) .ne. 0) then
c
c               The case is not one in which precipitation is to be
c               governed by instantaneous equilibrium. In that case,
c               let other controls apply.
c
                krzero = krzero + 1
                aaxxp = -axxp
                if (aaxxp .gt. xval) then
                  xval = aaxxp
                  nrzero = nrc
                  unam24 = ureac(nrc)
                  ilsign = 1
                  xtargv = -atgsar
                endif
              endif
            endif
          elseif (axx0 .lt. 0.) then
            if (axxp .gt. tolsar) then
c
c             A case has been found in which the cross over tolerance
c             is exceeded for the negative to positive case.
c
              if (jreac(nrc) .ne. 1) then
c
c               The case is not one in which any new forward reaction
c               is impossible because the reactant is exhausted.
c               That case may be ignored.
c
                krzero = krzero + 1
                aaxxp = axxp
                if (aaxxp .gt. xval) then
                  xval = aaxxp
                  nrzero = nrc
                  unam24 = ureac(nrc)
                  ilsign = -1
                  xtargv = atgsar
                endif
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
            xval0 = afrc0(nrzero)
            do n = 1,nord
              dxval0(n) = dafrc0(n,nrzero)
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
c
c       Print any occurrences of reactant affinity sign changes.
c
        if (krzero .gt. 0) then
          do nrc = 1,nrct
            if (jreac(nrc).eq.0 .or. jreac(nrc).eq. 1) then
              axx0 = afrc0(nrc)
              axxp = afrcp(nrc)
              if (abs(axxp) .le. tolsx) then
                if (axx0.gt.tolsx .and. dafrc0(1,nrc).lt.0.) then
                  if (nrk(2,nrc) .ne. 0) then
                    j2 = ilnobl(ureac(nrc))
                    write (noutpt,1010) ureac(nrc)(1:j2),xi1,delxi
 1010               format(/3x,"Taylor's series predict that the",
     $              ' affinity for reactant',/3x,a,' has decreased to',
     $              ' nearly zero at Xi= ',1pe12.5,',',/3x,'delxi= ',
     $              1pe12.5,'.')
                  endif
                elseif (axx0.lt.tolsx .and. dafrc0(1,nrc).gt.0.) then
                  j2 = ilnobl(ureac(nrc))
                  write (noutpt,1020) ureac(nrc)(1:j2),xi1,delxi
 1020             format(/3x,"Taylor's series predict that the",
     $            ' affinity for reactant',/3x,a,' has increased to',
     $            ' nearly zero at Xi= ',1pe12.5,',',/3x,'delxi= ',
     $            1pe12.5,'.')
                endif
              endif
            endif
          enddo
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
