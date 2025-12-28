      subroutine search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,
     $ nodbmx,nord,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,
     $ xtargv,xval0)
c
c     This subroutine finds the value of delxi at which occurs an event
c     whose type is described by the string in the usearch variable.
c
c     This subroutine is called by:
c
c       EQ6/chksar.f
c       EQ6/chksir.f
c       EQ6/chksrr.f
c       EQ6/chktmx.f
c       EQ6/chktpl.f
c       EQ6/chktpr.f
c       EQ6/ckawmn.f
c       EQ6/ckawmx.f
c       EQ6/ckehmn.f
c       EQ6/ckehmx.f
c       EQ6/cko2mn.f
c       EQ6/cko2mx.f
c       EQ6/ckphmn.f
c       EQ6/ckphmx.f
c       EQ6/fpbdpp.f
c       EQ6/fpbflo.f
c       EQ6/fpbnpp.f
c       EQ6/fpexrc.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       unam24 = name of entity, one of whose properties is involved in
c                  the search
c       ilsign = expected sign of the residual function at the
c                  left boundary (delxi = 0)
c       tolsx  = the convergence tolerance
c
c     Principal output:
c
c       ier     = error flag
c                   = 1  the search did not converge
c                   = 2  event being searched for appears to have
c                          occurred prior to the interval being
c                          searched
c                   = 3  event being searched for doesn't appear to
c                          occur in the interval being searched
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
      integer ier,ilsign,nord
c
      character*48 usearch
      character*24 unam24
c
      real*8 dxval0(nrd1mx)
c
      real*8 delxi,dlxmin,eps100,tolsx,xtargv,xval0
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer iter,j2,j3
c
      integer ilnobl
c
      real*8 x(2),y(2)
c
      real*8 adx,ares,ares0,resfnc,resx,slope,xleft,xnew,xright
c
c-----------------------------------------------------------------------
c
      ier = 0
c
c     The following local variables are signficant:
c
c       iter   = iteration counter
c       resfnc = convergence function
c
      iter = 0
      resfnc = 0.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      j2 = ilnobl(usearch)
      j3 = ilnobl(unam24)
c
      if (iodb(7) .ge. 1) then
        write (noutpt,1000) usearch(1:j2)
 1000   format(/3x,'--- Search for where ',a,' ---',/)
        if (j3 .gt. 0) write (noutpt,1010) unam24(1:j3)
 1010   format(9x,'The entity involved is ',a,'.',/)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Save the entering value of delxi in xright.
c
c       xleft  = the left boundary of the interval being examined
c       xright = the right boundary
c
      xleft = 0.
      xright = delxi
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check the left boundary point.
c
c       The current left point is ( x(1),y(1) ).
c       The current right point is ( x(2),y(2) ).
c
      x(1) = xleft
      delxi = xleft
c
      call sfncge(delxi,xval0,xtargv,dxval0,nord,nrd1mx,resx)
c
      y(1) = resx
      ares0 = abs(resx)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(7) .ge. 1) write (noutpt,1020) iter,x(1),resx,resfnc
 1020 format(3x,'iter= ',i2,' delxi= ',1pe11.4,' residual= ',e11.4,
     $ ' resfnc=',e11.4)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for a direct hit on the left boundary.
c
      if (ares0 .le. 0.) then
        delxi = dlxmin
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if ( (resx.lt.0 .and. ilsign.gt.0) .or.
     $  (resx.gt.0 .and. ilsign.lt.0) ) then
        write (noutpt,1030) usearch(1:j2)
        write (nttyo,1030) usearch(1:j2)
 1030   format(/' * Note - (EQ6/search) A search for where ',a,
     $  /7x,' indicates that the event has been stepped over.')
        if (j3 .gt. 0) then
          write (noutpt,1040) unam24(1:j3)
          write (nttyo,1040) unam24(1:j3)
 1040     format(7x,'The entity involved is ',a,'.')
        endif
        ier = 2
        delxi = dlxmin
        write (noutpt,1050) delxi
        write (nttyo,1050) delxi
 1050   format(7x,'The step size will be set to ',1pe11.4,'.')
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check the right boundary point.
c
  100 x(2) = xright
      delxi = xright
c
      call sfncge(delxi,xval0,xtargv,dxval0,nord,nrd1mx,resx)
c
      y(2) = resx
      ares = abs(resx)
c
      if (iodb(7) .ge. 1) write (noutpt,1020) iter,x(2),resx,resfnc
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for a direct hit on the right boundary.
c
      if (ares .le. 0.) then
        delxi = xright
        go to 999
      endif
c
      ares = min(ares0,ares)
c
      if ( (resx.lt.0 .and. ilsign.lt.0) .or.
     $  (resx.gt.0 .and. ilsign.gt.0) ) then
        write (noutpt,1120) usearch(1:j2)
        write (nttyo,1120) usearch(1:j2)
 1120   format(/' * Note - (EQ6/search) A search for where ',a,
     $  /7x,'indicates that the event being searched for does',
     $  /7x,'not take place in the interval being examined.')
        if (j3 .gt. 0) then
          write (noutpt,1040) unam24(1:j3)
          write (nttyo,1040) unam24(1:j3)
        endif
        write (noutpt,1130)
        write (nttyo,1130)
 1130   format(7x,'The step size will not be decreased.')
        ier = 3
        delxi = xright
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find a new point. Each new point will replace either the
c     left point or the right point, depending on the value of
c     the residual. At least one iteration is required. The
c     secant method is the principal algorithm for generating a
c     new point. However, if convergence is slow, interval halving
c     may occasionally be used instead.
c
  110 iter = iter + 1
      ares0 = ares
c
      if (iodb(7) .ge. 2) write (noutpt,1140) x(1),x(2)
 1140 format(7x,'Left point= ',1pe11.4,', Right point= ',e11.4)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iter.le.4 .or. resfnc.ge.0.5 .or. mod(iter,4).ne.0) then
c
c       Secant method.
c
        if (iodb(7) .ge. 1) write (noutpt,1150)
 1150   format(3x,'Secant method')
        slope =  (y(2) - y(1) )/( x(2) - x(1) )
        xnew = x(1) - (y(1)/slope)
      else
c
c       Interval halving method.
c
        if (iodb(7) .ge. 1) write (noutpt,1160)
 1160   format(3x,'Interval halving')
        xnew = 0.5*( x(1) + x(2) )
      endif
      delxi = xnew
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      call sfncge(delxi,xval0,xtargv,dxval0,nord,nrd1mx,resx)
c
      ares = abs(resx)
      resfnc = (ares0 - ares)/ares0
      if (iodb(7) .ge. 1) write (noutpt,1020) iter,xnew,resx,resfnc
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Note there are multiple convergence tests.
c
      if (ares .le. tolsx) go to 999
c
      if (xnew .gt. 0.) then
c
c       Check for underflow of the most recent correction.
c
        adx = abs((xnew - x(1))/xnew)
        if (adx .le. eps100) go to 999
      endif
c
      if (abs(xnew - x(1)) .le. dlxmin) then
c
c       Check for correction within the limit of the minimum
c       step size.
c
        delxi = max(delxi,dlxmin)
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iter .ge. 50) then
        delxi = max(x(1),dlxmin)
        write (noutpt,1170) usearch(1:j2)
        write (nttyo,1170) usearch(1:j2)
 1170   format(/' * Note - (EQ6/search) A search for where ',a,
     $  /7x,'has failed to converge within the maximum number of',
     $  ' iterations.')
        if (j3 .gt. 0) then
          write (noutpt,1040) unam24(1:j3)
          write (nttyo,1040) unam24(1:j3)
        endif
        write (noutpt,1180) delxi
        write (nttyo,1180) delxi
 1180   format(7x,'The step size will be set to a "safe" value of ',
     $  1pe11.4)
        ier = 1
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if ( (resx.lt.0 .and. ilsign.lt.0 ) .or.
     $  (resx.gt.0 .and. ilsign.gt.0) ) then
        x(1) = xnew
        y(1) = resx
      else
        x(2) = xnew
        y(2) = resx
      endif
c
      go to 110
c
  999 continue
      end
