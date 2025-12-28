      subroutine ldlxrd(delxi,dlxmin,dzvc0,fdzv0,iodb,ipndx1,kdim,
     $ km1,kmax,kxt,loph,nodbmx,nord,noutpt,nptmax,nrd1mx,nttyo,
     $ uphase,zklogu,zvec0)
c
c     This subroutine limits delxi when a phase is rapidly disappearing
c     from the equilibrium system. The particular mechanism here is
c     only designed to slow things down enough so that other mechanisms
c     can efficiently locate the phase disappearance boundary.
c     This mechanism is not employed to increase information density
c     around the boundary.
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
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer kmax,nodbmx,nptmax,nrd1mx
c
      integer noutpt,nttyo
c
      integer iodb(nodbmx),ipndx1(kmax)
c
      integer kdim,km1,kxt,nord
c
      character*24 uphase(nptmax)
c
      real*8 dzvc0(nrd1mx,kmax),fdzv0(nrd1mx,kmax),loph(nptmax),
     $ zvec0(kmax)
c
      real*8 delxi,dlxmin,zklogu
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j,j2,kcol,np
c
      integer ilnobl
c
      character*24 unamph
c
      real*8 dlxic,dxsv
c
c-----------------------------------------------------------------------
c
      dxsv = delxi
      unamph = 'Error                   '
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do kcol = km1,kxt
        if (fdzv0(1,kcol) .ge. 0.) go to 100
        do j = 1,nord
          if (dzvc0(j,kcol) .ge .0.) go to 100
        enddo
        np = ipndx1(kcol)
        if (loph(np) .le. zklogu) go to 100
c
c       Here dlxic is the value of delxi at which approximately
c       ninety per cent of the existing mass is destroyed.
c
        dlxic = -0.9*zvec0(kcol)/dzvc0(1,kcol)
        if (dlxic .lt. delxi) then
          delxi = dlxic
          delxi = max(delxi,dlxmin)
          unamph = uphase(np)
        endif
  100   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(5) .gt. 0) then
        if (delxi .lt. dxsv) then
          j2 = ilnobl(unamph)
          write (noutpt,1000) dxsv,delxi,unamph(1:j2)
          write (nttyo,1000) dxsv,delxi,unamph(1:j2)
 1000     format(/' * Note - (EQ6/ldlxrd) The step size has been cut',
     $    ' from ',1pe12.5,/7x,'to ',e12.5,' to limit the dissolution',
     $    ' of ',a,'.')
        endif
      endif
c
      end
