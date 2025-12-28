      subroutine chdxec(delxi,dlxmx0,dzvc0,iodb,kdim,kmax,km1,kxt,
     $ nodbmx,nord,noutpt,nrd1mx,qmin,qscon,scale,scalim,scnstd,
     $ scnsti,uzvec1,zvec0)
c
c     This subroutine chooses the step size (delxi) for economy and
c     super economy modes.
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
      integer kmax,nodbmx,nrd1mx
c
      integer noutpt
c
      integer iodb(nodbmx)
c
      integer kdim,km1,kxt,nord
c
      logical qmin,qscon
c
      character*48 uzvec1(kmax)
c
      real*8 dzvc0(nrd1mx,kmax),zvec0(kmax)
c
      real*8 delxi,dlxmx0,scale,scalim,scnstd,scnsti
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,j3,kcol,kcolcn
c
      integer ilnobl
c
      character*24 ubas,uzt,uzn
c
      real*8 aterm,d1,d2,d1sq,mxx0,rd1,rd2,ri1,ri2,rx,scfac,scnst,sx,
     $ termd,termi,thetad,thetai
c
c-----------------------------------------------------------------------
c
      data ubas /'basis variable          '/
c
c-----------------------------------------------------------------------
c
      if (qscon) then
c
c       Choose the step size for super economy mode.
c
        nord = 0
        scale = dlxmx0/delxi
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Choose step size for ordinary economy mode.
c
      if (nord .le. 0) then
        scale = dlxmx0/delxi
        go to 999
      endif
c
      uzt = 'None'
      uzn = ' '
      kcolcn = 0
      scfac = scale
c
      do kcol = 1,kdim
        qmin = kcol.ge.km1 .and. kcol.le.kxt
        d1 = dzvc0(1,kcol)
        d2 = dzvc0(2,kcol)
        mxx0 = zvec0(kcol)
        if (nord.le.1 .or. d2.eq.0.) then
c
c         First order bound.
c
          if (d1 .eq. 0.) go to 100
          if (d1 .gt. 0.) then
            scnst = scnsti
          else
            scnst = scnstd
          endif
          sx = (mxx0*scnst) / (d1*delxi)
          if (iodb(5) .ge. 2) then
            if (sx .lt. scale) then
              rx = sx*delxi
              write (noutpt,1000) uzvec1(kcol),rx
 1000         format(3x,'Element= ',a,/7x,'Root= ',1pe12.5)
              write (noutpt,1010) dzvc0(1,kcol)
 1010         format(7x,'Derivatives= ',4(1pe12.5,2x))
            endif
          endif
        else
c
c         Second order bound. There are two quadratic equations to be
c         solved, each of which has two roots. Only the smallest
c         positive root is of interest. If one root of an equation is
c         complex, the other will also be complex. If both roots are
c         real, they could be positive and negative, both positive, or
c         both negative.
c
          d1sq = d1*d1
          thetai = mxx0*scnsti
          thetad = mxx0*scnstd
          termi = 2.*d2*thetai
          termd = 2.*d2*thetad
c
c         Here ri1, ri2, rd1, and rd2 are values of Xi which are roots
c         of the equations.
c
          aterm = d1sq - termi
          if (aterm .ge. 0.) then
            ri1 = (d1 + sqrt(aterm))/d2
            if (ri1.le.0.) ri1 = 1.e+30
            ri2 = (d1 - sqrt(aterm))/d2
            if (ri2.le.0.) ri2=1.e+30
          else
            ri1 = 1.e+30
            ri2 = 1.e+30
          endif
c
          aterm = d1sq - termd
          if (aterm .ge. 0.) then
            rd1 = (d1 + sqrt(aterm))/d2
            if (rd1.le.0.) rd1 = 1.e+30
            rd2 = (d1 - sqrt(aterm))/d2
            if (rd2.le.0.) rd2 = 1.e+30
          else
            rd1 = 1.e+30
            rd2 = 1.e+30
          endif
c
          rx = min(ri1,ri2,rd1,rd2)
          if (rx .lt. 1.e+30) then
            sx = rx/delxi
          else
            sx = 1.e+30
          endif
c
          if (iodb(5) .ge. 2) then
            if (sx .lt. scale) then
              j2 = ilnobl(uzvec1(kcol))
              write (noutpt,1000) uzvec1(kcol)(1:j2),rx
              write (noutpt,1010) dzvc0(1,kcol),dzvc0(2,kcol)
            endif
          endif
        endif
c
        if (sx .lt. scfac) then
          scfac = sx
          kcolcn = kcol
          uzt(1:24) = ubas(1:24)
        endif
  100   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (kcolcn .gt. 0) uzn(1:24) = uzvec1(kcolcn)(1:24)
      if (iodb(5) .ge. 1) then
        write (noutpt,1020) scfac
 1020   format('   Economy mode: scale= ',f10.5)
        j2 = ilnobl(uzt)
        j3 = ilnobl(uzn)
        if (j3 .gt. 0) then
          write (noutpt,1030) scfac,uzt(1:j2),uzn(1:j3)
 1030     format(5x,'Constrained by ',a,1x,a)
        else
          write (noutpt,1040) uzt(1:j2)
 1040     format(5x,'Constrained by ',a)
        endif
      endif
      scale = min(scalim,scfac)
c
  999 continue
      end
