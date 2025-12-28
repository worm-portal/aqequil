      subroutine setup(coval,eh,ehfac,ier,irdxc3,itdsf3,jflag,mwtsp,
     $ narn1,nbaspd,nbt,nbtmax,noutpt,nstmax,nttyo,pe,rho,tdspkg,
     $ tdspl,uspec)
c
c     This subroutine converts input coval data which are not on the
c     molal concentration scale to that scale.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
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
      integer nbtmax,nstmax
c
      integer noutpt,nttyo
c
      integer jflag(nstmax),nbaspd(nbtmax)
c
      integer ier,irdxc3,itdsf3,narn1,nbt
c
      character(len=48) uspec(nstmax)
c
      real(8) coval(nbtmax),mwtsp(nstmax)
c
      real(8) eh,ehfac,pe,rho,tdspkg,tdspl
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jfl,jlen,nb,ns
c
      character(len=56) uspn56
c
      real(8) wfh2o
c
c-----------------------------------------------------------------------
c
      ier = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (itdsf3 .ge. 1) tdspkg = tdspl/rho
      wfh2o = (1000. - 0.001*tdspkg)/1000.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (irdxc3 .eq. -2) then
c
c       Convert input pe to Eh.
c
        eh = pe*ehfac
        irdxc3 = -1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do nb = 1,nbt
        ns = nbaspd(nb)
        jfl = jflag(ns)
        if (ns .eq. narn1) go to 100
        if (jfl .eq. 0) go to 100
        if (jfl .eq. 30) go to 100
        if (jfl .eq. 7) go to 100
        if (jfl.eq.16 .or. jfl.eq.17 .or. jfl.eq.18) go to 100
        if (jfl.eq.19 .or. jfl.eq.20 .or. jfl.eq.21) go to 100
        if (jfl.eq.22 .or. jfl.eq.23) go to 100
        if (jfl.eq.25 .or. jfl.eq.27) go to 100
        if (jfl .eq.- 1) go to 100
c
        if (jfl .eq. 1) then
c
c         Convert molarity to molality.
c
          coval(nb) = coval(nb)/(rho*wfh2o)
          jflag(ns) = 0
        elseif (jfl .eq. 2) then
c
c         Convert mg/L to molality.
c
          coval(nb) = coval(nb)*1.e-3/(mwtsp(ns)*rho*wfh2o)
          jflag(ns) = 0
        elseif (jfl .eq. 3) then
c
c         Convert mg/kg.sol to molality.
c
          coval(nb) = coval(nb)*1.e-3/(mwtsp(ns)*wfh2o)
          jflag(ns) = 0
c
        elseif (jfl .eq. 8) then
c
c         Convert alkalinity in eq/L to eq/kg.H2O.
c
          coval(nb) = coval(nb)/(rho*wfh2o)
          jflag(ns) = 7
        elseif (jfl .eq. 9) then
c
c         Convert alkalinity in eq/kg.sol to eq/kg.H2O.
c
          coval(nb) = coval(nb)/wfh2o
          jflag(ns) = 7
        elseif (jfl .eq. 10) then
c
c         Convert alkalinity in mg/L CaCO3 to eq/kg.H2O.
c
          coval(nb) = coval(nb)/(50000.*rho*wfh2o)
          jflag(ns) = 7
        elseif (jfl .eq. 11) then
c
c         Convert alkalinity in mg/L HCO3- to eq/kg.H2O.
c
          coval(nb) = coval(nb)/(60960*rho*wfh2o)
          jflag(ns) = 7
        else
c
c         No match.
c
c         Calling sequence substitutions:
c           uspec(ns) for unam48
c
          call fmspnx(jlen,uspec(ns),uspn56)
          write (noutpt,1000) jfl,uspn56(1:jlen)
          write (nttyo,1000) jfl,uspn56(1:jlen)
 1000     format(/' * Error - (EQ3NR/setup) An undefined jflag value',
     $    ' of ',i3,/7x,'was specified on the input file for the',
     $    ' species',/7x,a,'.')
          ier = 1
        endif
  100   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
