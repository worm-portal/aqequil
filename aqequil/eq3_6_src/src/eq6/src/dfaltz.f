      subroutine dfaltz(dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,
     $ dloplo,dloprn,dltplo,dltpll,dltprl,dltprn,dlxdmp,dlxmx0,
     $ dlxplo,dlxpll,dlxprl,dlxprn,iopt,itermx,ksplmx,ksppmx,kstpmx,
     $ net,noptmx,nordmx,nrct,noutpt,ntrymx,nttyo,prcinf,qecon,qscon,
     $ timmxi,tistti,tolbt,toldl,tolsat,tolxsf,tolxst,tolxsu,
     $ ximaxi,xistti)
c
c     This subroutine sets the defaults for various run parameters read
c     from the input file. In some cases, it forces the parameters to
c     take on certain values, or to fall in certain ranges.
c
c     This subroutine is called by:
c
c       EQ6/eq6.f
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
      integer noptmx
c
      integer noutpt,nttyo
c
      integer iopt(noptmx)
c
      integer itermx,ksplmx,ksppmx,kstpmx,net,nordmx,nrct,ntrymx
c
      logical qecon,qscon
c
      real*8 dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,
     $ dltplo,dltpll,dltprl,dltprn,dlxdmp,dlxmx0,dlxplo,dlxpll,
     $ dlxprl,dlxprn,prcinf,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,
     $ tolxst,tolxsu,ximaxi,xistti
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
c       None
c
c-----------------------------------------------------------------------
c
c     Iopt(1): Physical system model option.
c
      if (iopt(1) .lt. 0) then
        write (noutpt,1000) iopt(1)
        write (nttyo,1000) iopt(1)
 1000   format(/' * Note - (EQ6/dfaltz) The value of iopt(1) is ',i3,
     $  ',',/7x,'which is out of range. Resetting iopt(1) to 0.')
        iopt(1) = 0
      endif
c
      if (iopt(1) .gt. 2) then
        write (noutpt,1010) iopt(1)
        write (nttyo,1010) iopt(1)
 1010   format(/' * Note - (EQ6/dfaltz) The value of iopt(1) is ',i3,
     $  ',',/7x,'which is out of range. Resetting iopt(1) to 2.')
        iopt(1) = 2
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iopt(3): Phase boundary search option.
c
      if (iopt(3).gt.0 .and. iopt(1).ge.2) then
        write (noutpt,1020) iopt(3),iopt(1)
        write (nttyo,1020) iopt(3),iopt(1)
 1020   format(/' * Note - (EQ6/dfaltz) The value of iopt(3) is ',i3,
     $  ',',/7x,'which is not valid because iopt(1) is ',i3,'.',
     $  /7x,'Resetting iopt(3) = 0.')
        iopt(3) = 0
      endif
c
      if (iopt(3).gt.0 .and. iopt(2).ge.1) then
        write (noutpt,1030) iopt(3),iopt(2)
        write (nttyo,1030) iopt(3),iopt(2)
 1030   format(/' * Note - (EQ6/dfaltz) The value of iopt(3) is ',i3,
     $  ',',/7x,'which is not valid because iopt(2) is ',i3,'.',
     $  /7x,'Resetting iopt(3) = 0.')
        iopt(3) = 0
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum value of Xi.
c
      if (ximaxi .lt. xistti) then
        write (noutpt,1040)
        write (nttyo,1040)
 1040   format(/' * Note - (EQ6/dfaltz) The maximum value of reaction',
     $  /7x,'progress (ximaxi) is less than the starting value',
     $  ' (xistti).',/7x,'The maximum value is being reset to',
     $  ' infinity.')
        ximaxi = prcinf
      endif
      if (ximaxi .gt. prcinf) ximaxi = prcinf
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum value of time.
c
      if (timmxi .lt. tistti) then
        write (noutpt,1050)
        write (nttyo,1050)
 1050   format(/' * Note - (EQ6/dfaltz) The maximum value of time',
     $  /7x,'(timmxi) is less than the starting value (tistti).',
     $  /7x,'The maximum value is being reset to infinity.')
        timmxi = prcinf
      endif
      if (timmxi .gt. prcinf) timmxi = prcinf
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum number of steps.
c
      if (kstpmx .lt. 0) kstpmx = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Newton-Raphson beta tolerance.
c
      if (iopt(2) .le. 0) then
        if (tolbt .le. 0.) tolbt = 1.e-6
        if (tolbt .lt. 1.e-10) tolbt = 1.e-10
        if (tolbt .gt. 1.e-4) tolbt = 1.e-4
      else
        if (tolbt .le. 0.) tolbt = 1.e-8
        if (tolbt .lt. 1.e-10) tolbt = 1.e-10
        if (tolbt .gt. 1.e-6) tolbt = 1.e-6
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Newton-Raphson del tolerance.
c
      if (iopt(2) .le. 0) then
        if (toldl .le. 0.) toldl = 1.e-6
        if (toldl .lt. 1.e-10) toldl = 1.e-10
        if (toldl .gt. 1.e-4) toldl = 1.e-4
      else
        if (toldl .le. 0.) toldl = 1.e-8
        if (toldl .lt. 1.e-10) toldl = 1.e-10
        if (toldl .gt. 1.e-6) toldl = 1.e-6
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum finite-difference order.
c
      if (nordmx .eq. 0) nordmx = 6
      if (nordmx .lt. 1) nordmx = 1
      if (nordmx .gt. 8) nordmx = 8
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum number of Newton-Raphson iterations.
c
      if (itermx .le. 0) itermx = 200
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum number of tries to find the correct phase assemblage.
c
      if (ntrymx .le. 0) ntrymx = 100
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Ordinary search/find tolerance.
c
      if (tolxsf .le. 0.) tolxsf = 1.e-6
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Search tolerance on time (relative).
c
      if (tolxst .le. 0.) tolxst = 1.e-8
      if (tolxst .lt. 1.e-12) tolxst = 1.e-12
      if (tolxst .gt. 1.e-4) tolxst = 1.e-4
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Search tolerance on pH, Eh, log fO2, and activity of water
c     (linear).
c
      if (tolxsu .le. 0.) tolxsu = 1.e-5
      if (tolxsu .lt. 1.e-8) tolxsu = 1.e-8
      if (tolxsu .gt. 1.e-4) tolxsu = 1.e-4
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Saturation tolerance.
c
      if (tolsat .le. 0.) tolsat = 0.0005
      if (tolsat .lt. 0.00005) tolsat = 0.00005
      if (tolsat .gt. 0.05) tolsat = 0.05
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Limits on Newton-Raphson delta tolerance.
c
      if (tolxsf .lt. 1.e-10) toldl = 1.e-10
      if (tolxsf .gt. 1.e-2) toldl = 1.e-2
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Zero-order step size.
c
      if (dlxmx0 .le. 0.) then
        dlxmx0 = 1.e-8
        if (iopt(2) .ge. 1) dlxmx0 = 1.e-9
        if (iopt(1) .ge. 2) dlxmx0 = 1.e-9
        if (qecon) dlxmx0 = 1.e-6
        if (nrct .le. 0) dlxmx0 = 1.e-2
        if (qscon) then
          if (dlxprn.gt.0.) then
            dlxmx0 = dlxprn
          elseif (ximaxi .gt. xistti) then
            dlxmx0 = 0.1*(ximaxi - xistti)
          endif
        endif
      endif
c
      if (dlxmx0 .lt. 1.e-12) dlxmx0 = 1.e-12
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print interval in Xi.
c
      if (dlxprn .le. 0.) dlxprn = prcinf
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print interval in log Xi.
c
      if (dlxprl .le. 0.) dlxprl = 0.5
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print interval in time (seconds).
c
      if (dltprn .le. 0.) dltprn = prcinf
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print interval in log time.
c
      if (dltprl .le. 0.) dltprl = prcinf
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print interval in pH.
c
      if (dlhprn .le. 0.) dlhprn = prcinf
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print interval in Eh (v).
c
      if (dleprn .le. 0.) dleprn = prcinf
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print interval in log fO2.
c
      if (dloprn .le. 0.) dloprn = prcinf
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print interval in activity of water.
c
      if (dlaprn .le. 0.) dlaprn = prcinf
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print interval in steps.
c
      if (ksppmx .le. 0) ksppmx = 100
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Plot interval in Xi.
c
      if (dlxplo .le. 0.) dlxplo = prcinf
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Plot interval in log Xi.
c
      if (dlxpll .le. 0.) dlxpll = prcinf
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Plot interval in pH.
c
      if (dlhplo .le. 0.) dlhplo = prcinf
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Plot interval in Eh (v).
c
      if (dleplo .le. 0.) dleplo = prcinf
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Plot interval in log fO2.
c
      if (dloplo .le. 0.) dloplo = prcinf
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Plot interval in activity of water.
c
      if (dlaplo .le. 0.) dlaplo = prcinf
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Plot interval in time (seconds).
c
      if (dltplo .le. 0.) dltplo = prcinf
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Plot interval in log time.
c
      if (dltpll .le. 0.) dltpll = prcinf
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Plot interval in steps.
c
      if (ksplmx .le. 0.) ksplmx = 10000
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     PRS transfer interval in Xi.
c
      if (dlxdmp .le. 0.) then
        if (iopt(1) .le. 1) then
c
c         Have the titration model or the closed system model.
c
          dlxdmp = prcinf
        elseif (iopt(1) .eq. 2) then
c
c         Have the fluid-centered, flow-through open system model.
c
          dlxdmp = prcinf
          if (iopt(4).ge.1 .or. net.gt.0) then
            if (ximaxi .gt. 0.) then
              dlxdmp = ximaxi/25.
            else
              dlxdmp = 1.e-6
            endif
          endif
        endif
      endif
c
      if (iopt(1) .eq. 2) then
        if (iopt(4).ge.1 .or. net.gt.0) then
c
c         Have solid solutions or generic ion exchangers, or both.
c         Make sure the PRS transfer interval is no greater then the
c         linear plot interval.
c
          if (dlxdmp .gt. dlxplo) dlxdmp = dlxplo
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
