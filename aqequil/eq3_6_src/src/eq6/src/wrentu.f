      subroutine wrentu(actw,eh,fo2lg,iopg,iopt,kstep,nopgmx,noptmx,
     $ nttyo,ph,qredox,time1,xi1)
c
c     This subroutine writes entertainment for the user while the
c     run is underway. This output gives an idea of how the
c     calculation is progressing. In that respect, it is somewhat
c     like a progress bar.
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
      integer nopgmx,noptmx
c
      integer nttyo
c
      integer iopg(nopgmx),iopt(noptmx)
c
      integer kstep
c
      logical qredox
c
      real*8 actw,eh,fo2lg,ph,time1,xi1
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
c     None
c
c-----------------------------------------------------------------------
c
      if (iopt(2) .le. 0) then
        if (qredox) then
          if (iopg(1) .gt. 0) then
            write (nttyo,1000) kstep,xi1,ph,actw,fo2lg
 1000       format(3x,'step=',i5,', Xi=',1pe10.3,', pH=',0pf7.3,
     $      ', aw=',f6.3,', log fO2= ',f7.3)
          else
            write (nttyo,1010) kstep,xi1,ph,fo2lg
 1010       format(3x,'step=',i5,', Xi=',1pe10.3,', pH=',0pf7.3,
     $      ', log fO2= ',f7.3)
          endif
        else
          if (iopg(1) .gt. 0) then
            write (nttyo,1020) kstep,xi1,ph,actw
 1020       format(3x,'step=',i5,', Xi=',1pe10.3,', pH=',0pf7.3,
     $      ', aw=',f6.3)
          else
            write (nttyo,1030) kstep,xi1,ph
 1030       format(3x,'step=',i5,', Xi=',1pe10.3,', pH=',0pf7.3)
          endif
        endif
      else
        if (qredox) then
          if (iopg(1) .gt. 0) then
            write (nttyo,1100) kstep,time1,ph,actw,fo2lg
 1100       format(3x,'step=',i5,', time=',1pe10.3,' s, pH=',
     $      0pf7.3,', aw=',f6.3,', log fO2= ',f7.3)
          else
            write (nttyo,1110) kstep,time1,ph,fo2lg
 1110       format(3x,'step=',i5,', time=',1pe10.3,' s, pH=',
     $      0pf7.3,', log fO2= ',f7.3)
          endif
        else
          if (iopg(1) .gt. 0) then
            write (nttyo,1120) kstep,time1,ph,actw
 1120       format(3x,'step=',i5,', time=',1pe10.3,' s, pH=',
     $      0pf7.3,', aw=',f6.3)
          else
            write (nttyo,1130) kstep,time1,ph
 1130       format(3x,'step=',i5,', time=',1pe10.3,' s, pH=',
     $      0pf7.3)
          endif
        endif
      endif
c
      end
