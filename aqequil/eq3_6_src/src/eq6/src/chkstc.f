      subroutine chkstc(actw,awmax,awmin,eh,ehmax,ehmin,fo2lg,iopt,
     $ jreac,kstep,kstpmx,noptmx,noutpt,nrct,nrctmx,nttyo,o2max,o2min,
     $ ph,phmax,phmin,prcinf,qaft1,qcnpre,qcntmp,qconst,qredox,qstop,
     $ qvhfxi,qvlsow,timemx,time1,tolxst,tolxsu,ximax,xi1)
c
c     This subroutine checks for conditions which call for terminating
c     the current reaction path calculation.
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
      integer kstpmx,noptmx,nrctmx
c
      integer noutpt,nttyo
c
      integer iopt(noptmx),jreac(nrctmx)
c
      integer kstep,nrct
c
      logical qaft1,qcnpre,qcntmp,qconst,qredox,qstop,qvhfxi,qvlsow
c
      real*8 actw,awmax,awmin,eh,ehmax,ehmin,fo2lg,o2max,o2min,
     $ ph,phmax,phmin,prcinf,timemx,time1,tolxst,tolxsu,ximax,xi1
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer nrc,nrrct
c
      real*8 tres
c
c-----------------------------------------------------------------------
c
      qstop = .false.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (xi1 .ge. ximax) then
        write (noutpt,1000)
        write (nttyo,1000)
 1000   format(/3x,'Have reached the maximum value of Xi.')
        qstop = .true.
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(2) .ge. 1) then
        if (timemx .lt. prcinf) then
          tres = (time1 - timemx)/timemx
          if (abs(tres) .le. tolxst) then
            write (noutpt,1010)
            write (nttyo,1010)
 1010       format(/3x,'Have reached the maximum value of time.')
            qstop = .true.
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (phmin .gt. -prcinf) then
        if (abs(ph - phmin) .le. tolxsu) then
          write (noutpt,1020)
          write (nttyo,1020)
 1020     format(/3x,'Have reached the minimum value of pH.')
          qstop = .true.
        elseif (ph .lt. phmin) then
          write (noutpt,1030)
          write (nttyo,1030)
 1030     format(/3x,'The pH is below the minimum value.')
          qstop = .true.
        endif
      endif
c
      if (phmax .lt. prcinf) then
        if (abs(ph - phmax) .le. tolxsu) then
          write (noutpt,1040)
          write (nttyo,1040)
 1040     format(/3x,'Have reached the maximum value of pH.')
          qstop = .true.
        elseif (ph .gt. phmax) then
          write (noutpt,1050)
          write (nttyo,1050)
 1050     format(/3x,'The pH is above the maximum value.')
          qstop = .true.
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qredox) then
        if (ehmin .gt. -prcinf) then
          if (abs(eh - ehmin) .le. tolxsu) then
            write (noutpt,1060)
            write (nttyo,1060)
 1060       format(/3x,'Have reached the minimum value of Eh (v).')
            qstop = .true.
          elseif (eh .lt. ehmin) then
            write (noutpt,1070)
            write (nttyo,1070)
 1070       format(/3x,'The Eh (v) is below the minimum value.')
            qstop = .true.
          endif
        endif
c
        if (ehmax .lt. prcinf) then
          if (abs(eh - ehmax) .le. tolxsu) then
            write (noutpt,1080)
            write (nttyo,1080)
 1080       format(/3x,'Have reached the maximum value of Eh (v).')
            qstop = .true.
          elseif (eh .gt. ehmax) then
            write (noutpt,1090)
            write (nttyo,1090)
 1090       format(/3x,'The Eh (v) is above the maximum value.')
            qstop = .true.
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qredox) then
        if (o2min .gt. -prcinf) then
          if (abs(fo2lg - o2min) .le. tolxsu) then
            write (noutpt,1100)
            write (nttyo,1100)
 1100       format(/3x,'Have reached the minimum value of log fO2.')
            qstop = .true.
          elseif (fo2lg .lt. o2min) then
            write (noutpt,1110)
            write (nttyo,1110)
 1110       format(/3x,'The log fO2 is below the minimum value.')
            qstop = .true.
          endif
        endif
c
        if (o2max .lt. prcinf) then
          if (abs(fo2lg - o2max) .le. tolxsu) then
            write (noutpt,1120)
            write (nttyo,1120)
 1120       format(/3x,'Have reached the maximum value of log fO2.')
            qstop = .true.
          elseif (fo2lg .gt. o2max) then
            write (noutpt,1130)
            write (nttyo,1130)
 1130       format(/3x,'The log fO2 is above the maximum value.')
            qstop = .true.
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (awmin .gt. -prcinf) then
        if (abs(actw - awmin) .le. tolxsu) then
          write (noutpt,1140)
          write (nttyo,1140)
 1140     format(/3x,'Have reached the minimum value of the',
     $    ' activity of water.')
          qstop = .true.
        elseif (actw .lt. awmin) then
          write (noutpt,1150)
          write (nttyo,1150)
 1150     format(/3x,'The activity of water is below the minimum',
     $    ' value.')
          qstop = .true.
        endif
      endif
c
      if (awmax .lt. prcinf) then
        if (abs(actw - awmax) .le. tolxsu) then
          write (noutpt,1160)
          write (nttyo,1160)
 1160     format(/3x,'Have reached the maximum value of the',
     $    ' activity of water.')
          qstop = .true.
        elseif (actw .gt. awmax) then
          write (noutpt,1170)
          write (nttyo,1170)
 1170     format(/3x,'The activity of water is above the maximum',
     $    ' value.')
          qstop = .true.
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (kstep .ge. kstpmx) then
        write (noutpt,1200)
        write (nttyo,1200)
 1200   format(/3x,'Have done the maximum number of steps.')
        qstop = .true.
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qconst) then
        write (noutpt,1400)
        write (nttyo,1400)
 1400   format(/3x,'All irreversible reaction rates are now zero.')
        qstop=.true.
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qaft1) then
        write (noutpt,1410)
        write (nttyo,1410)
 1410   format(/3x,'The total affinity is repeatedly nearly zero.')
        qstop=.true.
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Test to see if any reactants remain active. If none, terminate
c     the run unless the temperature or pressure are changing.
c
      if (nrct .gt. 0) then
        if (qcntmp .and. qcnpre) then
          nrrct = 0
          do nrc = 1,nrct
            if (jreac(nrc) .eq. 0) nrrct = nrrct + 1
          enddo
          if (nrrct .le. 0) then
            write (noutpt,1420)
            write (nttyo,1420)
 1420       format(/3x,'Each reactant is now saturated or exhausted.')
            qstop = .true.
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qvlsow) then
        write (noutpt,1430)
        write (nttyo,1430)
 1430   format(/3x,'Have very nearly exhausted solvent water.')
        qstop=.true.
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qvhfxi) then
        write (noutpt,1440)
        write (nttyo,1440)
 1440   format(/3x,'Have extremely high ionic strength.')
        qstop=.true.
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
