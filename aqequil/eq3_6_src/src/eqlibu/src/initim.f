      subroutine initim(iexec0,jexec0,texec0,noutpt,nttyo,
     $ udate0,utime0)
c
c     This subroutine gets the time and date of the start of execution.
c     It also calculates the data in the form needed to compute the
c     actual run time at the end of execution (see EQLIBU/runtim.f).
c     This routine should be called at the start of execution.
c
c     This subroutine is called by:
c
c       EQPT/eqpt.f
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       noutpt = unit number of the output file
c       nttyo  = unit number of the screen file
c
c     Output:
c
c       texec0 = time at the start of execution, since midnight,
c                  in seconds
c       iexec0 = day of the year at the start of execution
c       jexec0 = the last two digits of the year at the start of
c                  execution
c       udate0 = the date at the start of execution
c       utime0 = the time at the start of execution
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer noutpt,nttyo
c
      integer iexec0,jexec0
c
      character*11 utime0
      character*9 udate0
c
      real*8 texec0
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ndays(12)
c
      integer mo,mm,dd,yy,hr,hs,mi,se
      integer j2,j3
c
      integer ilnobl
c
      logical qleapy
c
      character*8 umonth(12)
      character*8 ux
      character*3 um
c
c-----------------------------------------------------------------------
c
      data (umonth(mo), mo = 1,12) /'Jan','Feb','Mar','Apr','May',
     $ 'Jun','Jul','Aug','Sep','Oct','Nov','Dec'/
c
      data ndays(1) /31/
      data (ndays(mo), mo = 3,12) /31,30,31,30,31,31,30,31,30,31/
c
c-----------------------------------------------------------------------
c
c     Get the time and date.
c
c     Calling sequence substitutions:
c       udate0 for udate
c       utime0 for utime
c
      call timdat(udate0,utime0)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the data needed to calculate the actual run time (this
c     will be used by EQLIBU/runtim.f).
c
      ux = udate0(1:2)
      read (ux,1000,err=900) dd
 1000 format(i2)
c
      um = udate0(3:5)
      do mo = 1,12
        if (umonth(mo)(1:3) .eq. um(1:3)) then
          mm = mo
          go to 110
        endif
      enddo
c
      j2 = ilnobl(um)
      j3 = ilnobl(udate0)
      write (noutpt,1010) um(1:j2),udate0(1:j3)
      write (nttyo,1010) um(1:j2),udate0(1:j3)
 1010 format(/' * Error - (EQLIBU/initim) The string "',a,'"',
     $ " doesn't match the standard",/7x,'abbreviation for any',
     $ ' month. This was extracted from the string "',a,'",',
     $ /7x,'which should contain the date at the start of',
     $ ' execution. Check',/7x,'porting modifications in',
     $ ' EQLIBU/timdat.f.')
      stop
c
  110 if (udate0(8:9) .ne. '  ') then
        ux = udate0(6:9)
        read (ux,1020,err=900) yy
 1020   format(i4)
      else
        ux = udate0(6:7)
        read (ux,1000,err=900) yy
      endif
c
      ux = utime0(1:2)
      read (ux,1000,err=910) hr
c
      ux = utime0(4:5)
      read (ux,1000,err=910) mi
c
      ux = utime0(7:8)
      read (ux,1000,err=910) se
c
      ux = utime0(10:11)
      read (ux,1000,err=910) hs
c
      texec0 = hr*3600. + mi*60. + se*1. + hs*0.01
c
c     Get the number of days in February of the current year.
c
      ndays(2) = 28
      call tleapy(yy,qleapy)
      if (qleapy) ndays(2) = 29
c
      iexec0 = dd
      if (mm .gt. 1) then
        do mo = 1, mm - 1
          iexec0 = iexec0 + ndays(mo)
        enddo
      endif
      jexec0 = yy
c
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  900 j3 = ilnobl(udate0)
      write (noutpt,1030) udate0(1:j3)
      write (nttyo,1030) udate0(1:j3)
 1030 format(/' * Error - (EQLIBU/initim) The date at the start of',
     $ ' execution (udate0)'/7x,'must be a string in the format "',
     $ 'ddMmmccyy" (e.g., 23Sep1997) or',/7x,'"ddmmyy", (e.g.,',
     $ ' 23Sep97). The actual string "',a,'"'," isn't",/7x,'in',
     $ ' compliance with this requirement. Check porting',
     $ ' modifications.',/7x,'in EQLIBU/timdat.f.')
      stop
c
  910 j3 = ilnobl(utime0)
      write (noutpt,1040) utime0(1:j3)
      write (nttyo,1040) utime0(1:j3)
 1040 format(/' * Error - (EQLIBU/initim) The time at the start of',
     $ ' execution (utime0)'/7x,'must be a string in the format "',
     $ 'hh:mm:ss:hs" (e.g., 15:32:19:28).',7x,'The actual string "',
     $ a,'"'," isn't in compliance with this requirement.",
     $ /7x, 'Check porting modifications in EQLIBU/timdat.f.')
      stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
