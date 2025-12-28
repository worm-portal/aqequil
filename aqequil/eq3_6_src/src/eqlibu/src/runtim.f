      subroutine runtim(iexec0,jexec0,texec0,noutpt,nttyo,trun,
     $ tuser,tcpu,udate1,utime1)
c
c     This subroutine gets the time and date at the end of execution,
c     and also the run, user, and cpu times. The user and cpu times
c     are not defined on all machines. In such cases, zero values are
c     returned. In order to use this subroutine, EQLIBU/initim.f must
c     be called at the start of execution to get values for iexec0,
c     jexec0, and texec0. This routine should be called at the end of
c     execution.
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
c       iexec0 = day of the year at the start of execution
c       jexec0 = the last two digits of the year at the start of
c                  execution
c       texec0 = time at the start of execution, since midnight,
c                in seconds
c       noutpt = unit number of the output file
c       nttyo  = unit number of the screen file
c
c     Output:
c
c       tcpu   = cpu time, in seconds
c       trun   = run time, in seconds
c       tuser  = user time, in seconds
c       udate1 = the date at the end of execution
c       utime1 = the time at the end of execution
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
      character*11 utime1
      character*9 udate1
c
      real*8 texec0,trun,tuser,tcpu
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ndays(12)
c
      integer idays,iexec1,iy,jexec1,mo,mm,dd,yy,hr,hs,mi,se
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
      real*8 texec1
c
c-----------------------------------------------------------------------
c
com   BEGIN_UNIX_DEPENDENT_CODE (Optional)
com
c       real*4 timear(2),tx
com
c       real*8 etime
com
c       external etime
com
com   END_UNIX_DEPENDENT_CODE
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
c       udate1 for udate
c       utime1 for utime
c
      call timdat(udate1,utime1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the data needed to calculate the actual run time.
c
      ux = udate1(1:2)
      read (ux,1000,err=900) dd
 1000 format(i2)
c
      um = udate1(3:5)
      do mo = 1,12
        if (umonth(mo)(1:3) .eq. um(1:3)) then
          mm = mo
          go to 110
        endif
      enddo
c
      j2 = ilnobl(um)
      j3 = ilnobl(udate1)
      write (noutpt,1010) um(1:j2),udate1(1:j3)
      write (nttyo,1010) um(1:j2),udate1(1:j3)
 1010 format(/' * Error - (EQLIBU/runtim) The string "',a,'"',
     $ " doesn't match the standard",/7x,'abbreviation for any',
     $ ' month. This was extracted from the string "',a,'",',
     $ /7x,'which should contain the date at the end of',
     $ ' execution. Check',/7x,'porting modifications in',
     $ ' EQLIBU/timdat.f.')
      stop
c
  110 if (udate1(8:9) .ne. '  ') then
        ux = udate1(6:9)
        read (ux,1020,err=900) yy
 1020   format(i4)
      else
        ux = udate1(6:7)
        read (ux,1000,err=900) yy
      endif
c
      ux = utime1(1:2)
      read (ux,1000,err=910) hr
c
      ux = utime1(4:5)
      read (ux,1000,err=910) mi
c
      ux = utime1(7:8)
      read (ux,1000,err=910) se
c
      ux = utime1(10:11)
      read (ux,1000,err=910) hs
c
      texec1 = hr*3600. + mi*60. + se*1. + hs*0.01
c
c     Get the number of days in February of the current year.
c
      ndays(2) = 28
      call tleapy(yy,qleapy)
      if (qleapy) ndays(2) = 29
c
      iexec1 = dd
      if (mm .gt. 1) then
        do mo = 1, mm - 1
          iexec1 = iexec1 + ndays(mo)
        enddo
      endif
      jexec1 = yy
c
      idays = iexec1 - iexec0
      if (jexec1 .gt. jexec0) then
        do iy = jexec0, jexec1 - 1
          idays = idays + 365
          if (mod(iy,4) .eq. 4) then
            if (mod(iy,100) .ne. 0) then
              idays = idays + 1
            endif
          endif
        enddo
      endif
      trun = idays*86400. + (texec1 - texec0)
c
c     Normally, actual values are not returned for the user and cpu
c     times. These times are not defined on all systems. The user time
c     is usually not useful. The cpu time gives an idea of how much
c     time a code is actually carrying out calculations, as opposed to
c     executing other functions such as I/O.
c
      tuser = 0.
      tcpu = 0.
c
com   BEGIN_UNIX_DEPENDENT_CODE (Optional)
com
com     Activate this block only if you have a UNIX system and want
com     to obtain the user and cpu times on the output files.
com
com     Note: the variable declarations in the corresponding block
com     above must be activated in order to successfully activate
com     the present block.
com
com     The subroutine ETIME called below is a UNIX subroutine that
com     gets the elapsed user and cpu times, returning the former in
com     TIMEAR(1), the latter in TIMEAR(2). The array TIMEAR is real*4.
com     On some systems, the name of this subroutine may be "ETIME_"
com     instead of "ETIME".
com
c       tx = etime(timear)
c       tuser = timear(1)
c       tcpu = timear(2)
com
com     Note: the following statement doesn't really do anything except
com     cause the compiler not to complain that tx is not used.
com
c       timear(1) = tx
com
com   END_UNIX_DEPENDENT_CODE
c
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  900 j3 = ilnobl(udate1)
      write (noutpt,1030) udate1(1:j3)
      write (nttyo,1030) udate1(1:j3)
 1030 format(/' * Error - (EQLIBU/runtim) The date at the end of',
     $ ' execution (udate1)'/7x,'must be a string in the format "',
     $ 'ddMmmccyy" (e.g., 23Sep1997) or',/7x,'"ddmmyy", (e.g.,',
     $ ' 23Sep97). The actual string "',a,'"'," isn't",/7x,'in',
     $ ' compliance with this requirement. Check porting',
     $ ' modifications.',/7x,'in EQLIBU/timdat.f.')
      stop
c
  910 j3 = ilnobl(utime1)
      write (noutpt,1040) utime1(1:j3)
      write (nttyo,1040) utime1(1:j3)
 1040 format(/' * Error - (EQLIBU/runtim) The time at the end of',
     $ ' execution (utime1)'/7x,'must be a string in the format "',
     $ 'hh:mm:ss:hs" (e.g., 15:32:19:28).',7x,'The actual string "',
     $ a,'"'," isn't in compliance with this requirement.",
     $ /7x, 'Check porting modifications in EQLIBU/timdat.f.')
      stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
