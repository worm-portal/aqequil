      subroutine timdat(udate,utime)
c
c     This subroutine gets the date and clock time in ASCII for time and
c     date stamping on output files. The variable "udate" is returned in
c     the form "ddMmmccyy" (e.g., 23Sep1997) under Fortran 90, in the
c     form 'ddMmmyy" (e.g., 23Sep97) under Fortran 77. The variable
c     "utime" is returned in the form "hh:mm:ss:hs", though in some
c     cases the hundred-th second part ("hs") may be blank.
c
c     See also EQLIBU/initim.f and EQLIBU/runtim.f.
c
c     This subroutine is called by:
c
c       EQLIBU/initim.f
c       EQLIBU/runtim.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       None
c
c     Output:
c
c       udate  = string containing the date
c       utime  = string containing the clock time
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      character*11 utime
      character*9 udate
c
c-----------------------------------------------------------------------
c
        integer values(8)
c
        character*8 date
        character*10 time
        character*6 zone
        character*8 ux
        character*8 umonth(12)
c
        integer mo,mm
c
        data (umonth(mo), mo = 1,12) /'Jan','Feb','Mar','Apr','May',
     $  'Jun','Jul','Aug','Sep','Oct','Nov','Dec'/
c
c-----------------------------------------------------------------------
c
c       Note: date_and_time is a Fortran 90 intrinsic subroutine.
c
        call date_and_time(date,time,zone,values)
c
c       Note: date has the form ccyymmdd; udate has the form ddMmmccyy.
c
        udate(1:2) = date(7:8)
        ux = date(5:6)
        read (ux,'(i5)') mm
        udate(3:5) = umonth(mm)(1:3)
        udate(6:9) = date(1:4)
c
c       Note: time has the form hhmmss.sss (sss = milliseconds);
c       utime has the form hh:mm:ss:hs.
c
        utime(1:2) = time(1:2)
        utime(3:3) = ':'
        utime(4:5) = time(3:4)
        utime(6:6) = ':'
        utime(7:8) = time(5:6)
        utime(9:9) = ':'
        utime(10:11) = time(8:9)
c
c       The following is nonsense intended to keep the compiler from
c       complaining that zone and values are not used.
c
        ux = zone
        mm = values(1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
