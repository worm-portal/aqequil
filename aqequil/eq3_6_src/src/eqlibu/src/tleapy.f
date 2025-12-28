      subroutine tleapy(yy,qleapy)
c
c     This subroutine determines if the year yy is a leap year. The
c     logical flag qleapy is returned as true if this is the case.
c     Modern era rules apply:
c
c       1. Any year not evenly divisible by 4 is not a leap year.
c       2. If a year is evenly divisible by 100 (i.e., is a "century
c          year"), it is not a leap year, unless it is also evenly
c          divisible by 400.
c
c     Examples:
c
c        1. The year 1998 is not a leap year because it is not evenly
c           divisible by 4.
c        2. Century years like 1700, 1800, and 1900 are not leap years
c           because they are not evenly divisible by 400.
c        3. Century years like 1600 and 2000 are leap years because
c           they are evenly divisible by 400.
c
c     The modern era rules (part of the Gregorian calendar) first
c     applied in 1582. The first leap year of the Gregorian calendar
c     was therefore 1584. This subroutine will not return any prior
c     year as a leap year. Thus, 1580 is not a leap year.
c
c     This subroutine is called by:
c
c       EQLIB/initim.f
c       EQLIB/runtim.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       yy     = year, assumed to be an integer up to four digits
c
c     Output:
c
c       qleapy = logical flag, true if yy is a leap year.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer yy
c
      logical qleapy
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
c     None
c
c-----------------------------------------------------------------------
c
      qleapy = .false.
c
      if (yy .ge. 1584) then
c
c       Have a year sufficiently recent to possibly be a leap year.
c
        if (mod(yy,4) .eq. 0) then
c
c         Have a year evenly divisible by 4. It might be a leap year.
c
          if (mod(yy,100) .eq. 0) then
c
c           Have a century year. It might be a leap year.
c
            if (mod(yy,400) .eq. 0) then
c
c             Have a century year that is evenly divisible by 400.
c             This is a leap year.
c
              qleapy = .true.
            endif
          else
c
c           Have a year evenly divisible by 4, but it is not a
c           century year. This is a leap year.
c
            qleapy = .true.
          endif
        endif
      endif
c
      end
