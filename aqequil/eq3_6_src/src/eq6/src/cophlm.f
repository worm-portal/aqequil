      subroutine cophlm(actw,awmax,awmin,delxi,dlxmin,eh,ehmax,
     $ ehmin,fo2lg,iodb,nodbmx,noutpt,o2max,o2min,ph,phmax,phmin,
     $ qadjdx,qredox,tolxsu)
c
c     This subroutine checks for oversteps with regard to specified
c     minimum and maximum values for the pH, Eh, log fO2, and
c     activity of water.
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
      integer nodbmx
c
      integer noutpt
c
      integer iodb(nodbmx)
c
      logical qadjdx,qredox
c
      real*8 actw,awmax,awmin,delxi,dlxmin,eh,ehmax,ehmin,fo2lg,
     $ o2max,o2min,ph,phmax,phmin,tolxsu
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
c     None
c
c-----------------------------------------------------------------------
c
      qadjdx = .false.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Is the pH less than the requested minimum value?
c
      if ((ph - phmin) .lt. -tolxsu) then
        if (iodb(1) .gt. 0) then
          write (noutpt,1000) ph,phmin
 1000      format(/3x,'The pH is ',f7.4,', which is less than the',
     $     ' minimum value',/5x,'of ',f7.4,'. This overstep exceeds',
     $     ' the specified tolerance.',/)
        endif
c
c       Determine whether or not to reduce delxi and go back and try
c       again to better locate the event in question.
c
        call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)
        if (qadjdx) go to 999
      endif
c
c     Is the pH greater than the requested maximum value?
c
      if ((ph - phmax) .gt. tolxsu) then
        if (iodb(1) .gt. 0) then
          write (noutpt,1010) ph,phmax
 1010     format(/3x,'The pH is ',f7.4,', which is more than the',
     $    ' maximum value',/5x,'of ',f7.4,'. This overstep exceeds',
     $    ' the specified tolerance.',/)
        endif
c
c       Determine whether or not to reduce delxi and go back and try
c       again to better locate the event in question.
c
        call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)
        if (qadjdx) go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qredox) then
c
c       Is the Eh less than the requested minimum value?
c
        if ((eh - ehmin) .lt. -tolxsu) then
          if (iodb(1) .gt. 0) then
            write (noutpt,1020) eh,ehmin
 1020       format(/3x,'The Eh is ',f9.4,' v, which is less than the',
     $      ' minimum value',/5x,'of ',f9.4,'. This overstep exceeds',
     $      ' the specified tolerance.',/)
          endif
c
c         Determine whether or not to reduce delxi and go back and try
c         again to better locate the event in question.
c
          call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)
          if (qadjdx) go to 999
        endif
c
c       Is the Eh greater than the requested maximum value?
c
        if ((eh - ehmax) .gt. tolxsu) then
          if (iodb(1) .gt. 0) then
            write (noutpt,1040) eh,ehmax
 1040       format(/3x,'The Eh is ',f9.4,' v, which is more than the',
     $      ' maximum value',/5x,'of ',f9.4,'. This overstep exceeds',
     $      ' the specified tolerance.',/)
          endif
c
c         Determine whether or not to reduce delxi and go back and try
c         again to better locate the event in question.
c
          call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)
          if (qadjdx) go to 999
        endif
c
c       Is the log fO2 less than the requested minimum value?
c
        if ((fo2lg - o2min) .lt. -tolxsu) then
          if (iodb(1) .gt. 0) then
            write (noutpt,1050) fo2lg,o2min
 1050       format(/3x,'The log fO2 is ',f9.4,', which is less than',
     $      ' the minimum value',/5x,'of ',f9.4,'. This overstep',
     $      ' exceeds the specified tolerance.',/)
          endif
c
c         Determine whether or not to reduce delxi and go back and try
c         again to better locate the event in question.
c
          call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)
          if (qadjdx) go to 999
        endif
c
c       Is the Eh greater than the requested maximum value?
c
        if ((fo2lg - o2max) .gt. tolxsu) then
          if (iodb(1) .gt. 0) then
            write (noutpt,1060) fo2lg,o2max
 1060       format(/3x,'The log fO2 is ',f9.4,', which is more than',
     $      ' the maximum value',/5x,'of ',f9.4,'. This overstep',
     $      ' exceeds the specified tolerance.',/)
          endif
c
c         Determine whether or not to reduce delxi and go back and try
c         again to better locate the event in question.
c
          call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)
          if (qadjdx) go to 999
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Is the activity of water less than the requested minimum value?
c
      if ((actw - awmin) .lt. -tolxsu) then
        if (iodb(1) .gt. 0) then
          write (noutpt,1070) actw,awmin
 1070     format(/3x,'The activity of water is ',f6.4,', which is',
     $    ' less than the',/5x,'minimum value of ',f6.4,'. This',
     $    ' overstep exceeds the',/5x,'specified tolerance.',/)
        endif
c
c       Determine whether or not to reduce delxi and go back and try
c       again to better locate the event in question.
c
        call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)
        if (qadjdx) go to 999
      endif
c
c     Is the activity of water greater than the requested maximum
c     value?
c
      if ((actw - awmax) .gt. tolxsu) then
        if (iodb(1) .gt. 0) then
          write (noutpt,1080) actw,awmax
 1080     format(/3x,'The activity of water is ',f6.4,', which is',
     $    ' greater than the',/5x,'maximum value of ',f6.4,'. This',
     $    ' overstep exceeds the',/5x,'specified tolerance.',/)
        endif
c
c       Determine whether or not to reduce delxi and go back and try
c       again to better locate the event in question.
c
        call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)
        if (qadjdx) go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
