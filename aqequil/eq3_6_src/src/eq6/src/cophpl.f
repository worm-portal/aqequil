      subroutine cophpl(actw,aw0plo,aw1plo,delxi,dlxmin,eh,eh0plo,
     $ eh1plo,fo2lg,iodb,nodbmx,noutpt,o20plo,o21plo,ph,ph0plo,
     $ ph1plo,qadjdx,qredox,tolxsu)
c
c     This subroutine checks for oversteps with regard to currently
c     defined lesser and greater plot point values for the pH, Eh,
c     log fO2, and activity of water.
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
      real*8 actw,aw0plo,aw1plo,delxi,dlxmin,eh,eh0plo,eh1plo,fo2lg,
     $ o20plo,o21plo,ph,ph0plo,ph1plo,tolxsu
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
c     Is the pH less than the current lesser plot point value?
c
      if ((ph - ph0plo) .lt. -tolxsu) then
        if (iodb(1) .gt. 0) then
          write (noutpt,1000) ph,ph0plo
 1000     format(/3x,'The pH is ',f7.4,', which is less than the',
     $    ' current lesser plot point',/5x,'value of ',f7.4,'. This',
     $    ' overstep exceeds the specified tolerance.',/)
        endif
c
c       Determine whether or not to reduce delxi and go back and try
c       again to better locate the event in question.
c
        call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)
        if (qadjdx) go to 999
      endif
c
c     Is the pH greater than the current greater plot point value?
c
      if ((ph - ph1plo) .gt. tolxsu) then
        if (iodb(1) .gt. 0) then
          write (noutpt,1010) ph,ph1plo
 1010     format(/3x,'The pH is ',f7.4,', which is greater than the',
     $    ' current greater plot point',/5x,'value of ',f7.4,'. This',
     $    ' overstep exceeds the specified tolerance.',/)
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
c       Is the Eh less than the current lesser plot point value?
c
        if ((eh - eh0plo) .lt. -tolxsu) then
          if (iodb(1) .gt. 0) then
            write (noutpt,1020) eh,eh0plo
 1020       format(/3x,'The Eh is ',f7.4,' v, which is less than the',
     $      ' current lesser plot point',/5x,'value of ',f7.4,' v.',
     $      ' This overstep exceeds the specified tolerance.',/)
          endif
c
c         Determine whether or not to reduce delxi and go back and try
c         again to better locate the event in question.
c
          call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)
          if (qadjdx) go to 999
        endif
c
c       Is the Eh greater than the current greater plot point value?
c
        if ((eh - eh1plo) .gt. tolxsu) then
          if (iodb(1) .gt. 0) then
            write (noutpt,1030) eh,eh1plo
 1030        format(/3x,'The Eh is ',f7.4,' v, which is greater than',
     $       ' the current greater plot point',/5x,'value of ',f7.4,
     $       ' v. This overstep exceeds the specified tolerance.',/)
          endif
c
c         Determine whether or not to reduce delxi and go back and try
c         again to better locate the event in question.
c
          call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)
          if (qadjdx) go to 999
        endif
c
c       Is the log fO2 less than the current lesser plot point value?
c
        if ((fo2lg - o20plo) .lt. -tolxsu) then
          if (iodb(1) .gt. 0) then
            write (noutpt,1040) fo2lg,o20plo
 1040       format(/3x,'The log fO2 is ',f7.4,', which is less than',
     $      ' the current lesser plot point',/5x,'value of ',f7.4,
     $      '. This overstep exceeds the specified tolerance.',/)
          endif
c
c         Determine whether or not to reduce delxi and go back and try
c         again to better locate the event in question.
c
          call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)
          if (qadjdx) go to 999
        endif
c
c       Is the log fO2 greater than the current greater plot point
c       value?
c
        if ((fo2lg - o21plo) .gt. tolxsu) then
          if (iodb(1) .gt. 0) then
            write (noutpt,1050) fo2lg,o21plo
 1050        format(/3x,'The log fO2 is ',f7.4,', which is greater',
     $       ' than the current greater plot',/5x,'point value of ',
     $       f7.4,'. This overstep exceeds the specified tolerance.',/)
          endif
c
c         Determine whether or not to reduce delxi and go back and try
c         again to better locate the event in question.
c
          call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)
          if (qadjdx) go to 999
        endif
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Is the activity of water less than the current lesser plot
c     point value?
c
      if ((actw - aw0plo) .lt. -tolxsu) then
        if (iodb(1) .gt. 0) then
          write (noutpt,1060) actw,aw0plo
 1060     format(/3x,'The activity of water is ',f7.4,', which is less',
     $    ' than the current lesser',/5x,'plot point value of ',f7.4,
     $    '. This overstep exceeds the specified',/5x,'tolerance.',/)
        endif
c
c       Determine whether or not to reduce delxi and go back and try
c       again to better locate the event in question.
c
        call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)
        if (qadjdx) go to 999
      endif
c
c     Is the activity of water greater than the current greater plot
c     point value?
c
      if ((actw - aw1plo) .gt. tolxsu) then
        if (iodb(1) .gt. 0) then
          write (noutpt,1070) actw,aw1plo
 1070     format(/3x,'The activity of water is ',f7.4,', which is',
     $    ' greater than the current',/5x,'greater plot point value',
     $    ' of ',f7.4,'. This overstep exceeds the',/5x,'specified',
     $    ' tolerance.',/)
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
