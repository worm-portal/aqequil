      subroutine gpress(iopt,jpress,noptmx,noutpt,nptkmx,nttyo,presg,
     $ presh,press,pressb,time1,ptk,xi1)
c
c     This subroutine computes the pressure (press) as a function of
c     reaction progress (xi1) or time (time1).
c
c     This subroutine is called by:
c
c       EQ6/eq6.f
c       EQ6/tpadv.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       jpress = pressure tracking flag
c       ptk    = pressure tracking coefficients
c       time1  = time variable
c       xi1    = reaction progress variable
c
c     Principal output:
c
c       press  = pressure, bars
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer noptmx,nptkmx
c
      integer noutpt,nttyo
c
      integer iopt(noptmx)
      integer jpress
c
      real*8 ptk(nptkmx)
      real*8 presg,presh,press,pressb,time1,xi1
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      real*8 presmx
c
c-----------------------------------------------------------------------
c
c     Maximum allowed pressure (bars):
c
      data presmx /90000./
c
c-----------------------------------------------------------------------
c
      if (jpress .eq. 0) then
c
c       The pressure corresponds to the data file reference presssure
c       curve.
c
        press = presg
      elseif (jpress .eq. 1) then
c
c       The pressure corresponds to the 1.013-bar/steam-saturation
c       presssure curve.
c
        press = presh
      elseif (jpress .eq. 2) then
c
c       Constant pressure.
c
        press = pressb
      elseif (jpress .eq. 3) then
c
c       Linear tracking in reaction progress.
c
        press = pressb + ptk(1)*xi1
      elseif (jpress .eq. 4) then
c
c       Linear tracking in time.
c
        press = pressb + ptk(1)*time1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Stop if the pressure isn't greater than 0.
c
      if (press .le. 0) then
        write (noutpt,1010) press
        write (nttyo,1010) press
 1010   format(/' * Error - (EQ6/gpress) The pressure must be',
     $  /7x,'greater than zero. The current pressure is ',1pg12.5,
     $  ' bars.')
        stop
      endif
c
c     Stop if press is greater than presmx.
c
      if (press .gt. presmx) then
        write (noutpt,1020) press,presmx,xi1
        write (nttyo,1020) press,presmx,xi1
 1020   format(/' * Error - (EQ6/gpress) The calculated pressure',
     $  /7x,'is ',1pg12.5," C, greater than the code's built-in",
     $  /7x,'maximum value of ',g12.5,' bars. Other limits',
     $  /7x,'associated with the supporting data file may also apply.',
     $  /7x,'The current value of reaction progress is ',
     $  g12.5,' moles.')
        if (iopt(2) .gt. 0) then
          write (noutpt,1030) time1
          write (nttyo,1030) time1
 1030     format(7x,'The current time value is ',g12.5,' seconds.')
        endif
        stop
      endif
c
  999 continue
      end
