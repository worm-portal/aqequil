      subroutine gtemp(afcnst,al10,iopt,jtemp,noptmx,noutpt,nttkmx,
     $ nttyo,rconst,rtcnst,tempc,tempcb,tempk,time1,ttk,xi1)
c
c     This subroutine computes the temperature (tempc) as a function of
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
c       jtemp  = temperature tracking flag
c       ttk    = temperature tracking coefficients
c       time1  = time variable
c       xi1    = reaction progress variable
c
c     Principal output:
c
c       tempc  = temperature, C
c       tempk  = temperature, K
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer noptmx,nttkmx
c
      integer iopt(noptmx)
      integer jtemp,noutpt,nttyo
c
      real(8) ttk(nttkmx)
      real(8) afcnst,al10,rconst,rtcnst,tempc,tempcb,tempk,time1,xi1
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      real(8) tmpmin,tmpmax
c
c-----------------------------------------------------------------------
c
c     Allowed range of the temperature (C):
c
      data tmpmin,tmpmax /-273.15,1400./
c
c-----------------------------------------------------------------------
c
      if (jtemp .eq. 0) then
c
c       Constant temperature.
c
        tempc = tempcb
      elseif (jtemp .eq. 1) then
c
c       Linear tracking in Xi.
c
        tempc = tempcb + ttk(1)*xi1
      elseif (jtemp .eq. 2) then
c
c       Linear tracking in time.
c
        tempc = tempcb + ttk(1)*time1
      elseif (jtemp .eq. 3) then
c
c       Fluid mixing tracking.
c
        tempc = ( tempcb*ttk(1) + xi1*ttk(2) )/( xi1 + ttk(1) )
      endif
c
      tempk = tempc + 273.15
      rtcnst = 0.001*rconst*tempk
      afcnst = al10*rtcnst
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Stop if tempc is less than tmpmin or greater than tmpmax.
c
      if (tempc.lt.tmpmin .or. tempc.gt.tmpmax) then
        write (noutpt,1010) tempc,tmpmin,tmpmax,xi1
        write (nttyo,1010) tempc,tmpmin,tmpmax,xi1
 1010   format(/' * Error - (EQ6/gtemp) The calculated temperature',
     $  /7x,'is ',g10.3," C, outside the code's built-in allowed",
     $  /7x,'range of ',g10.3,'-',g10.3,' C. Other limits associated',
     $  /7x,'with the supporting data file may also apply.',
     $  /7x,'The current value of reaction progress is ',
     $  g12.5,' moles.')
        if (iopt(2) .gt. 0) then
          write (noutpt,1020) time1
          write (nttyo,1020) time1
 1020     format(7x,'The current time value is ',g10.3,' seconds.')
        endif
        stop
      endif
c
  999 continue
      end
