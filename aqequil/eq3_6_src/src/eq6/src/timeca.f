      subroutine timeca(deltim,delxi,drir0,iodb,nodbmx,nord,noutpt,
     $ nrd1mx,nttyo,prcinf,qriinf,rirec0,time0,time1)
c
c     This subroutine calculates the interval (deltim) of model time
c     that corresponds to the reaction progress interval (delxi).
c
c     This subroutine is called by:
c
c       EQ6/chksti.f
c       EQ6/eqshel.f
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
      integer nodbmx,nrd1mx
c
      integer noutpt,nttyo
c
      integer iodb(nodbmx)
      integer nord
c
      logical qriinf
c
      real*8 drir0(nrd1mx)
      real*8 deltim,delxi,prcinf,rirec0,time0,time1
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,n
c
      real*8 dxp
c
      real*8 fctrl
c
c-----------------------------------------------------------------------
c
      if (qriinf) then
        write (noutpt,1000)
        write (nttyo,1000)
 1000   format(/' * Note - (EQ6/timeca) Have encountered an infinite',
     $  /7x,'time interval.')
        time1 = prcinf
        deltim = prcinf
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      deltim = rirec0*delxi
      if (nord .gt. 0) then
        dxp = delxi
        do n = 1,nord
          i = n + 1
          dxp = dxp*delxi
          deltim = deltim + ( drir0(n)*dxp )/fctrl(i)
        enddo
      endif
      time1 = time0 + deltim
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (deltim.le.0. .and. delxi.gt.0.) then
        if (iodb(2) .gt. 0) then
          write (noutpt,1010) deltim,nord
          write (nttyo,1010) deltim,nord
 1010     format(/' * Note - (EQ6/timeca) Have encountered a',
     $    /7x,'time increment of ',1pe12.5,' seconds calculated for',
     $    /7x,'order ',i2,'.')
          write (noutpt,1020) delxi,rirec0
 1020     format(7x,'delxi= ',1pe12.5,/7x,'rirec0= ',e12.5)
          if (nord .gt. 0) write (noutpt,1030) (drir0(n), n = 1,nord)
 1030     format(7x,'drir0: ',/12x,1pe12.5,3x,e12.5,3x,e12.5,
     $    /12x,e12.5,3x,e12.5,3x,e12.5)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
