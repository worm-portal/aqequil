      subroutine chksti(akmat0,drer0,drir0,deltim,delxi,dlxmin,
     $ fdri0,fdrr0,iodb,jreac,kly,kmax,nodbmx,kord,nord,noutpt,
     $ npts,nrct,nrctmx,nrd1mx,nttyo,prcinf,qriinf,rirec0,smp100,
     $ time0,time1,xi0,xi1)
c
c     This subroutine checks the sign of the time increment.
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
      integer kmax,nodbmx,nrctmx,nrd1mx
c
      integer noutpt,nttyo
c
      integer iodb(nodbmx),jreac(nrctmx)
c
      integer kly,kord,nord,npts,nrct
c
      logical qriinf
c
      real*8 akmat0(nrd1mx,nrd1mx),drer0(nrd1mx,nrctmx),
     $ drir0(nrd1mx),fdri0(nrd1mx),fdrr0(nrd1mx,nrctmx)
c
      real*8 deltim,delxi,dlxmin,rirec0,prcinf,smp100,time0,time1,
     $ xi0,xi1
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      real*8 dx,dxsv
c
c-----------------------------------------------------------------------
c
  100 xi1 = xi0 + delxi
      call timeca(deltim,delxi,drir0,iodb,nodbmx,nord,noutpt,
     $ nrd1mx,nttyo,prcinf,qriinf,rirec0,time0,time1)
      if (deltim .le. smp100) then
        dxsv = delxi
        dx = 0.25*delxi
        delxi = max(dx,dlxmin)
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
          write (noutpt,1000) deltim
          write (nttyo,1000) deltim
 1000     format(/' * Warning - (EQ6/chksti) The calculated time',
     $    /7x,"increment is ',1pe11.4. This increment can't be zero",
     $    /7x,"or negative, and shouldn't be too close to such a",
     $    /7x,'condition. The finite difference representation of',
     $    /7x,'the inverse rate may be losing accuracy.')
        endif
        if (dxsv .gt. dlxmin) then
          if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1010) dxsv,delxi
            write (nttyo,1010) dxsv,delxi
 1010       format(7x,'The step size will be cut from ',1pe11.4,' to ',
     $      1pe11.4,'.')
          endif
        else
          nord = nord - 1
          kord = nord
          npts = nord + 1
          kly = 6
          if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1020) nord
            write (nttyo,1020) nord
 1020       format(7x,'The step size is already at the  minimum value.',
     $      /7x,'Cutting the order to ',i2,'.')
          endif
          call rderiv(akmat0,drer0,drir0,fdri0,fdrr0,jreac,nord,
     $    nrct,nrctmx,nrd1mx)
        endif
        go to 100
      endif
c
      end
