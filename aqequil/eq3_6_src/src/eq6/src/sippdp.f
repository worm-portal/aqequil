      subroutine sippdp(actw,aw0plo,aw0prn,aw1plo,aw1prn,dlaplo,
     $ dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,
     $ dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,dlxprl,
     $ dlxprn,eh,eh0plo,eh0prn,eh1plo,eh1prn,eps100,fo2lg,lprcin,
     $ o20plo,o20prn,o21plo,o21prn,ph,ph0plo,ph0prn,ph1plo,ph1prn,
     $ prcinf,qredox,tiplol,tiplot,tiprnl,tiprnt,tistsv,xidump,
     $ xiplol,xiplot,xiprnl,xiprnt,xistsv)
c
c     This subroutine advances the print, plot, and dump points.
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
      logical qredox
c
      real*8 actw,aw0plo,aw0prn,aw1plo,aw1prn,dlaplo,dlaprn,dleplo,
     $ dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,
     $ dlxdmp,dlxmx0,dlxpll,dlxplo,dlxprl,dlxprn,eh,eh0plo,eh0prn,
     $ eh1plo,eh1prn,eps100,fo2lg,lprcin,o20plo,o20prn,o21plo,o21prn,
     $ ph,ph0plo,ph0prn,ph1plo,ph1prn,prcinf,tiplol,tiplot,tiprnl,
     $ tiprnt,tistsv,xidump,xiplol,xiplot,xiprnl,xiprnt,xistsv
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n
c
      real*8 lloc,tsmall,t1day,xloc,yloc
c
cXX   real*8 t1sec,t1min,t1hr,t1year
c
      real*8 texp,tlg
c
c-----------------------------------------------------------------------
c
      xiprnt = prcinf
      xiprnl = prcinf
      xiplot = prcinf
      xiplol = prcinf
c
      xidump = prcinf
c
      tiprnt = prcinf
      tiprnl = prcinf
      tiplot = prcinf
      tiplol = prcinf
c
      ph0prn = -prcinf
      eh0prn = -prcinf
      o20prn = -prcinf
      aw0prn = -prcinf
      ph0plo = -prcinf
      eh0plo = -prcinf
      o20plo = -prcinf
      aw0plo = -prcinf
c
      ph1prn = prcinf
      eh1prn = prcinf
      o21prn = prcinf
      aw1prn = prcinf
      ph1plo = prcinf
      eh1plo = prcinf
      o21plo = prcinf
      aw1plo = prcinf
c
cXX   t1sec= 1.0
cXX   t1min= 60.
cXX   t1hr = 3600.
      t1day = 86400.
cXX   t1year = t1day*365.25
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the initial print point for the specified Xi interval.
c
      if (dlxprn .lt. prcinf) then
        n = int( xistsv/dlxprn ) + 1
        xiprnt = n*dlxprn
        if ((xiprnt + eps100) .le. xistsv) xiprnt = xiprnt + dlxprn
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the initial print point for the specified log Xi
c     interval.
c
      if (dlxprl .lt. lprcin) then
        if (xistsv .lt. dlxmx0) then
          xloc = tlg(dlxmx0)
          lloc = int(xloc)
          n = int( ( xloc - lloc )/dlxprl )
          yloc = lloc + n*dlxprl
          if ((yloc + eps100) .le. xloc) yloc = yloc + dlxprl
        else
          xloc = tlg(xistsv)
          lloc = int(xloc)
          n = int( ( xloc - lloc )/dlxprl ) + 1
          yloc = lloc + n*dlxprl
          if ((yloc + eps100) .le. xloc) yloc = yloc + dlxprl
        endif
        xiprnl = texp(yloc)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the initial plot point for the specified Xi interval.
c
      if (dlxplo .lt. prcinf) then
        n = int( xistsv/dlxplo ) + 1
        xiplot = n*dlxplo
        if ((xiplot + eps100) .le. xistsv) xiplot = xiplot + dlxplo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the initial plot point for the specified log Xi
c     interval.
c
      if (dlxpll .lt. lprcin) then
        if (xistsv .lt. dlxmx0) then
          xloc = tlg(dlxmx0)
          lloc = int(xloc)
          n = int( ( xloc - lloc )/dlxpll )
          yloc = lloc + n*dlxpll
          if ((yloc + eps100) .le. xloc) yloc = yloc + dlxpll
        else
          xloc = tlg(xistsv)
          lloc = int(xloc)
          n = int( ( xloc - lloc )/dlxpll ) + 1
          yloc = lloc + n*dlxpll
          if ((yloc + eps100) .le. xloc) yloc = yloc + dlxpll
        endif
        xiplol = texp(yloc)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the initial dump point for the specified Xi
c     interval.
c
      if (dlxdmp. lt. prcinf) then
        n = int( xistsv/dlxdmp ) + 1
        xidump = n*dlxdmp
        if ((xidump + eps100) .le. xistsv) xidump = xidump + dlxdmp
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the initial print point for the specified time interval.
c
      if (dltprn .lt. prcinf) then
        n = int( tistsv/dltprn ) + 1
        tiprnt = n*dltprn
        if ((tiprnt + eps100) .le. tistsv) tiprnt = tiprnt + dltprn
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the initial print point for the specified log time
c     interval.
c
      if (dltprl .lt. lprcin) then
        tsmall = t1day
        if (tistsv .lt. tsmall) then
          xloc = tlg(tsmall)
          lloc = int(xloc)
          n = int( ( xloc - lloc )/dltprl )
          yloc = lloc + n*dltprl
          if ((yloc + eps100) .le. xloc) yloc = yloc + dltprl
        else
          xloc = tlg(tistsv)
          lloc = int(xloc)
          n = int( ( xloc - lloc )/dltprl ) + 1
          yloc = lloc + n*dltprl
          if ((yloc + eps100) .le. xloc) yloc = yloc + dltprl
        endif
        tiprnl = texp(yloc)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the initial plot point for the specified time interval.
c
      if (dltplo .lt. prcinf) then
        n = int( tistsv/dltplo ) + 1
        tiplot = n*dltplo
        if ((tiplot + eps100) .le. tistsv) tiplot = tiplot + dltplo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the initial plot point for the specified log time
c     interval.
c
      if (dltpll .lt. lprcin) then
        tsmall = t1day
        if (tistsv .lt. tsmall) then
          xloc = tlg(tsmall)
          lloc = int(xloc)
          n = int( ( xloc - lloc )/dltpll )
          yloc = lloc + n*dltpll
          if ((yloc + eps100) .le. xloc) yloc = yloc + dltpll
        else
          xloc = tlg(tistsv)
          lloc = int(xloc)
          n = int( ( xloc - lloc )/dltpll ) + 1
          yloc = lloc + n*dltpll
          if ((yloc + eps100) .le. xloc) yloc = yloc + dltpll
        endif
        tiplol = texp(yloc)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the initial print points for the specified pH interval.
c     Note there are two such points, one lower and one higher than the
c     current value.
c
      if (dlhprn .lt. prcinf) then
        n = int( ph/dlhprn )
        ph0prn = n*dlhprn
        if ((ph0prn - eps100) .ge. ph) ph0prn = ph0prn - dlhprn
        ph1prn = (n + 1)*dlhprn
        if ((ph1prn + eps100) .le. ph) ph1prn = ph1prn + dlhprn
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the initial print points for the specified aw interval.
c     Note there are two such points, one lower and one higher than the
c     current value.
c
      if (dlaprn .lt. prcinf) then
        n = int( actw/dlaprn )
        aw0prn = n*dlaprn
        if ((aw0prn - eps100) .ge. actw) aw0prn = aw0prn - dlaprn
        aw1prn = (n + 1)*dlaprn
        if ((aw1prn + eps100) .le. actw) aw1prn = aw1prn + dlaprn
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
