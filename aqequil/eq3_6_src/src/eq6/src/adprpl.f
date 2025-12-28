      subroutine adprpl(actw,aw0plo,aw0prn,aw1plo,aw1prn,dlaplo,
     $ dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,
     $ dltplo,dltprl,dltprn,dlxpll,dlxplo,dlxprl,dlxprn,dlxdmp,
     $ eh,eh0plo,eh0prn,eh1plo,eh1prn,fo2lg,fxprpl,o20plo,o20prn,
     $ o21plo,o21prn,ph,ph0plo,ph0prn,ph1plo,ph1prn,qplaw0,qplaw1,
     $ qpleh0,qpleh1,qplolt,qplolx,qplott,qplotx,qplo20,qplo21,
     $ qplph0,qplph1,qpraw0,qpraw1,qpreh0,qpreh1,qprnlx,qprnlt,
     $ qprntx,qprntt,qpro20,qpro21,qprph0,qprph1,qredox,time1,
     $ tiplol,tiplot,tiprnl,tiprnt,xidump,xiplol,xiplot,xiprnl,
     $ xiprnt,xi1)
c
c     This subroutine adjusts print, plot, and PRS transfer points
c     as needed for the next step along the reaction path.
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
      logical qplaw0,qplaw1,qpleh0,qpleh1,qplolt,qplolx,qplott,qplotx,
     $ qplo20,qplo21,qplph0,qplph1,qpraw0,qpraw1,qpreh0,qpreh1,qprnlx,
     $ qprnlt,qprntx,qprntt,qpro20,qpro21,qprph0,qprph1,qredox
c
      real*8 actw,aw0plo,aw0prn,aw1plo,aw1prn,dlaplo,dlaprn,dleplo,
     $ dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,
     $ dlxpll,dlxplo,dlxprl,dlxprn,dlxdmp,eh,eh0plo,eh0prn,eh1plo,
     $ eh1prn,fo2lg,fxprpl,o20plo,o20prn,o21plo,o21prn,ph,ph0plo,
     $ ph0prn,ph1plo,ph1prn,time1,tiplol,tiplot,tiprnl,tiprnt,xidump,
     $ xiplol,xiplot,xiprnl,xiprnt,xi1
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      real*8 tolx,tx,xx
c
      real*8 texp,tlg
c
c-----------------------------------------------------------------------
c
c     Reset print points as necessary.
c
      if (qprntx) then
        tolx = fxprpl*dlxprn
  100   xiprnt = xiprnt + dlxprn
        if ((xi1 + tolx) .ge. xiprnt) go to 100
      endif
c
      if (qprntt) then
        tolx = fxprpl*dlxprn
  110   tiprnt = tiprnt + dltprn
        if ((time1 + tolx) .ge. tiprnt) go to 110
      endif
c
      if (qprnlx) then
        tolx = fxprpl*dlxprl
  120   xiprnl = texp(tlg(xiprnl) + dlxprl)
        xx = texp(tlg(xi1) + tolx)
        if (xx .ge. xiprnl) go to 120
      endif
c
      if (qprnlt) then
        tolx = fxprpl*dltprl
  130   tiprnl = texp(tlg(tiprnl) + dltprl)
        tx = texp(tlg(time1) + tolx)
        if (tx .ge. tiprnl) go to 130
      endif
c
      if (qprph0) then
        tolx = fxprpl*dlhprn
  140   ph0prn = ph0prn - dlhprn
        if ((ph - tolx) .le. ph0prn) go to 140
        ph1prn = ph0prn
  150   ph1prn = ph1prn + dlhprn
        if ((ph + tolx) .ge. ph1prn) go to 150
      endif
c
      if (qprph1) then
        tolx = fxprpl*dlhprn
  160   ph1prn = ph1prn + dlhprn
        if ((ph + tolx) .ge. ph1prn) go to 160
        ph0prn = ph1prn
  170   ph0prn = ph0prn - dlhprn
        if ((ph - tolx) .le. ph0prn) go to 170
      endif
c
      if (qredox) then
        if (qpreh0) then
          tolx = fxprpl*dleprn
  180     eh0prn = eh0prn - dleprn
          if ((eh - tolx) .le. eh0prn) go to 180
          eh1prn = eh0prn
  190     eh1prn = eh1prn + dleprn
          if ((eh + tolx) .ge. eh1prn) go to 190
        endif
c
        if (qpreh1) then
          tolx = fxprpl*dleprn
  200     eh1prn = eh1prn + dleprn
          if ((eh + tolx) .ge. eh1prn) go to 200
          eh0prn = eh1prn
  210     eh0prn = eh0prn - dleprn
          if ((eh - tolx) .le. eh0prn) go to 210
        endif
c
        if (qpro20) then
          tolx = fxprpl*dloprn
  220     o20prn = o20prn - dloprn
          if ((fo2lg - tolx) .le. o20prn) go to 220
          o21prn = o20prn
  230     o21prn = o21prn + dloprn
          if ((fo2lg + tolx) .ge. o21prn) go to 230
        endif
c
        if (qpro21) then
          tolx = fxprpl*dloprn
  240     o21prn = o21prn + dloprn
          if ((fo2lg + tolx) .ge. o21prn) go to 240
          o20prn = o21prn
  250     o20prn = o20prn - dloprn
          if ((fo2lg - tolx) .le. o20prn) go to 250
        endif
c
      endif
c
      if (qpraw0) then
        tolx = fxprpl*dlaprn
  260   aw0prn = aw0prn - dlaprn
        if ((actw - tolx) .le. aw0prn) go to 260
        aw1prn = aw0prn
  270   aw1prn = aw1prn + dlaprn
        if ((actw + tolx) .ge. aw1prn) go to 270
      endif
c
      if (qpraw1) then
        tolx = fxprpl*dlaprn
  280   aw1prn = aw1prn + dlaprn
        if ((actw + tolx) .ge. aw1prn) go to 280
        aw0prn = aw1prn
  290   aw0prn = aw0prn - dlaprn
        if ((actw - tolx) .le. aw0prn) go to 290
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Reset plot points as necessary.
c
      if (qplotx) then
        tolx = fxprpl*dlxplo
  400   xiplot = xiplot + dlxplo
        if ((xi1 + tolx)  .ge. xiplot) go to 400
      endif
c
      if (qplott) then
        tolx = fxprpl*dltplo
  410   tiplot = tiplot + dltplo
        if ((time1 + tolx)  .ge. tiplot) go to 410
      endif
c
      if (qplolx) then
        tolx = fxprpl*dlxpll
  420   xiplol = texp(tlg(xiplol) + dlxpll)
        xx = texp(tlg(xi1) + tolx)
        if (xx  .ge. xiplol) go to 420
      endif
c
      if (qplolt) then
        tolx = fxprpl*dltpll
  430   tiplol = texp(tlg(tiplol) + dltpll)
        tx = texp(tlg(time1) + tolx)
        if (tx .ge. tiplol) go to 430
      endif
c
      if (qplph0) then
        tolx = fxprpl*dlhplo
  440   ph0plo = ph0plo - dlhplo
        if ((ph - tolx) .le. ph0plo) go to 440
        ph1plo = ph0plo
  450   ph1plo = ph1plo + dlhplo
        if ((ph + tolx) .ge. ph1plo) go to 450
      endif
c
      if (qplph1) then
        tolx = fxprpl*dlhplo
  460   ph1plo = ph1plo + dlhplo
        if ((ph + tolx) .ge. ph1plo) go to 460
        ph0plo = ph1plo
  470   ph0plo = ph0plo - dlhplo
        if ((ph - tolx) .le. ph0plo) go to 470
      endif
c
      if (qredox) then
        if (qpleh0) then
          tolx = fxprpl*dleplo
  480     eh0plo = eh0plo - dleplo
          if ((eh - tolx) .le. eh0plo) go to 480
          eh1plo = eh0plo
  490     eh1plo = eh1plo + dleplo
          if ((eh + tolx) .ge. eh1plo) go to 490
        endif
c
        if (qpleh1) then
          tolx = fxprpl*dleplo
  500     eh1plo = eh1plo + dleplo
          if ((eh + tolx) .ge. eh1plo) go to 500
          eh0plo = eh1plo
  510     eh0plo = eh0plo - dleplo
          if ((eh - tolx) .le. eh0plo) go to 510
        endif
c
        if (qplo20) then
          tolx = fxprpl*dloplo
  520     o20plo = o20plo - dloplo
          if ((fo2lg - tolx) .le. o20plo) go to 520
          o21plo = o20plo
  530     o21plo = o21plo + dloplo
          if ((fo2lg + tolx) .ge. o21plo) go to 530
        endif
c
        if (qplo21) then
          tolx = fxprpl*dloplo
  540     o21plo = o21plo + dloplo
          if ((fo2lg + tolx) .ge. o21plo) go to 540
          o20plo = o21plo
  550     o20plo = o20plo - dloplo
          if ((fo2lg - tolx) .le. o20plo) go to 550
        endif
c
      endif
c
      if (qplaw0) then
        tolx = fxprpl*dlaplo
  560   aw0plo = aw0plo - dlaplo
        if ((actw - tolx) .le. aw0plo) go to 560
        aw1plo = aw0plo
  570   aw1plo = aw1plo + dlaplo
        if ((actw + tolx) .ge. aw1plo) go to 570
      endif
c
      if (qplaw1) then
        tolx = fxprpl*dlaplo
  580   aw1plo = aw1plo + dlaplo
        if ((actw + tolx) .ge. aw1plo) go to 580
        aw0plo = aw1plo
  590   aw0plo = aw0plo - dlaplo
        if ((actw - tolx) .le. aw0plo) go to 590
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Reset the PRS transfer point as necessary.
c
      if (xidump. le. xi1) then
        tolx = fxprpl*dlxdmp
  700   xidump = xidump + dlxdmp
        if ((xi1 + tolx)  .ge. xidump) go to 700
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
