      subroutine ndiagx(delmax,delvec,eps100,idelmx,iebal,iindx1,
     $ irdxc3,jflag,kcarb,kebal,khydr,kmax,ko2gaq,nbasp,nbtmax,nhydr,
     $ noutpt,nstmax,nttyo,screwd,uspec)
c
c     This subroutine attempts to generate diagnostics if hybrid
c     Newton-Raphson iteration fails.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
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
      integer kmax,nbtmax,nstmax
c
      integer noutpt,nttyo
c
      integer iindx1(kmax),jflag(nstmax),nbasp(nbtmax)
      integer idelmx,iebal,irdxc3,kcarb,kebal,khydr,ko2gaq,nhydr
c
      character(len=48) uspec(nstmax)
c
      real(8) delvec(kmax)
      real(8) delmax,eps100,screwd
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,nb,ns
c
      integer ilnobl
c
      logical qadtst,qdltst
c
      real(8) adx,xx
c
c-----------------------------------------------------------------------
c
c     Check to see if idelmx is zero. If so, diagnostics can not be
c     generated.
c
      if (idelmx .eq. 0) then
        write (noutpt,1000)
        write (nttyo,1000)
 1000   format(/' * Note - (EQ3NR/ndiagx) No crash diagnostics were',
     $  /7x,'generated. Look at the contents of the delvec and beta',
     $  /7x,'arrays in the iteration summary for clues.',/)
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check to see if the ion being adjusted for electrical balance
c     was crashing to zero. This implies that an ion of opposite
c     charge must be used to achieve electrical balance.
c
      qdltst = (delvec(idelmx) + screwd) .le. eps100
      if (idelmx.eq.kebal .and. qdltst) then
        write (noutpt,1010)
        write (nttyo,1010)
 1010   format(/' * Note - (EQ3NR/ndiagx) the ion adjusted for',
     $  /7x,'electrical balance is crashing to zero. Electrical',
     $  /7x,'balancing requires an ion of opposite charge.')
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check to see if the ion associated with alkalinity balance is
c     crashing to zero. This implies that non-carbonate alkalinity is
c     greater than or equal to carbonate alkalinity.
c
      if (idelmx.eq.kcarb .and. qdltst) then
        nb = iindx1(idelmx)
        ns = nbasp(nb)
        if (jflag(ns) .eq. 7) then
          j2 = ilnobl(uspec(ns)(1:24))
          write (noutpt,1020) uspec(ns)(1:j2)
          write (nttyo,1020) uspec(ns)(1:j2)
 1020     format(/' * Note - (EQ3NR/ndiagx) The basis species ',a,' is',
     $    /7x,'crashing to zero because the non-carbonate alkalinity',
     $    /7x,'exceeds the specified total alkalinity. The total',
     $    /7x,"alkalinity can't be used to constrain this species",
     $    /7x,'because of excessive interference from OH-, borate,',
     $    /7x,'phosphate, organic acids, metal hydroxy complexes, or',
     $    /7x,'the like. A total CO2 (IR) or ion chromatographic',
     $    /7x,'analysis is required instead.',/)
          go to 999
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check to see if fO2 is crashing. This could be due to
c     a bad combination of constraining the redox state by a
c     non-fO2 option and divergence to zero of an associated
c     ion that is constrained by electrical balance.
c
      qadtst = (delmax - screwd) .le. eps100
      if (idelmx.eq.ko2gaq .and. qadtst .and. iebal.gt.0) then
        if (irdxc3.eq.-1 .and. nbasp(iebal).eq.nhydr) then
          xx = 0.1*screwd
          adx = abs(delvec(khydr))
          if (adx .ge. xx) then
            write (noutpt,1030)
            write (nttyo,1030)
 1030       format(/' * Note - (EQ3NR/ndiagx) The fO2 is crashing,',
     $      /7x,'probably because a bad electrical balance constraint',
     $      /7x,'on H+ is causing the concentration of that ion',
     $      /7x,'to crash to zero.',/)
            go to 999
          endif
        endif
c
        if (irdxc3 .eq. 1) then
          write (noutpt,1040)
          write (nttyo,1040)
 1040     format(/' * Note - (EQ3NR/ndiagx) The fO2 is crashing,',
     $    /7x,'probably because of a bad constraint on one of the',
     $    /7x,'aqueous species appearing in the redox reaction that',
     $    /7x,'is being used to constrain the fO2.',/)
          go to 999
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
