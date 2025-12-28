      subroutine prtpct(conc,csts,ctb,iopr,jcsort,jflag,narn1,narn2,
     $ nbasp,nbaspd,nbt,nbtmax,nelect,nhydr,noprmx,no2gaq,noutpt,
     $ nstmax,nsts,nstsmx,nstsr,uspec)
c
c     This subroutine prints tables giving the percentages of species
c     making up solute mass totals in the aqueous solution. The
c     level of printing is controlled by the print control flag
c     iopr(6):
c
c        0 = Don't print any tables
c        1 = Print tables including 99% of all contributing species
c        2 = Print tables including all contributing species
c
c     This subroutine is called by:
c
c       EQ3NR/scripx.f
c       EQ6/scripz.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       conc   = array of species concentrations
c       csts   = array of stoichiometric mass balance factors
c       ctb    = array of total molalities of data file basis
c                  species
c       iopr   = array of print control options
c       jcsort = array of species indices, in order of increasing
c                  concentration, but with sorting restricted to within
c                  phase ranges
c       narn1  = start of the range of aqueous species
c       narn2  = end of the range of aqueous species
c       nbasp  = array of the species indices in the active basis set
c       nbaspd = array of the species indices in the data file basis
c                  set
c       nbt    = number of species in the basis set
c       nelect = index of the fictive species aqueous e-
c       nhydr  = index of the species aqueous H+
c       no2gaq = index of the fictive species aqueous O2(g)
c       nsts   = array of indices of species corresponding to the
c                  coefficients in the csts array
c       nstsr  = array giving the range in the csts and nsts arrays
c                  corresponding to a given species
c       nwater = index of the species liquid water
c       uspec  = array of species names
c
c     Principal output:
c
c       None
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer nbtmax,noprmx,nstmax,nstsmx
c
      integer noutpt
c
      integer iopr(noprmx),jcsort(nstmax),jflag(nstmax),nbasp(nbtmax),
     $ nbaspd(nbtmax),nsts(nstsmx),nstsr(2,nstmax)
      integer narn1,narn2,nbt,nelect,nhydr,no2gaq
c
      character*(48) uspec(nstmax)
c
      real*8 conc(nstmax),csts(nstsmx),ctb(nbtmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,nb,ns,ns1,ns2,nss,nssi
c
      integer ilnobl
c
      real*8 apct,cx,cxn,cxs,cxt,fx,px,pxs,pxt
      real*8 coefst
c
c-----------------------------------------------------------------------
c
      if (iopr(6) .le. -1) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopr(6) .le. 0) then
        write (noutpt,1000)
 1000   format(//6x,'--- Major Species by Contribution to Aqueous',
     $  ' Mass Balances ---',/)
      else
        write (noutpt,1010)
 1010   format(//6x,'--- Species by Contribution to Aqueous Mass',
     $  ' Balances ---',/)
      endif
c
      do nb = 1,nbt
        ns1 = nbaspd(nb)
        ns2 = nbasp(nb)
        if (ns1.lt.narn1 .or. ns1.gt.narn2) go to 110
        if (conc(ns2).le.0. .or. jflag(ns2).eq.30) go to 110
        if (ns1.eq.no2gaq .or. ns1.eq.nelect) go to 110
        if (ns1.eq.narn1 .or. ns1.eq.nhydr) go to 110
c
        j2 = ilnobl(uspec(ns1)(1:24))
        if (iopr(6) .eq. 0) then
          write (noutpt,1020) uspec(ns1)(1:j2)
 1020     format(/' Species Accounting for 99% or More of Aqueous ',a,
     $    //5x,'Species',19x,'Factor',4x,'Molality',5x,'Per Cent',/)
        else
          write (noutpt,1030) uspec(ns1)(1:j2)
 1030     format(/' Species Accounting for Total Aqueous ',a,
     $    //5x,'Species',19x,'Factor',4x,'Molality',5x,'Per Cent',/)
        endif
c
        cxt = ctb(nb)
        pxt = 100.
        if (cxt .gt. 0.) then
          cxn = pxt/cxt
        else
          cxn = 0.
        endif
        cxs = 0.
        pxs = 0.
c
        do nss = narn1,narn2
          nssi = narn2 - nss + narn1
          ns = jcsort(nssi)
          if (conc(ns) .gt. 0.) then
            fx = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
            if (fx .gt. 0.) then
              cx = fx*conc(ns)
              cxs = cxs + cx
              px = cx*cxn
              pxs = pxs + px
              apct = abs(px)
c
              if (apct .gt. 999.99) then
                write (noutpt,1040) uspec(ns),fx,conc(ns),px
 1040           format(3x,a24,3x,f6.2,3x,1pe11.4,3x,e10.3)
              elseif (apct .lt. 1.00) then
                write (noutpt,1040) uspec(ns),fx,conc(ns),px
              else
                write (noutpt,1050) uspec(ns),fx,conc(ns),px
 1050           format(3x,a24,3x,f6.2,3x,1pe11.4,3x,0pf7.2)
              endif
c
              if (iopr(6).le.0 .and. pxs.ge.99.) go to 105
            endif
          endif
        enddo
c
  105   write (noutpt,1060)
 1060   format(' - - - - - - - - - - - - - - - - - - - - - - - - - - -',
     $  ' - - - - -')
        if (iopr(6) .eq. 0) then
          write (noutpt,1070) cxs,pxs
 1070     format(3x,'Subtotal',28x,1pe11.4,3x,0pf7.2,/)
        else
          write (noutpt,1080) cxt,pxt
 1080     format(3x,'Total',31x,1pe11.4,3x,0pf7.2,/)
        endif
  110   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
