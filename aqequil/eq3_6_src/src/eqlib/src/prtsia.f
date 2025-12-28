      subroutine prtsia(affsd,jflagd,jflgi,jsflag,narn1,narn2,nbasp,
     $ nbaspd,nbt,nbtmax,ndrsd,ndrsmx,ndrsrd,nhydr,noutpt,nrdxsp,
     $ nstmax,sidrsp,uspec)
c
c     This subroutine prints tables of saturation indices and affinities
c     for reactions among aqueous species which are not constrained
c     to be at equilibrium.
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
c       affsd  = array of affinities for reactions corresponding to
c                  the various species, as these reactions were read
c                  from the data file
c       jflagd = array of species control flags, used with the 'd'
c                  set of reactions
c       jsflag = array of status flags for species
c       nbasp  = array of the indices of the species in the active
c                  basis set
c       nbaspd = array of the indices of the species in the data file
c                  basis set
c       nbt    = the number of species in the basis set
c       ndrsd  = array of the indices of the species corresponding to
c                  the reaction coefficients in the cdrsd array
c       ndrsrd = array giving the range in the cdrsd and ndrsd arrays
c                  corresponding to the reaction for a given species
c       nhydr  = index of the species aqueous H+
c       nrdxsp = index of the redox basis species
c       sidrsp = array of saturation indices for reactions
c                  corresponding to the various species, as these
c                  reactions were read from the data file
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
c     Calling sequence variable declarations.
c
      integer nbtmax,ndrsmx,nstmax
c
      integer noutpt
c
      integer jflagd(nstmax),jflgi(nbtmax),jsflag(nstmax),
     $ nbasp(nbtmax),nbaspd(nbtmax),ndrsd(ndrsmx),ndrsrd(2,nstmax)
      integer narn1,narn2,nbt,nhydr,nrdxsp
c
      character(len=48) uspec(nstmax)
c
      real(8) affsd(nstmax),sidrsp(nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,j3,kount,n,nb,nbb,nr1,nr2,nse,ns1,ns2,nt
c
      integer ilnobl
c
      character(len=24) ux24
c
      character(len=16) ux16a,ux16b
c
c-----------------------------------------------------------------------
c
      write (noutpt,1000)
 1000 format(//6x,'--- Saturation States of Aqueous Reactions Not',
     $ ' Fixed at Equilibrium ---')
      write (noutpt,1010)
 1010   format(/3x,'Reaction',27x,'Log Q/K',4x,'Affinity, kcal',/)
c
      kount = 0
      do nb = 1,nbt
        ns1 = nbaspd(nb)
        if (ns1.lt.narn1 .or. ns1.gt.narn2) go to 250
        if (jflagd(ns1) .eq. 27) go to 250
        if (jflagd(ns1) .eq. 30) go to 250
        if (jsflag(ns1) .gt. 0) go to 250
c
        nr1 = ndrsrd(1,ns1)
        nr2 = ndrsrd(2,ns1)
        nt = nr2 - nr1 + 1
        if (nt .lt. 2) go to 250
c
c       Find the appropriate matching species to complete the
c       description of the reaction.
c
        ns2 = 0
        do n = nr1,nr2
          nse = ndrsd(n)
          if (nse.ne.ns1 .and. nse.ne.nrdxsp
     $      .and. nse.ne.narn1 .and. nse.ne.nhydr) then
            ns2 = nse
            go to 130
          endif
        enddo
  130   continue
c
        if (ns2 .eq. 0) then
          do n = nr1,nr2
            nse = ndrsd(n)
            if (nse.ne.ns1 .and. nse.eq.narn1) then
              ns2 = narn1
              go to 140
            endif
          enddo
  140     continue
        endif
c
        if (ns2 .eq. 0) then
          do n = nr1,nr2
            nse = ndrsd(n)
            if (nse.ne.ns1 .and. nse.eq.nhydr) then
              ns2 = nhydr
              go to 150
            endif
          enddo
  150     continue
        endif
c
        if (ns2 .eq. 0) then
          do n = nr1,nr2
            nse = ndrsd(n)
            if (nse.ne.ns1 .and. nse.eq.nrdxsp) then
              ns2 = nrdxsp
              go to 160
            endif
          enddo
  160     continue
        endif
c
c       At this point, ns2 should not be zero.
c
        if (ns2 .eq. 0) then
          go to 250
        endif
c
        if (jsflag(ns2) .gt. 0) go to 250
        ux24 = uspec(ns2)
c
        kount = kount + 1
c
        j2 = ilnobl(uspec(ns1)(1:24))
        j3 = ilnobl(ux24)
c
  240   if ((j2 + j3) .gt. 32) then
          if (j2 .gt. j3) then
            j2 = j2 - 1
          else
            j3 = j3 - 1
          endif
          go to 240
        endif
c
        if (sidrsp(ns1).gt.-9999999.
     $   .and. sidrsp(ns1).lt.9999999.) then
          write (ux16a,'(f10.5)') sidrsp(ns1)
          write (ux16b,'(f10.5)') affsd(ns1)
        else
          ux16a = '    N/A   '
          ux16b = '    N/A   '
        endif
c
        write (noutpt,1020) uspec(ns1)(1:j2),ux24(1:j3),ux16a(1:10),
     $  ux16b(1:10)
 1020   format(1x,a,'/',a,t37,a,3x,a)
  250   continue
      enddo
c
      if (kount.le.0) write (noutpt,1030)
 1030 format(1x,'None')
c
      write (noutpt,1040)
 1040 format(1x)
c
      end
