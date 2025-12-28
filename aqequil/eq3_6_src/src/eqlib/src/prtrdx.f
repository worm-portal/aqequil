      subroutine prtrdx(ah,ahrc,cdrsd,eh,ehrc,fo2lg,fo2lrc,jflgi,
     $ jsflag,narn1,nbasp,nbaspd,nbt,nbtmax,ndrsd,ndrsmx,ndrsrd,
     $ nelect,nhydr,no2gaq,noutpt,nstmax,pe,perc,uspec)
c
c     This subroutine prints a table of the Eh, pe-, log fO2, and Ah for
c     the default redox constraint and each aqueous redox couple which
c     is not required to satisfy the default redox constraint.
c
c     Most of the data printed by this subroutine are calculated by
c     EQLIB/cdardx.f.
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
c       ah     = Default Ah, kcal
c       ahrc   = array of couple-specific Ah values, kcal
c       cdrsd  = coefficients of reactions in the 'd' set
c       eh     = Default Eh, volts
c       ehrc   = array of couple-specific Eh values, volts
c       fo2lg  = Default log fO2
c       fo2lrc = array of couple-specific log fO2 values
c       jsflag = array of species status flags
c       nbt    = the number of species in the basis set
c       ndrsd  = array of the indices of the species corresponding to
c                  the reaction coefficients in the cdrsd array
c       ndrsrd = array giving the range in the cdrsd and ndrsd arrays
c                  corresponding to the reaction for a given species
c       nelect = index of the fictive species aqueous e-
c       no2gaq = index of the fictive species aqueous O2(g)
c       pe     = Default pe-
c       perc   = array of couple-specific pe- values
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
      integer jflgi(nbtmax),jsflag(nstmax),nbasp(nbtmax),
     $ nbaspd(nbtmax),ndrsd(ndrsmx),ndrsrd(2,nstmax)
      integer narn1,nbt,nelect,nhydr,no2gaq
c
      character*48 uspec(nstmax)
c
      real*8 ahrc(nbtmax),cdrsd(ndrsmx),ehrc(nbtmax),
     $ fo2lrc(nbtmax),perc(nbtmax)
      real*8 ah,eh,fo2lg,pe
c
      real*8 coefdr
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,j3,n,nb,nr1,nr2,ns1,ns2,nse,nt
c
      integer ilnobl
c
      character*24 ux24
c
      real*8 cx1
c
c-----------------------------------------------------------------------
c
      write (noutpt,1000)
 1000 format(//16x,'--- Aqueous Redox Reactions ---',
     $ //3x,'Couple',27x,'Eh, volts',6x,'pe-',6x,'log fO2',3x,
     $ 'Ah, kcal',/)
      write (noutpt,1010) eh,pe,fo2lg,ah
 1010 format(1x,'DEFAULT',t37,f7.3,3x,1pe11.4,2x,0pf8.3,2x,f8.3)
c
      do nb = 1,nbt
        if (jflgi(nb) .eq. 30) go to 110
        if (jflgi(nb) .eq. 27) go to 110
        ns1 = nbaspd(nb)
        ns2 = nbasp(nb)
        if (jsflag(ns1) .ne. 0) go to 110
        if (jsflag(ns2) .ne. 0) go to 110
        nt = ndrsrd(2,ns1) - ndrsrd(1,ns1) + 1
        if (nt .lt. 2) go to 110
c
c       Test the reaction read from the data file to see if it contains
c       an entry for O2(g).
c
c       Calling sequence substitutions:
c         cdrsd for cdrs
c         ndrsd for ndrs
c         ndrsrd for ndrsr
c         no2gaq for nse
c         ns1 for ns
c
        cx1 = coefdr(cdrsd,ndrsd,ndrsmx,ndrsrd,no2gaq,ns1,nstmax)
        if (cx1 .ne. 0.) then
          ux24 = uspec(narn1)
          nr1 = ndrsrd(1,ns1)
          nr2 = ndrsrd(2,ns1)
          do n = nr1,nr2
            nse = ndrsd(n)
            if (nse .eq. 0) go to 110
            if (nse.ne.ns1 .and. nse.ne.no2gaq .and. nse.ne.nelect
     $      .and. nse.ne.narn1 .and. nse.ne.nhydr) ux24 = uspec(nse)
          enddo
c
          j2 = ilnobl(uspec(ns1)(1:24))
          j3 = ilnobl(ux24)
c
  100     if ((j2 + j3) .gt. 32) then
            if (j2 .gt. j3) then
              j2 = j2 - 1
            else
              j3 = j3 - 1
            endif
            go to 100
          endif
c
          write (noutpt,1020) uspec(ns1)(1:j2),ux24(1:j3),ehrc(nb),
     $    perc(nb),fo2lrc(nb),ahrc(nb)
 1020     format(1x,a,'/',a,t37,f7.3,3x,1pe11.4,2x,0pf8.3,2x,f8.3)
        endif
  110   continue
      enddo
      write (noutpt,1030)
 1030 format(/4x,'Couples required to satisfy the default redox',
     $ ' constraint are not listed.',/)
c
      end
