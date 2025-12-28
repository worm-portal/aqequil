      subroutine cdardx(actlg,actwlg,ah,ahrc,cdrsd,eh,ehfac,ehrc,
     $ farad,fo2lg,fo2lrc,jsflag,mosp,nbasp,nbaspd,nbt,nbtmax,ndrsd,
     $ ndrsmx,ndrsrd,no2gaq,nstmax,pe,perc,ph,xlke,xlksd)
c
c     This subroutine computes the default Eh, pe-, and Ah from the
c     default log fO2, and also computes the Eh, pe-, log fO2, and Ah
c     for each aqueous redox couple.
c
c     This subroutine is called by:
c
c       EQ3NR/scripx.f
c       EQ6/cdappl.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       actlg  = array of log activities of species
c       cdrsd  = coefficients of reactions in the 'd' set
c       jsflag = array of species status flags
c       nbt    = the number of species in the basis set
c       ndrsd  = array of the indices of the species corresponding to
c                  the reaction coefficients in the cdrsd array
c       ndrsrd = array giving the range in the cdrsd and ndrsd arrays
c                  corresponding to the reaction for a given species
c       no2gaq = index of the fictive species aqueous O2(g)
c       xlke   = the log K for the "Eh" reaction
c
c     Principal output:
c
c       fo2lrc = array of couple-specific log fO2 values
c       ehrc   = array of couple-specific Eh values
c       perc   = array of couple-specific pe values
c       ahrc   = array of couple-specific Ah values
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
      integer jsflag(nstmax),nbasp(nbtmax),nbaspd(nbtmax),
     $ ndrsd(ndrsmx),ndrsrd(2,nstmax)
      integer nbt,no2gaq
c
      real*8 actlg(nstmax),ahrc(nbtmax),cdrsd(ndrsmx),ehrc(nbtmax),
     $ fo2lrc(nbtmax),mosp(nstmax),perc(nbtmax),xlksd(nstmax)
      real*8 actwlg,ah,eh,ehfac,farad,fo2lg,pe,ph,xlke
c
      real*8 coefdr
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n,nb,nr1,nr2,ns1,ns2,nse,nt
c
      real*8 ahfac,cx1,ehx,fo2lgx
c
c-----------------------------------------------------------------------
c
      ahfac = 0.001*farad
c
      do nb = 1,nbt
        fo2lrc(nb) = fo2lg
        ehrc(nb) = eh
        perc(nb) = pe
        ahrc(nb) = ah
        ns1 = nbaspd(nb)
        ns2 = nbasp(nb)
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
        if (cx1 .eq. 0.) go to 110
c
        if (xlksd(ns1) .lt. 9999999.) then
          fo2lgx = xlksd(ns1)
          nr1 = ndrsrd(1,ns1)
          nr2 = ndrsrd(2,ns1)
          do n = nr1,nr2
            nse = ndrsd(n)
            if (nse .le. 0) go to 110
            if (nse .ne. no2gaq) then
              if (mosp(nse) .le. 0.) go to 110
              fo2lgx = fo2lgx - cdrsd(n)*actlg(nse)
            endif
          enddo
          fo2lgx = fo2lgx/cx1
          ehx = (ehfac/4.)*(fo2lgx - 4.*ph - 2.*actwlg - xlke)
c
          fo2lrc(nb) = fo2lgx
          ehrc(nb) = ehx
          perc(nb) = ehx/ehfac
          ahrc(nb) = ahfac*ehx
        endif
  110   continue
      enddo
c
      end
