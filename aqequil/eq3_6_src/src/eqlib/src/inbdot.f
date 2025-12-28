      subroutine inbdot(azeroa,insgfa,nad1,narn1a,narn2a,nata,
     $ nata_asv,nerr,noutpt,nsta_asv,nttyo,uspeca)
c
c     This subroutine reads from the data file the block of individual
c     species parameters for the B-dot model of aqueous species
c     activity coefficients. For each solute species, these parameters
c     are its hard core diameter (azeroa) and its neutral species flag
c     (insgfa). The latter is relevant only for electrically neutral
c     species. If it is has a value of 0, the log activity coefficient
c     is set to zero; if it has a value of -1, the log activity
c     coefficient is represented by the Drummond (1981) polynomial.
c     Water (the solvent) is not subject to control by the insgfa flag.
c
c     This subroutine is called by:
c
c       EQLIB/indata.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       noutpt = unit number of the output file
c       nttyo  = unit number of the screen file
c       uspeca = array of names of species read from the data file
c
c     Principal output:
c
c       azeroa = array of hard core diameters read from the data file
c       insgfa = array of flags for treating the activity coefficients
c                  of species that happen to be electrically neutral
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nata_asv,nsta_asv
c
      integer insgfa(nata_asv)
      integer nad1,narn1a,narn2a,nata,nerr,noutpt,nttyo
c
      character*48 uspeca(nsta_asv)
c
      real*8 azeroa(nata_asv)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer igx,j2,na,naz,ncount,ns
c
      integer ilnobl
c
      real*8 azdef,azx
c
      character*24 unam,unam1,unam2
      character*8 uendit
c
c-----------------------------------------------------------------------
c
      data uendit /'endit.  '/
c
c     The following variable is the default value for the hard core
c     diameter of an aqueous solute species. Note that hard core
c     diameters are given in units of 10**-8 cm.
c
      data azdef  /4.0/
c
c-----------------------------------------------------------------------
c
c     Initialize the azeroa and insgfa arrays to zero.
c
      do na = 1,nata_asv
        azeroa(na) = 0.
        insgfa(na) = 0
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the block contents.
c
      naz = 0
c
c     Read the species name and azeroa and insgfa data.
c
  110 read (nad1) unam,azx,igx
c
c     Test for the end of the data block.
c
      if (unam(1:8) .eq. uendit(1:8)) go to 120
c
c     Find the species index, ns.
c
c     Calling sequence substitutions:
c       narn1a for nrn1a
c       narn2a for nrn2a
c
      call srchn(narn1a,narn2a,ns,nsta_asv,unam,uspeca)
c
c     If not found among the loaded species, skip.
c
      if (ns .le. 0) then
        j2 = ilnobl(unam)
        write (noutpt,1020) unam(1:j2)
        write (nttyo,1020) unam(1:j2)
 1020   format(/' * Warning - (EQLIB/inbdot) Have B-dot model',
     $  ' parameters listed',/7x,'for an aqueous species named ',a,
     $  ', but there is',/7x,'no data block for such a species.')
        go to 110
      endif
c
c     Compute the aqueous species index na.
c
      na = ns - narn1a + 1
c
c     Test for duplicate input.
c
      if (azeroa(na) .gt. 0.) then
        write (noutpt,1030) unam(1:j2),azx,igx,azeroa(na),insgfa(na)
        write (nttyo,1030) unam(1:j2),azx,igx,azeroa(na),insgfa(na)
 1030   format(/' * Warning - (EQLIB/inbdot) Have a duplicate entry',
     $  ' on the data file for',/7x,'the B-dot model parameters of',
     $  ' the species ',a,'.',/7x,'The duplicate entry values are:',
     $  //9x,'azero= ',f7.2,', insgf= ',i3,
     $  //7x,'The first entry values are:',
     $  //9x,'azero= ',f7.2,', insgf= ',i3,
     $  //7x,'The first entry values will be used.')
        go to 110
      endif
c
      naz = naz + 1
      if (naz .gt. nata) then
c
        write (noutpt,1040) nata
        write (nttyo,1040) nata
 1040   format(/' * Error - (EQLIB/inbdot) There are more entries',
     $  /7x,'for species with  B-dot model parameters than there are',
     $  /7x,'aqueous species, which number ',i4,'.')
        nerr = nerr + 1
        go to 999
      endif
c
      azeroa(na) = azx
      insgfa(na) = igx
      go to 110
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for solute species with no entries.
c
  120 do na = 1,nata
        if (azeroa(na) .le. 0.) go to 140
      enddo
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Assign the default value if needed.
c
  140 write (noutpt,1050) azdef
 1050 format(/' * Note - (EQLIB/inbdot) The following aqueous species',
     $ ' have been assigned',/7x,'a default hard core diameter of',
     $ ' ',f6.3,' x 10**-8 cm:',/)
      ncount = 0
      do ns = narn1a,narn2a
        na = ns - narn1a + 1
        if (azeroa(na) .le. 0.) then
          azeroa(na) = azdef
          if (ns .ne. narn1a) then
            if (ncount .le. 0) then
              unam1 = uspeca(ns)(1:24)
              ncount = 1
            else
              unam2 = uspeca(ns)(1:24)
              j2 = ilnobl(unam2)
              write (noutpt,1060) unam1,unam2(1:j2)
 1060         format(9x,a24,3x,a)
              ncount = 0
            endif
          endif
        endif
      enddo
c
      if (ncount .eq. 1) then
        j2 = ilnobl(unam1)
        write (noutpt,1070) unam1(1:j2)
 1070   format(9x,a)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
