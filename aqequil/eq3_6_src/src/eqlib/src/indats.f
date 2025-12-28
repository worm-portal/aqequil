      subroutine indats(adhfsa,advfsa,axhfsa,axlksa,axvfsa,cdrsa,cdrsv,
     $ cessa,cessv,ipch,ipch_asv,ipcv,ipcv_asv,mwtspa,nad1,narxt,
     $ narx_asv,nata,nata_asv,nbta,nbta_asv,nbta1_asv,nbtafd,ncmpra,
     $ ncta,ncta_asv,ndrsa,ndrsa_asv,ndrsn,ndrsra,nerr,nessa,nessa_asv,
     $ nessn,nessra,ngta,ngta_asv,nlta,nlta_asv,nmta,nmta_asv,noutpt,
     $ np,npta_asv,ns,nsta_asv,ntprt,ntpr_asv,nttyo,uaqsln,ubasp,udrsv,
     $ uelema,uendit,uessv,uphasa,uphasv,uptgas,uptliq,uptsld,uptypa,
     $ usblkf,uspeca,vosp0a,zchara)
c
c     This subroutine reads a superblock of species blocks from the
c     supporting data file "data1".
c
c     This subroutine is called by:
c
c       EQLIB/indata.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       nad1   = unit number of the data file
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
      integer ipch_asv,ipcv_asv,narx_asv,nata_asv,nbta_asv,nbta1_asv,
     $ ncta_asv,ndrsa_asv,nessa_asv,ngta_asv,nlta_asv,nmta_asv,npta_asv,
     $ nsta_asv,ntpr_asv
c
      integer noutpt,nttyo
c
      integer narxt(ntpr_asv),ncmpra(2,npta_asv),ndrsa(ndrsa_asv),
     $ ndrsra(2,nsta_asv),nessa(nessa_asv),nessra(2,nsta_asv)
c
      integer ipch,ipcv,nad1,nata,nbta,nbtafd,ncta,ndrsn,nerr,
     $ nessn,ngta,nlta,nmta,np,ns,ntprt
c
      character*48 ubasp(nbta_asv),uspeca(nsta_asv)
      character*24 udrsv(nbta1_asv),uphasa(npta_asv),uptypa(npta_asv)
      character*24 uphasv,uaqsln,uptsld,uptliq,uptgas,usblkf
      character*8 uelema(ncta_asv),uessv(ncta_asv)
      character*8 uendit
c
c
      real*8 adhfsa(narx_asv,ntpr_asv,ipch_asv,nsta_asv),
     $ advfsa(narx_asv,ntpr_asv,ipcv_asv,nsta_asv),
     $ axhfsa(narx_asv,ntpr_asv,nsta_asv),
     $ axlksa(narx_asv,ntpr_asv,nsta_asv),
     $ axvfsa(narx_asv,ntpr_asv,nsta_asv),cdrsa(ndrsa_asv),
     $ cdrsv(nbta1_asv),cessa(nessa_asv),cessv(ncta_asv),
     $ mwtspa(nsta_asv),vosp0a(nsta_asv),zchara(nsta_asv)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ipc,j2,j3,j4,j5,j6,j7,n,nb,nc,nctsv,ndrstv,nf,ntpr
c
      integer ilnobl
c
      character*24 uelect,uhx,uspecv,ux24
c
c-----------------------------------------------------------------------
c
      data uelect /'e-                      '/
c
c-----------------------------------------------------------------------
c
c     The label below defines a loop for reading the blocks in the
c     superblock.
c
  100 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read a species block.
c
c     uspecv = species name
c     nctsv  = number of chemical elements comprising the species
c     ndrstv = number of reaction coefficients in the reaction
c              associated with the species
c
      read (nad1) uspecv,nctsv,ndrstv
      j2 = ilnobl(uspecv)
      j5 = ilnobl(usblkf)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for the end of the current superblock.
c
      if (uspecv(1:8) .eq. uendit(1:8)) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     If the species is a pure mineral or pure liquid, set the phase
c     name to match the species name.
c
      if (usblkf(1:24).eq.uptsld(1:24) .or.
     $ usblkf(1:24).eq.uptliq(1:24)) uphasv(1:24) = uspecv(1:24)
      j3 = ilnobl(uphasv)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Have another species.
c
      ns = ns + 1
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Load the species name array. Note that the last 24 characters
c     contain the name of the associated phase.
c
      uspeca(ns)(1:24) = uspecv(1:24)
      uspeca(ns)(25:48) = uphasv(1:24)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     mwtspa = molecular weight
c     zchara = electrical charge
c     vosp0a = molar volume (not read).
c
      read (nad1) mwtspa(ns),zchara(ns)
c
      vosp0a(ns) = 0.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     cessv  = chemical element stoichiometric coefficient
c     uessv  = name of corresponding chemical element.
c
      if (nctsv .gt. 0) read (nad1) (cessv(n),uessv(n), n = 1,nctsv)
c
c     Decode the elemental composition.
c
      nessra(1,ns) = nessn + 1
      if (nctsv .gt. 0) then
        do n = 1,nctsv
          do nc = 1,ncta
            if (uessv(n)(1:8) .eq. uelema(nc)(1:8)) then
              nessn = nessn + 1
              cessa(nessn) = cessv(n)
              nessa(nessn) = nc
              go to 120
            endif
          enddo
          j4 = ilnobl(uessv(n))
          write (noutpt,1010) uspecv(1:j2),uphasv(1:j3),usblkf(1:j5),
     $    uessv(n)(1:j4)
          write (nttyo,1010) uspecv(1:j2),uphasv(1:j3),usblkf(1:j5),
     $    uessv(n)(1:j4)
 1010     format(/' * Error - (EQLIB/indats) The species',
     $    /7x,a,' (',a,')',
     $    /7x,'appearing on the data file in the ',a,' superblock',
     $    /7x,'is listed as being composed of an unrecognized',
     $    /7x,'chemical element ',a,'.')
          nerr = nerr + 1
  120     continue
        enddo
        nessra(2,ns) = nessn
      else
c
c       Have a species comprised of no chemical elements. This is
c       permitted only for the following:
c         O2(g) in a phase other than "gas"
c         e-    in any phase
c
        if (uspecv(1:5) .eq. 'O2(g)') then
          if (usblkf(1:24) .eq. uptgas(1:24)) then
            write (noutpt,1020) uspecv(1:j2),uphasv(1:j3),usblkf(1:j5)
            write (nttyo,1020) uspecv(1:j2),uphasv(1:j3),usblkf(1:j5)
 1020       format(/' * Error - (EQLIB/indats) The species',
     $      /7x,a,' (',a,')',
     $      /7x,'appearing on the data file in the ',a,' superblock',
     $      /7x,'is comprised of no chemical elements.')
            nerr = nerr + 1
          endif
        elseif (uspecv(1:2) .eq. uelect(1:2)) then
          continue
        else
          write (noutpt,1020) uspecv(1:j2),uphasv(1:j3),usblkf(1:j5)
          write (nttyo,1020) uspecv(1:j2),uphasv(1:j3),usblkf(1:j5)
          nerr = nerr + 1
        endif
c
c       A species with no elements is treated as being comprised of
c       one "null" element.
c
  125   nessn = nessn + 1
        cessa(nessn) = 0.
        nessa(nessn) = 0
        nessra(2,ns) = nessn
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the associated reaction data, if any.
c
c     Read the reaction, if any.
c
      if (ndrstv .gt. 0) then
        read (nad1) (cdrsv(n),udrsv(n), n = 1,ndrstv)
      endif
c
c     Decode the associated reaction into a temporary format. The
c     species associated with the reaction is moved if necessary so
c     that it is the first species to appear in the reaction. Its
c     species index is stored in the ndrsa array. The other species
c     are presumed to be basis species. If one of these species has
c     not been already loaded into the ubasp array, it is loaded at
c     this time. The basis number of the species is stored in the
c     ndrsa array. It will be converted to the corresponding species
c     number after all species are loaded. This scheme permits the
c     loading of a reaction which is written in terms of a basis
c     species which itself has not yet been loaded.
c
      ndrsra(1,ns) = ndrsn + 1
      if (ndrstv .le. 0) then
c
c       A species with no reaction listed on the data file is treated
c       as a basis species. An associated reaction with one "null"
c       species is created. The actual reaction is formally treated
c       as one in which one mole of the associated species is destroyed
c       and one mole of the same species is created. This is termed
c       an "identity reaction." However, the reaction stored is not
c       the identity reaction but the one with the "null" species as
c       the only species.
c
        nbta = nbta + 1
        ubasp(nbta) = uspeca(ns)
        ndrsn = ndrsn + 1
        cdrsa(ndrsn) = 0.
        ndrsa(ndrsn) = 0
        ndrsra(2,ns) = ndrsn
        do ntpr = 1,ntpr_asv
          do n = 1,narx_asv
            axlksa(n,ntpr,ns) = 0.
          enddo
        enddo
      elseif (ndrstv .eq. 1) then
        write (noutpt,1030) uspecv(1:j2),uphasv(1:j3),usblkf(1:j5)
        write (nttyo,1030) uspecv(1:j2),uphasv(1:j3),usblkf(1:j5)
 1030   format(/' * Error - (EQLIB/indats) The species'
     $  /7x,a,' (',a,')',
     $  /7x,'appearing on the data file in the ',a,' superblock',
     $  /7x,'has a reaction with only one species in it.')
        nerr = nerr + 1
      else
c
c       Make sure that the species associated with the reaction is
c       the first species in the reaction.
c
        do n = 1,ndrstv
          if (udrsv(n)(1:24) .eq. uspecv(1:24)) then
            ndrsn = ndrsn + 1
            cdrsa(ndrsn) = cdrsv(n)
            ndrsa(ndrsn) = ns
            go to 160
          endif
        enddo
c
        write (noutpt,1040) uspecv(1:j2),uphasv(1:j3),usblkf(1:j5)
        write (nttyo,1040) uspecv(1:j2),uphasv(1:j3),usblkf(1:j5)
 1040   format(/' * Error - (EQLIB/indats) The species',
     $  /7x,a,' (',a,')',
     $  /7x,'appearing on the data file in the ',a,' superblock',
     $  /7x,'does not appear in its associated reaction.')
        j4 = ilnobl(udrsv(1))
        write (noutpt,1050) udrsv(1)(1:j4),uphasv(1:j3),uspecv(1:j2),
     $  uphasv(1:j3)
        write (nttyo,1050) udrsv(1)(1:j4),uphasv(1:j3),uspecv(1:j2),
     $  uphasv(1:j3)
 1050   format(/9x,'Is ',a,' (',a,') a doppelganger of',
     $  /12x,a,' (',a,')?')
        udrsv(1) = uspecv
        nerr = nerr + 1
c
  160   nf = 0
        do n = 1,ndrstv
          if (udrsv(n)(1:24).eq.uspecv(1:24) .and. nf.eq.0) then
            nf = 1
            if (ns .le. nbtafd) then
c
c             The species is in the range of aqueous species that are
c             intended to be basis species. Add it to the basis set.
c
              nbta = nbta + 1
              ubasp(nbta) = udrsv(n)
            endif
          else
            ndrsn = ndrsn + 1
            cdrsa(ndrsn) = cdrsv(n)
            do nb = 1,nbta
              if (ubasp(nb)(1:24) .eq. udrsv(n)(1:24)) then
                ndrsa(ndrsn) = nb
                go to 180
              endif
            enddo
c
c           The reaction is written in terms of another species that
c           is not yet in the basis set. Add it to that set.
c
            nbta = nbta + 1
            ubasp(nbta) = udrsv(n)
            ndrsa(ndrsn) = nbta
          endif
  180     continue
        enddo
        ndrsra(2,ns) = ndrsn
c
c       Read the interpolating polynomial coefficients for the log K
c       of the reaction.
c
        read (nad1) ux24
        j6 = ilnobl(ux24)
        uhx = 'Log K'
        j7 = ilnobl(uhx)
        if (uhx(1:j7) .ne. ux24(1:j6)) then
          write (noutpt,1052) ux24(1:j6),uhx(1:j7)
          write (nttyo,1052) ux24(1:j6),uhx(1:j7)
 1052     format(/' * Error - (EQLIB/indats) Have read the tag',
     $    ' string "',a,'"',/7x,'where "',a,'" was expected.')
        endif
        do ntpr = 1,ntprt
          read (nad1) (axlksa(n,ntpr,ns), n = 1,narxt(ntpr))
        enddo
      endif
c
      if (ipch .ge. 0) then
c
c       Read the interpolating polynomial coefficients for the enthalpy
c       function and its pressure derivatives. For a basis species, the
c       enthalpy function is the standard partial molar enthalpy of
c       formation; for all other species, the enthalpy function is
c       the standard partial molar enthalpy of reaction.
c
        read (nad1) ux24
        j6 = ilnobl(ux24)
        uhx = 'Enthalpy'
        j7 = ilnobl(uhx)
        if (uhx(1:j7) .ne. ux24(1:j6)) then
          write (noutpt,1052) ux24(1:j6),uhx(1:j7)
          write (nttyo,1052) ux24(1:j6),uhx(1:j7)
        endif
        do ntpr = 1,ntprt
          read (nad1) (axhfsa(n,ntpr,ns), n = 1,narxt(ntpr))
        enddo
        do ipc = 1,ipch
          read (nad1) ux24
          j6 = ilnobl(ux24)
          if (ipc .eq. 1) then
            uhx = 'dEnthalpy/dP'
          else
            uhx = 'd( )Enthalpy/dP( )'
            write (uhx(3:3),'(i1)') ipc,ipc
            write (uhx(17:17),'(i1)') ipc,ipc
          endif
          j7 = ilnobl(uhx)
          if (uhx(1:j7) .ne. ux24(1:j6)) then
            write (noutpt,1052) ux24(1:j6),uhx(1:j7)
            write (nttyo,1052) ux24(1:j6),uhx(1:j7)
          endif
          do ntpr = 1,ntprt
            read (nad1) (adhfsa(n,ntpr,ipc,ns), n = 1,narxt(ntpr))
          enddo
        enddo
      endif
c
      if (ipcv .ge. 0) then
c
c       Read the interpolating polynomial coefficients for the volume
c       function and its pressure derivatives. For a basis species, the
c       volume function is the standard partial molar enthalpy of
c       formation; for all other species, the volume function is
c       the standard partial molar volume of reaction.
c
        read (nad1) ux24
        j6 = ilnobl(ux24)
        uhx = 'Volume'
        j7 = ilnobl(uhx)
        if (uhx(1:j7) .ne. ux24(1:j6)) then
          write (noutpt,1052) ux24(1:j6),uhx(1:j7)
          write (nttyo,1052) ux24(1:j6),uhx(1:j7)
        endif
        do ntpr = 1,ntprt
          read (nad1) (axvfsa(n,ntpr,ns), n = 1,narxt(ntpr))
        enddo
        do ipc = 1,ipcv
          read (nad1) ux24
          j6 = ilnobl(ux24)
          if (ipc .eq. 1) then
            uhx = 'dVolume/dP'
          else
            uhx = 'd( )Volume/dP( )'
            write (uhx(3:3),'(i1)') ipc,ipc
            write (uhx(15:15),'(i1)') ipc,ipc
          endif
          j7 = ilnobl(uhx)
          if (uhx(1:j7) .ne. ux24(1:j6)) then
            write (noutpt,1052) ux24(1:j6),uhx(1:j7)
            write (nttyo,1052) ux24(1:j6),uhx(1:j7)
          endif
          do ntpr = 1,ntprt
            read (nad1) (advfsa(n,ntpr,ipc,ns), n = 1,narxt(ntpr))
          enddo
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (usblkf(1:24) .eq. uaqsln(1:24)) nata = nata + 1
c
      if (usblkf(1:24) .eq. uptsld(1:24)) then
        nmta = nmta + 1
        np = np + 1
        uphasa(np) = uspecv
        uptypa(np) = uptsld
        ncmpra(1,np) = ns
        ncmpra(2,np) = ns
      endif
c
      if (usblkf(1:24) .eq. uptliq(1:24)) then
        nlta = nlta + 1
        np = np + 1
        uphasa(np) = uspecv
        uptypa(np) = uptliq
        ncmpra(1,np) = ns
        ncmpra(2,np) = ns
      endif
c
      if (usblkf(1:24) .eq. uptgas(1:24)) ngta = ngta + 1
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Loop back to read another species block.
c
      go to 100
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
