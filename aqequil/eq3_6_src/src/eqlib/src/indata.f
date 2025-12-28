      subroutine indata(aadh,aadhh,aadhv,aaphi,abdh,abdhh,abdhv,
     $ abdot,abdoth,abdotv,adadhh,adadhv,adbdhh,adbdhv,adbdth,
     $ adbdtv,adhfe,adhfsa,advfe,advfsa,amua,aprehw,apresg,
     $ apxa,aslma,atwta,axhfe,axhfsa,axlke,axlksa,axvfe,axvfsa,
     $ azeroa,bpxa,cco2,cdrsa,cessa,eps100,iapxa_asv,iapxta,
     $ iaqsla,ibpxa_asv,ibpxta,ielam,igas,ikta_asv,insgfa,ipbt_asv,
     $ ipch,ipch_asv,ipcv,ipcv_asv,ixrn1a,ixrn2a,jpdblo,jpfc_asv,
     $ jptffl,jsola,mwtspa,nad1,nalpaa,napa_asv,napta,narn1a,narn2a,
     $ narx_asv,narxt,nata,nata_asv,nbaspa,nbta,nbta_asv,nbtafd,
     $ nbta1_asv,ncmpra,ncta,ncta_asv,ndrsa,ndrsa_asv,ndrsra,nessa,
     $ nessa_asv,nessra,ngrn1a,ngrn2a,ngta,ngta_asv,nlrn1a,nlrn2a,
     $ nlta,nlta_asv,nmrn1a,nmrn2a,nmta,nmta_asv,nmuta,nmuta_asv,
     $ nmuxa,noutpt,npta,npta_asv,nslta,nslta_asv,nslxa,nsta,
     $ nsta_asv,ntid_asv,ntitld,ntpr_asv,ntprt,nttyo,nxrn1a,nxrn2a,
     $ nxta,nxta_asv,qclnsa,palpaa,tdamax,tdamin,tempcu,ubasp,
     $ udakey,udatfi,uelema,uspeca,uphasa,uptypa,utitld,vosp0a,
     $ zchara)
c
c     This subroutine reads the supporting data file DATA1, starting
c     with its title. This subroutine should be called after subroutine
c     indath.f has been called to read the header section (which
c     mostly contains array size allocation data required to allocate
c     sufficient array space to store the data read by the present
c     subroutine.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       noutpt = unit number of the output file
c       nttyo  = unit number of the screen file
c
c     Principal output:
c
c       aadh   = array of polynomial coefficients for computing
c                  the Debye-Huckel A(gamma,10) parameter
c       aadhh  = array of polynomial coefficients for computing
c                  the Debye-Huckel A(H) parameter
c       aadhv  = array of polynomial coefficients for computing
c                  the Debye-Huckel A(V) parameter
c       aaphi  = array of polynomial coefficients for computing
c                  the Debye-Huckel A(phi) parameter
c       abdh   = array of polynomial coefficients for computing
c                  the Debye-Huckel B parameter
c       abdhh  = array of polynomial coefficients for computing
c                  the Debye-Huckel B(H) parameter
c       abdhv  = array of polynomial coefficients for computing
c                  the Debye-Huckel B(V) parameter
c       abdot  = array of polynomial coefficients for computing
c                  the Helgeson (1969) B-dot parameter
c       abdoth = array of polynomial coefficients for computing
c                  the B-dot(H) parameter
c       abdotv = array of polynomial coefficients for computing
c                  the B-dot(V) parameter
c       adadhh = array of polynomial coefficients for computing
c                  the the pressure derivatives of the Debye-Huckel
c                  A(H) parameter
c       adadhv = array of polynomial coefficients for computing
c                  the the pressure derivatives of the Debye-Huckel
c                  A(V) parameter
c       adbdhh = array of polynomial coefficients for computing
c                  the the pressure derivatives of the Debye-Huckel
c                  B(H) parameter
c       adbdhv = array of polynomial coefficients for computing
c                  the the pressure derivatives of the Debye-Huckel
c                  B(V) parameter
c       adbdth = array of polynomial coefficients for computing
c                  the the pressure derivatives of the B-dot(H)
c                  parameter
c       adbdtv = array of polynomial coefficients for computing
c                  the the pressure derivatives of the B-dot(V)
c                  parameter
c       adhfsa = array of polynomial coefficients for computing
c                  the pressure derivatives of the enthalpy of reaction
c                  for reactions read from the data file
c       adhfe  = array of polynomial coefficients for computing
c                  the pressure derivatives of the enthalpy of
c                  reaction for the "Eh" reaction
c       advfe  = array of polynomial coefficients for computing
c                  the pressure derivatives of the volume of
c                  reaction for the "Eh" reaction
c       advfsa = array of polynomial coefficients for computing
c                  the pressure derivatives of the volume of reaction
c                  for reactions read from the data file
c       amua   = coefficients for calculating Pitzer mu parameters
c                  as a function of temperature, read from the
c                  data file
c       aprehw = array of polynomial coefficients for computing
c                  the recommended pressure envelope half-width (bars)
c       apresg = array of polynomial coefficients for computing
c                  the pressure (bars) on the 1-atm:steam saturation
c                  curve
c       apxa   = array of interaction coefficients for computing
c                   activity coefficients in solid solutions
c       aslm   = coefficients for calculating Pitzer S-lambda(n)
c                  parameters as a function of temperature, read
c                  read from the data file
c       atwta  = array of atomic weights of the elements read from
c                  the data file
c       axhfe  = array of polynomial coefficients for computing
c                  the enthalpy of reaction for the "Eh" reaction
c       axhfsa = array of polynomial coefficients for computing
c                  the enthalpy of reaction for reactions read
c                  from the data file
c       axlke  = array of polynomial coefficients for computing
c                  the log K for the "Eh" reaction
c       axlksa = array of polynomial coefficients for computing
c                  the log K functions of reactions read from the
c                  data file
c       axvfe  = array of polynomial coefficients for computing
c                  the volume of reaction for the "Eh" reaction
c       axvfsa = array of polynomial coefficients for computing
c                  the volume of reaction for reactions read
c                  from the data file
c       azeroa = array of hard core diameters of the species read
c                  from the data file
c       bpxa   = array of site-mixing parameters for computing activity
c                  coefficients in solid solutions
c       cco2   = array of coefficients for computing the activity
c                  coefficient of CO2(aq) from the Drummond (1981)
c                  equation
c       cdrsa  = array of reaction coefficients for the reactions
c                  read from the data file
c       cessa  = array of elemental composition coefficients for the
c                  species read from the data file
c       iapxta = array of the numbers of non-zero interaction
c                  coefficients for computing activity coefficients in
c                  solid solutions
c       iaqsla = index of the aqueous solution phase, as read from
c                  the data file
c       ibpxta = array of the numbers of non-zero site-mixing
c                  paramters for computing activity coefficients in
c                  solid solutions
c       ielam  = flag indicating no use or use of E-lambda terms
c       igas   = index of the gas phase
c       insgfa = array of flags for treating the activity coefficients
c                  of species that happen to be electrically neutral,
c                  as read from the data file
c       ixrn1a = index of the first non-aqueous solution phase,
c                  as read from the data file
c       ixrn2a = index of the last non-aqueous solution phase,
c                  as read from the data file
c       jpdblo = Pitzer data block organization flag
c       jptffl = Pitzer parameter temperature function flag
c       jsola  = array of flags for the activity coefficient models
c                  to use for solid solutions, as read from the
c                  data file
c       mwtspa = array of the molecular weights of the species read
c                  from the data file
c       nalpaa = pointer array giving the index of the set of alpha
c                  coeffcients for a given set of S-lambda
c                  coefficients, as read from the data file
c       narn1a = start of the range of aqueous species, as read from
c                  the data file
c       narn1a = end of the range of aqueous species, as read from
c                  the data file
c       nbaspa = array of indices of species in basis set, as read
c                  from the data file
c       ncmpra = array giving the start and the end of the range of the
c                  species belonging to a given phase, as read from the
c                  data file
c       ndrsa  = array of indices of the species corresponding to
c                  the reaction coefficients in the cdrsa array
c       ndrsra = array giving the start and end of the range in the
c                  cdrsa and ndrsa arrays of the species in the reaction
c                  which belongs to a given species
c       nessa  = array of indices of the chemical elements
c                  corresponding to the composition coefficients in
c                  the cessa array
c       nessra = array giving the start and end of the range in the
c                  cessa and nessa arrays of the chemical elements
c                  comprising a given species
c       ngrn1a = index of the first species in the gas range, as read
c                  from the data file
c       ngrn2a = index of the last species in the gas range, as read
c                  from the data file
c       nlrn1a = index of the first species in the pure liquid range,
c                  as read from the data file
c       nlrn2a = index of the last species in the pure liquid range,
c                  as read from the data file
c       nmrn1a = index of the first species in the pure mineral range,
c                  as read from the data file
c       nmrn2a = index of the last species in the pure mineral range,
c                  as read from the data file
c       nmuta  = the number of species triplets for which mu
c                  coefficients are defined, as read from the data file
c       nmuxa  = indices of species in triplets for which mu
c                  coefficients are defined, as read from the data file
c       nslta  = number of species pairs for which S-lambda
c                  coefficients are defined, as read from the data file
c       nslxa  = indices of species in pairs for which S-lambda
c                  coefficients are defined, as read from the data file
c       nxrn1a = index of the first species in the solid solution range,
c                  as read from the data file
c       nxrn2a = index of the last species in the solid solution range,
c                  as read from the data file
c       palpaa = array of sets of alpha coefficients read from the
c                  data file
c       qclnsa = array of flags indicating if a species was created
c                  by cloning
c       tdamax = the nominal upper temperature limit of the data file
c       tdamin = the nominal lower temperature limit of the data file
c       ubasp  = array of names of the species in the basis set
c       uelema = array of names of the chemical elements read from
c                  the data file
c       uphasa = array of names of phases read from the data file
c       uptypa = array giving the type of a given phase (e.g.,
c                  liquid, solid, gas) read from the data file
c       uspeca = array of names of species read from the data file
c       utitld = the title read from the data file
c       vosp0a = array of standard state molar volumes of the species
c                  read from the data file
c       zchara = array of electrical charge numbers of the species
c                  read from the data file
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer iapxa_asv,ibpxa_asv,ikta_asv,ipbt_asv,ipch_asv,ipcv_asv,
     $ jpfc_asv,napa_asv,narx_asv,nata_asv,nbta_asv,nbta1_asv,ncta_asv,
     $ ndrsa_asv,nessa_asv,ngta_asv,nlta_asv,nmta_asv,nmuta_asv,
     $ npta_asv,nslta_asv,nsta_asv,ntid_asv,ntpr_asv,nxta_asv
c
      integer nad1,noutpt,nttyo
c
      integer iapxta(nxta_asv),ibpxta(nxta_asv),insgfa(nata_asv),
     $ jsola(nxta_asv),nalpaa(nslta_asv),narxt(ntpr_asv),
     $ nbaspa(nbta_asv),ncmpra(2,npta_asv),ndrsa(ndrsa_asv),
     $ ndrsra(2,nsta_asv),nessa(nessa_asv),nessra(2,nsta_asv),
     $ nmuxa(3,nmuta_asv),nslxa(2,nslta_asv)
c
      integer iaqsla,ielam,igas,ipch,ipcv,ixrn1a,ixrn2a,jpdblo,jptffl,
     $ napta,narn1a,narn2a,nata,nbta,nbtafd,ncta,ngrn1a,ngrn2a,ngta,
     $ nlrn1a,nlrn2a,nlta,nmrn1a,nmrn2a,nmta,nmuta,npta,nslta,nsta,
     $ ntitld,ntprt,nxrn1a,nxrn2a,nxta
c
      logical qclnsa(nsta_asv)
c
      character(len=8) uelema(ncta_asv)
      character(len=24) uphasa(npta_asv),uptypa(npta_asv)
      character(len=48) ubasp(nbta_asv),uspeca(nsta_asv)
      character(len=80) utitld(ntid_asv)
c
      character(len=8) udakey,udatfi
c
      real(8) tempcu(ntpr_asv)
c
      real(8) aadh(narx_asv,ntpr_asv),aadhh(narx_asv,ntpr_asv),
     $ aadhv(narx_asv,ntpr_asv),aaphi(narx_asv,ntpr_asv),
     $ abdh(narx_asv,ntpr_asv),abdhh(narx_asv,ntpr_asv),
     $ abdhv(narx_asv,ntpr_asv),abdot(narx_asv,ntpr_asv),
     $ abdoth(narx_asv,ntpr_asv),abdotv(narx_asv,ntpr_asv),
     $ adadhh(narx_asv,ntpr_asv,ipch_asv),
     $ adadhv(narx_asv,ntpr_asv,ipcv_asv),
     $ adbdhh(narx_asv,ntpr_asv,ipch_asv),
     $ adbdhv(narx_asv,ntpr_asv,ipcv_asv),
     $ adbdth(narx_asv,ntpr_asv,ipch_asv),
     $ adbdtv(narx_asv,ntpr_asv,ipcv_asv),
     $ adhfe(narx_asv,ntpr_asv,ipch_asv),
     $ advfe(narx_asv,ntpr_asv,ipcv_asv),
     $ adhfsa(narx_asv,ntpr_asv,ipch_asv,nsta_asv),
     $ advfsa(narx_asv,ntpr_asv,ipcv_asv,nsta_asv),
     $ aprehw(narx_asv,ntpr_asv),apresg(narx_asv,ntpr_asv),
     $ apxa(iapxa_asv,nxta_asv),atwta(ncta_asv),
     $ axhfe(narx_asv,ntpr_asv),axhfsa(narx_asv,ntpr_asv,nsta_asv),
     $ axlke(narx_asv,ntpr_asv),axlksa(narx_asv,ntpr_asv,nsta_asv),
     $ axvfe(narx_asv,ntpr_asv),axvfsa(narx_asv,ntpr_asv,nsta_asv)
c
      real(8) amua(jpfc_asv,nmuta_asv),
     $ aslma(jpfc_asv,0:ipbt_asv,nslta_asv),azeroa(nata_asv),
     $ bpxa(ibpxa_asv,nxta_asv),cco2(5),cdrsa(ndrsa_asv),
     $ cessa(nessa_asv),mwtspa(nsta_asv),palpaa(ipbt_asv,napa_asv),
     $ vosp0a(nsta_asv),zchara(nsta_asv)
c
      real(8) eps100
c
      real(8) tdamax,tdamin
c
c-----------------------------------------------------------------------
c
c     Local variables with global dimensioning.
c     None of these variables need to be SAVEd.
c
      character(len=24), dimension(:), allocatable :: udrsv
      character(len=8), dimension(:), allocatable :: uessv
c
      real(8), dimension(:), allocatable :: cdrsv,cessv
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ipc,j2,j3,n,nc,ndrsn,nerr,nessn,nmax,nn,np,ns,nsc,
     $ ntpr,nwatra
c
      integer ilnobl
c
      logical qx
c
      character(len=80) ux80
      character(len=56) ux56
      character(len=24) uphasv,uwater,uaqsln,uptsld,uptliq,uptgas,
     $ usblkf,ux24
      character(len=8) uendit,usedh,upitz,ustr,ustr2,ustr3,uterm
c
c-----------------------------------------------------------------------
c
      data uwater /'H2O                     '/
      data uaqsln /'Aqueous solution        '/
      data uptsld /'Solid                   '/
      data uptliq /'Liquid                  '/
      data uptgas /'Gas                     '/
c
      data uendit /'endit.  '/,uterm  /'+-------'/,
     $     usedh  /'SEDH    '/,upitz  /'Pitzer  '/
c
c-----------------------------------------------------------------------
c
c     Allocate scratch arrays.
c
      ALLOCATE(cessv(ncta_asv))
      ALLOCATE(uessv(ncta_asv))
c
      ALLOCATE(cdrsv(nbta1_asv))
      ALLOCATE(udrsv(nbta1_asv))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Zero some arrays.
c
      nmax = nbta_asv
      call initiz(nbaspa,nmax)
c
      nmax = ndrsa_asv
      call initiz(ndrsa,nmax)
      call initaz(cdrsa,nmax)
c
      nmax = nsta_asv
      call initaz(vosp0a,nmax)
c
      nmax = 2*nsta_asv
      call initiz(ndrsra,nmax)
c
      nmax = narx_asv*ntpr_asv*nsta_asv
      call initaz(axhfsa,nmax)
      call initaz(axlksa,nmax)
      call initaz(axvfsa,nmax)
c
c     SEDH arrays.
c
      nmax = nata_asv
      call initaz(azeroa,nmax)
      call initiz(insgfa,nmax)
c
c     Pitzer arrays.
c
      nmax = ipbt_asv*napa_asv
      call initaz(palpaa,nmax)
      nmax = jpfc_asv*(ipbt_asv + 1)*nslta_asv
      call initaz(aslma,nmax)
      nmax = jpfc_asv*nmuta_asv
      call initaz(amua,nmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Don't rewind the DATA1 file (nad1). It should be positioned
c     properly by the preceding call to subroutine indath.f.
c
      write (noutpt,1000)
      write (nttyo,1000)
 1000 format(/' Reading the rest of the DATA1 file ...')
c
c     Initialize the error counter.
c
      nerr = 0
c
c     Set the number of elements on the data file. This equal to
c     the dimensioned limit. Note that nbta, the number of basis
c     species, is determined in this subroutine as it is not usually
c     equal to the corresponding dimensioned limit nbta_asv.
c     In fact, nbta as returned by the subroutine may be subsequently
c     increased due to "promotions" implied by the contents of the
c     input file.
c
      ncta = ncta_asv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the data file title.
c
      ntitld = 0
      do n = 1,ntid_asv
        read (nad1) utitld(n)
        ntitld = ntitld + 1
        if (utitld(n)(1:8) .eq. uterm(1:8)) go to 140
      enddo
  140 continue
c
c     Write the first 5 lines of the data file title to the screen
c     file.
c
      write (nttyo,1060)
 1060 format(/3x,'The data file title is (first 5 lines maximum):',/)
c
      nn = min(ntitld,5)
      do n = 1,nn
        if (utitld(n)(1:8) .eq. uterm(1:8)) go to 150
        j2 = ilnobl(utitld(n))
        j2 = min(j2,74)
        write (nttyo,1070) utitld(n)(1:j2)
 1070   format(5x,a)
      enddo
  150 continue
c
c     Write the complete data file title, less terminator, if any,
c     to the output file.
c
      write (noutpt,1080)
 1080 format(/3x,'The data file title is:',/)
c
      do n = 1,ntitld
        if (utitld(n)(1:8) .eq. uterm(1:8)) go to 160
        j2 = ilnobl(utitld(n))
        j2 = min(j2,74)
        write (noutpt,1070) utitld(n)(1:j2)
      enddo
  160 continue
c
c     Search the data file title for the usual identifying three-letter
c     keystring (e.g., com, hmw, etc.).
c
      udatfi = 'current'
      do n = 1,ntitld
        ux80 = utitld(n)
        call lejust(ux80)
        if (ux80(1:6) .eq. 'data0.') then
          udatfi = ux80(7:9)
          go to 170
        endif
      enddo
  170 continue
c
      write (noutpt,1090)
      write (nttyo,1090)
 1090 format(/3x,'Continuing to read the DATA1 file ...')
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the actual number of temperature ranges in the temperature
c     grid.
c
      read (nad1) ux56
      read (nad1) ntprt
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the actual number of coefficients in each range.
c
      do ntpr = 1,ntprt
        read (nad1) ux56
        read (nad1) narxt(ntpr)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read a flag indicating the degree of representation of enthalpy
c     functions.
c
      read (nad1) ux56
      read (nad1) ipch
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read a flag indicating the degree of representation of volume
c     functions.
c
      read (nad1) ux56
      read (nad1) ipcv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (udakey(1:8) .eq. upitz(1:8)) then
c
c       Read the Pitzer data block organization flag.
c
        read (nad1) ux56
        read (nad1) jpdblo
c
c       Read the Pitzer parameter temperature function flag.
c
        read (nad1) ux56
        read (nad1) jptffl
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read chemical elements data.
c
      read (nad1) ux24
c
      do nc = 1,ncta
        read (nad1) uelema(nc),atwta(nc)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read nominal temperature limits (Celsius) for the data file.
c
      read (nad1) ux56
      read (nad1) tdamin,tdamax
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the maximum temperatures (tempcu) of the temperature ranges.
c     These are used to find the temperature range flag (ntpr) for any
c     specified temperature (see EQLIB/gntpr.f).
c
      read (nad1) ux56
      do n = 1,ntprt
        read (nad1) tempcu(n)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read interpolating polynomial coefficients for computing the
c     standard grid pressure (bars) as a function of temperature
c     (Celsius). Often this grid corresponds to 1.013 bar up to 100C
c     and the steam-liquid water saturation pressure at higher
c     temperatures.
c
c     Calling sequence substitutions:
c       apresg for arr
c
      call indatc(apresg,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ipcv .ge. 0) then
c
c       Read interpolating polynomial coefficients for computing the
c       recommended pressure envelope half width (bars).
c
c       Calling sequence substitutions:
c         aprehw for arr
c
       call indatc(aprehw,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (udakey(1:8) .eq. usedh(1:8)) then
c
c       The aqeuous species activity coefficient formalism is
c       consistent with the B-dot equation, the Davies equation,
c       or similar equations.
c
c       Read interpolating polynomial coefficients for the Debye-Huckel
c       A(gamma,10) parameter and related parameters.
c
c       Calling sequence substitutions:
c         aadh for arr
c
        call indatc(aadh,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)
c
        if (ipch .ge. 0) then
c
c         Calling sequence substitutions:
c           aadhh for arr
c
          call indatc(aadhh,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)
c
          do ipc = 1,ipch
c
c           Calling sequence substitutions:
c             adadhh for arr
c             ipch_asv for ipcx_asv
c
            call indatd(adadhh,ipc,ipch_asv,nad1,narxt,narx_asv,ntprt,
     $      ntpr_asv,ux24)
          enddo
        endif
c
        if (ipcv .ge. 0) then
c
c         Calling sequence substitutions:
c           aadhv for arr
c
          call indatc(aadhv,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)
c
          do ipc = 1,ipcv
c
c           Calling sequence substitutions:
c             adadhv for arr
c             ipcv_asv for ipcx_asv
c
            call indatd(adadhv,ipc,ipcv_asv,nad1,narxt,narx_asv,ntprt,
     $      ntpr_asv,ux24)
          enddo
        endif
c
c       Read interpolating polynomial coefficients for the Debye-Huckel
c       B(gamma) and related parameters.
c
c       Calling sequence substitutions:
c         abdh for arr
c
        call indatc(abdh,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)
c
        if (ipch .ge. 0) then
c
c         Calling sequence substitutions:
c           abdhh for arr
c
          call indatc(abdhh,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)
c
          do ipc = 1,ipch
c
c           Calling sequence substitutions:
c             adbdhh for arr
c             ipch_asv for ipcx_asv
c
            call indatd(adbdhh,ipc,ipch_asv,nad1,narxt,narx_asv,ntprt,
     $      ntpr_asv,ux24)
          enddo
        endif
c
        if (ipcv .ge. 0) then
c
c         Calling sequence substitutions:
c           abdhv for arr
c
          call indatc(abdhv,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)
c
          do ipc = 1,ipcv
c
c           Calling sequence substitutions:
c             adbdhv for arr
c             ipcv_asv for ipcx_asv
c
            call indatd(adbdhv,ipc,ipcv_asv,nad1,narxt,narx_asv,ntprt,
     $      ntpr_asv,ux24)
          enddo
        endif
c
c       Read interpolating polynomial coefficients for the B-dot
c       parameter and related parameters.
c
c       Calling sequence substitutions:
c         abdot for arr
c
        call indatc(abdot,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)
c
        if (ipch .ge. 0) then
c
c         Calling sequence substitutions:
c           abdoth for arr
c
          call indatc(abdoth,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)
c
          do ipc = 1,ipch
c
c           Calling sequence substitutions:
c             adbdth for arr
c             ipch_asv for ipcx_asv
c
            call indatd(adbdth,ipc,ipch_asv,nad1,narxt,narx_asv,ntprt,
     $      ntpr_asv,ux24)
          enddo
        endif
c
        if (ipcv .ge. 0) then
c
c         Calling sequence substitutions:
c           abdotv for arr
c
          call indatc(abdotv,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)
c
          do ipc = 1,ipcv
c
c           Calling sequence substitutions:
c             adbdtv for arr
c             ipcv_asv for ipcx_asv
c
            call indatd(adbdtv,ipc,ipcv_asv,nad1,narxt,narx_asv,ntprt,
     $      ntpr_asv,ux24)
          enddo
        endif
c
c       Read coefficients for the Drummond (1981) equation.
c       This equation represents the activity coeficient
c       of CO2(aq) as a function of temperature and ionic
c       strength in sodium chloride solutions.
c
        read (nad1) ux24
        read (nad1) (cco2(n), n = 1,5)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (udakey(1:8) .eq. upitz(1:8)) then
c
c       The aqeuous species activity coefficient formalism is
c       consistent with Pitzer's equations. Read interpolating
c       polynomial coefficients for Debye-Huckel A(phi) and
c       related parameters.
c
c       Calling sequence substitutions:
c         aaphi for arr
c
        call indatc(aaphi,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)
c
        if (ipch .ge. 0) then
c
c         Calling sequence substitutions:
c           aadhh for arr
c
          call indatc(aadhh,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)
c
          do ipc = 1,ipch
c
c           Calling sequence substitutions:
c             adadhh for arr
c             ipch_asv for ipcx_asv
c
            call indatd(adadhh,ipc,ipch_asv,nad1,narxt,narx_asv,ntprt,
     $      ntpr_asv,ux24)
          enddo
        endif
c
        if (ipcv .ge. 0) then
c
c         Calling sequence substitutions:
c           aadhv for arr
c
          call indatc(aadhv,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)
c
          do ipc = 1,ipcv
c
c           Calling sequence substitutions:
c             adadhv for arr
c             ipcv_asv for ipcx_asv
c
            call indatd(adadhv,ipc,ipcv_asv,nad1,narxt,narx_asv,ntprt,
     $      ntpr_asv,ux24)
          enddo
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read interpolating polynomial coefficients for the log K of
c     the "Eh" reaction:
c
c       2 H2O(l) = 2 O2(g) + 4 H+ + 4 e-
c
c     Calling sequence substitutions:
c       axlke for arr
c
      call indatc(axlke,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)
c
      if (ipch .ge. 0) then
c
c       Calling sequence substitutions:
c         axhfe for arr
c
        call indatc(axhfe,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)
c
        do ipc = 1,ipch
c
c         Calling sequence substitutions:
c           adhfe for arr
c           ipch_asv for ipcx_asv
c
          call indatd(adhfe,ipc,ipch_asv,nad1,narxt,narx_asv,ntprt,
     $    ntpr_asv,ux24)
        enddo
      endif
c
      if (ipcv .ge. 0) then
c
c       Calling sequence substitutions:
c         axvfe for arr
c
        call indatc(axvfe,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)
c
        do ipc = 1,ipcv
c
c         Calling sequence substitutions:
c           advfe for arr
c           ipcv_asv for ipcx_asv
c
          call indatd(advfe,ipc,ipcv_asv,nad1,narxt,narx_asv,ntprt,
     $    ntpr_asv,ux24)
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the phases and species and assocated data from the
c     data file. The following variables are used as counters:
c
c       np     = phase index
c       ns     = species index
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize the qclnsa array to .false.
c
      qx = .false.
c
      call initlv(qclnsa,nsta_asv,qx)
c
c     Initialize the following variables, all of which are used
c     below as counters.
c
c       nessn  = counter for the cessa and nessa arrays
c       ndrsn  = counter for the cdrsa and ndrsa arrays
c       ns     = species index
c       nbta   = number of basis species
c       nata   = number of aqueous species
c       ngta   = number of gas species
c       nmta   = number of pure mineral species
c       nlta   = number of pure liquid species
c       nxta   = number of solid solution phases
c
      nessn = 0
      ndrsn = 0
      ns = 0
      nbta = 0
      nata = 0
      ngta = 0
      nmta = 0
      nlta = 0
      nxta = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Process the aqueous species superblock. Set up the aqueous
c     solution as the first phase.
c
      np = 1
      uphasv(1:24) = uaqsln(1:24)
      iaqsla = np
      uphasa(np)(1:24) = uphasv(1:24)
      uptypa(np)(1:24) = uaqsln(1:24)
      usblkf(1:24) = uaqsln(1:24)
      ncmpra(1,1) = 1
      narn1a = 1
c
c     Read the superblock header.
c
      read (nad1) ustr,ustr2,ustr3
c
c     Read the blocks in the current superblock.
c
      call indats(adhfsa,advfsa,axhfsa,axlksa,axvfsa,cdrsa,cdrsv,
     $ cessa,cessv,ipch,ipch_asv,ipcv,ipcv_asv,mwtspa,nad1,narxt,
     $ narx_asv,nata,nata_asv,nbta,nbta_asv,nbta1_asv,nbtafd,ncmpra,
     $ ncta,ncta_asv,ndrsa,ndrsa_asv,ndrsn,ndrsra,nerr,nessa,nessa_asv,
     $ nessn,nessra,ngta,ngta_asv,nlta,nlta_asv,nmta,nmta_asv,noutpt,
     $ np,npta_asv,ns,nsta_asv,ntprt,ntpr_asv,nttyo,uaqsln,ubasp,udrsv,
     $ uelema,uendit,uessv,uphasa,uphasv,uptgas,uptliq,uptsld,uptypa,
     $ usblkf,uspeca,vosp0a,zchara)
c
      ncmpra(2,np) = ns
      narn2a = ns
c
c     Check to see that water is the first species in the phase
c     "aqueous solution". This is required so that all solute species
c     can be cleanly referenced in simple loop; that is, in a loop of
c     the form "do ns = narn1a + 1,narn2a".
c
      nwatra = 0
      if (uspeca(narn1a)(1:24) .eq. uwater(1:24)) nwatra = 1
c
      if (nwatra .eq. 0) then
        j3 = ilnobl(uwater)
        write (noutpt,1150) uwater(1:j3),uaqsln(1:j2)
        write (nttyo,1150) uwater(1:j3),uaqsln(1:j2)
 1150   format(/' * Error - (EQLIB/indata) The species ',
     $  /7x,a,' (',a,') must be listed on the',
     $  /7x,'data file as the first species in its phase.')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Process the pure minerals superblock. Each pure mineral is both
c     a species and a phase. The phase name variable uphasv will be set
c     to the name of each mineral as its block is read.
c
      usblkf(1:24) = uptsld(1:24)
      nmrn1a = ns + 1
c
c     Read the superblock header.
c
      read (nad1) ustr,ustr2,ustr3
c
c     Read the blocks in the current superblock.
c
      call indats(adhfsa,advfsa,axhfsa,axlksa,axvfsa,cdrsa,cdrsv,
     $ cessa,cessv,ipch,ipch_asv,ipcv,ipcv_asv,mwtspa,nad1,narxt,
     $ narx_asv,nata,nata_asv,nbta,nbta_asv,nbta1_asv,nbtafd,ncmpra,
     $ ncta,ncta_asv,ndrsa,ndrsa_asv,ndrsn,ndrsra,nerr,nessa,nessa_asv,
     $ nessn,nessra,ngta,ngta_asv,nlta,nlta_asv,nmta,nmta_asv,noutpt,
     $ np,npta_asv,ns,nsta_asv,ntprt,ntpr_asv,nttyo,uaqsln,ubasp,udrsv,
     $ uelema,uendit,uessv,uphasa,uphasv,uptgas,uptliq,uptsld,uptypa,
     $ usblkf,uspeca,vosp0a,zchara)
c
      nmrn2a = ns
      if (nmrn2a .le. narn2a) nmrn1a = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     There is presently no superblock corresponding to pure liquids.
c     Set up water as a pure liquid.
c
      usblkf(1:24) = uptliq(1:24)
      nlrn1a = 0
      nlrn2a = 0
      nlta = 0
c
      if (nwatra .gt. 0) then
        np = np + 1
        ns = ns + 1
        uphasv(1:24) = uspeca(nwatra)(1:24)
        uphasa(np)(1:24) = uphasv(1:24)
        uptypa(np)(1:24) = uptliq(1:24)
        ncmpra(1,np) = ns
        ncmpra(2,np) = ns
        nsc = 1
        call clones(axlksa,cdrsa,cessa,mwtspa,narx_asv,
     $  ndrsa,ndrsa_asv,ndrsn,ndrsra,nessa,nessa_asv,nessn,
     $  nessra,np,npta_asv,ns,nsc,nsta_asv,ntpr_asv,uphasa,
     $  uspeca,zchara)
        qclnsa(ns) = .true.
        nlrn1a = ns
        nlrn2a = ns
        nlta = 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Process the gas species superblock. Set up the gas phase
c     as the next phase.
c
      uphasv(1:24) = uptgas(1:24)
      usblkf(1:24) = uptgas(1:24)
      np = np + 1
      igas = np
      uphasa(np)(1:24) = uphasv(1:24)
      uptypa(np)(1:24) = uptgas(1:24)
      ncmpra(1,np) = ns + 1
      ngrn1a = ns + 1
c
c     Read the superblock header.
c
      read (nad1) ustr,ustr2,ustr3
c
c     Read the blocks in the current superblock.
c
      call indats(adhfsa,advfsa,axhfsa,axlksa,axvfsa,cdrsa,cdrsv,
     $ cessa,cessv,ipch,ipch_asv,ipcv,ipcv_asv,mwtspa,nad1,narxt,
     $ narx_asv,nata,nata_asv,nbta,nbta_asv,nbta1_asv,nbtafd,ncmpra,
     $ ncta,ncta_asv,ndrsa,ndrsa_asv,ndrsn,ndrsra,nerr,nessa,nessa_asv,
     $ nessn,nessra,ngta,ngta_asv,nlta,nlta_asv,nmta,nmta_asv,noutpt,
     $ np,npta_asv,ns,nsta_asv,ntprt,ntpr_asv,nttyo,uaqsln,ubasp,udrsv,
     $ uelema,uendit,uessv,uphasa,uphasv,uptgas,uptliq,uptsld,uptypa,
     $ usblkf,uspeca,vosp0a,zchara)
c
      ngrn2a = ns
      if (ngrn2a .le. nlrn2a) ngrn1a = 0
      ncmpra(2,np) = ns
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Process the solid solution superblock.
c
      read (nad1) ustr,ustr2,ustr3
      nxrn1a = ns + 1
      ixrn1a = np + 1
c
      call indatp(apxa,axlksa,bpxa,cdrsa,cessa,iapxa_asv,iapxta,
     $ ibpxa_asv,ibpxta,ikta_asv,jsola,mwtspa,nad1,narx_asv,ncmpra,
     $ ndrsa,ndrsa_asv,ndrsn,ndrsra,nerr,nessa,nessa_asv,nessn,
     $ nessra,nmrn1a,nmrn2a,noutpt,np,npta_asv,ns,nsta_asv,ntpr_asv,
     $ nttyo,nxta,nxta_asv,qclnsa,uendit,uspeca,uphasa,uptsld,
     $ uptypa,zchara)
c
      npta = np
      nsta = ns
      nxrn2a = ns
      ixrn2a = np
      if (nxrn2a .le. nlrn2a) then
        nxrn1a = 0
        nxrn2a = 0
        ixrn1a = 0
        ixrn2a = 0
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set up the nbaspa array and expand the elements of the ubasp array
c     to the full 48 characters. Convert basis species indices in the
c     ndrsa array to species indices.
c
      call dfbasp(nbaspa,nbta,nbta_asv,ndrsa,ndrsa_asv,ndrsra,
     $ nerr,noutpt,nsta,nsta_asv,nttyo,ubasp,uspeca)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (udakey(1:8) .eq. usedh(1:8)) then
c
c       The aqeuous species activity coefficient formalism is
c       consistent with the B-dot equation, the Davies equation,
c       or similar equations. Read the appropriate data for
c       this formalism.
c
        call inbdot(azeroa,insgfa,nad1,narn1a,narn2a,nata,
     $  nata_asv,nerr,noutpt,nsta_asv,nttyo,uspeca)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (udakey(1:8) .eq. upitz(1:8)) then
c
c       The aqeuous species activity coefficient formalism is
c       consistent with Pitzer's equations. Read the appropriate
c       data for this formalism.
c
        call inupt(amua,aslma,ielam,ipbt_asv,jpdblo,jpfc_asv,
     $  nad1,nalpaa,napa_asv,napta,narn1a,narn2a,nerr,nmuta,
     $  nmuta_asv,nmuxa,noutpt,nslta,nslta_asv,nslxa,nsta_asv,
     $  nttyo,palpaa,uspeca,zchara)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The following is a bit of nonsense so compiler warnings will
c     not be generated saying that ux80, ustr, ustr2, and ustr3 are
c     not used (they are all used to read unused data).
c
      ux80(1:56) = ux56
      ustr3 = ux80(1:8)
      ustr = ustr2
      ustr2 = ustr3
      ustr3 = ustr
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Deallocate the scratch arrays used by calls to EQLIB/indats.f.
c
      DEALLOCATE(cessv)
      DEALLOCATE(uessv)
c
      DEALLOCATE(cdrsv)
      DEALLOCATE(udrsv)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nerr .gt. 0) stop
c
      write (noutpt,1200)
      write (nttyo,1200)
 1200 format(/'   Done reading the DATA1 file.')
c
      end
