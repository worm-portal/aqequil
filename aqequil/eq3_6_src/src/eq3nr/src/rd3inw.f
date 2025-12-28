      subroutine rd3inw(cgexj,cgexpi,covali,ehi,egexsi,fo2lgi,iebal3,
     $ ietmax,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,jgext,jetmax,
     $ jflgi,jgexti,jpres3,kxmod,mwtges,nbti,nbtmax,ncmpri,net,neti,
     $ netmax,ngexti,ninpts,ngexrt,nobswt,nodbmx,nopgmx,noprmx,
     $ noptmx,noutpt,nprob,nsbswt,ntitl,ntitmx,nttyo,nxicmx,nxmdmx,
     $ nxmod,nxti,nxtimx,pei,press,qend,qgexsh,qrderr,rho,scamas,
     $ tdspkg,tdspl,tempc,tgexp,tolbt,toldl,tolspf,ucospi,uebal,
     $ ugexj,ugexji,ugexmo,ugexp,ugexpi,ugexr,ugexsi,umemi,uobsw,
     $ uredox,usbsw,usoli,uspeci,utitl,uhfgex,uvfgex,uxkgex,uxmod,
     $ xbari,xgexsi,xhfgex,xlkgex,xvfgex,xlkmod,zgexj)
c
c     This subroutine reads the EQ3NR input file in compact ("W")
c     format.
c
c     This subroutine is a near-clone of XCON3/rd3w8.f. The present
c     subroutine differs from the latter, in that in addition to the
c     pure read function, it writes an instant echo of what is read to
c     the output file. To maintain close consistency with XCON3/rd3w8.f,
c     this subroutine assigns no default values and only performs such
c     checking of what is read to ensure that what follows is readable.
c
c     The calling sequence of this subroutine is identical to that of
c     EQ3NR/rd3ind.f, XCON3/rd3w8.f, and XCON3/rd3d8.f.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ninpts = unit number of the stripped input file
c       noutpt = unit number of the output file
c       nttyo  = unit number of the screen file
c
c     Principal output:
c
c       qrderr = flag denoting a read error
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ietmax,jetmax,nbtmax,netmax,nodbmx,nopgmx,noprmx,noptmx,
     $ ntitmx,nxmdmx,nxicmx,nxtimx
c
      integer ninpts,noutpt,nttyo
c
      integer iodb(nodbmx),iopg(nopgmx),iopr(noprmx),iopt(noptmx),
     $ jgext(netmax),jgexti(netmax),jflgi(nbtmax),kxmod(nxmdmx),
     $ ncmpri(2,nxtimx),ngexrt(jetmax,netmax),ngexti(jetmax,netmax)
c
      integer iebal3,irdxc3,itdsf3,itermx,jpres3,nbti,net,neti,nobswt,
     $ nprob,nsbswt,ntitl,nxmod,nxti
c
      logical qend,qgexsh,qrderr
c
      character(len=80) utitl(ntitmx)
      character(len=56) ugexr(ietmax,jetmax,netmax)
      character(len=48) ucospi(nbtmax),uobsw(2,nbtmax),usbsw(2,nbtmax),
     $ uspeci(nbtmax),uxmod(nxmdmx)
      character(len=24) ugexmo(netmax),ugexp(netmax),ugexpi(netmax),
     $ ugexsi(ietmax,jetmax,netmax),umemi(nxicmx),usoli(nxtimx)
      character(len=24) uebal,uredox
      character(len=8) ugexj(jetmax,netmax),ugexji(jetmax,netmax),
     $ uhfgex(ietmax,jetmax,netmax),uvfgex(ietmax,jetmax,netmax),
     $ uxkgex(ietmax,jetmax,netmax)
c
      real(8) cgexj(jetmax,netmax),cgexpi(netmax),covali(nbtmax),
     $ egexsi(ietmax,jetmax,netmax),mwtges(netmax),tgexp(netmax),
     $ xbari(nxicmx),xhfgex(ietmax,jetmax,netmax),
     $ xlkgex(ietmax,jetmax,netmax),xvfgex(ietmax,jetmax,netmax),
     $ xlkmod(nxmdmx),zgexj(jetmax,netmax),xgexsi(ietmax,jetmax,netmax)
c
      real(8) ehi,fo2lgi,pei,press,rho,scamas,tdspkg,tdspl,tempc,tolbt,
     $ toldl,tolspf
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,iei,je,jei,j2,j3,j4,n,nbi,ne,nei,nn,nr1,nr2,nssoti,nxi,
     $ nxic
c
      integer ilnobl
c
      logical qgexef
c
      character(len=80) uline
      character(len=48) ux48
      character(len=8) uendit,ux8
c
c-----------------------------------------------------------------------
c
      data uendit /'endit.  '/
c
c-----------------------------------------------------------------------
c
c     qend   = .true if the end of the input file has been encountered
c     qrderr = .true if the current problem can't be read because of
c                a read format error or a dimensional overflow
c
      qend = .false.
      qrderr = .false.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Title.
c
c       utitl1 = Title
c       ntitl1 = number of lines in the title
c
c     Note: if there are exactly ntitpa lines in the title, the title
c     need not be terminated by an 'endit.'. The 'endit.' if present
c     is here not considered to be part of the title.
c
      n = 0
      read (ninpts,1000,end=100,err=990) uline
 1000 format(a80)
      go to 110
c
  100 qend = .true.
      go to 999
c
  110 write (noutpt,1010) nprob
      write (nttyo,1010) nprob
 1010 format(//' Reading problem ',i3,' from the input file ...',/)
c
      j2 = ilnobl(uline)
      j2 = min(j2,79)
      write (noutpt,1020) uline
 1020 format(1x,a)
c
      if (uline(1:8) .eq. uendit(1:8)) then
        write (noutpt,1030)
        write (nttyo,1030)
 1030   format(3x,'This input file has no title.')
        go to 120
      else
        n = 1
        utitl(n) = uline
      endif
c
      do nn = 2,ntitmx + 1
        read (ninpts,1000,err=990) uline
        j2 = ilnobl(uline)
        j2 = min(j2,79)
        write (noutpt,1020) uline(1:j2)
        if (uline(1:8) .eq. uendit(1:8)) go to 120
c
        n = n + 1
c
        if (n .gt. ntitmx) then
          write (nttyo,1015) ntitmx
 1015     format(/' * Error - (EQ3NR/rd3inw) Have too many lines in',
     $    /7x,'the title. The code is only dimensioned for ',i4,
     $    /7x,'lines. Reduce the size of the title or increase the',
     $    /7x,'dimensioning parameter ntitpa.')
          go to 990
        endif
c
        utitl(n) = uline
      enddo
  120 ntitl = n
c
      if (ntitl .gt. 0) then
c
c       Write the first 5 lines of the input problem title to the
c       screen file.
c
        write (nttyo,1040)
 1040   format(3x,'The input problem title is (first 5 lines',
     $   ' maximum):',/)
c
        i = min(ntitl,5)
        do n = 1,i
          j2 = ilnobl(utitl(n))
          j2 = min(j2,74)
          write (nttyo,1050) utitl(n)(1:j2)
 1050     format(5x,a)
        enddo
      endif
c
      write (nttyo,1060)
 1060 format(/3x,'Continuing to read the problem input ...')
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Special basis switches.
c
c       nsbswt = the number of special basis switches
c       usbsw(1,n) = the species to be switched from the strict
c                      basis to the auxiliary basis set in the
c                      in the n-th switch
c       usbsw(2,n) = the species to be switched from the auxiliary
c                      basis set to the strict basis set in the
c                      n-th switch
c
      write (noutpt,1100)
 1100 format(' * Special basis switches')
c
      read (ninpts,1110,err=990) nsbswt
 1110 format(12x,i3)
      write (noutpt,1120) nsbswt
 1120 format(5x,'nsbswt= ',i3)
c
      do n = 1,nsbswt
        read (ninpts,1130,err=990) usbsw(1,n)
 1130   format(9x,a48)
        j2 = ilnobl(usbsw(1,n))
        write (noutpt,1140) usbsw(1,n)(1:j2)
 1140   format(1x,'species= ',a)
c
        read (ninpts,1150,err=990) usbsw(2,n)
 1150   format(15x,a48)
        j2 = ilnobl(usbsw(2,n))
        write (noutpt,1160) usbsw(2,n)(1:j2)
 1160   format(3x,'switch with= ',a)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Temperature.
c
c       tempc  = temperature, C
c
      write (noutpt,1165)
 1165 format(' * General')
c
      read (ninpts,1170,err=990) tempc
 1170 format(12x,e12.5)
      write (noutpt,1180) tempc
 1180 format(6x,'tempc= ',1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Pressure.
c
c       jpres3 = pressure code:
c
c         0 = Compute the pressure from the temperauture, using the
c             data file's reference pressure curve:
c               press = func(data file, T)
c
c         1 = Compute the pressure from the temperature, using the
c             1.013-bar/steam-saturation curve:
c               press = built-in func(T)
c
c         2 = Constant pressure:
c               press = pressb
c
c       press  = pressure, bars
c
      read (ninpts,1190,err=990) jpres3,press
 1190 format(12x,i3,/12x,e12.5)
      write (noutpt,1200) jpres3,press
 1200 format(5x,'jpres3= ',i3,/6x,'press= ',1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Density.
c
c       rho    = density of aqueous solution in g/ml
c
      read (ninpts,1220,err=990) rho
 1220 format(12x,e12.5)
      write (noutpt,1222) rho
 1222 format(8x,'rho= ',1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Total dissolved solutes.
c
c       itdsf3 = TDS code:
c
c         0 = Use the input value in mg/kg.sol (tdspkg)
c         1 = Use the input value in mg/L (tdspl)
c
c       tdspkg = total dissolved solutes, mg/kg.sol
c       tdspl  = total dissolved solutes, mg/L
c
      read (ninpts,1224,err=990) itdsf3,tdspkg,tdspl
 1224 format(12x,i3,/12x,e12.5,12x,e12.5)
      write (noutpt,1227) itdsf3,tdspkg,tdspl
 1227 format(5x,'itdsf3= ',i3,/5x,'tdspkg= ',1pe12.5,5x,'tdspl= ',e12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Species to adjust for electrical balance.
c
c       iebal3 = electrical balancing code:
c
c         0 = No balancing
c         1 = Balance on the species whose name is entered as uredox
c
c       uebal = name of a basis species whose concentration is to be
c                 adjusted to achieve electrical balance
c
      read (ninpts,1500,err=990) iebal3,uebal
 1500 format(12x,i3,/12x,a24)
      j2 = ilnobl(uebal)
      write (noutpt,1510) iebal3,uebal(1:j2)
 1510 format(5x,'iebal3= ',i3/6x,'uebal= ',a)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Default redox constraint.
c
c       irdxc3 = default redox option switch:
c         -3 = Read a constraint on aqueous "O2(g)" in the regular
c                basis species constraint block
c         -2 = The pe is fixed at the value entered as "pei"
c         -1 = The Eh (V) is fixed at the value entered as "ehi"
c          0 = The log fO2 is fixed at the value entered as "fo2lgi"
c          1 = The log Fo2, Eh, pe, etc., is determined by a redox
c                couple, one member of which is defined by the name
c                of an auxiliary basis species entered as "uredox",
c                the member other of which is defined by its appearance
c                in the reaction for the former. Usually Fe++ is a
c                strict basis species, and Fe++ is in the auxiliary
c                basis set, so entering uredox = "Fe+++" specifies
c                the Fe+++/Fe++ couple. Concentration or other data
c                must be entered below for both members of the couple.
c                Couples involving non-aqueous phases can be defined
c                indirectly, using the jflgi = 25 option described
c                below. For example, to use the Hematite/Magnetite
c                couple, enter uredox = "Fe+++". In the basis species
c                constraint block, for Fe++, enter jflgi = 25 and
c                ucospi = "Magnetite", and for Fe+++, enter jflgi = 25
c                and ucospi = "Hematite".
c
c       fo2lgi = log fO2
c       ehi    = Eh, oxidation-reduction parameter
c       pei    = pe, electron activity
c       uredox = name (all 24 letters) of an auxiliary basis species
c                  that is part of a redox couple that will be used
c                  to  distribute other redox couples
c
c       Note: heterogeneous equilibria can be used also by setting
c         irdxc3 = 0, fep = 0., and jflgi(O2(g)) = 25
c
      read (ninpts,1232,err=990) irdxc3
 1232 format(12x,i3)
      write (noutpt,1234) irdxc3
 1234 format(5x,'irdxc3= ',i3)
      read (ninpts,1240,err=990) fo2lgi,ehi
 1240 format(12x,e12.5,12x,e12.5)
      write (noutpt,1242) fo2lgi,ehi
 1242 format(5x,'fo2lgi= ',1pe12.5,7x,'ehi= ',e12.5)
      read (ninpts,1250,err=990) pei,uredox
 1250 format(12x,e12.5,12x,a24)
      j2 = ilnobl(uredox)
      write (noutpt,1252) pei,uredox(1:j2)
 1252 format(8x,'pei= ',1pe12.5,4x,'uredox= ',a)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Aqueous basis species and associated constraints.
c
c       uspeci = name of an aqueous species in the basis set. If
c                  the species is not defined as a basis species on
c                  the input file, it will be treated as a basis species
c                  in the current run. The name is constructed in
c                  species/phase format (the first 24 characters
c                  comprise the species part,the second 24 characters,
c                  the name of the corresponding phase). Only the
c                  species part is input; the code provides the
c                  phase part.
c
c       jflgi  = flag variable defining the type of constraint placed
c                  on a basis species. It determines the meaning of the
c                  "covali" parameter.
c
c       covali = a floating point datum (concentration, activity, etc.)
c                  whose meaning is specified by the jflgi parameter.
c
c       ucospi = the name of a species/phase which must be designated
c                  to complete a specification for the jflgi = 17, 18,
c                  or 25 options. For jflgi = 17 or 18, the species is
c                  an aqueous ion. For jflgi = 25, the species
c                  specified must belong to a phase other than the
c                  aqueous solution. The name is constructed in
c                  species/phase format. For the jflgi = 17 or 18
c                  options, the phase part may be left blank, as the
c                  phase by implication is 'aqueous solution.'
c                  For the jflgi = 25 option, the phase part may also
c                  be left blank. The phase will then be taken as that
c                  which first produces a name match on the species
c                  part. If the phase is intended to be a solid or
c                  liquid solution, failure to specify the phase part
c                  will result in matching the pure component phase
c                  instead.
c
c       Summary of the jflgi conventions:
c
c         -1 = Not present in the system. This value is assigned
c                as a default to strict basis species.
c
c          0 = Total concentration, molal ("covali")
c                Note: enter jflgi = 0 for water and O2(g). The
c                "covali" parameter for O2(g) is always log fugacity.
c                A line for O2(g) is not entered here unless
c                irdxc3 = -3.
c
c          1 = Total concentration, molar ("covali")
c
c          2 = Total concentration, mg/L ("covali")
c
c          3 = Total concentration, mg/kg.sol ("covali")
c
c          7 = Total alkalinity, eq/kg.H2O ("covali")
c
c          8 = Total alkalinity, eq/L ("covali")
c
c          9 = Total alkalinity, eq/kg.sol ("covali")
c
c         10 = Total alkalinity, mg/L CaCO3 ("covali")
c
c         11 = Total alkalinity, mg/L HCO3- ("covali")
c
c         16 = Log activity (molal scale) ("covali")
c
c         17 = The function abs(zj) log ai +/- abs(zi) log aj,
c                where i and j are the ions specified by "uspeci"
c                and "ucospi". The sign of the second term is positive
c                for cation-anion and anion-cation combinations,
c                negative for like-charge combinations. To input the
c                pHCl defined as pH + pCl, enter H+ as "uspeci", -pHCl
c                for "covali", and Cl- for uphas1.
c
c         18 = Log mean activity of the neutral electrolyte composed
c                of the ions specified by "uspeci" and "ucospi".
c
c         19 = pX, or -log activity. This is very similar to the
c                jflgi= 16 option, except that the "covali" input
c                has the opposite sign. This can be used to enter
c                pH, but see jflag= 20 below.
c
c         20 = pH. This is very similar to the jflgi= 19 (pX) option.
c                It can only be applied to H+, however.
c
c         21 = pHCl. This is very similar to the jflgi= 17 (log
c                activity combination) option. It can only be applied
c                to H+ or Cl-, however. The other ion is then the
c                counterion. Note that the counterion is not specified
c                in the "ucospi" input.
c
c         25 = Heterogeneous equilibrium. Here "usp" contains the name
c                of the relevant species (first 24 characters) and
c                phase (last 24 characters). If the relevant species
c                is a gas species, its fugacity is specified by the
c                "covali" input. This is a dangerous option.
c
c         27 = Homogeneous equilibrium - the species must be in the
c                auxiliary basis when the calculation commences;
c                it and its own dependent species are treated as
c                depenent species of the corresponding species in the
c                strict basis, but are not limited by the total
c                concentration, if one is specified, for that
c                corresponding species. No "covali" input is required.
c                This is a dangerous option.
c
c                For example, suppose you have a sulfate analysis
c                of a water and you expect that sulfide species
c                may be significant and possibly in equilibrium with
c                sulfate. If you think the sulfate analysis
c                measured both sulfate and sulfide, enter the
c                analytical value as the total concentration of
c                sulfate and specify jflgi(HS-) = 30. If you think
c                the analysis measured only sulfate proper, then
c                enter jflgi(HS-) = 27 instead. If the redox state of
c                the water is very reducing you may get an infinite
c                or near-infinite amount of calculated sulfide.
c
c         16 = Log activity (molal scale) ("covali")
c
c         30 = Make non-basis. The species is treated as above,
c                but it and its own dependent species are constrained
c                by the total concentration entered for the
c                corresponding strict basis species. No "covali" input
c                is required.
c
      write (noutpt,2490)
 2490 format(' * Aqueous basis species')
c
      nbi = 0
c
      do nn = 1,nbtmax + 1
        read (ninpts,2500,err=990) ux8,ux48
 2500   format(a8,1x,a48)
        j2 = ilnobl(ux48)
c
        if (ux8(1:8) .eq. uendit(1:8)) then
          j3 = ilnobl(uendit)
          write (noutpt,2510) uendit(1:j3)
 2510     format(1x,a)
          go to 400
        endif
c
        nbi = nbi + 1
c
        if (nbi .gt. nbtmax) then
          write (noutpt,2520) nbtmax,ux48(1:j2)
          write (nttyo,2520) nbtmax,ux48(1:j2)
 2520     format(/' * Error - (EQ3NR/rd3inw) The number of basis',
     $    ' species read',/7x,'from the input file exceeded the',
     $    ' maximum of ',i3,' while',/7x,'trying to read data for',
     $    ' the species ',a,'.',/7x,'Increase the dimensioning',
     $    ' parameter nbtpar.')
          go to 990
        endif
c
        write (noutpt,2530) ux48(1:j2)
 2530   format(1x,'species= ',a)
        uspeci(nbi) = ux48
c
        read (ninpts,2540,err=990) jflgi(nbi),covali(nbi)
 2540   format(10x,i2,12x,e12.5)
        write (noutpt,2550) jflgi(nbi),covali(nbi)
 2550   format(4x,'jflgi= ',i2,4x,'covali= ',1pe12.5)
c
        if (jflgi(nbi).eq.17 .or. jflgi(nbi).eq.18
     $  .or. jflgi(nbi).eq.25) then
          read (ninpts,2570,err=990) ucospi(nbi)
 2570     format(10x,a48)
          j3 = ilnobl(ucospi(nbi))
          write (noutpt,2580) ucospi(nbi)(1:j3)
 2580     format(3x,'ucospi= ',a)
        endif
      enddo
c
  400 nbti = nbi
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Ion exchanger creation. This section reads the directives for
c     creating ion exchange phases and species and their associated
c     intrinsic (including thermodynamic) properties. The data is
c     essentially that which might be read from a supporting data file.
c     Scenario-specific data (e.g., a specified amount of exchanger or
c     a specific composition for an exchanger phase) are not included
c     in this section.
c
c       net    = the number of ion exchangers
c       ugexp  = array of ion exchange phase names
c       mwtges = array of molecular weights of substrates of ion
c                  exchange phases
c       ugexmo = array of strings identifying models for composing
c                  exact ion exchange species and corresponding
c                  reactions; examples include:
c                    'Vanselow' = Vanselow model
c                    'Gapon'    = Gapon model
c                    'Site-mixing' = the general site-mixing model
c       tgexp  = array of reference temperatures (C)  for the
c                  thermodynamic data for the reactions specified
c                  for the exchanger phases
c       jgext  = array of numbers of exchange sites in ion exchange
c                  phases
c       ugexj  = array of names specified for the sites on the various
c                 exchangers
c       cgexj  = array of numbers of moles of exchange sites per mole
c                  of substrate in ion exchange phases
c       zgexj  = array of electrical charges of exchange sites
c                  in ion exchange phases. The total charge number
c                  for a site is the product of this charge and
c                  and the number of moles of the site per mole
c                  of exchanger. The sign of the charge on a site
c                  is generally the opposite of that of the ions
c                  which exchange on it.
c       ngexrt = array of numbers of ion exchange reactions specified
c                  for a given site and exchanger
c       ugexr  = array of strings containing compact representations
c                  of the exchange reactions; e.g., 'Na+ = Ca++' for a
c                  reaction in which Na+ on the exchanger is replaced
c                  by Ca++. One may make specifications such as
c                  'Na+ = ' in which case the ion goes into solution
c                  leaving a bare substrate. All reactions are
c                  normalized to the exchange (or loss) of one
c                  equivalent. The exact form of the reaction is
c                  otherwise dependent on the mixing law specifed in
c                  the element of the ugexmo array for the current
c                  exchanger.
c       xlkgex = array of equilibrium constants or related data for
c                  the reaction specified in the above string
c       uxkgex = array of strings denoting the kind of data in the
c                  corresponding entry of the xlkgex array:
c                    ' '       = log K per equivalent
c                    'LogK/eq' = log K per equivalent
c                    'kcal/eq' = DeltaG0r, kcal, per equivalent
c                    'kJ/eq'   = DeltaG0r, kJ, per equivalent
c       xhfgex = array of enthalpy of reaction data for the reaction
c                  specified in the above string
c       uhfgex = array of strings denoting the kind of data in the
c                  corresponding entry of the xhfgex array:
c                    ' '       = DeltaH0r, kcal, per equivalent
c                    'kcal/eq' = DeltaH0r, kcal, per equivalent
c                    'kJ/eq'   = DeltaH0r, kJ, per equivalent
c       xvfgex = array of volume of reaction data for the reaction
c                  specified in the above string
c       uvfgex = array of strings denoting the kind of data in the
c                  corresponding entry of the xhfgex array:
c                    ' '      = DeltaV0r, cm3 per equivalent
c                    'cm3/eq' = DeltaV0r, cm3 per equivalent
c
      write (noutpt,1590)
 1590 format(' * Ion exchangers')
c
      read (ninpts,1595) ux8
 1595 format(12x,a8)
      call lejust(ux8)
      call locase(ux8)
      qgexsh = ux8(1:2).eq.'.t' .or. ux8(1:1).eq. 't'
      write (noutpt,1597)
 1597 format(5x,'qgexsh= ',l8)
c
      read (ninpts,1600,err=990) net
 1600 format(12x,i3)
      write (noutpt,1610) net
 1610 format(8x,'net= ',i3)
c
      if (net .gt. netmax) then
        write (noutpt,1620) netmax,net
        write (nttyo,1620) netmax,net
 1620   format(/' * Error - (EQ3NR/rd3inw) Have exceeded the maximum',
     $  ' number of ',i3,/7x,'generic ion exchange phases while',
     $  ' reading the data to create',/7x,'such phases. Increase',
     $  ' the dimensioning parameter netpar',/7x,'to at least ',i3,'.')
        go to 990
      endif
c
      do ne = 1,net
        read (ninpts,1630,err=990) ugexp(ne)
 1630   format(12x,a24)
        j2 = ilnobl(ugexp(ne))
        write (noutpt,1640) ugexp(ne)(1:j2)
 1640   format(6x,'ugexp= ',a)
c
        read (ninpts,1650,err=990) mwtges(ne)
 1650   format(12x,e12.5)
        write (noutpt,1660) mwtges(ne)
 1660   format(5x,'mwtges= ',1pe12.5)
c
        read (ninpts,1630,err=990) ugexmo(ne)
        j3 = ilnobl(ugexmo(ne))
        write (noutpt,1670) ugexmo(ne)(1:j3)
 1670   format(5x,'ugexmo= ',a)
c
        read (ninpts,1650,err=990) tgexp(ne)
        write (noutpt,1680) tgexp(ne)
 1680   format(6x,'tgexp= ',1pe12.5)
c
        read (ninpts,1690,err=990) jgext(ne)
 1690   format(12x,i3)
        write (noutpt,1700) jgext(ne)
 1700   format(6x,'jgext= ',i3)
c
        if (jgext(ne) .gt. jetmax) then
          write (noutpt,1710) jetmax,ugexp(ne)(1:j2),jgext(ne)
          write (nttyo,1710) jetmax,ugexp(ne)(1:j2),jgext(ne)
 1710     format(/' * Error - (EQ3NR/rd3inw) Have exceeded the maximum',
     $    ' number of ',i3,/7x,'exchange sites on a generic ion',
     $    ' exchange phase while reading',/7x,'the data to create',
     $    a,'. Increase the',/7x,'dimensioning parameter jetpar to',
     $    ' at least ',i3,'.')
          go to 990
        endif
c
        do je = 1,jgext(ne)
          read (ninpts,1730,err=990) ugexj(je,ne)
 1730     format(12x,a8)
          j3 = ilnobl(ugexj(je,ne))
          write (noutpt,1740) ugexj(je,ne)(1:j3)
 1740     format(6x,'ugexj= ',a)
c
          read (ninpts,1750,err=990) cgexj(je,ne),zgexj(je,ne)
 1750     format(12x,e12.5,12x,e12.5)
          write (noutpt,1760) cgexj(je,ne),zgexj(je,ne)
 1760     format(6x,'cgexj= ',1pe12.5,5x,'zgexj= ',e12.5)
c
          read (ninpts,1690,err=990) ngexrt(je,ne)
          write (noutpt,1810) ngexrt(je,ne)
 1810     format(5x,'ngexrt= ',i3)
c
          if (ngexrt(je,ne) .gt. ietmax) then
            write (noutpt,1820) netmax,ugexj(je,ne)(1:j3),
     $      ugexp(ne)(1:j2),ngexrt(je,ne)
            write (nttyo,1820) netmax,ugexj(je,ne)(1:j3),
     $      ugexp(ne)(1:j2),ngexrt(je,ne)
 1820       format(/' * Error - (EQ3NR/rd3inw) Have exceeded the',
     $      ' maximum number of ',i3,/7x,'reactions for a site',
     $      ' belonging to a generic ion exchange',/7x,'phase while',
     $      ' reading the data for site ',a,' of exchange phase',
     $      /7x,a,'. Increase the dimensioning parameter',
     $      /7x,'ietpar to at least ',i3,'.')
            go to 990
          endif
c
          do n = 1,ngexrt(je,ne)
            read (ninpts,1830,err=990) ugexr(n,je,ne)
 1830       format(12x,a56)
            j2 = ilnobl(ugexr(n,je,ne))
            write (noutpt,1840) ugexr(n,je,ne)(1:j2)
 1840       format(6x,'ugexr= ',a)
            read (ninpts,1850,err=990) xlkgex(n,je,ne),uxkgex(n,je,ne)
 1850       format(12x,e12.5,12x,a8)
            j4 = ilnobl(uxkgex(n,je,ne))
            write (noutpt,1860) xlkgex(n,je,ne),uxkgex(n,je,ne)(1:j4)
 1860       format(5x,'xlkgex= ',1pe12.5,5x,'units= ',a)
            read (ninpts,1850,err=990) xhfgex(n,je,ne),uhfgex(n,je,ne)
            j4 = ilnobl(uhfgex(n,je,ne))
            write (noutpt,1870) xhfgex(n,je,ne),uhfgex(n,je,ne)(1:j4)
 1870       format(5x,'xhfgex= ',1pe12.5,5x,'units= ',a)
            read (ninpts,1850,err=990) xvfgex(n,je,ne),uvfgex(n,je,ne)
            j4 = ilnobl(uvfgex(n,je,ne))
            write (noutpt,1880) xvfgex(n,je,ne),uvfgex(n,je,ne)(1:j4)
 1880       format(5x,'xvfgex= ',1pe12.5,5x,'units= ',a)
          enddo
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Specified generic ion exchanger compositions.
c
c       neti   = the number of generic exchangers for which compositions
c                  are entered
c       ugexpi = array of names of generic ion exchange phases for
c                  which compositions are entered
c       cgexpi = concentration (moles/kg H2O) of generic ion exchange
c                  phases for which compositions are entered
c       jgexti  = array of numbers of sites for the various generic
c                  exchanger phases
c       ugexji = array of names of exchange sites of generic ion
c                  exchange phases
c       ngexti  = array of numbers of exchange species for sites of
c                   the various generic exchanger phasess
c       ugexsi = array of names of species for sites of generic ion
c                  exchange phases
c       egexsi = array of equivalent fractions of ions on exchange sites
c                  in an ion exchange phase
c       xgexsi = array of moles fractions of ions on exchange sites
c                  in an ion exchange phase
c
      write (noutpt,1900)
 1900 format(' * Ion exchanger compositions')
c
      read (ninpts,1910,err=990) neti
 1910 format(12x,i3)
      write (noutpt,1920) neti
 1920 format(7x,'neti= ',i3)
      if (neti .le. 0) go to 250
c
      if (neti .gt. netmax) then
        write (noutpt,1930) netmax,neti
        write (nttyo,1930) netmax,neti
 1930   format(/' * Error - (EQ3NR/rd3inw) Have exceeded the maximum',
     $  ' number of ',i3,/7x,'generic ion exchange phases while',
     $  ' reading the data for concentrations',7x,'and compositions',
     $  ' of such phases. Increase the dimensioning',/7x,'parameter',
     $  ' netpar to at least ',i3,'.')
        go to 990
      endif
c
      do nei = 1,neti
        read (ninpts,1940,err=990) ugexpi(nei)
 1940   format(12x,a24)
        j2 = ilnobl(ugexpi(nei))
        write (noutpt,1950) ugexpi(nei)(1:j2)
 1950   format(5x,'ugexpi= ',a)
c
        read (ninpts,1960,err=990) cgexpi(nei)
 1960   format(12x,e12.5)
        write (noutpt,1970) cgexpi(nei)
 1970   format(5x,'cgexpi= ',1pe12.5)
c
        read (ninpts,1910,err=990) jgexti(nei)
        write (noutpt,1990) jgexti(nei)
 1990   format(5x,'jgexti= ',i3)
c
        if (jgexti(nei) .gt. jetmax) then
          write (noutpt,2000) jetmax,ugexpi(nei)(1:j2),jgexti(nei)
          write (nttyo,2000) jetmax,ugexpi(nei)(1:j2),jgexti(nei)
 2000     format(/' * Error - (EQ3NR/rd3inw) Have exceeded the maximum',
     $    ' number of ',i3,/7x,'exchange sites on a generic ion',
     $    ' exchange phase while reading',/7x,'the data for the',
     $    ' concentration and composition of',/7x,a,'. Increase the',
     $    ' dimensioning parameter jetpar',/7x,'to at least ',i3,'.')
          go to 990
        endif
c
c       Find the corresponding ne index.
c
        do ne = 1,net
          j2 = ilnobl(ugexp(ne))
          j3 = ilnobl(ugexpi(nei))
          if (ugexp(nei)(1:j3) .eq. ugexp(ne)(1:j2)) go to 240
        enddo
c
        write (noutpt,2002) ugexpi(nei)(1:j3)
        write (nttyo,2002) ugexpi(nei)(1:j3)
 2002   format(/' * Error - (EQ3NR/rd3inw) Data are present on the',
     $  ' input file for the',/7x,'generic ion exchange phase "',a,
     $  ', but this phase',/7x,"hasn't been previously defined.",
     $  " Can't determine the requisite model,",/7x,"therefore can't",
     $  ' finish reading the current data for this phase.')
        stop
c
  240   j2 = ilnobl(ugexmo(ne))
        if (ugexmo(ne)(1:j2) .eq. 'Gapon' .or.
     $    ugexmo(ne)(1:6) .eq. 'Gapon-' .or.
     $    ugexmo(ne)(1:j2) .eq. 'Vanselow' .or.
     $    ugexmo(ne)(1:9) .eq. 'Vanselow-') then
c
c         Input composition is described in terms of equivalent
c         fractions on the sites.
c
          qgexef = .true.
        else
c
c         Input composition is described in terms of mole fractions
c         on the sites.
c
          qgexef = .false.
        endif
c
        do jei = 1,jgexti(nei)
          read (ninpts,2010,err=990) ugexji(jei,nei)
 2010     format(12x,a8)
          j3 = ilnobl(ugexji(jei,nei))
          write (noutpt,2020) ugexji(jei,nei)(1:j3)
 2020     format(5x,'ugexji= ',a)
c
          read (ninpts,1910,err=990) ngexti(jei,nei)
          write (noutpt,2030) ngexti(jei,nei)
 2030     format(5x,'ngexti= ',i3)
c
          if (ngexti(jei,nei) .gt. ietmax) then
            write (noutpt,2040) ietmax,ugexji(jei,nei),ugexpi(nei),
     $      ngexti(jei,nei)
            write (nttyo,2040) ietmax,ugexji(jei,nei),ugexpi(nei),
     $      ngexti(jei,nei)
 2040       format(/' * Error - (EQ3NR/rd3inw) Have exceeded the',
     $      ' maximum number of ',i3,/7x,'species on an exchange',
     $      ' site of a generic ion exchange phase',/7x,'while reading',
     $      /7x,'the data for the composition of site ',a,' of',
     $      /7x,'exchange phase ',a,'. Increase the dimensioning',
     $      ' parameter ietpar to at least ',i3,'.')
            go to 990
          endif
c
          if (qgexef) then
c
c           The data are in equivalent fractions.
c
            do iei = 1,ngexti(jei,nei)
              read (ninpts,2050,err=990) ugexsi(iei,jei,nei),
     $        egexsi(iei,jei,nei)
 2050         format(12x,a24,10x,e12.5)
              write (noutpt,2052) ugexsi(iei,jei,nei),
     $        egexsi(iei,jei,nei)
 2052         format(5x,'ugexsi= ',a24,2x,'egexsi= ',1pe12.5)
            enddo
          else
c
c           The data are in mole fractions.
c
            do iei = 1,ngexti(jei,nei)
              read (ninpts,2050,err=990) ugexsi(iei,jei,nei),
     $        xgexsi(iei,jei,nei)
              write (noutpt,2060) ugexsi(iei,jei,nei),
     $        xgexsi(iei,jei,nei)
 2060         format(5x,'ugexsi= ',a24,2x,'xgexsi= ',1pe12.5)
            enddo
          endif
        enddo
      enddo
c
  250 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Specified solid solution compositions.
c
c       nxti   = the number of solid solutions for which compositions
c                  are entered
c       usoli  = array of names of solid solutions for which
c                  compositions are entered
c       ncmpri = range pointer array for solid solutions for which
c                  compositions are entered; defines the range in
c                  the umemi and xbari arrays which corresponds to the
c                  end-member components of a given solid solution
c       nssoti = number of end-member components in a given solid
c                  solution
c       umemi  = array of names of end-member components of solid
c                  solutions for which compositions are entered
c       xbari  = array of moles fractions of end-member components
c                  of solid solutions for which compositions are entered
c
      write (noutpt,2100)
 2100 format(' * Solid solution compositions')
c
      read (ninpts,2102,err=990) nxti
 2102 format(12x,i3)
      write (noutpt,2104) nxti
 2104 format(7x,'nxti= ',i3)
c
      if (nxti .gt. nxtimx) then
        write (noutpt,2107) nxtimx,nxti
        write (nttyo,2107) nxtimx,nxti
 2107   format(/' * Error - (EQ3NR/rd3inw) Have exceeded the maximum',
     $  /7x,'number of ',i5,' solid solutions for which',
     $  ' compositions',/7x,'may be specified on the input file.',
     $  ' Increase the dimensioning',/7x,'parameter nxtipa to at',
     $  ' least ',i3,'.')
        go to 990
      endif
c
      nxic = 0
c
      do nxi = 1,nxti
        read (ninpts,2120,err=990) usoli(nxi)
 2120   format(12x,a24)
        j2 = ilnobl(usoli(nxi))
        write (noutpt,2130) usoli(nxi)(1:j2)
 2130   format(6x,'usoli= ',a24)
c
        read (ninpts,2140,err=990) nssoti
 2140   format(12x,i3)
        write (noutpt,2142) nssoti
 2142   format(5x,'nssoti= ',i3)
c
        ncmpri(1,nxi) = nxic + 1
        ncmpri(2,nxi) = nxic + nssoti
        nr1 = ncmpri(1,nxi)
        nr2 = ncmpri(2,nxi)
c
        if (nr2 .gt. nxicmx) then
          write (noutpt,2150) nxicmx,usoli(nxi)(1:j2),nr2
          write (nttyo,2150) nxicmx,usoli(nxi)(1:j2),nr2
 2150     format(/' * Error - (EQ3NR/rd3inw) Have exceeded the',
     $    ' maximum number',/7x,'of ',i5,' solid solution species',
     $    ' for which mole fractions',/7x,'are specified on the',
     $    ' input file. This occurred while reading ',/7x,'data',
     $    ' for the solid solution ',a,'. Increase',/7x,'the',
     $    ' dimensioning parameter nxicpa to at least ',i3,'.')
          go to 990
        endif
c
        do nxic = nr1,nr2
          read (ninpts,2160,err=990) umemi(nxic),xbari(nxic)
 2160     format(12x,a24,12x,e12.5)
          write (noutpt,2170) umemi(nxic),xbari(nxic)
 2170     format(6x,'umemi= ',a24,5x,'xbari= ',1pe12.5)
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Number of nxmod options.
c
c       nxmod = the number of suppressed/altered species/reactions
c
      write (noutpt,2330)
 2330 format(' * Alter/suppress options')
c
      read (ninpts,2340,err=990) nxmod
 2340 format(12x,i3)
      write (noutpt,2350) nxmod
 2350 format(6x,'nxmod= ',i3)
c
      if (nxmod .gt. nxmdmx) then
        write (noutpt,2360) nxmdmx
        write (nttyo,2360) nxmdmx
 2360   format(/' * Error - (EQ3NR/rd3inw) Have too many nxmod',
     $  /7x,'alter/suppress options. The code is only dimensioned',
     $  /7x,'for ',i3,' such options. Reduce the number of such',
     $  /7x,'options or increase the dimensioning parameter nxmdpa.')
        go to 990
      endif
c
c     Nxmod options.
c
c       uxmod = the name (all 48 letters) of the species. If the
c                 phase part is not given, the option is applied to
c                 every species for which the species part of its name
c                 generates a match.
c       kxmod  = alter/suppress code
c         -1 = the species is suppressed
c          0 = the log K is replaced by xlkmod
c          1 = the log K is augmented by xlkmod
c          2 = same as kxmod=1, but xlkmod is input in units of kcal
c                per mole of the associated species
c       xlkmod = log K value alteration function as defined above
c
      do n = 1,nxmod
        read (ninpts,2370,err=990) uxmod(n),kxmod(n),xlkmod(n)
 2370   format(12x,a48,/12x,i2,22x,e12.5)
        j2 = ilnobl(uxmod(n))
        write (noutpt,2380) uxmod(n)(1:j2),kxmod(n),xlkmod(n)
 2380   format(4x,'species= ',a,/5x,'option= ',i2,14x,
     $  'xlkmod= ',1pe12.5)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iopt model option switches.
c     Note: iopt(1) = iopt1, etc.
c
c       iopt(1) - iopt(3) = Used only by EQ6
c
c       iopt(4) = Solid solution products:
c          0 = Solid solutions are ignored
c          1 = Solid solutions are permitted
c          2 = Solid solutions are permitted and compositions to
c                use in computing saturation indices are read
c
c       iopt(5) - iopt(7) = Used only by EQ6
c
c       iopt(8) = Not used
c
c       iopt(9) - iopt(10) = Used only by EQ6
c
c       iopt(11) = Auto basis switching, in pre-N-R optimization:
c          0 = Off
c          1 = On
c
c       iopt(12) - iopt(13) Used only by EQ6
c
c       iopt(14) = Not used
c
c       iopt(15) - iopt(16) = Used only by EQ6
c
c       iopt(17) = pickup file options:
c         -1 = Don't write a pickup file
c          0 = Write a pickup file at the end of the run
c
c       iopt(18) = Used only by EQ6
c
c       iopt(19) = Advanced EQ3NR pickup file options:
c          0 = Write a normal EQ3NR pickup file
c          1 = Write an EQ6 input file with Quartz dissolving, using
c                a relative rate law
c          2 = Write an EQ6 input file with Albite dissolving, using
c                a TST rate law
c          3 = Write an EQ6 input file to mix two fluids appearing
c                on a single EQ3NR input file; the first fluid on that
c                file is set up as a special reactant called "Fluid 2"
c                that will be added to the second fluid on that
c                file, which in effect becomes "Fluid 1". Set iopt(19)
c                to 3 only for the second fluid.
c
c       iopt(20) = Used only by EQ6
c
      write (noutpt,1365)
 1365 format(' * Iopt, iopg, iopr, and iodb options')
c
      write (noutpt,1370)
 1370 format(' *',15x,'1    2    3    4    5    6    7    8    9   10')
      read (ninpts,1380,err=990) (iopt(i), i = 1,20)
 1380 format(12x,10i5)
      write (noutpt,1390) (iopt(i), i = 1,20)
 1390 format(3x,'iopt1-10= ',10i5,/2x,'iopt11-20= ',10i5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iopg activity coefficient option switches.
c     Note: iopg(1) = iopg1, etc.
c
c       iopg(1) = Model for aqueous species activity coefficients:
c         -1 = The  Davies equation
c          0 = The B-dot equation
c          1 = Pitzer's  equations
c          2 = HC + DH equations
c
c       iopg(2) = Rescaling of aqueous ionic activity coefficients for
c                   consistency with a desired pH scale:
c         -1 = "Internal" pH scale; no scaling (e.g., the "Davies"
c                scale if iopg(1) = -1, the "B-dot" scale if
c                iopg(1) = 0, or the Pitzer scale if iopg(1) = 1)
c          0 = The NBS pH scale (log gamma(Cl-) is defined by the
c                Bates-Guggenheim equation)
c          1 = The Mesmer pH scale (log gamma(H+) = 0)
c
c       iopg(3) = iopg(20) = Not used
c
      read (ninpts,1410,err=990) (iopg(i), i = 1,20)
 1410 format(12x,10i5)
      write (noutpt,1420) (iopg(i), i = 1,20)
 1420 format(3x,'iopg1-10= ',10i5,/2x,'iopg11-20= ',10i5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iopr print option switches.
c     Note: iopr(1) = iopr1, etc.
c
c       iopr(1) = List the names of all species read from the
c                   supporting data file:
c          0 = Don't list
c          1 = List
c
c       iopr(2) = List all reactions:
c          0 = Don't list
c          1 = List all reactions (this can be quite lengthy)
c          2 = Also print the log K values
c          3 = Also print the coefficients of the interpolating
c                polynomials
c
c       iopr(3) = List the hard core diameters of the aqueous species:
c          0 = Don't list
c          1 = List
c
c       iopr(4) = Print a table at each print point of the
c                   concentrations, activities, and activity
c                   coefficients of the aqueous species:
c         -3  = Omit species with molalities < 1.e-8
c         -2 =  Omit species with molalities < 1.e-12
c         -1 =  Omit species with molalities < 1.e-20
c          0 =  Omit species with molalities < 1.e-100
c          1 =  Include all species
c
c       iopr(5) = Print a table at each print point of the cation/H+
c                   activity ratios, anion-H+ activity products, and
c                   neutral species activities:
c          0 = Don't print
c          1 = Print cation/H+ activity ratios only
c          2 = Print cation/H+ activity ratios and anion-H+ activity
c                products only
c          3 = Print cation/H+ activity ratios, anion-H+ activity
c                products, and neutral species activities
c
c       iopr(6) = At each print point, print a table of the percentage
c                    contributions for each aqueous mass balance total:
c         -1 = Don't print any tables
c          0 = Print tables including 99% of all contributing species
c          1 = Print tables including all contributing species
c
c       iopr(7) = Print a table at each print point of the saturation
c                   indices and affinities of the various non-aqueous
c                   phases:
c         -1 = Don't print
c          0 = Print for those phases not undersaturated by
c                more than 10 kcal
c          1 = Print for all phases
c
c       iopr(8) = Print a table at each print point of the fugacities
c                   of the gas species:
c         -1 = Don't print
c          0 = Print
c          1 = Print
c
c       iopr(9) = Print a table at each print of the mean molal
c                   activity coefficients:
c         -1 = Don't print
c          0 = Don't print
c          1 = Print
c
c       iopr(10) = Print a tabulation at the start of running the
c                    current problem of the Pitzer interaction
c                    coefficients:
c          0 = Don't print
c          1 = Print a summary of the names of the species present and
c                the number of Pitzer interaction coefficients
c          2 = Print a summary of the names of the species present and
c                the number of Pitzer interaction coefficients
c
c       iopr(11) - iopr(16) = Not used
c
c       iopr(17) = pickup file format:
c          0 = Use the same format ("D" or "W") as the input file
c          1 = Use "W" format
c          2 = Use "D" format
c
c       iopr(18) - iopr(20) = Not used
c
      read (ninpts,1430,err=990) (iopr(i), i = 1,20)
 1430 format(12x,10i5)
      write (noutpt,1440) (iopr(i), i = 1,20)
 1440 format(3x,'iopr1-10= ',10i5,/2x,'iopr11-20= ',10i5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iodb debugging print option switches.
c     Note: iodb(1) = iodb1, etc.
c
c       iodb(1) = General diagnostic messages:
c          0 = Don't print
c          1 = Print Level 1 diagnostic messages
c          2 = Print Level 1 and Level 2 diagnostic messages
c
c       iodb(2)  = Used only by EQ6
c
c       iodb(3) = Pre-Newton-Raphson optimization:
c          0 = Don't print
c          1 = Print summary information
c          2 = Print detailed information
c          3 = Print more detailed information
c          4 = Also print changes to activity coefficients
c
c       iodb(4) = Newton-Raphson iteration:
c          0 = Don't print
c          1 = Print summary information
c          2 = Print detailed information including the residual and
c                correction vectors
c          3 = Also print the Jacobian matrix
c          4 = Also print changes to activity coefficients
c
c       iodb(5)  = Used only by EQ6
c
c       iodb(6) = Hypothetical affinity calculations (for solid
c                   solutions, etc.):
c          0 = Don't print
c          1 = Print summary information
c          2 = Print detailed information
c
c       iodb(7)  = Used only by EQ6
c
c       iodb(8) -iodb(20) = Not used
c
      read (ninpts,1470,err=990) (iodb(i), i = 1,20)
 1470 format(12x,10i5)
      write (noutpt,1480) (iodb(i), i = 1,20)
 1480 format(3x,'iodb1-10= ',10i5,/2x,'iodb11-20= ',10i5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Numerical parameters.
c
c       tolbt = convergence tolerance on residual magnitude
c                 (in Newton-Raphson iteration)
c       toldl = convergence tolerance on correction magnitude
c                 (in Newton-Raphson iteration)
c       itermx = maximum number of Newton-Raphson iterations
c
      write (noutpt,1270)
 1270 format(' * Numerical parameters')
c
      read (ninpts,1280,err=990) tolbt,toldl
 1280 format(12x,e12.5,12x,e12.5)
      write (noutpt,1290) tolbt,toldl
 1290 format(6x,'tolbt= ',1pe12.5,5x,'toldl= ',e12.5)
c
      read (ninpts,1300,err=990) itermx
 1300 format(12x,i3)
      write (noutpt,1310) itermx
 1310 format(5x,'itermx= ',i3)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Ordinary basis switches.
c
c       nobswt = the number of ordinary basis switches
c       uobsw(1,n) = the species to be switched from the basis
c                      set in the n-th switch
c       uobsw(2,n) = the species to be switched into the basis
c                      set in the n-th switch
c
      write (noutpt,2700)
 2700 format(' * Ordinary basis switches')
c
      read (ninpts,1110,err=990) nobswt
      write (noutpt,2710) nobswt
 2710 format(5x,'nobswt= ',i3)
c
      do n = 1,nobswt
        read (ninpts,1130,err=990) uobsw(1,n)
        j2 = ilnobl(uobsw(1,n))
        write (noutpt,1140) uobsw(1,n)(1:j2)
        read (ninpts,1150,err=990) uobsw(2,n)
        j2 = ilnobl(uobsw(2,n))
        write (noutpt,1160) uobsw(2,n)(1:j2)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Saturation flag tolerance.
c
c       tolspf = saturation flag tolerance for minerals and
c                  various other phases, in kcal. This parameter
c                  controls the printing of the "SATD" and "SSATD"
c                  flags in the saturation state tables. It does not
c                  have any other effect in EQ3NR.
c
      write (noutpt,1292)
 1292 format(' * Saturation flag tolerance')
c
      read (ninpts,1294,err=990) tolspf
 1294 format(12x,e12.5)
      write (noutpt,1296) tolspf
 1296 format(5x,'tolspf= ',1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Scale factor for the mass of aqueous solution to write
c     on the pickup file.
c
c       scamas = scale factor (the default is 1.0, which corresponds
c                  to a mass of aqueous solution containing 1 kg
c                  of solvent water)
c
      write (noutpt,1315)
 1315 format(' * Aqueous phase scale factor')
c
      read (ninpts,1320,err=990) scamas
 1320 format(12x,e12.5)
      write (noutpt,1330) scamas
 1330 format(5x,'scamas= ',1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,3000) nprob
      write (nttyo,3000) nprob
 3000 format(/'   Done reading problem ',i3,'.',/)
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  990 qrderr = .true.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
