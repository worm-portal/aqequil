      subroutine intexi(al10,axhfs,axlks,axvfs,cegexs,cess,cdrs,
     $ cgexj,cpgexs,egexjf,iern1,iern2,ietmax,jern1,jern2,jetmax,
     $ jflag,jgext,jpflag,jsflag,jsitex,kern1,kern2,ketmax,kgexsa,
     $ mwtges,mwtsp,narn1,narn2,narxmx,narxt,nbasp,nbt,nbtmax,ncmpr,
     $ ndrs,ndrsmx,ndrsr,ness,nessmx,nessr,netmax,net,nern1,nern2,
     $ ngexro,ngexrt,ngexsa,ngexso,ngext,noutpt,nphasx,npt,nptmax,
     $ nst,nstmax,ntprt,ntprmx,nttyo,nvetmx,rconst,tgexp,ugexj,
     $ ugexmo,ugexmv,ugexp,ugexr,ugexs,ugexsr,uhfgex,uphase,uspec,
     $ uvfgex,uxkgex,xhfgex,xlkgex,xvfgex,zchar,zgexj)
c
c     This subroutine inteprets data read from the input file for the
c     generic ion exchange model. It composes the necessary generic
c     exchange phases and species and sets up the associated reactions
c     and thermodynamic properties.
c
c     "Ion exchange" may be strictly or loosely interpreted here,
c     depending on the models chosen. In the case of the Gapon and
c     Vanselow models, the interpretation is strict. For example, a
c     site must be negatively or positively charged, only ions of
c     the opposite sign may occupy it, and the formal exchange
c     capacity is exactly 100% utilized, save for a negligible
c     trace of unfilled positions. The end-member exchanger species
c     are electrically neutral, save for the trace unfilled site
c     species. In contrast, the "Site-mixing" model allows more
c     flexibility, including charged exchanger end-members, and
c     overloading or underloading of the formal exchange capacities.
c
c     Each exchanger species is specific to a given site (type of
c     site), and these species are organized in the species list such
c     that those belonging to a given site are listed contiguously.
c     An exchanger with only one site is treated in the same manner
c     as a phase composed of constituent species. In the case of an
c     exchanger with more than one site, each site is treated as a
c     kind of subphase, composed of its own constituent species.
c
c     To illustrate, an exchange phase with substrate Z and two exchange
c     sites is represented by:
c
c       (Na+,K+)3(Fe+++,Al+++)2(Z)
c
c     Here a pair of parentheses represents a site. Z represents the
c     substrate. The site it occupies is not an exchange site. There is
c     always one mole of substrate occupying this site. The first
c     exchange site, here called site A, contains a mixture of Na+
c     and K+ ions. There are 3 moles of this site per mole of substrate.
c     This site has an intrinsic charge (charge when the site is not
c     occupied) of -1. Thus, 3 moles of monovalent cations per mole of
c     substrate are required to attain electrical neutrality on this
c     site. The second exchange site, here called site B, contains a
c     mixture of Fe+++ and Al+++ ions. There are 2 moles of this site
c     per mole of substrate, and it has an intrinsic charge of -3.
c
c     The species for this example may be written as:
c
c       site A:
c
c         (Na+)3()2()
c         (K+)3()2()
c
c       site B:
c
c         ()3(Fe+++)2()
c         ()3(Al+++)2()
c
c       substrate site:
c
c         ()3()2(Z)
c
c     Here an empty pair of parentheses denotes an empty site.
c     However, an empty site is a definite part of the species.
c     An empty site makes no contribution to the molecular weight of
c     of such a species. Thus, the molecular weight of (Na+)3()2()
c     is 3 times the atomic weight of Na. The mass in grams of the
c     complete exchanger phase is the sum of the masses of all of
c     the constitutent species. The number of moles of the exchanger
c     phase is equal to the number of moles of the substrate.
c
c     One could conceivably define an exchanger with a non-unit
c     stoichiometry for the substrate component; e.g., something
c     like:
c
c       (Na+,K+)3(Fe+++,Al+++)2(Z')2
c
c     where Z' is a redefined substrate for the previous example.
c     However, there do not seem to be any advantages to introducing
c     or even allowing such a complexity. One can always define a
c     stoichiometric with unit substrate, in this example by taking
c     Z = 2*Z'.
c
c     The activity of such a species in the ideal case (implicitly in
c     the site mixing sense) is the mole fraction on the site for which
c     the species is defined, raised to the power equal to the
c     stoichiometric factor for the site. In the case of (Na+)3()2(),
c     the relation is:
c
c       a{(Na+)3()2()} = x{(Na+)3()2()}**3.
c
c     The origin of this relationship is as follows. Write a pseudo-
c     reaction in which the (Na+)3()2() species goes to a form of
c     itself, but with the stoichiometry redefined so that there is
c     one mole of the site of interest per mole of the exchanger
c     species:
c
c       (Na+)3()2() = 3(Na+)()[2/3]()[1/3]
c
c     Because this is a pseudo-reaction, reflecting only a change
c     in components, there is no change in the Gibbs energy, enthalpy,
c     or entropy. Alternatively, the equilibrium constant is unity.
c     Thus,
c
c       a{(Na+)3()2()} = a{(Na+)()[2/3]()[1/3]}**3
c
c     In the ideal case, it is clear that:
c
c       a{(Na+)()[2/3]()[1/3]} = x{(Na+)()[2/3]()[1/3]}
c
c     Thus,
c
c       a{(Na+)3()2()} = x{(Na+)()[2/3]()[1/3]}**3
c
c     The mole fractions of the Na+ species are not affected by the
c     stoichiometric arbitrariness in how the species are defined.
c     Hence:
c
c       x{(Na+)3()2()} = x{(Na+)()[2/3]()[1/3]}
c
c     Subsitution of this into the equaiton immediately above then
c     yields:
c
c       a{(Na+)3()2()} = x{(Na+)3()2()}**3.
c
c     This explains the origin of the site number appearing as an
c     exponent in the calculation of the activity of a component.
c
c     The above kinds of species are not end-members. End-members for
c     the above exchanger phase would be:
c
c       (Na+)3(Fe+++)2(Z)
c       (Na+)3(Al+++)2(Z)
c       (K+)3(Fe+++)2(Z)
c       (K+)3(Al+++)2(Z)
c
c     Note that for example:
c
c       a((Na+)3(Fe+++)2(Z)) = a((Na+)3()2())**3 a(()3(Fe+++)2())**2
c
c     Note that a factor of a(Z)**1 is missing from the right hand
c     side. This is because this factor has a fixed value of 1 by
c     definition. This would be true even if we allowed for the
c     formal possibility of an arbitrary number of moles of substrate
c     in the exchanger phase, because there is no mixing on the
c     substrate site. Thus, a(Z) = 1.
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
c       ugexj  = array of exchanger site names
c       ugexmo = array of strings identifying exchange models, which
c                  are used for composing ion exchange species and
c                  corresponding reactions; examples include:
c                    'Gapon'    = Gapon model
c                    'Vanselow' = Vanselow model
c                    'Site-mixing' = the general site-mixing model
c       ugexmv = array of valid strings identifying exchange models
c       ugexp  = array of exchanger phase names
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
c       ugexs  = array of exchange ions names
c       ugexsr = array of exchange ions names extracted from the
c                  ugexr array
c
c     Principal input/output:
c
c       uphase = array of phase names (extended to include the names
c                  of generic exchange phases)
c       uspec  = array of species names (extended to include the names
c                  of species belonging to generic exchange phases)
c
c     Principal output:
c
c       cegexs = array of coefficients giving the number of equivalents
c                  per mole of each exchanger species.
c       cpgexs = array of coefficients giving the number of moles of
c                  exchanger substrate per mole of each exchanger
c                  species.
c       egexjf = array of formal exchange capacities of the sites of
c                  the exchangers, defined in terms of equivalents per
c                  mole of exchanger substrate. The formal exchange
c                  capacity of an exchanger phase is the sum of these
c                  for all its exchange sites.
c       jern1  = array giving the start of the range in the species
c                  list corresponding to species in the je-th site
c                  of the ne-th exchanger.
c       jern2  = array giving the end of the range in the species
c                  list corresponding to species in the je-th site
c                  of the ne-th exchanger.
c       jgext  = array giving the number of exchange sites in each of
c                  the ion exchange phases
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ietmax,jetmax,ketmax,narxmx,nbtmax,ndrsmx,nessmx,netmax,
     $ nptmax,nstmax,ntprmx,nvetmx
c
      integer noutpt,nttyo
c
      integer jern1(jetmax,netmax),jern2(jetmax,netmax),jflag(nstmax),
     $ jgext(netmax),jpflag(nptmax),jsflag(nstmax),jsitex(nstmax),
     $ kern1(netmax),kern2(netmax),kgexsa(ketmax,netmax),narxt(ntprmx),
     $ nbasp(nbtmax),ncmpr(2,nptmax),ndrs(ndrsmx),ndrsr(2,nstmax),
     $ ness(nessmx),nessr(2,nstmax),ngexro(ietmax,jetmax,netmax),
     $ ngexrt(jetmax,netmax),ngext(jetmax,netmax),
     $ ngexsa(ietmax,jetmax,netmax),ngexso(ietmax,jetmax,netmax),
     $ nphasx(nstmax)
c
      integer iern1,iern2,narn1,narn2,nbt,net,nern1,nern2,npt,nst,ntprt
c
      character*56 ugexr(ietmax,jetmax,netmax)
      character*48 uspec(nstmax)
      character*24 ugexmo(netmax),ugexmv(nvetmx),ugexp(netmax),
     $ ugexs(ietmax,jetmax,netmax),ugexsr(2,ietmax,jetmax,netmax),
     $ uphase(nptmax)
      character*8 ugexj(jetmax,netmax),uhfgex(ietmax,jetmax,netmax),
     $ uvfgex(ietmax,jetmax,netmax),uxkgex(ietmax,jetmax,netmax)
c
      real*8 axhfs(narxmx,ntprmx,nstmax),axlks(narxmx,ntprmx,nstmax),
     $ axvfs(narxmx,ntprmx,nstmax),cegexs(ietmax,jetmax,netmax),
     $ cess(nessmx),cdrs(ndrsmx),cgexj(jetmax,netmax),
     $ cpgexs(ietmax,jetmax,netmax),egexjf(jetmax,netmax),
     $ mwtges(netmax),mwtsp(nstmax),tgexp(netmax),xhfgex(ietmax,
     $ jetmax,netmax),xlkgex(ietmax,jetmax,netmax),
     $ xvfgex(ietmax,jetmax,netmax),zchar(nstmax),zgexj(jetmax,netmax)
c
      real*8 al10,rconst
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ie,iee,iej,ieo,itot,j,je,jee,jj,jj2,jj3,j2,j3,j4,j5,j6,
     $ k,ke,kee,ke1,n,nb,nd2,ne,nee,nerr,nn,np,nrf1,nrf2,nr1,nr2,ns,
     $ nsba,nse,nsl,nspect,nss,nsse,nve,nvet
c
      integer ilnobl,nbasis
c
      logical qmoerr,qsperr
c
      character*56 ustr56
      character*24 ugexpd,ustr1,ustr1e,ustr2
      character*8 ugexjd
c
      real afhx,afvx,afxb,afx0,arcnst,bfx,cfactr,cgx,cfxi,cfxz,cx,
     $ efx,tfxb,tfx0,xx,zprod,zpx,zsi,zxj,zxjt
c
c-----------------------------------------------------------------------
c
c     Initialize some constants.
c
      nerr = 0
      arcnst = -al10*0.001*rconst
      tfx0 = 273.15
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Count the number of valid exchange model options.
c
      nvet = 0
      do nve = 1,nvetmx
        if (ugexmv(nve)(1:6) .ne. 'ERROR ') nvet = nvet + 1
      enddo
c
      if (net.gt.0 .and. nvet.le.0) then
        write (noutpt,1000)
        write (nttyo,1000)
 1000   format (/' * Error - (EQLIB/intexi) Programming error trap:',
     $  /7x,'There are no programmed valid strings representing',
     $  /7x,'exchange models for the generic ion phases. Check the',
     $  /7x,'data statement which initializes the ugexmv array. This',
     $  /7x,'is located in the EQLIB INCLUDE file eqlo8d.h. If the',
     $  /7x,'INCLUDE files have been "pre-stuffed" in your source',
     $  /7x,'code, this data statement may be found in the main',
     $  /7x,'programs for EQ3NR and EQ6.')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Pass 1. Loop on exchanger phases. Check the input. Provide
c     defaults as necessary. build a list of exchanger species, and
c     map exchange ions to their aqueous species counterparts.
c     The exchange phases and species will be created (in the
c     sense of being added to the general lists of phases and species)
c     below in pass 2.
c
      do ne = 1,net
        call lejust(ugexp(ne))
        j4 = ilnobl(ugexp(ne))
c
        if (j4 .le. 0) then
c
c         Have a blank exchanger name. Assign a default name.
c
          call adgexp(ne,noutpt,nttyo,ugexpd)
          ugexp(ne) = ugexpd
          j4 = ilnobl(ugexp(ne))
        endif
c
c       Make sure that the phase name is unique among all exchangers.
c
        do nee = 1,ne - 1
          if (ugexp(ne)(1:24) .eq. ugexp(nee)(1:24)) then
            write (noutpt,1010) ugexp(ne)(1:j4),ne,nee
            write (nttyo,1010) ugexp(ne)(1:j4),ne,nee
 1010       format (/' * Error - (EQLIB/intexi) "',a,'" is the name',
     $      ' specified',/7x,'for exchanger phases ',i3,' and ',i3,
     $      '. Exchanger phase names,',/7x,'including any assigned',
     $      ' defaults, must be unique.')
           nerr = nerr + 1
          endif
        enddo
c
c       Check the prescribed exchange model.
c
        call lejust(ugexmo(ne))
c
c       Provide a default.
c
        if (ugexmo(ne)(1:3) .eq. '   ') then
          ugexmo(ne) = ugexmv(1)
          write (noutpt,1012) ugexp(ne)(1:j4),ugexmo(ne)(1:j2)
          write (nttyo,1012) ugexp(ne)(1:j4),ugexmo(ne)(1:j2)
 1012     format (/' * Warning - (EQLIB/intexi) No exchange model',
     $    ' was species for',/7x,'the exchanger phase ',a,'. The ',a,
     $    ' model has been',/7x,'assigned as a default.')
        endif
c
        j5 = ilnobl(ugexmo(ne))
        do nve = 1,nvetmx
          j6 = ilnobl(ugexmv(nve))
          if (ugexmo(ne)(1:j5) .eq. ugexmv(nve)(1:j6)) go to 130
        enddo
c
        write (noutpt,1020) ugexmo(ne)(1:j2),ugexp(ne)(1:j4)
        write (nttyo,1020) ugexmo(ne)(1:j2),ugexp(ne)(1:j4)
 1020   format (/" * Error - (EQLIB/intexi) Don't recognize the",
     $  ' exchange model type',/7x,'"',a,'" specified for the',
     $  ' generic exchanger phase ',a,'.',/7x,'Valid choices',
     $  ' include the following:',/)
c
        do nve = 1,nvetmx
          if (ugexmv(nve)(1:6) .ne. 'ERROR ') then
            j6 = ilnobl(ugexmv(nve))
            write (noutpt,1022) ugexmv(nve)(1:j6)
            write (nttyo,1022) ugexmv(nve)(1:j6)
 1022       format(9x,a)
          endif
        enddo
c
        j6 = ilnobl(ugexmv(1))
        write (noutpt,1024) ugexmv(1)(1:j6)
        write (nttyo,1024) ugexmv(1)(1:j6)
 1024   format(/7x,'A blank input defaults to "',a,'".')
        nerr = nerr + 1
c
  130   if (mwtges(ne) .le. 0) then
          mwtges(ne) = 100.
          write (noutpt,1030) ugexp(ne)(1:j4)
          write (nttyo,1030) ugexp(ne)(1:j4)
 1030     format (/" * Warning - (EQLIB/intexi) Don't have a valid",
     $    ' input for the',/7x,'molecular weight of the substrate for',
     $    ' exchange phase'/7x,a,'. Assigning a default value of',
     $    ' 100 grams/mol.')
        endif
c
        if (tgexp(ne) .eq. 0.) then
          tgexp(ne) = 25.0
          write (noutpt,1040) ugexp(ne)(1:j4)
          write (nttyo,1040) ugexp(ne)(1:j4)
 1040     format (/" * Warning - (EQLIB/intexi) Don't have a valid",
     $    ' input for the',/7x,'temperature reference required for',
     $    ' the thermodynamic data for',/7x,'exchange phase ',a,
     $    '. Assigning a default value of 25C.')
        endif
        cfactr = 1./(arcnst*(tgexp(ne) + 273.15))
c
c       Loop on sites.
c
        do je = 1,jgext(ne)
          call lejust(ugexj(je,ne))
          j3 = ilnobl(ugexj(je,ne))
c
c         Calculate the formally declared exchange capacity (in
c         equivalents per mole of exchanger substrate) of the
c         current site. The sign of this is opposite to that of the
c         charge on the site itself.
c
          egexjf(je,ne) = -cgexj(je,ne)*zgexj(je,ne)
c
c         Check the name of the current site.
c
          if (j3 .eq. 0) then
c
c           No name was given. Assign one (e.g., "S(1)" to site 1,
c           "S(2)" to site 2).
c
            call adgexj(je,noutpt,nttyo,ugexjd)
            ugexj(je,ne) = ugexjd
            j3 = ilnobl(ugexj(je,ne))
          endif
c
c         Make sure that the site name is unique among all sites for
c         the current exchanger phase.
c
          do jee = 1,je - 1
            if (ugexj(je,ne)(1:8) .eq. ugexj(jee,nee)(1:8)) then
              write (noutpt,1060) ugexj(je,ne)(1:j3),je,jee,
     $        ugexp(ne)(1:j4)
              write (nttyo,1060) ugexj(je,ne)(1:j3),je,jee,
     $        ugexp(ne)(1:j4)
 1060         format (/' * Error - (EQLIB/intexi) "',a,'" is the name',
     $        ' of',/7x,'sites ',i3,' and ',i3,' of exchanger phase ',
     $        a,'.',/7x,'Site names, including any assigned',
     $        ' defaults, must be unique',/7x,'for a given',
     $        ' exchanger phase.')
              nerr = nerr + 1
            endif
          enddo
c
c         Check the electrical charge of the site against the specified
c         exchange model.
c
          if (ugexmo(ne)(1:5).eq.'Gapon' .or.
     $      ugexmo(ne)(1:6).eq.'Gapon-' .or.
     $      ugexmo(ne)(1:8).eq.'Vanselow' .or.
     $      ugexmo(ne)(1:9).eq.'Vanselow-') then
            if (zgexj(je,ne) .eq. 0.) then
              write (noutpt,1070) ugexj(je,ne)(1:j3),ugexp(ne)(1:j4),
     $        ugexmo(ne)(1:j5)
              write (nttyo,1070) ugexj(je,ne)(1:j3),ugexp(ne)(1:j4),
     $        ugexmo(ne)(1:j5)
 1070         format (/' * Error - (EQLIB/intexi) The electrical',
     $        ' charge of site ',a,''/7x,'of exchanger phase ',
     $        a,' is zero. This is not valid',/7x,'for the ',a,
     $        ' exchange model.')
              nerr = nerr + 1
            endif
          endif
c
c         Loop on condensed reactions read from the input file. Map
c         these to species corresponding to ions on the sites.
c
          qsperr = .false.
          do n = 1,ngexrt(je,ne)
            k = index(ugexr(n,je,ne),' == ')
            if (k .gt. 0) then
              ustr56 = ugexr(n,je,ne)(1:k - 1)
              call lejust(ustr56)
              ustr1 = ustr56(1:24)
              ustr56 = ugexr(n,je,ne)(k + 4:56)
              call lejust(ustr56)
              ustr2 = ustr56(1:24)
            else
              k = index(ugexr(n,je,ne),' = ')
              if (k .gt. 0) then
                ustr56 = ugexr(n,je,ne)(1:k - 1)
                call lejust(ustr56)
                ustr1 = ustr56(1:24)
                ustr56 = ugexr(n,je,ne)(k + 3:56)
                call lejust(ustr56)
                ustr2 = ustr56(1:24)
              else
                k = index(ugexr(n,je,ne),'==')
                if (k .gt. 0) then
                  ustr56 = ugexr(n,je,ne)(1:k - 1)
                  call lejust(ustr56)
                  ustr1 = ustr56(1:24)
                  ustr56 = ugexr(n,je,ne)(k + 2:56)
                  call lejust(ustr56)
                  ustr2 = ustr56(1:24)
                else
                  k = index(ugexr(n,je,ne),'=')
                  if (k .gt. 0) then
                    ustr56 = ugexr(n,je,ne)(1:k - 1)
                    call lejust(ustr56)
                    ustr1 = ustr56(1:24)
                    ustr56 = ugexr(n,je,ne)(k + 1:56)
                    call lejust(ustr56)
                    ustr2 = ustr56(1:24)
                  else
                    j2 = ilnobl(ugexr(n,je,ne))
                    write (noutpt,1080) ugexr(n,je,ne)(1:j2),
     $              ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
                    write (nttyo,1080) ugexr(n,je,ne)(1:j2),
     $              ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
 1080               format (/" * Error - (EQLIB/intexi) Can't decode",
     $              ' the condensed reaction string',/7x,'"',a,'",',
     $              ' which is given for site ',a,' of exchange',
     $              /7x,'phase ',a,'. The two species in the string',
     $              ' must be',/7x,'separated by " == ", " = ", "==",',
     $              ' or "=".')
                    ustr1 = 'Error'
                    ustr2 = 'Error'
                    nerr = nerr + 1
                    qsperr = .true.
                  endif
                endif
              endif
            endif
c
c           Prettify the strings in the ugexr array. Use "__" to
c           indicate a bare exchange site.
c
            if (ustr1(1:3) .eq. '   ') ustr1 = '__'
            if (ustr1(1:3) .eq. '_  ') ustr1 = '__'
            if (ustr2(1:3) .eq. '   ') ustr2 = '__'
            if (ustr2(1:3) .eq. '_  ') ustr2 = '__'
            ugexsr(1,n,je,ne) = ustr1
            ugexsr(2,n,je,ne) = ustr2
            j2 = ilnobl(ustr1)
            ugexr(n,je,ne)(j2 + 1:j2 + 3) = ' = '
            ugexr(n,je,ne)(j2 + 4:56) = ustr2
c
c           Make sure the two species in the reaction are not identical.
c
            if (ustr1(1:6) .ne. 'Error ') then
              if (ustr1(1:24) .eq. ustr2(1:24)) then
                j2 = ilnobl(ugexr(n,je,ne))
                if (ustr1(1:3) .eq. '__ ') then
                  write (noutpt,1090) ugexr(n,je,ne)(1:j2),
     $            ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
                  write (nttyo,1090) ugexr(n,je,ne)(1:j2),
     $            ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
 1090             format (/' * Error - (EQLIB/intexi) The condensed',
     $            ' exchange reaction "',a,'"',/7x,'for site ',a,
     $            ' of exchange phase ',a,' is an',/7x,'identity',
     $            " reaction. Such a reaction shouldn't appear on the",
     $            /7x,'input file.')
                else
                  write (noutpt,1100) ugexr(n,je,ne)(1:j2),
     $            ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
                  write (nttyo,1100) ugexr(n,je,ne)(1:j2),
     $            ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
 1100             format (/' * Error - (EQLIB/intexi) The condensed',
     $            ' exchange reaction "',a,'"',/7x,'for site ',a,
     $            ' of exchange phase ',a,' is invalid.',/7x,'The two',
     $            " species in the reaction can't be identical.")
                endif
                nerr = nerr + 1
                qsperr = .true.
              endif
            endif
          enddo
c
          if (qsperr) go to 200
c
c         Add a condensed reaction which involves one of the
c         exchangeable species and the bare site species, if none
c         such was included on the input.
c
          do n = 1,ngexrt(je,ne)
            ustr1 = ugexsr(1,n,je,ne)
            ustr2 = ugexsr(2,n,je,ne)
            if (ustr1(1:3) .eq. '__ ') go to 140
            if (ustr2(1:3) .eq. '__ ') go to 140
          enddo
c
          n = ngexrt(je,ne) + 1
          ngexrt(je,ne) = n
          ustr1 = ugexsr(1,1,je,ne)
          ustr2 = '__'
          ugexsr(1,n,je,ne) = ustr1
          ugexsr(2,n,je,ne) = ustr2
          j2 = ilnobl(ustr1)
          ugexr(n,je,ne) = ustr1
          ugexr(n,je,ne)(j2 + 1:j2 + 3) = ' = '
          ugexr(n,je,ne)(j2 + 4:56) = ustr2
          xlkgex(n,je,ne) = -12.0
          xhfgex(n,je,ne) = 0.
          xvfgex(n,je,ne) = 0.
          uxkgex(n,je,ne) = 'LogK/eq'
          uhfgex(n,je,ne) = 'kcal/eq'
          uvfgex(n,je,ne) = 'cm3/eq'
  140     continue
c
c         Add an identity reaction for the bare site species.
c         This is done to simplify the indexing relations between
c         exchanger species and corresponding reactions by making
c         them one to one. This simplification is important, because
c         the species must be created in a special order so that the
c         corresponding reaction properties for dissocation to the
c         bare site species can be calculated from the specified
c         exchange data without having to resort to matrix equations.
c
          n = ngexrt(je,ne) + 1
          ngexrt(je,ne) = n
          ustr1 = '__'
          ustr2 = '__'
          ugexsr(1,n,je,ne) = ustr1
          ugexsr(2,n,je,ne) = ustr2
          j2 = ilnobl(ustr1)
          ugexr(n,je,ne) = ustr1
          ugexr(n,je,ne)(j2 + 1:j2 + 3) = ' = '
          ugexr(n,je,ne)(j2 + 4:56) = ustr2
          xlkgex(n,je,ne) = 0.
          xhfgex(n,je,ne) = 0.
          xvfgex(n,je,ne) = 0.
          uxkgex(n,je,ne) = 'LogK/eq'
          uhfgex(n,je,ne) = 'kcal/eq'
          uvfgex(n,je,ne) = 'cm3/eq'
c
c         Put the names of the ions on the current site
c         in the ugexs array.
c
          nspect = 2
          ustr1 = ugexsr(1,1,je,ne)
          ustr2 = ugexsr(2,1,je,ne)
          ugexs(1,je,ne) = ustr1
          ugexs(2,je,ne) = ustr2
c
c         The last reaction is the identity reaction for the bare site
c         species, and it is guaranteed that a previous reaction
c         involves this species, so looping from 2,ngexrt(je,ne) - 1
c         is sufficient.
c
          do n = 2,ngexrt(je,ne) - 1
            do i = 1,2
              ustr1 = ugexsr(i,n,je,ne)
              do nn = 1,nspect
                if (ustr1(1:24) .eq. ugexs(nn,je,ne)) go to 150
              enddo
              nspect = nspect + 1
              ugexs(nspect,je,ne) = ustr1
  150         continue
            enddo
          enddo
c
          if (nspect .ne. ngexrt(je,ne)) then
            write (noutpt,1110) nspect,ngexrt(je,ne),
     $      ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
            write (nttyo,1110) nspect,ngexrt(je,ne),
     $      ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
 1110       format (/' * Error - (EQLIB/intexi) Have ',i3,' species',
     $      ' and ',i3,' reactions for',/7x,'site ',a,' of exchange',
     $      ' phase ',a,' after automatic',/7x,'additions to deal with',
     $      ' the bare site species. The number',/7x,'of species',
     $      ' (including the bare site species) must be equal to',
     $      /7x,'the number of reactions. Check the input. If you have',
     $      ' listed for',/7x,'a site for example two exchange',
     $      ' reactions such as "Na+ = K+" and',/7x,'"Rb+ = Cs+,"',
     $      ' you must include a third one to complete the linkage,',
     $       /7x,'such as "Na+ = Rb+".')
            nerr = nerr + 1
          endif
c
c         Map each exchange species to its aqueous species counterpart
c         (the exchangeable species).
c
          do ie = 1,nspect
            if (ugexs(ie,je,ne)(1:3) .eq. '__ ') then
              ngexsa(ie,je,ne) = 0
              go to 170
            endif
            j2 = ilnobl(ugexs(ie,je,ne))
c
            do nss = narn1,narn2
              if (ugexs(ie,je,ne)(1:24) .eq. uspec(nss)(1:24)) then
                ngexsa(ie,je,ne) = nss
                go to 170
              endif
            enddo
c
            write (noutpt,1130) ugexs(ie,je,ne)(1:j2),
     $      ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
            write (nttyo,1130) ugexs(ie,je,ne)(1:j2),
     $      ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
 1130       format (/' * Error - (EQLIB/intexi) The exchange species ',
     $      a,/7x,'specified for site ',a,' of exchange phase ',a,
     $      /7x,'does not correspond to any aqueous species read from',
     $      ' the data file.')
            nerr = nerr + 1
            qsperr = .true.
c
  170       continue
          enddo
c
          if (qsperr) go to 200
c
c         Test the electrical charges of the exchangeable species
c         and the charge of the site against the exchanger model.
c
          if (ugexmo(ne)(1:5).eq.'Gapon' .or.
     $      ugexmo(ne)(1:6).eq.'Gapon-' .or.
     $      ugexmo(ne)(1:8).eq.'Vanselow' .or.
     $      ugexmo(ne)(1:9).eq.'Vanselow-') then
            do ie = 1,nspect
              nss = ngexsa(ie,je,ne)
              if (nss .gt. 0) then
                zpx = zchar(nss)*zgexj(je,ne)
                if (zpx .ge. 0.) then
                  write (noutpt,1150) ugexs(ie,je,ne)(1:j2),
     $            zchar(nss),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4),
     $            ugexmo(ne)(1:j5)
                  write (nttyo,1150) ugexs(ie,je,ne)(1:j2),
     $            zchar(nss),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4),
     $            ugexmo(ne)(1:j5)
 1150             format (/' * Error - (EQLIB/intexi) The electrical',
     $            ' charge of ',a,/7x,'is ',f5.1,'. This species',
     $            " can't be exchanged onto site ",a,/7x,'of exchanger',
     $            ' phase ',a,' because the site',/7x,'has the same',
     $            ' charge sign or zero charge. Opposite charge',
     $            /7x,'signs are required for the ',a,' exchange',
     $            ' model.')
                  nerr = nerr + 1
                endif
              endif
            enddo
          endif
c
c         Convert the input thermodynamic parameters for exchange or
c         dissociation reactions to standard units.
c
          do ie = 1,nspect
            call lejust(uxkgex(ie,je,ne))
            if (uxkgex(ie,je,ne)(1:3) .eq. '   ') then
              uxkgex(ie,je,ne) = 'LogK/eq'
            elseif (uxkgex(ie,je,ne)(1:8) .eq. 'kcal/eq ') then
              xlkgex(ie,je,ne) = cfactr*xlkgex(ie,je,ne)
              uxkgex(ie,je,ne) = 'LogK/eq'
            elseif (uxkgex(ie,je,ne)(1:8) .eq. 'kJ/eq   ') then
              xlkgex(ie,je,ne) = xlkgex(ie,je,ne)/4.184
              xlkgex(ie,je,ne) = cfactr*xlkgex(ie,je,ne)
              uxkgex(ie,je,ne) = 'LogK/eq'
            elseif (uxkgex(ie,je,ne)(1:8) .ne. 'LogK/eq ') then
              j2 = ilnobl(ugexr(ie,je,ne))
              j5 = ilnobl(uxkgex(ie,je,ne))
              write (noutpt,1170) uxkgex(ie,je,ne)(1:j5),
     $        ugexr(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
              write (nttyo,1170) uxkgex(ie,je,ne)(1:j5),
     $        ugexr(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
 1170         format (/" * Error - (EQLIB/intexi) Don't recognize the",
     $        ' string "',a,'", which is',/7x,'given to define the',
     $        ' units for the xlkgex (log K) input for the',
     $        /7x,'reaction "',a,'" on site ',a,' of exchanger phase',
     $        /7x,a,'. The string must be blank (defaults to',
     $        ' "LogK/eq"),',/7x,'"LogK/eq", "kcal/eq", or "kJ/eq".')
              nerr = nerr + 1
            endif
c
            call lejust(uhfgex(ie,je,ne))
            if (uhfgex(ie,je,ne)(1:3) .eq. '   ') then
              uhfgex(ie,je,ne) = 'kcal/eq'
            elseif (uhfgex(ie,je,ne)(1:8) .eq. 'kJ/eq   ') then
              xhfgex(ie,je,ne) = xhfgex(ie,je,ne)/4.184
              uhfgex(ie,je,ne) = 'kcal/eq'
            elseif (uhfgex(ie,je,ne)(1:8) .ne. 'kcal/eq ') then
              j2 = ilnobl(ugexr(ie,je,ne))
              j5 = ilnobl(uhfgex(ie,je,ne))
              write (noutpt,1180) uhfgex(ie,je,ne)(1:j5),
     $        ugexr(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
              write (nttyo,1180) uhfgex(ie,je,ne)(1:j5),
     $        ugexr(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
 1180         format (/" * Error - (EQLIB/intexi) Don't recognize the",
     $        ' string "',a,'", which is',/7x,'given to define the',
     $        ' units for the xhfgex (enthalpy) input for the',
     $        /7x,'reaction "',a,'" on site ',a,' of exchanger phase',
     $        /7x,a,'. The string must be blank (defaults to',
     $        ' "kcal/eq"),',/7x,'"kcal/eq", or "kJ/eq".')
              nerr = nerr + 1
            endif
c
            call lejust(uvfgex(ie,je,ne))
            if (uvfgex(ie,je,ne)(1:3) .eq. '   ') then
              uvfgex(ie,je,ne) = 'cm3/eq'
            elseif (uvfgex(ie,je,ne)(1:8) .ne. 'cm3/eq  ') then
              j2 = ilnobl(ugexr(ie,je,ne))
              j5 = ilnobl(uvfgex(ie,je,ne))
              write (noutpt,1190) uvfgex(ie,je,ne)(1:j5),
     $        ugexr(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
              write (nttyo,1190) uvfgex(ie,je,ne)(1:j5),
     $        ugexr(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
 1190         format (/" * Error - (EQLIB/intexi) Don't recognize the",
     $        ' string "',a,'", which is',/7x,'given to define the',
     $        ' units for the xvfgex (volume) input for the',
     $        /7x,'reaction "',a,' on site ',a,' of exchanger phase',
     $        /7x,a,'. The string must be blank (defaults to',
     $        ' "cm3/eq")',/7x,'or "cm3/eq".')
              nerr = nerr + 1
            endif
          enddo
c
  200     continue
        enddo
      enddo
c
      if (nerr .gt. 0) stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      np = npt
      ns = nst
      nern1 = nst + 1
      iern1 = npt + 1
c
c     Pass 2. Loop on exchanger phases. Formally create exchanger phases
c     and species, and comute the corresponding reactions and
c     thermodynamic data for such species.
c
      do ne = 1,net
        j4 = ilnobl(ugexp(ne))
c
        if (npt .ge. nptmax) then
          write (noutpt,1240) nptmax,ugexp(ne)(1:j2)
          write (nttyo,1240) nptmax,ugexp(ne)(1:j2)
 1240     format (/' * Error - (EQLIB/intexi) The maximum ',i3,
     $    ' phases would be',/7x,'exceeded by creating the exchange',
     $    ' phase ',a,'. Increase',/7x,'the dimensioning parameter',
     $    ' nptpar.')
          stop
        endif
c
c       Add the current exchange phase to the list of phases.
c
        np = np + 1
        npt = np
        jpflag(np) = 0
        uphase(np) = ugexp(ne)
        ncmpr(1,np) = ns + 1
        tfxb = tgexp(ne)
c
c       Loop on sites.
c
        nsba = 0
        do je = 1,jgext(ne)
          j3 = ilnobl(ugexj(je,ne))
c
          nspect = ngexrt(je,ne)
c
c         Change the condensed reactions if necessary so that any
c         association reactions are converted to dissociation
c         reactions (e.g., for Na+, from "__ = Na+" to "Na+ = __").
c         Do not modify the last reaction, which is the identity
c         reaction for the bare site species.
c
          do ie = 1,nspect - 1
            if (ugexsr(1,ie,je,ne)(1:3) .eq. '__ ') then
              xlkgex(ie,je,ne) = -xlkgex(ie,je,ne)
              xhfgex(ie,je,ne) = -xhfgex(ie,je,ne)
              xvfgex(ie,je,ne) = -xvfgex(ie,je,ne)
              ugexsr(1,ie,je,ne)(1:24) = ugexsr(2,ie,je,ne)(1:24)
              ugexsr(2,ie,je,ne)(1:24) = '__ '
            endif
          enddo
c
c         Now convert all exchange reactions in this set to
c         dissociation reactions. To give an idea of what is happening,
c         a single exchange reaction "Ca++ = Na+" read from the input
c         file would now be expanded to the following condensed
c         reaction set:
c
c           1. Ca++ = Na+
c           2. Ca++ = __
c           3. __ = __
c
c         The corresponding species (corresponding in the set of having
c         matching indices) would be:
c
c           1. Ca++
c           2. Na+
c           3. __
c
c         Note that there is no formal correspondence in that species 2
c         (Na+) does not appear in reaction 2 (Ca++ = __). Also, the
c         base site species ("__") need not appear last, depending
c         on what condensed reactions were read from the input file.
c         The identity reaction for this species is guaranteed to
c         appear last, as it can't be read from the input file, and
c         is created in this position by the current subroutine in the
c         "Pass 1" section above.
c
c         The first step is to convert the set of existing condensed
c         reactions into a set of linearly independent dissociation
c         reactions. In the example, the set of reactions becomes:
c
c           1. Na+ = __
c           2. Ca++ = __
c           3. __ = __
c
c         Note that it is guaranteed that the last reaction is the
c         identity reaction for the bare site species and that at
c         least one of the remaining reactions is a dissociation
c         reaction. Tests performed above also guarantee that the
c         desired conversion is possible. Basically, the process
c         requires going through the reactions in an order which
c         insures that for each exchange reaction there currently
c         exists a dissociation reaction for one of the two species.
c
c         First, loop through the reactions, and mark all those which
c         are in the desired form (the identity reaction for the bare
c         site species and all reactions in dissociation format) by
c         setting a value of 1 in the corresponding element of the
c         ngexro array. It is guaranteed that there is at least one
c         reaction in dissociation format. Then repeatedly loop
c         through the remaining reactions, each time finding one to
c         transform from exchange to dissociation format by combining
c         it with a reaction in the processed set. Mark it as
c         transformed using the ngexro array. Repeated looping
c         is required to insure that all such reactions can be
c         transformed. Stop when no reactions remain in exchange
c         format.
c
c         The ngexro array will later be used to create a mapping
c         between the species and the reactions.
c
          do iee = 1,nspect
            ngexro(iee,je,ne) = 0
          enddo
c
          itot = 1
          ngexro(nspect,je,ne) = 1
c
          do ie = 1,nspect - 1
            ustr1 = ugexsr(1,ie,je,ne)
            ustr2 = ugexsr(2,ie,je,ne)
            if (ustr2(1:3) .eq. '__ ') then
              itot = itot + 1
              ngexro(ie,je,ne) = 1
            endif
          enddo
c
  210     if (itot .eq. nspect) go to 220
          do ie = 1,nspect - 1
            if (ngexro(ie,je,ne) .eq. 0) then
              ustr1 = ugexsr(1,ie,je,ne)
              ustr2 = ugexsr(2,ie,je,ne)
              do iee = 1,nspect - 1
                if (ngexro(iee,je,ne) .eq. 1) then
                  ustr1e = ugexsr(1,iee,je,ne)
                  if (ustr1e(1:24) .eq. ustr2(1:24)) then
                    itot = itot + 1
                    ngexro(ie,je,ne) = 1
                    ugexsr(1,ie,je,ne) = ustr1
                    ugexsr(2,ie,je,ne) = '__'
                    j2 = ilnobl(ustr1)
                    ugexr(ie,je,ne) = ustr1
                    ugexr(ie,je,ne)(j2 + 1:j2 + 3) = ' = '
                    ugexr(ie,je,ne)(j2 + 4:56) = '__'
                    xx = xlkgex(ie,je,ne) + xlkgex(iee,je,ne)
                    xlkgex(ie,je,ne) = xx
                    xx = xhfgex(ie,je,ne) + xhfgex(iee,je,ne)
                    xhfgex(ie,je,ne) = xx
                    xx = xvfgex(ie,je,ne) + xvfgex(iee,je,ne)
                    xvfgex(ie,je,ne) = xx
                    go to 210
                  elseif (ustr1e(1:24) .eq. ustr1(1:24)) then
                    itot = itot + 1
                    ngexro(ie,je,ne) = 1
                    ugexsr(1,ie,je,ne) = ustr2
                    ugexsr(2,ie,je,ne) = '__'
                    j2 = ilnobl(ustr2)
                    ugexr(ie,je,ne) = ustr2
                    ugexr(ie,je,ne)(j2 + 1:j2 + 3) = ' = '
                    ugexr(ie,je,ne)(j2 + 4:56) = '__'
                    xx = -xlkgex(ie,je,ne) + xlkgex(iee,je,ne)
                    xlkgex(ie,je,ne) = xx
                    xx = -xhfgex(ie,je,ne) + xhfgex(iee,je,ne)
                    xhfgex(ie,je,ne) = xx
                    xx = -xvfgex(ie,je,ne) + xvfgex(iee,je,ne)
                    xvfgex(ie,je,ne) = xx
                    go to 210
                  endif
                endif
              enddo
            endif
          enddo
c
          write (noutpt,1250) ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
          write (nttyo,1250) ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
 1250     format(/' * Error - (EQLIB/intexi) Programming error trap:',
     $    " Couldn't transform",/7x,'one or more  condensed reactions',
     $    ' for site ',a,' of exchange',/7x,'phase ',a,' from',
     $    ' exchange format to dissociation',/7x,'format. Check the',
     $    ' responsible coding.')
          stop
c
  220     continue
c
c         Now set up the ngexro array as a pointer array giving the
c         index of the reaction corresponding to a given species.
c
          itot = 0
          do ie = 1,nspect
            ngexro(ie,je,ne) = 0
            ustr1 = ugexs(ie,je,ne)
            do iee = 1,nspect
              if (ustr1(1:24) .eq. ugexsr(1,iee,je,ne)(1:24)) then
                itot = itot + 1
                ngexro(ie,je,ne) = iee
                go to 230
              endif
            enddo
c
            write (noutpt,1260) ugexs(ie,je,ne)(1:j2),
     $      ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
            write (nttyo,1260) ugexs(ie,je,ne)(1:j2),
     $      ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
 1260       format(/' * Error - (EQLIB/intexi) Programming error trap:',
     $      " Couldn't find",/7x,'a condensed identity or dissociation',
     $      ' reaction for ',a,/7x,'on site ',a,' of exchange phase ',
     $      a,'. Check the',/7x,'responsible coding.')
            stop
c
  230       continue
          enddo
c
c         Order the species for creation, using the ngexso array.
c         Use the existing order, except that the bare site species
c         is created first.
c
          itot = 0
          do ie = 1,nspect
            if (ugexs(ie,je,ne)(1:3) .eq. '__ ') then
              itot = itot + 1
              ngexso(itot,je,ne) = ie
              go to 240
            endif
          enddo
  240     continue
          do ie = 1,nspect
            if (ugexs(ie,je,ne)(1:3) .ne. '__ ') then
              itot = itot + 1
              ngexso(itot,je,ne) = ie
            endif
          enddo
c
c         Create the species for the current site.
c
          jern1(je,ne) = ns + 1
c
c         Loop on ions on sites. Note that the bare exchanger species
c         will be created first.
c
          do ieo = 1,nspect
            ie = ngexso(ieo,je,ne)
            j2 = ilnobl(ugexs(ie,je,ne))
c
            if (nst .ge. nstmax) then
              write (noutpt,1290) nstmax,ugexs(ie,je,ne)(1:j2),
     $        ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
              write (nttyo,1290) nstmax,ugexs(ie,je,ne)(1:j2),
     $        ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
 1290         format (/' * Error - (EQLIB/intexi) The maximum ',i4,
     $        ' species would be',/7x,'exceeded by creating an',
     $        ' exchange species for ',a,/7x,'on site ',a,
     $        ' of exchange phase ',a,'. Increase the',
     $        /7x,'dimensioning parameter nstpar.')
              stop
            endif
c
            ns = ns + 1
            nsl = nst
            nst = ns
            jsflag(ns) = 0
c
            jsitex(ns) = je
            nphasx(ns) = np
c
c           If the current species is the bare one for the current
c           site, add it to the set of strict basis species.
c
            nss = ngexsa(ie,je,ne)
            if (nss .eq. 0) then
              nbt = nbt + 1
              if (nbt .gt. nbtmax) then
                write (noutpt,1300) nbtmax,ugexs(ie,je,ne)(1:j2),
     $          ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
                write (nttyo,1300) nbtmax,ugexs(ie,je,ne)(1:j2),
     $          ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
 1300           format (/' * Error - (EQLIB/intexi) The maximum,',
     $          i3,' basis species would be',/7x,'exceeded by',
     $          ' adding a species for ',a,' on site ',a,/7x,'of',
     $          ' exchange phase ',a,'. Increase the dimensioning',
     $          /7x,'parameter nbtpar.')
                stop
              else
                nbasp(nbt) = ns
                jflag(ns) = 0
                nsba = ns
              endif
            endif
c
c           Get the index of the corresponding aqueous species
c           (the exchangeable species). Compute the stoichiometric
c           relationship between the exchanger species and the
c           corresponding aqueous species. This relationship depends
c           on the specified exchange model.
c
            nss = ngexsa(ie,je,ne)
            if (nss .gt. 0) then
              zsi = zchar(nss)
            else
              zsi = 0.
            endif
            zxj = zgexj(je,ne)
            zprod = zsi*zxj
c
c           Note: cfxi below is the reaction coefficient for the
c           exchangeable ion in the dissociation reaction for the
c           corresponding exchanger species. The dissociation reaction
c           is always written such that one mole of exchanger species
c           is dissociated.
c
            j5 = ilnobl(ugexmo(ne))
c
            if (nss .gt. 0) then
              if (ugexmo(ne)(1:5).eq.'Gapon' .or.
     $          ugexmo(ne)(1:6).eq.'Gapon-' .or.
     $          ugexmo(ne)(1:8).eq.'Vanselow' .or.
     $          ugexmo(ne)(1:9).eq.'Vanselow-') then
c
c               Check electrical charges of the exchangeable ion and
c               the exchange site.
c
                qmoerr = .false.
                if (zsi .eq. 0.) then
                  write (noutpt,1320) ugexs(ie,je,ne)(1:j2),
     $            ugexj(je,ne)(1:j3),ugexp(ne)(1:j4),ugexmo(ne)(1:j5)
                  write (nttyo,1320) ugexs(ie,je,ne)(1:j2),
     $            ugexj(je,ne)(1:j3),ugexp(ne)(1:j4),ugexmo(ne)(1:j5)
 1320             format (/" * Error - (EQLIB/intexi) Can't create",
     $            ' a generic ion exchange species',/7x,'for ',a,
     $            ' on site ',a,' of exchange phase',/7x,a,
     $            ' for the ',a,' model because',/7x,'the',
     $            ' exchangeable ion has no electrical charge.')
                  nerr = nerr + 1
                  qmoerr = .true.
                endif
c
                if (zxj .eq. 0.) then
                  write (noutpt,1330) ugexs(ie,je,ne)(1:j2),
     $            ugexj(je,ne)(1:j3),ugexp(ne)(1:j4),ugexmo(ne)(1:j5)
                  write (nttyo,1330) ugexs(ie,je,ne)(1:j2),
     $            ugexj(je,ne)(1:j3),ugexp(ne)(1:j4),ugexmo(ne)(1:j5)
 1330             format (/" * Error - (EQLIB/intexi) Can't create",
     $            ' a generic ion exchange species',/7x,'for ',a,
     $            ' on site ',a,' of exchange phase',/7x,a,
     $            ' for the ',a,' model because',/7x,'the',
     $            ' the exchange site has no electrical charge.')
                  nerr = nerr + 1
                  qmoerr = .true.
                endif
c
                if (qmoerr) go to 270
c
                if (zprod .gt. 0.) then
                  write (noutpt,1340) ugexs(ie,je,ne)(1:j2),
     $            ugexj(je,ne)(1:j3),ugexp(ne)(1:j4),ugexmo(ne)(1:j5)
                  write (nttyo,1340) ugexs(ie,je,ne)(1:j2),
     $            ugexj(je,ne)(1:j3),ugexp(ne)(1:j4),ugexmo(ne)(1:j5)
 1340             format (/" * Error - (EQLIB/intexi) Can't create",
     $            ' a generic ion exchange species',/7x,'for ',a,
     $            ' on site ',a,' of exchange phase',/7x,a,' for the ',
     $            a,' model because the exchangeable',/7x,'ion and',
     $            ' the exchange site have the same electrical',
     $            ' charge sign.')
                  nerr = nerr + 1
                  go to 270
                endif
              endif
            endif
c
            cgx = cgexj(je,ne)
            zxjt = cgx*zxj
c
c           cfxi = the number of moles of exchangeable ion per
c                    mole of the exchanger species
c           cfxz = the number of moles of exchanger substrate per
c                    mole of exchanger species
c
            if (nss .gt. 0) then
c
c             Have other than the bare site species.
c
              if (ugexmo(ne)(1:j5).eq.'Gapon' .or.
     $          ugexmo(ne)(1:6).eq.'Gapon-') then
c
c               Gapon (Gapon-?) model.
c
c                 For:
c
c                   cgexj(je,ne) = 1
c                   zgexj(je,ne) = -1
c                   "Na+ = Ca++"
c
c                 the species are:
c
c                   Na-Z and Ca[1/2]-Z
c
                cfxi = -cgx*zxj/zsi
                cfxz = -cfxi*zsi/zxjt
c
              elseif (ugexmo(ne)(1:j5).eq.'Vanselow' .or.
     $          ugexmo(ne)(1:9).eq.'Vanselow-') then
c
c               Vanselow (Vanselow-?) model.
c
c                 For:
c
c                   cgexj(je,ne) = 1
c                   zgexj(je,ne) = -1
c                   "Na+ = Ca++"
c
c                 the species are:
c
c                   Na-Z and Ca-Z2
c
                cfxi = -cgx*zxj*abs(zsi)/zsi
                cfxz = -cfxi*zsi/zxjt
c
              elseif (ugexmo(ne)(1:j5) .eq. 'Site-mixing') then
c
c               Site-mixing model.
c
                cfxi = cgx
                cfxz = 1.0
c
              else
                write (noutpt,1350) ugexmo(ne)(1:j5),ugexp(ne)(1:j4)
                write (nttyo,1350) ugexmo(ne)(1:j5),ugexp(ne)(1:j4)
 1350           format (/' * Error - (EQLIB/intexi.f) Programming',
     $          " error trap: Don't recognize the",/7x,'exchange',
     $          ' model type"',a,'" specified for the exchanger phase',
     $          /7x,a,'. Valid choices include "Gapon" and "Vanselow".',
     $          /7x,'Invalid choices should have been trapped above in',
     $          ' this subroutine.',/7x,'Check the coding at the',
     $          ' present point. Now trying to compute',/7x,'factors',
     $          ' for the compositions and reactions of species',
     $          ' belonging',/7x,'to this exchanger phase.')
                stop
              endif
            else
c
c             Have the bare site species.
c
              cfxi = cgx
              cfxz = 1.0
            endif
c
c           Here efx is the number of equivalents of exchangeable
c           ion appearing in the dissociation reaction.
c
            if (nss .gt. 0) then
c
c             Have other than the bare site species.
c
              zchar(ns) = cfxi*zsi - zxjt
              mwtsp(ns) = cfxi*mwtsp(nss)
              efx = abs(cfxi*zsi)
            else
c
c             Have the bare site species.
c
              zchar(ns) = zxjt
              mwtsp(ns) = 0.
              efx = 0.
            endif
c
c           Compose the full name of the exchanger species. The phase
c           part of the name is the exchange phase name. The species
c           part is a concatenation of the species part of the name of
c           the ion released in an exchange reaction (the aqueous
c           ion) and the site name. A blank space is included in the
c           concatenation between the two parts. The concatenation is
c           trimmed as necessary to fit into 24 characters.
c           Approximately two characters of the exchange ion name are
c           trimmed for each character of the site name. The full name
c           of the species should should then look like
c           "Na+ S(1)                Exchanger(A)". This algorithm
c           matches that in EQLIB/intgex.f. As formatted by
c           EQLIBU/fmspnm.f or EQLIBU/fmspnx.f, the name in this
c           example would appear as "Na+ S(1) Exchanger(A)".
c
c           Warning: the stoichiometry of the species can't be deduced
c           from the name alone. The the number of formula units of
c           exchangeable Na+ in site S(1) in the above example could be
c           1.0, 2.0, or some other number. The same is true for the
c           number of moles of the site itself per mole of exchanger.
c           Thus, "Na+ S(1) Exchanger(A)" has one stoichiometry for
c           the Vanselow model, but a different one for the Gapon model.
c
            jj2 = j2
            jj3 = j3
            nd2 = 0
  250       jj = jj2 + jj3 + 1
            if (jj .gt. 24) then
              if (nd2 .lt. 2) then
                jj2 = jj2 - 1
                nd2 = nd2 + 1
                go to 250
              else
                jj3 = jj3 - 1
                nd2 = 0
                go to 250
              endif
            endif
            uspec(ns) = ' '
            uspec(ns)(1:jj2) = ugexs(ie,je,ne)(1:jj2)
            uspec(ns)(jj2 + 1:jj2 + 1) = ' '
            uspec(ns)(jj2 + 2:jj) = ugexj(je,ne)(1:jj3)
            uspec(ns)(25:48) = ugexp(ne)(1:24)
c
c           Set up the elemental composition.
c
            nrf1 = nessr(2,nsl) + 1
            if (nss .gt. 0) then
              nr1 = nessr(1,nss)
              nr2 = nessr(2,nss)
              nrf2 = nrf1 + nr2 - nr1
            else
              nrf2 = nrf1
            endif
c
            if (nrf2 .gt. nessmx) then
              write (noutpt,1380) nessmx,ugexs(ie,je,ne)(1:j2),
     $        ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
              write (nttyo,1380) nessmx,ugexs(ie,je,ne)(1:j2),
     $        ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
 1380         format (/' * Error - (EQLIB/intexi) The maximum ',i5,
     $        ' entries in the',/7x,'cess/ness arrays would be',
     $        ' exceeded by creating the',/7x,'exchange species',
     $        ' for ',a,' on site',/7x,a,' of exchange phase ',a,'.',
     $        ' Increase the',/7x,'dimensioning parameter nesspa.')
              stop
            endif
c
            nessr(1,ns) = nrf1
            nessr(2,ns) = nrf2
            if (nss .gt. 0) then
              k = nr1 - 1
              do n = nrf1,nrf2
                k = k + 1
                cess(n) = cfxi*cess(k)
                ness(n) = ness(k)
              enddo
            else
              n = nrf1
              cess(n) = 0.
              ness(n) = 0
            endif
c
c           Set up the corresponding reaction.
c
            nrf1 = ndrsr(2,nsl) + 1
            if (nss .gt. 0) then
c
c             Case of a non-bare site species.
c
              nrf2 = nrf1 + 2
            else
c
c             Case of the bare site species.
c
              nrf2 = nrf1
            endif
c
            if (nrf2 .gt. ndrsmx) then
              write (noutpt,1390) ndrsmx,ugexs(ie,je,ne)(1:j2),
     $        ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
              write (nttyo,1390) ndrsmx,ugexs(ie,je,ne)(1:j2),
     $        ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
 1390         format (/' * Error - (EQLIB/intexi) The maximum ',i5,
     $        ' entries in the',/7x,'cdrs/ndrs arrays would be',
     $        ' exceeded by creating the',/7x,'exchange species',
     $        ' for ',a,'on site',/7x,a,' of exchange phase ',a,'.',
     $        ' Increase the',/7x,'dimensioning parameter ndrspa.')
              stop
            endif
c
            ndrsr(1,ns) = nrf1
            ndrsr(2,ns) = nrf2
            if (nss .gt. 0) then
c
c             Case of a non-bare site species. Create a reaction in
c             which the species dissociates to the bare site species.
c
              n = nrf1
              cdrs(n) = -1.0
              ndrs(n) = ns
              n = n + 1
              cdrs(n) = cfxi
              ndrs(n) = nss
              n = n + 1
              cdrs(n) = cfxz
              ndrs(n) = nsba
            else
c
c             Case of the bare site species. Create a null reaction.
c             The bare site species then becomes a strict basis
c             species, though it doesn't actually formally correspond
c             to a chemical element.
c
              n = nrf1
              cdrs(n) = 0.
              ndrs(n) = 0
            endif
c
            if (nss .eq. 0) go to 270
c
c           Set up to the thermodynamic properties of the reaction
c           just created. This reaction should match one of the
c           existing condensed reactions.
c
            iee = ngexro(ie,je,ne)
c
c           Calculate the thermodynamic properties of the reaction.
c
            do j = 1,ntprt
              if (efx .ne. 0.) then
                afxb = efx*xlkgex(iee,je,ne)
                afhx = efx*xhfgex(iee,je,ne)
                afvx = efx*xvfgex(iee,je,ne)
              else
                afxb = xlkgex(iee,je,ne)
                afhx = xhfgex(iee,je,ne)
                afvx = xvfgex(iee,je,ne)
              endif
c
c             Here afxb is the log K value and tfxb the temperature (C)
c             at the base temperature. Use the van't Hoff relation to
c             find afx0, the log K at 0C. Approximate the van't Hoff
c             relation across the total temperature range using the
c             standard power series.
c
              bfx = -afhx/(arcnst*tfx0)
              afx0 = afxb + bfx*(tfxb/(tfxb + tfx0))
              axlks(1,j,ns) = afx0
              xx = -bfx
              do i = 2,narxt(j)
                xx = -xx/tfx0
                axlks(i,j,ns) = xx
              enddo
c
              axhfs(1,j,ns) = afhx
              axvfs(1,j,ns) = afvx
              do i = 2,narxt(j)
                axhfs(i,j,ns) = 0.
                axvfs(i,j,ns) = 0.
              enddo
            enddo
c
c           Is the aqueous species which goes onto the substrate a
c           basis species? If not, the reaction must be rewritten
c           in terms of equivalent basis species.
c
            nrf1 = ndrsr(1,ns)
            nrf2 = ndrsr(2,ns)
            do n = nrf1 + 1,nrf2
              nse = ndrs(n)
              if (nse.ge.narn1 .and. nse.le.narn2) then
c
c               Calling sequence substitutions:
c                 nse for ns
c
                nb = nbasis(nbasp,nbt,nbtmax,nse)
                if (nb .eq. 0) then
                  nbt = nbt + 1
                  if (nbt .gt. nbtmax) then
                    j6 = ilnobl(uspec(nse))
                    write (noutpt,1400) nbtmax,uspec(nse)(1:j6),
     $              ugexs(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),
     $              ugexp(ne)(1:j4)
                    write (nttyo,1400) nbtmax,uspec(nse)(1:j6),
     $              ugexs(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),
     $              ugexp(ne)(1:j4)
 1400               format (/' * Error - (EQLIB/intexi) The maximum,',
     $              i3,' basis species would be',/7x,'exceeded by',
     $              ' adding ',a,' in order to',/7x,'accomodate the',
     $              ' required setup for species ',a,/7x,'on site ',a,
     $              ' of exchange phase ',a,'.',/7x,'Increase the',
     $              ' dimensioning parameter nbtpar.')
                    stop
                  else
                    nbasp(nbt) = nse
                    jflag(nse) = 30
                  endif
                endif
              endif
            enddo
c
  270       continue
          enddo
c
c         Modify the ngexsa array if necessary so that the species
c         index of the bare site species appears in the first position
c         for the current site of the current exchange phase.
c
          do ie = 1,nspect
            nss = ngexsa(ie,je,ne)
            if (nss .eq. 0) then
              do i = 1,ie - 1
                iee = ie - i
                ngexsa(iee + 1,je,ne) = ngexsa(iee,je,ne)
              enddo
              ngexsa(1,je,ne) = 0
              go to 280
            endif
          enddo
  280     continue
c
          jern2(je,ne) = nst
          ngext(je,ne) = jern2(je,ne) - jern1(je,ne) + 1
        enddo
        ncmpr(2,np) = nst
      enddo
c
      nern2 = nst
      iern2 = npt
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Pass 3. Loop on exchanger phases. Change the ordering of the
c     condensed reactions and associated data to match that of the
c     corresponding created species. Then eliminate the identity
c     reaction for the bare site species. This puts these data
c     in the desired form for writing on a pickup file.
c
c     In the example discussed above (input of a single exchange
c     reaction, Ca++ = Na+), the expanded set of condensed reactions
c     was:
c
c           1. Na+ = __
c           2. Ca++ = __
c           3. __ = __
c
c     The actual species created were processed in the order
c     __, Ca++, Na+. The bare site species is always processed first.
c     The other species are processed in order of appearance in the
c     condensed reactions read from the input file. Thus, the order of
c     creation of the corresponding reactions from the condensed
c     counterparts was:
c
c           1. (3. __ = __)
c           2. (2. Ca++ = __)
c           3. (1. Na+ = __)
c
c     Here the set of condensed reactions is changed to:
c
c           1. Ca++ = __
c           2. Na+ = __
c
c     As input, this is expanded to:
c
c           1. Ca++ = __
c           2. Na+ = __
c           3. __ = __
c
c     The actual species are then created in the same order as before,
c     (__, Ca++, Na+).
c
c     Note that less setup work is required for the case of the input
c     of two dissociation reactions than for the case of the input of
c     the single equivalent exchange reaction.
c
      do ne = 1,net
c
c       Loop on sites.
c
        do je = 1,jgext(ne)
          nspect = ngexrt(je,ne)
c
c         Re-arrange the condensed reactions so that their order
c         matches that of the local species list. This eliminates
c         the need for the ngexro pointer array to find the
c         condensed reaction for a given species in the local
c         species list.
c
          do ie = 1,nspect
            iee = ngexro(ie,je,ne)
            if (iee .ne. ie) then
c
c             Exchange positions.
c
              ustr56 = ugexr(iee,je,ne)
              ugexr(iee,je,ne) = ugexr(ie,je,ne)
              ugexr(ie,je,ne) = ustr56
              xx = xlkgex(iee,je,ne)
              xlkgex(iee,je,ne) = xlkgex(ie,je,ne)
              xlkgex(ie,je,ne) = xx
              xx = xhfgex(iee,je,ne)
              xhfgex(iee,je,ne) = xhfgex(ie,je,ne)
              xhfgex(ie,je,ne) = xx
              xx = xvfgex(iee,je,ne)
              xvfgex(iee,je,ne) = xvfgex(ie,je,ne)
              xvfgex(ie,je,ne) = xx
              do iej = 1,nspect
                if (ngexro(iej,je,ne) .eq. ie) go to 310
              enddo
  310         ngexro(iej,je,ne) = iee
              ngexro(ie,je,ne) = ie
            endif
          enddo
c
          do ie = 1,nspect
            iee = ngexso(ie,je,ne)
            if (iee .ne. ie) then
c
c             Exchange positions.
c
              ustr56 = ugexr(iee,je,ne)
              ugexr(iee,je,ne) = ugexr(ie,je,ne)
              ugexr(ie,je,ne) = ustr56
              xx = xlkgex(iee,je,ne)
              xlkgex(iee,je,ne) = xlkgex(ie,je,ne)
              xlkgex(ie,je,ne) = xx
              xx = xhfgex(iee,je,ne)
              xhfgex(iee,je,ne) = xhfgex(ie,je,ne)
              xhfgex(ie,je,ne) = xx
              xx = xvfgex(iee,je,ne)
              xvfgex(iee,je,ne) = xvfgex(ie,je,ne)
              xvfgex(ie,je,ne) = xx
              do iej = 1,nspect
                if (ngexso(iej,je,ne) .eq. ie) go to 320
              enddo
  320         ngexso(iej,je,ne) = iee
              ngexso(ie,je,ne) = ie
            endif
          enddo
c
c         Remove the identity reaction from the set of condensed
c         reactions.
c
          do ie = 1,nspect - 1
            iee = ie + 1
            ugexr(ie,je,ne) = ugexr(iee,je,ne)
            xlkgex(ie,je,ne) = xlkgex(iee,je,ne)
            xhfgex(ie,je,ne) = xhfgex(iee,je,ne)
            xvfgex(ie,je,ne) = xvfgex(iee,je,ne)
          enddo
          ie = nspect
          ugexr(ie,je,ne) = '__ = __ '
          xlkgex(ie,je,ne) = 0.
          xhfgex(ie,je,ne) = 0.
          xvfgex(ie,je,ne) = 0.
          nspect = nspect - 1
          ngexrt(je,ne) = nspect
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the cegexs array. This is an array of coefficients
c     giving the number of equivalents per mole of each exchanger
c     species. This supports the calculation of equivalent fractions.
c
      do ne = 1,net
        do je = 1,jgext(ne)
          ns = jern1(je,ne) - 1
          do ie = 1,ngext(je,ne)
            ns = ns + 1
            if (uspec(ns)(1:3) .eq. '__ ') then
            else
c
c             Have the bare site species.
c
              cegexs(ie,je,ne) = 0.
c
c             Have other than the bare site species.
c
              nr1 = ndrsr(1,ns)
              nr2 = ndrsr(2,ns)
              cx = 0.
              do n = nr1 + 1,nr2
                nse = ndrs(n)
c
c               There are contributions only from the aqueous
c               species. There is none from the bare site species,
c               which also appears in the reaction.
c
                if (nse.ge.narn1 .and. nse.le.narn2) then
                  cx = cx + cdrs(n)*zchar(nse)
                endif
              enddo
              cegexs(ie,je,ne) = cx
            endif
          enddo
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the cpgexs array. This is an array of coefficients
c     giving the number of moles of exchanger substrate (Z) per mole of
c     exchanger species (e.g., Na-Z, Ca-Z2). This supports the
c     calculation of the number of moles of the exchanger phase.
c
      do ne = 1,net
        j5 = ilnobl(ugexmo(ne))
        if (ugexmo(ne)(1:j5).eq.'Gapon' .or.
     $    ugexmo(ne)(1:6).eq.'Gapon-' .or.
     $    ugexmo(ne)(1:j5) .eq. 'Site-mixing') then
c
c         Gapon (Gapon-?) or Site-mixing model.
c
          do je = 1,jgext(ne)
            do ie = 1,ngext(je,ne)
              cpgexs(ie,je,ne) = 1.0
            enddo
          enddo
c
        elseif (ugexmo(ne)(1:j5).eq.'Vanselow' .or.
     $    ugexmo(ne)(1:9).eq.'Vanselow-') then
c
c         Vanselow (Vanselow-?) model.
c
          do je = 1,jgext(ne)
            do ie = 1,ngext(je,ne)
              cpgexs(ie,je,ne) = cegexs(ie,je,ne)/egexjf(je,ne)
            enddo
          enddo
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set up the kgexsa, kern1, and kern2 pointer arrays. The kgexsa
c     array contains the indices of the exchange species (e.g., Na+).
c     The latter two arrays point to the start and end of the range
c     in kgexsa corresponding to a given generic ion exchanger phase.
c     Here kgexsa(ke,ne) is the index of the ke-th species exchanging
c     on the ne-th generic ion exchanger phase.
c
      ke = 0
      do ne = 1,net
        ke1 = ke + 1
        kern1(ne) = ke1
c
        do je = 1,jgext(ne)
          do ie = 1,ngext(je,ne)
            nss = ngexsa(ie,je,ne)
c
c           Has this exchange species already been loaded because
c           because it appears in another site?
c
            do kee = ke1,ke
              nsse = kgexsa(kee,ne)
              if (nsse .eq. nss) go to 400
            enddo
c
c           Have an exchange species which has not already been
c           loaded. Load it.
c
            ke = ke + 1
            kgexsa(ke,ne) = nss
  400       continue
          enddo
        enddo
c
        kern2(ne) = ke
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  990 if (nerr .gt. 0) stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
