      subroutine cmpdat(adhfs,adhfsd,advfs,advfsd,amu,amua,apx,
     $ apxa,aslm,aslma,atwt,atwta,axhfs,axhfsd,axlks,axlksd,axvfs,
     $ axvfsd,azero,azeroa,bpx,bpxa,cdrs,cdrsd,cess,cessa,iapxmx,
     $ iapxt,iapxta,iaqsla,iaqsln,ibpxmx,ibpxt,ibpxta,ilrn1,ilrn2,
     $ imrn1,imrn2,ipbtmx,ipch,ipchmx,ipcv,ipcvmx,insgf,insgfa,
     $ iopg,ixrn1,ixrn1a,ixrn2,ixrn2a,jflag,jpfcmx,jpflag,jsflag,
     $ jsitex,jsol,jsola,mwtsp,mwtspa,nalpaa,nalpha,napmax,napt,
     $ napta,narn1,narn1a,narn2,narn2a,narxmx,nat,natmax,nbasp,
     $ nbaspd,nbmap,nbt,nbtd,nbti,nbtmax,nchlor,ncmap,ncmpr,ncmpra,
     $ nct,ncta,nctmax,ndecsp,ndrs,ndrsd,ndrsmx,ndrsr,ndrsrd,ness,
     $ nessa,nessmx,nessr,nessra,ngrn1,ngrn1a,ngrn2,ngrn2a,ngt,
     $ nlrn1,nlrn1a,nlrn2,nlrn2a,nlt,nmrn1,nmrn1a,nmrn2,nmrn2a,nmt,
     $ nmut,nmuta,nmutmx,nmux,nmuxa,nopgmx,nslt,nsltmx,nslta,nslx,
     $ nslxa,nphasx,npt,npta,nptmax,nsmap,nsta,nst,nstmax,ntf1,ntf1a,
     $ ntf1mx,ntf1t,ntf1ta,ntf2,ntf2a,ntf2mx,ntf2t,ntf2ta,ntprmx,
     $ nxrn1,nxrn1a,nxrn2,nxrn2a,nxt,nxtmax,palpaa,palpha,qchlor,
     $ tf1,tf1a,tf2,tf2a,uelem,uelema,uphasa,uphase,uptypa,uptype,
     $ uspec,uspeca,vosp0,vosp0a,zchar,zchara)
c
c     This sburoutine compresses the data read from the data file. The
c     objective is to write data arrays (and associated variables)
c     which exclude data belonging to phases and species that are not
c     required by the current problem. Such phases are marked by
c     the condition jpflag = 2. If a phase is not needed, none of its
c     constituent species are needed. If a phase is needed, any
c     constituent species that are not needed are marked by the
c     condition jsflag = 2.
c
c     The data read from the data file are stored in arrays (and a
c     few important variables) which are characterized by an 'a'
c     suffix, but which otherwise have the same names as the working
c     arrays (and associated variables) which are created by this
c     compression. Thus, the number of species read from the data
c     file is nsta. The number remaining after compression is nst.
c     The array of reaction coefficients read from the data file is
c     cdrsd. The array created by the compression is cdrs.
c
c     During compression, the basis set is compressed from nbtd elements
c     in the nbaspd array to nbt elements in the nbasp array. Here nbasp
c     and nbt are basically used as workspace. At the end of this
c     subroutine, the nbaspd is set equal to the compressed nbasp array
c     and nbtd is set equal to nbt.
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
      integer iapxmx,ibpxmx,ipbtmx,ipchmx,ipcvmx,jpfcmx,napmax,narxmx,
     $ natmax,nbtmax,nctmax,ndrsmx,nessmx,nmutmx,nopgmx,nptmax,nsltmx,
     $ nstmax,ntf1mx,ntf2mx,ntprmx,nxtmax
c
      integer iapxt(nxtmax),iapxta(nxtmax),ibpxt(nxtmax),ibpxta(nxtmax),
     $ insgf(natmax),insgfa(natmax),iopg(nopgmx),jflag(nstmax),
     $ jpflag(nptmax),jsflag(nstmax),jsitex(nstmax),jsol(nxtmax),
     $ jsola(nxtmax),nalpaa(nsltmx),nalpha(nsltmx),nbasp(nbtmax),
     $ nbaspd(nbtmax),nbmap(nbtmax),ncmap(nctmax),ncmpr(2,nptmax),
     $ ncmpra(2,nptmax),ndecsp(nbtmax),ndrs(ndrsmx),ndrsd(ndrsmx),
     $ ndrsr(2,nstmax),ndrsrd(2,nstmax),ness(nessmx),nessa(nessmx),
     $ nessr(2,nstmax),nessra(2,nstmax),nmux(3,nmutmx),nmuxa(3,nmutmx),
     $ nphasx(nstmax),nslx(2,nsltmx),nslxa(2,nsltmx),nsmap(nstmax),
     $ ntf1(ntf1mx),ntf1a(ntf1mx),ntf2(ntf2mx),ntf2a(ntf2mx)
c
      integer iaqsla,iaqsln,ilrn1,ilrn2,imrn1,imrn2,ipch,ipcv,ixrn1,
     $ ixrn1a,ixrn2,ixrn2a,napt,napta,narn1,narn1a,narn2,narn2a,nat,nbt,
     $ nbtd,nbti,nchlor,nct,ncta,ngrn1,ngrn1a,ngrn2,ngrn2a,ngt,nlrn1,
     $ nlrn1a,nlrn2,nlrn2a,nlt,nmrn1,nmrn1a,nmrn2,nmrn2a,nmt,nmut,nmuta,
     $ nslt,nslta,npt,npta,nsta,nst,ntf1t,ntf1ta,ntf2t,ntf2ta,nxrn1,
     $ nxrn1a,nxrn2,nxrn2a,nxt
c
      logical qchlor
c
      character*8 uelem(nctmax),uelema(nctmax)
      character*48 uspec(nstmax),uspeca(nstmax)
      character*24 uphasa(nptmax),uphase(nptmax),uptypa(nptmax),
     $ uptype(nptmax)
c
      real*8 adhfs(narxmx,ntprmx,ipchmx,nstmax),
     $ adhfsd(narxmx,ntprmx,ipchmx,nstmax),
     $ advfs(narxmx,ntprmx,ipcvmx,nstmax),
     $ advfsd(narxmx,ntprmx,ipcvmx,nstmax),apx(iapxmx,nxtmax),
     $ apxa(iapxmx,nxtmax),atwt(nctmax),atwta(nctmax),
     $ axhfs(narxmx,ntprmx,nstmax),axhfsd(narxmx,ntprmx,nstmax),
     $ axlks(narxmx,ntprmx,nstmax),axlksd(narxmx,ntprmx,nstmax),
     $ axvfs(narxmx,ntprmx,nstmax),axvfsd(narxmx,ntprmx,nstmax),
     $ azero(natmax),azeroa(natmax)
c
      real*8 amu(jpfcmx,nmutmx),amua(jpfcmx,nmutmx),
     $ aslm(jpfcmx,0:ipbtmx,nsltmx),aslma(jpfcmx,0:ipbtmx,nsltmx),
     $ bpx(ibpxmx,nxtmax),bpxa(ibpxmx,nxtmax),cdrs(ndrsmx),
     $ cdrsd(ndrsmx),cess(nessmx),cessa(nessmx),mwtsp(nstmax),
     $ mwtspa(nstmax),palpaa(ipbtmx,napmax),palpha(ipbtmx,napmax),
     $ tf1(ntf1mx),tf1a(ntf1mx),tf2(ntf2mx),tf2a(ntf2mx),vosp0(nstmax),
     $ vosp0a(nstmax),zchar(nstmax),zchara(nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ipc,j,n,na,naa,nap,napa,nb,nba,nbi,nc,nca,nd,ne,nmu,
     $ nmua,nqt,nsa,np,npa,nrn1a,nrn2a,nr1,nr1a,nr2,nr2a,ns,nsl,nsla,
     $ ns1,ns1a,ns2,ns2a,ns3,ns3a,nx,nxa
c
c-----------------------------------------------------------------------
c
c     As long as the aqueous phase itself is not eliminated by
c     compression, keep the chloride ion from being eliminated
c     even if it is not present in the system to be modeled.
c     This allows the interpretation of NBS pH, which requires
c     the activity coefficient of the chloride ion.
c
      qchlor = .false.
      if (jpflag(iaqsla) .le. 0) then
        if (jsflag(nchlor) .eq. 2) then
          jsflag(nchlor) = 0
          qchlor = .true.
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The nsmap array is a pointer array. Here nsmap(nsa) = ns,
c     where nsa is the index of a species before compression and ns
c     is that of the same species after compression. If the species
c     is eliminated by the compression, nsmap(nsa) = 0.
c
      do nsa = 1,nsta
        nsmap(nsa) = 0
      enddo
c
      np = 0
      ns = 0
      ne = 0
      nd = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do npa = 1,npta
        if (jpflag(npa) .lt. 2) then
          np = np + 1
          jpflag(np) = jpflag(npa)
          uphase(np) = uphasa(npa)
          uptype(np) = uptypa(npa)
          nrn1a = ncmpra(1,npa)
          nrn2a = ncmpra(2,npa)
c
c         Loop over all the species belonging to this phase.
c         At least one species should have jsflag .lt. 2.
c
          ncmpr(1,np) = ns + 1
          do nsa = nrn1a,nrn2a
            if (jsflag(nsa) .lt. 2) then
              ns = ns + 1
              nsmap(nsa) = ns
              jsflag(ns) = jsflag(nsa)
              uspec(ns) = uspeca(nsa)
              mwtsp(ns) = mwtspa(nsa)
              zchar(ns) = zchara(nsa)
              vosp0(ns) = vosp0a(nsa)
              jsitex(ns) = 1
              nphasx(ns) = np
c
              nr1a = nessra(1,nsa)
              nr2a = nessra(2,nsa)
              nessr(1,ns) = ne + 1
              do n = nr1a,nr2a
                ne = ne + 1
                ness(ne) = nessa(n)
                cess(ne) = cessa(n)
              enddo
              nessr(2,ns) = ne
c
              nr1a = ndrsrd(1,nsa)
              nr2a = ndrsrd(2,nsa)
              ndrsr(1,ns) = nd + 1
              do n = nr1a,nr2a
                nd = nd + 1
                ndrs(nd) = ndrsd(n)
                cdrs(nd) = cdrsd(n)
              enddo
              ndrsr(2,ns) = nd
c
              do j = 1,ntprmx
                do i = 1,narxmx
                  axlks(i,j,ns) = axlksd(i,j,nsa)
                enddo
              enddo
c
              if (ipch .ge. 0) then
                do j = 1,ntprmx
                  do i = 1,narxmx
                    axhfs(i,j,ns) = axhfsd(i,j,nsa)
                  enddo
                enddo
                do ipc = 1,ipch
                  do j = 1,ntprmx
                    do i = 1,narxmx
                      adhfs(i,j,ipc,ns) = adhfsd(i,j,ipc,nsa)
                    enddo
                  enddo
                enddo
              endif
c
              if (ipcv .ge. 0) then
                do j = 1,ntprmx
                  do i = 1,narxmx
                    axvfs(i,j,ns) = axvfsd(i,j,nsa)
                  enddo
                enddo
                do ipc = 1,ipcv
                  do j = 1,ntprmx
                    do i = 1,narxmx
                      advfs(i,j,ipc,ns) = advfsd(i,j,ipc,nsa)
                    enddo
                  enddo
                enddo
              endif
c
            endif
          enddo
c
          ncmpr(2,np) = ns
        endif
      enddo
      npt = np
      nst = ns
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Suppress the aqueous chloride ion if it is not present.
c
      if (qchlor) then
        ns = nsmap(nchlor)
        jsflag(ns) = 2
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Update the jflag array.
c
      do nsa = 1,nsta
        ns = nsmap(nsa)
        if (ns .gt. 0) jflag(ns) = jflag(nsa)
      enddo
c
      do ns = nst + 1,nsta
        jflag(ns) = 0
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Update the species indices stored in the ndrs array.
c
      do n = 1,nd
        nsa = ndrs(n)
        if (nsa .gt. 0.) ndrs(n) = nsmap(nsa)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Update species range pointers.
c
c     Calling sequence substitutions:
c       narn1 for nlim1
c       narn2 for nlim2
c       narn1a for nlim1a
c       narn2a for nlim2a
c       nat for ntot
c
      call cnvndx(narn1,narn1a,narn2,narn2a,nsmap,nstmax,nat)
c
c     Calling sequence substitutions:
c       nmrn1 for nlim1
c       nmrn2 for nlim2
c       nmrn1a for nlim1a
c       nmrn2a for nlim2a
c       nmt for ntot
c
      call cnvndx(nmrn1,nmrn1a,nmrn2,nmrn2a,nsmap,nstmax,nmt)
c
c     Calling sequence substitutions:
c       nlrn1 for nlim1
c       nlrn2 for nlim2
c       nlrn1a for nlim1a
c       nlrn2a for nlim2a
c       nlt for ntot
c
      call cnvndx(nlrn1,nlrn1a,nlrn2,nlrn2a,nsmap,nstmax,nlt)
c
c     Calling sequence substitutions:
c       ngrn1 for nlim1
c       ngrn2 for nlim2
c       ngrn1a for nlim1a
c       ngrn2a for nlim2a
c       ngt for ntot
c
      call cnvndx(ngrn1,ngrn1a,ngrn2,ngrn2a,nsmap,nstmax,ngt)
c
c     Note- nxt is the number of solid solution phases, not the
c     number of solid solution species. Hence the use of nqt below
c     for the number of species in solid solution phases. This
c     variable (nqt) is not used elsewhere in this code.
c
c     Calling sequence substitutions:
c       nxrn1 for nlim1
c       nxrn2 for nlim2
c       nxrn1a for nlim1a
c       nxrn2a for nlim2a
c       nqt for ntot
c
      call cnvndx(nxrn1,nxrn1a,nxrn2,nxrn2a,nsmap,nstmax,nqt)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the phase index of the aqueous solution.
c
      iaqsln = nphasx(narn1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Update phase range pointers for pure liquids, pure solids,
c     and solid solutions.
c
      if (nlrn1 .le. 0) then
        ilrn1 = 0
        ilrn2 = -1
      else
        ilrn1 = nphasx(nlrn1)
        ilrn2 = nphasx(nlrn2)
      endif
c
      if (nmrn1 .le. 0) then
        imrn1 = 0
        imrn2 = -1
      else
        imrn1 = nphasx(nmrn1)
        imrn2 = nphasx(nmrn2)
      endif
c
      if (nxrn1 .le. 0) then
        nxt = 0
        ixrn1 = 0
        ixrn2 = -1
      else
        ixrn1 = nphasx(nxrn1)
        ixrn2 = nphasx(nxrn2)
        nxt = ixrn2 - ixrn1 + 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compress the basis set. Use the nbmap array as a pointer array
c     analogous to the nsmap array.
c
      nb = 0
      do nba = 1,nbtd
        nbmap(nba) = 0
        nsa = nbaspd(nba)
        ns = nsmap(nsa)
        if (ns .gt. 0) then
          nb = nb + 1
          nbmap(nba) = nb
          nbasp(nb) = ns
        endif
      enddo
      nbt = nb
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Reset the ndecsp basis pointer array.
c
      do nbi = 1,nbti
        nba = ndecsp(nbi)
        if (nba .gt. 0) ndecsp(nbi) = nbmap(nba)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compress the set of chemical elements. Use the ncmap
c     array as a pointer array analogous to the nsmap array.
c
      nc = 0
      do nca = 1,ncta
        ncmap(nca) = 0
        do nba = 1,nbtd
          nsa = nbaspd(nba)
          ns = nsmap(nsa)
          if (ns .gt. 0) then
            nr1 = nessra(1,nsa)
            nr2 = nessra(2,nsa)
            do n = nr1,nr2
              if (nessa(n) .eq. nca) then
                nc = nc + 1
                ncmap(nca) = nc
                uelem(nc) = uelema(nca)
                atwt(nc) = atwta(nca)
                go to 460
              endif
            enddo
          endif
        enddo
  460   continue
      enddo
      nct = nc
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Update the element indices stored in the ness array.
c
      do n = 1,ne
        nca = ness(n)
        if (nca .gt. 0) ness(n) = ncmap(nca)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Redefine the nbaspd array so that it matches the compressed
c     nbasp array.
c
      do nb = 1,nbtd
        nbaspd(nb) = 0
      enddo
c
      do nb = 1,nbt
        nbaspd(nb) = nbasp(nb)
      enddo
      nbtd = nbt
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Create compressed activity coefficient data for solution
c     phases other than the aqueous solution.
c
cXX   Note- the scheme below is not entirely satisfactory, as it
cXX   does not account for the effects of compression of component
cXX   species. The problem is that the apx-apxa arrays do not have
cXX   the data explcitly tied to the corresponding pairs or larger
cXX   groups of species. Revisit this issue when solid solutions
cXX   are redone.
c
      nxa = 0
      nx = 0
      do npa = ixrn1a,ixrn2a
        nxa = npa - ixrn1a + 1
        if (jpflag(npa) .lt. 2) then
          nx = nx + 1
          jsol(nx) = jsola(nxa)
c
          iapxt(nx) = iapxta(nxa)
          do i = 1,iapxta(nx)
            apx(i,nx) = apxa(i,nxa)
          enddo
c
          ibpxt(nx) = ibpxta(nxa)
          do i = 1,ibpxta(nx)
            bpx(i,nx) = bpxa(i,nxa)
          enddo
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Create compressed arrays for titration factors.
c
      ntf1t = 0
      do n = 1,ntf1ta
        nsa = ntf1a(n)
        ns = nsmap(nsa)
        if (ns .gt. 0) then
          ntf1t = ntf1t + 1
          ntf1(ntf1t) = ns
          tf1(ntf1t) = tf1a(n)
        endif
      enddo
c
      ntf2t = 0
      do n = 1,ntf2ta
        nsa = ntf2a(n)
        ns = nsmap(nsa)
        if (ns .gt. 0) then
          ntf2t = ntf2t + 1
          ntf2(ntf2t) = ns
          tf2(ntf2t) = tf2a(n)
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Create compressed arrays for hard core diameters, neutral
c     solute species flags, and hydration numbers.
c
      na = 0
      do nsa = narn1a,narn2a
        ns = nsmap(nsa)
        if (ns .gt. 0) then
          naa = nsa - narn1a + 1
          na = na + 1
          azero(na) = azeroa(naa)
          insgf(na) = insgfa(naa)
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopg(1) .eq. 1) then
c
c       Create compressed arrays for the parameters belonging to
c       Pitzer's equations.
c
        nsl = 0
        do nsla = 1,nslta
          ns1a = nslxa(1,nsla)
          ns2a = nslxa(2,nsla)
          ns1 = nsmap(ns1a)
          ns2 = nsmap(ns2a)
          n = ns1*ns2
          if (n .gt. 0) then
            nsl = nsl + 1
c
            nalpha(nsl) = nalpaa(nsla)
c
            nslx(1,nsl) = ns1
            nslx(2,nsl) = ns2
            do i = 0,ipbtmx
              do j = 1,jpfcmx
                aslm(j,i,nsl) = aslma(j,i,nsla)
              enddo
            enddo
          endif
        enddo
        nslt = nsl
c
c       Copy, but don't compress, the Pitzer alpha parameters.
c
        do napa = 1,napta
          nap = napa
          do i = 1,ipbtmx
            palpha(i,nap) = palpaa(i,napa)
          enddo
        enddo
        napt = napta
c
        nmu = 0
        do nmua = 1,nmuta
          ns1a = nmuxa(1,nmua)
          ns2a = nmuxa(2,nmua)
          ns3a = nmuxa(3,nmua)
          ns1 = nsmap(ns1a)
          ns2 = nsmap(ns2a)
          ns3 = nsmap(ns3a)
          n = ns1*ns2*ns3
          if (n .gt. 0) then
            nmu = nmu + 1
            nmux(1,nmu) = ns1
            nmux(2,nmu) = ns2
            nmux(3,nmu) = ns3
            do j = 1,jpfcmx
              amu(j,nmu) = amua(j,nmua)
            enddo
          endif
        enddo
        nmut = nmu
      endif
c
c
      end
