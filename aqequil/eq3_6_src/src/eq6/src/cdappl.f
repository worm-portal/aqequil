      subroutine cdappl(acflg,acfw,acfwlg,actlg,actw,actwlg,adwipp,
     $ afcnst,affpd,affsd,ah,ahrc,alk,alk1,alk2,alki,atwt,bdwipp,
     $ cdrsd,cess,conc,csts,ctb,cteaq,dvoso,dwoso,eh,ehfac,ehrc,
     $ eps100,farad,fdpe0,fdse0,fjest,fo2lg,fo2lrc,fxist,iaqsln,
     $ iemop0,iemos0,iern1,iern2,ifrn1,ifrn2,ilrn1,ilrn2,imrn1,
     $ imrn2,iopt,ixrn1,ixrn2,jcode,jcsort,jern1,jflag,jflagd,
     $ jgext,jpflag,jsflag,jssort,modr,moph,mophg,mophj,mopht,morr,
     $ mosp,mospg,mospj,mospt,mprph,mprsp,mrgers,mrmlra,mtb,
     $ mtbaq,mte,mteaq,mwtges,mwtrc,mwtsp,narn1,narn2,nat,nbasp,
     $ nbaspd,nbt,nchlor,ncmpe0,ncmpr,nct,ndrsd,ndrsrd,nelect,
     $ nern1,nern2,nert,ness,nessr,net,nfrn1,nfrn2,ngext,ngrn1,
     $ ngrn2,ngt,nhydr,nlrn1,nlrn2,nlt,nmrn1,nmrn2,nmrt,nmt,no2gaq,
     $ npchk,npet,npet0,npt,npts,nrct,nrndex,nst,nsts,nstsr,ntf1,
     $ ntf1t,ntf2,ntf2t,nxridx,nxrn1,nxrn2,nxrt,nxt,osc,oscst,
     $ omega,pe,perc,ph,phmes,ppmwb,ppmwe,qriinf,qxknph,rreacn,
     $ rreac1,rxbar,sfcar,sidrsp,sidrph,sigmam,sigmst,tdays,tempc,
     $ tf1,tf2,thours,time1,tmins,tyears,uphase,uspec,vodrt,voph,
     $ vophg,vophj,vopht,vosoct,vosp,vospg,vospj,vospt,vosp0,
     $ vreac,wfh2o,wkgwi,wodr,wodrt,woph,wophg,wophj,wopht,
     $ worr,worrt,wosoct,wosp,wospg,wospj,wospt,xbar,xlke,
     $ xlksd,zchcu6,zchsq2)
c
c     This subroutine computes various data pertaining to the modeled
c     system for printing and plotting at the current point of
c     reaction progress (xi1).
c
c     This subroutine is called by:
c
c       EQ6/path.f
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
      include 'eqlib/eqlpar.h'
      include 'eqlib/eqldv.h'
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer iemop0(npetmx),iemos0(nsetmx),iopt(noptmx),
     $ jcode(nrctmx),jcsort(nstmax),jern1(jetmax,netmax),
     $ jflag(nstmax),jflagd(nstmax),jgext(netmax),jpflag(nptmax),
     $ jsflag(nstmax),jssort(nstmax),nbasp(nbtmax),nbaspd(nbtmax),
     $ ncmpe0(2,npetmx),ncmpr(2,nptmax),ndrsd(ndrsmx),ndrsrd(2,nstmax),
     $ ness(nessmx),nessr(2,nstmax),ngext(jetmax,netmax),npchk(nptmax),
     $ nrndex(nrctmx),nsts(nstsmx),nstsr(2,nstmax),ntf1(ntf1mx),
     $ ntf2(ntf2mx),nxridx(nrctmx)
c
      integer iern1,iern2,ifrn1,ifrn2,ilrn1,ilrn2,imrn1,imrn2,
     $ ixrn1,ixrn2
c
      integer nat,nbt,nct,net,nlt,nmt,nst,npt,ngt,nxt
c
      integer narn1,narn2,nern1,nern2,nfrn1,nfrn2,ngrn1,ngrn2,
     $ nlrn1,nlrn2,nmrn1,nmrn2,nxrn1,nxrn2
c
      integer iaqsln,nchlor,nelect,nert,nhydr,nmrt,no2gaq,npet,npet0,
     $ npts,nrct,ntf1t,ntf2t,nxrt
c
      logical qxknph(nptmax)
c
      logical qriinf
c
      character(len=48) uspec(nstmax)
      character(len=24) uphase(nptmax)
c
      real*8 acflg(nstmax),actlg(nstmax),affpd(nptmax),affsd(nstmax),
     $ ahrc(nbtmax),atwt(nctmax),cdrsd(ndrsmx),cess(nessmx),
     $ conc(nstmax),csts(nstsmx),ctb(nbtmax),cteaq(nctmax),ehrc(nbtmax),
     $ fdpe0(nordmx,npetmx),fdse0(nordmx,nsetmx),fo2lrc(nbtmax)
c
      real*8 modr(nrctmx),moph(nptmax),mophg(nptmax),mophj(nptmax),
     $ mopht(nptmax),morr(nrctmx),mosp(nstmax),mospg(nstmax),
     $ mospj(nstmax),mospt(nstmax),mprph(nptmax),mprsp(nstmax),
     $ mrgers(ietmax,jetmax,nertmx),mtb(nbtmax),mtbaq(nbtmax),
     $ mte(nctmax),mteaq(nctmax),mwtges(netmax),mwtrc(nrctmx),
     $ mwtsp(nstmax)
c
      real*8 perc(nbtmax),ppmwb(nbtmax),ppmwe(nctmax),
     $ rreacn(nrctmx),rreac1(nrctmx),rxbar(iktmax,nxrtmx),
     $ sfcar(nrctmx),sidrph(nptmax),sidrsp(nstmax),tf1(ntf1mx),
     $ tf2(ntf2mx),voph(nptmax),vophg(nptmax),vophj(nptmax),
     $ vopht(nptmax),vosp(nstmax),vospg(nstmax),vospj(nstmax),
     $ vospt(nstmax),vosp0(nstmax),vreac(nrctmx),wodr(nrctmx),
     $ woph(nptmax),wophg(nptmax),wophj(nptmax),wopht(nptmax),
     $ worr(nrctmx),wosp(nstmax),wospg(nstmax),wospj(nstmax),
     $ wospt(nstmax),xbar(nstmax),xlksd(nstmax),zchcu6(nstmax),
     $ zchsq2(nstmax)
c
      real*8 acfw,acfwlg,actw,actwlg,adwipp,afcnst,ah,alk,alki,
     $ alk1,alk2,bdwipp,dvoso,dwoso,eh,ehfac,eps100,farad,fjest,
     $ fo2lg,fxist,mrmlra,osc,oscst,omega,pe,ph,phmes,sigmam,
     $ sigmst,tdays,tempc,thours,time1,tmins,tyears,vosoct,vodrt,
     $ wfh2o,wkgwi,wosoct,wodrt,worrt,xlke
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ie,je,ik,n,nb,nc,ne,ner,np,npe,np1,ns,nse,nrc,nr1,nr2,
     $ nss,nxr
c
      logical qrct
c
      real*8 mophx,mospx,osfac,vophx,vospx,wophx,wospx
c
      real*8 texp,tlg
c
c-----------------------------------------------------------------------
c
c     Compute the model time in various units.
c
      if (iopt(2) .ge. 1) then
        if (qriinf) then
          tmins = 1.e+38
          thours = 1.e+38
          tdays = 1.e+38
          tyears = 1.e+38
        else
          tmins = time1/60.
          thours = tmins/60.
          tdays = thours/24.
          tyears = tdays/365.25
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute some data for reactants and corresponding reactions.
c
      vodrt = 0.
      wodrt = 0.
      worrt = 0.
c
      do nrc = 1,nrct
        worr(nrc) = mwtrc(nrc)*morr(nrc)
        wodr(nrc) = mwtrc(nrc)*modr(nrc)
        worrt = worrt + worr(nrc)
        wodrt = wodrt + wodr(nrc)
        vodrt = vodrt + modr(nrc)*vreac(nrc)
      enddo
c
      if (iopt(2) .ge. 1) then
        do nrc = 1,nrct
          rreacn(nrc) = 0.
          if (sfcar(nrc) .gt. 0.) rreacn(nrc) = rreac1(nrc)/sfcar(nrc)
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute data for a summary of the elemental composition of the
c     aqueous phase.
c
      call initaz(mte,nctmax)
      call initaz(mteaq,nctmax)
c
      do nss = 1,nst
        ns = jssort(nss)
        if (mosp(ns) .gt. 0.) then
          nr1 = nessr(1,ns)
          nr2 = nessr(2,ns)
          do n = nr1,nr2
            nc = ness(n)
            if (nc .gt. 0) mte(nc) = mte(nc) + cess(n)*mosp(ns)
          enddo
        endif
      enddo
c
      do nss = narn1,narn2
        ns = jcsort(nss)
        if (ns .ne. no2gaq) then
          nr1 = nessr(1,ns)
          nr2 = nessr(2,ns)
          do n = nr1,nr2
            nc = ness(n)
            if (nc .gt. 0) mteaq(nc) = mteaq(nc) + cess(n)*mosp(ns)
          enddo
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the elemental composition of the aqueous solution
c     (molalities and ppm: mg/kg.sol).
c
      call initaz(cteaq,nctmax)
      call initaz(ppmwe,nctmax)
c
      do nc = 1,nct
        cteaq(nc) = wkgwi*mteaq(nc)
        ppmwe(nc) = 1000.*cteaq(nc)*atwt(nc)*wfh2o
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the mass and concentration totals in the aqueous phase
c     in terms of data file basis species.
c
      call initaz(mtbaq,nbtmax)
      call initaz(ctb,nbtmax)
      call initaz(ppmwb,nbtmax)
c
      do nss = narn1,narn2
        ns = jcsort(nss)
        if (mosp(ns) .ne. 0.) then
          nr1 = nstsr(1,ns)
          nr2 = nstsr(2,ns)
          do n = nr1,nr2
            nb = nsts(n)
            mtbaq(nb) = mtbaq(nb) + csts(n)*mosp(ns)
          enddo
        endif
      enddo
c
      do nb = 1,nbt
        ctb(nb) = wkgwi*mtbaq(nb)
        ns = nbaspd(nb)
        ppmwb(nb) = 1000.*ctb(nb)*mwtsp(ns)*wfh2o
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute various aqueous phase parameters.
c
      acfwlg = acflg(narn1)
      acfw = texp(acfwlg)
c
c     Compute the sum of the stoichiometric molalities.
c
      sigmst = 0.
      do nb = 1,nbt
        ns = nbaspd(nb)
c
c       The 'abs' comes into play below only when ns = nhydr or nhydx
c       and the value is negative (e.g., a negative ctb for H+ is
c       interpreted as a positive value for OH-).
c
        if (ns.ne.no2gaq .and. ns.ne.narn1)
     $  sigmst = sigmst + abs(ctb(nb))
      enddo
c
c     Compute the true and stoichiometric osmotic coefficients.
c
      osfac = -omega * log(actw)
      osc = 0.
      if (sigmam .gt. eps100) osc = osfac/sigmam
      oscst = 0.
      if (sigmst .gt. eps100) oscst = osfac/sigmst
c
c     Compute the stoichiometric ionic strength.
c
c     Calling sequence substitution:
c       fxist for fxistc
c
      call cfxist(ctb,fxist,nbaspd,nbt,nbtmax,nstmax,zchsq2)
c
c     Compute the stoichiometric ionic asymmetry.
c
c     Calling sequence substitution:
c       fjest for fjestc
c
      call cfjest(ctb,fjest,nbaspd,nbt,nbtmax,nstmax,zchcu6)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute alkalinity parameters.
c
c     Compute the HCO3-CO3-OH total alkalinity.
c
c     Calling sequence substitutions:
c       alk1 for alkc
c       ntf1 for ntfx
c       ntf1mx for ntfxmx
c       ntf1t for ntfxt
c       tf1 for tfx
c
      call calk(alk1,conc,nstmax,ntf1,ntf1mx,ntf1t,tf1)
c
c     Compute the extended total alkalinity.
c
c     Calling sequence substitutions:
c       alk2 for alkc
c       ntf2 for ntfx
c       ntf2mx for ntfxmx
c       ntf2t for ntfxt
c       tf2 for tfx
c
      call calk(alk2,conc,nstmax,ntf2,ntf2mx,ntf2t,tf2)
c
      alki = 0.
      alk = alk2
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute data for redox reactions that are not constrained to be
c     at equilibrium.
c
      call cdardx(actlg,actwlg,ah,ahrc,cdrsd,eh,ehfac,ehrc,
     $ farad,fo2lg,fo2lrc,jsflag,mosp,nbasp,nbaspd,nbt,nbtmax,ndrsd,
     $ ndrsmx,ndrsrd,no2gaq,nstmax,pe,perc,ph,xlke,xlksd)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate affinities and saturation indices using the 'd' set
c     of reactions.
c
      call gaffsd(actlg,afcnst,affpd,affsd,cdrsd,jflagd,jpflag,
     $ ncmpr,ndrsd,ndrsmx,ndrsrd,npt,nptmax,nst,nstmax,qxknph,sidrph,
     $ sidrsp,uphase,uspec,xbar,xlksd)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      vosoct = 0.
      wosoct = 0.
c
      call initaz(voph,nptmax)
      call initaz(woph,nptmax)
      call initaz(vosp,nstmax)
      call initaz(wosp,nstmax)
c
      call initaz(mopht,nptmax)
      call initaz(vopht,nptmax)
      call initaz(wopht,nptmax)
      call initaz(mospt,nstmax)
      call initaz(vospt,nstmax)
      call initaz(wospt,nstmax)
c
      call initaz(mophg,nptmax)
      call initaz(vophg,nptmax)
      call initaz(wophg,nptmax)
      call initaz(mospg,nstmax)
      call initaz(vospg,nstmax)
      call initaz(wospg,nstmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(1).eq.0 .or. iopt(1).eq.1) then
c
c       Compute data which describe the solid phases in the ES.
c
        if (npet .gt. 1) then
          do np = 1,npt
            if (np .ne. iaqsln) then
              if (npchk(np).eq.0 .and. jpflag(np).eq.-1) then
                nr1 = ncmpr(1,np)
                nr2 = ncmpr(2,np)
                vophx = 0.
                wophx = 0.
                do ns = nr1,nr2
                  wospx = mosp(ns)*mwtsp(ns)
                  wosp(ns) = wospx
                  wophx = wophx + wospx
                  vospx = mosp(ns)*vosp0(ns)
                  vosp(ns) = vospx
                  vophx = vophx + vospx
                enddo
                if (np.ge.iern1 .and. np.le.iern2) then
                  ne = np - iern1 + 1
                  wophx = wophx + moph(np)*mwtges(ne)
                endif
                woph(np) = wophx
                voph(np) = vophx
                wosoct = wosoct + wophx
                vosoct = vosoct + vophx
              endif
            endif
          enddo
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(1) .eq. 2) then
c
c       Compute the product mineral assemblage which is instantaneous
c       with respect to Xi (for the fluid-centered flow-through system
c       only).
c
        call initaz(mophj,nptmax)
        call initaz(vophj,nptmax)
        call initaz(wophj,nptmax)
        call initaz(mospj,nstmax)
        call initaz(vospj,nstmax)
        call initaz(wospj,nstmax)
c
        if (npts .gt. 0) then
c
c         Use the two-point derivative.
c
          do npe = 1,npet0
            np = iemop0(npe)
            mophj(np) = fdpe0(1,npe)
            if (mophj(np) .lt. 0.) mophj(np) = 0.
            wophx = 0.
            vophx = 0.
c
            nr1 = ncmpe0(1,npe)
            nr2 = ncmpe0(2,npe)
            do nse = nr1,nr2
              ns = iemos0(nse)
              mospx = fdse0(1,nse)
              mospj(ns) = mospx
              wospx = mospx*mwtsp(ns)
              wospj(ns) = wospx
              wophx = wophx + wospx
              vospx = mospx*vosp0(ns)
              vospj(ns) = vospx
              vophx = vophx + vospx
            enddo
            if (np.ge.iern1 .and. np.le.iern2) then
              ne = np - iern1 + 1
              wophx = wophx + mophj(np)*mwtges(ne)
            endif
            wophj(np) = wophx
            vophj(np) = vophx
          enddo
        endif
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(1) .eq. 2) then
c
c       Compute data which describe the solid phases in the system
c       ES + PRS.
c
        do np = 1,npt
          if (np .ne. iaqsln) then
            if (jpflag(np) .lt. 2) then
              mopht(np) = moph(np) + mprph(np)
              nr1 = ncmpr(1,np)
              nr2 = ncmpr(2,np)
              vophx = 0.
              wophx = 0.
              do ns = nr1,nr2
                mospx = mosp(ns) + mprsp(ns)
                mospt(ns) = mospx
                wospx = mospx*mwtsp(ns)
                wospt(ns) = wospx
                wophx = wophx + wospx
                vospx = mospx*vosp0(ns)
                vospt(ns) = vospx
                vophx = vophx + vospx
              enddo
              if (np.ge.iern1 .and. np.le.iern2) then
                ne = np - iern1 + 1
                wophx = wophx + mopht(np)*mwtges(ne)
              endif
              wopht(np) = wophx
              vopht(np) = vophx
              wosoct = wosoct + wophx
              vosoct = vosoct + vophx
            endif
          endif
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if ((nmrt + nxrt + nert) .gt. 0) then
c
c       Print a grand summary of the solid phases in the system
c       ES + PRS + reactants.
c
        do np = imrn1,imrn2
          mophx = moph(np) + mprph(np)
          qrct = .false.
          if (nmrt .gt. 0) then
            do nrc = 1,nrct
              if (jcode(nrc) .eq. 0) then
                np1 = nrndex(nrc)
                if (np1 .eq. np) then
                  mophx = mophx + morr(nrc)
                  qrct = .true.
                  go to 607
                endif
              endif
            enddo
          endif
c
  607     if (mophx.ne.0. .or. qrct) then
            mophg(np) = mophx
c
            ns = ncmpr(1,np)
            mospx = mophx
            mospg(ns) = mospx
            wospx = mospx*mwtsp(ns)
            wospg(ns) = wospx
            vospx = mospx*vosp0(ns)
            vospg(ns) = vospx
c
            wophx = wospx
            vophx = vospx
            wophg(np) = wophx
            vophg(np) = vophx
          endif
        enddo
c
        do np = ifrn1,ifrn2
          mophx = moph(np) + mprph(np)
          qrct = .false.
          if (nmrt .gt. 0) then
            do nrc = 1,nrct
              if (jcode(nrc) .eq. 0) then
                np1 = nrndex(nrc)
                if (np1 .eq. np) then
                  mophx = mophx + morr(nrc)
                  qrct = .true.
                  go to 609
                endif
              endif
            enddo
          endif
c
  609     if (mophx.ne.0. .or. qrct) then
            mophg(np) = mophx
c
            ns = ncmpr(1,np)
            mospx = mophx
            mospg(ns) = mospx
            wospx = mospx*mwtsp(ns)
            wospg(ns) = wospx
            vospx = mospx*vosp0(ns)
            vospg(ns) = vospx
c
            wophx = wospx
            vophx = vospx
            wophg(np) = wophx
            vophg(np) = vophx
          endif
        enddo
c
        if (iopt(4) .gt. 0) then
          do np = ixrn1,ixrn2
            if (jpflag(np) .lt. 2) then
              mophx = moph(np) + mprph(np)
              nr1 = ncmpr(1,np)
              nr2 = ncmpr(2,np)
              qrct = .false.
              nrc = 0
c
              if (nxrt .gt. 0) then
                do nrc = 1,nrct
                  if (jcode(nrc) .eq. 1) then
                    np1 = nrndex(nrc)
                    if (np1 .eq. np) then
                      mophx = mophx + morr(nrc)
                      qrct = .true.
                      go to 612
                    endif
                  endif
                enddo
              endif
c
  612         if (mophx.ne.0. .or. qrct) then
                mophg(np) = mophx
                if (qrct) nxr = nxridx(nrc)
                ik = 0
                vophx = 0.
                wophx = 0.
                do ns = nr1,nr2
                  mospx = mosp(ns) + mprsp(ns)
                  ik = ik + 1
                  if (qrct) mospx = mospx + morr(nrc)*rxbar(ik,nxr)
                  mospg(ns) = mospx
                  wospx = mospx*mwtsp(ns)
                  wospg(ns) = wospx
                  wophx = wophx + wospx
                  vospx = mospx*vosp0(ns)
                  vospg(ns) = vospx
                  vophx = vophx + vospx
                enddo
                vophg(np) = vophx
                wophg(np) = wophx
              endif
            endif
          enddo
        endif
c
        do np = iern1,iern2
          if (jpflag(np) .lt. 2) then
            mophx = moph(np) + mprph(np)
            qrct = .false.
            nrc = 0
c
            if (nert .gt. 0) then
              do nrc = 1,nrct
                if (jcode(nrc) .eq. 5) then
                  np1 = nrndex(nrc)
                  if (np1 .eq. np) then
                    mophx = mophx + morr(nrc)
                    qrct = .true.
                    go to 614
                  endif
                endif
              enddo
            endif
c
  614       if (mophx.ne.0. .or. qrct) then
              mophg(np) = mophx
              if (qrct) ner = nxridx(nrc)
              vophx = 0.
              wophx = 0.
              ne = np - iern1 + 1
              do je = 1,jgext(ne)
                ns = jern1(je,ne) - 1
                do ie = 1,ngext(je,ne)
                  ns = ns + 1
                  mospx = mosp(ns) + mprsp(ns)
                  if (qrct) mospx = mospx + morr(nrc)*mrgers(ie,je,ner)
                  mospg(ns) = mospx
                  wospx = mospx*mwtsp(ns)
                  wospg(ns) = wospx
                  wophx = wophx + wospx
                  vospx = mospx*vosp0(ns)
                  vospg(ns) = vospx
                  vophx = vophx + vospx
                enddo
              enddo
              wophx = wophx + mophg(np)*mwtges(ne)
c
              vophg(np) = vophx
              wophg(np) = wophx
            endif
          endif
        enddo
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      dvoso = vosoct - vodrt
      dwoso = wosoct - wodrt
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
