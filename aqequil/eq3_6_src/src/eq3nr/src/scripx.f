      subroutine scripx(abar,acflg,act,actlg,adh,afcnst,affpd,
     $ affsd,ahrc,alki,apx,atwt,a3bar,a3bars,bpx,cdrsd,cegexs,
     $ cess,conc,conclg,coval,csts,ctb,cteaq,egexjc,egexjf,egexpa,
     $ egexpc,egexs,egexw,ehfac,ehrc,eps100,eh,farad,fje,fo2,fo2lg,
     $ fo2lrc,fugac,fugalg,fxi,iapxmx,ibpxmx,iebal,iern1,iern2,ietmax,
     $ igas,iktmax,iopg,iopr,iopt,ilrn1,ilrn2,imrn1,imrn2,ixrn1,
     $ ixrn2,jcsort,jern1,jern2,jetmax,jflag,jflagd,jflgi,jfleba,jgext,
     $ jgsort,jpflag,jsflag,jsol,jsomax,kern1,kern2,ketmax,kgexsa,
     $ mlmrra,mrmlra,moph,mosp,mte,mteaq,mwtsp,narn1,narn2,natmax,
     $ nbasp,nbaspd,nbt,nbtmax,nchlor,ncmpr,nct,nctmax,ndrsd,ndrsmx,
     $ ndrsrd,nelect,ness,nessmx,nessr,net,neti,netmax,ngexpi,ngexsa,
     $ ngext,ngrn1,ngrn2,ngt,ngtmax,nhydr,nhydx,nopgmx,noprmx,noptmx,
     $ noutpt,no2gaq,npnxp,npt,nptmax,nrdxsp,nst,nstmax,nsts,nstsmx,
     $ nstsr,ntf1,ntf1mx,ntf1t,ntf2,ntf2mx,ntf2t,nttyo,nxrn1,nxrn2,
     $ nxt,nxti,nxtimx,nxtmax,omega,pe,perc,ppmwe,qrho,qxknph,rho,
     $ rhoc,rhowc,sidrph,sidrsp,sigmam,sigzi,tdsglw,tdspkc,tdspkg,
     $ tdspl,tdsplc,tempc,tf1,tf2,tolspf,uelem,ugexj,ugexmo,uphase,
     $ uspec,uxtype,vosol,wfac,wfh2o,wftds,wkgwi,woh2o,wosol,wotds,
     $ xbar,xbarlg,xbarw,xbrwlg,xgexw,xlke,xlksd,zchar,zchcu6,zchsq2)
c
c     This subroutine writes a description of the computed aqueous
c     solution model on the output file.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
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
      integer iapxmx,ibpxmx,ietmax,iktmax,jetmax,jsomax,ketmax,natmax,
     $ nbtmax,nctmax,ndrsmx,nessmx,netmax,ngtmax,nopgmx,noprmx,noptmx,
     $ nptmax,nstmax,nstsmx,ntf1mx,ntf2mx,nxtimx,nxtmax
c
      integer noutpt,nttyo
c
      integer iopg(nopgmx),iopr(noprmx),iopt(noptmx),jcsort(nstmax),
     $ jern1(jetmax,netmax),jern2(jetmax,netmax),jgsort(ngtmax),
     $ jflag(nstmax),jflagd(nstmax),jflgi(nbtmax),jgext(netmax),
     $ jpflag(nptmax),jsflag(nstmax),jsol(nxtmax),kern1(netmax),
     $ kern2(netmax),kgexsa(ketmax,netmax),nbasp(nbtmax),nbaspd(nbtmax),
     $ ncmpr(2,nptmax),ndrsd(ndrsmx),ndrsrd(2,nstmax),ness(nessmx),
     $ nessr(2,nstmax),ngexpi(netmax),ngexsa(ietmax,jetmax,netmax),
     $ ngext(jetmax,netmax),npnxp(nxtimx),nsts(nstsmx),nstsr(2,nstmax),
     $ ntf1(ntf1mx),ntf2(ntf2mx)
c
      integer iebal,iern1,iern2,igas,ilrn1,ilrn2,imrn1,imrn2,ixrn1,
     $ ixrn2,jfleba,narn1,narn2,nbt,nchlor,nct,nelect,net,neti,ngrn1,
     $ ngrn2,ngt,nhydr,nhydx,no2gaq,npt,nrdxsp,nst,ntf1t,ntf2t,nxrn1,
     $ nxrn2,nxt,nxti
c
      logical qxknph(nptmax)
c
      logical qrho
c
      character(len=48) uspec(nstmax)
      character(len=32) uxtype(jsomax)
      character(len=24) ugexmo(netmax),uphase(nptmax)
      character(len=8) uelem(nctmax),ugexj(jetmax,netmax)
c
      real(8) acflg(nstmax),act(nstmax),actlg(nstmax),ahrc(nbtmax),
     $ atwt(nctmax),a3bars(natmax),affpd(nptmax),affsd(nstmax),
     $ apx(iapxmx,nxtmax),bpx(ibpxmx,nxtmax),cdrsd(ndrsmx),
     $ cegexs(ietmax,jetmax,netmax),cess(nessmx),conc(nstmax),
     $ conclg(nstmax),coval(nbtmax),csts(nstsmx),ctb(nbtmax),
     $ cteaq(nctmax),egexjc(jetmax,netmax),egexjf(jetmax,netmax),
     $ egexpa(netmax),egexpc(netmax),egexs(ietmax,jetmax,netmax),
     $ egexw(ketmax,netmax),ehrc(nbtmax),fo2lrc(nbtmax),
     $ fugac(ngtmax),fugalg(ngtmax)
c
      real(8) moph(nptmax),mosp(nstmax),mte(nctmax),mteaq(nctmax),
     $ mwtsp(nstmax),perc(nbtmax),ppmwe(nctmax),sidrph(nptmax),
     $ sidrsp(nstmax),tf1(ntf1mx),tf2(ntf2mx),wfac(iktmax,nxtmax),
     $ xbar(nstmax),xbarlg(nstmax),xgexw(ketmax,netmax),xlksd(nstmax),
     $ zchar(nstmax),zchcu6(nstmax),zchsq2(nstmax)
c
      real(8) abar,adh,afcnst,alki,a3bar,ehfac,eps100,eh,farad,
     $ fje,fo2,fo2lg,fxi,omega,pe,rho,sigmam,sigzi,tdspkg,tdspl,
     $ tempc,tolspf,vosol,wfh2o,wftds,wkgwi,woh2o,wosol,wotds,
     $ xbarw,xbrwlg,xlke
c
      real(8) mlmrra,mrmlra,rhoc,rhowc,tdsglw,tdspkc,tdsplc
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ier,j2,nb,ncount,nei,np,nr1,nr2,ns,nss,ns1,ns2,nt,nx,nxi
c
      integer ilnobl
c
      logical qphcl,qredox
c
      character(len=8) uadj,ufinal,uinput,ux1,ux2,ux3
c
      real(8) acfw,acfwlg,alk1,alk2,actw,actwlg,ah,ahmes,ahnbs,axc,axd,
     $ axx,cte,ctebal,cx,cxv,cxw,ehmes,ehnbs,fjest,fxist,msigzm,osc,
     $ oscst,osfac,pch,pemes,penbs,ph,phc,phcl,phd,phnbs,phmes,phx,
     $ ppmv,ppmw,pxc,pxd,pxx,pzmean,pztot,sanion,sigmst,sigza,sigzc,
     $ sigzm,sx
c
      real(8) coefst,texp,tlg
c
c-----------------------------------------------------------------------
c
      data uinput /'Input   '/,ufinal /'Final   '/,uadj   /'Adj     '/
c
c-----------------------------------------------------------------------
c
      write (noutpt,1000)
 1000 format(/' - - - - - - - - - - - - - - - - - - - - - - - - - - -',
     $ ' - - - - - - - - - - - -',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
cXX   Begin jsol values with 0 or 1? Need to work something out with
cXX   the data base group as to which. Wait until solid solutions are
cXX   redone.
c
      do nx = 1,nxt
        if (jsol(nx) .eq. 0) jsol(nx) = 1
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print a table of the elemental composition of the aqueous
c     solution.
c
      call prteca(cteaq,mrmlra,nct,nctmax,noutpt,ppmwe,
     $ qrho,rho,uelem)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute and save the calculated total concentration of the
c     species adjusted for electrical balance.
c
      if (iebal .gt. 0) then
        ctebal = 0.
        nb = iebal
        do nss = narn1,narn2
          ns = jcsort(nss)
          if (jsflag(ns) .le. 0) then
            cx = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
            ctebal = ctebal + cx*conc(ns)
          endif
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print a table of the numerical composition of the aqueous
c     solution.
c
      call prtnca(ctb,jflag,jsflag,mrmlra,mwtsp,narn1,narn2,
     $ nbasp,nbaspd,nbt,nbtmax,noutpt,nstmax,qrho,rho,uspec,wfh2o)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute and print a table of the sensible composition of the
c     aqueous solution.
c
      call prtsca(ctb,jflag,jsflag,mrmlra,mwtsp,narn1,narn2,
     $ nbasp,nbaspd,nbt,nbtmax,nelect,nhydr,nhydx,noutpt,no2gaq,
     $ nstmax,qrho,rho,uspec,wfh2o)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute and print various aqueous solution parameters.
c
      actwlg = actlg(narn1)
      actw = texp(actwlg)
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
      call prtvpa(abar,acfw,acfwlg,actw,actwlg,a3bar,fje,
     $ fjest,fo2,fo2lg,fxi,fxist,iopg,mlmrra,mrmlra,nopgmx,noutpt,
     $ osc,oscst,qrho,rhoc,rhowc,sigmam,sigmst,tdsglw,tdspkc,
     $ tdsplc,vosol,wfh2o,wftds,woh2o,wosol,wotds,xbarw,xbrwlg)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute and print pH, Eh, and pe-, all with reference to
c     appropriate pH scales. Also compute and print the pHCl.
c
      qredox = .true.
c
      call gpheh(acflg,actlg,actwlg,adh,ah,ahmes,ahnbs,conc,
     $ eh,ehfac,ehmes,ehnbs,farad,fo2lg,fxi,iopg,mrmlra,nchlor,nhydr,
     $ nopgmx,noutpt,nstmax,nttyo,pch,pe,pemes,penbs,ph,phcl,phmes,
     $ phnbs,qphcl,qredox,qrho,xlke)
c
      call prpheh(ah,ahmes,ahnbs,eh,ehmes,ehnbs,iopg,
     $ nopgmx,noutpt,pch,pe,pemes,penbs,ph,phcl,phmes,phnbs,
     $ qphcl,qredox,qrho)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print a table of computed alkalinity parameters.
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
      call prtalk(alki,alk1,alk2,mrmlra,noutpt,ntf1t,ntf2t,
     $ qrho,rho,tempc,wfh2o)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate and print the electrical balance and the cation and
c     anion contributions.
c
      call gszm(conc,jcsort,narn1,narn2,nstmax,sigza,sigzc,
     $ sigzi,sigzm,zchar)
c
      msigzm = 0.5*sigzm
      sanion = -sigza
c
      write (noutpt,1250)
 1250 format(/11x,'--- Electrical Balance Totals ---',//
     $ 34x,'eq/kg.H2O',/)
      write (noutpt,1260) sigzc,sanion,sigzm,msigzm,sigzi
 1260 format(8x,'Sigma(mz) cations= ',3x,1pe17.10,
     $ /9x,'Sigma(mz) anions= ',3x,1pe17.10,
     $ /13x,'Total charge= ',3x,1pe17.10,
     $ /14x,'Mean charge= ',3x,1pe17.10,
     $ /9x,'Charge imbalance= ',3x,1pe17.10)
      sx = 100.*sigzi
      pztot = sx/sigzm
      pzmean = sx/msigzm
      write (noutpt,1270) pztot,pzmean
 1270 format(//9x,'The electrical imbalance is:',
     $ //11x,f9.4,' per cent of the total charge',
     $ /11x,f9.4,' per cent of the mean charge',/)
c
c     Print results of electrical balancing.
c
      if (iebal .gt. 0) then
        ns = nbaspd(iebal)
        j2 = ilnobl(uspec(ns)(1:24))
        write (noutpt,1280) uspec(ns)(1:j2)
 1280   format(/6x,'--- Electrical Balancing on ',a,' ---')
        ux1 = uinput
        ux2 = ufinal
        ux3 = uadj
        if (jfleba .eq. 0) then
          write (noutpt,1290)
 1290     format(/17x,'mg/L',10x,'mg/kg.sol',7x,'Molality',/)
          cte = coval(iebal)
          cxw = 1000.*cte*mwtsp(ns)*wfh2o
          cxv = cxw*rho
          write (noutpt,1300) ux1,cxv,cxw,cte
 1300     format(5x,a6,f13.4,3x,f13.4,3x,1pe17.10)
          cte = ctebal
          cxw = 1000.*cte*mwtsp(ns)*wfh2o
          cxv = cxw*rho
          write (noutpt,1300) ux2,cxv,cxw,cte
          cte = ctebal - coval(iebal)
          cxw = 1000.*cte*mwtsp(ns)*wfh2o
          cxv = cxw*rho
          write (noutpt,1300) ux3,cxv,cxw,cte
        elseif (jfleba .eq. 16) then
          write (noutpt,1310)
 1310     format(/3x,'        Log Activity',/)
          axx = coval(iebal)
          axc = actlg(ns)
          write (noutpt,1320) ux1,axx
          write (noutpt,1320) ux2,axc
 1320     format(5x,a6,f13.4)
          axd = axc - axx
          write (noutpt,1320) ux3,axd
        elseif (jfleba .eq. 19) then
          write (noutpt,1312)
 1312     format(/3x,'            pX',/)
          pxx = coval(iebal)
          pxc = -actlg(ns)
          write (noutpt,1320) ux1,pxx
          write (noutpt,1320) ux2,pxc
          pxd = pxc - pxx
          write (noutpt,1320) ux3,pxd
        elseif (jfleba .eq. 20) then
          write (noutpt,1314)
 1314     format(/3x,'            pH',/)
          phx = coval(iebal)
          phc = -actlg(ns)
          write (noutpt,1320) ux1,phx
          write (noutpt,1320) ux2,phc
          phd = phc - phx
          write (noutpt,1320) ux3,phd
        elseif (jfleba .eq. 21) then
          write (noutpt,1315)
 1315     format(/3x,'            pHCl',/)
          phx = coval(iebal)
          phc = -actlg(nhydr)-actlg(nchlor)
          write (noutpt,1320) ux1,phx
          write (noutpt,1320) ux2,phc
          phd = phc - phx
          write (noutpt,1320) ux3,phd
        elseif (jfleba .eq. 22) then
          write (noutpt,1316)
 1316     format(/3x,'            pmH',/)
          phx = coval(iebal)
          phc = -conclg(ns)
          write (noutpt,1320) ux1,phx
          write (noutpt,1320) ux2,phc
          phd = phc - phx
          write (noutpt,1320) ux3,phd
        elseif (jfleba .eq. 23) then
          write (noutpt,1318)
 1318     format(/3x,'            pxH',/)
          pxx = coval(iebal)
          pxc = -conclg(ns)
          write (noutpt,1320) ux1,pxx
          write (noutpt,1320) ux2,pxc
          pxd = pxc - pxx
          write (noutpt,1320) ux3,pxd
        endif
        write (noutpt,1322)
 1322   format(1x)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopr(5) .gt. 0) then
c
c       Compute and print activity ratios of aqueous species.
c
        call prtacr(actlg,iopr,jsflag,nbaspd,nbt,nbtmax,nelect,
     $  nhydr,noprmx,no2gaq,noutpt,nstmax,uspec,zchar)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print the aqueous species distribution.
c
      call prtaqs(acflg,actlg,conc,conclg,iopr,jcsort,narn1,
     $ narn2,noprmx,noutpt,nstmax,uspec)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopr(9) .gt. 0) then
c
c       Compute and print the mean ionic activities and activity
c       coefficients.
c
        call prtmip(acflg,actlg,conclg,ctb,nbaspd,nbt,nbtmax,
     $  nelect,nhydr,nhydx,noutpt,nstmax,uspec,zchar)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print the contributions of aqueous species to each mass balance.
c
      call prtpct(conc,csts,ctb,iopr,jcsort,jflag,narn1,narn2,
     $ nbasp,nbaspd,nbt,nbtmax,nelect,nhydr,noprmx,no2gaq,noutpt,
     $ nstmax,nsts,nstsmx,nstsr,uspec)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute and print the state of redox reactions that are not
c     constrained to be at equilibrium.
c
      call cdardx(actlg,actwlg,ah,ahrc,cdrsd,eh,ehfac,ehrc,
     $ farad,fo2lg,fo2lrc,jsflag,mosp,nbasp,nbaspd,nbt,nbtmax,ndrsd,
     $ ndrsmx,ndrsrd,no2gaq,nstmax,pe,perc,ph,xlke,xlksd)
c
      call prtrdx(ah,ahrc,cdrsd,eh,ehrc,fo2lg,fo2lrc,jflgi,
     $ jsflag,narn1,nbasp,nbaspd,nbt,nbtmax,ndrsd,ndrsmx,ndrsrd,
     $ nelect,nhydr,no2gaq,noutpt,nstmax,pe,perc,uspec)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the saturation index for the reaction associated to
c     each species. Do this for the reactions as they were written
c     on the data file.
c
      if (nxrn1 .gt. 0) then
        do ns = nxrn1,nxrn2
          if (actlg(ns) .le. -99999.) actlg(ns) = 0.0
        enddo
      endif
c
c     Calculate affinities and saturation indices using the 'd' set
c     of reactions.
c
      call gaffsd(actlg,afcnst,affpd,affsd,cdrsd,jflagd,jpflag,
     $ ncmpr,ndrsd,ndrsmx,ndrsrd,npt,nptmax,nst,nstmax,qxknph,sidrph,
     $ sidrsp,uphase,uspec,xbar,xlksd)
c
c     Compute and print saturation index and affinity tables for
c     reactions in the aqueous phase not constrained to be at
c     equilibrium.
c
      call prtsia(affsd,jflagd,jflgi,jsflag,narn1,narn2,nbasp,
     $ nbaspd,nbt,nbtmax,ndrsd,ndrsmx,ndrsrd,nhydr,noutpt,nrdxsp,
     $ nstmax,sidrsp,uspec)
c
      if (iopr(7) .ge. 0) then
c
c       Print saturation index and affinity tables for the
c       various non-aqueous phases.
c
        call prtsat(affpd,iern1,iern2,ilrn1,ilrn2,imrn1,imrn2,
     $  iopr,iopt,ixrn1,ixrn2,jpflag,noutpt,noprmx,noptmx,nptmax,
     $  sidrph,tolspf,uphase)
      endif
c
      if (iopt(4).ge.1 .and. nxti.gt.0) then
c
c       Print affinities and saturation indices for solid solutions
c       for which compositions were entered on the input file.
c
        write (noutpt,1600)
 1600   format(//16x,'--- Saturation States of Input Solid Solutions',
     $   ' ---',/)
c
        ncount = 0
        do nxi = 1,nxti
          np = npnxp(nxi)
          ncount = ncount + 1
          call prtsso(acflg,actlg,affpd,affsd,ixrn1,ixrn2,
     $    jsol,jsomax,ncmpr,noutpt,np,nptmax,nstmax,nxtmax,sidrph,
     $    sidrsp,tolspf,uspec,uphase,uxtype,xbar,xbarlg)
        enddo
        if (ncount .le. 0) write (noutpt,1710)
 1710   format(1x,'None')
      endif
c
      if (iopt(4).ge.1 .and. ixrn2.ge.ixrn1) then
c
c       Compute hypothetical saturation states for solid solutions.
c
        write (noutpt,1720)
 1720   format(/16x,'--- Saturation States of Hypothetical Solid',
     $  ' Solutions ---',/)
c
        ncount = 0
        do np = ixrn1,ixrn2
          nr1 = ncmpr(1,np)
          nr2 = ncmpr(2,np)
          nt = nr2 - nr1 + 1
          if (nt .gt. 1) then
c
            ier = 0
            ncount = ncount + 1
c
c           Calling sequence substitutions:
c             affpd for affp
c             affsd for affs
c             cdrsd for cdrs
c             ndrsd for ndrs
c             ndrsrd for ndrsr
c             xlksd for xlks
c
            call hpsat(acflg,act,actlg,afcnst,affpd,affsd,apx,bpx,
     $      cdrsd,eps100,iapxmx,ibpxmx,ier,iktmax,ixrn1,jflag,jpflag,
     $      jsflag,jsol,ncmpr,ndrsd,ndrsmx,ndrsrd,noutpt,np,nptmax,
     $      nstmax,nttyo,nxrn1,nxrn2,nxtmax,sidrsp,sidrph,uphase,
     $      uspec,wfac,xbar,xbarlg,xlksd)
c
c           Check to see if the hypothetical affinity calculation
c           converged.
c
            if (ier .le. 0) then
              qxknph(np) = .true.
            else
              qxknph(np) = .false.
              j2 = ilnobl(uphase(np))
              write (noutpt,1730) uphase(np)(1:j2)
              write (nttyo,1730) uphase(np)(1:j2)
 1730         format(/' * Warning - (EQ3NR/scripz) The hypothetical',
     $        /7x,'affinity calculation failed for ',a,'.')
              go to 310
            endif
c
            call prtsso(acflg,actlg,affpd,affsd,ixrn1,ixrn2,
     $      jsol,jsomax,ncmpr,noutpt,np,nptmax,nstmax,nxtmax,sidrph,
     $      sidrsp,tolspf,uspec,uphase,uxtype,xbar,xbarlg)
          endif
  310    continue
        enddo
        if (ncount .le. 0) write (noutpt,1740)
 1740   format(1x,'None')
      endif
c
cXX   Temporarily suppress these tables by setting neti = 0.
cXX   Maybe reactivate them in the future.
      neti = 0
cXX
      if (neti .gt. 0) then
c
c       Print affinities and saturation indices for generic ion
c       exchangers for which compositions were entered on the input
c       file.
c
        write (noutpt,1750)
 1750   format(//11x,'--- Saturation States of Input Generic Ion',
     $   ' Exchangers ---',/)
c
        ncount = 0
        do nei = 1,neti
          np = ngexpi(nei)
          ncount = ncount + 1
c
          call prtgex(acflg,actlg,affpd,affsd,cegexs,conc,
     $    egexjc,egexjf,egexpa,egexpc,egexs,egexw,iern1,iern2,ietmax,
     $    jern1,jern2,jetmax,jgext,kern1,kern2,ketmax,kgexsa,moph,
     $    mosp,netmax,ngexsa,ngext,noutpt,np,nptmax,nstmax,sidrph,
     $    sidrsp,tolspf,ugexj,ugexmo,uspec,uphase,xbar,xbarlg,xgexw,
     $    wkgwi)
c
        enddo
        if (ncount .le. 0) write (noutpt,1760)
 1760   format(1x,'None')
      endif
c
      if (net .gt. 0) then
c
c       Print compositions of generic ion exchangers.
c
        write (noutpt,1770)
 1770   format(/11x,'--- Generic Ion Exchangers ---',/)
c
        ncount = 0
        do np = iern1,iern2
          nr1 = ncmpr(1,np)
          nr2 = ncmpr(2,np)
          nt = nr2 - nr1 + 1
          if (nt .gt. 1) then
c
            ier = 0
            ncount = ncount + 1
c
            call prtgex(acflg,actlg,affpd,affsd,cegexs,conc,
     $      egexjc,egexjf,egexpa,egexpc,egexs,egexw,iern1,iern2,ietmax,
     $      jern1,jern2,jetmax,jgext,kern1,kern2,ketmax,kgexsa,moph,
     $      mosp,netmax,ngexsa,ngext,noutpt,np,nptmax,nstmax,sidrph,
     $      sidrsp,tolspf,ugexj,ugexmo,uspec,uphase,xbar,xbarlg,xgexw,
     $      wkgwi)
c
          endif
        enddo
        if (ncount .le. 0) write (noutpt,1790)
 1790   format(1x,'None')
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopr(8) .ge. 0) then
c
c       Print a table of equilibrium fugacities.
c
        call prtfug(jgsort,fugac,fugalg,jsflag,ngrn1,ngt,ngtmax,
     $  noutpt,nstmax,uspec)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,2000)
 2000 format(/' - - - - - - - - - - - - - - - - - - - - - - - - - - -',
     $ ' - - - - - - - - - - - -',/)
c
      end
