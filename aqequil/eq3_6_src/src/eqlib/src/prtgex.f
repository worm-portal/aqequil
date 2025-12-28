      subroutine prtgex(acflg,actlg,affpd,affsd,cegexs,conc,
     $ egexjc,egexjf,egexpa,egexpc,egexs,egexw,iern1,iern2,ietmax,
     $ jern1,jern2,jetmax,jgext,kern1,kern2,ketmax,kgexsa,moph,
     $ mosp,netmax,ngexsa,ngext,noutpt,np,nptmax,nstmax,sidrph,
     $ sidrsp,tolspf,ugexj,ugexmo,uspec,uphase,xbar,xbarlg,xgexw,
     $ wkgwi)
c
c     This subroutine prints tables describing the state and properties
c     of the np-th phase (a generic ion exchanger). It is analogous to
c     EQLIB/prtsso.f, which prints the same tables for solid solution
c     phases.
c
c     This subroutine is called by:
c
c       EQ3NR/scripx.f
c       EQ6/scripz.f
c
c----------------------------------------------------------------------
c
c     Principal input:
c
c       acflg  = array of log activity coefficients
c       actlg  = array of log activities of species
c       affpd  = array of phase affinities
c       affsd  = array of species affinities
c       sidrph = array of phase saturation indices
c       sidrsp = array of species saturation indices
c       tolspf = saturation print flag tolerance, used to flag those
c                  phases which are close to saturation
c       uphase = array of phase names
c       xbar   = array of mole fractions of species
c       xbarlg = array of log mole fractions of species
c
c     Principal output:
c
c       None
c
c----------------------------------------------------------------------
c
      implicit none
c
c----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ietmax,jetmax,ketmax,netmax,nptmax,nstmax
c
      integer noutpt
c
      integer jern1(jetmax,netmax),jern2(jetmax,netmax),jgext(netmax),
     $ kern1(netmax),kern2(netmax),kgexsa(ketmax,netmax),
     $ ngexsa(ietmax,jetmax,netmax),ngext(jetmax,netmax)
c
      integer iern1,iern2,np
c
      character*48 uspec(nstmax)
      character*24 ugexmo(netmax),uphase(nptmax)
      character*8 ugexj(jetmax,netmax)
c
      real*8 acflg(nstmax),actlg(nstmax),affpd(nptmax),affsd(nstmax),
     $ cegexs(ietmax,jetmax,netmax),conc(nstmax),egexjc(jetmax,netmax),
     $ egexjf(jetmax,netmax),egexpa(netmax),egexpc(netmax),
     $ egexs(ietmax,jetmax,netmax),egexw(ketmax,netmax),moph(nptmax),
     $ mosp(nstmax),sidrph(nptmax),sidrsp(nstmax),xbar(nstmax),
     $ xbarlg(nstmax),xgexw(ketmax,netmax)
c
      real*8 tolspf,wkgwi
c
c----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ie,je,j2,j3,j4,ke,nr1,nr2,ne,ns,nss
c
      integer ilnobl
c
      character*24 ux24
c
      real*8 aafx,afx,cxp,elx,ex,expkg,mxp,xlx,xx
c
      real*8 tlg
c
c----------------------------------------------------------------------
c
      mxp = moph(np)
      j2 = ilnobl(uphase(np))
c
      if (np.lt.iern1 .or. np.gt. iern2) then
        write (noutpt,1000) uphase(np)(1:j2)
        write (noutpt,1000) uphase(np)(1:j2)
 1000   format(/' * Error - (EQLIB/prtgex) Programming error trap:',
     $  /7x,"Can't write an output table for ",a," because",
     $  /7x,"it isn't a generic ion exchanger.")
        stop
      endif
c
      ne = np - iern1 + 1
c
      write (noutpt,1010) uphase(np)(1:j2)
 1010 format(/16x,'--- ',a,' ---',/)
c
      cxp = mxp*wkgwi
      j3 = ilnobl(ugexmo(ne))
      write (noutpt,1020) ugexmo(ne)(1:j3),mxp,cxp
 1020 format(/5x,'Exchanger model= ',a,
     $ /5x,'Amount of substrate (Z)= ',1pe11.4,
     $ ' mol',/5x,'Concentration of substrate= ',e11.4)
c
      if (egexpc(ne) .ne. 0.) then
        ex = egexpc(ne)*mxp
        expkg = ex*wkgwi
        write (noutpt,1024) egexpc(ne),ex,expkg
 1024   format(/5x,'Apparent cation exchange capacity= ',1pe11.4,
     $  ' eq/mol',/5x,'Total exchangeable cations= ',e11.4,' eq',
     $  /5x,'Total exchangeable cations= ',e11.4,' eq/kg.H2O')
      endif
c
      if (egexpa(ne) .ne. 0.) then
        ex = egexpa(ne)*mxp
        expkg = ex*wkgwi
        write (noutpt,1030) egexpa(ne),ex,expkg
 1030   format(/5x,'Apparent anion exchange capacity= ',1pe11.4,
     $  ' eq/mol',/5x,'Total exchangeable anions= ',e11.4,' eq',
     $  /5x,'Total exchangeable anions= ',e11.4,' eq/kg.H2O')
      endif
c
      do je = 1,jgext(ne)
        j4 = ilnobl(ugexj(je,ne))
        ex = egexjc(je,ne)*mxp
        expkg = ex*wkgwi
        write (noutpt,1040) ugexj(je,ne)(1:j4),egexjf(je,ne),
     $  egexjc(je,ne),ex,expkg
 1040   format(/7x,'Site: ',a,//9x,'Formal exchange capacity= ',
     $  1pe11.4,' eq/mol',/9x,'Apparent exchange capacity= ',
     $  e11.4,' eq/mol',/9x,'Total exchangeable ions= ',
     $  e11.4,' eq',/9x,'Total exchangeable ions= ',
     $  e11.4,' eq/kg.H2O',/)
c
c       Print equivalents and equivalents/kg.H2O.
c
        write (noutpt,1070)
 1070   format(/13x,'Component',19x,'eq',10x,'eq/kg.H2O',/)
c
        ns = jern1(je,ne) - 1
        do ie = 1,ngext(je,ne)
          ns = ns + 1
          ex = cegexs(ie,je,ne)*mosp(ns)
          expkg = ex*wkgwi
          nss = ngexsa(ie,je,ne)
          if (nss .le. 0) then
            ux24 = '__ '
          else
            ux24 = uspec(nss)
          endif
          if (nss .gt. 0) write (noutpt,1080) ux24,ex,expkg
 1080     format(10x,a24,3x,1pe11.4,3x,e11.4)
        enddo
c
c       Print equivalent fractions.
c
        write (noutpt,1110)
 1110   format(//13x,'Component',20x,'e',11x,'Log e',/)
c
        ns = jern1(je,ne) - 1
        do ie = 1,ngext(je,ne)
          ns = ns + 1
          ex = egexs(ie,je,ne)
          nss = ngexsa(ie,je,ne)
          if (nss .gt. 0) then
            ux24 = uspec(nss)
            elx = tlg(ex)
            write (noutpt,1130) ux24,ex,elx
 1130       format(10x,a24,3x,1pe11.4,3x,0pf9.4)
          endif
        enddo
        write (noutpt,1090)
 1090   format(1x)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1140)
 1140 format(/6x,'Whole-Phase Equivalent Fractions',
     $ //4x,'Component',20x,'e',11x,'Log e',/)
c
      do ke = kern1(ne),kern2(ne)
        nss = kgexsa(ke,ne)
        if (nss .gt. 0) then
          ex = egexw(ke,ne)
          elx = tlg(ex)
          write (noutpt,1160) uspec(nss),ex,elx
 1160     format(1x,a24,3x,1pe11.4,3x,0pf9.4)
        endif
      enddo
c
      write (noutpt,1170)
 1170 format(//6x,'Whole-Phase Mole Fractions',
     $ //4x,'Component',20x,'x',11x,'Log x',/)
c
      do ke = kern1(ne),kern2(ne)
        nss = kgexsa(ke,ne)
        if (nss .gt. 0) then
          xx = xgexw(ke,ne)
          xlx = tlg(xx)
          write (noutpt,1160) uspec(nss),xx,xlx
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print mole fractions, activity coefficients, and activities.
c
      write (noutpt,1180)
 1180 format(//6x,'Thermodynamic Summary',
     $ //4x,'Component',20x,'x',11x,'Log x',3x,'Log lambda',2x,
     $ 'Log activity',/)
c
      do je = 1,jgext(ne)
        nr1 = jern1(je,ne)
        nr2 = jern2(je,ne)
        do ns = nr1,nr2
          if (xbar(ns) .gt. 0.) then
            write (noutpt,1190) uspec(ns),xbar(ns),xbarlg(ns),
     $      acflg(ns),actlg(ns)
 1190       format(1x,a24,3x,1pe11.4,3(3x,0pf9.4))
          endif
        enddo
        write (noutpt,1090)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print saturation states and affinities.
c
      write (noutpt,1200)
 1200 format(/4x,'Exchanger',21x,'Log Q/K',9x,'Aff, kcal',4x,'State',/)
      afx = affpd(np)
      aafx = abs(afx)
      if (aafx .le. tolspf) then
        write (noutpt,1210) uphase(np),sidrph(np),affpd(np)
 1210   format(1x,a24,2x,2(3x,f13.4),3x,'SATD')
      elseif (afx .gt. tolspf) then
        write (noutpt,1220) uphase(np),sidrph(np),affpd(np)
 1220   format(1x,a24,2x,2(3x,f13.4),3x,'SSATD')
      else
        write (noutpt,1230) uphase(np),sidrph(np),affpd(np)
 1230   format(1x,a24,2x,2(3x,f13.4))
      endif
      write (noutpt,1090)
c
      do je = 1,jgext(ne)
        nr1 = jern1(je,ne)
        nr2 = jern2(je,ne)
        do ns = nr1,nr2
          if (xbar(ns) .gt. 0.) then
            afx = affsd(ns)
            aafx = abs(afx)
            if (aafx .le. tolspf) then
              write (noutpt,1240) uspec(ns),sidrsp(ns),affsd(ns)
 1240         format(3x,a24,2(3x,f13.4),3x,'SATD')
            elseif (afx .gt. tolspf) then
              write (noutpt,1250) uspec(ns),sidrsp(ns),affsd(ns)
 1250         format(3x,a24,2(3x,f13.4),3x,'SSATD')
            else
              write (noutpt,1260) uspec(ns),sidrsp(ns),affsd(ns)
 1260         format(3x,a24,2(3x,f13.4))
            endif
          endif
        enddo
        write (noutpt,1090)
      enddo
c
      end
