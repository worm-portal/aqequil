      subroutine intrct(cbsr,cbsri,cesr,cesri,egers,egersi,ibsrti,
     $ iern1,iern2,iesrti,ietmax,igerti,iktmax,imrn1,imrn2,ixrn1,
     $ ixrn2,ixrti,jcode,jetmax,jgerti,jgext,narn1,narn2,nbaspd,
     $ nbt,nbtmax,nbt1mx,ncmpr,nct,nctmax,nertmx,netmax,ngexsa,
     $ ngext,ngrn1,ngrn2,noutpt,nptmax,nrct,nrctmx,nrndex,nsrtmx,
     $ nstmax,nttyo,nxridx,nxrtmx,rxbar,rxbari,ubsri,ucxri,uelem,
     $ uesri,ugerji,ugermo,ugersi,ugexj,ugexmo,uphase,ureac,uspec,
     $ xgers,xgersi)
c
c     This subroutine assigns phase or species indices corresponding to
c     reactants.
c
c     This subroutine is called by:
c
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
c       jcode  = flag which determines the type of reactant:
c
c                  0 =  Mineral
c                  1 =  Solid solution
c                  2 =  Special reactant
c                  3 =  Aqueous species
c                  4 =  Gas
c                  5 =  Generic ion exchanger
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ietmax,iktmax,jetmax,nbtmax,nbt1mx,nctmax,nertmx,netmax,
     $ nptmax,nrctmx,nsrtmx,nstmax,nxrtmx
c
      integer ibsrti(nsrtmx),iesrti(nsrtmx),igerti(jetmax,nertmx),
     $ ixrti(nxrtmx),jcode(nrctmx),jgerti(nertmx),jgext(netmax),
     $ ncmpr(2,nptmax),nbaspd(nbtmax),ngexsa(ietmax,jetmax,netmax),
     $ ngext(jetmax,netmax),nrndex(nrctmx),nxridx(nrctmx)
c
      integer iern1,iern2,imrn1,imrn2,ixrn1,ixrn2,narn1,narn2,nbt,nct,
     $ ngrn1,ngrn2,noutpt,nrct,nttyo
c
      character*48 uspec(nstmax)
      character*24 ubsri(nbt1mx,nsrtmx),ucxri(iktmax,nxrtmx),
     $ ugermo(nertmx),ugersi(ietmax,jetmax,nertmx),ugexmo(netmax),
     $ uphase(nptmax),ureac(nrctmx)
      character*8 uelem(nctmax),uesri(nctmax,nsrtmx),
     $ ugerji(jetmax,nertmx),ugexj(jetmax,netmax)
c
      real*8 cesr(nctmax,nsrtmx),cesri(nctmax,nsrtmx),
     $ cbsr(nbt1mx,nsrtmx),cbsri(nbt1mx,nsrtmx),
     $ egers(ietmax,jetmax,netmax),egersi(ietmax,jetmax,nertmx),
     $ rxbar(iktmax,nxrtmx),rxbari(iktmax,nxrtmx),
     $ xgers(ietmax,jetmax,netmax),xgersi(ietmax,jetmax,nertmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ie,iei,ik,iki,je,jei,j2,j3,j4,j5,j6,nb,nbi,nc,nci,ne,ner,
     $ nerr,np,nrc,nr1,nr2,ns,nsr,nss,nt,nxr
c
      integer ilnobl
c
      character*24 urn,ux
      character*24 uaq,umn,ugs,uss,uge
      character*16 ux16
      character*8 uxd,ux8a,ux8b
c
      real*8 ex,exc,exo,rx,rxo,rxsum
c
c-----------------------------------------------------------------------
c
      data uaq    /'aqueous species         '/,
     $     umn    /'minerals                '/,
     $     ugs    /'gases                   '/,
     $     uss    /'solid solutions         '/,
     $     uge    /'generic ion exchangers  '/
c
c-----------------------------------------------------------------------
c
      if (nrct .eq. 0) go to 999
c
      nerr = 0
      nsr = 0
      nxr = 0
      ner = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do nrc = 1,nrct
        urn = ureac(nrc)
c
        if (jcode(nrc) .eq. 0) then
c
c         Pure minerals.
c
          do np = imrn1,imrn2
           if (urn(1:24) .eq. uphase(np)(1:24)) then
             nrndex(nrc) = np
             go to 200
           endif
          enddo
c
          j2 = ilnobl(urn)
          j3 = ilnobl(umn)
          write (noutpt,1000) urn(1:j2),umn(1:j3)
          write (nttyo,1000) urn(1:j2),umn(1:j3)
 1000     format(/' * Error - (EQ6/intrct) The reactant ',a,
     $    " isn't among",/7x,'the loaded ',a,'.')
          nerr = nerr + 1
c
        elseif (jcode(nrc) .eq. 1) then
c
c         Solid solutions.
c
          nxr = nxr + 1
          do np = ixrn1,ixrn2
            if (urn(1:24) .eq. uphase(np)(1:24)) then
              nrndex(nrc) = np
              nxridx(nrc) = nxr
              nr1 = ncmpr(1,np)
              nr2 = ncmpr(2,np)
              nt = nr2 - nr1 + 1
c
c             Make sure that the mole fractions sum to unity.
c
              rxsum = 0.
              do iki = 1,ixrti(nxr)
                rxsum = rxsum + rxbari(iki,nxr)
              enddo
c
              if (abs(rxsum - 1.0) .gt. 1.e-3) then
                j3 = ilnobl(urn)
                write (ux16,'(g11.4)') rxsum
                call lejust(ux16)
                j4 = ilnobl(ux16)
                write (noutpt,1002) urn(1:j3),ux16(1:j4)
                write (nttyo,1002) urn(1:j3),ux16(1:j4)
 1002           format(/' * Warning - (EQ6/intrct) The mole fractions',
     $          ' given for the composition',/7x,'of the solid',
     $          ' solution reactant ',a,' sum to',/7x,a,', not unity.',
     $          ' They will be renormalized.')
c
c               Renormalize the mole fractions.
c
                write (noutpt,1004)
                write (nttyo,1004)
 1004           format(/5x,'Renormalizing the Input Mole Fractions',
     $          //9x,'Component',19x,'Old x',9x,'New x',/)
                do iki = 1,ixrti(nxr)
                  rxo = rxbari(iki,nxr)
                  rx = rxo/rxsum
                  write (noutpt,1006) ucxri(iki,nxr),rxo,rx
                  write (nttyo,1006) ucxri(iki,nxr),rxo,rx
 1006             format(7x,a24,3x,1pe11.4,3x,e11.4)
                  rxbari(iki,nxr) = rx
                enddo
              endif
c
c             Decode the solid solution composition.
c
              do iki = 1,ixrti(nxr)
                ux = ucxri(iki,nxr)
                ns = nr1 - 1
                do ik = 1,nt
                  ns = ns + 1
                  if (uspec(ns)(1:24) .eq. ux(1:24)) then
                    rxbar(ik,nxr) = rxbari(iki,nxr)
                    go to 100
                  endif
                enddo
c
                j2 = ilnobl(ux)
                j3 = ilnobl(urn)
                write (noutpt,1010) ux(1:j2),urn(1:j3)
                write (nttyo,1010) ux(1:j2),urn(1:j3)
 1010           format(/' * Error - (EQ6/intrct) The species ',a,
     $          " doesn't appear",/7x,'on the data file as a',
     $          ' component of ',a,'.')
                nerr = nerr + 1
c
  100           continue
              enddo
              go to 200
            endif
          enddo
c
          j2 = ilnobl(urn)
          j3 = ilnobl(uss)
          write (noutpt,1000) urn(1:j2),uss(1:j3)
          write (nttyo,1000) urn(1:j2),uss(1:j3)
          nerr = nerr + 1
c
        elseif (jcode(nrc) .eq. 2) then
c
c         Special reactants.
c
          nsr = nsr + 1
          nrndex(nrc) = nsr
          nxridx(nrc) = nsr
c
c         Decode the chemical composition.
c
          do nci = 1,iesrti(nsr)
            uxd = uesri(nci,nsr)
c
            do nc = 1,nct
              if (uelem(nc)(1:8) .eq. uxd(1:8)) then
                cesr(nc,nsr) = cesri(nci,nsr)
                go to 110
              endif
            enddo
c
            j2 = ilnobl(uxd)
            j3 = ilnobl(urn)
            write (noutpt,1020) uxd(1:j2),urn(1:j3)
            write (nttyo,1020) uxd(1:j2),urn(1:j3)
 1020       format(/' * Error - (EQ6/intrct) The element ',a,
     $      ' in special reactant',/7x,a," isn't on the data file.")
            nerr = nerr + 1
c
  110       continue
          enddo
c
c         Decode the chemical reaction.
c
          ux = ubsri(1,nsr)
          if (ux(1:24) .eq. ureac(nrc)(1:24)) then
            cbsr(nbt1mx,nsr) = cbsri(1,nsr)
          else
            j2 = ilnobl(ux)
            j3 = ilnobl(urn)
            write (noutpt,1030) ux(1:j2),urn(1:j3)
            write (nttyo,1030) ux(1:j2),urn(1:j3)
 1030       format(' * Error - (EQ6/intrct) The species ',a,' appears',
     $      ' first',/7x,'in the reaction associated with the special',
     $      ' reactant',/7x,a,'. The special reactant itself must',
     $      ' appear first.')
            nerr = nerr + 1
          endif
c
          do nbi = 2,ibsrti(nsr)
            ux = ubsri(nbi,nsr)
c
            do nb = 1,nbt
              ns = nbaspd(nb)
              if (uspec(ns)(1:24) .eq. ux(1:24)) then
                cbsr(nb,nsr) = cbsri(nbi,nsr)
                go to 120
              endif
            enddo
c
            j2 = ilnobl(ux)
            j3 = ilnobl(urn)
            write (noutpt,1040) ux(1:j2),urn(1:j3)
            write (nttyo,1040) ux(1:j2),urn(1:j3)
 1040       format(' * Error - (EQ6/intrct) The species ',a,' which',
     $      ' appears',/7x,'in the reaction associated with the',
     $      ' special reactant',/7x,a," isn't in the data file',
     $      ' basis set.")
            nerr = nerr + 1
  120       continue
          enddo
c
        elseif (jcode(nrc) .eq. 3) then
c
c         Aqueous species.
c
          do ns = narn1,narn2
            if (urn(1:24) .eq. uspec(ns)(1:24)) then
              nrndex(nrc) = ns
              go to 200
            endif
          enddo
c
          j2 = ilnobl(urn)
          j3 = ilnobl(uaq)
          write (noutpt,1000) urn(1:j2),uaq(1:j3)
          write (nttyo,1000) urn(1:j2),uaq(1:j3)
          nerr = nerr + 1
c
        elseif (jcode(nrc) .eq. 4) then
c
c         Gases.
c
          do ns = ngrn1,ngrn2
            if (urn(1:24) .eq. uspec(ns)(1:24)) then
              nrndex(nrc) = ns
              go to 200
            endif
          enddo
c
          j2 = ilnobl(urn)
          j3 = ilnobl(ugs)
          write (noutpt,1000) urn(1:j2),ugs(1:j3)
          write (nttyo,1000) urn(1:j2),ugs(1:j3)
          nerr = nerr + 1
c
        elseif (jcode(nrc) .eq. 5) then
c
c         Generic ion exchangers.
c
          ner = ner + 1
          do np = iern1,iern2
            if (urn(1:24) .eq. uphase(np)(1:24)) then
              nxridx(nrc) = ner
              nrndex(nrc) = np
              ne = np - iern1 + 1
c
              if (ugermo(ner)(1:24) .ne. ugexmo(ne)(1:24)) then
                j2 = ilnobl(urn)
                j5 = ilnobl(ugermo(ner))
                j6 = ilnobl(ugexmo(ne))
                write (noutpt,1042) urn(1:j2),ugermo(ner)(1:j5),
     $          ugexmo(ne)(1:j6)
                write (nttyo,1042) urn(1:j2),ugermo(ner)(1:j5),
     $          ugexmo(ne)(1:j6)
 1042           format(/' * Error - (EQ6/intrct) The generic exchanger',
     $          ' phase reactant',/7x,a,' has the specified exchange',
     $          ' model',/7x,'"',a,'", but the defined exchange model',
     $          ' for this',/7x,'exchanger phase is "',a,'".')
                nerr = nerr + 1
              endif
c
              if (jgerti(ner) .ne. jgext(ne)) then
                j2 = ilnobl(urn)
                write (ux8a,'(i5)') jgerti(ner)
                call lejust(ux8a)
                j5 = ilnobl(ux8a)
                write (ux8b,'(i5)') jgext(ne)
                call lejust(ux8b)
                j6 = ilnobl(ux8b)
                write (noutpt,1044) urn(1:j2),ux8a(1:j5),ux8b(1:j6)
                write (nttyo,1044) urn(1:j2),ux8a(1:j5),ux8b(1:j6)
 1044           format(/' * Error - (EQ6/intrct) The generic exchanger',
     $          ' phase reactant',/7x,a,' has a specified composition',
     $          ' including',/7x,a,' exchange site(s), but this phase',
     $          ' actually has ',a,' exchange',/7x,'site(s). The',
     $          ' number of such sites must match.')
                nerr = nerr + 1
              endif
c
              do jei = 1,jgerti(ner)
c
c               Make sure that the exchange fractions on the current
c               site sum to unity.
c
                exc = 0.
                do iei = 1,igerti(jei,ner)
                  exc = exc + egersi(iei,jei,ner)
                enddo
c
                if (abs(exc - 1.0) .gt. 1.e-3) then
                  j3 = ilnobl(urn)
                  j4 = ilnobl(ugerji(jei,ner))
                  write (ux16,'(g11.4)') exc
                  call lejust(ux16)
                  j5 = ilnobl(ux16)
                  write (noutpt,1050) ugerji(1,j5),urn(1:j3),ux16(1:j5)
                  write (nttyo,1050) ugerji(1,j5),urn(1:j3),ux16(1:j5)
 1050             format(/' * Warning - (EQ6/intrct) The exchange',
     $            ' fractions given for the composition',/7x,'of',
     $            ' exchange site ',a,' of the generic ion exchanger',
     $            /7x,'reactant ',a,' sum to ',a,', not unity.',
     $            /7x,'They will be renormalized.')
c
c                 Renormalize the exchange fractions.
c
                  write (noutpt,1060)
                  write (nttyo,1060)
 1060             format(/5x,'Renormalizing the Input Exchange',
     $            ' Fractions',//9x,'Component',19x,'Old e',9x,
     $            'New e',/)
                  do iei = 1,igerti(jei,ner)
                    exo = egersi(iei,jei,ner)
                    ex = exo/exc
                    write (noutpt,1006) ugersi(iei,jei,ner),exo,ex
                    write (nttyo,1006) ugersi(iei,jei,ner),exo,ex
                    egersi(iei,jei,ner) = ex
                  enddo
                endif
c
c               Decode the composition of the current exchange site.
c
                do je = 1,jgext(ne)
                  if (ugexj(je,ne)(1:8) .eq. ugerji(jei,ner)(1:8))
     $            go to 150
                enddo
c
                j4 = ilnobl(ugerji(jei,ner))
                write (noutpt,1070) ugerji(jei,ner)(1:j4),urn(1:j3)
                write (nttyo,1070) ugerji(jei,ner)(1:j4),urn(1:j3)
 1070           format(/' * Error - (EQ6/intrct) The exchange site ',
     $          a,' specified for',/7x,'the generic ion exchanger',
     $          ' reactant ',a,'does not',/7x,'match any site',
     $          ' defined for that exchanger.')
                nerr = nerr + 1
                go to 170
c
  150           do iei = 1,igerti(jei,ner)
                  ux = ugersi(iei,jei,ner)
c
c                 Match with the exchange species name. This contains
c                 the site name.
c
                  ns = ncmpr(1,np)
                  do je = 1,jgext(ne)
                    do ie = 1,ngext(je,ne)
                      ns = ns + 1
                      if (uspec(ns)(1:24) .eq. ux(1:24)) then
                        egers(ie,je,ne) = egersi(iei,jei,ner)
                        go to 160
                      endif
                    enddo
                  enddo
c
c                 No match was found above. Match with the name of
c                 the aqueous species which exchanges. This does
c                 not contain the site name.
c
                  do je = 1,jgext(ne)
                    do ie = 1,ngext(je,ne)
                      nss = ngexsa(ie,je,ne)
                      if (nss .gt. 0) then
                        if (uspec(nss)(1:24) .eq. ux(1:24)) then
                          egers(ie,je,ne) = egersi(iei,jei,ner)
                          go to 160
                        endif
                      endif
                    enddo
                  enddo
c
                  j2 = ilnobl(ux)
                  j3 = ilnobl(urn)
                  write (noutpt,1080) ux(1:j2),urn(1:j3)
                  write (nttyo,1080) ux(1:j2),urn(1:j3)
 1080             format(/' * Error - (EQ6/intrct) The species ',
     $            a,' specified for',/7x,'the generic ion exchanger',
     $            ' reactant ',a,' does not',/7x,'match any species',
     $            ' defined for that exchanger.')
                  nerr = nerr + 1
c
  160             continue
                enddo
  170           continue
              enddo
              go to 200
            endif
          enddo
c
          j2 = ilnobl(urn)
          j3 = ilnobl(uge)
          write (noutpt,1000) urn(1:j2),uge(1:j3)
          write (nttyo,1000) urn(1:j2),uge(1:j3)
          nerr = nerr + 1
c
        endif
c
  200   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nerr .gt. 0) stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
