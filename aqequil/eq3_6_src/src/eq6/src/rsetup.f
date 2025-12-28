      subroutine rsetup(atwt,cbsr,cesr,iern1,ietmax,iindx1,iktmax,
     $ jcode,jern1,jern2,jetmax,jgext,kbt,kmax,mwtges,mwtrc,mwtsp,
     $ nbaspd,nbt,nbtmax,nbt1mx,ncmpr,nct,nctmax,nertmx,netmax,ngext,
     $ noutpt,nptmax,nrct,nrctmx,nrndex,nsrtmx,nstmax,nsts,nstsmx,
     $ nstsr,nttyo,nxridx,nxrtmx,rxbar,ureac,uspec,vosp0,vreac,xgers)
c
c     This subroutine assigns the molecular weights and molar volumes
c     of the reactants. Volumes are treated only for solid reactants
c     for later comparison with the volume of solid products.
c
c        nrct = total number of reactants
c        ureac = name of reactant
c        mwtrc = molecular weight of reactant
c        vreac = molar volume of (solid) reactant
c        jcode = reactant type code
c                0   Mineral
c                1   Solid solution
c                2   Special reactant
c                3   Aqueous species
c                4   Gas
c                5   Generic ion exchanger
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
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ietmax,iktmax,jetmax,kmax,nbtmax,nbt1mx,nctmax,nertmx,
     $ netmax,nptmax,nrctmx,nsrtmx,nstmax,nstsmx,nxrtmx
c
      integer noutpt,nttyo
c
      integer iindx1(kmax),jcode(nrctmx),jern1(jetmax,netmax),
     $ jern2(jetmax,netmax),jgext(netmax),nbaspd(nbtmax),
     $ ncmpr(2,nptmax),ngext(jetmax,netmax),nrndex(nrctmx),
     $ nsts(nstsmx),nstsr(2,nstmax),nxridx(nrctmx)
c
      integer iern1,kbt,nbt,nct,nrct
c
      character*48 uspec(nstmax)
      character*24 ureac(nrctmx)
c
      real*8 atwt(nctmax),cbsr(nbt1mx,nsrtmx),cesr(nctmax,nsrtmx),
     $ mwtges(netmax),mwtrc(nrctmx),mwtsp(nstmax),
     $ rxbar(iktmax,nxrtmx),vosp0(nstmax),vreac(nrctmx),
     $ xgers(ietmax,jetmax,nertmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ie,ik,je,j2,jlen,kcol,n,nb,nbb,nc,ne,ner,np,nrc,nrn1,nrn2,
     $ nr1,nr2,ns,nsr,nss,nxr
c
      integer ilnobl
c
      logical qstop
c
      character*56 uspn56
c
c-----------------------------------------------------------------------
c
      qstop = .false.
      do nrc = 1,nrct
        if (jcode(nrc) .eq. 0) then
c
c         Pure mineral reactants.
c
          np = nrndex(nrc)
          ns = ncmpr(1,np)
          mwtrc(nrc) = mwtsp(ns)
          vreac(nrc) = vosp0(ns)
c
          nr1 = nstsr(1,ns)
          nr2 = nstsr(1,ns)
          do n = nr1,nr2
            nb = nsts(n)
            do kcol = 1,kbt
              nbb = iindx1(kcol)
              if (nbb .eq. nb) go to 100
            enddo
c
            nss = nbaspd(nb)
            j2 = ilnobl(ureac(nrc))
c
c           Calling sequence substitutions:
c             uspec(nss) for unam48
c
            call fmspnx(jlen,uspec(nss),uspn56)
c
            write (noutpt,2100) ureac(nrc)(1:j2),uspn56(1:jlen)
            write (nttyo,2100) ureac(nrc)(1:j2),uspn56(1:jlen)
 2100       format(/' * Error - (EQ6/rsetup) The reactant ',a,' is',
     $      ' composed of',/7x,'the basis species ',a,', which',
     $      /7x,'is not in the active basis set.')
            qstop = .true.
  100       continue
          enddo
c
        elseif (jcode(nrc) .eq. 1) then
c
c         Solid solution reactants.
c
          np = nrndex(nrc)
          nrn1 = ncmpr(1,np)
          nrn2 = ncmpr(2,np)
          nxr = nxridx(nrc)
          mwtrc(nrc) = 0.
          vreac(nrc) = 0.
c
          ik = 0
          do ns = nrn1,nrn2
            ik = ik + 1
            mwtrc(nrc) = mwtrc(nrc) + rxbar(ik,nxr)*mwtsp(ns)
            vreac(nrc) = vreac(nrc) + rxbar(ik,nxr)*vosp0(ns)
c
            nr1 = nstsr(1,ns)
            nr2 = nstsr(1,ns)
            do n = nr1,nr2
              nb = nsts(n)
              do kcol = 1,kbt
                nbb = iindx1(kcol)
                if (nbb .eq. nb) go to 110
              enddo
c
              nss = nbaspd(nb)
              j2 = ilnobl(ureac(nrc))
c
c             Calling sequence substitutions:
c               uspec(nss) for unam48
c
              call fmspnx(jlen,uspec(nss),uspn56)
c
              write (noutpt,2100) ureac(nrc)(1:j2),uspn56(1:jlen)
              write (nttyo,2100) ureac(nrc)(1:j2),uspn56(1:jlen)
              qstop = .true.
  110         continue
            enddo
          enddo
c
        elseif (jcode(nrc) .eq. 2) then
c
c         Special reactants.
c
          mwtrc(nrc) = 0.
          nsr = nrndex(nrc)
          do nc = 1,nct
            if (cesr(nc,nsr) .ne. 0.)
     $      mwtrc(nrc) = mwtrc(nrc) + atwt(nc)*cesr(nc,nsr)
          enddo
          vreac(nrc) = 0.
c
          do nb = 1,nbt
            if (cbsr(nb,nsr) .ne. 0.) then
              do kcol = 1,kbt
                nbb = iindx1(kcol)
                if (nbb .eq. nb) go to 120
              enddo
c
              nss = nbaspd(nb)
              j2 = ilnobl(ureac(nrc))
c
c             Calling sequence substitutions:
c               uspec(nss) for unam48
c
              call fmspnx(jlen,uspec(nss),uspn56)
c
              write (noutpt,2100) ureac(nrc)(1:j2),uspn56(1:jlen)
              write (nttyo,2100) ureac(nrc)(1:j2),uspn56(1:jlen)
              qstop = .true.
            endif
  120       continue
          enddo
c
        elseif (jcode(nrc).eq.3 .or. jcode(nrc).eq.4) then
c
c         Aqueous species reactants and gas reactants.
c
          ns = nrndex(nrc)
          mwtrc(nrc) = mwtsp(ns)
          vreac(nrc) = 0.
c
          nr1 = nstsr(1,ns)
          nr2 = nstsr(1,ns)
          do n = nr1,nr2
            nb = nsts(n)
            do kcol = 1,kbt
              nbb = iindx1(kcol)
              if (nbb .eq. nb) go to 130
            enddo
c
            nss = nbaspd(nb)
            j2 = ilnobl(ureac(nrc))
c	
c           Calling sequence substitutions:
c             uspec(nss) for unam48
c
            call fmspnx(jlen,uspec(nss),uspn56)
c
            write (noutpt,2100) ureac(nrc)(1:j2),uspn56(1:jlen)
            write (nttyo,2100) ureac(nrc)(1:j2),uspn56(1:jlen)
            qstop = .true.
  130       continue
          enddo
c
        elseif (jcode(nrc) .eq. 5) then
c
c         Generic ion exchangers.
c
          np = nrndex(nrc)
          ner = nxridx(nrc)
          ne = np - iern1 + 1
c
          mwtrc(nrc) = mwtges(ne)
          vreac(nrc) = 0.
c
          do je = 1,jgext(ne)
            ns = jern1(je,ne) - 1
            do ie = 1,ngext(je,ne)
              ns = ns + 1
              mwtrc(nrc) = mwtrc(nrc) + xgers(ie,je,ner)*mwtsp(ns)
            enddo
          enddo
c
          do je = 1,jgext(ne)
            do ns = jern1(je,ne),jern2(je,ne)
              nr1 = nstsr(1,ns)
              nr2 = nstsr(1,ns)
              do n = nr1,nr2
                nb = nsts(n)
                do kcol = 1,kbt
                  nbb = iindx1(kcol)
                  if (nbb .eq. nb) go to 140
                enddo
c
                nss = nbaspd(nb)
                j2 = ilnobl(ureac(nrc))
c
c               Calling sequence substitutions:
c                 uspec(nss) for unam48
c
                call fmspnx(jlen,uspec(nss),uspn56)
c
                write (noutpt,2100) ureac(nrc)(1:j2),uspn56(1:jlen)
                write (nttyo,2100) ureac(nrc)(1:j2),uspn56(1:jlen)
                qstop = .true.
c
  140           continue
              enddo
            enddo
          enddo
        endif
      enddo
c
      if (qstop) stop
c
      end
