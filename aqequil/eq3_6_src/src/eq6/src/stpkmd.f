      subroutine stpkmd(cbsri,cbsr1,cdac,cesri,cesr1,csigma,eact,
     $ fkrc,hact,iact,ibsrti,ibsrt1,iesrti,iesrt1,iktmax,imchmx,
     $ imech,iopt,ixrti,jcode,jreac,modr,morr,morrw1,nbt1mx,nctmax,
     $ ndact,ndctmx,noptmx,noutpt,nprob,nrct,nrctmx,nrk,nsk,nsrt,
     $ nsrtmx,ntitl1,ntitmx,nttyo,nxrt,nxrtmx,rkb,rxbari,sfcar,
     $ ssfcar,trkb,ubsri,ubsr1,ucxri,udac,uesri,uesr1,ureac,ureac1,
     $ utitl1,vreac)
c
c     This subroutine sets up certain variables and arrays for writing
c     the top half of an EQ6 input file when one of the advanced pickup
c     file options is selected. Presently there is only one such option,
c     which is to replace the set of reactants by a special reactant
c     defined as the aqueous solution at the last point of reaction
c     progress. This subroutine is somewhat analogous to EQ3NR/stpk36.f,
c     which sets up certain variables and arrays for the advanced
c     pickup file options for EQ3NR.
c
c     This subroutine should only be called if iopt(1) is greater than
c     zero.
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
c     Calling sequence variable declarations.
c
      integer iktmax,imchmx,nbt1mx,nctmax,ndctmx,noptmx,nsrtmx,ntitmx,
     $ nrctmx,nxrtmx
c
      integer noutpt,nttyo
c
      integer iact(imchmx,2,nrctmx),ibsrti(nsrtmx),iesrti(nsrtmx),
     $ imech(2,nrctmx),iopt(noptmx),ixrti(nxrtmx),jcode(nrctmx),
     $ jreac(nrctmx),ndact(imchmx,2,nrctmx),nrk(2,nrctmx),nsk(nrctmx)
c
      integer ibsrt1,iesrt1,nprob,nrct,nsrt,ntitl1,nxrt
c
      character*80 utitl1(ntitmx)
      character*24 ubsri(nbt1mx,nsrtmx),ubsr1(nbt1mx),
     $ ucxri(iktmax,nxrtmx),udac(ndctmx,imchmx,2,nrctmx),ureac(nrctmx)
      character*8 uesri(nctmax,nsrtmx),uesr1(nctmax)
c
      character*24 ureac1
c
      real*8 cbsri(nbt1mx,nsrtmx),cbsr1(nbt1mx),
     $ cdac(ndctmx,imchmx,2,nrctmx),cesri(nctmax,nsrtmx),cesr1(nctmax),
     $ csigma(imchmx,2,nrctmx),eact(imchmx,2,nrctmx),fkrc(nrctmx),
     $ hact(imchmx,2,nrctmx),modr(nrctmx),morr(nrctmx),
     $ rkb(imchmx,2,nrctmx),rxbari(iktmax,nxrtmx),sfcar(nrctmx),
     $ ssfcar(nrctmx),trkb(imchmx,2,nrctmx),vreac(nrctmx)
c
      real*8 morrw1
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n,nmax
c
c-----------------------------------------------------------------------
c
      if (iopt(20) .le. 0) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set up the initial part of the new main title.
c
      call initcb(utitl1,ntitmx)
      ntitl1 = 25
      utitl1(1) = 'EQ6 input file name= sample.6i'
      utitl1(2) = 'Description= "Sample"'
      utitl1(3) = 'Version level= 8.0'
      utitl1(4) = 'Revised mm/dd/yy    Revisor= Username'
      utitl1(6) = '  This is a sample EQ6 input file written as an EQ6 P
     $ICKUP file.'
      utitl1(7) = 'The EQ6 input file used to generate this output is id
     $entified in'
      utitl1(8) = 'the second title given below.'
      utitl1(10) = '  You are expected to modify this EQ6 input file to
     $meet your own needs.'
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Delete all current reactants.
c
      nrct = 0
      call initiz(jcode,nrctmx)
      call initiz(jreac,nrctmx)
      call initiz(nsk,nrctmx)
c
      nsrt = 0
      call initiz(ibsrti,nsrtmx)
      call initiz(iesrti,nsrtmx)
c
      nmax = nbt1mx*nsrtmx
      call initaz(cbsri,nmax)
      call initcb(ubsri,nmax)
c
      nmax = nctmax*nsrtmx
      call initaz(cesri,nmax)
      call initcb(uesri,nmax)
c
      call initiz(ixrti,nxrtmx)
c
      nxrt = 0
      nmax = iktmax*nxrtmx
      call initaz(rxbari,nmax)
      call initcb(ucxri,nmax)
c
      call initcb(ureac,nrctmx)
      call initaz(fkrc,nrctmx)
      call initaz(sfcar,nrctmx)
      call initaz(ssfcar,nrctmx)
c
      nmax = 2*nrctmx
      call initiz(imech,nmax)
      call initiz(nrk,nmax)
c
      nmax = imchmx*2*nrctmx
      call initaz(csigma,nmax)
      call initaz(rkb,nmax)
      call initaz(trkb,nmax)
      call initaz(eact,nmax)
      call initaz(hact,nmax)
      call initiz(iact,nmax)
      call initiz(ndact,nmax)
c
      nmax = ndctmx*imchmx*2*nrctmx
      call initaz(cdac,nmax)
      call initcb(udac,nmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write an EQ6 input file in which "Fluid 2" is the sole reactant
c     added to the equilibrium system (ES) containing "Fluid 1" and
c     possible other phases. "Fluid 2" is defined as a special
c     reactant and appears in the top half of the file. It is the
c     aqueous solution from the end of the first reaction-path problem
c     of one or more such problems included on the current EQ6 input
c     file. "Fluid 1" is the aqueous solution from the end of the
c     current such problem, and the whole ES described is that at
c     the end of this problem.
c
c     If this option is chosen for the first (perhaps the only)
c     problem on the EQ6 input file, the pickup file produced
c     is set up for mixing of "Fluid 2" with an ES containing an
c     identical aqueous solution. This would normally be not very
c     useful as a whole, thought the top and bottom halves might be
c     separately useful.
c
      utitl1(12) = '  This sample file has "Fluid 2" mixing into an equi
     $librium system (ES)'
      utitl1(13) = 'which contains "Fluid 1". "Fluid 2" is the aqueous s
     $olution at the end'
      utitl1(14) = 'of the first EQ6 problem. "Fluid 1" is the aqueous s
     $olution at the'
      utitl1(15) = 'end of the current EQ6 problem, and the ES is the co
     $rresponding ES.'
c
      if (nprob .eq. 1) then
        utitl1(17) = '  Warning: This EQ6 input file is set up to mix Fl
     $"uid 2" with an ES'
        utitl1(18) = 'containing the identical aqueous solution. You may
     $ wish to change'
        utitl1(19) = 'the bottom half of this file so that this fluid mi
     $xes with something else.'
      endif
c
c     Note: "Fluid 2" as a special reactant is defined in EQ6/path.f.
c
      nrct = 1
      nsrt = 1
      ureac(1) = ureac1
      jcode(1) = 2
      jreac(1) = 0
      morr(1) = morrw1
      modr(1) = 0.
      nrk(1,1) = 1
      rkb(1,1,1) = 1.0
c
      vreac(1) = 1.0
c
      iesrti(1) = iesrt1
      do n = 1,iesrt1
        uesri(n,1) = uesr1(n)
        cesri(n,1) = cesr1(n)
      enddo
c
      ibsrti(1) = ibsrt1
      do n = 1,ibsrt1
        ubsri(n,1) = ubsr1(n)
        cbsri(n,1) = cbsr1(n)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
