      subroutine prtvpa(abar,acfw,acfwlg,actw,actwlg,a3bar,fje,
     $ fjest,fo2,fo2lg,fxi,fxist,iopg,mlmrra,mrmlra,nopgmx,noutpt,
     $ osc,oscst,qrho,rhoc,rhowc,sigmam,sigmst,tdsglw,tdspkc,
     $ tdsplc,vosol,wfh2o,wftds,woh2o,wosol,wotds,xbarw,xbrwlg)
c
c     This subroutine prints a table of various aqueous solution
c     parameters.
c
c     This subroutine is called by:
c
c       EQ3NR/scripx.f
c       EQ6/scripz.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       abar   = average hard core diameter
c       acfw   = activity coefficient of water
c       acfwlg = log activity coefficient of water
c       actw   = activity of water
c       actwlg = log activity of water
c       a3bar  = average cubed hard core diameter
c       fje    = the ionic asymmetry (the 3rd-order electrostatic
c                  moment function J)
c       fjest  = stoichiometric ionic asymmetry
c       fo2    = oxygen fugacity
c       fo2lg  = log oxygen fugacity
c       fxi    = the ionic strength (the 2nd-order electrostatic
c                  moment function I)
c       fxist  = stoichiometric ionic strength
c       mlmrra = molality/molarity ratio
c       mrmlra = molarity/molality ratio
c       osc    = osmotic coefficient
c       oscst  = stoichiometric osmotic coefficient
c       qrho   = .true. if a solution density value is available
c       rhoc   = the solution density (g/mL)
c       rhowc  = the solution density (g/L)
c       sigmam = sum of solute molalities
c       sigmst = sum of stoichiometric solute molalities
c       tdsglw = TDS (total dissolved solutes) (g/L)
c       tdspkc = TDS (total dissolved solutes) (mg/kg.sol)
c       tdsplc = TDS (total dissolved solutes) (mg/L)
c       vosol  = Volume (L) of aqueous solution
c       wfh2o  = weight (mass) fraction of water in aqueous solution
c       wftds  = weight (mass) fraction of solutes in aqueous solution
c       woh2o  = weight (mass) of solvent water
c       wosol  = weight (mass) of aqueous solution
c       wotds  = weight (mass) of total dissolved solutes
c       xbarw  = mole fraction of water in aqueous solution
c       xbrwlg = log mole fraction of water in aqueous solution
c
c     Principal output:
c
c       None
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nopgmx
c
      integer noutpt
c
      integer iopg(nopgmx)
c
      logical qrho
c
      real(8) abar,acfw,acfwlg,actw,actwlg,a3bar,fje,fjest,fo2,fo2lg,
     $ fxi,fxist,osc,oscst,sigmam,sigmst,vosol,wfh2o,wftds,woh2o,wosol,
     $ wotds,xbarw,xbrwlg
c
      real(8) mlmrra,mrmlra,rhoc,rhowc,tdsglw,tdspkc,tdsplc
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j1,j2,j3
c
      integer ilnobl
c
      character(len=16) ux16a,ux16b,ux16c
c
c-----------------------------------------------------------------------
c
      if (fo2lg .gt. -99999.) then
        write (ux16a,'(1pg12.5)') fo2
        call lejust(ux16a)
        j1 = ilnobl(ux16a)
        write (ux16b,'(1pg12.5)') fo2lg
        call lejust(ux16b)
        j2 = ilnobl(ux16b)
        write (noutpt,1000) ux16a(1:j1),ux16b(1:j2)
 1000   format(//20x,'Oxygen fugacity= ',a,' bars',
     $ /16x,'Log oxygen fugacity= ',a)
      else
        write (noutpt,1010)
 1010   format(/8x,'Oxygen fugacity is not defined.')
      endif
c
      write (ux16a,'(1pg12.5)') actw
      call lejust(ux16a)
      j1 = ilnobl(ux16a)
      write (ux16b,'(1pg12.5)') actwlg
      call lejust(ux16b)
      j2 = ilnobl(ux16b)
      write (noutpt,1020) ux16a(1:j1),ux16b(1:j2)
 1020 format(/18x,'Activity of water= ',a,
     $ /14x,'Log activity of water= ',a)
c
      write (ux16a,'(1pg12.5)') xbarw
      call lejust(ux16a)
      j1 = ilnobl(ux16a)
      write (ux16b,'(1pg12.5)') xbrwlg
      call lejust(ux16b)
      j2 = ilnobl(ux16b)
      write (noutpt,1030) ux16a(1:j1),ux16b(1:j2)
 1030 format(/13x,'Mole fraction of water= ',a,
     $ /9x,'Log mole fraction of water= ',a)
c
      write (ux16a,'(1pg12.5)') acfw
      call lejust(ux16a)
      j1 = ilnobl(ux16a)
      write (ux16b,'(1pg12.5)') acfwlg
      call lejust(ux16b)
      j2 = ilnobl(ux16b)
      write (noutpt,1040) ux16a(1:j1),ux16b(1:j2)
 1040 format(/6x,'Activity coefficient of water= ',a,
     $ /2x,'Log activity coefficient of water= ',a)
c
      write (ux16a,'(1pg12.5)') osc
      call lejust(ux16a)
      j1 = ilnobl(ux16a)
      write (ux16b,'(1pg12.5)') oscst
      call lejust(ux16b)
      j2 = ilnobl(ux16b)
      write (noutpt,1050) ux16a(1:j1),ux16b(1:j2)
 1050 format(/16x,'Osmotic coefficient= ',a,
     $ /1x,'Stoichiometric osmotic coefficient= ',a)
c
      write (ux16a,'(1pg12.5)') sigmam
      call lejust(ux16a)
      j1 = ilnobl(ux16a)
      write (ux16b,'(1pg12.5)') sigmst
      call lejust(ux16b)
      j2 = ilnobl(ux16b)
      write (noutpt,1060) ux16a(1:j1),ux16b(1:j2)
 1060 format(/18x,'Sum of molalities= ',a,
     $ /3x,'Sum of stoichiometric molalities= ',a)
c
      write (ux16a,'(1pg12.5)') fxi
      call lejust(ux16a)
      j1 = ilnobl(ux16a)
      write (ux16b,'(1pg12.5)') fxist
      call lejust(ux16b)
      j2 = ilnobl(ux16b)
      write (noutpt,1070) ux16a(1:j1),ux16b(1:j2)
 1070 format(/17x,'Ionic strength (I)= ',a,' molal',
     $ /6x,'Stoichiometric ionic strength= ',a,' molal')
c
      write (ux16a,'(1pg12.5)') fje
      call lejust(ux16a)
      j1 = ilnobl(ux16a)
      write (ux16b,'(1pg12.5)') fjest
      call lejust(ux16b)
      j2 = ilnobl(ux16b)
      write (noutpt,1080) ux16a(1:j1),ux16b(1:j2)
 1080 format(/16x,'Ionic asymmetry (J)= ',a,' molal',
     $ /5x,'Stoichiometric ionic asymmetry= ',a,' molal')
c
      if (iopg(1) .eq. 2) then
        write (ux16a,'(1pg12.5)') abar
        call lejust(ux16a)
        j1 = ilnobl(ux16a)
        write (ux16b,'(1pg12.5)') a3bar
        call lejust(ux16b)
        j2 = ilnobl(ux16b)
        write (noutpt,1090) ux16a(1:j1),ux16b(1:j2)
 1090   format(/31x,'Abar = ',a,/30x,'A3bar = ',a)
      endif
c
      write (ux16a,'(1pg12.5)') woh2o
      call lejust(ux16a)
      j1 = ilnobl(ux16a)
      write (ux16b,'(1pg12.5)') wotds
      call lejust(ux16b)
      j2 = ilnobl(ux16b)
      write (ux16c,'(1pg12.5)') wosol
      call lejust(ux16c)
      j3 = ilnobl(ux16c)
      write (noutpt,1100) ux16a(1:j1),ux16b(1:j2),ux16c(1:j3)
 1100 format(/23x,'Solvent mass= ',a,' g',
     $ /17x,'Solutes (TDS) mass= ',a,' g',
     $ /14x,'Aqueous solution mass= ',a,' g')
c
      if (qrho) then
        write (ux16a,'(1pg12.5)') vosol
        call lejust(ux16a)
        j1 = ilnobl(ux16a)
        write (noutpt,1110) ux16a(1:j1)
 1110   format(/12x,'Aqueous solution volume= ',a,' L')
      endif
c
      write (ux16a,'(1pg12.5)') wfh2o
      call lejust(ux16a)
      j1 = ilnobl(ux16a)
      write (ux16b,'(1pg12.5)') wftds
      call lejust(ux16b)
      j2 = ilnobl(ux16b)
      write (noutpt,1120) ux16a(1:j1),ux16b(1:j2)
 1120 format(/19x,'Solvent fraction= ',a,' kg.H2O/kg.sol',
     $ /20x,'Solute fraction= ',a,' kg.tds/kg.sol')
c
      write (ux16a,'(1pg12.5)') tdspkc
      call lejust(ux16a)
      j1 = ilnobl(ux16a)
      write (noutpt,1140) ux16a(1:j1)
 1140 format(/6x,'Total dissolved solutes (TDS)= ',a,' mg/kg.sol')
c
      if (qrho) then
        write (ux16a,'(1pg12.5)') tdsplc
        call lejust(ux16a)
        j1 = ilnobl(ux16a)
        write (ux16b,'(1pg12.5)') tdsglw
        call lejust(ux16b)
        j2 = ilnobl(ux16b)
        write (noutpt,1210) ux16a(1:j1),ux16b(1:j2)
 1210   format(32x,'TDS= ',a,' mg/L',/32x,'TDS= ',a,' g/L')
c
        write (ux16a,'(1pg12.5)') rhoc
        call lejust(ux16a)
        j1 = ilnobl(ux16a)
        write (ux16b,'(1pg12.5)') rhowc
        call lejust(ux16b)
        j2 = ilnobl(ux16b)
        write (noutpt,1220) ux16a(1:j1),ux16b(1:j2)
 1220   format(/19x,'Solution density= ',a,' g/mL',
     $  /19x,'Solution density= ',a,' g/L')
c
        write (ux16a,'(1pg12.5)') mrmlra
        call lejust(ux16a)
        j1 = ilnobl(ux16a)
        write (ux16b,'(1pg12.5)') mlmrra
        call lejust(ux16b)
        j2 = ilnobl(ux16b)
        write (noutpt,1240) ux16a(1:j1),ux16b(1:j2)
 1240   format(/18x,'Molarity/molality= ',a,' kg.H2O/L',
     $  /18x,'Molality/molarity= ',a,' L/kg.H2O')
      endif
c
      write (noutpt,1300)
 1300 format(1x)
c
      end
