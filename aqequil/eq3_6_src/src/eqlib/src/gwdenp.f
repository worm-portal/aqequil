      subroutine gwdenp(adwipp,bdwipp,jcsort,mlmrra,mosp,mrmlra,
     $ mwtsp,narn1,narn2,nstmax,qdwipp,rhoc,rhowc,tdsgks,tdsglw,
     $ tdspkc,tdsplc,tempc,vosol,wfh2o,wftds,wkgwi,woh2o,
     $ wosol,wotds)
c
c     This subroutine gets the weights of solute, total dissolved
c     solutes, and aqueous solution, and obtains the density of
c     the aqueous solution by evaluating the WIPP brine density
c     model. It also computes some related parameters, such
c     as the solution volume and the molarity/molality ratio.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/path.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       adwipp = 'a' constant in the WIPP brine density model
c       bdwipp = 'b' constant in the WIPP brine density model
c       jcsort = array of species indices in order of increasing
c                  concentration, organized into blocks, one per phase
c       mosp   = numbers of moles of species
c       mwtsp  = molecular weights of species
c       mwtsp  = molecular weights of species
c
c     Principal output:
c
c       mlmrra = molality/molarity ratio
c       mrmlra = molarity/molality ratio
c       qdwipp = .true. if results from the WIPP brine density model
c                  are available
c       rhoc   = the solution density (g/mL), calculated from the WIPP
c                  brine density model
c       rhowc  = the solution density (g/L), calculated from the WIPP
c                  brine density model
c       tdsglw = TDS (total dissolved solutes) (g/L), calculated
c                  from the WIPP brine density model
c       tdsgks = TDS (total dissolved solutes) (g/kg.sol)
c       tdspkc = TDS (total dissolved solutes) (mg/kg.sol)
c       tdsplc = TDS (total dissolved solutes) (mg/L), calculated
c                  from the WIPP brine density model
c       vosol  = Volume (L) of aqueous solution
c       wfh2o  = weight (mass) fraction of water in aqueous solution
c       wftds  = weight (mass) fraction of solutes in aqueous solution
c       woh2o  = weight (mass) of solvent water
c       wosol  = weight (mass) of aqueous solution
c       wotds  = weight (mass) of total dissolved solutes
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nstmax
c
      integer jcsort(nstmax)
c
      integer narn1,narn2
c
      logical qdwipp
c
      real(8) mosp(nstmax),mwtsp(nstmax)
c
      real(8) adwipp,bdwipp,mlmrra,mrmlra,rhoc,rhowc,tdsgks,tdsglw,
     $ tdspkc,tdsplc,tempc,vosol,wfh2o,wftds,wkgwi,woh2o,wosol,wotds
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
c     None
c
c-----------------------------------------------------------------------
c
c     Calculate the weights (actually, the masses) of total dissolved
c     solutes, the solvent, and the solution, and also some related
c     quantities.
c
c     Get the total dissolved solute mass (wotds, g).
c
      call ctds(jcsort,mosp,mwtsp,narn1,narn2,nstmax,wotds)
c
c     Get the solvent mass (fixed at 1000 grams in EQ3NR, variable
c     in EQ6). The following works for either code.
c
      woh2o = mosp(narn1)*mwtsp(narn1)
c
c     Get the mass of aqueous solution.
c
      wosol =  woh2o + wotds
c
c     Get the weight fraction of solvent and the weight fraction
c     of total dissolved solutes.
c
      wfh2o = woh2o/wosol
      wftds = wotds/wosol
c
c     Get the inverse of the weight of solvent, in kg.
c
      wkgwi = 1000./woh2o
c
c     Get the TDS (g/kg.solution)
c
      tdsgks = 1000.*wftds
c
c     Get the TDS (mg/kg.solution)
c
      tdspkc = 1000.*tdsgks
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the solution density (rhowc) from the WIPP brine
c     density model if the temperature lies in the range 20-30 C.
c     This model is a fit to NaCl solution data having the form:
c
c       density (g/L) = a + b x TDS (g/L)
c
c     where a = 1000.96 and b = 0.639963. Note that this formulation
c     requires knowing TDS in g/L (tdsglw), which prevents a direct
c     evaluation using the above form because the density is not
c     available in the appropriate units. Also calculate the molarity/
c     molality conversion ratio (mrmlra, molarity = mrmlra x molality)
c     and its inverse.
c
      qdwipp = .false.
      rhowc = 0.
      vosol = 0.
      tdsglw = 0.
      mrmlra = 0.
      mlmrra = 0.
      rhoc = 0.
      tdsplc = 0.
c
      if (tempc.ge.20. .and. tempc.le.30.) then
c
c       Note that wfh2o corresponds to Nw, the "solvent fraction"
c       in the density calculations writeup. Here wftds corresponds
c       to "(1 - Nw)".
c
        qdwipp = .true.
        rhowc = adwipp/(1.0 - bdwipp*wftds)
        vosol = wosol/rhowc
        tdsglw = rhowc*wftds
        mrmlra = 0.001*wfh2o*rhowc
        if (mrmlra .gt. 0.) then
          mlmrra = 1./mrmlra
        else
          mlmrra = 1.e+38
        endif
        rhoc = 0.001*rhowc
        tdsplc = 1000.*tdsglw
      endif
c
  999 continue
      end
