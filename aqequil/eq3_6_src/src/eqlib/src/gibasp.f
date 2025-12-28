      subroutine gibasp(cgexj,cjbasp,iern1,ixbasp,jern1,jern2,
     $ jetmax,jgext,narn1,narn2,nbasp,nbt,nbtmax,nern1,nern2,
     $ netmax,nphasx,nstmax)
c
c     This subroutine sets up the ixbasp and cjbasp arrays. The former
c     is a flag array, each member of which denotes whether the
c     thermodynamic activity of the corresponding basis species is
c     defined in terms of molality (= 0) or mole fraction (= 1). The
c     cjbasp array contains any site stoichiometric factor associated
c     with a given basis species. If there is no such factor, then
c     an element of cjbasp may be assigned a value of unity.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ3NR/arrset.f
c       EQ6/eqphas.f
c       EQ6/eqshel.f
c       EQ6/optmzr.f
c       EQ6/path.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       cgexj  = array of site stoichiometric factors for generic ion
c                  exchanger species
c       jgext  = array giving the number of sites for a generic ion
c                  exchanger species.
c       nbasp  = array containing the species indices of the basis
c                  species
c       nbt    = the number of basis species
c       nphasx = array giving the phase to which a species belongs
c
c     Principal output:
c
c       cjbasp = array of site stoichiometric factors corresponding to
c                  sites occupied by basis species, if any
c       ixbasp = array of flag switches indicating whether the
c                  thermodynamic activity of a basis species is based
c                  on molality (= 0) or mole fraction (= 1)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer jetmax,nbtmax,netmax,nstmax
c
      integer ixbasp(nbtmax),jern1(jetmax,netmax),jern2(jetmax,netmax),
     $ jgext(netmax),nbasp(nbtmax),nphasx(nstmax)
c
      integer iern1,narn1,narn2,nbt,nern1,nern2
c
      real*8 cjbasp(nbtmax),cgexj(jetmax,netmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer je,nb,ne,np,nrr1,nrr2,ns,nss
c
c-----------------------------------------------------------------------
c
      do nb = 1,nbt
        ixbasp(nb) = 0
        cjbasp(nb) = 1.0
c
        ns = nbasp(nb)
        if (ns .eq. narn1) then
          ixbasp(nb) = 1
        elseif (ns.ge.nern1 .and. ns.le.nern2) then
          np = nphasx(ns)
          ne = np - iern1 + 1
          ixbasp(nb) = 1
          do je = 1,jgext(ne)
            nrr1 = jern1(je,ne)
            nrr2 = jern2(je,ne)
            do nss = nrr1,nrr2
              if (ns .eq. nss) then
                cjbasp(nb) = cgexj(je,ne)
              endif
            enddo
          enddo
        endif
      enddo
c
  999 continue
      end
