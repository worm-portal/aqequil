      integer function jesite(jern1,jern2,jetmax,jgext,ne,netmax,ns)
c
c     This subroutine returns the index of the site which contains
c     the ns-th species, a species of the ne-th generic ion exchanger
c     phase. If a match is not found, a zero value is returned.
c
c     This subroutine is called by:
c
c       None
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       jern1  = array giving the starts of the species ranges of sites
c                of generic ion exchanger phases
c       jern2  = array giving the ends of the species ranges of sites
c                of generic ion exchanger phases
c       jgext  = array of the number of sites in generic ion exchanger
c                  phases
c       ne     = the index of the desired generic ion exchanger phase
c       ns     = the index of the desired species
c
c     Output:
c
c       jesite = the index of the site containing the species whose
c                  index is ns, and which is a species in the ne-th
c                  generic ion exchanger phase
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer jetmax,netmax
c
      integer jern1(jetmax,netmax),jern2(jetmax,netmax),jgext(netmax)
      integer ne,ns
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer je,nrr1,nrr2
c
c-----------------------------------------------------------------------
c
      jesite = 0
      do je = 1,jgext(ne)
        nrr1 = jern1(je,ne)
        nrr2 = jern2(je,ne)
        if (ns .le. nrr2) then
          if (ns .ge. nrr1) then
            jesite = je
            go to 999
          endif
        endif
      enddo
c
  999 continue
      end
