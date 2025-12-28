      subroutine srchn(nrn1a,nrn2a,ns,nstmax,unam,uspeca)
c
c     This subroutine matches the species name unam with the
c     corresponding entry in the nrn1a-th through nrn2a-th range of
c     the species name array uspeca. Only the first 24 characters are
c     compared (uspeca has 48 characters, unam only 24). This
c     subroutine returns the species index ns. If there is no match,
c     ns is returned with a value of 0.
c
c     This subroutine is called by:
c
c       EQLIB/inbdot.f
c       EQLIB/inupt.f
c       EQLIB/srchne.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       nrn1a  = start of the range of species to search
c       nrn2a  = end of the range of species to search
c       unam   = name of the species whose index is to be found
c       uspeca = array of species names
c
c     Principal output:
c
c       ns     = index of the species whose name is unam
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
      integer nrn1a,nrn2a,ns
c
      character*48 uspeca(nstmax)
      character*24 unam
c
c----------------------------------------------------------------------
c
c     Local variable declarations.
c
c
c----------------------------------------------------------------------
c
      do ns = nrn1a,nrn2a
        if (uspeca(ns)(1:24) .eq. unam(1:24)) go to 999
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     No match was found.
c
      ns = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
