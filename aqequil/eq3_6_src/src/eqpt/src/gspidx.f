      subroutine gspidx(ier,n,nat,natmax,uaqsp,unams)
c
c     Get the index (n) of the aqueous species whose name is unams.
c
c     This subroutine is called by:
c
c       EQPT/rdpz2.f
c       EQPT/rdpz3.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       nat    = the number of aqueous species
c       uaqsp  = array of names of aqueous species
c       unams  = the name of an aqueous species whose index n is
c                  desired
c
c     Principal output:
c
c       ier    = error flag
c       n      = the index of the species unams in the uaqsp array
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer natmax
c
      integer ier,n,nat
c
      character*24 uaqsp(natmax)
      character*24 unams
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
c
c-----------------------------------------------------------------------
c
      ier = 0
      do n = 1,nat
        if (unams(1:24) .eq. uaqsp(n)(1:24)) then
          go to 100
        endif
      enddo
c
c     The species was not found. Set the error flag.
c
      ier = 1
      n = 0
c
  100 continue
c
      end
