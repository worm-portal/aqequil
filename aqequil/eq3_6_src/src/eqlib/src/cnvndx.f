      subroutine cnvndx(nlim1,nlim1a,nlim2,nlim2a,nsmap,nstmax,ntot)
c
c     This subroutine finds the first and last species in a reduced
c     range corresponding to an original set defined by the species
c     index limits nlim1a, nlim2a. Note that this is not a
c     straightforward mapping of such as "nlim1 = nsmap(nlim1a),"
c     because the species whose original index is nlim1a may have been
c     eliminated in the reduction.
c
c     This subroutine is called by:
c
c       EQLIB/cmpdat.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       nlim1a = starting index of the range for the original set
c       nlim2a = ending index of the range for the original set
c       nsmap  = pointer array mapping species indices from the
c                  original set to the reduced set
c
c     Principal output:
c
c       nlim1  = starting index of the range for the reduced set
c       nlim2  = ending index of the range for the reduced set
c       ntot   = number of species in the range for the reduced set
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
      integer nsmap(nstmax)
      integer nlim1,nlim1a,nlim2,nlim2a,ntot
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j
c
c-----------------------------------------------------------------------
c
c     Search for the first species that is in the list.
c
      do j = nlim1a,nlim2a,1
        if (nsmap(j) .ne. 0) go to 120
      enddo
c
c     None found, error.
c
      nlim1 = 0
      nlim2 = -1
      ntot = 0
      go to 999
c
 120  continue
c
c     Save the first species.
c
      nlim1 = nsmap(j)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Search for last species that is in the list.
c
      do j = nlim2a,nlim1a,-1
         if (nsmap(j) .ne. 0) go to 140
      enddo
c
 140  continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Save the last species.
c
      nlim2 = nsmap(j)
      ntot = nlim2 - nlim1 + 1
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
