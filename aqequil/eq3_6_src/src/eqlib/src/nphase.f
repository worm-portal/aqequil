      integer function nphase(ncmpr,npt,nptmax,ns)
c
c     This subroutine returns the index of the phase which contains
c     the ns-th species. If a match is not found, a zero value is
c     returned.
c
c     This subroutine is called by:
c
c       EQ6/mincsp.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       ncmpr  = array giving the start and the end of the range of the
c                species belonging to a given phase
c       npt    = the number of phases
c       ns     = the index of the desired species
c
c     Output:
c
c       nphase = the index of the phase containing the ns-th species
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nptmax
c
      integer ncmpr(2,nptmax)
      integer npt,ns
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer np,nr1,nr2
c
c-----------------------------------------------------------------------
c
      nphase = 0
      do np = 1,npt
        nr1 = ncmpr(1,np)
        nr2 = ncmpr(2,np)
        if (ns .le. nr2) then
          if (ns .ge. nr1) then
            nphase = np
            go to 999
          endif
        endif
      enddo
c
  999 continue
      end
