      subroutine nrxnck(nbtmx1,ndrsts,nentri,nerr,qduprs,udrsi)
c
c     Check the associated raection of a species to ensure that
c     each species name appearing in the reaction is unique.
c
c     This subroutine is called by:
c
c       EQPT/rxnsck.f
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
      integer nbtmx1
c
      integer nentri(nbtmx1)
c
      integer ndrsts,nerr
c
      logical qduprs
c
      character(len=24) udrsi(nbtmx1)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      character*24 unam
c
      integer i,j,k,ncount
c
c-----------------------------------------------------------------------
c
      qduprs = .false.
c
c     Check each name for uniqueness.
c
      do i = 1,ndrsts
        unam = udrsi(i)
        ncount = 1
        if (unam(1:7) .ne. '<blank>') then
          do j = i + 1,ndrsts
            if (unam(1:24) .eq. udrsi(j)(1:24)) then
c
c             Have found a duplicate entry moving down the list from
c             the entry now being tested for uniqueness. For the first
c             such duplicate entry only, make sure that there is not a
c             duplicate entry moving up the list from the entry now
c             being tested. If there is such a entry, then the
c             duplication in question has been noted previously and
c             should not be noted again.
c
              if (ncount .eq. 1) then
                do k = 1,i - 1
                  if (unam(1:24) .eq. udrsi(k)(1:24)) then
c
c                   Have found a duplicate entry preceding the entry
c                   now being tested for uniqueness.
c
                    go to 100
                  endif
                enddo
              endif
c
c             Have found a duplication that has not been previously
c             noted.
c
              ncount = ncount + 1
            endif
          enddo
  100     continue
        endif
c
        nentri(i) = ncount
        if (ncount .gt. 1) qduprs = .true.
      enddo
c
      end
