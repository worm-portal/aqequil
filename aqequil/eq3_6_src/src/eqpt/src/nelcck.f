      subroutine nelcck(nctmax,ncts,nentei,nerr,qdupes,uessi)
c
c     Check the elemental composition of a species to ensure that
c     each chemical element name appearing in the composition is unique.
c
c     This subroutine is called by:
c
c       EQPT/elesck.f
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
      integer nctmax
c
      integer nentei(nctmax)
c
      integer ncts,nerr
c
      logical qdupes
c
      character(len=8) uessi(nctmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      character*8 unam
c
      integer i,j,k,ncount
c
c-----------------------------------------------------------------------
c
      qdupes = .false.
c
c     Check each name for uniqueness.
c
      do i = 1,ncts
        unam = uessi(i)
        ncount = 1
        if (unam(1:7) .ne. '<blank>') then
          do j = i + 1,ncts
            if (unam(1:8) .eq. uessi(j)(1:8)) then
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
                  if (unam(1:8) .eq. uessi(k)(1:8)) then
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
        nentei(i) = ncount
        if (ncount .gt. 1) qdupes = .true.
      enddo
c
      end
