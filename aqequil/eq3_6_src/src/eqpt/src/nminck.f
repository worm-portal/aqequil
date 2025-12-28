      subroutine nminck(nerr,nmt,nmtmax,noutpt,nttyo,uminsp)
c
c     Check the names of pure minerals for uniqueness.
c
c     This subroutine is called by:
c
c       EQPT/eqpt.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       uminsp = array containing the names of pure minerals
c
c     Principal output:
c
c       nerr   = incremented error counter
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer noutpt,nttyo
c
      integer nmtmax
c
      integer nmt,nerr
c
      character(len=24) uminsp(nmtmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      character*24 unam
      character*8 ux8
c
      integer i,ilist,j,j2,j3,k,ncount
c
      integer ilnobl
c
c-----------------------------------------------------------------------
c
c     Note: ilist = the number of species with duplicate species blocks
c     present on the data file.
c
      ilist = 0
c
c     Check each name for uniqueness.
c
      do i = 1,nmt
        unam = uminsp(i)
        ncount = 1
        if (unam(1:7) .ne. '<blank>') then
          do j = i + 1,nmt
            if (unam(1:24) .eq. uminsp(j)(1:24)) then
c
c             Have found a duplicate block moving down the list from
c             the block now being tested for uniqueness. For the first
c             such duplicate block only, make sure that there is not a
c             duplicate block moving up the list from the block now
c             being tested. If there is such a block, then the
c             duplication in question has been noted previously and
c             should not be noted again.
c
              if (ncount .eq. 1) then
                do k = 1,i - 1
                  if (unam(1:24) .eq. uminsp(k)(1:24)) then
c
c                   Have found a duplicate block preceding the block
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
        if (ncount .gt. 1) then
c
c         Have a duplicate block or blocks.
c
          ilist = ilist + 1
c
          if (ilist .eq. 1) then
c
c           Write error message header.
c
            write (noutpt,1000)
            write (nttyo,1000)
 1000       format(/' * Error - (EQPT/nminck) The following pure',
     $      ' minerals have multiple',/7x,'species blocks:',/)
          endif
c
          write (ux8,'(i5)') ncount
          call lejust(ux8)
          j3 = ilnobl(ux8)
          j2 = ilnobl(unam)
          write (noutpt,1010) unam(1:j2),ux8(1:j3)
          write (nttyo,1010) unam(1:j2),ux8(1:j3)
 1010     format(9x,a,' has ',a,' blocks')
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      nerr = nerr + ilist
c
      end
