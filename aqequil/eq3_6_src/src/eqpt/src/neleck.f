      subroutine neleck(nct,nctmax,nerr,noutpt,nttyo,uelem)
c
c     Check the names of chemical elements for uniqueness.
c
c     This subroutine is called by:
c
c       EQPT/eqpt.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
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
      integer nctmax
c
      integer noutpt,nttyo
c
      integer nct,nerr
c
      character(len=8) uelem(nctmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      character*8 unam,ux8
c
      integer i,ilist,j,j2,j3,k,ncount
c
      integer ilnobl
c
c-----------------------------------------------------------------------
c
      ilist = 0
c
c     Check each name for uniqueness.
c
      do i = 1,nct
        unam = uelem(i)
        ncount = 1
        if (unam(1:7) .ne. '<blank>') then
          do j = i + 1,nct
            if (unam(1:8) .eq. uelem(j)(1:8)) then
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
                  if (unam(1:8) .eq. uelem(k)(1:8)) then
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
        if (ncount .gt. 1) then
c
c         Have a duplicate entry or entries.
c
          ilist = ilist + 1
c
          if (ilist .eq. 1) then
c
c           Write error message header.
c
            write (noutpt,1000)
            write (nttyo,1000)
 1000       format(/' * Error - (EQPT/neleck) The following chemical',
     $      ' elements have multiple',/7x,'entries in the elements',
     $      ' block:',/)
          endif
c
          write(ux8,'(i5)') ncount
          call lejust(ux8)
          j3 = ilnobl(ux8)
          j2 = ilnobl(unam)
          write (noutpt,1010) unam(1:j2),ux8(1:j3)
          write (nttyo,1010) unam(1:j2),ux8(1:j3)
 1010     format(9x,a,' has ',a,' entries')
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      nerr = nerr + ilist
c
      end
