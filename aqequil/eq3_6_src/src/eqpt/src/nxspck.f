      subroutine nxspck(iktmax,issot,nerr,noutpt,nttyo,nxt,nxtmax,
     $ ussoph,ussosp)
c
c     Check the names of solid solution end-members for uniqueness
c     within each solid solution.
c
c     This subroutine is called by:
c
c       EQPT/eqpt.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ussoph = array containing the names of solid solutions
c       ussosp = array containing the names of solid solution
c                  end-members
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
      integer iktmax,nxtmax
c
      integer nxt,nerr
c
      integer issot(nxtmax)
c
      character(len=24) ussoph(nxtmax)
      character(len=24) ussosp(iktmax,nxtmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      character*24 unam
      character*8 ux8
c
      integer i,ikt,ilist,j,j2,j3,j4,k,nx,ncount
c
      integer ilnobl
c
c-----------------------------------------------------------------------
c
c     Loop over all solid solution phases.
c
      do nx = 1,nxt
c
        ikt = issot(nx)
c
c       Note: ilist = the number of end-members that are duplicated
c       for a given solid solution.
c
        ilist = 0
c
c       Check each name for uniqueness.
c
        do i = 1,ikt
          unam = ussosp(i,nx)
          ncount = 1
          if (unam(1:7) .ne. '<blank>') then
            do j = i + 1,ikt
              if (unam(1:24) .eq. ussosp(j,nx)(1:24)) then
c
c               Have found a duplicate entry moving down the list from
c               the entry now being tested for uniqueness. For the first
c               such duplicate entry only, make sure that there is not a
c               duplicate entry moving up the list from the entry now
c               being tested. If there is such a entry, then the
c               duplication in question has been noted previously and
c               should not be noted again.
c
                if (ncount .eq. 1) then
                  do k = 1,i - 1
                    if (unam(1:24) .eq. ussosp(k,nx)(1:24)) then
c
c                     Have found a duplicate entry preceding the entry
c                     now being tested for uniqueness.
c
                      go to 100
                    endif
                  enddo
                endif
c
c               Have found a duplication that has not been previously
c               noted.
c
                ncount = ncount + 1
              endif
            enddo
  100       continue
          endif
c
          if (ncount .gt. 1) then
c
c           Have a duplicate entry or entries.
c
            ilist = ilist + 1
c
            if (ilist .eq. 1) then
c
c             Write error message header.
c
              j2 = ilnobl(ussoph(nx))
              write (noutpt,1000) ussoph(nx)(1:j2)
              write (nttyo,1000) ussoph(nx)(1:j2)
 1000         format(/' * Error - (EQPT/nxspck) The following',
     $        ' end-members are listed more than',/7x,'once for solid',
     $        ' solution ',a,':',/)
            endif
c
            write (ux8,'(i5)') ncount
            call lejust(ux8)
            j4 = ilnobl(ux8)
            j3 = ilnobl(unam)
            write (noutpt,1010) unam(1:j3),ux8(1:j4)
            write (nttyo,1010) unam(1:j3),ux8(1:j4)
 1010       format(9x,a,' has ',a,' entries')
          endif
        enddo
c
        nerr = nerr + ilist
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
