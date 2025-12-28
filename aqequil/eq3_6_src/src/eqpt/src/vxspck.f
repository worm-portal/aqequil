      subroutine vxspck(iktmax,issot,nerr,nmt,nmtmax,noutpt,nttyo,
     $ nxt,nxtmax,uminsp,ussoph,ussosp)
c
c     Validate the names of solid solution end-members. Make sure that
c     these names appear on the list of pure minerals.
c
c     This subroutine is called by:
c
c       EQPT/eqpt.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       nxt    = the number of solid solutions
c       issot  = array of numbers of end-members of solid solutions
c       ussoph = array of names of solid solutions
c       ussosp = array of names of end-members of solid solutions
c       uminsp = array of pure mineral names
c
c     Principal output:
c
c       nerr   = error counter
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer iktmax,nmtmax,nxtmax
c
      integer noutpt,nttyo
c
      integer nerr,nmt,nxt
c
      integer issot(nxtmax)
c
      character(len=24) uminsp(nmtmax),ussoph(nxtmax),
     $ ussosp(iktmax,nxtmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      character*24 unam
c
      integer i,ikt,j2,j3,kerr,nm,nx
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
        kerr = 0
c
c       Check that each end-member appears on the list of pure minerals.
c
        do i = 1,ikt
          unam = ussosp(i,nx)
          if (unam(1:7) .ne. '<blank>') then
c
c           Check that each end-member appears on the list of
c           pure minerals.
c
            do nm = 1,nmt
              if (unam(1:24) .eq. uminsp(nm)(1:24)) go to 100
            enddo
c
c           Did not find this end-member on the list of pure minerals.
c
            kerr = kerr + 1
            nerr = nerr + 1
c
            if (kerr .eq. 1) then
c
c             Write a header for the solid solution.
c
              j2 = ilnobl(ussoph(nx))
              write (noutpt,1000) ussoph(nx)(1:j2)
              write (nttyo,1000) ussoph(nx)(1:j2)
 1000         format(/' * Error - (EQPT/vxspck) The following',
     $        ' end-members of solid solution',/7x,a,' are not present',
     $        ' on the data file as pure minerals:',/)
            endif
c
            j3 = ilnobl(unam)
            write (noutpt,1010) unam(1:j3)
            write (nttyo,1010) unam(1:j3)
 1010       format(9x,a)
c
  100       continue
          endif
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
