      subroutine vazpck(nat,natmax,nazt,naztmx,noutpt,nttyo,
     $ uaqsp,uazp)
c
c     Validate the aqueous species names used to specify hard core
c     diameters. Write a note if such a name does not appear on the
c     main list of aqueous species.
c
c     This subroutine is called by:
c
c       EQPT/eqpt.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       nat    = the number of aqueous species
c       nazt   = the number of specified hard core diameters
c       uaqsp  = array of names of aqueous species
c       uazp   = array of aqueous species names used to specify
c                  hard core diamters on the data file
c
c     Principal output:
c
c       None
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer natmax,naztmx
c
      integer noutpt,nttyo
c
      integer nat,nazt
c
      character*24 uaqsp(natmax),uazp(naztmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      character*24 unam
c
      integer j2,na,naz
c
      integer ilnobl
c
c-----------------------------------------------------------------------
c
c     Loop over all aqueous species names used to specify hard core
c     diameters.
c
      do naz = 1,nazt
        unam = uazp(naz)
        if (unam(1:7) .ne. '<blank>') then
c
c         Check that each such name appears on the main list of
c         aqueous species.
c
          do na = 1,nat
            if (unam(1:24) .eq. uaqsp(na)(1:24)) go to 100
          enddo
c
c         Did not find this name on the main list of aqueous species.
c
          j2 = ilnobl(unam)
          write (noutpt,1000) unam(1:j2)
          write (nttyo,1000) unam(1:j2)
 1000     format(/' * Note - (EQPT/vazpck) A hard core diameter is',
     $    ' specified on the',/7x,'data file for ',a,', but there is',
     $    ' no species block for this',/7x,'in the aqueous species',
     $    ' superblock. Therefore, this species will',/7x,'not appear',
     $    ' in any computed model and the hard core diameter value',
     $    /7x,'will not be used. It is suggested to remove the',
     $    ' hard core diameter',/7x,'in question or to add the',
     $    ' corresponding species block.')
c
  100     continue
        endif
      enddo
c
      end
