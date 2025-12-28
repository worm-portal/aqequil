      subroutine srch22(jpair,unam1,unam2,upair,npx2mx,npx2t)
c
c     This subroutine searches for the species pair corresponding
c     to unam1, unam2 or unam2, unam1 in the upair array, and returns
c     the index jpair of that pair in that aray. The variable jpair
c     is returned as 0 if no match is found.
c
c     This suboutine is called by:
c
c       EQPT/tprca.f
c       EQPT/tprn2.f
c       EQPT/tprnn.f
c       EQPT/tprnc.f
c       EQPT/tprna.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       npx2t  = number of species pairs in the upair array
c       unam1  = the name of the first aqueous species in a pair
c       unam2  = the name of the second aqueous species in a pair
c       upair  = array containing the names of species in a pair
c
c     Principal output:
c
c       jpair  = the index of the pair unam1, unam2 or unam2, unam1 in
c                  the upair array, if that pair is present; else 0
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer npx2mx
c
      integer jpair,npx2t
c
      character*24 unam1,unam2
      character*24 upair(2,npx2mx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j
c
c-----------------------------------------------------------------------
c
c     Search the upair array for an element matching (unam1, unam2).
c
      do j = 1,npx2t
        if (unam1(1:24) .eq. upair(1,j)(1:24)) then
c
c         Look for unam2 in the matching position in the upair array.
c
          if (unam2(1:24) .eq. upair(2,j)(1:24)) then
            jpair = j
            go to 999
          endif
        endif
      enddo
c
c     Search the upair array for an element matching (unam2, unam1).
c
      do j = 1,npx2t
        if (unam2(1:24) .eq. upair(1,j)(1:24)) then
c
c         Look for unam1 in the matching position in the upair array.
c
          if (unam1(1:24) .eq. upair(2,j)(1:24)) then
            jpair = j
            go to 999
          endif
        endif
      enddo
c
c     Did not find the pair in the upair array.
c
      jpair = 0
c
  999 continue
      end
