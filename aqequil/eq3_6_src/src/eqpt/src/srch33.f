      subroutine srch33(jtripl,unam1,unam2,unam3,utripl,npx3mx,npx3t)
c
c     This subroutine searches for the species triplet corresponding
c     to unam1, unam2, unam3 in the utripl array, and returns the
c     index jtripl of that triplet in that array. The variable jtripl
c     is returned as 0 if no match is found. The species names in
c     the triplet array must be in the same order to generate a match.
c
c     This subroutine is called by:
c
c       EQPT/wrpz3.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       npx3t  = the number of aqueous species triplets in the utripl
c                  array
c       unam1  = the name of the first aqueous species in a triplet
c       unam2  = the name of the second aqueous species in a triplet
c       unam3  = the name of the third aqueous species in a triplet
c       utripl = array of aqueous species triplets
c
c     Principal output:
c
c       jtripl = the index of the triplet unam1, unam2, unam3 in the
c                  utripl array, if that triplet is present; else 0
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer npx3mx
c
      integer jtripl,npx3t
c
      character(len=24) unam1,unam2,unam3
      character(len=24) utripl(3,npx3mx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j
c
c-----------------------------------------------------------------------
c
      do j = 1,npx3t
        if (unam1(1:24) .eq. utripl(1,j)(1:24)) then
          if (unam2(1:24) .eq. utripl(2,j)(1:24)) then
            if (unam3(1:24) .eq. utripl(3,j)(1:24)) then
c
c             Found in the triplet in the utripl array.
c
              jtripl = j
              go to 999
            endif
          endif
        endif
      enddo
c
c     Did not find the triplet in the utripl array.
c
      jtripl = 0
c
  999 continue
      end
