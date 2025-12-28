      subroutine gspion(narn1,narn2,nchlor,nelect,nhydr,nhydx,noutpt,
     $ no2gaq,nstmax,nttyo,uspec)
c
c     This subroutine finds the species indices of H+, OH-, Cl-,
c     aqueous O2(g), and aqueous e-.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c      Principal input:
c
c        narn1  = start of the range of aqueous species; also the
c                   species index of solvent water
c        narn2  = end of the range of aqueous species
c        uspec  = array of species names
c
c      Principal output:
c
c        nchlor = species index of the Cl- ion
c        nhydr  = species index of the H+ ion
c        nhydx  = species index of the OH- ion
c        no2gaq = species index of the fictive species, aqueous O2(g)
c        nelect = species index of the fictive species, aqueous e-
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nstmax
c
      integer noutpt,nttyo
c
      integer narn1,narn2,nchlor,nelect,nhydr,nhydx,no2gaq
c
      character(len=48) uspec(nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,nerr,ns
c
      integer ilnobl
c
      character(len=8) uhydr,uhydx,uchlor,uo2gaq,uelect
c
c-----------------------------------------------------------------------
c
      data uhydr  /'H+      '/,uhydx  /'OH-     '/,
     $     uchlor /'Cl-     '/,uo2gaq /'O2(g)   '/,
     $     uelect /'e-      '/
c
c-----------------------------------------------------------------------
c
      nhydr = 0
      nchlor = 0
      nhydx = 0
      no2gaq = 0
      nelect = 0
      nerr = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find H+.
c
      do ns = narn1,narn2
        if (uspec(ns)(1:8) .eq. uhydr(1:8)) then
          nhydr = ns
          go to 100
        endif
      enddo
      j2 = ilnobl(uhydr)
      write (noutpt,1000) uhydr(1:j2)
      write (nttyo,1000) uhydr(1:j2)
 1000 format(/' * Error - (EQLIB/gspion) The species ',a," isn't",
     $ ' present in the model',/7x,'system, as is required. Check to',
     $ ' see that this species is not missing',/7x,'from the',
     $ ' supporting data file. If it is there, and the present input',
     $ /7x,'file is an EQ6 input file, check to see that this file',
     $ ' was created',/7x,'by a previous EQ3NR or EQ6 run using the',
     $ ' same data file.',/7x,'Otherwise,the input file may have been',
     $ ' corrupted.')
      nerr = nerr + 1
c
  100 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find OH-.
c
      do ns = narn1,narn2
        if (uspec(ns)(1:8) .eq. uhydx(1:8)) then
          nhydx = ns
          go to 110
        endif
      enddo
      j2 = ilnobl(uhydx)
      write (noutpt,1000) uhydx(1:j2)
      write (nttyo,1000) uhydx(1:j2)
      nerr = nerr + 1
c
  110 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find Cl-.
c
      do ns = narn1,narn2
        if (uspec(ns)(1:8) .eq. uchlor(1:8)) then
          nchlor = ns
          go to 120
        endif
      enddo
      j2 = ilnobl(uchlor)
      write (noutpt,1000) uchlor(1:j2)
      write (nttyo,1000) uchlor(1:j2)
      nerr = nerr + 1
c
  120 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find the fictive, aqueous O2(g).
c
      do ns = narn1,narn2
        if (uspec(ns)(1:5) .eq. uo2gaq(1:5)) then
          no2gaq = ns
          go to 130
        endif
      enddo
      j2 = ilnobl(uo2gaq)
      write (noutpt,1000) uo2gaq(1:j2)
      write (nttyo,1000) uo2gaq(1:j2)
      write (noutpt,1010)
      write (nttyo,1010)
 1010 format(/7x,'Note- The species O2(g) (the fictive aqueous',
     $ ' species, not',/7x,'the real gas species) is required to be',
     $ ' present even if a problem',/7x,'has no redox aspect. The',
     $ ' species is then inactive.')
      nerr = nerr + 1
c
  130 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find the fictive, aqueous e-.
c
      do ns = narn1,narn2
        if (uspec(ns)(1:2) .eq. uelect(1:2)) then
          nelect = ns
          go to 140
        endif
      enddo
c
c     Note: it is currently not an error for this species to not be
c     present. In future development, it may serve as an alternative
c     to the fictive, aqueous O2(g).
c
  140 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nerr .gt. 0) stop
c
      end
