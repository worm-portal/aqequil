      subroutine inttfx(narn1a,narn2a,noutpt,nsta_asv,ntfxa,ntfxmx,
     $ ntfxta,nttyo,tfxa,uspeca,utfxxd)
c
c     This subroutine sets up arrays for handling alkalinity
c     coefficients Species names used to tag alkalinity coefficients
c     are matched against the aqueous species names read from the data
c     file. If a match is not found, the corresponding alkalinity
c     coefficient is ignored.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       uspeca = array of names of species read from the data file
c       utfxxd = array of string entries containing the names of
c                  of species and corresponding alkalinity factors
c
c     Principal output:
c
c       ntfxa  = pointer array containing the indices of species
c                  corresponding to the titration factors in the
c                  tfxa array
c       ntfxta = number of species with alkalinity factors
c       tfxa   = array of titration factors
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
      integer nsta_asv,ntfxmx
c
      integer ntfxa(ntfxmx)
      integer narn1a,narn2a,ntfxta
c
      character*48 uspeca(nsta_asv)
      character*32 utfxxd(ntfxmx)
c
      real*8 tfxa(ntfxmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,j,j2,n,nerr,nn,ns
c
      integer ilnobl
c
      real*8 tfx
c
      character*32 uentry,uscr1,uscr2
      character*24 unam
      character*8 ux8
c
c-----------------------------------------------------------------------
c
      nerr = 0
c
      ntfxta = 0
      do nn = 1,ntfxmx
c
c       Decode the current entry in the pseudo-data file. The
c       content of an entry is a string of the form: '"CaCO3(aq)", 2.0'.
c       Double-quote marks surround the name. A comma separates the
c       name from the number. All blanks are optional.
c
        uentry = utfxxd(nn)
        i = index(uentry,'"')
c
        if (i .le. 0) then
          j2 = ilnobl(uentry)
          write (noutpt,1000) uentry(1:j2)
          write (nttyo,1000) uentry(1:j2)
 1000     format(/' * Error - (EQ3NR/inttfx) Encountered an error',
     $    ' while reading a',/7x,'pseudo-data file of alkalinity',
     $    ' factors. The following entry does not',/7x,'contain an',
     $    ' initial double-quote mark:',/9x,"'",a,"'")
          nerr = nerr + 1
        endif
        if (nerr .gt. 4) stop
c
c       Okay, found the first double-quote mark.
c
        uscr1 = uentry(i + 1:32)
        call lejust(uscr1)
        j = index(uscr1,'"')
c
        if (j .le. 0) then
          j2 = ilnobl(uentry)
          write (noutpt,1010) uentry(1:j2)
          write (nttyo,1010) uentry(1:j2)
 1010     format(/' * Error - (EQ3NR/inttfx) Encountered an error',
     $    ' while reading a',/7x,'pseudo-data file of alkalinity',
     $    ' factors. The following entry does not',/7x,'contain a',
     $    ' second double-quote mark:',/9x,"'",a,"'")
          nerr = nerr + 1
        endif
        if (nerr .gt. 4) stop
c
c       Okay, found the second double-quote mark.
c
        if (j .eq. 1) then
          j2 = ilnobl(uentry)
          write (noutpt,1020) uentry(1:j2)
          write (nttyo,1020) uentry(1:j2)
 1020     format(/' * Error - (EQ3NR/inttfx) Encountered an error',
     $    ' while reading a',/7x,'pseudo-data file of alkalinity',
     $    ' factors. The following entry does not',/7x,'contain a',
     $    ' species name between the double-quote marks:',
     $     /9x,"'",a,"'")
          nerr = nerr + 1
        endif
        if (nerr .gt. 4) stop
c
        unam = uscr1(1:j - 1)
c
c       Test for the end of data.
c
        if (unam(1:6).eq.'endit.') go to 150
c
        uscr2 = uscr1(j:32)
        i = index(uscr2,',')
c
        if (i .le. 0) then
          j2 = ilnobl(uentry)
          write (noutpt,1030) uentry(1:j2)
          write (nttyo,1030) uentry(1:j2)
 1030     format(/' * Error - (EQ3NR/inttfx) Encountered an error',
     $    ' while reading a',/7x,'pseudo-data file of alkalinity',
     $    ' factors. The following entry does not',/7x,'contain a',
     $    ' comma after the species name:',/9x,"'",a,"'")
          nerr = nerr + 1
        endif
        if (nerr .gt. 4) stop
c
        uscr1 = uscr2(i + 1:32)
        call lejust(uscr1)
        ux8 = uscr1(1:8)
c
c       Read the titration factor.
c
        read (ux8,'(f5.1)',err=100,end=100) tfx
        go to 110
c
  100   j2 = ilnobl(uentry)
        write (noutpt,1040) uentry(1:j2)
        write (nttyo,1040) uentry(1:j2)
 1040   format(/' * Error - (EQ3NR/inttfx) Encountered a read error',
     $  ' while reading a',/7x,'pseudo-data file of alkalinity',
     $  ' factors. The following entry does not',/7x,'contain a',
     $  ' readable titration number:',/9x,"'",a,"'")
        nerr = nerr + 1
        if (nerr .gt. 4) stop
c
  110   continue
c
c       Search for a matching species name.
c
c       Calling sequence substitutions:
c         narn1a for nrn1a
c         narn2a for nrn2a
c
        call srchn(narn1a,narn2a,ns,nsta_asv,unam,uspeca)
c
        if (ns .gt. 0) then
c
c         If a previous match was found, skip.
c
          do n = 1,ntfxta
            if (ns .eq. ntfxa(n)) go to 140
          enddo
c
          ntfxta = ntfxta + 1
          ntfxa(ntfxta) = ns
          tfxa(ntfxta) = tfx
  140     continue
        endif
      enddo
  150 if (nerr .gt. 0) stop
c
      end
