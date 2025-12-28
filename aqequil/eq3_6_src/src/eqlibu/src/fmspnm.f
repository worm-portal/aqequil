      subroutine fmspnm(jlen,unam48,uspn56)
c
c     This subroutine formats a 48-character species name (unam48) into
c     a 56-character string (uspn56) so that the phase part of the
c     name (in the field composed of the second 24 characters) appears
c     in parentheses following the actual species part (in the field
c     composed of the first 24 characters). For example, ignoring
c     trailing blanks, "Albite                  Plagioclase     " is
c     reformatted as "Albite (Plagioclase)". However, the phase
c     part (including the surrounding parentheses) is omitted from
c     the formatted string if the phase name is not present in the
c     48-character string or if the phase name is identical to the
c     species name. Thus, the returned string for pure Albite is
c     'Albite', not 'Albite (Albite)'.
c
c     The formatted string is used mostly in error, warning, and note
c     messages written by EQ3NR and EQ6.
c
c     This subroutine is similar to fmspnx.f. However, that subroutine
c     omits the phase part if the phase name is 'Aqueous solution'.
c
c     This subroutine is called by:
c
c       Any
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       unam48 = the species name (unformatted)
c
c     Principal output:
c
c       jlen   = the length of the formatted species name
c       uspn56 = the formatted species name
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer jlen
c
      character*48 unam48
      character*56 uspn56
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2
c
      integer ilnobl
c
      character*24 uphnam
c
c-----------------------------------------------------------------------
c
      uspn56 = unam48(1:24)
      jlen = ilnobl(uspn56)
c
      uphnam = unam48(25:48)
      j2 = ilnobl(uphnam)
      if (uphnam(1:24) .ne. unam48(1:24)) then
        if (j2 .gt. 0) then
          uspn56(jlen + 1:jlen + 2) = ' ('
          uspn56(jlen + 3:jlen + 2 + j2) = uphnam(1:j2)
          jlen = jlen + 3 + j2
          uspn56(jlen:jlen) = ')'
        endif
      endif
c
      end
