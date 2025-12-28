      subroutine srchne(nrn1a,nrn2a,ns,nsta_asv,unam,uspeca)
c
c     This subroutine matches the species name unam with the
c     corresponding entry in the nrn1a-th through nrn2a-th range of
c     the species name array uspeca. Only the first 24 characters are
c     compared (uspeca has 48 characters, unam only 24). This
c     subroutine returns the species index ns. If there is no match,
c     the species name is converted to all lower case and the code again
c     searches for a match. If none is found, the species name is
c     coverted to all upper case and the code again searches for a
c     match. If no match is found, ns is returned with a value of 0.
c
c     This subroutine is much like EQLIB/srchn.f, except that if no
c     match is found, the above-noted case conversions are tried in
c     an effort to find a match.
c
c     This subroutine is called by:
c
c       None
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       nrn1a  = start of the range of species to search
c       nrn2a  = end of the range of species to search
c       unam   = name of the species whose index is to be found
c       uspeca = array of species names
c
c     Principal output:
c
c       ns     = index of the species whose name is unam
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nsta_asv
c
      integer nrn1a,nrn2a,ns
c
      character*48 uspeca(nsta_asv)
      character*24 unam
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer k,nfa,nfacap,nfdif,nfz,nfzcap,nlettr
c
      character*24 unamlc,unamuc
      character*1 ulettr
c
c-----------------------------------------------------------------------
c
c     Determine the ranges of upper and lower case letters.
c
      nfa = ichar('a')
      nfz = ichar('z')
      nfacap = ichar('A')
      nfzcap = ichar('Z')
      nfdif = ichar('A') - ichar('a')
c
c     Try the species name as is.
c
      call srchn(nrn1a,nrn2a,ns,nsta_asv,unam,uspeca)
c
      if (ns .le. 0) then
c
c       Convert to lower case.
c
        unamlc = unam
        do k = 1,24
          ulettr = unam(k:k)
          nlettr = ichar(ulettr)
          if (nlettr.ge.nfacap .and. nlettr.le.nfzcap)
     $    unamlc(k:k) = char(nlettr - nfdif)
        enddo
c
c       Test to see if this is different. If so, try it.
c
        if (unamlc .ne. unam) then
c
c         Calling sequence substitutions:
c           unamlc for unam
c
          call srchn(nrn1a,nrn2a,ns,nsta_asv,unamlc,uspeca)
        endif
      endif
c
      if (ns .le. 0) then
c
c       Convert to upper case.
c
        unamuc = unam
        do k = 1,24
          ulettr = unam(k:k)
          nlettr = ichar(ulettr)
          if (nlettr.ge.nfa .and. nlettr.le.nfz)
     $    unamuc(k:k) = char(nlettr + nfdif)
        enddo
c
c       Test to see if this is different. If so, try it.
c
        if (unamuc. ne. unam) then
c
c         Calling sequence substitutions:
c           unamuc for unam
c
          call srchn(nrn1a,nrn2a,ns,nsta_asv,unamuc,uspeca)
        endif
      endif
c
      end
