      subroutine dfbasp(nbaspa,nbta,nbta_asv,ndrsa,ndrsa_asv,ndrsra,
     $ nerr,noutpt,nsta,nsta_asv,nttyo,ubasp,uspeca)
c
c     This subroutine decodes the ubasp array (list of basis species
c     names built while reading the data file). It puts their
c     species indices in the nbaspa array. The basis indices that
c     were put in the ndrsa array are replaced by the corresponding
c     species indices. The working basis set may be subsequently
c     expanded when the input file is interpreted, as any species
c     on the input file that is referenced as a basis species is
c     added to that set if not already in it.
c
c     This subroutine is called by:
c
c       EQLIB/indata.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c     Principal output:
c
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nbta_asv,ndrsa_asv,nsta_asv
c
      integer nbaspa(nbta_asv),ndrsa(ndrsa_asv),ndrsra(2,nsta_asv)
      integer nbta,nerr,noutpt,nsta,nttyo
c
      character*48 ubasp(nbta_asv),uspeca(nsta_asv)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jlen,n,nb,ns,nr1,nr1p1,nr2,nt
c
      character*56 uspn56
c
c-----------------------------------------------------------------------
c
c     Set up the nbaspa array and expand the elements of the ubasp array
c     to the full 48 characters.
c
      do nb = 1,nbta
        do ns = 1,nsta
          if (uspeca(ns)(1:24) .eq. ubasp(nb)(1:24)) then
            nbaspa(nb) = ns
            ubasp(nb)(1:48) = uspeca(ns)(1:48)
            go to 105
          endif
        enddo
c
        nbaspa(nb) = 0
c
c       Calling sequence substitutions:
c         ubasp(nb) for unam48
c
        call fmspnm(jlen,ubasp(nb),uspn56)
        write (noutpt,1000) uspn56(1:jlen)
        write (nttyo,1000) uspn56(1:jlen)
 1000   format(/" * Error - (EQLIB/dfbasp) Couldn't find a basis",
     $  ' species named',/7x,a,' among the species read from the',
     $  ' data file.',/7x,'It must be referenced erroneously in',
     $  ' the associated reaction',/7x,'for another species.')
        nerr = nerr + 1
  105   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Convert basis species indices in the ndrsa array to species
c     indices.
c
      do ns = 1,nsta
        nt = ndrsra(2,ns) - ndrsra(1,ns) + 1
        if (nt .ge. 2) then
          nr1 = ndrsra(1,ns)
          nr1p1 = nr1 + 1
          nr2 = ndrsra(2,ns)
          do n = nr1p1,nr2
            nb = ndrsa(n)
            ndrsa(n) = nbaspa(nb)
          enddo
        endif
      enddo
c
      end
