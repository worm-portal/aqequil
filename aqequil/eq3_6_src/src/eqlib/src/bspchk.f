      subroutine bspchk(jsflag,nbaspd,nbtd,nbtmax,ndrsd,ndrsmx,ndrsrd,
     $ noutpt,nrdxsp,nstmax,nttyo,uspeca)
c
c     This subroutine looks at each active auxiliary basis species.
c     It prints a warning if any other species in the corresponding
c     dissociation reaction is not present.
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
c       jsflag = array of status flags for species
c       nbaspd = array of indices of 'data file' basis species
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
      integer nbtmax,ndrsmx,nstmax
c
      integer noutpt,nttyo
c
      integer jsflag(nstmax),nbaspd(nbtmax),ndrsd(ndrsmx),
     $ ndrsrd(2,nstmax)
      integer nbtd,nrdxsp
c
      character(len=48) uspeca(nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jlen,jlene,n,nb,nr1,nr2,ns,nse,nt1
c
      character(len=56) uspe56,uspn56
c
c-----------------------------------------------------------------------
c
      do nb = 1,nbtd
        ns = nbaspd(nb)
        if (jsflag(ns) .le. 0) then
          nr1 = ndrsrd(1,ns)
          nr2 = ndrsrd(2,ns)
          nt1 = nr2 - nr1 + 1
          if (nt1 .ge. 2) then
            do n = nr1,nr2
              nse = ndrsd(n)
              if (jsflag(nse) .gt. 0) then
                if (nse .ne. nrdxsp) then
c
c                 Calling sequence substitutions:
c                   jlene for jlen
c                   uspeca(nse) for unam48
c                   uspe56 for uspn56
c
                  call fmspnx(jlene,uspeca(nse),uspe56)
c
c                 Calling sequence substitutions:
c                   uspeca(ns) for unam48
c
                  call fmspnx(jlen,uspeca(ns),uspn56)
c
                  write (noutpt,1000) uspn56(1:jlen),uspe56(1:jlene)
                  write (nttyo,1000) uspn56(1:jlen),uspe56(1:jlene)
 1000             format(/' The auxiliary basis species ',a,
     $            ' is active even though',/3x,'it is a dependent',
     $            ' species of ',a,', which is not',/3x,'present',
     $            ' in the current model. This detached auxililary',
     $            ' basis species',/3x,'will behave much like a',
     $            ' strict basis species.')
                endif
              endif
            enddo
          endif
        endif
      enddo
  110 continue
c
      end
