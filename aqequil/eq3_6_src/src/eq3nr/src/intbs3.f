      subroutine intbs3(covali,ier,jflag,jflgi,narn1a,narn2a,nbaspd,
     $ nbtd,nbti,nbtmax,ndrsrd,ndecsp,noutpt,nrdxsp,nsta,nstmax,
     $ nttyo,uspeca,uspeci)
c
c     This subroutine interprets the basis species listed on the input
c     file. It sets up the jflag arrays.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
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
      integer nbtmax,nstmax
c
      integer jflag(nstmax),jflgi(nbtmax),nbaspd(nbtmax),
     $ ndecsp(nbtmax),ndrsrd(2,nstmax)
c
      integer ier,narn1a,narn2a,nbtd,nbti,noutpt,nrdxsp,nsta,nt,nttyo
c
      character(len=48) uspeca(nstmax),uspeci(nbtmax)
c
      real(8) covali(nbtmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,nb,nbi,ns
c
      integer ilnobl,nbasis
c
c-----------------------------------------------------------------------
c
      ier = 0
c
      do ns = 1,nsta
        jflag(ns) = 0
      enddo
c
      do ns = narn1a,narn2a
        jflag(ns) = 30
      enddo
c
      do nb = 1,nbtd
        ns = nbaspd(nb)
        nt = ndrsrd(2,ns) - ndrsrd(1,ns) + 1
        if (nt .lt. 2) jflag(ns) = -1
      enddo
c
      jflag(narn1a) = 0
      if (nrdxsp .gt. 0) jflag(nrdxsp) = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do nbi = 1,nbti
        j2 = ilnobl(uspeci(nbi)(1:24))
        do ns = 1,nsta
          if (uspeca(ns)(1:24) .eq. uspeci(nbi)(1:24)) then
            if (ndecsp(nbi) .gt. 0) then
              write (noutpt,1000) uspeci(nbi)(1:j2)
              write (nttyo,1000) uspeci(nbi)(1:j2)
 1000         format(/' * Error - (EQ3NR/intbs3) More than one',
     $        ' constraint was specified',/7x,'on the input file for',
     $        /7x,'the species ',a,'.')
              ier = ier + 1
            else
              if (jflgi(nbi).le.15 .and. covali(nbi).le.0.)
     $        jflgi(nbi) = -1
              jflag(ns) = jflgi(nbi)
c
c             Calling sequence substitutions:
c               nbaspd for nbasp
c               nbtd for nbt
c
              nb = nbasis(nbaspd,nbtd,nbtmax,ns)
              if (nb .eq. 0) then
c
c               Add a species to the working basis set.
c
                nbtd = nbtd + 1
                if (nbtd .gt. nbtmax) then
                  write (noutpt,1005) nbtmax,uspeci(nbi)(1:j2)
                  write (nttyo,1005) nbtmax,uspeci(nbi)(1:j2)
 1005             format(/' * Error - (EQ3NR/intbs3) The maximum ',i7,
     $            ' basis species have been exceeded',/7x,'while',
     $            ' trying to load ',a,'. Increase the dimensioning',
     $            /7x,'parameter nbtpar.')
                  ier = ier + 1
                  go to 115
                endif
                nb = nbtd
                nbaspd(nb) = ns
              endif
              ndecsp(nbi) = nb
            endif
            go to 115
          endif
        enddo
c
        write (noutpt,1010) uspeci(nbi)(1:j2)
        write (nttyo,1010) uspeci(nbi)(1:j2)
 1010   format(/' * Error - (EQ3NR/intbs3) The species ',a,' is',
     $  /7x,"on the input file, but it isn't on the supporting",
     $  ' data file.')
        ier = ier + 1
c
  115   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
