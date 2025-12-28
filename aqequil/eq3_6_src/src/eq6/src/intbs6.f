      subroutine intbs6(jflag,jflgi,kmax,narn1a,narn2a,nbaspd,nbtd,
     $ nbti,nbtmax,ndrsrd,ndecsp,noutpt,nsta,nstmax,nttyo,uspeca,
     $ ubmtbi)
c
c     This subroutine interprets the data file basis species read from
c     the input file. "Data file" species to be created according to
c     directives read from the input file are ignored. This subroutine
c     also sets up the jflag array.
c
c     This subroutine is called by:
c
c       EQ6/eq6.f
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
      integer kmax,nbtmax,nstmax
c
      integer jflag(nstmax),jflgi(nbtmax),nbaspd(nbtmax),
     $ ndecsp(nbtmax),ndrsrd(2,nstmax)
c
      integer narn1a,narn2a,nbtd,nbti,noutpt,nsta,nttyo
c
      character*56 uspn56
      character*48 uspeca(nstmax),ubmtbi(nbtmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jlen,j2,nb,nbi,nerr,ns,nt
c
      integer ilnobl,nbasis
c
      character*8 ux8
c
c-----------------------------------------------------------------------
c
      nerr = 0
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
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do nbi = 1,nbti
        do ns = 1,nsta
          if (uspeca(ns)(1:48) .eq. ubmtbi(nbi)(1:48)) then
            if (jflgi(nbi).ne.0 .and. jflgi(nbi).ne.30) then
c
c             Calling sequence substitutions:
c               ubmtbi(nbi) for unam48
c
              call fmspnx(jlen,ubmtbi(nbi),uspn56)
              write (ux8,'(i5)') jflgi(nbi)
              call lejust(ux8)
              j2 = ilnobl(ux8)
              write (noutpt,1000) uspn56(1:jlen),ux8(1:j2)
              write (nttyo,1000) uspn56(1:jlen),ux8(1:j2)
 1000         format(/' * Error - (EQ6/intbs6) The species ',a,
     $        /7x,'has an has illegal jflag value of ',a,' on the',
     $        ' input file',/7x,'(the only legal values are 0 and 30).')
              nerr = nerr + 1
            endif
            jflag(ns) = jflgi(nbi)
c
c           Calling sequence substitutions:
c             nbaspd for nbasp
c             nbtd for nbt
c
            nb = nbasis(nbaspd,nbtd,nbtmax,ns)
            if (nb .eq. 0) then
c
c             Add a species to the working basis set.
c
              nbtd = nbtd + 1
              if (nbtd .gt. nbtmax) then
c
c               Calling sequence substitutions:
c                 ubmtbi(nbi) for unam48
c
                call fmspnx(jlen,ubmtbi(nbi),uspn56)
                write (ux8,'(i5)') nbtmax
                call lejust(ux8)
                j2 = ilnobl(ux8)
                write (noutpt,1005) ux8(1:j2),uspn56(1:jlen)
                write (nttyo,1005) ux8(1:j2),uspn56(1:jlen)
 1005           format(/' * Error - (EQ6/intbs6) The maximum ',a,
     $          /7x,'basis species have been exceeded in',
     $          ' interpreting the',/7x,'input file while processing',
     $          ' data file basis species',/7x,a,'. Increase the',
     $          ' dimensioning parameter nbtpar.')
                nerr = nerr + 1
                go to 115
              endif
              nb = nbtd
              nbaspd(nb) = ns
            endif
            ndecsp(nbi) = nb
            go to 115
          endif
        enddo
c
c       Calling sequence substitutions:
c         ubmtbi(nbi) for unam48
c
        call fmspnx(jlen,ubmtbi(nbi),uspn56)
        write (noutpt,1010) uspn56(1:jlen)
        write (nttyo,1010) uspn56(1:jlen)
 1010   format(/' * Note - (EQ6/intbs6) The species "',a,'"',
     $  /7x,'appears as a basis species on the input file, but it',
     $  " wasn't",/7x,'read from the data file. If it is a species',
     $  ' to be created by the code,',/7x,'such as a generic ion',
     $  ' exchanger species, there is no problem.')
c
  115   continue
      enddo
c
      if (nerr .gt. 0) stop
c
      end
