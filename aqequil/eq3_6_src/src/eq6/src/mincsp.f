      subroutine mincsp(cdrsd,jpflag,jsflag,nbaspd,nbtd,nbtmax,ncmpra,
     $ ndrsd,ndrsmx,ndrsrd,nmrn1a,nmrn2a,noutpt,npta,nptmax,nstmax,
     $ nttyo,nxopex,nxopmx,nxopt,nxpemx,uspeca,uxcat,uxopex,uxopt)
c
c     This subroutine executes the mineral subset-selection suppression
c     (nxopt) options.
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
      integer nbtmax,ndrsmx,nptmax,nstmax,nxopmx,nxpemx
c
      integer jpflag(nptmax),jsflag(nstmax),ncmpra(2,nptmax),
     $ nbaspd(nbtmax),ndrsd(ndrsmx),ndrsrd(2,nstmax)
c
      integer nbtd,nmrn1a,nmrn2a,noutpt,npta,nttyo,nxopex,nxopt
c
      character(len=48) uspeca(nstmax)
      character(len=24) uxcat(nxopmx),uxopex(nxpemx)
      character(len=8) uxopt(nxopmx)
c
      real*8 cdrsd(ndrsmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,j,jlen,j2,j3,n,nb,nerr,np,np1m1,ns,nse
c
      integer ilnobl,nphase
c
      character(len=56) uspn56
      character(len=24) ux24
      character(len=8) ux8
c
      real*8 cx
c
      real*8 coefdr
c
c-----------------------------------------------------------------------
c
      if (nxopt .le. 0) go to 999
c
      nerr = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      i = 0
      do n = 1,nxopt
        ux8 = uxopt(n)
        call locase(ux8)
        if (ux8(1:4) .eq. 'all ') i = i + 1
      enddo
c
      if (i .gt. 1) then
        write (noutpt,1010)
        write (nttyo,1010)
 1010   format(/' * Error - (EQ6/mincsp) The "all" mineral subset-',
     $  'selection suppression',/7x,'option has been specified more',
     $  ' than once.')
        nerr = nerr + 1
      endif
c
      do n = 1,nxopt
        ux8 = uxopt(n)
        call locase(ux8)
        if (ux8(1:4) .eq. 'all ') then
          j2 = ilnobl(uxcat(n))
          if (j2 .gt. 0) then
            write (noutpt,1020) uxcat(n)(1:j2)
            write (nttyo,1020) uxcat(n)(1:j2)
 1020       format(/' * Error - (EQ6/mincsp) A category may not be',
     $      ' specified for an',/7x,'"all" mineral subset-selection',
     $      ' suppression option. Here "',a,'"',/7x,'is specified',
     $      ' as defining a category.')
            nerr = nerr + 1
          endif
        endif
      enddo
c
      do n = 1,nxopt
        ux8 = uxopt(n)
        call locase(ux8)
        if (ux8(1:8).eq.'alwith ' .or. ux8(1:8).eq.'allwith') then
          ux24 = uxcat(n)
          j2 = ilnobl(ux24)
          if (j2 .le. 0) then
            write (noutpt,1030)
            write (nttyo,1030)
 1030       format(/' * Error - (EQ6/mincsp) A category must be',
     $      ' specified for each',/7x,'"allwith" mineral subset-',
     $      'selection suppression option. Have an',/7x,'instance',
     $      ' where this is lacking.')
            nerr = nerr + 1
          else
            do j = n + 1,nxopt
              ux8 = uxopt(j)
              call locase(ux8)
              if (ux8(1:8).eq.'alwith ' .or. ux8(1:8).eq.'allwith') then
                if (uxcat(j)(1:24) .eq. ux24(1:24)) then
                  write (noutpt,1040) ux24(1:j2)
                  write (nttyo,1040) ux24(1:j2)
 1040             format(/' * Error - (EQ6/mincsp) The string "',a,'"',
     $            ' is repeated in defining',/7x,'a category for an',
     $            ' "allwith" mineral subset-selection suppression',
     $            /7x,'option. It may only be so used once.')
                  nerr = nerr + 1
                endif
              endif
            enddo
          endif
        endif
      enddo
c
      if (nerr .gt. 0) stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      ns = nmrn1a
c
c     Calling sequence substitutions:
c       ncmpra for ncmpr
c       npta for npt
c
      np1m1 = nphase(ncmpra,npta,nptmax,ns) - 1
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do n = 1,nxopt
        ux8 = uxopt(n)
        call locase(ux8)
        if (ux8(1:4) .eq. 'all ') then
          write (noutpt,1110)
          write (nttyo,1110)
 1110     format(/' * Note - (EQ6/mincsp) Executing the "all" mineral',
     $    ' subset-selection',/7x,'suppression option.')
c
          if (nxopex .gt. 0) then
            write (noutpt,1120)
            write (nttyo,1120)
 1120       format(/9x,'Minerals that are exceptions are:',/)
          endif
c
          np = np1m1
          do ns = nmrn1a,nmrn2a
            np = np + 1
            if (jsflag(ns) .le. 0) then
              do i = 1,nxopex
                if (uspeca(ns)(1:24) .eq. uxopex(i)(1:24)) then
c
c                 Calling sequence substitutions:
c                   uspeca(ns) for unam48
c
                  call fmspnm(jlen,uspeca(ns),uspn56)
                  write (noutpt,1130) uspn56(1:jlen)
                  write (nttyo,1130) uspn56(1:jlen)
 1130             format(11x,a)
                  go to 200
                endif
              enddo
              jsflag(ns) = 1
              jpflag(np) = 1
            endif
  200       continue
          enddo
c
        elseif (ux8(1:8).eq.'alwith ' .or. ux8(1:8).eq.'allwith') then
          do nb = 1,nbtd
            nse = nbaspd(nb)
            if (uspeca(nse)(1:24) .eq. uxcat(n)(1:24)) then
              j3 = ilnobl(uxcat(n))
              write (noutpt,1140) uxcat(n)(1:j3)
              write (nttyo,1140) uxcat(n)(1:j3)
 1140         format(/' * Note - (EQ6/mincsp) Executing the "alwith ',
     $        a,'" mineral',/7x,'subset-selection suppression option.')
c
              if (nxopex .gt. 0) then
                write (noutpt,1120)
                write (nttyo,1120)
              endif
c
              np = np1m1
              do ns = nmrn1a,nmrn2a
                np = np + 1
                if (jsflag(ns) .le. 0) then
c
c                 Calling sequence substitutions:
c                   cdrsd for cdrs
c                   ndrsd for ndrs
c                   ndrsrd for ndrsr
c
                  cx =  coefdr(cdrsd,ndrsd,ndrsmx,ndrsrd,nse,ns,nstmax)
                  if (cx .ne. 0.) then
                    if (nxopex .gt. 0) then
                      do i = 1,nxopex
                        if (uspeca(ns)(1:24) .eq. uxopex(i)(1:24)) then
c
c                         Calling sequence substitutions:
c                           uspeca(ns) for unam48
c
                          call fmspnm(jlen,uspeca(ns),uspn56)
                          write (noutpt,1130) uspn56(1:jlen)
                          write (nttyo,1130) uspn56(1:jlen)
                          go to 210
                        endif
                      enddo
                    endif
                    jsflag(ns) = 1
                    jpflag(np) = 1
                  endif
                endif
  210           continue
              enddo
              go to 220
            endif
          enddo
c
          j3 = ilnobl(uxcat(n))
          write (noutpt,1150) uxcat(n)(1:j3)
          write (nttyo,1150) uxcat(n)(1:j3)
 1150     format(" * Error - (EQ6/mincsp) Don't recognize",
     $    ' "',a,'" as',/7x,'an argument to the "alwith" mineral',
     $    ' subset-selection',/7x,'suppression option.')
          nerr = nerr + 1
c
  220     continue
        else
c
          j2 = ilnobl(uxopt(n))
          j3 = ilnobl(uxcat(n))
          write (noutpt,1160) uxopt(n)(1:j2),uxcat(n)(1:j3)
          write (nttyo,1160) uxopt(n)(1:j2),uxcat(n)(1:j3)
 1160     format(/" * Error - (EQ6/mincsp) Don't recognize",
     $    ' "',a,1x,a,'"',/7x,'as a mineral subset-selection ',
     $    ' suppression option.')
          nerr = nerr + 1
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1170)
      write (nttyo,1170)
 1170 format(1x)
c
      if (nerr .gt. 0) stop
c
cr * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
