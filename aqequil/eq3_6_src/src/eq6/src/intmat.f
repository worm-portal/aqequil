      subroutine intmat(iaqsln,iindx1,ipndx1,kbt,kdim,kelect,khydr,
     $ khydx,kmax,km1,kmt,ko2gaq,kwater,kx1,kxt,narn1,narn2,nbasp,
     $ nbt,nbti,nbtmax,ncmpr,nelect,nern1,nern2,nhydr,nhydx,nobswt,
     $ noutpt,no2gaq,nphasx,npt,nptmax,nstmax,nttyo,qloffg,ubmtbi,
     $ ufixf,uobsw,uphase,uspec,uzveci,uzvec1,zvclgi,zvclg1,zvec1)
c
c     This routine interprets matrix variables read from the input
c     file. It constructs the master variable array zvclgi and
c     builds iindx1, the master variable index array, deleting any
c     left over fixed fugacity phases. Note that iindx1(kcol) gives
c     the basis index for kcol = 1,kbt. For kcol = km1,kxt, it gives
c     the species index. this subroutine also constructs the ipndx1
c     array. Note that ipndx1(kcol) gives the phase index for
c     kcol = 1,kxt.
c
c     This routine also resolves any issues between the matrix variable
c     input (which may imply instances of ordinary basis switching) and
c     the ordinary basis switching input.
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
c       qloffg = logical flag, = .true. if left over fixed fugacity
c                  phases are in the system
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer kmax,nbtmax,nptmax,nstmax
c
      integer noutpt,nttyo
c
      integer iindx1(kmax),ipndx1(kmax),nbasp(nbtmax),ncmpr(2,nptmax),
     $ nphasx(nstmax)
c
      integer iaqsln,kbt,kdim,kelect,khydr,khydx,km1,kmt,ko2gaq,
     $ kwater,kx1,kxt,narn1,narn2,nbt,nbti,nelect,nern1,nern2,nhydr,
     $ nhydx,nobswt,no2gaq,npt
c
      logical qloffg
c
      character(len=48) ubmtbi(nbtmax),uobsw(2,nbtmax),uspec(nstmax),
     $ uzveci(kmax),uzvec1(kmax)
      character(len=24) uphase(nptmax)
      character(len=8) ufixf
c
      real(8) zvclgi(kmax),zvclg1(kmax),zvec1(kmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jlena,jlenb,jlenc,jlend,jlen2,j2,kbts,kcol,km1s,kmts,
     $ krow,kx1s,kxts,n,nb,nbi,nerr,nn,np,nr1,nr2,ns,ns1,ns2,nt
c
      integer ilnobl
c
      character(len=56) uspa56,uspb56,uspc56,uspd56,uspn56
      character(len=24) ux24
      character(len=8) ux8
c
      real(8) lxx
c
      real(8) texp
c
c-----------------------------------------------------------------------
c
      nerr = 0
      qloffg = .false.
      km1s = km1
      kmts = kmt
      kx1s = kx1
      kxts = kxt
      kbts = kbt
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Aqueous species and generic ion exchanger species.
c
      kcol = 0
      do krow = 1,kbts
        do nb = 1,nbt
          ns = nbasp(nb)
c
          if (ns.ge.nern1 .and. ns.le.nern2) ns = ns + 1
c
          if (uzveci(krow)(1:48) .eq. uspec(ns)(1:48)) then
            kcol = kcol + 1
            iindx1(kcol) = nb
            if (ns.ge.narn1 .and. ns.le.narn2) then
c
c             Have an aqueous species.
c
              ipndx1(kcol) = iaqsln
            elseif (ns.ge.nern1 .and. ns.le.nern2) then
c
c             Have a generic ion exchanger species.
c
              ipndx1(kcol) = nphasx(ns)
            endif
            uzvec1(kcol) = uspec(ns)
            lxx = zvclgi(krow)
            zvclg1(kcol) = lxx
            zvec1(kcol) = texp(lxx)
            go to 100
          endif
        enddo
c
c       Did not find the listed active basis species in the
c       "data file" set. It may be a species that is to be
c       switched in by an ordinary basis switch.
c
        if (krow .le. nbti) then
          nbi = krow
          do ns2 = narn1,narn2
            if (uzveci(krow)(1:48) .eq. uspec(ns2)(1:48)) then
c
c             Calling sequence substitutions:
c               jlenb for jlen2
c               uspb56 for uspn56
c
                call fmspnx(jlenb,uzveci(krow),uspb56)
c
              do nb = 1,nbt
                ns1 = nbasp(nb)
                if (ubmtbi(nbi)(1:48) .eq. uspec(ns1)(1:48)) then
c
c                 A basis switch is implied, with the species whose
c                 index is ns2 being switched in for the one whose
c                 index is ns1. The actual switch will be made later.
c                 For the moment, put the name of the species to be
c                 switched out (the "data file" basis species) in the
c                 current position in the uzvec1 array.
c
c                 Calling sequence substitutions:
c                   jlena for jlen2
c                   uspa56 for uspn56
c
                  call fmspnx(jlena,ubmtbi(nbi),uspa56)
c
                  kcol = kcol + 1
                  iindx1(kcol) = nb
                  ipndx1(kcol) = iaqsln
                  uzvec1(kcol) = uspec(ns1)
                  lxx = zvclgi(krow)
                  zvclg1(kcol) = lxx
                  zvec1(kcol) = texp(lxx)
c
c                 Check for explicit basis switch directive.
c
                  do n = 1,nobswt
                    if (uobsw(1,n)(1:48) .eq. uspec(ns1)(1:48)) then
                      if (uobsw(2,n)(1:48) .eq. uspec(ns2)(1:48)) then
c
c                       The n-th ordinary basis switching directive
c                       covers this pair. This is the normally expected
c                       case.
c
                        go to 100
                      else
c
c                       The first species is set to be switched with
c                       a different second species.
c
c                       Calling sequence substitutions:
c                         jlend for jlen2
c                         uspd56 for uspn56
c
                        call fmspnx(jlend,uobsw(2,n),uspd56)
c
                        write (noutpt,1010) uspb56(1:jlenb),
     $                  uspa56(1:jlena),uspd56(1:jlend)
                        write (nttyo,1010) uspb56(1:jlenb),
     $                  uspa56(1:jlena),uspd56(1:jlend)
 1010                   format(/' * Error - (EQ6/intmat) The species ',
     $                  a,' must be in',/7x,'the active basis set for',
     $                  ' the system described on the input',
     $                  /7x,'file. This implies it should replace the',
     $                  ' "data file" basis',/7x,'species ',a,
     $                  '. However, an explicit replacement by',
     $                  ' the species',/7x,a,' is called for by an',
     $                  ' ordinary basis switch specified',
     $                  /7x,'on the input file.')
                        nerr = nerr + 1
                        go to 100
c
                      endif
                    elseif (uobsw(2,n)(1:48) .eq. uspec(ns2)(1:48)) then
c
c                     The second species is set to be switched with a
c                     different first species.
c
c                       Calling sequence substitutions:
c                         jlenc for jlen2
c                         uspc56 for uspn56
c
                        call fmspnx(jlenc,uobsw(1,n),uspc56)
c
                        write (noutpt,1020) uspb56(1:jlenb),
     $                  uspa56(1:jlena),uspc56(1:jlenc)
                        write (nttyo,1020) uspb56(1:jlenb),
     $                  uspa56(1:jlena),uspc56(1:jlenc)
 1020                   format(/' * Error - (EQ6/intmat) The species ',
     $                  a,' must be in',/7x,'the active basis set for',
     $                  ' the system described on the input',
     $                  /7x,'file. This implies it should replace the',
     $                  ' "data file" basis',/7x,'species ',a,
     $                  '. However, an explicit replacing of',
     $                  ' the species',/7x,a,' is called for by an',
     $                  ' ordinary basis switch specified',
     $                  /7x,'on the input file.')
                        nerr = nerr + 1
                        go to 100
c
                    endif
                  enddo
c
c                 The implied ordinary basis switch does not match
c                 any specified switch. Create one.
c
                  write (noutpt,1030) uspa56(1:jlena),
     $            uspb56(1:jlenb)
                  write (nttyo,1030) uspa56(1:jlena),
     $            uspb56(1:jlenb)
 1030             format(/' * Note - (EQ6/intmat) An implied',
     $            ' ordinary basis switch',/7x,'replacing ',a,
     $            ' with ',a,' was found.',/7x,'An explicit switch'
     $            ' will be constructed.')
c
                  nn = nobswt + 1
c
                  if (nn .gt. nbtmax) then
                    write (ux8,'(i8)') nbtmax
                    call lejust(ux8)
                    j2 = ilnobl(ux8)
                    write (noutpt,1040) ux8(1:j2)
                    write (nttyo,1040) ux8(1:j2)
 1040               format(/' * Error - (EQ6/intmat) Cannot construct',
     $              ' an explicit ordinary',/7x,'basis switch because',
     $              ' this would exceed the limit of ',a,' switches',
     $              /7x,'(which is equal to the dimensioned limit',
     $              ' for basis species).')
                    nerr = nerr + 1
                  endif
c
                  nobswt = nn
                  uobsw(1,nn) = uspec(ns1)
                  uobsw(2,nn) = uspec(ns2)
                  go to 100
c
                endif
              enddo
            endif
          enddo
        endif
c
c       Cannot identify the active basis species specified in the
c       uzveci array on the input file.
c
        write (noutpt,1060) uspb56(1:jlenb)
        write (nttyo,1060) uspb56(1:jlenb)
 1060   format(/' * Error - (EQ6/intmat) The species ',a,' is in the',
     $  /7x,'active basis set for the system described on the input',
     $  ' file,',/7x,"but it isn't on the data file and wasn't created",
     $  ' by',/7x,'an input file directive.')
        nerr = nerr + 1
c
  100   continue
      enddo
      kbt = kcol
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Pure phases.
c
      km1 = kcol + 1
      if (kmts .lt. km1s) go to 120
      qloffg = .false.
c
      do krow = km1s,kmts
c
        do np = 1,npt
          nr1 = ncmpr(1,np)
          nr2 = ncmpr(2,np)
          nt = nr2 - nr1 + 1
          if (nt .eq. 1) then
            ns = nr1
            if (uzveci(krow)(1:48) .eq. uspec(ns)(1:48)) then
              kcol = kcol + 1
              iindx1(kcol) = ns
              ipndx1(kcol) = np
              uzvec1(kcol) = uspec(ns)
              lxx = zvclgi(krow)
              zvclg1(kcol) = lxx
              zvec1(kcol) = texp(lxx)
              go to 110
            endif
          endif
        enddo
c
        if (uzveci(krow)(1:5) .eq. ufixf(1:5)) then
          ux24 = uzveci(krow)(6:24)
          j2 = ilnobl(ux24)
          write (noutpt,1070) ux24(1:j2)
          write (nttyo,1070) ux24(1:j2)
 1070     format(/' * Note - (EQ6/intmat) A left over fixed fugacity',
     $    ' phase for ',a,/7x,'is in the system described on the',
     $    ' input file. It will be purged.')
          qloffg = .true.
        else
          call fmspnm(jlen2,uzveci(krow),uspn56)
          write (noutpt,1080) uspn56(1:jlen2)
          write (nttyo,1080) uspn56(1:jlen2)
 1080     format(/' * Error - (EQ6/intmat) The species ',a,
     $    /7x,'is in the system described on the input file, but it',
     $    " isn't",/7x,"on the data file and wasn't created",
     $  ' following an input file directive.')
          nerr = nerr + 1
        endif
c
  110   continue
      enddo
c
  120 kmt = kcol
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Non-aqueous solutions.
c
      kx1 = kcol + 1
c
      do krow = kx1s,kxts
c
        do np = 1,npt
          nr1 = ncmpr(1,np)
          nr2 = ncmpr(2,np)
          nt = nr2 - nr1 + 1
          if (nt.gt.1 .and.
     $     uzveci(krow)(25:48).eq.uphase(np)(1:24)) go to 130
        enddo
c
        call fmspnm(jlen2,uzveci(krow),uspn56)
        write (noutpt,1080) uspn56(1:jlen2)
        write (nttyo,1080) uspn56(1:jlen2)
        nerr = nerr + 1
        go to 140
c
  130   do ns = nr1,nr2
          if (uzveci(krow)(1:24) .eq. uspec(ns)(1:24)) then
            kcol = kcol + 1
            iindx1(kcol) = ns
            ipndx1(kcol) = np
            uzvec1(kcol) = uspec(ns)
            lxx = zvclgi(krow)
            zvclg1(kcol) = lxx
            zvec1(kcol) = texp(lxx)
            go to 140
          endif
        enddo
c
        call fmspnm(jlen2,uzveci(krow),uspn56)
        write (noutpt,1080) uspn56(1:jlen2)
        write (nttyo,1080) uspn56(1:jlen2)
        nerr = nerr + 1
c
  140   continue
      enddo
c
      kxt = kcol
      kdim = kxt
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      kwater = 0
      if (narn1 .ne. 0) then
        do kcol = 1,kbt
          nb = iindx1(kcol)
          ns = nbasp(nb)
          if (ns .eq. narn1) then
            kwater = kcol
            go to 160
          endif
        enddo
c
        call fmspnx(jlen2,uspec(narn1),uspn56)
        write (noutpt,1200) uspn56(1:jlen2)
        write (nttyo,1200) uspn56(1:jlen2)
 1200   format(/' * Error - (EQ6/intmat) The species ',a," isn't",
     $  /7x,'on the input file, as is required.')
        nerr = nerr + 1
      endif
c
  160 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      khydr = 0
      if (nhydr .ne. 0) then
        do kcol = 1,kbt
          nb = iindx1(kcol)
          ns = nbasp(nb)
          if (ns .eq. nhydr) then
            khydr = kcol
            go to 170
          endif
        enddo
c
        if (nhydx .eq. 0) then
          call fmspnx(jlen2,uspec(nhydr),uspn56)
          write (noutpt,1200) uspn56(1:jlen2)
          write (nttyo,1200) uspn56(1:jlen2)
          nerr = nerr + 1
        endif
      endif
c
  170 continue
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      khydx = 0
      if (nhydx .ne. 0) then
        do kcol = 1,kbt
          nb = iindx1(kcol)
          ns = nbasp(nb)
          if (ns .eq. nhydx) then
            khydx = kcol
            go to 180
          endif
        enddo
c
        if (nhydr .eq. 0) then
          call fmspnx(jlen2,uspec(nhydx),uspn56)
          write (noutpt,1200) uspn56(1:jlen2)
          write (nttyo,1200) uspn56(1:jlen2)
          nerr = nerr + 1
        endif
      endif
c
  180 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      ko2gaq = 0
      if (no2gaq .ne. 0) then
        do kcol = 1,kbt
          nb = iindx1(kcol)
          ns = nbasp(nb)
          if (ns .eq. no2gaq) then
            ko2gaq = kcol
            go to 190
          endif
        enddo
c
        if (nelect .eq. 0) then
          call fmspnm(jlen2,uspec(no2gaq),uspn56)
          write (noutpt,1200) uspn56(1:jlen2)
          write (nttyo,1200) uspn56(1:jlen2)
          nerr = nerr + 1
        endif
      endif
c
  190 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      kelect = 0
      if (nelect .ne. 0) then
        do kcol = 1,kbt
          nb = iindx1(kcol)
          ns = nbasp(nb)
          if (ns .eq. nelect) then
            kelect = kcol
            go to 200
          endif
        enddo
c
        if (no2gaq .eq. 0) then
          call fmspnm(jlen2,uspec(nelect),uspn56)
          write (noutpt,1200) uspn56(1:jlen2)
          write (nttyo,1200) uspn56(1:jlen2)
          nerr = nerr + 1
        endif
      endif
c
  200 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nerr .gt. 0) stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
