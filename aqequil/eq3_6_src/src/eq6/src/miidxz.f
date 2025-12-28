      subroutine miidxz(ier,iindx1,ipndx1,jpflag,jsflag,kbt,kdim,
     $ kmax,km1,kmt,kx1,kxt,losp,ncmpr,noutpt,npt,nptmax,nstmax,
     $ nttyo,uspec,uzvec1,zvclg1,zvec1)
c
c     This subroutine modifies the matrix indexing originally read in
c     from the input file. Modification occurs whenever a phase
c     boundary is crossed.
c
c     This subroutine is called by:
c
c       EQ6/dumpdp.f
c       EQ6/eqcalc.f
c       EQ6/setffg.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c     Principal output:
c
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
      integer kmax,nptmax,nstmax
c
      integer iindx1(kmax),ipndx1(kmax),jpflag(nptmax),jsflag(nstmax),
     $ ncmpr(2,nptmax)
c
      integer ier,kbt,kdim,km1,kmt,kx1,kxt,noutpt,npt,nttyo
c
      character*48 uspec(nstmax),uzvec1(kmax)
c
      real*8 losp(nstmax),zvclg1(kmax),zvec1(kmax)
c
c-----------------------------------------------------------------------
c
c     Local sequence variable declarations.
c
      integer jlen,kcol,nr1,nr2,np,ns,nt
c
      character*56 uspn56
c
      real*8 lxx
c
      real*8 texp
c
c-----------------------------------------------------------------------
c
      ier = 0
c
      kcol = kbt
c
c     Do non-aqueous phases with only one species.
c
      do np = 1,npt
        if (jpflag(np) .eq. -1) then
          nr1 = ncmpr(1,np)
          nr2 = ncmpr(2,np)
          nt = nr2 - nr1 + 1
          if (nt .eq. 1) then
            ns = nr1
            kcol = kcol + 1
            if (kcol .le. kmax) then
              iindx1(kcol) = ns
              ipndx1(kcol) = np
              uzvec1(kcol) = uspec(ns)
              jsflag(ns) = -1
              lxx = losp(ns)
              zvclg1(kcol) = lxx
              zvec1(kcol) = texp(lxx)
            else
c
c             Calling sequence substitutions:
c               uspec(ns) for unam48
c
              call fmspnm(jlen,uspec(ns),uspn56)
c
              write (noutpt,1000) kmax,uspn56(1:jlen)
              write (nttyo,1000) kmax,uspn56(1:jlen)
 1000         format(/' * Warning - (EQ6/miidxz) Have exceeded the',
     $        ' maximum ',i4,' elements of the',/7x,'iindx1 array',
     $        ' while trying to add ',a,/7x,'as a member of the',
     $        ' equilibrium system. Increase the dimensioning',
     $        /7x,'parameter kpar.')
              ier = 1
            endif
          endif
        endif
      enddo
      km1 = kbt + 1
      kmt = kcol
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Do non-aqueous phases with more than one species.
c
      do np = 2,npt
        if (jpflag(np) .eq. -1) then
          nr1 = ncmpr(1,np)
          nr2 = ncmpr(2,np)
          nt = nr2 - nr1 + 1
          if (nt .ge. 2) then
            do ns = nr1,nr2
              if (jsflag(ns) .eq. -1) then
                kcol = kcol + 1
                if (kcol .le. kmax) then
                  iindx1(kcol) = ns
                  ipndx1(kcol) = np
                  uzvec1(kcol) = uspec(ns)
                  jsflag(ns) = -1
                  lxx = losp(ns)
                  zvclg1(kcol) = lxx
                  zvec1(kcol) = texp(lxx)
                else
c
c                 Calling sequence substitutions:
c                   uspec(ns) for unam48
c
                  call fmspnm(jlen,uspec(ns),uspn56)
c
                  write (noutpt,1000) uspn56(1:jlen)
                  write (nttyo,1000) uspn56(1:jlen)
                  stop
                endif
              endif
            enddo
          endif
        endif
      enddo
c
      kx1 = kmt + 1
      kxt = kcol
      kdim = kcol
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
