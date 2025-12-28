      subroutine alters(afcnst,apresg,axlks,cdrs,kxmod,narxmx,narxt,
     $ ndrs,ndrsmx,ndrsr,noutpt,npt,nptmax,nst,nstmax,ntpr,ntprmx,
     $ nttyo,nxmdmx,nxmod,tempc,uphase,uspec,uxmod,xlkmod)
c
c     This subroutine alters the log K values of reactions for the
c     destruction of the specified species after they are read from
c     the data file. This is done before the reactions are rewritten
c     by any basis switching. The reactions are identified on the
c     input file by the names (uxmod) of the associated species.
c     The log K values (actually here the corresponding interpolating
c     polynomial coefficients) are altered according the following
c     options:
c
c       kxmod = 0   Replace the log K value by xlkmod
c       kxmod = 1   Augment the log K value by xlkmod
c       kxmod = 2   Destablize the associated species by xlkmod kcal
c                    (at the temperature at the start of the run)
c
c     This subroutine maps any kxmod = 0 or 2 options into the
c     equivalent kxmod = 1 option. An unaltered pickup file (written
c     by EQ3NR or EQ6) therefore may have only the kxmod = 1 alter
c     option in addition to the kxmod = -1 suppression option (which
c     is handled by EQLIB/supprs.f). The effect of this is that if
c     temperature in EQ6 is a function of reaction progress, the
c     log K at any temperature is augmented by a constant amount.
c
c     Note that uxmod is a 48-character name with a 24-character
c     species part followed by a 24-character phase part. If the phase
c     part is blank, which would normally be the case, name matching
c     is done only on the species part. That way, for example,
c     directing the code to change the log K of 'albite' means that
c     the code makes the change for 'albite'  in 'albite,' 'albite' in
c     'plagioclase,' and 'albite' in (Na,K)-feldspar.'
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
c       axlks  = array of polynomial coefficients for computing log K
c                  values (altered by the options carried out by this
c                  subroutine)
c       kxmod  = array of flag variables defining the type of nxmod
c                  options
c       nxmod  = the number of 'nxmod' alter/suppress options
c       uxmod  = array of names of species whose thermodynamic data
c                  are to be altered
c       xlkmod = array of values used by the 'nxmod' options to replace
c                  or augment the the log K values
c
c     Principal output:
c
c       axlks  = array of polynomial coefficients for computing log K
c                  values
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer narxmx,ndrsmx,nptmax,nstmax,ntprmx,nxmdmx
c
      integer noutpt,nttyo
c
      integer kxmod(nxmdmx),narxt(ntprmx),ndrs(ndrsmx),ndrsr(2,nstmax)
      integer npt,nst,ntpr,nxmod
c
      character*48 uxmod(nxmdmx),uspec(nstmax)
      character*24 uphase(nptmax)
c
      real*8 axlks(narxmx,ntprmx,nstmax),apresg(narxmx,ntprmx),
     $ cdrs(ndrsmx),xlkmod(nxmdmx)
      real*8 afcnst,tempc
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j,jlen,j2,n,nchar,nerr,nhits,nn,np,nr1,nr2,ns,nt
c
      integer ilnobl
c
      character*56 uspn56
      character*48 unam48
      character*24 ublk24
      character*8 ufix
c
      real*8 cds,pgrid,xlkdel,xlknew,xlkold
c
c-----------------------------------------------------------------------
c
      data ublk24 /'                        '/
      data ufix   /'fix_    '/
c
c-----------------------------------------------------------------------
c
      nerr = 0
      if (nxmod .le. 0) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the grid pressure.
c
c     Calling sequence substitutions:
c       apresg for arr
c       pgrid for prop
c
      call evdat2(apresg,narxmx,narxt,ntpr,ntprmx,pgrid,tempc)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do n = 1,nxmod
        if (kxmod(n) .lt. 0) go to 130
        unam48 = uxmod(n)
        xlkdel = xlkmod(n)
c
        nchar = 24
        if (unam48(25:48) .ne. ublk24(1:24)) nchar = 48
c
        nhits = 0
        do ns = 1,nst
          if (unam48(1:nchar) .ne. uspec(ns)(1:nchar)) go to 120
c
          nhits = nhits + 1
          nr1 = ndrsr(1,ns)
          nr2 = ndrsr(2,ns)
          nt = nr2 - nr1 + 1
          if (nt .lt. 2) then
            call fmspnm(jlen,unam48,uspn56)
            write (noutpt,1005) uspn56(1:jlen)
            write (nttyo,1005) uspn56(1:jlen)
 1005       format(/' * Error - (EQLIB/alters) The species ',a,
     $      /7x,'is in the strict basis set, so it can not be affected',
     $      /7x,'by the  specifed nxmod alter option.')
            nerr = nerr + 1
            if (nchar .eq. 48) go to 130
            go to 120
          endif
c
c         Calling sequence substitutions:
c           axlks for arr
c           ns for k
c           nstmax for nmax
c           xlkold for prop
c
          call evdat3(axlks,ns,nstmax,narxmx,narxt,ntpr,ntprmx,xlkold,
     $    tempc)
c
c         Map the kxmod = 0 and kxmod = 2 options to the kxmod = 1
c         option.
c
          if (kxmod(n) .eq. 0) then
            kxmod(n) = 1
            xlkdel = xlkdel - xlkold
            xlkmod(n) = xlkdel
          elseif (kxmod(n) .eq. 2) then
            do nn = nr1,nr2
              if (ndrs(nn) .eq. ns) then
                cds = cdrs(nn)
                go to 110
              endif
            enddo
  110       if (cds .le. 0) then
              call fmspnm(jlen,unam48,uspn56)
              write (noutpt,1010) uspn56(1:jlen)
              write (nttyo,1010) uspn56(1:jlen)
 1010         format(/" * Error - (EQLIB/alters) Couldn't find a",
     $        ' non-zero reaction coefficient',/7x,'for the species ',
     $        a,', which is',/7x,'specified in an nxmod alter option.',
     $        " Therefore can't alter the",/7x,'corresponding',
     $        ' equilibrium constant.')
              nerr = nerr + 1
              if (nchar .eq. 48) go to 130
              go to 120
            endif
            xlkdel = -cds*xlkdel/afcnst
            kxmod(n) = 1
            xlkmod(n) = xlkdel
          endif
c
c         Change the polynomial coefficients.
c
          do j = 1,ntprmx
            axlks(1,j,ns) = axlks(1,j,ns) + xlkdel
          enddo
c
c         Calling sequence substitutions:
c           noutpt for nf
c
          call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns,nstmax,uspec)
c
c         Calling sequence substitutions:
c           nttyo for nf
c
          call prreac(cdrs,ndrs,ndrsmx,ndrsr,nttyo,ns,nstmax,uspec)
          xlknew = xlkold + xlkdel
          write (noutpt,1020) tempc,pgrid,xlkold,xlknew
          write (nttyo,1020) tempc,pgrid,xlkold,xlknew
 1020     format(/7x,'The log K of the above reaction at ',f10.3,'C',
     $    /7x,'and ',f10.3,' bars was changed from',
     $    f10.4,' to ',f10.4,'.',/)
          if (nchar .eq. 48) go to 130
c
  120     continue
        enddo
c
        if (nhits .le. 0) then
          call fmspnm(jlen,unam48,uspn56)
          write (noutpt,1025) uspn56(1:jlen)
          write (nttyo,1025) uspn56(1:jlen)
 1025     format(/" * Error- (EQLIB/alters) Can't find the species",
     $    /7x,a,', which is specified in an nxmod alter option.',
     $    ' Check to make sure',/7x,'that a species of this name',
     $    ' appears on the supporting data file.')
c
          j2 = ilnobl(unam48(1:24))
          do np = 1,npt
            if (unam48(1:24) .eq. uphase(np)(1:24)) then
              write (noutpt,1030) unam48(1:j2)
              write (nttyo,1030) unam48(1:j2)
 1030         format(/' * Note - (EQLIB/alters) The entity "',a,'"',
     $        ' specified',/7x,'in an nxmod alter option is a phase,',
     $        ' not a species. Such an option applies only to species.')
              go to 130
            endif
          enddo
        endif
c
c       Check to see that the species is not a fictive fugacity
c       fixing species.
c
        if (unam48(1:4) .eq. ufix(1:4)) then
          call fmspnm(jlen,unam48,uspn56)
          write (noutpt,1035) uspn56(1:jlen)
          write (nttyo,1035) uspn56(1:jlen)
 1035     format(/' * Note - (EQLIB/alters) The species "',a,
     $    '" specified',/7x,'in an nxmod alter option comprises a',
     $    ' fictive fugacity-fixing.',/7x,' phase. Such an option',
     $    " can't be applied to this kind of species.")
        endif
c
  130   continue
      enddo
c
      if (nerr .gt. 0) stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
