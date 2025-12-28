      subroutine bdmlx(narn1,narn2,natmax,nmut,nmutmx,nmux,nmxi,
     $ nmxmax,nmxx,noutpt,nttyo)
c
c     This subroutine builds the nmxi and nmxx arrays. These are
c     pointer arrays used in connection with the mu parts of
c     Pitzer's equations.
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
c       narn1  = start of range of aqueous species
c       narn2  = end of range of aqueous species
c       natmax = the maximum number of aqeuous species
c       nmux   = array of aqueous species indices defining triplets
c                  of species for which mu data are present
c       nmutmx = the maximum number of triplets of aqueous species
c                  for which Pitzer interaction parameters are
c                  defined; the second dimension of the nmux array
c       nmxmax = the second dimension of the nmxx pointer array
c
c     Principal output:
c
c       nmxi   = range pointer array into the nmxx array:
c                  nmxi(1,na) and nmxi(2,na) are the first and last
c                  values of the second subscript (nmx) of the nmxx
c                  array for entries pertaining to the na-th
c                  aqueous species (ns-th species)
c       nmxx   = pointer array:
c                  nmxx(1,nmx) = the species index of the second
c                  species in the nmu-th triplet, nmxx(2,nmx) is the
c                  species index of the third species, and
c                  nmxx(3,nmx) = nmu
c
c     Note: the usage of the nmxi and nmxx arrays is illustrated by
c     the following pseudo-code to evaluate SUM(jk) mu(ijk)m(j)m(k),
c     where i, j, and k are species indices and na is the aqueous
c     species index of the i-th species:
c
c       sum = 0.
c       i = na + narn1 - 1
c       do nmx = nmxi(1,na),nmxi(2,na)
c         j = nmxx(1,nmx)
c         k = nmxx(2,nmx)
c         nmu = nmxx(3,nmx)
c         sum = sum + mu(nmu)*m(j)*m(k)
c       enddo
c       sum = 2.*sum
c
c     This sum is actually evaluated in EQLIBG/gmdsm.f.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer natmax,nmutmx,nmxmax
c
      integer nmux(3,nmutmx),nmxi(2,natmax),nmxx(3,nmxmax)
c
      integer narn1,narn2,nmut,noutpt,nttyo
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n,na,nmu,nmx,ns
c
c-----------------------------------------------------------------------
c
      nmx = 1
c
c     The following loop assumes that water is the narn1-th species.
c
      do ns = narn1 + 1,narn2
c
        na = ns - narn1 + 1
c
c       Set the beginning of the range for the current species.
c
        nmxi(1,na) = nmx
c
c       Search column 1 of the nmux array for the na-th aqueous species.
c
        do nmu = 1,nmut
          if (nmux(1,nmu) .eq. ns) then
c
c           Found one. Get the other two aqueous species indices and the
c           nmx index of the triplet.
c
            nmxx(1,nmx) = nmux(2,nmu)
            nmxx(2,nmx) = nmux(3,nmu)
            nmxx(3,nmx) = nmu
            nmx = nmx + 1
            if (nmx .gt. nmxmax) then
              n = 3*nmxmax
              write (noutpt,1010) nmxmax,n
              write (nttyo,1010) nmxmax,n
 1010         format(/' * Error - (EQLIBG/bdmlx) Have overflow of the',
     $        /7x,'mu pointer array nmxx. Increase the dimensioning',
     $        /7x,'parameter nmxpar from ',i5,' to no more than ',
     $        i5,'.')
              stop
            endif
          endif
        enddo
c
c       Search column 2.
c
        do nmu = 1,nmut
          if (nmux(2,nmu) .eq. ns) then
c
c           Found one. Skip if the same species is also in the first
c           column.
c
            if (ns .ne. nmux(1,nmu)) then
              nmxx(1,nmx) = nmux(1,nmu)
              nmxx(2,nmx) = nmux(3,nmu)
              nmxx(3,nmx) = nmu
              nmx = nmx + 1
              if (nmx .gt. nmxmax) then
                n = 3*nmxmax
                write (noutpt,1010) nmxmax,n
                write (nttyo,1010) nmxmax,n
                stop
              endif
            endif
          endif
        enddo
c
c       Search column 3.
c
        do nmu = 1,nmut
          if (nmux(3,nmu) .eq. ns) then
c
c           Found one. Skip if the species is also in column 1 or
c           column 2.
c
            if (ns.ne.nmux(1,nmu) .and. ns.ne.nmux(2,nmu)) then
              nmxx(1,nmx) = nmux(1,nmu)
              nmxx(2,nmx) = nmux(2,nmu)
              nmxx(3,nmx) = nmu
              nmx = nmx + 1
              if (nmx .gt. nmxmax) then
                n = 3*nmxmax
                write (noutpt,1010) nmxmax,n
                write (nttyo,1010) nmxmax,n
              stop
            endif
          endif
        endif
        enddo
c
c       Set the end of the range for the current species.
c
        nmxi(2,na) = nmx - 1
      enddo
c
      end
