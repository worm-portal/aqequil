      subroutine bdslx(narn1,narn2,natmax,noutpt,nslt,nsltmx,
     $ nslx,nsxi,nsxx,nsxmax,nttyo)
c
c     This subroutine builds the nsxi and nsxx arrays. These are
c     pointer arrays used in connection with the S-lambda
c     parts of Pitzer's equations.
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
c       narn1  = start of aqueous species range
c       narn2  = end of aqueous species range
c       natmax = the maximum number of aqeuous species
c       nslx   = array of aqueous species indices defining pairs
c                  of species for which S-lambda data are present
c       nsltmx = the maximum number of pairs of aqueous species
c                  for which Pitzer interaction parameters are
c                  defined; the second dimension of the nslx array
c       nsxmax = the second dimension of the nsxx pointer array
c
c
c     Principal output:
c
c       nsxi   = range pointer array into the nsxx array:
c                  nsxi(1,na) and nsxi(2,na) are the first and last
c                  values of the second subscript (nsx) of the nsxx
c                  array for entries pertaining to the na-th
c                  aqueous species (ns-th species)
c       nsxt   = the total number of entries in the nsxx array
c       nsxx   = pointer array:
c                  nsxx(1,nsx) = the species index of the second
c                  species in the nsl-th pair, and nsxx(2,nsx) = nsl
c
c     Note: the usage of the nsxi and nsxx arrays is illustrated by
c     the following pseudo-code to evaluate the sum:
c     SUM(j) S-lambda(ij)m(j)
c
c       sum = 0.
c       do nsx = nsxi(1,i),nsxi(2,i)
c         j = nsxx(1,nsx)
c         nsl = nsxx(2,nsx)
c         sum = sum + S-lambda(nsl)*m(j)
c       enddo
c
c     This sum is actually evaluated in EQLIBG/gsgsm.f.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer natmax,nsltmx,nsxmax
c
      integer nslx(2,nsltmx),nsxi(2,natmax),nsxx(2,nsxmax)
c
      integer narn1,narn2,noutpt,nslt,nttyo
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n,na,nsx,ns,nsl
c
c-----------------------------------------------------------------------
c
      nsx = 1
c
c     The following loop assumes that water is the narn1-th species.
c
      do ns = narn1 + 1,narn2
        na = ns - narn1 + 1
c
c       Set the beginning of the range for the current species.
c
        nsxi(1,na) = nsx
c
c       Search column 1 of the nslx array for the ns-th species.
c
        do nsl = 1,nslt
          if (nslx(1,nsl) .eq. ns) then
c
c           Found one. Get the other aqueous species index and the
c           nsl index of the pair.
c
            nsxx(1,nsx) = nslx(2,nsl)
            nsxx(2,nsx) = nsl
            nsx = nsx + 1
            if (nsx .gt. nsxmax) then
              n = 2*nsltmx
              write (noutpt,1010) nsxmax,n
              write (nttyo,1010) nsxmax,n
 1010         format(/' * Error - (EQLIBG/bdslx) Have overflow of the',
     $        /7x,'S-lambda pointer array nsxx. Increase the ',
     $        'dimensioning',
     $        /7x,'parameter nsxpar from ',i5,' to no more than ',i5,
     $        '.')
              stop
            endif
          endif
        enddo
c
c       Search column 2.
c
        do nsl = 1,nslt
          if (nslx(2,nsl) .eq. ns) then
c
c           Found one. Skip if the species is also in the first column.
c
            if (ns .ne. nslx(1,nsl)) then
              nsxx(1,nsx) = nslx(1,nsl)
              nsxx(2,nsx) = nsl
              nsx = nsx + 1
              if (nsx .gt. nsxmax) then
                n = 2*nsltmx
                write (noutpt,1010) nsxmax,n
                write (nttyo,1010) nsxmax,n
                stop
              endif
            endif
          endif
        enddo
c
c       Set the end of the range for the current species.
c
        nsxi(2,na) = nsx - 1
      enddo
c
      end
