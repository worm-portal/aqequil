      subroutine ptzchk(narn1,narn2,natmax,nmxi,noutpt,nstmax,nsxi,
     $ nttyo,uspec)
c
c     This subroutine checks each aqueous species in the current model
c     to see if it has any S-lambda or mu Pitzer coefficients. If
c     not, a warning is printed. If a species has at least one
c     coefficient of either type but no S-lambda coefficients (i.e.,
c     has only one or more mu coefficients), a warning is printed.
c     No warning is printed if a species has one or more S-lambda
c     coefficients but no mu coefficients.
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
c       uspec  = names of all species
c       nmxi   = range pointer array into the nmxx array:
c                  nmxi(1,na) and nmxi(2,na) are the first and last
c                  values of the second subscript (nmx) of the nmxx
c                  array for entries pertaining to the na-th
c                  aqueous species (ns-th species)
c       nsxi   = range pointer array into the nsxx array:
c                  nsxi(1,na) and nsxi(2,na) are the first and last
c                  values of the second subscript (nsx) of the nsxx
c                  array for entries pertaining to the na-th
c                  aqueous species (ns-th species)
c
c     Principal output:
c
c        None
c
c     Not used here, but referenced above:
c
c       nmxx   = pointer array:
c                  nmxx(1,nmx) = the species index of the second
c                  species in the nmu-th triplet, nmxx(2,nmx) is the
c                  species index of the third species, and
c                  nmxx(3,nmx) = nmu
c       nsxx   = pointer array:
c                  nsxx(1,nsx) = the species index of the second
c                  species in the nsl-th pair, where
c                  nsxx(2,nsx) = nsl
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer natmax,nstmax
c
      integer nmxi(2,natmax),nsxi(2,natmax)
      integer narn1,narn2,noutpt,nttyo
c
      character*48 uspec(nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,islt,imut,na,ns
      integer ilnobl
c
c-----------------------------------------------------------------------
c
c     The following loop assumes that water is the narn1-th species.
c
      do ns = narn1 + 1,narn2
        na = ns - narn1 + 1
c
        islt = nsxi(2,na) - nsxi(1,na) + 1
        imut = nmxi(2,na) - nmxi(1,na) + 1
c
        if (islt.eq.0 .and. imut.eq.0) then
          if (uspec(ns)(1:6).ne.'O2(g) ' .and.
     $      uspec(ns)(1:3).ne.'e- ') then
            j2 = ilnobl(uspec(ns)(1:24))
            write (noutpt,2000) uspec(ns)(1:j2)
            write (nttyo ,2000) uspec(ns)(1:j2)
 2000       format(/' * Warning - (EQLIBG/ptzchk) The species ',a,
     $      /7x,"has no Pitzer interaction coefficients.")
          endif
        elseif (islt .eq. 0) then
          j2 = ilnobl(uspec(ns)(1:24))
          write (noutpt,2010) uspec(ns)(1:j2)
          write (nttyo ,2010) uspec(ns)(1:j2)
 2010     format(/' * Warning - (EQLIBG/ptzchk) The species ',a,
     $    /7x,"has no Pitzer S-lambda interaction coefficients.")
        endif
c
      enddo
c
      end
