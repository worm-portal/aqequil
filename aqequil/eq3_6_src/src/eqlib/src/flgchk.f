      subroutine flgchk(jpflag,jsflag,ncmpra,npta,nptmax,nstmax,qclnsa)
c
c     This subroutine recalculates the jpflag array to insure
c     consistency with the jsflag array. The jpflag value of a phase
c     can't be lower than the lowest value of any of the species which
c     comprise that phase. For example, if all of the species of
c     a phase are not present (jsflag = 2), the phase can not be
c     present (jpflag = 2). Also, if a solution phase (one with more
c     than one component species) has only one species that is active
c     (jsflag .lt. 2) and that species was created by cloning, the
c     solution phase is purged (jpflag is set to 2). This subroutine
c     should be called after calling other subroutines which may purge
c     or suppress species. EQLIB/flgset.f is an example of such.
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
c       jpflag = array of status flags for phases
c       jsflag = array of status flags for species
c       ncmpra = array giving the range of indices of species
c                  belonging to a given phase ('a' data set)
c       npta   = the number of phaes in the 'a' data set.
c       qclnsa = array of logical flags indicating if a species was
c                  created by cloning.
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
      integer nptmax,nstmax
c
      integer jpflag(nptmax),jsflag(nstmax),ncmpra(2,nptmax)
      integer npta
c
      logical qclnsa(nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ijs0,ijs1,jsflmn,njs0,npa,nr1,nr2,nsa,nt
c
c-----------------------------------------------------------------------
c
      do npa = 1,npta
      if (jpflag(npa) .lt. 2) then
        nr1 = ncmpra(1,npa)
        nr2 = ncmpra(2,npa)
        jsflmn = 2
        do nsa = nr1,nr2
          jsflmn = min(jsflag(nsa),jsflmn)
        enddo
        jpflag(npa) = max(jpflag(npa),jsflmn)
      endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do npa = 1,npta
        if (jpflag(npa) .lt. 2) then
          nr1 = ncmpra(1,npa)
          nr2 = ncmpra(2,npa)
          nt = nr2 - nr1 + 1
          if (nt .ge. 2) then
c
c           Note the following local variables:
c
c             njs0 = counter for non-cloned species with jsflag = 0
c             ijs0 = counter for cloned and non-cloned species with
c                    jsflag = 0
c             ijs1 = counter for cloned and non-cloned species with
c                    jsflag = 1
c
            njs0 = 0
            ijs1 = 0
            ijs0 = 0
            do nsa = nr1,nr2
              if (jsflag(nsa) .eq. 0) then
                ijs0 = ijs0 + 1
                if (.not.qclnsa(nsa)) njs0 = njs0 + 1
              elseif (jsflag(nsa) .eq. 1) then
                ijs1 = ijs1 + 1
              endif
            enddo
            if (njs0.le.0 .and. ijs0.lt.2) then
              jpflag(npa) = 1
              if (ijs1 .le. 0) jpflag(npa) = 2
            endif
          endif
        endif
      enddo
c
      end
