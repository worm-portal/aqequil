      subroutine rxnsck(nbtmx1,cdrsi,nct,ndrsts,nentri,nerr,noutpt,
     $ ns,nsb,nttyo,qblkrs,qzerrs,udrsi,usblkf,uspec)
c
c     This subroutine conducts tests on the associated reaction
c     specified for a species. It detects any blank and duplicate
c     species names in the reaction and any zero-valued reaction
c     coefficients. It checks to see if the first species in the
c     reaction is the one formally associated with it, and that
c     The corresponding coefficient is negative.
c
c     Note: any blank names are replaced by the string '<blank>'.
c
c     This subroutine is called by:
c
c       eqpt/pcraq.f
c       eqpt/pcrsg.f
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nbtmx1
c
      integer noutpt,nttyo
c
      integer nentri(nbtmx1)
c
      integer nct,ndrsts,nerr,ns,nsb
c
      logical qblkrs,qzerrs
c
      character(len=24) udrsi(nbtmx1),uspec(nbtmx1)
      character(len=24) usblkf
c
      real(8) cdrsi(nbtmx1)
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer j2,j3,j4,n,ncount
c
      integer ilnobl
c
      logical qduprs
c
      character*24 ux24
      character*8 ux8
c
c-----------------------------------------------------------------------
c
c     Check for blank element names.
c
      qblkrs = .false.
      ncount = 0
      do n = 1,ndrsts
        j3 = ilnobl(udrsi(n))
        if (j3 .le. 0) then
          ncount = ncount + 1
          udrsi(n) = '<blank>'
        endif
      enddo
c
      if (ncount .gt. 0) then
        j2 = ilnobl(uspec(ns))
        j4 = ilnobl(usblkf)
        write (noutpt,1000) uspec(ns)(1:j2),usblkf(1:j4)
        write (nttyo,1000) uspec(ns)(1:j2),usblkf(1:j4)
 1000   format(/' * Error - (eqpt/rxnsck) The species ',a,
     $  ' appearing',/7x,'on the data file in the ',a,' superblock',
     $  ' has a specified',/7x,'reaction with blank species names',
     $  ' in the following positions:',/)
        do n = 1,ndrsts
          if (udrsi(n)(1:7) .eq. '<blank>') then
            write (ux8,'(i5)') n
            call lejust(ux8)
            j4 = ilnobl(ux8)
            write (noutpt,1010) ux8(1:j4)
            write (nttyo,1010) ux8(1:j4)
 1010       format(9x,a)
          endif
        enddo
        qblkrs = .true.
        nerr = nerr + 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for a species appearing more than once in the reaction
c
      call nrxnck(nbtmx1,ndrsts,nentri,nerr,qduprs,udrsi)
c
      if (qduprs) then
c
c       Make an exception for the reaction O2(g) = O2(g) when
c       in the gases superblock and O2(g) is the redox species.
c       Here the O2(g) on the left-hand side is an actual gas
c       species, while that on the right is actually the
c       fictive aqueous species.
c
        if (usblkf(1:5) .eq. 'gases') then
          if (uspec(ns)(1:6) .eq. 'O2(g) ') then
            if (nsb .eq. (nct + 1)) then
              if (uspec(nsb)(1:6) .eq. 'O2(g) ') then
                if (ndrsts .eq. 2) then
                  if (udrsi(1)(1:6) .eq. 'O2(g) ') then
                    if (udrsi(2)(1:6) .eq. 'O2(g) ') then
                      qduprs = .false.
                      nentri(1) = 1
                      nentri(2) = 1
                    endif
                  endif
                endif
              endif
            endif
          endif
        endif
      endif
c
      if (qduprs)  then
        j2 = ilnobl(uspec(ns))
        j4 = ilnobl(usblkf)
        write (noutpt,1040) uspec(ns)(1:j2),usblkf(1:j4)
        write (nttyo,1040) uspec(ns)(1:j2),usblkf(1:j4)
 1040   format(/' * Error - (eqpt/rxnsck) The species ',a,
     $  ' appearing',/7x,'on the data file in the ',a,' superblock',
     $  ' has a specified',/7x,'reaction for which there is more',
     $  ' than one entry for the following',/7x,'species(s):',/)
        do n = 1,ndrsts
          if (nentri(n) .gt. 1) then
            write (ux8,'(i5)') nentri(n)
            call lejust(ux8)
            j3 = ilnobl(ux8)
            j4 = ilnobl(udrsi(n))
            write (noutpt,1050) udrsi(n)(1:j4),ux8(1:j3)
            write (nttyo,1050) udrsi(n)(1:j4),ux8(1:j3)
 1050       format(9x,a,' (',a,' entries)')
          endif
        enddo
        nerr = nerr + 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for zero-valued stoichiometric coefficients.
c
      qzerrs = .false.
      do n = 1,ndrsts
        if (cdrsi(n) .eq. 0.) qzerrs = .true.
      enddo
c
      if (qzerrs) then
        j2 = ilnobl(uspec(ns))
        j4 = ilnobl(usblkf)
        write (noutpt,1080) uspec(ns)(1:j2),usblkf(1:j4)
        write (nttyo,1080) uspec(ns)(1:j2),usblkf(1:j4)
 1080     format(/' * Error - (eqpt/rxnsck) The species ',a,
     $    ' appearing',/7x,'on the data file in the ',a,' superblock',
     $    ' has an associated',/7x,'reaction with zero-valued',
     $    ' coefficients for the following species:',/)
        do n = 1,ndrsts
          if (cdrsi(n) .eq. 0.) then
            ux24 = udrsi(n)
            j3 = ilnobl(ux24)
            write (noutpt,1090) ux24(1:j3)
            write (nttyo,1090) ux24(1:j3)
 1090       format(9x,a)
          endif
        enddo
        nerr = nerr + 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c       Check the name of the species as it appears in the associated
c       reaction. A species must appear first in its associated
c       reaction.
c
        if (ndrsts .ge. 2) then
          if (uspec(ns)(1:24) .ne. udrsi(1)(1:24)) then
            j2 = ilnobl(uspec(ns))
            j4 = ilnobl(usblkf)
            ux24 = udrsi(1)
            j3 = ilnobl(ux24)
            write (noutpt,1210) uspec(ns)(1:j2),usblkf(1:j4),
     $      ux24(1:j3),uspec(ns)(1:j2)
            write (nttyo,1210) uspec(ns)(1:j2),usblkf(1:j4),
     $      ux24(1:j3),uspec(ns)(1:j2)
 1210       format(/' * Error - (eqpt/rxnsck) The species ',a,
     $      ' appearing',/7x,'on the data file in the ',a,' superblock',
     $      ' does not appear',/7x,'in its associated reaction.',
     $      ' Should "',a,'" in the',/7x,'reaction be "',a,'", or',
     $      ' vice versa?')
            nerr = nerr + 1
          else
c
c           Check to see that the associated reaction coefficient is
c           negative.
c
            if (cdrsi(1) .ge. 0.) then
              j2 = ilnobl(uspec(ns))
              j4 = ilnobl(usblkf)
              write (noutpt,1220) uspec(ns)(1:j2),usblkf(1:j4),cdrsi(1)
              write (nttyo,1220) uspec(ns)(1:j2),usblkf(1:j4),cdrsi(1)
 1220         format(/' * Error - (eqpt/rxnsck) The species ',a,
     $        ' appearing',/7x,'on the data file in the ',a,
     $        ' superblock does not have a',/7x,'negative reaction',
     $        ' coefficient, as is required. The coefficient',
     $        /7x,'has a value of ',1pe12.5,'.')
              nerr = nerr + 1
            endif
          endif
        endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 return
      end
