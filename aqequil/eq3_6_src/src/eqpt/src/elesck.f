      subroutine elesck(cessi,nbtmx1,nctmax,ncts,nentei,nerr,noutpt,
     $ ns,nttyo,qblkes,qzeres,uessi,usblkf,uspec)
c
c     This subroutine conducts tests on the elemental composition
c     specified for a species. It detects any blank and duplicate
c     element names in the composition and any zero-valued
c     stoichiometric coefficients.
c
c     Note: any blank names are replaced by the string '<blank>'.
c
c     This subroutine is called by:
c
c       EQPT/pcraq.f
c       EQPT/pcrsg.f
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nbtmx1,nctmax
c
      integer noutpt,nttyo
c
      integer nentei(nctmax)
c
      integer ncts,nerr,ns
c
      logical qblkes,qzeres
c
      character(len=24) uspec(nbtmx1)
      character(len=24) usblkf
      character(len=8) uessi(nctmax)
c
      real(8) cessi(nctmax)
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer j2,j3,j4,n,ncount
c
      integer ilnobl
c
      logical qdupes
c
      character*8 ux8
c
c-----------------------------------------------------------------------
c
c     Check for blank element names.
c
      qblkes = .false.
      ncount = 0
      do n = 1,ncts
        j3 = ilnobl(uessi(n))
        if (j3 .le. 0) then
          ncount = ncount + 1
          uessi(n) = '<blank>'
        endif
      enddo
c
      if (ncount .gt. 0) then
        j2 = ilnobl(uspec(ns))
        j4 = ilnobl(usblkf)
        write (noutpt,1000) uspec(ns)(1:j2),usblkf(1:j4)
        write (nttyo,1000) uspec(ns)(1:j2),usblkf(1:j4)
 1000   format(/' * Error - (EQPT/elesck) The species ',a,
     $  ' appearing',/7x,'on the data file in the ',a,' superblock',
     $  ' has a specified',/7x,'composition with blank chemical',
     $  ' element names',/7x,'in the following positions:',/)
        do n = 1,ncts
          if (uessi(n)(1:7) .eq. '<blank>') then
            write (ux8,'(i5)') n
            call lejust(ux8)
            j4 = ilnobl(ux8)
            write (noutpt,1010) ux8(1:j4)
            write (nttyo,1010) ux8(1:j4)
 1010       format(9x,a)
          endif
        enddo
        qblkes = .true.
        nerr = nerr + 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for a chemical element appearing more than once in the
c     composition.
c
      call nelcck(nctmax,ncts,nentei,nerr,qdupes,uessi)
c
      if (qdupes)  then
        j2 = ilnobl(uspec(ns))
        j4 = ilnobl(usblkf)
        write (noutpt,1040) uspec(ns)(1:j2),usblkf(1:j4)
        write (nttyo,1040) uspec(ns)(1:j2),usblkf(1:j4)
 1040   format(/' * Error - (EQPT/elesck) The species ',a,
     $  ' appearing',/7x,'on the data file in the ',a,' superblock',
     $  ' has a specified',/7x,'elemental composition for which',
     $  ' there is more than one entry',/7x,'for the following',
     $  ' element(s):',/)
        do n = 1,ncts
          if (nentei(n) .gt. 1) then
            write (ux8,'(i5)') nentei(n)
            call lejust(ux8)
            j3 = ilnobl(ux8)
            j4 = ilnobl(uessi(n))
            write (noutpt,1050) uessi(n)(1:j4),ux8(1:j3)
            write (nttyo,1050) uessi(n)(1:j4),ux8(1:j3)
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
      qzeres = .false.
      do n = 1,ncts
        if (cessi(n) .eq. 0.) qzeres = .true.
      enddo
c
      if (qzeres) then
        j2 = ilnobl(uspec(ns))
        j4 = ilnobl(usblkf)
        write (noutpt,1080) uspec(ns)(1:j2),usblkf(1:j4)
        write (nttyo,1080) uspec(ns)(1:j2),usblkf(1:j4)
 1080   format(/' * Error - (EQPT/elesck) The species ',a,
     $  ' appearing',/7x,'on the data file in the ',a,' superblock',
     $  ' has an elemental',/7x,'composition with zero-valued',
     $  ' coefficients for the following elements:',/)
        do n = 1,ncts
          if (cessi(n) .eq. 0.) then
            ux8 = uessi(n)
            j3 = ilnobl(ux8)
            write (noutpt,1090) ux8(1:j3)
            write (nttyo,1090) ux8(1:j3)
 1090       format(9x,a)
          endif
        enddo
        nerr = nerr + 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 return
      end
