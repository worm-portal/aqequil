      subroutine parsln(nfldmx,nfldt,nlchmx,ufield,uline1,ulscr)
c
c     This subroutine parses the input string uline1, based on the
c     delimiter '|'. Character strings found between delimiters are
c     copied into the output array ufield. The elements of this array
c     are not justified. EQLIBU/parslj.f can be used to call this
c     subroutine and follow up the parsing with left-justification.
c
c     The scheme is exemplified as follows:
c
c     |  field 1  |  field 2   |  field3  |
c
c     where the first '|' is in column 1 or
c
c        field 1  |  field 2   |  field3  |
c
c     The first field may or may not have a delimiter to its left.
c     Each field must have a delimiter to its right. A field to
c     the right of the last delimiter, if any, is ignored.
c
c     This subroutine is called by:
c
c       EQLIBU/parslj.f
c       XCON3/rd3d7.f
c       XCON6/rd6d7.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       nfldmx = maximum number of fields
c       nlchmx = the dimension of the ufield array
c       uline1 = the line to be parsed
c       ulscr  = scratch character variable
c
c     Output:
c
c       nfldt  = actual number of fields
c       ufield = array of fields
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
c     Note: the character length of the variables uline1, ulstr, and
c     ufield is actually nlchmx. The Sun Fortran compiler does not
c     allow this to be directly specified.
c
      integer nfldmx,nlchmx
c
      integer nfldt
c
      character*(*) uline1,ulscr
      character*(*) ufield(nfldmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,icht,i2,nfld
c
c-----------------------------------------------------------------------
c
      icht = 0
      nfldt = 0
c
      do i = 1,nfldmx
        ufield(i) = ' '
      enddo
c
c     Find the first delimiter. This is normally expected to appear in
c     column 1. The first field then begins in column 2. If the first
c     delimiter appears in a column beyond the first, the first
c     field begins in column 1.
c
      i  = index(uline1,'|')
      if (i .eq. 1) then
        icht = 1
      else
        icht = 0
      endif
c
      do nfld = 1,nfldmx
        ulscr = uline1(icht + 1:nlchmx)
        i = index(ulscr,'|')
        if (i .eq. 0) then
c
c         Have a last field with no delimiter to its right. Ignore it.
c
          go to 999
        else
          i2 = i - 1
        endif
        nfldt = nfldt + 1
        if (nfldt .le. nfldmx) then
          if (i2 .gt. 0) then
            ufield(nfld)(1:i2) = ulscr(1:i2)
          endif
        endif
        icht = icht + i2 + 1
        if (icht .ge. nlchmx) go to 999
      enddo
c
  999 continue
c
      end
