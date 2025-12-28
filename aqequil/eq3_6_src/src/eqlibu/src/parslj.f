      subroutine parslj(nfldmx,nfldt,nlchmx,ufield,uline1,ulscr)
c
c     This subroutine parses the input string uline1, based on the
c     delimiter '|'. Character strings found between delimiters are
c     copied into the output array ufield. The elements of this
c     array are then left-justified.
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
c       EQLIBU/rdd1l.f
c       EQLIBU/rdd1lh.f
c       EQLIBU/rdd2l.f
c       EQLIBU/rdd2lh.f
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
      integer i
c
c-----------------------------------------------------------------------
c
c     Do a simple parse.
c
      call parsln(nfldmx,nfldt,nlchmx,ufield,uline1,ulscr)
c
c     Left-justify the elements of the ufield array.
c
      do i = 1,nfldt
        call lejust(ufield(i))
      enddo
c
      end
