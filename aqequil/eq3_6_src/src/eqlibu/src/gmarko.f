      subroutine gmarko(nfldmx,nfldt,nmark,nttyo,ufield,uline1)
c
c     This subroutine determines the choice for an option. This choice
c     is marked with an asterisk on a line containing a field for
c     each choice.
c
c     This subroutine is called by:
c
c       Any
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       nfldmx = the maximum number of fields on line uline1
c       nfldt  = the number of fields on line uline1
c       nttyo  = the unit number of the screen file
c       ufield = the array of field contents in uline1
c       uline1 = the line containing the option being examined
c
c     Output:
c
c       nmark  = index of the first marked field on line uline1
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nfldmx
c
      integer nfldt,nmark,nttyo
c
      character*(*) ufield(nfldmx)
      character*(*) uline1
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j,j2,n,ncount
c
      integer ilnobl
c
c-----------------------------------------------------------------------
c
c     Find the first option marked with an asterisk.
c
      ncount = 0
      nmark = 0
      do n = 2,nfldt
        j = index(ufield(n),'*')
        if (j .gt. 0) then
          ncount = ncount + 1
          if (ncount .eq. 1) nmark = n
        endif
      enddo
c
      j2 = ilnobl(uline1)
      j2 = min(j2,70)
c
      if (ncount .eq. 0) then
        write (nttyo,1030) uline1(1:j2)
 1030   format(/' * Warning - (EQLIBU/gmarko) None of the options'
     $  /7x,'is marked with an asterisk on the line beginning with',
     $  /7x,'"',a,'".')
      endif
c
      if (ncount .gt. 1) then
        write (nttyo,1040) uline1(1:j2)
 1040   format(/' * Warning - (EQLIBU/gmarko) More than one of the'
     $  /7x,'options is marked with an asterisk on the line',
     $  ' beginning with',
     $  /7x,'"',a,'".')
      endif
c
      end
