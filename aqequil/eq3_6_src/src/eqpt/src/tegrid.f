      subroutine tegrid(itgenf,nacdpr,narxt,nerr,noutpt,ntprmx,
     $ ntprt,nttyo,nwarn,ustrgr)
c
c     This subroutine tests the contents of a "log K" temperature
c     grid for sparse contents in any range. The response to an
c     instance of sparseness is determined by the variable
c     itgenf (-1 = ignore, 0 = warn, 1 = error). This subroutine
c     need not be called if itgenf = -1.
c
c     This suboutine is called by:
c
c       EQPT/rdpar.f
c       EQPT/pcraq.f
c       EQPT/pcrsg.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       itgenf = integer flag controlling the response to a detected
c                  incidence of sparseness (-1 = ignore, 0 = warn,
c                  1 = error)
c       nacdpr = array giving the number of points containing actual
c                  data in a temperature range for the "log K" grid
c                  currently being examined
c       narxt  = array giving the number of points (actual data plus
c                  no data) in that range
c       ntprt  = the number of ranges in the current "log K" grid
c       ustrgr = string describing what the grid represents
c
c     Principal output:
c
c       nerr   = cumulative error counter
c       nwarn  = cumulative warning counter
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ntprmx
c
      integer noutpt,nttyo
c
      integer itgenf,nerr,ntprt,nwarn
c
      integer nacdpr(ntprmx),narxt(ntprmx)
c
      character(len=56) ustrgr
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,j3,j4,ntpr
c
      logical qx
c
      integer ilnobl
c
      character(len=8) ux8a,ux8b
c
c-----------------------------------------------------------------------
c
      do ntpr = 1,ntprt
c
        if (narxt(ntpr) .le. 3) then
c
c         For a range with up to three grid points, at least two
c         must contain actual data. Note that a range can not
c         be composed of fewer than two grid points.
c
          qx = nacdpr(ntpr) .lt. 2
        else
c
c         Otherwise, at least three points containing actual data
c         are required.
c
          qx = nacdpr(ntpr) .lt. 3
        endif
c
        if (qx .and. itgenf.ge.0) then
          write (ux8a,'(i5)') nacdpr(ntpr)
          call lejust(ux8a)
          j2 = ilnobl(ux8a)
          write (ux8b,'(i5)') ntpr
          call lejust(ux8b)
          j3 = ilnobl(ux8b)
          call lejust(ustrgr)
          j4 = ilnobl(ustrgr)
          if (itgenf .eq. 0) then
            nwarn = nwarn + 1
            if (nacdpr(ntpr) .gt. 0) then
              write (noutpt,1000) ux8a(1:j2),ux8b(1:j3),ustrgr(1:j4)
              write (nttyo,1000) ux8a(1:j2),ux8b(1:j3),ustrgr(1:j4)
 1000         format(/' * Warning - (EQPT/tegrid) The number of grid',
     $        ' points containing actual',/7x,'data is only ',a,' in',
     $        ' range ',a,' of the temperature grid representing',
     $        /7x,a,'.')
            else
              write (noutpt,1010) ux8b(1:j3),ustrgr(1:j4)
              write (nttyo,1010) ux8b(1:j3),ustrgr(1:j4)
 1010         format(/' * Warning - (EQPT/tegrid) There are no grid',
     $        ' points containing actual',/7x,'data in range ',a,
     $        ' of the temperature grid representing',/7x,a,'.')
            endif
          else
            nerr = nerr + 1
            if (nacdpr(ntpr) .gt. 0) then
              write (noutpt,1030) ux8a(1:j2),ux8b(1:j3),ustrgr(1:j4)
              write (nttyo,1030) ux8a(1:j2),ux8b(1:j3),ustrgr(1:j4)
 1030         format(/' * Error - (EQPT/tegrid) The number of grid',
     $        ' points containing actual',/7x,'data is only ',a,' in',
     $        ' range ',a,' of the temperature grid representing',
     $        /7x,a,'.')
            else
              write (noutpt,1040) ux8b(1:j3),ustrgr(1:j4)
              write (nttyo,1040) ux8b(1:j3),ustrgr(1:j4)
 1040         format(/' * Error - (EQPT/tegrid) There are no grid',
     $        ' points containing actual',/7x,'data in range ',a,
     $        ' of the temperature grid representing',/7x,a,'.')
            endif
          endif
        endif
      enddo
c
  999 continue
      end
