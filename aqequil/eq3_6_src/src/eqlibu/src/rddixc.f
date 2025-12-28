      subroutine rddixc(nxcon,uoldvd,unewf,unewv)
c
c     This subroutine reads user input from the IXCON control file used
c     by XCON3 and XCON6. The IXCON file should look as follows:
c
c     1 (Column 1) (This line is not part of the IXCON file.)
c     IXCON
c
c       This file contains the user-controlled options for
c       the EQ3/6 input file converters XCON3 and XCON6.
c       The code will scan the old input file to determine
c       its format. It will also try to determine its
c       version level. If it can't do so, the version level
c       must be set using the default in this control file.
c
c     OLD INPUT FILE VERSION LEVEL (DEFAULT ONLY)
c      | | 6.0 (versions 6.0-6.1)
c      |X| 7.0 (versions 7.0-7.1)
c      | | 7.2 (version 7.2)
c      | | 8.0 (version 8.0)
c     NEW INPUT FILE FORMAT
c      | | W   (compact)
c      |X| D   (menu-style)
c     NEW INPUT FILE VERSION LEVEL
c      | | 6.0 (versions 6.0-6.1)
c      | | 7.0 (versions 7.0-7.1)
c      |X| 7.2 (version 7.2)
c      | | 8.0 (version 8.0)
c     END
c     1 (Column 1) (This line is not part of the IXCON file.)
c
c
c     Note that the '|' characters defining the checkboxes must be
c     in columns 2 and 4. Select a choice with a non-blank character.
c     The first choice will be selected. Secondary choices will be
c     ignored.
c
c     This subroutine is called by:
c
c       XCON3/xcon3.f
c       XCON6/xcon6.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       nxcon  = the unit number of the IXCON file
c
c     Output:
c
c       uoldvd = the default version level of the old input file
c       unewf  = the desired format ("W" or "D") of the new input file
c       unewvd = the desired version level of the new input file
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nxcon
c
      character*3 uoldvd,unewv
      character*1 unewf
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ii,j
c
      character*80 uline
c
c-----------------------------------------------------------------------
c
      rewind(nxcon)
      do i = 1,1000
        read (nxcon,1000,end=999) uline
 1000   format(a80)
        j = index(uline,'OLD INPUT FILE VERSION LEVEL')
        if (j .gt. 0) then
          do ii = 1,20
            read (nxcon,1000,end=999) uline
            if (uline(2:2).ne.'|' .or. uline(4:4).ne.'|') go to 120
            if (uline(3:3) .ne. ' ') then
              uoldvd = uline(6:8)
              go to 120
            endif
          enddo
          go to 120
        endif
      enddo
  120 continue
c
      rewind(nxcon)
      do i = 1,1000
        read (nxcon,1000,end=999) uline
        j = index(uline,'NEW INPUT FILE FORMAT')
        if (j .gt. 0) then
          do ii = 1,20
            read (nxcon,1000,end=999) uline
            if (uline(2:2).ne.'|' .or. uline(4:4).ne.'|') go to 150
            if (uline(3:3) .ne. ' ') then
              unewf = uline(6:6)
              go to 150
            endif
          enddo
          go to 150
        endif
      enddo
  150 continue
c
      rewind(nxcon)
      do i = 1,1000
        read (nxcon,1000,end=999) uline
        j = index(uline,'NEW INPUT FILE VERSION LEVEL')
        if (j .gt. 0) then
          do ii = 1,20
            read (nxcon,1000,end=999) uline
            if (uline(2:2).ne.'|' .or. uline(4:4).ne.'|') go to 180
            if (uline(3:3) .ne. ' ') then
              unewv = uline(6:8)
              go to 180
            endif
          enddo
          go to 180
        endif
      enddo
  180 continue
c
  999 continue
      end
