      subroutine rdwele(atwt,nch,nco,nct,nctmax,ndata1,ndat0s,ndat1f,
     $ nerr,noutpt,nslist,nttyo,uelem)
c
c     This suboutine reads the chemical elements block from the
c     DATA0 file.
c
c     This suboutine is called by:
c
c       EQPT/eqpt.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       nct    = the number of chemical elements
c       ndat0s = unit number of the stripped DATA0 file
c
c     Principal output:
c
c       atwt   = array of atomic weghts
c       uelem  = array of chemical element names
c       nerr   = cumulative error counter
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nctmax
c
      integer ndata1,ndat0s,ndat1f,nerr,noutpt,nslist,nttyo
c
      integer nch,nco,nct
c
      character(len=8) uelem(nctmax)
c
      real(8) atwt(nctmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,j3,j4,nc
c
      integer ilnobl
c
      character(len=80) uline
      character(len=72) uterm,utermc
      character(len=24) ux24
      character(len=8) ux8a,ux8b
c
c-----------------------------------------------------------------------
c
      uterm(1:48) = '+-----------------------------------------------'
      uterm(49:72) = '------------------------'
      utermc = uterm
      utermc(1:1) = '*'
c
c     Skip to the 'elements' block.
c
  110 read(ndat0s,1000,end=120,err=120) uline
      if (uline(1:8) .ne. 'elements') go to 110
 1000 format(a)
      go to 130
c
  120 write(noutpt,1010)
      write(nttyo,1010)
 1010 format(/' * Error - (EQPT/eqpt) End-of-file hit or other read',
     $ /7x,'error occurred while searching for the elements block.')
      stop
c
  130 continue
c
c     Skip the delimiter line below the "elements" header.
c
      read (ndat0s,1000,end=990,err=995) uline
c
      ux24 = 'elements'
      j2 = ilnobl(ux24)
      write (ndata1) ux24
      write (ndat1f,1020) ux24(1:j2)
 1020 format(a)
      write (ndat1f,1020) utermc(1:72)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the chemical elements data.
c
      write (noutpt,1030)
      write (nslist,1030)
 1030 format(/)
      do nc = 1,nct
        read (ndat0s,1000,end=990,err=995) uline
        read (uline,1040,err=995) uelem(nc),atwt(nc)
 1040   format(a8,f10.5,5x,a8,5x,f10.5)
c
c       Check for a blank element name.
c
        j2 = ilnobl(uelem(nc))
        if (j2 .le. 0) then
          write (ux8a,'(i5)') nc
          call lejust(ux8a)
          j3 = ilnobl(ux8a)
          write (noutpt,1050) ux8a(1:j3)
          write (nttyo,1050) ux8a(1:j3)
 1050     format(/' * Error - (EQPT/rdwele) Have encountered a blank',
     $    ' name for chemical',/7x,'element number ',a,' in the list',
     $    ' of elements on the data file.')
          if (nc .gt. 1) then
            ux8b = uelem(nc - 1)
            if (ux8b(1:7) .ne. '<blank>') then
              j4 = ilnobl(ux8b)
              write (noutpt,1060) ux8b(1:j4)
              write (nttyo,1060) ux8b(1:j4)
 1060         format(7x,'This follows the element ',a,'.')
            endif
          endif
          stop
        endif
c
c       Check for a zero atomic weight.
c
        if (atwt(nc) .le. 0.) then
          j2 = ilnobl(uelem(nc))
          write (noutpt,1070) uelem(nc)(1:j2)
          write (nttyo,1070) uelem(nc)(1:j2)
 1070     format(/' * Error - (EQPT/rdwele) Have encountered a zero',
     $    ' atomic weight for ',a,/7x,'on the data file.')
          nerr = nerr + 1
        endif
c
        write (ndata1) uelem(nc),atwt(nc)
        write (ndat1f,1040) uelem(nc),atwt(nc)
        write (noutpt,1080) uelem(nc),atwt(nc)
        write (nslist,1080) uelem(nc),atwt(nc)
 1080   format(' element = ',a8,', atwt = ',f10.5)
      enddo
c
      write (ndat1f,1020) utermc(1:72)
c
c     Skip the ending delimiter.
c
      read (ndat0s,1000,end=990,err=995) uline
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find the indices of O and H.
c
      nerr = 0
      nco = 0
      do nc = 1,nct
        if (uelem(nc)(1:2) .eq. 'O ') then
          nco = nc
          go to 140
        endif
      enddo
c
      write (noutpt,1140)
      write (nttyo,1140)
 1140 format(/' * Error - (EQPT/rdwele) The element O is not present',
     $ ' on the data file',/7x,'as is required.')
      nerr = nerr + 1
c
  140 continue
c
      nch = 0
      do nc = 1,nct
        if (uelem(nc)(1:2) .eq. 'H ') then
          nch = nc
          go to 150
        endif
      enddo
c
      write (noutpt,1150)
      write (nttyo,1150)
 1150 format(/' * Error - (EQPT/rdwele) The element H is not present',
     $ ' on the data file',/7x,'as is required.')
      nerr = nerr + 1
c
  150 continue
c
      if (nco .gt. 0) then
        if (nco .ne. 1) then
          j2 = ilnobl(uelem(1))
          write (noutpt,1160) uelem(1)(1:j2)
          write (nttyo,1160) uelem(1)(1:j2)
 1160     format(/' * Error - (EQPT/rdwele) The first chemical element',
     $    ' on the data file',/7x,'is ',a,'. The first chemical',
     $    ' element must be O. Note also that',/7x,'the corresponding',
     $    ' strict basis species must be H2O. This is an',
     $    /7x,'idiosyncrasy of EQ3/6.')
          nerr = nerr + 1
        endif
      endif
c
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  990 write (noutpt,2000)
      write (nttyo,2000)
 2000 format(/' * Error - (EQPT/rdwele) Unexpectedly encountered',
     $ /7x,'end-of-file while reading the DATA0 file.')
c
      write (noutpt,2010) nc
      write (nttyo,2010) nc
 2010 format(7x,'The value of the local block line counter is ',i3,'.')
      j2 = ilnobl(uline)
      if (j2 .gt. 0) then
        j2 = min(j2,70)
        write (noutpt,2030) uline(1:j2)
        write (nttyo,2030) uline(1:j2)
 2030   format(7x,'The last line read was the following:',
     $  /7x,'"',a,'"')
      endif
      stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  995 write (noutpt,2040)
      write (nttyo,2040)
 2040 format(/' * Error - (EQPT/rdwele) Encountered a read format',
     $ /7x,'error while reading the DATA0 file.')
      write (noutpt,2010) nc
      write (nttyo,2010) nc
      j2 = ilnobl(uline)
      if (j2 .gt. 0) then
        j2 = min(j2,70)
        write (noutpt,2030) uline(1:j2)
        write (nttyo,2030) uline(1:j2)
      endif
      stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
