      subroutine dscramc(ntabx,ntabs,nllnmx,ulinex)
c
c     This routine descrambles a file of tables whose lines are
c     interspersed, but which are marked by an id string that
c     apears first on any such line. The id string may contain up
c     to 8 characters and is followed by a comma.
c
c     The present routine is intended to support the creation of
c     .csv file to be opened in Excel or other spreadsheet or plotting
c     software, rather than a simple text file (such as the original
c     EQ6 TAB file).
c
c     The id string conceptually consists of eight characters.
c     Trailing blanks may be omitted on a line read by this
c     subroutine. The id string starts with a six-character tag
c     that identifies the corresponding table. The remaining two
c     characters may be used for other purposes. The string "vh'
c     identifies a variable header. The id string for a "Table A"
c     might be "A       ".
c
c     The contents of the scrambled file are copied to the descrambled
c     file as the descrambling takes place. The descrambled file must
c     already be open.
c
c     This subroutine is called by:
c
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       ntabx  = unit number of the scrambled file
c       ntabs  = unit number of the descrambled file (here the
c                TAB scratch file)
c       nllnmx = the maximum character length of a line on the
c                TABX, TABS, and TAB files
c       ulinex = variable holding a line, including the id string
c
c     Output:
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
      integer ntabx,ntabs
c
      integer nllnmx
c
      character(len=nllnmx) ulinex
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j1,j2,j,jc,n,ntable,ntablx
c
      integer, dimension(:), allocatable :: iftabl,iltabl
c
      logical qtablf,qvhfix
c
      character(len=nllnmx), dimension(:), allocatable :: ulvhtb
c
      character(len=8), dimension(:), allocatable :: utable,utablx
c
      character(len=8) uidstr,uidstx,utx,ux8
c
c-----------------------------------------------------------------------
c
c     Determine the six-character table identifiers (utable)
c     that are present on the scrambled file.
c
      ntablx = 10
      ALLOCATE(utable(ntablx))
      do n = 1,ntablx
        utable(n) = ' '
      enddo
c
      rewind ntabx
      ntable = 0
c
      uidstr = ' '
      read (ntabx,'(a8)',end=110) ux8
      jc = index(ux8,',')
      if (jc .gt. 0) then
        uidstr(1:jc - 1) = ux8(1:jc - 1)
      else
        uidstr(1:8) = ux8(1:8)
      endif
c
      ntable = 1
      utable(1)(1:6) = uidstr(1:6)
c
  100 uidstr = ' '
      read (ntabx,'(a8)',end=110) ux8
      jc = index(ux8,',')
      if (jc .gt. 0) then
        uidstr(1:jc - 1) = ux8(1:jc - 1)
      else
        uidstr(1:8) = ux8(1:8)
      endif
c
c     Is the current line's id string in the list?
c
      do n = 1,ntable
        if (uidstr(1:6) .eq. utable(n)(1:6)) go to 100
      enddo
c
c     Does the list of table identifiers need to be redimensioned?
c
      if (ntable .eq. ntablx) then
        ALLOCATE(utablx(ntablx))
        do n = 1,ntable
          utablx(n) = utable(n)
        enddo
        do n = ntable + 1,ntablx
          utablx(n) = ' '
        enddo
        DEALLOCATE(utable)
        ntablx = ntablx + 10
        ALLOCATE(utable(ntablx))
        do n = 1,ntable
          utable(n) = utablx(n)
        enddo
        do n = ntable + 1,ntablx
          utable(n) = ' '
        enddo
        DEALLOCATE(utablx)
      endif
c
c     Add the new table identifier to the list.
c
      ntable = ntable + 1
      utable(ntable)(1:6) = uidstr(1:6)
      go to 100
c
  110 continue
c
c     Organize the scrambled file so that all the lines belonging
c     to a table are contiguous. Descramble to a scratch file (TABS).
c     In the process, find the last instance of a variable header
c     in each table, if there is one. A variable header instance is
c     marked by the string 'vh' in the last two characters of the
c     id string.
c
      ALLOCATE(ulvhtb(ntable))
      do n = 1,ntable
        ulvhtb(n) = 'None'
      enddo
c
      rewind ntabs
      do n = 1,ntable
c
c       Loop over table designators.
c
        utx = utable(n)
        rewind ntabx
c
  140   uidstr = ' '
        read (ntabx,'(a)',end=150) ulinex
c
        jc = index(ulinex,',')
        if (jc .gt. 0) then
          uidstr(1:jc - 1) = ulinex(1:jc - 1)
        else
          uidstr(1:8) = ulinex(1:8)
        endif
c
        if (uidstr(1:6) .eq. utx(1:6)) then
          j2 = len_trim(ulinex)
          write (ntabs,'(a)') ulinex(1:j2)
c
          if (uidstr(7:8) .eq. 'vh') then
            ux8 = uidstr
            ux8(7:8) = '  '
            j1 = len_trim(ux8)
            jc = index(ulinex,",")
            ulvhtb(n)(1:j1) = ux8(1:j1)
            ulvhtb(n)(j1 + 1:nllnmx) = ulinex(jc:nllnmx)
          endif
c
        endif
        go to 140
  150   continue
      enddo
c
c     Copy the partially descrambled file back over the
c     original descrambled file. After this has been done,
c     the scrambled file will still be suitable for extension
c     in a possible subsequent EQ6 run.
c
      rewind ntabx
      rewind ntabs
c
  170 read (ntabs,'(a)',end=180) ulinex
      j2 = len_trim(ulinex)
      write (ntabx,'(a)') ulinex(1:j2)
      go to 170
c
  180 continue
c
c     The last instance of a possible variable header line
c     has already been found. It is the only instance that
c     is desired, and it will be placed in the relative position
c     of the first instance on the descrambled file.
c
c     A variable header line defines quantities associated with
c     some of the table columns, but whose contents are built up
c     during the run and cannot be known in advance. For example,
c     a table might give the number of moles of precipitated
c     minerals. The identities of these minerals cannot be known
c     until the run is complete. A first instance is written
c     in the desired relative position for the full variable header.
c     This will usually be incomplete. As more minerals are
c     precipitated, updated instances of the variable header will
c     be written. Only the last one is complete (as far as the
c     current run is concerned).
c
c     Loop over the tables.
c
      rewind ntabs
      do n = 1,ntable
c
c       Loop over table designators.
c
        qtablf = .false.
        qvhfix = .false.
        utx = utable(n)
        rewind ntabx
  240   uidstr = ' '
        read (ntabx,'(a)',end=250) ulinex
c
        jc = index(ulinex,',')
        if (jc .gt. 0) then
          uidstr(1:jc - 1) = ulinex(1:jc - 1)
        else
          uidstr(1:8) = ulinex(1:8)
        endif
c
        if (uidstr(1:6) .eq. utx(1:6)) then
          qtablf = .true.
          if (uidstr(7:8) .eq. 'vh') then
            if (.not.qvhfix) then
              qvhfix = .true.
              ulinex = ulvhtb(n)
            else
              go to 240
            endif
          endif
          j2 = len_trim(ulinex)
          write (ntabs,'(a)') ulinex(1:j2)
        elseif (qtablf) then
          go to 250
        endif
        go to 240
  250   continue
      enddo
c
      DEALLOCATE(utable)
      DEALLOCATE(ulvhtb)
c
      end
