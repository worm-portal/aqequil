      subroutine tpraz(nat,natmax,nazt,naztmx,ncvaz,nerr,noutpt,nttyo,
     $ pcvaz,qpdaz,uaqsp,uazp)
c
c     Test and process the hard core diameter (azero) and neutral
c     activity coefficient flag (insgf) data (i.e., the 'bdot' data)
c     read from the DATA0 file. Find and flag errors, such as duplication
c     of data (e.g., two entries for the same aqueous species).
c
c     Check the coverage of entered data against all aqueous solute
c     species present on the data file.
c
c     This subroutine is called by:
c
c       EQPT/eqpt.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       nat    = the number of aqueous species
c       nazt   = the number of specified hard core diameters
c       uaqsp  = array of names of aqueous species
c       uazp   = array of aqueous species names used to specify
c                  hard core diamters on the data file
c
c     Principal output:
c
c       nerr   = error counter
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer natmax,naztmx
c
      integer noutpt,nttyo
c
      integer nat,nazt,ncvaz,nerr
c
      logical qpdaz(natmax)
c
      character*24 uaqsp(natmax),uazp(naztmx)
c
      real*8 pcvaz
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jaz,j2,j3,j4,j5,na,naz,ncount,ndupl,nlistl,nn,nodatc
c
      integer ilnobl
c
      character*24 unam
      character*8 ux8
c
c-----------------------------------------------------------------------
c
c     Limit on the list of aqueous species for which no data
c     were found.
c
      data nlistl / 20 /
c
c-----------------------------------------------------------------------
c
c     Initialize the qpdaz array. A value of .false. for a species
c     means that no azero or insgf data for this species were read
c     from the data file. Solvent water (the first aqueous species)
c     does not have such data, so set the flag to .true. for this
c     species.
c
      qpdaz(1) = .true.
      do na = 2,nat
        qpdaz(na) = .false.
      enddo
c
c     Check the entered azero and insgf data for aqueous species.
c
      nodatc = 0
c
      do na = 2,nat
c
c       Search for uaqsp(na) in the uazp array. That array contains
c       the aqueous species names for which azero and insgf data
c       were read from the data file.
c
        unam = uaqsp(na)
        do naz = 1,nazt
          if (unam(1:24) .eq. uazp(naz)(1:24)) go to 100
        enddo
        naz = 0
  100   continue
c
        if (naz .gt. 0) then
c
c         Have found an entry.
c
          qpdaz(na) = .true.
c
c         Check for duplicate data sets.
c
          ndupl = 0
          do jaz = naz + 1,nazt
            if (unam(1:24) .eq. uazp(jaz)(1:24)) then
              ndupl = ndupl + 1
            endif
          enddo
c
          if (ndupl .gt. 0) then
            j2 = ilnobl(unam)
            if (ndupl .eq. 1) then
              write (noutpt,1010) unam(1:j2)
              write (nttyo,1010) unam(1:j2)
 1010         format(/' * Error - (EQPT/tpraz) Have found a duplicate',
     $        ' entry on the DATA0 file',/7x,'for the hard core',
     $        ' diameter of ',a,'.')
            else
              ux8 = ' '
              write (ux8,'(i5)') ndupl
              call lejust(ux8)
              j5 = ilnobl(ux8)
              write (noutpt,1020) ux8(1:j5),unam(1:j2)
              write (nttyo,1020) ux8(1:j5),unam(1:j2)
 1020         format(/' * Error - (EQPT/tpraz) Have found ',a,
     $        ' duplicate entries on the DATA0 file',
     $        /7x,'for the hard core diameter of ',a,'.')
            endif
            nerr = nerr + ndupl
          endif
c
        else
c
c         No data entry was found on the DATA0 file.
c         Note that qpdaz(na) is left with a value of .false.
c
          nodatc = nodatc + 1
        endif
c
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nodatc .gt. 0) then
        write (noutpt,1040)
        write (nttyo,1040)
 1040   format(/' * Warning - (EQPT/tpraz) Did not find a hard core',
     $  ' diameter entry on the',/7x,'DATA0 file for any of the',
     $  ' following aqueous species:',/)
c
        ncount = 0
        do na = 1,nat
          if (.not.qpdaz(na)) then
            ncount = ncount + 1
            unam = uaqsp(na)
            j4 = ilnobl(unam)
            write (noutpt,1050) unam(1:j4)
            write (nttyo,1050) unam(1:j4)
 1050       format(9x,a)
            if (ncount .eq. nlistl) go to 200
          endif
        enddo
  200   continue
c
        nn = nodatc - ncount
        if (nn .gt. 0) then
          write (ux8,'(i5)') nn
          call lejust(ux8)
          j3 = ilnobl(ux8)
          write (noutpt,1060) ux8(1:j3)
          write (nttyo,1060) ux8(1:j3)
 1060     format(/9x,'plus ',a,' others')
        endif
        write (noutpt,1070)
        write (nttyo,1070)
 1070   format(1x)
      endif
c
      ncvaz = nat - nodatc - 1
      pcvaz = (100.*ncvaz)/float(nat - 1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
