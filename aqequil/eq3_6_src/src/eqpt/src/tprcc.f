      subroutine tprcc(alamcc,atheta,iccpr,ipbtmx,jpfcmx,natmax,
     $ nccpr,ncvcc,nerr,noutpt,npx3mx,nthdt,nttyo,nwarn,pcvcc,
     $ qpdcc,uaqsp,uthdtr)
c
c     Test and process the Pitzer data for cc' (cation, different
c     cation) pairs read from the DATA0 file. Find and flag errors,
c     such as duplication of data (e.g., two data blocks for the same
c     cc' pair). Calculate the conventional primitive Pitzer parameters
c     from the observable compound parameters read from the data file:
c
c       theta(cc') -> lambda(cc')
c
c     Check the coverage of entered Pitzer data against all possible
c     cation-distinct cation pairs that can be composed of the aqueous
c     cations present on the data file.
c
c     This subroutine is complementary to tpraa.f, which does the
c     same thing for aa' (anion, different anion) pairs.
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
c       uaqsp  = array of names of aqueous species
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
      integer ipbtmx,jpfcmx,natmax,nccpr,npx3mx
c
      integer noutpt,nttyo
c
      integer iccpr(2,nccpr)
c
      integer ncvcc,nerr,nthdt,nwarn
c
      logical qpdcc(nccpr)
c
      character(len=24) uaqsp(natmax),uthdtr(3,npx3mx)
c
      real*8 alamcc(jpfcmx,0:ipbtmx,nccpr),atheta(jpfcmx,npx3mx)
c
      real*8 pcvcc
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,j,jth,jthpr,j2,j3,j4,j5,n,ncount,ndupl,nlistl,
     $ nn,nodatc
c
      integer ilnobl
c
      logical qmatch,qzerov
c
      character*56 ustr56
      character*24 unam1,unam2
      character*8 ux8
c
c-----------------------------------------------------------------------
c
c     Limit on the list of pairs for which no data were found.
c
      data nlistl / 20 /
c
c-----------------------------------------------------------------------
c
c     Initialize the data arrays.
c
      do n = 1,nccpr
        qpdcc(n) = .false.
        do j = 1,jpfcmx
          do i = 0,ipbtmx
            alamcc(j,i,n) = 0.
          enddo
        enddo
      enddo
c
c     Check the entered data for cc' pairs.
c
      nodatc = 0
c
      do n = 1,nccpr
        i = iccpr(1,n)
        j = iccpr(2,n)
        unam1 = uaqsp(i)
        unam2 = uaqsp(j)
c
        jthpr = 0
        ndupl = 0
c
c       Search for unam1, unam2 in the species triplet blocks.
c       This is actually done by searching the uthdtr array,
c       not the utripl array.
c
        do jth = 1,nthdt
          if (unam1(1:24) .eq. uthdtr(1,jth)(1:24)) then
            if (unam2(1:24) .eq. uthdtr(2,jth)(1:24)) then
              jthpr = jth
              go to 110
            endif
          endif
        enddo
  110   continue
c
        if (jthpr .gt. 0) then
c
c         Have found an entry in the triplet blocks.
c
          qpdcc(n) = .true.
c
c         Store the data.
c
          do j = 1,jpfcmx
            alamcc(j,0,n) = atheta(j,jthpr)
            do i = 1,ipbtmx
              alamcc(j,i,n) = 0.
            enddo
          enddo
c
c         Search for duplicates in the species triplets blocks.
c         Ignore duplications in the triplet blocks if the
c         values are all zeros (take this to mean "no data input")
c         or if the values all match the first-encountered
c         values (exact duplication).
c
          do jth = jthpr + 1,nthdt
            if (unam1(1:24) .eq. uthdtr(1,jth)(1:24)) then
              if (unam2(1:24) .eq. uthdtr(2,jth)(1:24)) then
                qzerov = .true.
                do j = 1,jpfcmx
                  if (atheta(j,jth) .ne. 0.) qzerov = .false.
                enddo
                qmatch = .true.
                do j = 1,jpfcmx
                  if ((atheta(j,jth) - atheta(j,jthpr)) .gt. 1.e-6)
     $            qmatch = .false.
                enddo
                if (.not.qzerov .and. .not.qmatch)
     $          ndupl = ndupl + 1
              endif
            endif
          enddo
c
          if (ndupl .gt. 0) then
            j2 = ilnobl(unam1)
            j3 = ilnobl(unam2)
            if (ndupl .eq. 1) then
              write (noutpt,1010) unam1(1:j2),unam2(1:j3)
              write (nttyo,1010) unam1(1:j2),unam2(1:j3)
 1010         format(/' * Error - (EQPT/tprcc) Have found a',
     $        ' duplicate data block on the DATA0 file',
     $        /7x,"for the cc' pair ",a,', ',a,'.')
            else
              ux8 = ' '
              write (ux8,'(i5)') ndupl
              call lejust(ux8)
              j5 = ilnobl(ux8)
              write (noutpt,1020) ux8(1:j5),unam1(1:j2),unam2(1:j3)
              write (nttyo,1020) ux8(1:j5),unam1(1:j2),unam2(1:j3)
 1020         format(/' * Error - (EQPT/tprcc) Have found ',a,
     $        ' duplicate data blocks on the DATA0 file',
     $        /7x,"for the cc' pair ",a,', ',a,'.')
            endif
            nerr = nerr + ndupl
          endif
        endif
c
        if (jthpr.le.0) then
c
c         No data block was found on the DATA0 file.
c         Note that qpdcc(n) is left with a value of .false.
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
 1040   format(/' * Warning - (EQPT/tprcc) Did not find a data',
     $  ' block on the DATA0 file',/7x,'for any of the following',
     $  " cc' pairs:",/)
c
        ncount = 0
        do n = 1,nccpr
          if (.not.qpdcc(n)) then
            ncount = ncount + 1
            i = iccpr(1,n)
            j = iccpr(2,n)
            unam1 = uaqsp(i)
            unam2 = uaqsp(j)
            j2 = ilnobl(unam1)
            j3 = ilnobl(unam2)
            ustr56 = unam1(1:j2) // ', ' // unam2(1:j3)
            j4 = ilnobl(ustr56)
            write (noutpt,1050) ustr56(1:j4)
            write (nttyo,1050) ustr56(1:j4)
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
        nwarn = nwarn + 1
      endif
c
      ncvcc = nccpr - nodatc
      pcvcc = (100.*ncvcc)/float(nccpr)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
c
      end
