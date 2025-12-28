      subroutine thetck(athetx,jpfcmx,na,nc,nerr,nn,noutpt,nttyo,
     $ n1,n2,n3,unam1,unam2,unam3)
c
c     This suboutine checks data that were read from the S-theta
c     fields of the DATA0 file.
c
c     This suboutine is called by:
c
c       EQPT/rdpz3.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       na     = number of anions in a species triplet
c       nc     = number of cations in a species triplet
c       nn     = number of neutral species in a species triplet
c       n1     = index of the first species in a triplet
c       n2     = index of the second species in a triplet
c       n3     = index of the third species in a triplet
c       unam1  = first name in a species triplet
c       unam2  = second name in a species triplet
c       unam3  = third name in a species triplet
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
      integer jpfcmx
c
      integer noutpt,nttyo
c
      integer na,nc,nn,nerr,n1,n2,n3
c
      character*24 unam1,unam2,unam3
c
      real*8 athetx(jpfcmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j,j2,j3,j4
c
      integer ilnobl
c
      logical qzethx
c
c-----------------------------------------------------------------------
c
c     Check for erroneous data in the thetax data.
c
c
      qzethx = .true.
      do j = 1,jpfcmx
        if (athetx(j) .ne. 0.) qzethx = .false.
      enddo
c
      if (.not.qzethx) then
c
        if (nn .eq. 3) then
c
c         Note: the unique neutral is assumed to be the third
c         species in the triplet.
c
          if (n1.eq.n2 .and. n3.ne.n1) then
c
c           Have an nnn' combination.
c
            j2 = ilnobl(unam1)
            j3 = ilnobl(unam2)
            j4 = ilnobl(unam3)
            write (noutpt,1400) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1400) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1400       format(/' * Error - (EQPT/thetck) Have found an illegal',
     $      ' data block for the',/7x,'species triplet ',a,', ',a,
     $      ', ',a,' among the blocks for',/7x,'"mixture parameters".',
     $      ' Have non-zero input for lambda data in the',
     $      /7x,'"theta" field. Only mu data may be specified (in',
     $      ' the "psi" field)',/7x,'in a block of this type. This',
     $      " is an nnn' combination. Enter the",/7x,'corresponding',
     $      " lambda(nn') data",' in a "single-salt parameters" block,',
     $      /7x,'using the "beta0" field. Leave the'," mu(nnn') data",
     $      ' in the present',/7x,'"mixture parameters" block (in the',
     $      ' "psi" field). Enter the',/7x,"corresponding mu(n'n'n)",
     $      ' data in a second "mixture parameters" block',
     $      /7x,'(in the "psi" field).')
            nerr = nerr + 1
          endif
        endif
c
        if (nn.eq.1 .and. nc.eq.1 .and. na.eq.1) then
c
c         Have an nca combination.
c
          j2 = ilnobl(unam1)
          j3 = ilnobl(unam2)
          j4 = ilnobl(unam3)
          write (noutpt,1410) unam1(1:j2),unam2(1:j3),unam3(1:j4)
          write (nttyo,1410) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1410     format(/' * Error - (EQPT/thetck) Have found an illegal',
     $    ' data block for the',/7x,'species triplet ',a,', ',a,
     $    ', ',a,' among the blocks for',/7x,'"mixture parameters".',
     $    ' Have non-zero input for lambda data in the',/7x,'"theta"',
     $    ' field. Only zeta data may be specified (in the "psi"',
     $    ' field)',/7x,'in a block of this type. This is an nca',
     $    ' combination. Enter the',/7x,'corresponding lambda(nc)',
     $    ' and lambda(na) data in separate',/7x,'"single-salt',
     $    ' parameters" blocks, using the "beta0" fields.',
     $    /7x,'Leave the zeta(nca) data in the present "mixture',
     $    ' parameters"',/7x,'block (in the "psi" field).')
          nerr = nerr + 1
        endif
c
      endif
c
      end
