      subroutine inupt(amua,aslma,ielam,ipbt_asv,jpdblo,jpfc_asv,
     $ nad1,nalpaa,napa_asv,napta,narn1a,narn2a,nerr,nmuta,
     $ nmuta_asv,nmuxa,noutpt,nslta,nslta_asv,nslxa,nsta_asv,
     $ nttyo,palpaa,uspeca,zchara)
c
c     This subroutine reads the parameters for Pitzer's equations from
c     the data file.
c
c     This subroutine is called by:
c
c       EQLIB/indata.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       nad1      = unit number of the data1 file
c       ipbt_asv  = dimensioning variable, the number of Pitzer
c                     alpha coefficients for any cation-anion pairnts
c       jpdblo    = the Pitzer data block organization flag:
c                     -1 = classical (pre-version 8)
c                      0 = new
c       jpfc_asv  = dimensioning variable, the number of coefficients
c                     in a Pitzer parameter temperature functionientsies
c       nmuta_asv = dimensioning variable, the maximum number of species
c                     triplets for which mu coefficients are defined
c       nslta_asv = dimensioning variable, the maximum number of species
c                     pairs for which S-lambda coefficients are defined
c       zchara    = array of electrical charge numbers of the species
c                     read from the data file
c
c     Principal output:
c
c       amua   = coefficients for calculating Pitzer mu parameters
c                  as a function of temperature
c       aslma  = coefficients for calculating Pitzer S-lamnda(n)
c                  parameters as a function of temperature
c       nalpaa = pointer array giving the index of the set of alpha
c                  coeffcients for a given set of S-lambda coefficients
c       nmuta  = number of species triplets for which mu
c                  coefficients are defined
c       nmuxa  = indices of species in triplets for which mu
c                  coefficients are defined
c       nslta  = number of species pairs for which S-lambda
c                  coefficients are defined
c       nslxa  = indices of species in pairs for which S-lambda
c                  coefficients are defined
c       palpaa = array of sets of alpha coefficients
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ipbt_asv,jpfc_asv,napa_asv,nmuta_asv,nslta_asv,nsta_asv
c
      integer noutpt,nttyo
c
      integer nalpaa(nslta_asv),nmuxa(3,nmuta_asv),nslxa(2,nslta_asv)
c
      integer ielam,jpdblo,nad1,napta,narn1a,narn2a,nerr,nmuta,nslta
c
      character(len=48) uspeca(nsta_asv)
c
      real(8) amua(jpfc_asv,nmuta_asv),
     $ aslma(jpfc_asv,0:ipbt_asv,nslta_asv),palpaa(ipbt_asv,napa_asv),
     $ zchara(nsta_asv)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,j,ja,j2,nblock,nmu,ns,nsl
c
      integer ilnobl
c
      logical qx
c
      character(len=80) uline
      character(len=24) unam1,unam2,unam3
      character(len=16) ux16
      character(len=8) uelam,uendit,ux8
c
      real(8) aax,z1,z2
c
c-----------------------------------------------------------------------
c
c     Local variable declarations with global dimensioning.
c     These need not be SAVEd.
c
      real(8), dimension(:), allocatable :: alphai
c
c-----------------------------------------------------------------------
c
      data uendit /'endit.  '/
c
c-----------------------------------------------------------------------
c
      ALLOCATE(alphai(ipbt_asv))
c
c-----------------------------------------------------------------------
c
      if (jpdblo .eq. -1) then
c
        if (ipbt_asv .ne. 2) then
          write (ux8,'(i5)') ipbt_asv
          call lejust(ux8)
          j2 = ilnobl(ux8)
          write (noutpt,1000) ux8(1:j2)
          write (nttyo,1000) ux8(1:j2)
 1000     format(/' * Error - (EQLIB/inupt) Have an illegal value',
     $    ' of ',a,' for the array',/7x,'allocation size variable',
     $    ' ipbt_asv (number of Pitzer alpha parameters',/7x,'for ',
     $    ' any cation-anion or like pair. For the classical Pitzer',
     $    /7x,'data block organization present here, this variable',
     $    ' must have',/7x,'a value of 2.')
          stop
        endif
      endif
c
c-----------------------------------------------------------------------
c
      if (jpdblo .eq. -1) then
c
c       In the "classical" Pitzer data block organization, the E-lambda
c       (E-theta) flag may be "on" or "off", as specified on the data
c       file. Read and decode this flag.
c
        uelam = ' '
        read (nad1) uline
        uelam(1:8) = uline(17:24)
        if (uelam(1:3) .eq. 'off') then
          ielam = -1
        elseif (uelam(1:3) .eq. 'on ') then
          ielam = 0
        else
          j2 = ilnobl(uelam)
          write (noutpt,1010) uelam(1:j2)
          write (nttyo,1010) uelam(1:j2)
 1010     format(/' * Error - (EQLIB/inupt) Have read an unrecognized',
     $    /7x,'value of "',a,'" for the E-lambda (E-theta) flag from',
     $    /7x,'the data file. Allowed values are "on" and "off".')
          stop
        endif
      else
c
c       In the "new" Pitzer data block organization, the E-lambda
c       flag is always "on".
c
        ielam = 0
      endif
c
c-----------------------------------------------------------------------
c
c     Read the species pairs for which S-lambda coefficients are
c     defined. These data are comprised in two superblocks.
c
      ja = 0
      nsl = 0
      do nblock = 1,2
c
c       Read a block.
c       Read the names of the species in a pair.
c
   15   read (nad1) unam1,unam2
c
c       Test for end of superblock.
c
        if (unam1(1:6) .eq. uendit(1:6)) go to 130
c
c       Increment the S-lambda species pair counter.
c
        nsl = nsl + 1
c
c       Calling sequence substitutions:
c         narn1a for nrn1a
c         narn2a for nrn2a
c         unam1 for unam
c
        call srchn(narn1a,narn2a,ns,nsta_asv,unam1,uspeca)
        nslxa(1,nsl) = ns
        z1 = zchara(ns)
c
c       Calling sequence substitutions:
c         narn1a for nrn1a
c         narn2a for nrn2a
c         unam2 for unam
c
        call srchn(narn1a,narn2a,ns,nsta_asv,unam2,uspeca)
        nslxa(2,nsl) = ns
        z2 = zchara(ns)
c
c       Read the temperature function coefficients for calculating
c       Pitzer S-lambda(n) parameters. Read also the associated
c       Pitzer alpha parameters.
c
        if (jpdblo .eq. -1) then
c
c         Classical Pitzer data block organization.
c
c         Read the lowest-order coefficients, which are here
c         the 25C values of the corresponding interaction
c         coefficients.
c
          read (nad1) (aslma(1,i,nsl), i = 0,2)
c
c         Read the associated alpha coefficients.
c
          read (nad1) (alphai(i), i = 1,2)
c
c         Read the higher-order coefficients, which here are
c         first and second-order temperature derivatives.
c
          do i = 0,2
            read (nad1) aslma(2,i,nsl),aslma(3,i,nsl)
          enddo
        else
c
c         New Pitzer data block organization.
c
c         Read the associated alpha coefficients.
c
          if ((z1*z2) .lt. 0.) then
c
c           Read a larger set of data for cation-anion pairs.
c
            do i = 1,ipbt_asv
              read (nad1) ux16,alphai(i)
            enddo
c
            do i = 0,ipbt_asv
              read (nad1) ux16
              do j = 1,jpfc_asv
                read (nad1) aslma(j,i,nsl)
              enddo
            enddo
          else
c
c           Read a more limited set of data for all other pair types.
c
            read (nad1) ux16
            do j = 1,jpfc_asv
              read (nad1) aslma(j,0,nsl)
            enddo
          endif
        endif
c
c       Read the block terminator.
c
        read (nad1) uline(1:72)
c
c       Test for alpha coefficient set already in the palpaa array.
c
        do j = 1,ja
          qx = .true.
          do i = 1,ipbt_asv
            aax = abs(palpaa(i,j) - alphai(i))
            if (aax .gt. 1.e-12) then
              qx = .false.
              go to 110
            endif
          enddo
  110     continue
c
          if (qx) then
c
c           Found the index to an existing alpha pair. Store this
c           index, then go process the data for another species pair.
c
            nalpaa(nsl) = j
            go to 15
          endif
        enddo
c
c       Not in the set; add it and store the index for this species
c       pair.
c
        ja = ja + 1
c
        do i = 1,ipbt_asv
          palpaa(i,ja) = alphai(i)
        enddo
        nalpaa(nsl) = ja
c
c       Go process the data for another species pair.
c
        go to 15
c
  130   continue
      enddo
c
      nslta = nsl
      napta = ja
c
      DEALLOCATE(alphai)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the species triplets for which mu coefficients are
c     defined. These data are comprised in two superblocks.
c
      nmu = 0
      do nblock = 1,2
c
c       Read a block.
c       Read the names of the species in a triplet.
c
  150   read (nad1) unam1,unam2,unam3
c
c       Test for end of superblock.
c
        if (unam1(1:6) .eq. uendit(1:6)) go to 200
c
c       Increment the mu species triplet counter.
c
        nmu = nmu + 1
c
c       Calling sequence substitutions:
c         narn1a for nrn1a
c         narn2a for nrn2a
c         unam1 for unam
c
        call srchn(narn1a,narn2a,ns,nsta_asv,unam1,uspeca)
        nmuxa(1,nmu) = ns
c
c       Calling sequence substitutions:
c         narn1a for nrn1a
c         narn2a for nrn2a
c         unam2 for unam
c
        call srchn(narn1a,narn2a,ns,nsta_asv,unam2,uspeca)
        nmuxa(2,nmu) = ns
c
c       Calling sequence substitutions:
c         narn1a for nrn1a
c         narn2a for nrn2a
c         unam3 for unam
c
        call srchn(narn1a,narn2a,ns,nsta_asv,unam3,uspeca)
        nmuxa(3,nmu) = ns
c
c       Read the coefficients for calculating Pitzer mu
c       parameters as a function of temperature.
c
        if (jpdblo .eq. -1) then
c
c         Classical Pitzer data block organization.
c
          read (nad1) amua(1,nmu),amua(2,nmu),amua(3,nmu)
        else
c
c         New Pitzer data block organization.
c
          read (nad1) ux16
          do j = 1,jpfc_asv
            read (nad1) amua(j,nmu)
          enddo
        endif
c
c       Read the block terminator.
c
        read (nad1) uline(1:72)
c
        go to 150
c
  200   continue
      enddo
c
      nmuta = nmu
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
