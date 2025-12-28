      subroutine indath(ikta_asv,ipbt_asv,ipch_asv,ipcv_asv,
     $ jpfc_asv,nad1,napa_asv,narx_asv,nata_asv,nbta_asv,
     $ ncta_asv,ngta_asv,nlta_asv,nmta_asv,npta_asv,nmuta_asv,
     $ noutpt,nslta_asv,nsta_asv,ntid_asv,ntpr_asv,nttyo,
     $ nxta_asv,udakey)
c
c     This subroutine reads the header section of the data1 file. This
c     section consists of a record containing the string 'data1' (to
c     ensure that the file is indeed a data1 file), a record containing
c     the keystring for the activity coefficient model for aqueous
c     species to which this data file corresponds, and the array
c     allocation size variables required to allocate sufficient array
c     space to store the rest of the data on this data file.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       nad1   = unit number of the data1 file
c       noutpt = unit number of the output file
c       nttyo  = unit number of the screen file
c
c     Principal output:
c
c       nbta_asv  = the number of basis species on the data file
c       ncta_asv  = the number of chemical elements on the data file
c       ntid_asv  = the number of lines in the data file title
c
c       ikta_asv  = the maximum number of end-member component species
c                     in any solid solution on the data file
c       nata_asv  = the number of aqueous species on the data file
c       ngta_asv  = the number of gas species on the data file
c       nlta_asv  = the number of pure liquids on the data file
c       nmta_asv  = the number of pure minerals on the data file
c       npta_asv  = The number of phases of all types on the data file
c       nsta_asv  = the number of species of all types on the data file
c       nxta_asv  = the number of solid-solution phases on the data
c                     file
c
c       napa_asv  = the number of distinct sets of Pitzer alpha values
c       nmuta_asv = the number of triplets of ions on the data file
c                     for which distinct Pitzer mu coefficients
c                     are defined
c       nslta_asv = the number of pairs of ions on the data file
c                     for which distinct Pitzer S-lambda coefficients
c                     are defined
c
c       ipch_asv = the maximum order for pressure corrections to
c                    enthalpy functions
c       ipcv_asv = the maximum order for pressure corrections to
c                    volume functions; the maximum order for pressure
c                    corrections to log K and other Gibbs-energy-based
c                    functions is one greater than this
c       ipbt_asv  = the maximum number of Pitzer alpha parameters
c                     for any species pair
c       jpfc_asv  = the number of coefficients in the Pitzer parameter
c                    temperature function
c       narx_asv = maximum number of coefficients per temperature range
c       ntpr_asv = number of temperature ranges
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nad1,noutpt,nttyo
c
      integer ikta_asv,ipbt_asv,napa_asv,nata_asv,nbta_asv,ncta_asv,
     $ ngta_asv,nlta_asv,nmta_asv,npta_asv,nmuta_asv,nslta_asv,nsta_asv,
     $ ntid_asv,nxta_asv
c
      integer ipch_asv,ipcv_asv,jpfc_asv,narx_asv,ntpr_asv
c
      character(len=8) udakey
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      character(len=56) ux56
      character(len=8) ux8
c
c-----------------------------------------------------------------------
c
      write (noutpt,1000)
      write (nttyo,1000)
 1000 format(/' Reading the data1 file header section ...')
c
c     Rewind the data file. The unit number is nad1.
c
      rewind nad1
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the file header string.
c
      read (nad1) ux8
c
      if (ux8(1:5) .ne. 'data1') then
        write (noutpt,1010) ux8
        write (nttyo,1010) ux8
 1010   format(/' * Error - (EQLIB/indath) Have wrong file header',
     $  ' "',a5,'"',/7x,'on the data1 file. The first five characters',
     $  ' must be "data1".')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read key string indicating the type of option supported by the
c     data file for treating activity coefficients of aqueous species.
c     This will be used to test the model specified on the input file
c     for consistency with the data file. That test is made in
c     EQLIB/cdakey.f.
c
      read (nad1) udakey
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the array dimension required for chemical elements.
c
      read (nad1) ux56
      read (nad1) ncta_asv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the array dimension required for basis species.
c
      read (nad1) ux56
      read (nad1) nbta_asv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the array dimension required for total species.
c
      read (nad1) ux56
      read (nad1) nsta_asv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the array dimension required for total phases.
c
      read (nad1) ux56
      read (nad1) npta_asv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the array dimension required for aqueous species.
c
      read (nad1) ux56
      read (nad1) nata_asv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the array dimension required for pure minerals.
c
      read (nad1) ux56
      read (nad1) nmta_asv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the array dimension required for pure liquids.
c
      read (nad1) ux56
      read (nad1) nlta_asv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the array dimension required for gas species.
c
      read (nad1) ux56
      read (nad1) ngta_asv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the array dimension required for solid-solution phases.
c
      read (nad1) ux56
      read (nad1) nxta_asv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the array dimension required for solid-solution components
c     (maximum number per solid solution).
c
      read (nad1) ux56
      read (nad1) ikta_asv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the array dimension required for the data file title.
c
      read (nad1) ux56
      read (nad1) ntid_asv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the array dimension required for the temperature ranges in
c     the logK temperature grid.
c
      read (nad1) ux56
      read (nad1) ntpr_asv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the array dimension required for the points in a temperature
c     range on the logK temperature grid.
c
      read (nad1) ux56
      read (nad1) narx_asv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the array dimension required for the dH/dP order represented
c     on the data file.
c
      read (nad1) ux56
      read (nad1) ipch_asv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the array dimension required for the dV/dP order represented
c     on the data file.
c
      read (nad1) ux56
      read (nad1) ipcv_asv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the array dimensions for parameters needed for activity
c     coefficient models.
c
      if (udakey(1:8) .eq. 'Pitzer  ') then
c
c       Read the array dimension required for the sets of Pitzer
c       alpha coefficients.
c
        read (nad1) ux56
        read (nad1) napa_asv
c
c       Read the array dimension required for the set of Pitzer
c       S-lambda coefficients.
c
        read (nad1) ux56
        read (nad1) nslta_asv
c
c       Read the array dimension required for the set of Pitzer mu
c       coefficients.
c
        read (nad1) ux56
        read (nad1) nmuta_asv
c
c       Read the number of coefficients for the Pitzer parameter
c       temperature function.
c
        read (nad1) ux56
        read (nad1) jpfc_asv
      else
c
c       Minimum values.
c
        napa_asv = 1
        nslta_asv = 1
        nmuta_asv = 1
        jpfc_asv = 1
      endif
c
c     Set the number of parameters in a Pitzer alpha set.
c
      if (udakey(1:8) .eq. 'Pitzer  ') then
        ipbt_asv = 2
      else
        ipbt_asv = 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The following is a bit of nonsense so compiler warninings will
c     not be generated that ux56 is not used.
c
      ux8 = ux56(1:8)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
