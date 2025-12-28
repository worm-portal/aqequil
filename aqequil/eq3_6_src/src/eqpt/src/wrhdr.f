      subroutine wrhdr(ikt_asv,ipch_asv,ipcv_asv,jpfc_asv,nap_asv,
     $ narx_asv,nat_asv,nbt_asv,nct_asv,ndata1,ndat1f,ngt_asv,
     $ nlat_asv,nlt_asv,nmt_asv,nmut_asv,npt_asv,nst_asv,ntid_asv,
     $ ntpr_asv,nxt_asv,uakey)
c
c     This suboutine writes a header on the DATA1 and DATA1F files.
c     The header includes the string "data1', followed by the keystring
c     for the type of aqueous species activity coefficient model, and
c     the dimensioning parameters required to read the rest of the data
c     on the DATA1 file. The dimensioning parameters are used by EQ3NR
c     and EQ6.
c
c     This suboutine is called by:
c
c       EQPT/eqpt.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ndata1    = the unit number of the DATA1 file
c       ndat1f    = the unit number of the DATA1F file
c       uakey     = the keystring identifying the type of model for the
c                     aqueous species activity coefficients
c       ikt_asv  = the maximum number of end-member component species
c                    in any solid solution on the data file
c       jpfc_asv = the number of terms (coefficients) in the
c                    temperature function used to represent Pitzer
c                    interaction parameters
c       nap_asv  = maximum number of distinct sets of Pitzer alpha
c                    parameters
c       nat_asv  = the number of aqueous species on the data file
c       nbt_asv  = the number of basis species on the data file
c       nct_asv  = the number of chemical elements on the data file
c       ngt_asv  = the number of gas species on the data file
c       nlt_asv  = the number of pure liquid species on the data file
c       nmt_asv  = the number of gas species on the data file
c       npt_asv  = The number of phases of all types on the data file
c       nlat_asv = the number of members in the set of Pitzer lambda
c                    coefficients
c       nmut_asv = the number of members in the set of Pitzer mu
c                    coefficients
c       nst_asv  = the number of species of all types on the data file
c       ntid_asv = the number of lines in the data file title
c       nxt_asv  = the number of solid-solution phases on the data file
c
c     Principal output:
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
      integer ndata1,ndat1f
c
      integer ikt_asv,ipch_asv,ipcv_asv,jpfc_asv,nap_asv,narx_asv,
     $ nat_asv,nbt_asv,nct_asv,ngt_asv,nlat_asv,nlt_asv,nmt_asv,
     $ nmut_asv,npt_asv,nst_asv,ntid_asv,ntpr_asv,nxt_asv
c
      character(len=8) uakey
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,j3
c
      integer ilnobl
c
      character(len=72) uterm,utermc
      character(len=56) ux56
      character(len=8) udat1
c
c-----------------------------------------------------------------------
c
      uterm(1:48) = '+-----------------------------------------------'
      uterm(49:72) = '------------------------'
      utermc = uterm
      utermc(1:1) = '*'
c
c     Write 'data1' at the top of the DATA1 and DATA1F files.
c
      udat1 = 'data1'
      j2 = ilnobl(udat1)
c
      write (ndata1) udat1
      write (ndat1f,1000) udat1(1:j2)
 1000 format(a)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Next write the keystring for the type of model for the activity
c     coefficients of the aqueous species.
c
      j3 = ilnobl(uakey)
c
      write (ndata1) uakey
      write (ndat1f,1000) uakey(1:j3)
      write (ndat1f,1000) utermc(1:72)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the minimum dimensioning parameters required by the current
c     data file. This will be used by EQ3NR and EQ6 to allocate array
c     sizes for arrays which will contain the data read from this data
c     file.
c
c     Write the standard grid parameters on the output files.
c
      ux56 = 'Dim: chemical elements'
      write (ndata1) ux56
      j2 = ilnobl(ux56)
      write (ndat1f,1000) ux56(1:j2)
      write (ndata1) nct_asv
      write (ndat1f,1010) nct_asv
 1010 format(i5)
c
      ux56 = 'Dim: basis species'
      write (ndata1) ux56
      j2 = ilnobl(ux56)
      write (ndat1f,1000) ux56(1:j2)
      write (ndata1) nbt_asv
      write (ndat1f,1010) nbt_asv
c
      ux56 = 'Dim: total species of all types'
      write (ndata1) ux56
      j2 = ilnobl(ux56)
      write (ndat1f,1000) ux56(1:j2)
      write (ndata1) nst_asv
      write (ndat1f,1010) nst_asv
c
      ux56 = 'Dim: total phases'
      write (ndata1) ux56
      j2 = ilnobl(ux56)
      write (ndat1f,1000) ux56(1:j2)
      write (ndata1) npt_asv
      write (ndat1f,1010) npt_asv
c
      ux56 = 'Dim: aqueous species'
      write (ndata1) ux56
      j2 = ilnobl(ux56)
      write (ndat1f,1000) ux56(1:j2)
      write (ndata1) nat_asv
      write (ndat1f,1010) nat_asv
c
      ux56 = 'Dim: pure minerals'
      write (ndata1) ux56
      j2 = ilnobl(ux56)
      write (ndat1f,1000) ux56(1:j2)
      write (ndata1) nmt_asv
      write (ndat1f,1010) nmt_asv
c
      ux56 = 'Dim: pure liquids'
      write (ndata1) ux56
      j2 = ilnobl(ux56)
      write (ndat1f,1000) ux56(1:j2)
      write (ndata1) nlt_asv
      write (ndat1f,1010) nlt_asv
c
      ux56 = 'Dim: gas species'
      write (ndata1) ux56
      j2 = ilnobl(ux56)
      write (ndat1f,1000) ux56(1:j2)
      write (ndata1) ngt_asv
      write (ndat1f,1010) ngt_asv
c
      ux56 = 'Dim: solid solutions'
      write (ndata1) ux56
      j2 = ilnobl(ux56)
      write (ndat1f,1000) ux56(1:j2)
      write (ndata1) nxt_asv
      write (ndat1f,1010) nxt_asv
c
      ux56 = 'Dim: max. components in a solid solution'
      write (ndata1) ux56
      j2 = ilnobl(ux56)
      write (ndat1f,1000) ux56(1:j2)
      write (ndata1) ikt_asv
      write (ndat1f,1010) ikt_asv
c
      ux56 = 'Dim: lines in the data file title'
      write (ndata1) ux56
      j2 = ilnobl(ux56)
      write (ndat1f,1000) ux56(1:j2)
      write (ndata1) ntid_asv
      write (ndat1f,1010) ntid_asv
c
      ux56 = 'Dim: ranges in the logK temperature grid'
      write (ndata1) ux56
      j2 = ilnobl(ux56)
      write (ndat1f,1000) ux56(1:j2)
      write (ndata1) ntpr_asv
      write (ndat1f,1010) ntpr_asv
c
      ux56 = 'Dim: max. points in a logK temperature grid range'
      write (ndata1) ux56
      j2 = ilnobl(ux56)
      write (ndat1f,1000) ux56(1:j2)
      write (ndata1) narx_asv
      write (ndat1f,1010) narx_asv
c
      ux56 = 'Dim: dH/dP order'
      write (ndata1) ux56
      j2 = ilnobl(ux56)
      write (ndat1f,1000) ux56(1:j2)
      write (ndata1) ipch_asv
      write (ndat1f,1010) ipch_asv
c
      ux56 = 'Dim: dV/dP order'
      write (ndata1) ux56
      j2 = ilnobl(ux56)
      write (ndat1f,1000) ux56(1:j2)
      write (ndata1) ipcv_asv
      write (ndat1f,1010) ipcv_asv
c
      if (uakey(1:8) .eq. 'Pitzer  ') then
        ux56 = 'Dim: distinct sets of Pitzer alpha parameters'
        write (ndata1) ux56
        j2 = ilnobl(ux56)
        write (ndat1f,1000) ux56(1:j2)
        write (ndata1) nap_asv
        write (ndat1f,1010) nap_asv
c
        ux56 = 'Dim: Pitzer lambda parameters'
        write (ndata1) ux56
        j2 = ilnobl(ux56)
        write (ndat1f,1000) ux56(1:j2)
        write (ndata1) nlat_asv
        write (ndat1f,1010) nlat_asv
c
        ux56 = 'Dim: Pitzer mu parameters'
        write (ndata1) ux56
        j2 = ilnobl(ux56)
        write (ndat1f,1000) ux56(1:j2)
        write (ndata1) nmut_asv
        write (ndat1f,1010) nmut_asv
c
        ux56 = 'Dim: coefficients in the Pitzer parameter'
     $  // ' temp. function'
        write (ndata1) ux56
        j2 = ilnobl(ux56)
        write (ndat1f,1000) ux56(1:j2)
        write (ndata1) jpfc_asv
        write (ndat1f,1010) jpfc_asv
      endif
c
      write (ndat1f,1000) utermc(1:72)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
