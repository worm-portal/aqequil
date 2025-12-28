c eqlpar.h
c
c     Explanations of the dimensioning parameters and equivalent
c     variables may be found in the file eqldef.h.
c
c     Some parameters are also defined and set in the include files
c     eqlj8.h and eqlo8.h. These are involved in the menu-style ("D")
c     input format for Version 8.
c
c-----------------------------------------------------------------------
c
c     Dimensioning parameters and equivalent variables used in both
c     EQ3NR and EQ6 modes (fixed/variable mass of solvent water).
c
c       Primary parameters:
c
          integer iet_par,jet_par,jso_par,net_par,nodb_par,nopg_par,
     $    nopr_par,nopt_par,ntit_par,nxmd_par
c
          parameter(iet_par = 10,jet_par = 4,jso_par = 12,
     $    net_par = 12,nodb_par = 20,nopg_par = 20,nopr_par = 20,
     $    nopt_par = 20,ntit_par = 200,nxmd_par = 100)
c
          integer iapxa_par,ibpxa_par
c
          parameter(iapxa_par = 30,ibpxa_par = 10)
c
c       Secondary parameters:
c
          integer ket_par
c
          parameter(ket_par = iet_par*net_par)
c
c-----------------------------------------------------------------------
c
c     Dimensioning parameters and equivalent variables used in only
c     EQ3NR mode (fixed mass of solvent water).
c
c       Primary parameters:
c
          integer njf_par,nxti_par
c
          parameter(njf_par = 30,nxti_par = 50)
c
c-----------------------------------------------------------------------
c
c     Dimensioning parameters and equivalent variables used in only
c     EQ6 mode (variable mass of solvent water).
c
c       Primary parameters:
c
          integer imch_par,ndct_par,nert_par,nffg_par,nprp_par,
     $    nptk_par,nrct_par,nsrt_par,nttk_par,nxop_par,nxpe_par,
     $    nxrt_par,nssc_par
c
          parameter(imch_par = 4,ndct_par = 4,nert_par = 40,
     $    nffg_par = 20,nprp_par = 100,nptk_par = 3,nrct_par = 40,
     $    nsrt_par = 40,nttk_par = 3,nxop_par = 40,nxpe_par = 100,
     $    nxrt_par = 40,nssc_par = 20)
c
c       Secondary parameters:
c
          integer nprs_par
c
          parameter(nprs_par = 3*nprp_par)
c
c
c     End of INCLUDE file eqlpar.h
c-----------------------------------------------------------------------
