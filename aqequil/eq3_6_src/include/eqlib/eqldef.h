c eqldef.h
c
c     Definitions of dimensioning paramaters and equivalent variables.
c     The parameters themselves are set in the include file eqlpar.h.
c
c     Some parameters are also defined and set in the include files
c     eqlj8.h and eqlo8.h. These are involved in the menu-style ("D")
c     input format for Version 8.
c
c-----------------------------------------------------------------------
c
c     Those used in both EQ3NR and EQ6 modes (fixed/variable mass of
c     solvent water).
c
c
c       Primary parameters/variables:
c
c         iapxpa = iapxmx   Leading dimension of the apx array, which
c                           contains solid solution activity
c                           coefficient parameters other than
c                           site-mixing paramters
c
c         ibpxpa = ibpxmx   Leading dimension of the bpx array,
c                           which contains solid solution activity
c                           coefficient parameters of the site-mixing
c                           variety
c
c         ietpar = ietmax   Maximum number of species in a distinct
c                           site of a generic ion exchange phase,
c                           including the bare site species
c
c         iktpar = iktmax   Maximum number of end-member components
c                           (species) in a solid solution
c
c         ipchpa = ipchmx   Maximum order of pressure corrections
c                             to enthalpy and volume functions
c
c         ipcvpa = ipcvmx   Maximum order of pressure corrections
c                             to enthalpy and volume functions
c
c         jetpar = jetmax   Maximum number of distinct sites of
c                           a generic ion exchange phase
c
c         jodbpa = jodbmx   Maximum number of choices per iodb
c                             debugging print option switch
c
c         jopgpa = jopgmx   Maximum number of choices per iopg activity
c                             coefficient option switch
c
c         joprpa = joprmx   Maximum number of choices per iopr print
c                             option switch
c
c         joptpa = joptmx   Maximum number of choices per iopt option
c                             switch
c
c         jsopar = jsomax   Maximum number of solid solution mixing
c                           laws
c
c         nappar = napmax   Number of distinct sets of Pitzer
c                           alpha coefficients
c
c         narxpa = narxmx   Maximum number of interpolating
c                           coefficients used to describe the
c                           temperature dependence of a quantity
c                           (such as the log K of a reaction)
c                           in any given temperature range
c
c         natpar = natmax   Maximum number of aqueous species
c
c         naztpa = naztmx   Maximum value of the max norm of the
c                           electrical charge numbers of the aqueous
c                           species (nazmmx = -naztmx); used with
c                           Pitzer's equations
c
c         nbtpar = nbtmax   Maximum number of species in a basis set
c
c         nctpar = nctmax   Maximum number of chemical elements
c
c         netpar = netmax   Maximum number of generic ion exchange
c                           phases
c
c         ngtpar = ngtmax   Maximum number of gas species
c
c         nltpar = nltmax   Maximum number of pure liquid
c                           species/phases
c
c         nmtpar = nmtmax   Maximum number of pure mineral
c                           species/phases
c
c         nmutpa = nmutmx   Number of Pitzer mu coefficients
c
c         nodbpa = nodbmx   Number of iodb debugging print option
c                           switches
c
c         nopgpa = nopgmx   Number of iopg activity coefficient option
c                           switches
c
c         noprpa = noprmx   Number of iopr print control option
c                           switches
c
c         noptpa = noptmx   Number of iopt model option switches
c
c         nptpar = nptmax   Maximum number of phases
c
c         nsltpa = nsltmx   Number of Pitzer S-lambda function
c                           coefficients
c
c         nstpar = nstmax   Maximum number of species of any kind
c
c         ntfxpa = ntfxmx   Maximum number of aqueous species which
c                           contribute to the working definition of
c                           total alkalinity (HCO3-CO3-OH or extended)
c
c         ntf1pa = ntf1mx   Maximum number of aqueous species which
c                           contribute to total HCO3-CO3-OH alkalinity
c
c         ntf2pa = ntf2mx   Maximum number of aqueous species which
c                           contribute to the total extended alkalinity,
c                           defined as the HCO3-CO3-OH alkalinity with
c                           contributions from other species such as
c                           acetate, other organics, borate, silicate,
c                           phosphate, and transition metal-hydroxy
c                           complexes
c
c         ntitpa = ntitmx   Maximum number of lines per title on an
c                           input file or data file
c
c         ntprpa = ntprmx   Maximum number of temperature ranges for
c                           describing the temperature dependence of
c                           a quantity (such as the log K of a
c                           reaction) by means of interpolating
c                           polynomials
c
c         nvet_par = nvetmx   Maximum number of valid exchange models
c                           for generic ion exchange phases
c
c         nxmdpa = nxmdmx   Maximum number of species/reactions that
c                           can be affected by the suppress/alter
c                           options
c
c         nxtpar = nxtmax   Maximum number of solid solution phases
c
c
c       Secondary parameters/variables:
c
c         kpar   =  kmax    Maximum dimension of the Jacobian matrix
c                           (WARNING - could set equal to nbtpar in
c                           EQ3NR, but need to set equal to 2*nbtpar
c                           in EQ6)
c
c         ndrspa = ndrsmx   size of array for storing reaction
c                           coefficients, expected average
c                           number*nstpar
c
c         nesspa = nessmx   size of array for storing elemental
c                           composition coefficients, expected
c                           average number*nstpar
c
c         nmxpar = nmxmax   Number of representations of Pitzer mu
c                           coefficients in the nmxx array; can't
c                           be greater than 3*nmutpa
c
c         nstspa = nstsmx   Size of array for storing stoichiometric
c                           mass balance coefficients, expected
c                           average number*nstpar
c
c         nsxpar = nsxmax   Number of representations of sets of
c                           Pitzer S-lambda function coefficients
c                           in the nsxx array; can't be greater than
c                           2*nsltpa
c
c-----------------------------------------------------------------------
c
c     Those used in only EQ3NR mode (fixed mass of solvent water).
c
c
c       Primary parameters/variables:
c
c         njfpar = njfmax   Largest jflag option value
c
c         nxtipa = nxtimx   Maximum number of non-aqueous phases for
c                           which compositions can be read from the
c                           EQ3NR input file
c
c
c       Secondary parameters/variables:
c
c         nbtpa1 = nbtmx1   Maximum number of basis species + 1
c
c         nxicpa = nxicmx   Maximum number of non-aqueous phase
c                           species for which mole fractions can be
c                           read from the EQ3NR input file
c
c
c-----------------------------------------------------------------------
c
c     Those used in only EQ6 mode (variable mass of solvent water).
c
c
c       Primary parameters/variables:
c
c         imchpa = imchmx   Maximum number of terms in a rate law
c                           expression for a chemical reaction
c
c         ndctpa = ndctmx   Maximum number of species whose activities
c                             can appear in any term in a rate law
c                             expression for a chemical reaction
c
c         nertpa = nertmx   Maximum number of generic ion exchanger
c                           irreversible reactants/reactions
c
c         nffgpa = nffgmx   Maximum number of gases whose fugacities
c                           may be fixed
c
c         nordpa = nordmx   Maximum order of finite difference
c                           expressions or the equivalent truncated
c                           Taylor's series
c
c         nprspa = nprsmx   Maximum number of species in the
c                           physically removed system that can be
c                           read from the input file
c
c         nptkpa = nttkmx   Maximum number of pressure tracking
c                           coefficients
c
c         nrctpa = nrctmx   Maximum number of specified irreversible
c                           reactants/reactions
c
c         nttkpa = nttkmx   Maximum number of temperature tracking
c                           coefficients
c
c         nsrtpa = nsrtmx   Maximum number of special irreversible
c                           reactants/reactions that can be defined
c                           on the input file
c
c         nxoppa = nxopmx   Maximum number of subset-selection
c                           options for suppressing mineral phases
c
c         nxpepa = nxpemx   Maximum number of specified exceptions
c                           to the mineral subset-selection
c                           suppression options
c
c         nxrtpa = nxrtmx   Maximum number of solid solution
c                           irreversible reactants/reactions
c
c         nsscpa = nsscmx = Maximum number of setscrew parameters
c
c
c       Secondary parameters/variables:
c
c         npetpa = npetmx   Maximum number of phases in the ES
c
c         nrdpa1 = nrdmx1   Maximum order + 1
c
c
c     End of INCLUDE file eqldef.h
c-----------------------------------------------------------------------
