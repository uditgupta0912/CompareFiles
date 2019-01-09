!This source file contains modules with all the global variables.
!*********************************************************************************************************
module global_all
!This module contains variables of TYPE 1. These are global variables not specific to any reactor model.
!*********************************************************************************************************
  implicit none

  integer  ::   irxtr, nruns, ktt=1, iflow !Reactor type, number of runs, number of energy balances
  integer  ::   kgmax, ksmax, kbmax, kmax ! # of gas, surface, bulk, and total species
  integer  ::   ngrxn, nsrxn, nrxn  ! # of gas, surface, and total reactions
  integer  ::   ngphases, nsphases, nbphases, nphases  ! # of phases (gas, surface, bulk, total)
  integer, allocatable            ::   site_type_lastID(:) !Identity of inerts and site types
  integer, allocatable            ::   site_type_firstID(:) !Identity of first species in each phase (not bulk)
  integer, allocatable :: SpecTot(:)  !Number of species in each phase
  character(len=15), allocatable :: PhaseNames(:)
  integer  ::   nelm  !Number of elements in the surface mechanism
  integer  ::   leniwk, lenrwk, lencwk    ! CHEMKIN work array dimensions
  integer  ::   lenisk, lenrsk, lencsk    ! SURFACE-CHEMKIN work array dimensions
  integer  ::   maxtp !Number of temperatures for fitting NASA polynomials
  double precision, parameter  ::  patm=1.01325e6        ! atm --> dyne/cm2
  integer, allocatable     ::   ickwrk(:), iskwrk(:), imcwrk(:)  ! Integer work arrays
  double precision, allocatable   ::   rckwrk(:), rskwrk(:), rmcwrk(:)  ! Real work arrays
  double precision, allocatable   ::   rckwrk_save(:), rskwrk_save(:)   ! Saved real work arrays
  double precision, allocatable   ::   sden(:), wt(:)      ! Site density and molecular weights
  integer, allocatable     ::   stoich_matrix(:,:) !Stoichiometric coefficients
  integer, allocatable :: stoich_matrix_sparse(:,:) !Sparse representation of the stoichiometric matrix
  character(len=80), allocatable :: rxn_str(:)
  double precision, allocatable   ::   rka(:), rb(:), re(:)    ! Surface: A, Beta, Eact
  double precision, allocatable   ::   rkag(:), rbg(:), reg(:)    ! Gas: A, Beta, Eact
  double precision, allocatable   ::   A6(:,:), A7(:,:) !Original thermo constants
  double precision, allocatable   ::   DH_save(:), EA_save(:) !Original DH, EA for /surface/ reactions
  integer, allocatable     ::   kocc(:)          ! Site occupancy
  double precision, allocatable      ::      cov_matrix(:,:), thresh_cov(:,:), StatpQ(:), cov_matrix2(:,:), cov_matrix3(:,:), thresh_cov2(:,:)
  double precision, allocatable      ::      omega(:) !proximity factor
  double precision      ::      Rgas_cgs, Rgas_mks, kcal_to_J, kcal_to_erg, pi
  data                          Rgas_cgs/8.31451d7/
  data                          Rgas_mks/8.31451d0/
  data                          kcal_to_J/4184.0d0/
  data                          kcal_to_erg/4184.0d7/
  data                          pi/3.141592653589793d0/
  double precision, parameter :: Rgas_kcal=1.987E-3
  double precision, parameter :: Rgas_cal=1.987
  character(len=16), allocatable   ::   cckwrk(:), cskwrk(:)      ! CHEMKIN character arrays
  character(len=16), allocatable   ::   knams(:)        ! Names of species
  character(len=16), allocatable   ::   enams(:)        ! Names of elements
  integer, allocatable :: elem_comp(:,:)  !Number of each element in each species
  integer               ::     IrPAR,IrKTHM, NCP2, NCP2T,i_tube_1               !integers needed to edit work arrays

!The following variables are for use with the BEP option
  double precision, allocatable      ::   BEP_def(:,:), BEP_def_save(:,:)
  integer, allocatable               ::   BEP_rxn(:,:)
  integer               ::   nBEP
  logical, allocatable :: lEA_user(:)

!The following variables/parameters are for the scaling relation code
  integer :: nscale !The actual number of scaling relations
  integer :: natoms_scale !The number of binding energy descriptors
  !The array scale_bind_mode contains the binding modes for the scaling relations (positions 2:nscale+1).
  !The first position is a flag indicating whether the scaling relations should be used for that
  !species. The first column of scale_atom contains the actual atomic species index, while the second column
  !tells how the atomic binding energies map onto the correlations (that is, which binding energy is used with
  !which correlation). The code assumes that the reference and target binding energies must map in the same fashion.
  integer, allocatable :: scale_bind_mode(:,:)
  !Contains the slopes for the scaling relations and the species IDs for the corresponding
  !atoms. Also a saved version for restoring after a sensitivity analysis.
  double precision, allocatable :: scale_slope(:,:), scale_slope_save(:,:)
  double precision, allocatable :: scale_ref(:) !The mapping of the reference energies to the atom IDs
  double precision, allocatable :: scale_targ(:,:) !The current binding energy coordinate
  !The following is used to approximately correct the adsorbate interactions for
  !changes in atomic binding energy.
  double precision, allocatable :: cov_factor(:)
  logical :: cov_factor_adjust
  integer :: nBE_coords !The total number of binding energy coordinates
  integer :: BE_coord_idx !The current binding energy coordinate
  double precision, allocatable :: BE_coords(:,:) !Array for storing the BE coordinates

  !This is the current revision number of the code. Update this for every revision!
  integer, parameter :: revno=116

!*********************************************************************************************************
end module global_all
!*********************************************************************************************************

!*********************************************************************************************************
module file_handles
!This module contains parameter declarations for the input/output LU numbers in a centralized location.
!This is done to reduce/eliminate the chances that a file handle will be mistakenly used more than once.
!The output file handles should be numbered sequentially to facilitate closing all of them with minimal
!effort.
!*********************************************************************************************************

  implicit none

  !Input files
  integer, parameter :: linc=25, linksk=26, linkmc=35, lout=6    ! Link file numbers
  integer, parameter :: itubeinp=300, iTflowinp=301, iomegainp=302, iEAginp=303, iEAsinp=304
  integer, parameter :: itubemoleinp=305, itolinp=306, iBEPinp=307, iScaleinp=308, iDOEinp=309
  integer, parameter :: itubeCOVinp=310, iStatpQinp=311, iSAinp=312, iGAinp=313
  integer, parameter :: iSAstateinp=314, iPRNGstateinp=315
  character(len=*), parameter :: cklink='INP.d/cklink', sklink='INP.d/sklink'
  integer, parameter :: imvtHinp=316, imvtSinp=317, imvninp=318

  !Output files
  integer, parameter :: verout=500
  integer, parameter :: itranmatout=501
  integer, parameter :: ifiletuberpa=502, ifiletrpa=503, irpavisout=504  !RPA files
  integer, parameter :: ifiletubesen=505, isenSout=506, isenHout=507  !SA files
  integer, parameter :: itubeconvout=508, ituberestartout=509 !Conversion/restart composition
  integer, parameter :: irxnratess=543 !Steady state reaction rates
    !Thermodynamic and kinetic parameters
  integer, parameter :: iDHRTout=510, iDSRout=511, iDGRTout=512, ikout=513, ibetaout=514, iEARTout=515, iEQKCout=516
  integer, parameter :: iHRTfout=517, iSRfout=518
  integer, parameter :: ielemcomp=544 !Elemental compositions
    !Kinetic analysis
  integer, parameter :: ituberateout=519
  integer, parameter :: icovtraout=520, igasmoletraout=522, igasmasstraout=524, isdottraout=526, igdottraout=528, imassbaltraout=530
  integer, parameter :: icovssout=521, igasmolessout=523, igasmassssout=525, isdotssout=527, igdotssout=529, imassbalssout=545
    !Global sensitivity analysis
  integer, parameter :: iHfpertout=531, iSfpertout=532, iApertout=533, iBpertout=534
  integer, parameter :: iEApertout=535, iBEPpertout=536, iLSRpertout=537
  integer, parameter :: iSAstateout=538, iPRNGstateout=539
  integer, parameter :: iwSAout=540, irSAout=541, iStoichout=542

  !These constants are used to rapidly close all open output files in a single do-loop
  integer, parameter :: output_first=verout
  integer, parameter :: output_last=ielemcomp

  !These constants are used to rapidly flush all open SA output files in a do-loop
  !Make sure the GSA output file unit numbers are numbered sequentially!
  integer, parameter :: SAout_first=iHfpertout
  integer, parameter :: SAout_last=irSAout

!*********************************************************************************************************
end module file_handles
!*********************************************************************************************************

!*********************************************************************************************************
module global_tube1
!This module contains variables of TYPE 2. These are global variables specific to the reactor model.
!*********************************************************************************************************
  implicit none
  integer  ::  neqns
  integer  ::  indext, indexrho
  integer     ::  nnodes, ncalls, ntdec
  double precision :: ndec
  logical ::  lomega
  integer     ::  iiso
  logical :: lstp, lBEP, lEA, lDOE, lStatpQ, NonNeg, MultiInput, lGA
  logical :: lcov, HideSAstdout, ltol
  integer :: mu, ml, iSolver
  logical :: gas_chem, surf_chem, ltra
  double precision  ::  p, tin, abyv, trise
  double precision  ::  text, aextbyv, htc,ramp
  integer     ::  iScale
  integer     ::  mrpa, isenbrute, isen_gas, isen_surf,itpd
  logical     :: lsenbrute
  !mrpa = 1 for surface chemistry, 2 for gas chemistry, 3 for both, and 4 for RPA at a specified T
  logical :: verbose_rpa  !Should we write the RPA at *every* reactor node (T) or just the exit (F)
  LOGICAL     ::  trpa_called = .FALSE.
  LOGICAL :: SA_called = .FALSE.
  logical :: inco_first_call=.true.
  double precision  ::  rdiv, delS_SA, delH_SA, trpa
  double precision, allocatable :: tvec(:)
  double precision, allocatable  ::  atol(:), rtol(:)
  integer     ::  liw, lrw, lwp
  double precision  ::  Tref_beta, omega_default
  INTEGER :: MARI_index      !Used for controlling which surface coverage is printed to the screen
                             !during solution of ODEs
  INTEGER :: conv_index      !Specifies which species is the 'reactant' for calculating conversion

  !Allocatable arrays for the ID numbers of species and reactions in the input files.
  !This allows for input files to be specified using only the non-trivial entries and
  !independent of order. This reduces/eliminates the odds of misassigning input values.
  integer, allocatable :: iEAg_rxn(:), iEAs_rxn(:)
  integer, allocatable :: iMole_spec(:), iTol_spec(:)

  !Allocatable arrays for storing run-specific parameters
  double precision, allocatable :: TPFSV(:,:) !/T/emperature, /P/ressure, /F/lowrate, /S/urface-to-/V/olume ratio conditions
  double precision, allocatable :: atol_all(:,:), rtol_all(:,:) !Tolerances
  double precision, allocatable :: EAg_all(:,:), EAs_all(:,:) !Activation energies
  double precision, allocatable :: act_all(:,:) !Compositions

!*********************************************************************************************************
end module global_tube1
!*********************************************************************************************************

!*********************************************************************************************************
module update_tube2
!This module contains variables of TYPE 3. These are updating variables specific to the reactor model.
!*********************************************************************************************************
  implicit none

  double precision, allocatable  ::  actin(:), act(:)
  double precision, allocatable  ::  xin(:), x(:)
  double precision, allocatable  ::  yin(:), y(:)
  double precision, allocatable  ::  win(:), wentr(:), wexit(:)
  double precision :: ventr, vexit
  double precision  ::  fluxmass, velo, rhoin
  double precision, allocatable :: restart_save(:,:)
  double precision, allocatable  ::  gdot(:)
  double precision, allocatable  ::  sdot(:), sitdot(:)
  double precision, allocatable  ::  rop(:), ropfwd(:)
  double precision, allocatable  ::  rpa(:,:), rpag(:,:)
  double precision, allocatable  ::  rxn_rates_ss(:,:,:,:)  !Steady state reaction rates (net, forward)
  double precision, allocatable  ::  parasurf(:), paragas(:)
  double precision, allocatable :: FIM_SA_matrix(:,:)
  integer, allocatable     ::  iwork(:)
  double precision, allocatable  ::  rwork(:)
  integer     ::  idid, jac, restarts
  integer :: idid_restart_count_max
  integer     ::  info(20)
  double precision  ::  zout, ttout, zs, tts, zouts, ttouts, rtemp, rttemp
  integer     ::  Nouter, Ntouter, Ninner, Ntinner
  double precision :: tprint, zprint
  integer     ::  itube_restart, i_tube
  double precision  ::  fluxmass_inlet
  double precision  ::  rlen, rtime, abstol, reltol
  type        :: perturbations ! This is a structure for storing data about a single perturbed parameter
      integer          :: param    ! Paramter type flag
      integer          :: id       ! Species/reaction/correlation ID
      double precision :: value(2) ! Amount of perturbation
  end type
  type        :: runs          ! This is a structure for storing multiple sets of perturbed parameters
      integer          :: num_param
      type(perturbations), allocatable :: experiment(:)
  end type
  integer, allocatable :: thermo_basis(:)  ! Allocated in inco_DOE
  double precision, allocatable :: trans_mat(:,:) !Transformation matrix for gas basis to surface species

!*********************************************************************************************************
end module update_tube2
!*********************************************************************************************************

!*********************************************************************************************************
module globalSA
!This module contains variables related to the global sensitivity analysis.
!*********************************************************************************************************

  implicit none

  !Parameter definitions for possible cases
  integer, parameter :: max_SA_param_types=7
    !NOTE: thermochemistry should come first as activation energies may be estimated from it!!
    !H = enthalpy, S = entropy, A = pre-exponential, B = beta, EA = act energy
    !BEP = BEP correlation, LSR = LSR correlation
    !The maximum index value below should be /.le./ the maximum number of
    !parameters max_SA_param_types.
  integer, parameter :: globalSA_H=1, globalSA_S=2, globalSA_A=3, globalSA_B=4, globalSA_EA=5
  integer, parameter :: globalSA_BEP=6, globalSA_LSR=7

  !General
  integer :: nSAruns, nSAruns_old, nBEcoords_old, nexp_old
  logical :: perturb_param(max_SA_param_types)=.false.
  logical :: lDFTH_SA=.false., lGAH_SA=.false., lLSR_SA=.false. !Enthalpy flags
  logical :: lDFTS_SA=.false., lGAS_SA=.false.  !Entropy flags
  logical :: lDFTEA_SA=.false., lBEP_SA=.false. !Activation energy flags
  logical :: lmvtH_SA=.false., lmvtS_SA=.false., lmvn_EA=.false.  !Multivariate distribution flags
  logical :: fix_EA, rel_pert !Fix E_A/E_TS and relative perturbations on/off
  logical :: model_solve
  logical, allocatable :: GA_species(:) !Specifies whether a species thermochemistry is from GA or not
  integer :: iter_write_freq, write_precision
  logical :: SA_restart=.false., restart_on_failure, flush_buffer
  logical, allocatable :: spec_pert(:), rxn_pert(:)
  integer :: nparams_total, nspec_params, nrxn_params
  !Keep track of which reactions are affected by each species. The first column
  !is the number of affected reactions. Subsequent columns are the integer
  !indices for the affected reactions.
  integer, allocatable :: spec_rxn_list(:,:)
  !The following arrays are for perturbing species energies in a thermodynamically
  !consistent manner. We can perturb a single species property and then propagate
  !that perturbation to other species by assuming that reaction properties are
  !conserved. We can also perturb reaction properties and back-propagate those
  !to the species. The species and reactions included in the perturbation are
  !controlled by the spec_pert and rxn_pert vectors. The rxn_pert_mat array is
  !used for calculating species property perturbations from the rxn property
  !perturbations. The spec_pert_mat is used for calculating direct perturbations
  !to species properties from a perturbation to a single species.
  double precision, allocatable :: rxn_pert_mat(:,:)
  integer, allocatable :: rxn_pert_mat_ID(:)
  double precision, allocatable :: spec_pert_mat(:,:)
  integer, allocatable :: spec_pert_mat_ID(:)
  logical :: pert_ref_spec, pert_rxn_prop
  integer :: ref_spec

  !Derived type for perturbations. This is intended to simplify the allocation
  !of memory and the storage of perturbation values.
  !Parameter definition for locations in the pert_val data member. We use
  !parameters to make the code more readable/less error prone.
  integer, parameter :: H_DFT_loc=1, H_GA_loc=2, H_LSR_loc=3
  integer, parameter :: S_DFT_loc=1, S_GA_loc=2
  integer, parameter :: Preexp_loc=1
  integer, parameter :: Beta_loc=1
  integer, parameter :: EA_DFT_loc=1, EA_BEP_loc=2
  integer, parameter :: BEP_slope_loc=1, BEP_int_loc=2
  integer, parameter :: LSR_slope_loc=1, LSR_int_loc=2
  !Type definition
  type :: SApert
    integer :: par_type !Type of parameter (e.g., enthalpy, entropy, etc.)
    integer :: par_number !Species/reaction/correlation number
    double precision :: pert_val(3) !The amount of perturbation from different sources
    double precision :: base_val(2) !The base/nominal value of the parameter
  end type
  !Array of perturbations
  type(SApert), allocatable :: SApert_data(:,:)

  !Vectors/arrays containing perturbation distributions
  double precision :: pert_local_dist(2,max_SA_param_types)
  double precision, allocatable :: Hspec_pert_dist(:,:)
  double precision, allocatable :: Sspec_pert_dist(:,:)
  double precision, allocatable :: Preexp_pert_dist(:)
  double precision, allocatable :: Beta_pert_dist(:)
  double precision, allocatable :: EA_pert_dist(:)
  double precision, allocatable :: EA_mvn_pert_dist(:,:,:)
  double precision, allocatable :: BEP_coeff_pert_dist(:,:)
  double precision, allocatable :: LSR_coeff_pert_dist(:,:)
  double precision, allocatable :: LSR_est_pert_dist(:)
  double precision, allocatable :: GA_Hest_pert_dist(:)
  double precision, allocatable :: GA_Sest_pert_dist(:)
  integer :: mvt_dof(max_SA_param_types)=0  !Degrees of freedom
  integer :: nmvn_EA_series !Number of EA mvn homologous series
  !The following array is where the mapping between reactions and homologous
  !series is kept for the mvn EA perturbations.
  integer, allocatable :: EA_mvn_pert_info(:,:)
  !The following array keeps track of the mapping between the MVT arrays
  !and the perturbed species.
  integer, allocatable :: Spec_pert_mvt_ID(:)

  !Output arrays
  !Sensitivity coefficients/indices
  double precision, allocatable :: NSC_matrix(:,:)  !Normalized SC for gas phase species only + conversion
  double precision, allocatable :: elem_effect_matrix(:,:,:)
  double precision, allocatable :: first_order_SI_matrix(:,:,:)
  double precision, allocatable :: total_SI_matrix(:,:,:)
  double precision, allocatable :: var_SI(:)
  double precision, allocatable :: w_var_SI(:,:,:)
  double precision, allocatable :: w_SA(:,:)  !Store solution vectors for each iteration
  double precision, allocatable :: rop_SA(:,:)  !Rates at each iteration

!*********************************************************************************************************
end module globalSA
!*********************************************************************************************************
