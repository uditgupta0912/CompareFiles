!*********************************************************************************************************
module input
!This module contains input related subroutines. The species and reactions included in the perturbation are
!controlled by the spec_pert and rxn_pert vectors.
!*********************************************************************************************************

  implicit none

contains

!*********************************************************************************************************
subroutine inco_tube
! This subroutine gets the program options, etc. contained in tube.inp. It is
! read once at the beginning of the program and never used again as the values
! read here are constant.
!*********************************************************************************************************

  use global_all; use global_tube1; use update_tube2; use file_handles; use globalSA

  implicit none

  integer     ::  i, ierr
  logical     ::  liso
  double precision  ::  press
  double precision :: logdt, logtout, omega_tmp
  character(len=80)   ::  rxn_str_tmp
  integer :: mem_stat, mem_stat_total=0
  integer :: nlines
  integer :: RxnID, NT, NVAL
  double precision :: RVAL(1)
  logical :: kerr
  character(len=80) :: MARI_str, conv_str
  integer :: SpecID
  integer :: maxl
!-------------------------------------------------------------------------------------------------------

!======================================================================================================
! Section 1 -- Digest tube.inp
!======================================================================================================
!.........Open input file:
  open(unit=itubeinp, file='INP.d/tube.inp', status='old', action='read')
!-------------------------------------------------------------------------------------------------------

  !Reactor type, number of runs, and whether each run has unique parameters
  call skip_comment(itubeinp)
  read(itubeinp,*) irxtr, nruns, MultiInput

  !Read operating conditions:
  if (MultiInput) then
    call skip_comment(itubeinp)
    read(itubeinp,*) lstp
  else
    call skip_comment(itubeinp)
    read(itubeinp,*) lstp,tin,press,velo,abyv,trise
  endif

  !More code options: gas/surface chemistry, isothermality, and reset reactor conditions
  !to be consistent with the chosen options.
  call skip_comment(itubeinp)
  read(itubeinp,*)liso,itpd
  if (liso) then
    iiso=0  !Zero here will make the non-isothermal term drop out of the energy balance
  else
    iiso=1
  end if
  if (itpd==0) then
    trpa_called=.true.
  else
    if (itpd==1) irxtr=0  !Must use a /UHV/ reactor for a UHV TPD
    ltra=.true.  !TPD experiments are inherently transient -- transient behavior is desired
  end if
  if (irxtr<=1) velo=0.0D0  !Batch and UHV reactors should have no flow in/out

  !Heat transfer
  call skip_comment(itubeinp)
  read(itubeinp,*)text,aextbyv,htc,ramp

  !Specify which species should be (1) printed during output and (2) is the primary reactant
  !Use name/phase/ construction to simplify entry.
  call skip_comment(itubeinp)
  read(itubeinp,*) MARI_str, conv_str

  !Convert species names into Species ID's
  call sksnum(MARI_str,0,LOUT,knams,kmax,PhaseNames,nphases,SpecTot,SpecID,NT,NVAL,RVAL(1),KERR)
  MARI_index=SpecID
  call sksnum(conv_str,0,LOUT,knams,kmax,PhaseNames,nphases,SpecTot,SpecID,NT,NVAL,RVAL(1),KERR)
  conv_index=SpecID

  if (MARI_index==0 .or. conv_index==0) then
    write(*,*) 'Invalid species specified for the MARI or the reactant. Stopping.'
    STOP
  end if

  !Options for controlling the number of reactors, the time steps, and how much
  !output is desired
  call skip_comment(itubeinp)
  read(itubeinp,*) rlen, nnodes, ttout, rtime, ntdec, ltra
  if (irxtr/=3) nnodes=1  !If not a PFR, only 1 node
  if (irxtr>1) then
    iflow=1 !CSTR and PFR are flow reactors
  else
    iflow=0 !Batch/UHV are not flow reactors
  end if

  !Solver parameters
  call skip_comment(itubeinp)
  read(itubeinp,*) ltol,abstol,reltol,NonNeg, idid_restart_count_max
  call skip_comment(itubeinp)
  read(itubeinp,*) iSolver, mu, ml

  !Read the option flags for semi-empirical techniques
  call skip_comment(itubeinp)
  read(itubeinp,*) lcov,lStatpQ,lBEP,iScale,lEA,lomega,omega_default, Tref_beta

  !Reaction path and sensitivity analysis
  call skip_comment(itubeinp)
  read(itubeinp,*)mrpa,verbose_rpa,trpa,lsenbrute,lDOE

  close(itubeinp)

!======================================================================================================
! Section 2 -- Digest omega.inp if needed (or just set to default)
!======================================================================================================

  if (lomega) then  !Omega correction on, use default plus values from file
    omega=omega_default
    open(unit=iomegainp, file='INP.d/omega.inp', status='unknown',action='read')  !proximity factor
    call skip_comment(iomegainp)
    read(iomegainp,*)nlines
    call skip_comment(iomegainp)
    do i=1,nlines
      read(iomegainp,*)omega_tmp,rxn_str_tmp
        call skcomp(rxn_str_tmp,rxn_str,nrxn,RxnID,NT)
        if (RxnID==0) then !Reaction not found
          write(*,*) 'Invalid reaction specified in omega.inp on line ', i, '. Stopping.'
          STOP
        end if
      omega(RxnID)=omega_tmp
    end do
    close(iomegainp)
  else !omega correction off
    omega=0.0D0
  end if

!======================================================================================================
! Section 3 -- Initialize some miscellaneous stuff
!======================================================================================================

  !Set up the time vector
  logdt=1/float(ntdec)
  logtout=log10(ttout)
  ndec=log10(rtime)-logtout  !Total number of decades of magnitude spanned
  ncalls=ceiling(ndec*ntdec) !Want to round UP to overshoot the stop time
  if (.not. allocated(tvec)) allocate(tvec(ncalls+1), stat=ierr)
  if (.not. allocated(tvec)) then
    write(*,*) 'Time vector not allocated. Stopping.'
    STOP
  end if

  tvec=0.
  tvec(ncalls+1)=rtime
  do i=2,ncalls
    tvec(i)=10**logtout
    logtout=logtout+logdt
  end do

  !Keep track of which reactions have user specied activation energies so
  !that if BEPs are used which might overwrite them an error is printed
  if (.not. allocated(lEA_user)) allocate(lEA_user(ngrxn+nsrxn))
  lEA_user=.false.

  !Allocate memory for solver arrays & set problem sizes
  liw=40+neqns
  select case (iSolver)
  case (1) !Krylov solver, DASPK
    maxl=max(5,neqns)
    lwp=(2*ml+mu+1)*neqns+2*((neqns/(ml+mu+1))+1) !Preconditioner real storage
    lrw=50+10*(neqns)+(neqns)**2+maxl+(maxl+3)*neqns+(maxl+3)*maxl+lwp
  case (0) !Direct solver, DASPK
    lrw=50+9*(neqns)+(neqns)**2
  case default  !Error condition
    write(*,*) 'Invalid solver selection specified. Using direct method, DASPK.'
    lrw=50+9*(neqns)+(neqns)**2
  end select
  if (.not. allocated(atol)) allocate(atol(neqns), stat=ierr)
  if (.not. allocated(atol)) mem_stat_total=mem_stat_total+1
  if (.not. allocated(rtol)) allocate(rtol(neqns), stat=ierr)
  if (.not. allocated(rtol)) mem_stat_total=mem_stat_total+1
  if (.not. allocated(iwork)) allocate(iwork(liw), stat=ierr)
  if (.not. allocated(iwork)) mem_stat_total=mem_stat_total+1
  if (.not. allocated(rwork)) allocate(rwork(lrw), stat=ierr)
  if (.not. allocated(rwork)) mem_stat_total=mem_stat_total+1
  if (mem_stat_total>0) then
    write(*,*) 'Errors initializing storage arrays for DASPK. Stopping.'
    STOP
  end if

  !Initialize solver options
  info=0  !Use defaults unless otherwise specified
  info(2)=1 !Specifies vector tolerances
  if (NonNeg) info(10)=2  !Specifies that Non-negativity constraint is ON during integration
  if (iSolver==1) then
    info(12)=1  !This turns the Krylov iterative solver /on/
    info(15)=1  !This specifies that there is a subroutine to calculate the preconditioner
    iwork(27)=lwp
    iwork(28)=neqns
  end if

  !Initialize the BEP routines
  if (lBEP) then
    call inco_BEP
    call zero_BEP_EA
  else
    if (.not. allocated(BEP_rxn)) allocate(BEP_rxn(nrxn,1))
    BEP_rxn=0 !This ensures that we can easily tell which rxns use BEPs and which don't
  end if

  !Initialize scaling relations
  if (iScale>0) then
    call inco_Scale
  else
    nBE_coords=1
  end if

  !Read the sensitivity analysis file
  if (lsenbrute) call inco_tube_SA

  !Initialize the coverage parameters
  if (lcov) then  !Coverage effects (constant values)
    call inco_cov
  else
    cov_matrix=0.0
  end if
  cov_factor=1.0  !This always starts at 1 and is adjusted as needed
  if (iScale<3) cov_factor_adjust=.false. !Only adjust interaction parameters for maps

  !Apply the StatpQ correction
  if (lStatpQ) then !Approximate corrections for T to constant BE
    call inco_StatpQ
  else
    StatpQ=0.0
  end if

  !If DOE or global SA perturbation is enabled, allocate memory and save the original correlations
  if (lDOE .or. lsenbrute) then
    if (perturb_param(globalSA_BEP)) then
      if (.not. allocated(BEP_def_save)) allocate(BEP_def_save(nBEP,6), stat=mem_stat)
      if (.not. allocated(BEP_def_save)) then
        write(*,*) 'Failure to allocate BEP_def_save. Stopping.'
        STOP
      end if
      BEP_def_save=BEP_def
    end if
    if (perturb_param(globalSA_LSR)) then
      if (.not. allocated(scale_slope_save)) allocate(scale_slope_save(nscale,2), stat=mem_stat)
      if (.not. allocated(scale_slope_save)) then
        write(*,*) 'Failure to allocate scale_slope_save. Stopping.'
        STOP
      end if
      scale_slope_save=scale_slope
    end if
  end if

!======================================================================================================
! Section 4 -- Write chosen options to screen and sort out possibly conflicting options
!======================================================================================================
!-------------------------------------------------------------------------------------------------------
  !Write operating conditions:
  select case (irxtr)
  case (0)
    write(*,*) 'UHV OR MOLECULAR BEAM REACTOR'
    write(*,*) '  VELOCITY RESET TO 0 AUTOMATICALLY'
    write(*,*) '  TOF CALCULATIONS IN tube_conv.out ARE DUMMY VALUES'
  case (1)
    write(*,*)'BATCH REACTOR'
    write(*,*) '  VELOCITY RESET TO 0 AUTOMATICALLY'
    write(*,*) '  TOF CALCULATIONS IN tube_conv.out ARE DUMMY VALUES'
    write(*,*) '  TRANSIENT RESULTS AUTOMATICALLY SAVED'
    ltra=.true.
  case (2)
    write(*,*)'CONTINUOUS STIRRED TANK REACTOR'
  case (3)
    write(*,*)'PLUG FLOW REACTOR'
  case default
    write(*,*)'INVALID REACTOR SELECTION. STOPPING'
    STOP
  end select
  if (.not. gas_chem) then
    write(*,*)'NO GAS CHEMISTRY'
  else
    write(*,*)'GAS CHEMISTRY PRESENT'
  endif
  if (.not. surf_chem) then
    write(*,*)'NO SURFACE CHEMISTRY'
  else
    write(*,*)'SURFACE CHEMISTRY PRESENT'
  endif
  if (itpd>0) then
    iiso=1  !TPD specification overrides isothermality flag
    liso=.false.
    write(*,*) 'TPD SPECIFICATION OVERRIDES ISOTHERMAL SPECIFICATION.'
  end if
  if (liso) then
    write(*,*)'ISOTHERMAL RUN'
  else
    write(*,*)'NON-ISOTHERMAL RUN'
    if (itpd>0) then
      write(*,*) 'TPD -- TRANSIENT RESULTS AUTOMATICALLY SAVED'
      ltra=.true.
      if (itpd==1) then
        write(*,*) 'UHV RUN, PRESSURE SET TO 1E-9 TORR. COMPOSITION SET TO 100% INERT.'
        write(*,*) 'UHV RUN, UHV REACTOR MODEL USED.'
      else
        write(*,*) 'HIGH PRESSURE RUN'
      end if
    else if (htc/=0.0) then
      write(*,*)'EXTERNAL HEAT TRANSFER'
    else
      write(*,*)'ADIABATIC RUN'
    endif
  endif
  if (lomega) write(*,*) 'OMEGA CORRECTION ON'
  if (lStatpQ) write(*,*) 'STATPQ T CORRECTION APPLIED TO BINDING ENERGIES'
  if (lBEP) write(*,*) 'BEP RELATIONS USED TO ESTIMATE BARRIERS'
  if (iScale>0) write(*,*) 'SCALING RELATIONS USED TO ESTIMATE HEATS OF FORMATION'
  if (iScale==3) write(*,*) 'SCALING RELATIONS USED TO PERFORM CATALYST SCREENING'
  if (lDOE) write(*,*) 'PARAMETER PERTURBATIONS IN EFFECT AFTER INITIAL RUN'
  if (lEA) write(*,*) 'EXTERNAL USER-SPECIFIED ACTIVATION ENERGIES'
  if (lstp) write(*,*)'INLET VELOCITY CHANGES WITH TEMPERATURE AND PRESSURE'
  if (.not. ltra) write(*,*)'!!NOTE!! TRANSIENT RESULTS NOT SAVED'
  select case (mrpa)
  case (0)
    write(*,*)'REACTION PATH ANALYSIS OFF'
    verbose_rpa=.false. !No RPA implies no verbose RPA
  case (1)
    write(*,*)'REACTION PATH ANALYSIS ON'
    trpa_called=.true. !No TPD RPA should be performed
  case (2)
    write(*,*) 'TEMPERATURE SPECIFIC TPD REACTION PATH ANALYSIS REQUESTED'
  case default
    mrpa=0
    verbose_rpa=.false. !No RPA implies no verbose RPA
    write(*,*) 'INVALID REACTION PATH ANALYSIS OPTION SELECTED'
    write(*,*) 'REACTION PATH ANALYSIS NOT WRITTEN'
  end select
  if (mrpa>0 .and. verbose_rpa .and. irxtr==3) then
    write(*,*) 'REACTION PATH ANALYSIS PERFORMED AT EVERY REACTOR NODE'
  else if (mrpa>0) then
    write(*,*) 'REACTION PATH ANALYSIS PERFORMED AT REACTOR EXIT'
    verbose_rpa=.false. !Batch/CST reactors only have RPA at the exit
  end if
  if (lDOE .and. lsenbrute) then
    write(*,*) 'DESIGN OF EXPERIMENTS AND SENSITIVITY ANALYSIS RUNS INCOMPATIBLE. STOPPING.'
    STOP
  end if
  if (.not. lsenbrute) isenbrute=0
  select case (isenbrute)
  case (0)
    write(*,*)'BRUTE FORCE SENSITIVITY ANALYSIS OFF'
  case (1)
    write(*,*) 'APPROXIMATE FIM-BASED SENSITIVITY ANALYSIS ON PRE-EXPONENTIALS'
  case (2)
    write(*,*) 'LOCAL BRUTE FORCE SENSITIVITY ANALYSIS SELECTED'
    if (rel_pert) then
      write(*,*) 'RELATIVE PERTURBATIONS TO PARAMETER VALUES USED'
    else
      write(*,*) 'ABSOLUTE PERTURBATIONS TO PARAMETER VALUES USED'
    end if
    if (count(spec_pert)>0) then
      if (fix_EA) then
        write(*,*) 'ACTIVATION ENERGIES FIXED WHEN SPECIES PROPERTIES PERTURBED'
      else
        write(*,*) 'TRANSITION STATE ENERGIES FIXED WHEN SPECIES PROPERTIES PERTURBED'
      end if
      end if
  case (3)
    write(*,*) 'VARIANCE-BASED GLOBAL SENSITIVITY ANALYSIS SELECTED'
    if (count(spec_pert)>0) then
      if (fix_EA) then
        write(*,*) 'ACTIVATION ENERGIES FIXED WHEN SPECIES PROPERTIES PERTURBED'
      else
        write(*,*) 'TRANSITION STATE ENERGIES FIXED WHEN SPECIES PROPERTIES PERTURBED'
      end if
    end if
  case (4)
    write(*,*) 'DERIVATIVE-BASED GLOBAL SENSITIVITY ANALYSIS SELECTED'
    if (count(spec_pert)>0) then
      if (fix_EA) then
        write(*,*) 'ACTIVATION ENERGIES FIXED WHEN SPECIES PROPERTIES PERTURBED'
      else
        write(*,*) 'TRANSITION STATE ENERGIES FIXED WHEN SPECIES PROPERTIES PERTURBED'
      end if
    end if
  case (5)
    write(*,*) 'UNCERTAINTY QUANTIFICATION ONLY'
    if (count(spec_pert)>0) then
      if (fix_EA) then
        write(*,*) 'ACTIVATION ENERGIES FIXED WHEN SPECIES PROPERTIES PERTURBED'
      else
        write(*,*) 'TRANSITION STATE ENERGIES FIXED WHEN SPECIES PROPERTIES PERTURBED'
      end if
    end if
  end select
  if (lsenbrute .and. isenbrute>2) then
    if (.not. model_solve) write(*,*) 'NO MODEL SOLUTION -- RANDOM SAMPLING OF DESIGN POINTS ONLY'
  end if
  if (lcov) then
    write(*,*) 'COVERAGE EFFECTS ON'
    if (cov_factor_adjust) then
      write(*,*) '  COVERAGE EFFECT PARAMETERS ADJUSTED FOR ATOMIC BINDING ENERGIES'
    else
      write(*,*) '  COVERAGE EFFECT PARAMETERS NOT ADJUSTED FOR ATOMIC BINDING ENERGIES'
    end if
  end if
  if (nruns > 1) then
    write(*,*)'MULTIPLE SETS OF REACTOR INPUT CONDITIONS'
    if (.not. MultiInput) then
      if (lDOE) write(*,*) 'DOE PERTURBATIONS BASED ON A SINGLE PARAMETER SET SPECIFIED.'
      if (trise/=0) write(*,*) 'TEMPERATURE RAMP SPECIFIED.'
      if (.not. lDOE .and. trise==0) then
        write(*,*) 'TEMPERATURE RISE OF 0 SPECIFIED FOR NON-DOE RUN. STOPPING.'
        STOP
      end if
    end if
  else
    write(*,*)'SINGLE SET OF REACTOR INPUT CONDITIONS'
  endif
  if (Tref_beta == 0) then
    write(*,*)'T_ref FOR BETA TERM OF KINETIC RATE CONSTANTS = 300 K'
  else
    write(*,*)'T_ref FOR BETA TERM OF KINETIC RATE CONSTANTS =   1 K'
  endif

  write(*,*)'************ MultiSite Variable Check ************'
  write(*,*)'number of site types:', nsphases
  write(*,*)'site type IDs:       ', site_type_lastID(2:nsphases+1)
!-------------------------------------------------------------------------------------------------------
  if (itpd==1) press=1.32E-12  !1E-9 torr (typical UHV pressure)
  p=press*patm      ! atm --> dyne/cm2
!-------------------------------------------------------------------------------------------------------

!*********************************************************************************************************
end subroutine inco_tube
!*********************************************************************************************************

!*********************************************************************************************************
subroutine get_tube_conditions
!This subroutine is a wrapper which calls all the other subroutines to get
!run-specific parameters and store them in global arrays.
!*********************************************************************************************************

  implicit none

  call inco_Tflow  !Operating conditions (T, P, etc.)
  call inco_tube_mole  !Feed composition
  call inco_EA !User-specified activation energies
  call inco_ddaspk !Solver parameters, tolerances, etc.

!*********************************************************************************************************
end subroutine get_tube_conditions
!*********************************************************************************************************

!*********************************************************************************************************
subroutine update_tube_conditions
!This subroutine extracts the run-specific parameters from the global arrays.
!*********************************************************************************************************

  use global_all; use global_tube1; use update_tube2; use file_handles

  implicit none

  integer :: i, j, nrows
  double precision  ::  addo, t
  double precision :: EA_tmp

  !Feed condition
  tin=TPFSV(i_tube,1)
  p=TPFSV(i_tube,2)
  velo=TPFSV(i_tube,3)
  abyv=TPFSV(i_tube,4)
  t=tin

  !Check to see if we need to adjust for STP
  if (lstp) velo=velo*(t/273.15)*(patm/p)  ! Account for T,P effect on velocity

  !Calculate the density and the mass flux
  call ckrhoy(p,t,yin,ickwrk,rckwrk,rhoin)   ! Inlet density
  fluxmass=rhoin*velo         ! Mass flux

  !Feed composition
  !The rest of this is useful for converting the act_all values into the actual compositions...
  actin(1:kgmax+ksmax)=0.0D0  !Initialize to zero and set only non-zero entries

  !Extract the input compositiion
  nrows=size(act_all,1)
  if (itube_restart/=0) nrows=nrows-1
  do i=1,nrows
    actin(iMole_spec(i))=act_all(i,i_tube)
  end do

  !Correct for UHV TPD
  if (itpd==1) then
    forall (i=1:kgmax-1) actin(i)=0.0D0
    actin(kgmax)=1.0D0
  end if

  !Normalize compositions for all phases to have activities summing to 1
  !Activities are mole fractions for gas & coverages for each type of surface site
  do_999: do i=1,nsphases+1  ! Loop over all surface site types
      addo=sum(actin(site_type_firstID(i):site_type_lastID(i))) ! sums the initial values of the (i-1)th site
      if (abs(1-addo)>1E-4) then
        write(*,*) 'Warning: possible error in composition in phase ', &
          PhaseNames(i),' due to renormalization. Continuing.'
      end if
      actin(site_type_firstID(i):site_type_lastID(i))=actin(site_type_firstID(i):site_type_lastID(i))/addo ! normalize coverages for (i+1)th site
  end do do_999
  actin(kmax-kbmax+1:kmax)=1.0D0 !Bulk species have activity of 1

  act=actin
  x(1:kgmax)=act(1:kgmax)
  xin=x
  call ckxty(x,ickwrk,rckwrk,y)
  call ckxty(xin,ickwrk,rckwrk,yin)
  call ckrhoy(p,t,yin,ickwrk,rckwrk,rhoin) ! Inlet density
  win(1:kgmax)=yin
  forall (i=kgmax+1:kmax-kbmax) win(i)=actin(i)
  win(indext)=t !Temperature
  win(indexrho)=rhoin !Density

  fluxmass=rhoin*velo       ! Mass flux
  if (itube_restart/=0) fluxmass_inlet=act_all(nrows,i_tube)

  !Activation energies
  if (lEA) then
    if (gas_chem) then
      nrows=size(iEAg_rxn,1)
      do i=1,nrows
        j=iEAg_rxn(i)
        EA_tmp=EAg_all(i,i_tube)*tin !Use Tin, not T since the EA values are divided by RTin
        call ckreex(-j,rckwrk,EA_tmp)
      end do
    end if
    if (surf_chem) then
      nrows=size(iEAs_rxn,1)
      do i=1,nrows
        j=iEAs_rxn(i)
        EA_tmp=EAs_all(i,i_tube)*tin
        call skreex(-j,iskwrk,rskwrk,EA_tmp)
      end do
    end if
  end if

  !Solver tolerances
  if (ltol) then
    nrows=size(iTol_spec,1)
    do i=1,nrows
      atol(iTol_spec(i))=atol_all(i,i_tube)
      rtol(iTol_spec(i))=rtol_all(i,i_tube)
    end do
  end if

!*********************************************************************************************************
end subroutine update_tube_conditions
!*********************************************************************************************************

!*********************************************************************************************************
subroutine inco_Tflow
!This subroutine is responsible for reading conditions from T_flow.
!*********************************************************************************************************

  use global_all; use file_handles; use global_tube1; use update_tube2

  implicit none

  integer :: i, ierr
  double precision  ::  press
  character(len=80) :: crap

  if (.not. allocated(TPFSV)) allocate(TPFSV(nruns,4), stat=ierr)
  if (.not. allocated(TPFSV)) then
    write(*,*) 'Feed condition array not allocated. Stopping.'
    STOP
  end if

  if (.not. MultiInput) then !Assign the tube.inp values to the TPFSV array
    do i=1,nruns
      TPFSV(i,1)=tin+(i-1)*trise
    end do
    TPFSV(:,2)=p  !Already in dyne/cm^2 (converted in inco_tube)
    TPFSV(:,3)=velo
    TPFSV(:,4)=abyv
  else  !Read from T_flow.inp
    open(unit=iTflowinp, file='INP.d/T_flow.inp', status='unknown', action='read')
    call skip_comment(iTflowinp)
    do i=1,nruns
      read(iTflowinp,*) TPFSV(i,:) !tin,press,velo,abyv -- Current operating conditions
    end do
    TPFSV(:,2)=TPFSV(:,2)*patm      ! atm --> dyne/cm2
    close(iTflowinp)
  end if

!*********************************************************************************************************
end subroutine inco_Tflow
!*********************************************************************************************************

!*********************************************************************************************************
subroutine inco_tube_mole
!This subroutine is responsible for reading conditions from T_flow.
!*********************************************************************************************************

  use global_all; use file_handles; use global_tube1; use update_tube2

  implicit none

  integer :: i, ii, nlines, ierr
  character(len=80) :: crap
  double precision :: act_local

  !Read inlet mole fractions and coverages (ignore bulk species)
  open(unit=itubemoleinp, file='INP.d/tube_mole.inp', status='unknown', action='read')
  call skip_comment(itubemoleinp)
  read(itubemoleinp,*)itube_restart
  call skip_comment(itubemoleinp)
  read(itubemoleinp,*)nlines
  call skip_comment(itubemoleinp)

  !Check for the restart and adjust the number of lines read
  if (itube_restart>0) nlines=nlines-1

  !Initialize the vector mapping species ID numbers to input file lines
  if (.not. allocated(iMole_spec)) allocate(iMole_spec(nlines), stat=ierr)
  if (.not. allocated(iMole_spec)) then
    write(*,*) 'Species ID vector not allocated. Stopping.'
    STOP
  end if
  call spec_vec_init(nlines,itubemoleinp,iMole_spec)

  !Now allocate the memory for the stored feed composition array
  if (itube_restart>0) nlines=nlines+1 !Account for mass flux for restart
  if (.not. allocated(act_all)) allocate(act_all(nlines,nruns), stat=ierr)
  if (.not. allocated(act_all)) then
    write(*,*) 'Feed composition array not allocated. Stopping.'
    STOP
  end if

  !Note: if itube_restart>0, then the last row in act_all is actually the mass flux!
  do ii=1,nlines
    if (MultiInput) then
      read(itubemoleinp,*)crap, act_all(ii,:)
    else
      read(itubemoleinp,*)crap, act_local
      act_all(ii,:)=act_local
    end if
  end do
  close(itubemoleinp)

!*********************************************************************************************************
end subroutine inco_tube_mole
!*********************************************************************************************************

!*********************************************************************************************************
subroutine inco_EA
!This subroutine gets user-specified activation energies from EA?.inp
!*********************************************************************************************************

  use global_all; use global_tube1; use update_tube2; use file_handles

  implicit none

  integer :: irxn, jrxn, rxn_count, i
  character(len=80) :: crap

  if (.not. lEA) return

  !Read in user specified activation energies for gas reactions and reset
  !the corresponding value in the work arrays.
  if (gas_chem) then
    open(unit=iEAginp, file='INP.d/EAg.inp', status='unknown', action='read')
    call skip_comment(iEAginp)
    read(iEAginp,*) rxn_count
    call skip_comment(iEAginp)

    !Initialize the vector mapping reaction ID numbers to input file lines
    if (.not. allocated(iEAg_rxn)) allocate(iEAg_rxn(rxn_count))
    call rxn_vec_init(rxn_count,iEAginp,iEAg_rxn)

    !Initialize the storage arrays for the EAg values
    if (.not. allocated(EAg_all)) allocate(EAg_all(rxn_count,nruns))
    do irxn=1,rxn_count
      if (MultiInput) then
        read(iEAginp,*) crap, EAg_all(irxn,:)
      else
        read(iEAginp,*) crap, EAg_all(irxn,1)
        EAg_all(irxn,2:nruns)=EAg_all(irxn,1)
      end if
      jrxn=iEAg_rxn(irxn)
      lEA_user(jrxn)=.true.
    end do
    close(iEAginp)
  end if

  !Read in user specified activation energies for surface reactions and reset
  !the corresponding value in the work arrays.
  if (surf_chem) then
    open(unit=iEAsinp, file='INP.d/EAs.inp', status='unknown', action='read')
    call skip_comment(iEAsinp)
    read(iEAsinp,*) rxn_count
    call skip_comment(iEAsinp)

    !Initialize the vector mapping reaction ID numbers to input file lines
    if (.not. allocated(iEAs_rxn)) allocate(iEAs_rxn(rxn_count))
    call rxn_vec_init(rxn_count,iEAsinp,iEAs_rxn)

    !Initialize the storage arrays for the EAs values
    if (.not. allocated(EAs_all)) allocate(EAs_all(rxn_count,nruns))
    !Accounts for the fact that the rxn_str vector has both gas and surface
    !reaction ID numbers but the SK work arrays have only surface reactions.
    iEAs_rxn=iEAs_rxn-ngrxn

    do irxn=1,rxn_count
      if (MultiInput) then
        read(iEAsinp,*) crap, EAs_all(irxn,:)
      else
        read(iEAsinp,*) crap, EAs_all(irxn,1)
        EAs_all(irxn,2:nruns)=EAs_all(irxn,1)
      end if
      jrxn=iEAs_rxn(irxn)
      lEA_user(ngrxn+jrxn)=.true.
    end do
    close(iEAsinp)
  end if

!*********************************************************************************************************
end subroutine inco_EA
!*********************************************************************************************************

!*********************************************************************************************************
subroutine inco_ddaspk
!This subroutine gets initial parameters for ddaspk as well as allocates memory for arrays necessary
!to solving the problem.
!*********************************************************************************************************

  use global_all; use global_tube1; use update_tube2; use file_handles

  implicit none

  double precision :: tol_tmp(2*nruns)
  integer :: i, j, ierr, nlines, mem_stat_total=0
  character(len=80)   ::  crap
  logical :: lenergy

!-------------------------------------------------------------------------------------------------------

  !Set solver tolerances
  atol=abstol; rtol=reltol

  if (.not. ltol) return  !Only read the rest of the input for granular tolerances

  open(unit=itolinp, file='INP.d/tol.inp', status='old', action='read')
  call skip_comment(itolinp)
  read(itolinp,*) nlines

  if (nlines<=0) then !There's nothing else to do here
    close(itolinp)
    return
  end if

  call skip_comment(itolinp)
  read(itolinp,*) lenergy
  call skip_comment(itolinp)

  !Initialize storage and species ID arrays if needed
  if (.not. allocated(iTol_spec)) allocate(iTol_spec(nlines), stat=ierr)
  if (.not. allocated(iTol_spec)) then
    write(*,*) 'Species ID vector not allocated. Stopping.'
    STOP
  end if
  if (.not. allocated(atol_all)) allocate(atol_all(nlines,nruns), stat=ierr)
  if (.not. allocated(atol_all)) then
    write(*,*) 'Absolute tolerance storage vector not allocated. Stopping.'
    STOP
  end if
  if (.not. allocated(rtol_all)) allocate(rtol_all(nlines,nruns), stat=ierr)
  if (.not. allocated(rtol_all)) then
    write(*,*) 'Absolute tolerance storage vector not allocated. Stopping.'
    STOP
  end if
  if (lenergy) then !Last line is energy balance
    call spec_vec_init(nlines-1,itolinp,iTol_spec)
    iTol_spec(nlines)=neqns
  else
    call spec_vec_init(nlines,itolinp,iTol_spec)
  end if

  !Read in tolerances
  do j=1,nlines !Assume that the /last/ line is the energy balance
    if (MultiInput) then
      read(itolinp,*)crap,tol_tmp(:)
      atol_all(j,:)=tol_tmp(1:2*nruns-1:2) !Odd elements
      rtol_all(j,:)=tol_tmp(2:2*nruns:2) !Even elements
    else
      read(itolinp,*)crap,tol_tmp(1:2)
      atol_all(j,:)=tol_tmp(1)
      rtol_all(j,:)=tol_tmp(2)
    end if
  end do
  close(itolinp)

!*********************************************************************************************************
end subroutine inco_ddaspk
!*********************************************************************************************************

!*********************************************************************************************************
subroutine inco_BEP
!This subroutine reads BEP.inp and stores the necessary variables to global before being used in skBEP
!*********************************************************************************************************
  use global_all; use global_tube1; use update_tube2; use file_handles; use globalSA

  implicit none

  integer :: del, i, mem_stat, mem_stat_total=0, ierr, nlines, j
  integer :: RxnID, NT, BEP_rxn_tmp(nrxn,3)
  character(len=80) :: crap, rxn_str_tmp

!-------------------------------------------------------------------------------------------------------
  write(*,*)'************************************************'
!.........Open input files:
  open(unit=iBEPinp, file='INP.d/BEP.inp', status='old', action='read')
!-------------------------------------------------------------------------------------------------------
  call skip_comment(iBEPinp)
  read(iBEPinp,*) nBEP                                              !read number of correlations

  !Allocate the necessary storage space
  if (.not. allocated(BEP_def)) allocate(BEP_def(nBEP,6), stat=mem_stat)
  if (.not. allocated(BEP_def)) mem_stat_total=mem_stat_total+1
  if (.not. allocated(BEP_rxn)) allocate(BEP_rxn(nrxn,2), stat=mem_stat)
  if (.not. allocated(BEP_rxn)) mem_stat_total=mem_stat_total+1
  if (mem_stat_total>0) then
    write(*,*) 'Memory allocation failure in inco_BEP, stopping'
    STOP
  end if
  BEP_def=0.0D0
  BEP_rxn=0

  call skip_comment(iBEPinp)
  do i=1,nBEP
    !read in correlation info (type, m, b, direction)
    read(iBEPinp,*) del, (BEP_def(i,j), j=1,6)
    if ((BEP_def(i,1)<-1) .OR. BEP_def(i,1)>1) then
      write(*,*) 'Invalid correlation type in BEP ',i,'. Stopping'
      STOP
    end if
    if (int(abs(BEP_def(i,4)))/=1) then
      write(*,*) 'Invalid correlation direction in BEP',i,'. Stopping'
      STOP
    end if
  end do
  !Normalize by R so that energies are stored internally in K
  BEP_def(:,3)=BEP_def(:,3)/Rgas_kcal
  BEP_def(:,5:6)=BEP_def(:,5:6)/Rgas_kcal

  call skip_comment(iBEPinp)
  read(iBEPinp,*) nlines
  call skip_comment(iBEPinp)
  call BEP_rxn_read(nlines,nBEP,iBEPinp,BEP_rxn_tmp(1:nlines,:))

  !Assign the data
  do i=1,nlines
    BEP_rxn(BEP_rxn_tmp(i,3),:)=BEP_rxn_tmp(i,1:2)
  end do
  close(iBEPinp)

!*********************************************************************************************************
end subroutine inco_BEP
!*********************************************************************************************************

!*********************************************************************************************************
subroutine zero_BEP_EA
!This subroutine reinitializes the activation energies for reactions
!using BEP correlations to 0.
!*********************************************************************************************************

  use global_all;

  implicit none

  integer :: i
  double precision :: zero

  zero=0.0
  if (nBEP==0) return

  !Go through the reactions resetting original activation energies of reactions
  !using BEP correlations to zero
  do i=1,ngrxn
    if (BEP_rxn(i,1)/=0) then  !This reaction does not use a BEP
      call ckreex(-i,rckwrk,zero)
      reg(i)=0.0
    end if
  end do
  do i=1,nsrxn
    if (BEP_rxn(ngrxn+i,1)/=0) then  !This reaction does not use a BEP
      call skreex(-i,iskwrk,rskwrk,zero)
      re(i)=0.0
    end if
  end do

!*********************************************************************************************************
end subroutine zero_BEP_EA
!*********************************************************************************************************

!*********************************************************************************************************
subroutine inco_Scale
!This subroutine reads Scale.inp and stores the necessary variables to global before being used in skScale
!*********************************************************************************************************
  use global_all; use global_tube1; use update_tube2; use file_handles; use globalSA

  implicit none

  double precision, allocatable :: BE_range(:,:)
  logical, allocatable :: single_point(:)
  integer, allocatable :: spec(:)
  integer :: i, grid_type, nspec, skip
  character(len=80) :: spec_name
  double precision :: BE_delta
  !This maps the species IDs to entries in the scale_ref array
  integer :: specID(kmax)

!-------------------------------------------------------------------------------------------------------

  open(unit=iScaleinp, file='INP.d/Scale.inp',status='old', action='read')

  !Read the number of BE descriptors, their ranges, and references
  call skip_comment(iScaleinp)
  read(iScaleinp,*) natoms_scale
  !Allocate memory
  if (.not. allocated(scale_ref)) allocate(scale_ref(ksmax))
  if (.not. allocated(scale_targ)) allocate(scale_targ(natoms_scale,2))
  if (.not. allocated(BE_range)) allocate(BE_range(2,natoms_scale))
  if (.not. allocated(single_point)) allocate(single_point(natoms_scale))
  if (.not. allocated(spec)) allocate(spec(natoms_scale))
  call spec_vec_init(natoms_scale,iScaleinp,spec)  !This gets the atom names
  scale_targ(:,1)=dble(spec)  !Save the species ID numbers
  scale_ref=1.0

  !This gets the data
  specID=0
  do i=1,natoms_scale
    specID(spec(i))=i
    read(iScaleinp,*) spec_name, BE_range(:,i), scale_ref(spec(i)-kgmax), &
      single_point(i)
    !The following ensures that a single point calculation for an atomic
    !binding energy will result in a range restricted to the lower value.
    if (single_point(i) .or. iScale<3) BE_range(2,i)=BE_range(1,i)
  end do
  deallocate(spec)  !Deallocate so we can reallocate later

  !Read the grid options and setup the grid
  call skip_comment(iScaleinp)
  if (all(single_point) .or. iScale<3) then !Single point only
    nBE_coords=1
    if (.not. allocated(BE_coords)) allocate(BE_coords(natoms_scale,nBE_coords))
    BE_coords(:,1)=BE_range(1,:)
    read(iScaleinp,*) grid_type !Skip grid type
    if (grid_type/=4) then
      read(iScaleinp,*) !Skip grid parameters
    else
      read(iScaleinp,*) nBE_coords
      do i=1,nBE_coords
        read(iScaleinp,*) !Skip specified coordinates
      end do
      nBE_coords=1
    end if
  else  !Binding energy map
    read(iScaleinp,*) grid_type
    if (grid_type==0) then  !Default
      if (natoms_scale>3) then  !More than 3 dimensions -- use Sobol'
        grid_type=3
      else if (natoms_scale==1) then
        grid_type=1 !Use rectangular grid for 1-D case
      else  !Use hex-grid for 2-D/3-D cases
        grid_type=2
      end if
    else if (grid_type==2) then
      if (natoms_scale>3) then
        write(*,*) 'Fatal error: Hex grid not implemented for more than three dimensions. Stopping.'
        STOP
      end if
      if (natoms_scale==1) grid_type=1 !Switch to rectangular grid
    end if
    select case (grid_type)
    case (1, 2)  !Rectangular, Hex
      read(iScaleinp,*) BE_delta  !BE increment
      if (grid_type==1) then
        call make_rect_grid(BE_delta,BE_range)
      else
        call make_hex_grid(BE_delta,BE_range)
      end if
    case (3)  !Sobol'
      read(iScaleinp,*) nBE_coords, skip  !Number of binding energy coordinates
      call make_sobol_grid(BE_range,skip)
    case (4)  !User-specified
      read(iScaleinp,*) nBE_coords
      if (.not. allocated(BE_coords)) allocate(BE_coords(natoms_scale,nBE_coords))
      do i=1,nBE_coords
        read(iScaleinp,*) BE_coords(:,i)
      end do
    end select
  end if
  BE_coords=BE_coords/Rgas_kcal !Convert to dimensionless units

  !Read the scaling relation definitions and the corresponding atoms
  call skip_comment(iScaleinp)
  read (iScaleinp,*) nscale !Total number of scaling relations
  if (.not. allocated(spec)) allocate(spec(nscale))
  call spec_vec_init(nscale,iScaleinp,spec)  !This gets the atom names
  if (.not. allocated(scale_slope)) allocate(scale_slope(nscale,2))
  do i=1,nscale
    scale_slope(i,1)=dble(specID(spec(i)))
    read(iScaleinp,*) spec_name, scale_slope(i,2)
  end do
  deallocate(spec)

  !Read the binding modes and the reference binding energies
  call skip_comment(iScaleinp)
  read (iScaleinp,*) nspec !Total number of species
  if (.not. allocated(spec)) allocate(spec(nspec))
  call spec_vec_init(nspec,iScaleinp,spec)  !This gets the atom names
  if (.not. allocated(scale_bind_mode)) allocate(scale_bind_mode(nscale+1,ksmax))
  scale_bind_mode=0
  do i=1,nspec
    scale_bind_mode(1,spec(i)-kgmax)=1
    read(iScaleinp,*) spec_name, scale_bind_mode(2:nscale+1,spec(i)-kgmax), &
      scale_ref(spec(i)-kgmax)
  end do
  deallocate(spec)
  scale_ref=scale_ref/Rgas_kcal !Normalize to be in units of K

  !Read the uncertainty distribution
  call skip_comment(iScaleinp)
  if (.not. allocated(LSR_est_pert_dist)) allocate(LSR_est_pert_dist(2))
  read(iScaleinp,*) LSR_est_pert_dist
  LSR_est_pert_dist=LSR_est_pert_dist/Rgas_kcal

  close(iScaleinp)

  !Deallocate temporary memory
  deallocate(BE_range)
  deallocate(single_point)

!*********************************************************************************************************
end subroutine inco_Scale
!*********************************************************************************************************

!*********************************************************************************************************
subroutine inco_DOE(DOE_data)
! This subroutine reads DOE.inp and stores the necessary data to arrays.
! It also tests the independence of the gas phase basis set and stops the code if an error occurs.
!*********************************************************************************************************
  use global_all; use global_tube1; use update_tube2; use file_handles

  implicit none

  integer :: i, j
  integer :: num_basis, num_param, ierr, kktot, dummy, inx, job=11, trans_mat_ind = 1
  double precision :: det(2), rcond
  type(runs), intent(out) :: DOE_data(nruns-1)
  character(len=80) :: gas_basis_str(nelm), basis_elem_str(nelm)
  integer :: basis_spec(nelm), basis_elem(nelm), pert_spec(ksmax)
  integer :: nbasis_spec, nbasis_elem, npert_spec
  character(len=80) :: name_str
  double precision :: scratch_val(2)
  integer :: SpecID, NT, NVAL, RxnID
  double precision :: RVAL
  logical :: kerr

!-------------------------------------------------------------------------------------------------------

  !Initialize species info
  basis_spec=0
  basis_elem=0
  pert_spec=0

  ! Read in and store data

  open(unit=iDOEinp, file='INP.d/DOE.inp',status='old', action='read')
  call skip_comment(iDOEinp)
  read(iDOEinp,*) num_basis ! Number of basis species

  if (num_basis>nelm) then
    write(*,*) 'More basis species specified than elements. Aborting.'
    STOP
  end if

  call skip_comment(iDOEinp)
  read(iDOEinp,*) (gas_basis_str(i), i=1, num_basis) ! Read in basis species names

  !Need to populate basis_spec array with species ID numbers
  do i=1,num_basis
    name_str=trim(adjustl(gas_basis_str(i)))//'/GAS/' !Assume that these are gas-phase species
    call sksnum(name_str,0,LOUT,knams,kmax,PhaseNames,nphases,SpecTot,SpecID,NT,NVAL,RVAL,KERR)
    if (SpecID/=0) then
      basis_spec(i)=SpecID
    else
      name_str=trim(adjustl(gas_basis_str(i)))
      write(*,'(A,A16,A)') 'Invalid gas-phase basis species ', name_str,' supplied. Stopping.'
      STOP
    end if
  end do
  nbasis_spec=count(basis_spec/=0)

  !Allocate and initialize thermo_basis
  if (.not. allocated(thermo_basis)) allocate(thermo_basis(num_basis))
  thermo_basis=basis_spec(1:num_basis)

  call skip_comment(iDOEinp)

  do i = 1, nruns-1    ! Loop over number of runs-1
    call skip_comment(iDOEinp)
    read(iDOEinp,*) num_param
    if (allocated(DOE_data(i)%experiment)) deallocate(DOE_data(i)%experiment)
    allocate(DOE_data(i)%experiment(num_param))
    DOE_data(i)%num_param = num_param
    do j = 1,num_param               ! Loop over number of parameters for the current run
      read(iDOEinp,'(A)') name_str  !Get the species/reaction string (this will be processed later)
      name_str=trim(adjustl(name_str))  !Get rid of extra spaces
      read(iDOEinp,*) dummy              ! Read the first entry of the line to a temp var.
      backspace(iDOEinp)                 ! Return to the beginning of the line
      select case (dummy)
      case (1:2)  !This is a species -- Enthalpy (1), Entropy (2)
        call sksnum(name_str,0,LOUT,knams,kmax,PhaseNames,nphases,SpecTot,SpecID,NT,NVAL,RVAL,KERR)
        if (SpecID==0) then
          write(*,'(A,A33,A)') 'Invalid species ', name_str, ' specified in DOE.inp on run ', i, '. Stopping.'
          STOP
        end if
        if (any(SpecID==site_type_lastID(2:nsphases+1))) then
          write(*,'(A,A33,A)') 'Vacant site type ', name_str, ' specified in DOE.inp on run ', i, '. Stopping.'
          STOP
        end if
        DOE_Data(i)%experiment(j)%id=SpecID
        read(iDOEinp,*) DOE_Data(i)%experiment(j)%param, DOE_Data(i)%experiment(j)%value(1)
        if (dummy==1) DOE_Data(i)%experiment(j)%value(1)=DOE_Data(i)%experiment(j)%value(1)/Rgas_kcal
        if (dummy==2) DOE_Data(i)%experiment(j)%value(1)=DOE_Data(i)%experiment(j)%value(1)/Rgas_cal
      case (3:5)  !This is a reaction -- Pre-exponential (3), Beta (4), Activation energy (5)
        call skcomp(name_str,rxn_str,nrxn,RxnID,NT)
        if (RxnID==0) then !Reaction not found
          write(*,'(A,A80,/,A,I3,A)') 'Invalid reaction ', name_str, ' specified in DOE.inp on run ', i, '. Stopping.'
          STOP
        end if
        DOE_Data(i)%experiment(j)%id=RxnID
        read(iDOEinp,*) DOE_Data(i)%experiment(j)%param, DOE_Data(i)%experiment(j)%value(1)
        if (dummy==5) DOE_Data(i)%experiment(j)%value(1)=DOE_Data(i)%experiment(j)%value(1)/Rgas_kcal
      case (6)  !Linear scaling relation
        if (iScale==0) then !LSRs not turned on
          write(*,*) 'Error in DOE.inp, run ', i, '. Linear scaling relations not turned on in tube.inp.'
          write(*,*) 'Please change iScale in tube.inp or correct the perturbation type.'
          write(*,*) 'Stopping.'
          STOP
        end if
        read(name_str,*,iostat=ierr) SpecID !Internal read to convert the string to an integer
        if (ierr/=0) then !A conversion error occured
          write(*,*) 'Error in DOE.inp on run ', i, '.'
          write(*,*) 'Conversion of the correlation number ', name_str, ' to an integer failed.'
          write(*,*) 'Stopping.'
          STOP
        end if
        if (SpecID<1 .or. SpecID>nscale) then
          write(*,*) 'An invalid correlation number ', SpecID, ' was specified in DOE.inp on run ', i, '.'
          write(*,*) 'Please specify an integer between 1 and ', nscale, '.'
          write(*,*) 'Stopping.'
          STOP
        end if
        DOE_Data(i)%experiment(j)%id=SpecID
        read(iDOEinp,*) DOE_Data(i)%experiment(j)%param, DOE_Data(i)%experiment(j)%value(1)
      case (7)  !Free energy correlation
        if (.not. lBEP) then
          write(*,*) 'Error in DOE.inp, run ', i, '. BEP correlations not turned on in tube.inp.'
          write(*,*) 'Please change lBEP in tube.inp or correct the perturbation type.'
          write(*,*) 'Stopping.'
          STOP
        end if
        read(name_str,*,iostat=ierr) RxnID !Internal read to convert the string to an integer
        if (ierr/=0) then !A conversion error occured
          write(*,*) 'Error in DOE.inp on run ', i, '.'
          write(*,*) 'Conversion of the correlation number ', name_str, ' to an integer failed.'
          write(*,*) 'Stopping.'
          STOP
        end if
        if (RxnID<1 .or. RxnID>nBEP) then
          write(*,*) 'An invalid correlation number ', RxnID, ' was specified in DOE.inp on run ', i, '.'
          write(*,*) 'Please specify an integer between 1 and ', nBEP, '.'
          write(*,*) 'Stopping.'
          STOP
        end if
        DOE_Data(i)%experiment(j)%id=RxnID
        read(iDOEinp,*) DOE_Data(i)%experiment(j)%param, DOE_Data(i)%experiment(j)%value(1), &
          DOE_Data(i)%experiment(j)%value(2)
        DOE_Data(i)%experiment(j)%value(2)=DOE_Data(i)%experiment(j)%value(2)/Rgas_kcal
      case default  !Error case
        write(*,*) 'Invalid perturbation option of ', dummy, ' chosen on run ', i, '. Stopping.'
        STOP
      end select
    end do
  end do

  close(iDOEinp) ! End of reading in DOE input
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Check that the gas-phase species form an independent basis; if so, store data

  ! Create a similar matrix containing surface species molecular info (exclude vacancies, bulk)
  inx=1
  do i=1,ksmax
    if (any(site_type_lastID==i+kgmax)) cycle
    pert_spec(inx)=i+kgmax
    inx=inx+1
  end do
  npert_spec=count(pert_spec/=0)

  !Search through the composition matrix finding the elements which are used
  !either in the basis species or the perturbed species.
  inx=1
  do i=1,nelm
    if (any(elem_comp(i,basis_spec(1:nbasis_spec))/=0) .or. &
      any(elem_comp(i,pert_spec(1:npert_spec))/=0)) then
      !Add to the list of basis elements
      basis_elem(inx)=i
      inx=inx+1
    end if
  end do
  nbasis_elem=count(basis_elem/=0)

  ! Determine whether matrix is nonsingular (and therefore invertible)
  if (nbasis_elem/=num_basis) then
    write(*,*)'Error - basis matrix is not square - choose different basis set.'
    stop
  end if

  ! Compute transformation matrix
  allocate(trans_mat(npert_spec,num_basis),stat=ierr)
  call find_thermo_basis(nbasis_spec,basis_spec(1:nbasis_spec),nbasis_elem,&
    basis_elem(1:nbasis_elem),npert_spec,pert_spec(1:npert_spec),trans_mat)

  ! Write out transformation matrix to file
  open(unit=itranmatout,file='OUT.d/trans_matrix.out',status='replace',action='write')
  write(itranmatout,'(17X,10(1X,A16))') (knams(basis_spec(i)),i=1,num_basis)
  trans_mat_ind = 1
  do i=1,ksmax
    if (any(site_type_lastID==i+kgmax)) cycle
    write(itranmatout,'(1X,A16,10(1X,F16.4))')knams(kgmax+i),trans_mat(trans_mat_ind,:)
    trans_mat_ind = trans_mat_ind + 1
  end do
  close(itranmatout)

!*********************************************************************************************************
end subroutine inco_DOE
!*********************************************************************************************************

!*********************************************************************************************************
subroutine inco_cov
!This subroutine reads in coverage effects in form of matrix of effect on binding energy.
!*********************************************************************************************************

  use global_all; use global_tube1; use update_tube2; use file_handles

  implicit none

  double precision               ::      RVAL(5),Nlines1
  integer                        ::      NCov,NLines,ca,cb,NEXP,spec1TOT,spec2TOT,NT,NVAL,i
  integer                        ::      Nspec,Nspec2
  character(len=80)              ::      LINE
  character(len=80)              ::      species
  logical                        ::      KERR
  integer :: cov_model
!---------------------------------------------------------------------------------------------
  !cov_matrix is a matrix which lists each surface species' dependence on other surface species
  cov_matrix=0.0D0
  cov_matrix2=0.0D0
  cov_matrix3=0.0D0
  !thresh_cov is a matrix that gives the threshold below which coverage effects
  !are ignored
  thresh_cov=1.0D0
  thresh_cov2=1.0D0

  !Number of values per line
  NEXP=0

  !write(*,*)cov_matrix
  open(unit=itubeCOVinp, file='INP.d/tube_COV.inp', status='unknown', action='read')
  call skip_comment(itubeCOVinp)
  read(itubeCOVinp,*)NLines1 !number of total species with coverage effects
  Nlines=nint(Nlines1)
  call skip_comment(itubeCOVinp)
  read(itubeCOVinp,*) cov_factor_adjust, cov_model
  call skip_comment(itubeCOVinp)
  !----opening lines read
  !-------------------------

  do_c1: do ca=1,NLines  !loop over species that are affected and store the interaction parameters
    read(itubeCOVinp,'(a)')LINE
    !Internal read to get number of following lines and the species
    read(line,*) species, NCov
    call SKSNUM(species,NEXP,LOUT,knams,kmax,PhaseNames,nphases,SpecTot,Nspec,NT,NVAL,RVAL(1),KERR)
    if (Nspec-kgmax<1 .or. Nspec-kgmax>ksmax) then
      write(*,*) 'Non-surface species in line: '
      write(*,*) line
      write(*,*) 'in tube_COV.inp. Stopping.'
      close(itubeCOVinp)
      STOP
    else if (kerr) then
      write(*,*) 'Species name specified in line: '
      write(*,*) line
      write(*,*) 'in tube_COV.inp does not appear in mechanism. Stopping.'
      close(itubeCOVinp)
      STOP
    end if
    do_c2: do cb=1,NCov  !loop over coverage effects
      !Read the line
      read(itubeCOVinp,'(a)')LINE
      !Internal read to extract the species
      read(line,*) species
      !Identify the species in it
      call SKSNUM(species,NEXP,LOUT,knams,kmax,Phasenames,nphases,spectot,Nspec2,NT,NVAL,RVAL(1),KERR)
      if (Nspec2-kgmax<1 .or. Nspec2-kgmax>ksmax) then
        write(*,*) 'Non-surface species in line: '
        write(*,*) line
        write(*,*) 'in tube_COV.inp. Stopping.'
        close(itubeCOVinp)
        STOP
      else if (kerr) then
        write(*,*) 'Species name specified in line: '
        write(*,*) line
        write(*,*) 'in tube_COV.inp does not appear in mechanism. Stopping.'
        close(itubeCOVinp)
        STOP
      end if
      !Internal read to parse the line for its numerical entries
      RVAL=0.0D0  !Set all parameters to zero initially
      select case (cov_model)
      case (1)
        read(line,*) species, RVAL(1)
        cov_matrix(Nspec2-kgmax,Nspec-kgmax)=RVAL(1)
      case (3)
        read(line,*) species, RVAL(1), RVAL(2), RVAL(3)
        cov_matrix(Nspec2-kgmax,Nspec-kgmax)=RVAL(1)
        thresh_cov(Nspec2-kgmax,Nspec-kgmax)=RVAL(2)
        cov_matrix2(Nspec2-kgmax,Nspec-kgmax)=RVAL(3)
      case (5)
        read(line,*) species, RVAL(1), RVAL(2), RVAL(3), RVAL(4), RVAL(5)
        cov_matrix(Nspec2-kgmax,Nspec-kgmax)=RVAL(1)
        thresh_cov(Nspec2-kgmax,Nspec-kgmax)=RVAL(2)
        cov_matrix2(Nspec2-kgmax,Nspec-kgmax)=RVAL(3)
        thresh_cov2(Nspec2-kgmax,Nspec-kgmax)=RVAL(4)
        cov_matrix3(Nspec2-kgmax,Nspec-kgmax)=RVAL(5)
      end select
    end do do_c2
    call skip_comment(itubeCOVinp)
  end do do_c1

  !Convert interaction parameters from kcal/mol/ML to K/ML
  cov_matrix=cov_matrix/Rgas_kcal
  cov_matrix2=cov_matrix2/Rgas_kcal
  cov_matrix3=cov_matrix3/Rgas_kcal

  close(itubeCOVinp)

!*********************************************************************************************************
end subroutine inco_cov
!*********************************************************************************************************

!*********************************************************************************************************
subroutine inco_StatpQ
!This subroutine reads in the StatpQ values.
!*********************************************************************************************************

  use global_all; use global_tube1; use update_tube2; use file_handles

  implicit none

  double precision               ::      RVAL(1),Nlines1
  integer                        ::      NCov,NLines,ca,cb,NEXP,spec1TOT,spec2TOT,NT,NVAL,i
  integer                        ::      Nspec,Nspec2
  character(len=80)              ::      LINE
  character(len=80)              ::      crap2
  logical                        ::      KERR
!---------------------------------------------------------------------------------------------
  !StatpQ is an array which gives the approximate T dependence of the binding energy
  StatpQ=0.0D0

  NEXP=1
  open(unit=iStatpQinp, file='INP.d/StatpQ.inp', status='unknown',action='read')
  call skip_comment(iStatpQinp)
  read(iStatpQinp,*) nlines
  call skip_comment(iStatpQinp)
  do i=1,nlines
    read(iStatpQinp,'(A)')LINE
    call SKSNUM(LINE,NEXP,LOUT,knams,kmax,PhaseNames,nphases,SpecTot,Nspec,NT,NVAL,RVAL(1),KERR)
    if (kerr) then
      write(*,*) 'Invalid species:'
      write(*,*) line
      write(*,*) 'in StatpQ.inp. Stopping.'
      close(iStatpQinp)
      STOP
    end if
    StatpQ(Nspec)=RVAL(1)
  end do
  close(iStatpQinp)

!*********************************************************************************************************
end subroutine inco_StatpQ
!*********************************************************************************************************

!*********************************************************************************************************
subroutine spec_vec_init(nlines, file_LU, SpecID_out)
!This subroutine initializes the species ID vectors for the various species-related
!input files. It assumes that the file position marker is in the correct location.
!*********************************************************************************************************

  use global_all; use global_tube1; use update_tube2; use file_handles

  implicit none

  integer, intent(in) :: nlines, file_LU
  integer, intent(out) :: SpecID_out(nlines)
  integer :: ierr, i, SpecID, NT, NVAL
  logical :: kerr
  double precision :: RVAL(1)
  character(len=33) :: species

  do i=1,nlines
    read(file_LU,*) species
    call sksnum(species,0,LOUT,knams,kmax,PhaseNames,nphases,SpecTot,SpecID,NT,NVAL,RVAL(1),KERR)
    if (kerr) then
      write(*,*) 'Processing error in species entry: ', species, '. Stopping.'
      STOP
    end if
    SpecID_out(i)=SpecID
  end do

  !Rewind to the starting point
  do i=1,nlines
    backspace(file_LU)
  end do

!*********************************************************************************************************
end subroutine spec_vec_init
!*********************************************************************************************************

!*********************************************************************************************************
subroutine rxn_vec_init(nlines, file_LU, RxnID_out)
!This subroutine initializes the reaction ID vectors for the various reaction-related
!input files. It assumes that the file position marker is in the correct location.
!*********************************************************************************************************

  use global_all; use global_tube1; use update_tube2; use file_handles

  implicit none

  integer, intent(in) :: nlines, file_LU
  integer, intent(out) :: RxnID_out(nlines)
  integer :: i, RxnID, NT
  character(len=80) :: rxn

  do i=1,nlines
    read(file_LU,*) rxn
    call skcomp(rxn,rxn_str,nrxn,RxnID,NT)
    if (RxnID==0) then !Reaction not found
      write(*,*) 'Invalid reaction specified in EA*.inp: '
      write(*,*) rxn
      write(*,*) 'Stopping.'
      STOP
    end if
    RxnID_out(i)=RxnID
  end do

  !Rewind to the starting point
  do i=1,nlines
    backspace(file_LU)
  end do

!*********************************************************************************************************
end subroutine rxn_vec_init
!*********************************************************************************************************

!*********************************************************************************************************
subroutine BEP_rxn_read(nlines,nseries,file_LU,RxnID)
!This subroutine reads the BEP-style reaction specifications in BEP.inp
!and mvn_EA.inp. It also does initial error checking.
!*********************************************************************************************************

  use global_all; use file_handles

  implicit none

  integer, intent(in) :: nlines, file_LU, nseries
  integer, intent(out) :: RxnID(nlines,3)
  integer :: i, RID, NT
  character(len=80) :: rxn(nlines), file_name

  !Set the file name for error handling later
  if (file_LU==iBEPinp) then
    file_name='BEP.inp'
  else if (file_LU==iSAinp) then
    file_name='SA.inp'
  end if

  !Read the data
  do i=1,nlines
    read(file_LU,*) RxnID(i,:2), rxn(i)
    call skcomp(rxn(i),rxn_str,nrxn,RID,NT)
    RxnID(i,3)=RID
  end do

  !Now check the RxnID array for correctness & consistency
  !Test 1: homologous series flags -- a 0 means a dummy entry
  if (any(RxnID(1:nlines,1)<0) .or. any(RxnID(1:nlines,1)>nseries)) then
    write(*,*) 'Invalid homologous series flag(s) detected in reaction(s):'
    do i=1,nlines
      if (RxnID(i,1)<0 .or. RxnID(i,1)>nseries) write(*,*) '  ', rxn(i)
    end do
    write(*,*) 'in ', trim(adjustl(file_name)), '. Stopping'
    STOP
  end if
  !Test 2: Reaction direction flags -- ignore dummy entries
  if (any((abs(RxnID(1:nlines,2))/=1) .and. (RxnID(1:nlines,1)/=0))) then
    write(*,*) 'Invalid reaction direction flag(s) detected in reaction(s):'
    do i=1,nlines
      if ((abs(RxnID(i,2))/=1) .and. (RxnID(i,1)/=0)) write(*,*) '  ', rxn(i)
    end do
    write(*,*) 'in ', trim(adjustl(file_name)), '. Stopping'
    STOP
  end if
  !Test 3: Reaction not found -- ignore dummy entries
  if (any((RxnID(1:nlines,3)==0) .and. (RxnID(1:nlines,1)/=0))) then
    write(*,*) 'Invalid reaction string(s) detected in reaction(s):'
    do i=1,nlines
      if ((RxnID(i,3)==0) .and. (RxnID(i,1)/=0)) write(*,*) '  ', rxn(i)
    end do
    write(*,*) 'in ', trim(adjustl(file_name)), '. Stopping'
    STOP
  end if
  !Test 4: User-specifiec activation energies -- only applies to BEPs
  if (file_LU==iBEPinp) then
    if (any(lEA_user(RxnID(1:nlines,3)) .and. RxnID(1:nlines,1)/=0)) then
      write(*,*) 'Both BEP and user-specified activation energies specified for reaction(s)'
      do i=1,nlines
        if ((lEA_user(RxnID(i,3))) .and. RxnID(i,1)/=0) write(*,*) '  ', rxn(i)
      end do
      write(*,*) 'The BEP value will overwrite the value specified in EA?.inp'
    end if
  end if

!*********************************************************************************************************
end subroutine BEP_rxn_read
!*********************************************************************************************************

!*********************************************************************************************************
subroutine inco_tube_SA
!This subroutine will read in the parameters associated with the sensitivity
!analysis given in SA.inp. It also allocates the necessary memory for the
!associated arrays.
!*********************************************************************************************************

  use file_handles; use globalSA; use mt19937_module; use global_all; use global_tube1
  use update_tube2

  implicit none

  !General values
  integer(4) :: seed
  integer :: nSAparams, SA_param_type, nlines
  integer :: i, j, k, correlation, NT, NVAL, Nspec, RxnID
  integer :: nrxn_list(kgmax+ksmax)
  integer :: basis_spec(nelm), pert_spec(kgmax+ksmax), basis_elem(nelm)
  integer :: nbasis_spec, npert_spec, nbasis_elem, nrxn_pert
  double precision :: temp_val(2)
  character(len=80) :: buffer
  logical :: kerr
  integer :: stoich_local(kmax-kbmax,nrxn)  !Work array for stoichiometric matrix
  integer :: RxnInfo(nrxn,3)
  integer :: skip

  open(unit=iSAinp, file='INP.d/SA.inp',status='old',action='read')
  call skip_comment(iSAinp)
  read(iSAinp,*) isenbrute
  !Decide what type (if any) of SA is performed. If the SA is not performed
  !(includes invalid options) or is FIM-based, then the rest of the file is
  !irrelevant and can be safely ignored. In that case, we close the file and
  !bail out of the routine. Note: The FIM-based SA is performed on all
  !reactions, so its array size is known a priori.
  if (isenbrute<-2) then
    model_solve=.false.
    isenbrute=-isenbrute
  else
    model_solve=.true.
  end if
  select case (isenbrute)
  case default  !Invalid option
    write(*,*) 'Warning: Invalid SA type selected. No SA performed.'
    lsenbrute=.false. !Later in inco_tube this will set isenbrute=0
    close(iSAinp)
    return
  case (0)  !No SA
    close(iSAinp)
    lsenbrute=.false.
    return
  case (1)  !FIM SA
    if (.not. allocated(FIM_SA_matrix)) allocate(FIM_SA_matrix(nrxn,2))
    if (.not. allocated(FIM_SA_matrix)) then
      write(*,*) 'FIM sensitivity analysis matrix not allocated. FIM SA not performed.'
      lsenbrute=.false.
    else
      FIM_SA_matrix=0.0
    end if
    call skip_comment(iSAinp)
    read(iSAinp,*) write_precision
    close(iSAinp)
    return
  case (2)  !Brute force LSA
    call skip_comment(iSAinp)
    read(iSAinp,*) write_precision, fix_EA
    SA_restart=.false.
    restart_on_failure=.false.
    nSAruns=1
    iter_write_freq=1
  case (3:5)  !Brute force GSA and UQ
    backspace(iSAinp)
    read(iSAinp,*) skip, seed, SA_restart, nSAruns, iter_write_freq, restart_on_failure, flush_buffer
    if (seed/=0) then
      call sgrnd(seed) !If non-default seed requested, initialize the generator
    end if
    if (iter_write_freq<=0) then
      write(*,*) 'Warning: A nonpositive number was selected for output frequency. Resetting to 1.'
      iter_write_freq=1
    end if
    call skip_comment(iSAinp)
    read(iSAinp,*) write_precision, fix_EA
  end select

  call skip_comment(iSAinp)
  read(iSAinp,*) rel_pert
  !Relative perturbations are only used for LSA (2) and derivative-based GSA (4)
  if (isenbrute/=2 .and. isenbrute/=4) rel_pert=.false.

  call skip_comment(iSAinp)
  read(iSAinp,*) nSAparams

  !Read in each of the data blocks and allocate arrays as needed. If a parameter
  !type is duplicated, a fatal error is thrown and the program terminates. All
  !enthalpies and entropies are divided by R.
  do i=1,nSAparams
    call skip_comment(iSAinp)
    read(iSAinp,*) SA_param_type
    if (SA_param_type<1 .or. SA_param_type>max_SA_param_types) then
      write(*,*) 'Invalid parameter type ', SA_param_type, ' selected in SA.inp. Stopping.'
      STOP
    end if
    if (.not. perturb_param(SA_param_type)) then
      perturb_param(SA_param_type)=.true. !We have decided to perturb this
      select case (SA_param_type)
      case (globalSA_H) !Enthalpies

        read(iSAinp,*) lDFTH_SA, lGAH_SA, lLSR_SA
        read(iSAinp,*) pert_ref_spec, pert_rxn_prop, lmvtH_SA

        if (isenbrute==2) then
          lDFTH_SA=.true.
          lGAH_SA=.false.
          lLSR_SA=.false.
          lmvtH_SA=.false.
          pert_ref_spec=.false.
          pert_rxn_prop=.false.
        end if

        if (lLSR_SA .and. iScale==0) then
          write(*,*) 'Scaling relations not enabled in tube.inp but LSR perturbations requested in SA.inp. Stopping.'
          STOP
        end if

        if (lmvtH_SA) then
          lDFTH_SA=.true.
          pert_ref_spec=.false.
          pert_rxn_prop=.false.
          read(iSAinp,*) pert_local_dist(:,globalSA_H)
          read(iSAinp,*)  !Skip the lines
          read(iSAinp,*)
          if (.not. rel_pert) pert_local_dist(:,globalSA_H)=&
            pert_local_dist(:,globalSA_H)/Rgas_kcal
          call mvt_init(globalSA_H) !Subroutine to read the MVT input file
        else if (lDFTH_SA) then
          !In Hspec_pert_dist, the first row is for species energies, while the
          !second is for reaction energies.
          if (.not. allocated(Hspec_pert_dist)) allocate(Hspec_pert_dist(4,2))
          read(iSAinp,*) pert_local_dist(:,globalSA_H)
          read(iSAinp,*) Hspec_pert_dist(:,1)
          read(iSAinp,*) Hspec_pert_dist(:,2)
          if (.not. rel_pert) then
            pert_local_dist(:,globalSA_H)=&
              pert_local_dist(:,globalSA_H)/Rgas_kcal
            Hspec_pert_dist=Hspec_pert_dist/Rgas_kcal
          end if
        else
          read(iSAinp,*)  !Skip the lines
          read(iSAinp,*)
          read(iSAinp,*)
        end if

      case (globalSA_S) !Entropies

        read(iSAinp,*) lDFTS_SA, lGAS_SA
        read(iSAinp,*) pert_ref_spec, pert_rxn_prop, lmvtS_SA

        if (isenbrute==2) then  !Local SA, use DFT perturbation only
          lDFTS_SA=.true.
          lGAS_SA=.false.
          lmvtS_SA=.false.
          pert_ref_spec=.false.
          pert_rxn_prop=.false.
        end if

        if (lmvtS_SA) then
          lDFTS_SA=.true.
          pert_ref_spec=.false.
          pert_rxn_prop=.false.
          read(iSAinp,*) pert_local_dist(:,globalSA_S)
          read(iSAinp,*)  !Skip the lines
          read(iSAinp,*)
          if (.not. rel_pert) pert_local_dist(:,globalSA_S)=&
            pert_local_dist(:,globalSA_S)/Rgas_cal
          call mvt_init(globalSA_S) !Subroutine to read the MVT input file
        else if (lDFTS_SA) then
          !In Sspec_pert_dist, the first row is for species entropies, while the
          !second is for reaction entropies.
          if (.not. allocated(Sspec_pert_dist)) allocate(Sspec_pert_dist(4,2))
          read(iSAinp,*) pert_local_dist(:,globalSA_S)
          read(iSAinp,*) Sspec_pert_dist(:,1)
          read(iSAinp,*) Sspec_pert_dist(:,2)
          if (.not. rel_pert) then
            pert_local_dist(:,globalSA_S)=&
              pert_local_dist(:,globalSA_S)/Rgas_cal
            Sspec_pert_dist=Sspec_pert_dist/Rgas_cal
          end if
        else
          read(iSAinp,*)  !Skip the lines
          read(iSAinp,*)
        end if

      case (globalSA_A) !Pre-exponentials

        if (.not. allocated(Preexp_pert_dist)) allocate(Preexp_pert_dist(4))
        read(iSAinp,*) pert_local_dist(:,globalSA_A)
        read(iSAinp,*) Preexp_pert_dist

      case (globalSA_B) !Betas

        if (.not. allocated(Beta_pert_dist)) allocate(Beta_pert_dist(4))
        read(iSAinp,*) pert_local_dist(:,globalSA_B)
        read(iSAinp,*) Beta_pert_dist

      case (globalSA_EA) !Activation barriers

        read(iSAinp,*) lDFTEA_SA, lBEP_SA, lmvn_EA
        if (isenbrute==2) then  !Local SA, use DFT perturbation only
          lDFTEA_SA=.true.
          lBEP_SA=.false.
          lmvn_EA=.false.
        end if

        if (lBEP_SA .and. .not. lBEP) then
          write(*,*) 'BEP correlations not used in tube.inp but perturbations specified in SA.inp. Stopping.'
          STOP
        end if

        if (lmvn_EA) then
          lDFTEA_SA=.true.
          pert_ref_spec=.false.
          pert_rxn_prop=.false.
          read(iSAinp,*) pert_local_dist(:,globalSA_EA)
          read(iSAinp,*)  !Skip the line
          if (.not. rel_pert) pert_local_dist(:,globalSA_EA)=&
            pert_local_dist(:,globalSA_EA)/Rgas_kcal
          call mvn_init()
        else if (lDFTEA_SA) then
          if (.not. allocated(EA_pert_dist)) allocate(EA_pert_dist(4))
          read(iSAinp,*) pert_local_dist(:,globalSA_EA)
          read(iSAinp,*) EA_pert_dist
          if (.not. rel_pert) then
            pert_local_dist(:,globalSA_EA)=&
              pert_local_dist(:,globalSA_EA)/Rgas_kcal
            EA_pert_dist=EA_pert_dist/Rgas_kcal
          end if
        else
          read(iSAinp,*)  !Skip the line
          read(iSAinp,*)
        end if

      case (globalSA_BEP) !BEP coefficients

        if (.not. lBEP) then
          write(*,*) 'BEP correlations not used in tube.inp but perturbations specified in SA.inp. Stopping.'
          STOP
        end if
        if (.not. allocated(BEP_coeff_pert_dist)) allocate(BEP_coeff_pert_dist(nBEP,4))
        do j=1,nBEP
          call skip_comment(iSAinp)
          read(iSAinp,*) correlation
          if (correlation<1 .or. correlation>nBEP) then
            write(*,*) 'Invalid BEP correlation number ', correlation, ' specified in SA.inp. Stopping.'
            STOP
          end if
          read(iSAinp,*) BEP_coeff_pert_dist(correlation,1:2)
          read(iSAinp,*) BEP_coeff_pert_dist(correlation,3:4)
        end do
        if (.not. rel_pert) BEP_coeff_pert_dist(:,3:4)=BEP_coeff_pert_dist(:,3:4)/Rgas_kcal

      case (globalSA_LSR) !LSR coefficients

        if (iScale==0) then
          write(*,*) 'Linear scaling relations not used in tube.inp but perturbations specified in SA.inp. Stopping.'
          STOP
        end if
        if (.not. allocated(LSR_coeff_pert_dist)) allocate(LSR_coeff_pert_dist(nscale,4))
        do j=1,nscale
          call skip_comment(iSAinp)
          read(iSAinp,*) correlation
          if (correlation<1 .or. correlation>nscale) then
            write(*,*) 'Invalid LSR correlation number ', correlation, ' specified in SA.inp. Stopping.'
            STOP
          end if
          read(iSAinp,*) LSR_coeff_pert_dist(correlation,1:2)
          read(iSAinp,*) LSR_coeff_pert_dist(correlation,3:4)
        end do
        if (.not. rel_pert) &
          LSR_coeff_pert_dist(correlation,3:4)=LSR_coeff_pert_dist(correlation,3:4)/Rgas_kcal

      end select
    else
      write(*,*) 'Duplicate perturbations for parameter type ', SA_param_type, ' detected in SA.inp. Stopping.'
      close(iSAinp)
      STOP
    end if
  end do
  !Check for mutual exclusivity of certain parameter types
  !Species properties
  if ((perturb_param(globalSA_H) .or. perturb_param(globalSA_S)) .and. perturb_param(globalSA_LSR)) then
    write(*,*) 'Per-species and LSR property perturbations are mutually exclusive. Stopping.'
    close(iSAinp)
    STOP
  end if
  !Reaction properties
  if ((perturb_param(globalSA_A) .or. perturb_param(globalSA_B) .or. perturb_param(globalSA_EA)) &
    .and. perturb_param(globalSA_BEP)) then
    write(*,*) 'Per-reaction and BEP property perturbations are mutually exclusive. Stopping.'
    close(iSAinp)
    STOP
  end if

  !Calculate the total number of parameters needed
  nparams_total=0

  !Now read in the species/reactions which should be perturbed
  !Read in the number of perturbed species
  call skip_comment(iSAinp)
  read(iSAinp,*) nlines !Number of species to perturb

  if (.not.allocated(spec_pert)) allocate(spec_pert(kgmax+ksmax))
  if (.not.allocated(spec_pert)) then
    write(*,*) 'Array of perturbed species not allocated. Stopping.'
    close(iSAinp)
    STOP
  end if
  if (nlines<-1 .or. nlines>kgmax+ksmax) then
    write(*,*) 'Invalid number of species specified. Stopping.'
    close(iSAinp)
    STOP
  end if
  if (.not. any(perturb_param((/globalSA_H,globalSA_S/)))) then
    !No species parameters perturbed, so we can skip this block.
    do i=1,nlines
      read(iSAinp,*)
    end do
    spec_pert=.false.
  elseif (nlines==0) then
    spec_pert=.false. !No species parameters perturbed
    if (perturb_param(globalSA_H) .or. perturb_param(globalSA_S)) &
      write(*,*) 'Species perturbation requested in SA.inp but no species&
        & names specified. No species perturbation performed.'
    perturb_param(globalSA_H)=.false.
    perturb_param(globalSA_S)=.false.
  else if (nlines==-1) then !Perturb everything but the free sites/inert
    spec_pert=.true.
    spec_pert(site_type_lastID(1:nsphases+1))=.false.
    !Initialize the array of species IDs
    if (.not. allocated(spec_pert_mvt_ID)) &
      allocate(spec_pert_mvt_ID(count(spec_pert)))
    j=1
    do i=1,kgmax+ksmax
      if (.not. spec_pert(i)) cycle
      spec_pert_mvt_ID(j)=i
      j=j+1
    end do
  else
    spec_pert=.false.
    if (.not. allocated(spec_pert_mvt_ID)) allocate(spec_pert_mvt_ID(nlines))
    j=0
    do i=1,nlines
      read(iSAinp,'(A)') buffer
      call SKSNUM(buffer,0,LOUT,knams,kmax,PhaseNames,nphases,SpecTot,Nspec,NT,NVAL,temp_val(1),KERR)
      if (Nspec==0) then !Species not found
        write(*,*) 'Invalid species specified in SA.inp. Stopping.'
        write(*,*) buffer
        close(iSAinp)
        STOP
      else if (Nspec>kgmax+ksmax) then  !Bulk species
        write(*,*) 'Warning: Bulk species specified in SA.inp. Species ignored.'
      else if (any(site_type_LastID(1:nsphases+1)==Nspec)) then !Vacant site
        write(*,*) 'Warning: Inert species or vacant site specified in SA.inp. Species ignored.'
      else  !Valid species
        j=j+1
        spec_pert(Nspec)=.true.
        spec_pert_mvt_ID(j)=Nspec
      end if
    end do
    !If the number of valid species is less than the number of entries,
    !reallocate the list of species ID numbers
    if (j<nlines) then
      pert_spec(1:j)=spec_pert_mvt_ID(1:j)
      deallocate(spec_pert_mvt_ID)
      allocate(spec_pert_mvt_ID(j))
      spec_pert_mvt_ID=pert_spec(1:j)
    end if
  end if
  !Reset the spec_pert vector if no species properties perturbed
  if (count(perturb_param((/globalSA_H,globalSA_S/)))==0) then
    spec_pert=.false.
    pert_rxn_prop=.false.
    pert_ref_spec=.false.
  end if
  nspec_params=count(spec_pert) !This counts the number of perturbed species

  !Check for size consistency of MVT array(s) with number of perturbed species
  !and with each other
  !  Case 1: both H and S from MVT
  !  Case 2: H from MVT -- S either from SA.inp or ignored
  !  Case 3: S from MVT -- H either from SA.inp or ignored
  if (lmvtH_SA .and. lmvtS_SA) then !Case 1
    if (any(shape(Hspec_pert_dist)/=shape(Sspec_pert_dist))) then
      write(*,*) 'Enthalpy and entropy MVT distribution arrays are of&
        & different shapes. Stopping.'
      STOP
    end if
    if (size(Hspec_pert_dist,2)/=nspec_params) then
      write(*,*) 'The number of perturbed species does not match the number&
        & of variables in the MVT distribution arrays. Stopping.'
      STOP
    end if
  else if (lmvtH_SA) then           !Case 2
    if (size(Hspec_pert_dist,2)/=nspec_params) then
      write(*,*) 'The number of perturbed species does not match the number&
        & of variables in the enthalpy MVT distribution array. Stopping.'
      STOP
    end if
  else if (lmvtS_SA) then           !Case 3
    if (size(Sspec_pert_dist,2)/=nspec_params) then
      write(*,*) 'The number of perturbed species does not match the number&
        & of variables in the entropy MVT distribution array. Stopping.'
      STOP
    end if
  end if

  call skip_comment(iSAinp)
  read(iSAinp,*) nlines !Number of reactions to perturb
  if (.not.allocated(rxn_pert)) allocate(rxn_pert(nrxn))
  if (.not.allocated(rxn_pert)) then
    write(*,*) 'Array of perturbed reactions not allocated. Stopping.'
    close(iSAinp)
    STOP
  end if
  if (nlines<-1 .or. nlines>nrxn) then
    write(*,*) 'Invalid number of reactions specified. Stopping.'
    close(iSAinp)
    STOP
  end if
  if (lmvn_EA) then
    rxn_pert=.false.
    if (.not. allocated(EA_mvn_pert_info)) allocate(EA_mvn_pert_info(nrxn,2))
    if (.not.allocated(EA_mvn_pert_info)) then
      write(*,*) 'Array of perturbed reaction info for MVN EA sampling not allocated. Stopping.'
      close(iSAinp)
      STOP
    end if

    !Read the BEP-style reactions
    RxnInfo=0
    call BEP_rxn_read(nlines,nmvn_EA_series,iSAinp,RxnInfo(1:nlines,:))

    !Assign the information to the permanent arrays
    do i=1,nlines
      if (RxnInfo(i,1)==0) then !Ignore this reaction
        rxn_pert(RxnInfo(i,3))=.false.
        EA_mvn_pert_info(RxnInfo(i,3),:)=0
      else
        rxn_pert(RxnInfo(i,3))=.true.
        EA_mvn_pert_info(RxnInfo(i,3),:)=RxnInfo(i,:2)
      end if
    end do
  elseif (.not. any(perturb_param((/globalSA_A,globalSA_B,globalSA_EA/)))) then
    !The information in this block does not need to be processed because
    !no reaction parameters are perturbed.
    do i=1,nlines
      read(iSAinp,*)  !Skip this block
    end do
    rxn_pert=.false.
  else
    if (nlines==0) then
      if (any(perturb_param((/globalSA_A,globalSA_B,globalSA_EA/)))) then
        write(*,*) 'Perturbation of reaction parameters requested but no &
          &reactions specified in SA.inp. No reaction parameters perturbed.'
      end if
      rxn_pert=.false.  !No reaction parameters perturbed
      perturb_param(globalSA_A)=.false.
      perturb_param(globalSA_B)=.false.
      perturb_param(globalSA_EA)=.false.
    else if (nlines==-1) then
      rxn_pert=.true.
    else
      rxn_pert=.false.
      do i=1,nlines
        read(iSAinp,'(A)') buffer
        call skcomp(buffer,rxn_str,nrxn,RxnID,NT)
        if (RxnID==0) then !Reaction not found
          write(*,*) 'Invalid reaction specified in SA.inp. Stopping.'
          close(iSAinp)
          STOP
        end if
        rxn_pert(RxnID)=.true.
      end do
    end if
  end if
  !Reset the rxn_pert vector if no reaction properties perturbed directly or
  !indirectly
  if (count(perturb_param((/globalSA_A,globalSA_B,globalSA_EA/)))==0 &
    .and. .not. pert_rxn_prop) rxn_pert=.false.
  !Count the number of reaction properties which are perturbed. If only
  !indirect perturbations, set nrxn_params to 0. We start with the base
  !assumption of direct perturbation and correct for the case where only
  !indirect perturbations are made. This assumption also works for the no
  !perturbation (of either type) case.
  nrxn_params=count(rxn_pert) !This counts the number of perturbed reactions
  !Check for no direct perturbations with indirect perturbations
  if (count(perturb_param((/globalSA_A,globalSA_B,globalSA_EA/)))==0 &
    .and. pert_rxn_prop) nrxn_params=0

  !If species properties require thermodynamically consistent perturbations,
  !read the rest of the file.
  basis_spec=0
  if (pert_ref_spec .or. pert_rxn_prop) then
    call skip_comment(iSAinp)
    read(iSAinp,*) nlines !Number of basis species
    if (nlines>nelm) then
      write(*,*) 'The number of basis species must not be larger than the &
        &number of elements. Stopping.'
      STOP
    end if
    do i=1,nlines
      read(iSAinp,'(A)') buffer
      call SKSNUM(buffer,0,LOUT,knams,kmax,PhaseNames,nphases,SpecTot,Nspec,NT,NVAL,temp_val(1),KERR)
      if (Nspec==0) then !Species not found
        write(*,*) 'Invalid species specified in SA.inp. Stopping.'
        write(*,*) buffer
        close(iSAinp)
        STOP
      else if (Nspec>kgmax+ksmax) then  !Bulk species
        write(*,*) 'Error: Bulk species specified in SA.inp basis set. Stopping.'
        STOP
        close(iSAinp)
      else if (any(site_type_LastID(1:nsphases+1)==Nspec)) then !Vacant site
        write(*,*) 'Error: Inert species or vacant site specified in SA.inp basis set. Stopping.'
        STOP
        close(iSAinp)
      else  !Valid species
        basis_spec(i)=Nspec
      end if
    end do
    nbasis_spec=count(basis_spec/=0)
  end if

  !Close the file; we're done with it
  close(iSAinp)

  !Check for special treatment of basis species and revise number of species
  !parameters
  if (pert_rxn_prop .or. pert_ref_spec) then  !Set all basis species to .false.
    nspec_params=nspec_params-count(spec_pert(basis_spec(1:nbasis_spec)))
    spec_pert(basis_spec(1:nbasis_spec))=.false.
  end if
  if (pert_ref_spec) then !Correct for being able to perturb the first basis species
    nspec_params=nspec_params+1
    spec_pert(basis_spec(1))=.true.
  end if

  !Now we have the total number of species/reaction/correlations perturbed for
  !each parameter type. Add them together to obtain the total number of
  !parameters.
  nparams_total=count(perturb_param((/globalSA_H,globalSA_S/)))*nspec_params+&
    count(perturb_param((/globalSA_A,globalSA_B,globalSA_EA/)))*nrxn_params
  if (lBEP .and. perturb_param(globalSA_BEP)) nparams_total=nparams_total+2*nBEP
  if (iScale>0 .and. perturb_param(globalSA_LSR)) nparams_total=nparams_total+2*nscale

  !Allocate space for the sensitivity coefficient matrices
  if (.not. allocated(w_SA)) allocate(w_SA(nparams_total+2,neqns))
  if (.not. allocated(w_SA)) then
    write(*,*) 'Storage for SA mass solution vectors not allocated. Sensitivity analysis not performed.'
    lsenbrute=.false.
    return
  end if
  if (.not. allocated(rop_SA)) allocate(rop_SA(nparams_total+2,2*nrxn))
  if (.not. allocated(rop_SA)) then
    write(*,*) 'Storage for SA reaction rate vectors not allocated. Sensitivity analysis not performed.'
    lsenbrute=.false.
    return
  end if
  if (isenbrute==2) then  !Local SA
    if (.not. allocated(NSC_matrix)) allocate(NSC_matrix(nparams_total,kgmax))
    if (.not. (allocated(NSC_matrix))) then
      write(*,*) 'Local SA matrix not allocated. Sensitivity analysis not performed.'
      lsenbrute=.false.
      return
    end if
    NSC_matrix=0.0
  else  !Global SA
    if (.not. allocated(elem_effect_matrix)) allocate(elem_effect_matrix(nparams_total,kgmax,nBE_coords))
    if (.not. allocated(first_order_SI_matrix)) allocate(first_order_SI_matrix(nparams_total,kgmax,nBE_coords))
    if (.not. allocated(total_SI_matrix)) allocate(total_SI_matrix(nparams_total,kgmax,nBE_coords))
    if (.not. (allocated(elem_effect_matrix) .and. allocated(first_order_SI_matrix) &
      .and. allocated(total_SI_matrix))) then
      write(*,*) 'Global SA matrices not allocated. Sensitivity analysis not performed.'
      lsenbrute=.false.
      return
    end if
    if (.not. allocated(w_var_SI)) allocate(w_var_SI(nSAruns/iter_write_freq,kgmax,nBE_coords))
    if (.not. allocated(var_SI)) allocate(var_SI(kgmax))
    if (.not. (allocated(w_var_SI) .and. allocated(var_SI))) then
      write(*,*) 'Storage for SI mass solution vectors not allocated. Sensitivity analysis not performed.'
      lsenbrute=.false.
      return
    end if
    elem_effect_matrix=0.0
    first_order_SI_matrix=0.0
    total_SI_matrix=0.0
    w_var_SI=0.0
    var_SI=0.0
  end if

  !Now allocate memory for the perturbation array
  if (.not. allocated(SApert_data)) allocate(SApert_data(nparams_total,2))

  !Initialize the perturbation array
  j=1 !Current index of the perturbation
  do i=1,max_SA_param_types
    if (.not. perturb_param(i)) cycle
    select case (i)
    case (globalSA_H)
      do k=1,kgmax+ksmax
        if (.not. spec_pert(k)) cycle
        !The SApert(...) construction is a fast way to initialize the members of
        !the derived data type.
        SApert_data(j,:)=SApert(globalSA_H,k,(/0.0,0.0,0.0/),(/0.0,0.0/))
        j=j+1
      end do
    case (globalSA_S)
      do k=1,kgmax+ksmax
        if (.not. spec_pert(k)) cycle
        SApert_data(j,:)=SApert(globalSA_S,k,(/0.0,0.0,0.0/),(/0.0,0.0/))
        j=j+1
      end do
    case (globalSA_A)
      do k=1,nrxn
        if (.not. rxn_pert(k)) cycle
        SApert_data(j,:)=SApert(globalSA_A,k,(/0.0,0.0,0.0/),(/0.0,0.0/))
        j=j+1
      end do
    case (globalSA_B)
      do k=1,nrxn
        if (.not. rxn_pert(k)) cycle
        SApert_data(j,:)=SApert(globalSA_B,k,(/0.0,0.0,0.0/),(/0.0,0.0/))
        j=j+1
      end do
    case (globalSA_EA)
      do k=1,nrxn
        if (.not. rxn_pert(k)) cycle
        SApert_data(j,:)=SApert(globalSA_EA,k,(/0.0,0.0,0.0/),(/0.0,0.0/))
        j=j+1
      end do
    case (globalSA_BEP)
      do k=1,nBEP
        !Assign the slope as +k and the intercept as -k
        SApert_data(j,:)=SApert(globalSA_BEP,k,(/0.0,0.0,0.0/),(/0.0,0.0/))
        SApert_data(j+1,:)=SApert(globalSA_BEP,-k,(/0.0,0.0,0.0/),(/0.0,0.0/))
        j=j+2
      end do
    case (globalSA_LSR)
      do k=1,nscale
        !Assign the slope as +k and the intercept as -k
        SApert_data(j,:)=SApert(globalSA_LSR,k,(/0.0,0.0,0.0/),(/0.0,0.0/))
        SApert_data(j+1,:)=SApert(globalSA_LSR,-k,(/0.0,0.0,0.0/),(/0.0,0.0/))
        j=j+2
      end do
    end select
  end do

  !Initialize the species/reaction index list if required by fix_EA
  if (.not. fix_EA) then

    !Count the number of reactions for gas and surface mechanisms for each species
    nrxn_list=count(stoich_matrix(1:kgmax+ksmax,:)<0 .and. &
      spread(spec_pert,2,nrxn),2)

    !Allocate memory for the list
    if (.not. allocated(spec_rxn_list)) &
      allocate(spec_rxn_list(kgmax+ksmax,maxval(nrxn_list)+1))
    spec_rxn_list=0

    !Assign the number of reactions to the index list
    spec_rxn_list(:,1)=nrxn_list

    !Assign reaction numbers for mechanism
    do i=1,nrxn
      do j=1,kgmax+ksmax

        !First check to see if this species/reaction combination is unused, and
        !if so, skip to the next one
        if (stoich_matrix(j,i)>=0 .or. spec_rxn_list(j,1)==0) cycle

        !The following loop finds the first available location for the next
        !reaction index for this reactant species
        k=1 !Column 1 is the number of entries
        do
          k=k+1
          if (spec_rxn_list(j,k)==0 .or. k==maxval(nrxn_list)+1) exit
        end do

        !Now assign reaction number to list
        spec_rxn_list(j,k)=i
      end do
    end do

  end if

  !We are done with the LSA stuff now, so we can quit if this is a LSA run as
  !what remains is only used in the GSA
  if (isenbrute==2) return

  !Determine which species are actually important to the perturbed mechanism
  if (pert_ref_spec .or. pert_rxn_prop) then
    ref_spec=basis_spec(1)  !First basis species is reference species
    !Find the perturbed species for the transformation matrix
    j=1
    pert_spec=0
    do i=1,kgmax+ksmax
      !Don't add any species which are either not perturbed or part of the basis
      !set.
      if ((.not. spec_pert(i)) .or. any(basis_spec==i)) cycle
      pert_spec(j)=i
      j=j+1
    end do
    npert_spec=count(pert_spec/=0)

    !Find the basis elements
    j=1
    basis_elem=0
    do i=1,nelm
      if (any(elem_comp(i,basis_spec(1:nbasis_spec))/=0) .or. &
        any(elem_comp(i,pert_spec(1:npert_spec))/=0)) then
        !Add this to the set of basis elements
        basis_elem(i)=j
        j=j+1
      end if
    end do
    nbasis_elem=count(basis_elem/=0)

    !Save the species IDs
    if (.not. allocated(spec_pert_mat_ID)) allocate(spec_pert_mat_ID(npert_spec))
    spec_pert_mat_ID=pert_spec(1:npert_spec)
  end if

  !Finish initializing the species perturbation arrays
  if (pert_ref_spec) then
    !Calculate the transformation matrix
    if (.not. allocated(spec_pert_mat)) &
      allocate(spec_pert_mat(npert_spec,nbasis_spec))
    spec_pert_mat=0.0
    call find_thermo_basis(nbasis_spec,basis_spec(1:nbasis_spec),nbasis_elem,&
      basis_elem(1:nbasis_elem),npert_spec,pert_spec(1:npert_spec),spec_pert_mat)

    ! Write out transformation matrix to file
    open(unit=itranmatout,file='OUT.d/trans_matrix.out',status='replace',action='write')
    write(itranmatout,'(17X,10(1X,A16))') &
      (knams(basis_spec(i)),i=1,nbasis_spec)
    do i=1,npert_spec
      write(itranmatout,'(1X,A16,10(1X,F16.4))') &
        knams(pert_spec(i)),spec_pert_mat(i,:)
    end do
    close(itranmatout)
  end if

  !Finish initializing the reaction property perturbation matrix
  if (pert_rxn_prop) then

    !Initialize working array, throwing away bulk species
    stoich_local=stoich_matrix(1:kmax-kbmax,:)

    !Throw out other, unused species
    do i=kgmax+ksmax,1,-1
      !The following test will ensure no basis species (including the first!)
      !are perturbed using reaction properties. It also ensures no species
      !containing an element /not/ in the basis set is discarded. It should
      !ensure the resulting matrix is full rank. We have to go backwards through
      !the species so that we don't accidentally skip any.
      if (.not. any(pert_spec==i) .or. any(basis_spec==i)) then
        !Throw away this species/reaction info with an end-off array shift
        stoich_local(i:,:)=eoshift(stoich_local(i:,:),shift=1,dim=1)
      end if
    end do
    nrxn_pert=count(rxn_pert)

    !Throw out unused/unperturbed reactions and store the IDs of the kept ones
    if (.not. allocated(rxn_pert_mat_ID)) &
      allocate(rxn_pert_mat_ID(nrxn_pert))
    rxn_pert_mat_ID=0
    j=1
    do i=1,nrxn
      if (.not. rxn_pert(i)) then
        stoich_local(:,i:)=eoshift(stoich_local(:,i:),shift=1,dim=2)
      else
        rxn_pert_mat_ID(j)=i
        j=j+1
      end if
    end do

    !Now calculate the transformation matrix for the least squares regression
    !of the species property perturbations.
    if (.not. allocated(rxn_pert_mat)) &
      allocate(rxn_pert_mat(npert_spec,nrxn_pert))
    rxn_pert_mat=0.0
    call find_rxn_pert_mat(npert_spec,nrxn_pert,&
      transpose(stoich_local(1:npert_spec,1:nrxn_pert)),rxn_pert_mat)

  end if

  !Read in group additivity related information if applicable
  if (lGAH_SA .or. lGAS_SA) then  !GA-based species thermochemistry used
    open(unit=iGAinp, file='INP.d/GA.inp', status='unknown', action='read')
    call skip_comment(iGAinp)
    if (lGAH_SA .and. .not. allocated(GA_Hest_pert_dist)) allocate(GA_Hest_pert_dist(4))
    if (lGAS_SA .and. .not. allocated(GA_Sest_pert_dist)) allocate(GA_Sest_pert_dist(4))
    if (lGAH_SA) then
      read(iGAinp,*) GA_Hest_pert_dist
      GA_Hest_pert_dist=GA_Hest_pert_dist/Rgas_kcal
    else
      read(iGAinp,*)  !Skip the line
    end if
    call skip_comment(iGAinp)
    if (lGAS_SA) then
      read(iGAinp,*) GA_Sest_pert_dist
      GA_Sest_pert_dist=GA_Sest_pert_dist/Rgas_cal
    else
      read(iGAinp,*)  !Skip the line
    end if
    call skip_comment(iGAinp)
    read(iGAinp,*) nlines
    if (.not. allocated(GA_species)) allocate(GA_species(kgmax+ksmax))
    GA_species=.false.
    call skip_comment(iGAinp)
    do j=1,nlines
      read(iGAinp,'(A)') buffer
      call SKSNUM(buffer,0,LOUT,knams,kmax,PhaseNames,nphases,SpecTot,Nspec,NT,NVAL,temp_val(1),KERR)
      if (Nspec==0) then !Species not found
        write(*,*) 'Invalid species specified in GA.inp. Stopping.'
        close(iGAinp)
        STOP
      else if (Nspec>kgmax+ksmax) then  !Bulk species
        write(*,*) 'Warning: Bulk species specified in GA.inp. Species ignored.'
      else if (any(site_type_lastID==Nspec)) then  !Inert/vacant site
        write(*,*) 'Warning: Inert species or vacant site specified in GA.inp. Species ignored.'
      else
        GA_species(Nspec)=.true.
      end if
    end do
    close(iGAinp)
  end if

  !If applicable, read in the saved statistics and PRNG state
  if (SA_restart) then
    inquire(file='INP.d/PRNG_save.inp',exist=kerr)
    inquire(file='OUT.d/tube_w_sen.out',exist=kerr)
    inquire(file='OUT.d/tube_rate_sen.out',exist=kerr)
    if (.not. kerr) then
      write(*,*) 'Warning: SA restart requested but saved state files do not exist.'
      write(*,*) 'Starting instead with default values.'
      SA_restart=.false.
    else
      write(*,*) 'Restarting SA with PRNG state from previous runs.'
      write(*,*) 'New results will be appended to previous results.'
      open(unit=iPRNGstateinp,file='INP.d/PRNG_save.inp',action='read')
      call mtget(iPRNGstateinp,'f')  !Read the PRNG state
      close(iPRNGstateinp)
    end if
  end if

!*********************************************************************************************************
end subroutine inco_tube_SA
!*********************************************************************************************************

!*********************************************************************************************************
subroutine skip_comment(LU)
!This subroutine will move the position marker in an open file to the last
!exclamation point-delimited comment line in a block of such comments.
!*********************************************************************************************************

  implicit none

  integer, intent(in) :: LU !File handle
  integer :: ierr
  character(len=80) :: line !Buffer

  do
    read(LU,*,iostat=ierr) line
    line=adjustl(line)
    if (line(1:1)/='!' .or. ierr<0) then
      backspace(LU)
      exit
    end if
  end do

!*********************************************************************************************************
end subroutine skip_comment
!*********************************************************************************************************

!*********************************************************************************************************
subroutine make_rect_grid(BE_delta, BE_range)
!This subroutine makes an isotropic rectangular grid in an arbitrary number of
!dimensions.
!*********************************************************************************************************

  use global_all

  implicit none

  double precision, intent(in) :: BE_delta, BE_range(2,natoms_scale)
  integer :: nscale_iter(natoms_scale)
  integer :: i, j
  integer :: scale_divisor(natoms_scale), scale_counter(natoms_scale)
  double precision :: BE_origin(natoms_scale)

  !Calculate the number of mesh points in each direction and the BE origin
  !The BE origin is shifted so that any left-over distance is slpit evenly
  !between the low and high ends of the range
  do i=1,natoms_scale
    nscale_iter(i)=floor((BE_range(2,i)-BE_range(1,i))/BE_delta)+1
    BE_origin(i)=(BE_range(2,i)-BE_range(1,i)-(nscale_iter(i)-1)*BE_delta)/2+BE_range(1,i)
  end do

  !Calculate the rollover points for an n-dimensional 'odometer'
  scale_counter=1
  scale_divisor(natoms_scale)=1
  do i=1,natoms_scale-1
    scale_divisor(i)=product(nscale_iter(i+1:natoms_scale))
  end do
  nBE_coords=product(nscale_iter)  !Total number of iterations

  !Allocate memory for the BE coordinate array
  if (.not. allocated(BE_coords)) allocate(BE_coords(natoms_scale,nBE_coords))

  !Initialize the BE coordinate array
  do j=1,nBE_coords
    BE_coords(:,j)=BE_origin+(scale_counter-1)*BE_delta
    !Update counter
    do i=1,natoms_scale
      if (mod(j,scale_divisor(i))==0) then !Increment odometer
        scale_counter(i)=scale_counter(i)+1
        if (scale_counter(i)>nscale_iter(i)) scale_counter(i)=1 !Digit overflow, reset to 1
      end if
    end do
  end do

!*********************************************************************************************************
end subroutine make_rect_grid
!*********************************************************************************************************

!*********************************************************************************************************
subroutine make_hex_grid(BE_delta, BE_range)
!This subroutine makes an isotropic hexagonol grid in up to three dimensions.
!*********************************************************************************************************

  use global_all

  implicit none

  double precision, intent(in) :: BE_delta, BE_range(2,natoms_scale)
  integer :: i, j, k
  integer :: nscale_iter(natoms_scale)
  double precision :: BE_origin(natoms_scale)
  double precision :: BE_spacing(natoms_scale)
  double precision :: BE_offset(natoms_scale)
  integer :: BE_no

  !Initialize the mesh spacings and number of points in each direction
  if (natoms_scale==2) then
    BE_spacing=(/1.,sqrt(3.)/2/)*BE_delta  !The 1 is to prevent divide by zero on the next line
  else
    BE_spacing=(/1.,sqrt(3.)/2,sqrt(2./3.)/)*BE_delta
  end if
  nscale_iter=floor((BE_range(2,:)-BE_range(1,:))/BE_spacing)+1
  !General case for BE origin
  BE_origin(:)=(BE_range(2,:)-BE_range(1,:)-(nscale_iter(:)-1)*BE_spacing)/2+BE_range(1,:)

  !Count the total number of points in the mesh
  nBE_coords=0
  do i=1,nscale_iter(1)
    do j=1,nscale_iter(2)
      !The last point on an even row is skipped
      if (i==nscale_iter(1) .and. mod(j,2)==0) cycle
      nBE_coords=nBE_coords+1
    end do
  end do
  if (natoms_scale==3) nBE_coords=nBE_coords*nscale_iter(3)

  !Allocate memory
  if (.not. allocated(BE_coords)) allocate(BE_coords(natoms_scale,nBE_coords))

  !Assign binding energies to first layer
  BE_no=0
  do i=1,nscale_iter(1)
    do j=1,nscale_iter(2)
      !The last point on an even row is skipped
      if (i==nscale_iter(1) .and. mod(j,2)==0) cycle
      BE_no=BE_no+1
      if (natoms_scale==2) then !2-D case
        if (mod(j,2)==0) then !Even row
          BE_coords(:,BE_no)=BE_origin+((/i,j/)-(/0.5,1./))*BE_spacing
        else  !Odd row
          BE_coords(:,BE_no)=BE_origin+((/i,j/)-1)*BE_spacing
        end if
      else  !3-D case
        if (mod(j,2)==0) then !Even row
          BE_coords(:,BE_no)=BE_origin+((/i,j,1/)-(/0.5,1.,1./))*BE_spacing
        else  !Odd row
          BE_coords(:,BE_no)=BE_origin+((/i,j,1/)-1)*BE_spacing
        end if
      end if
    end do
  end do

  !Return from subroutine if 2-D case or 1-layer 3-D because we are now done
  if (natoms_scale==2) return
  if (nscale_iter(3)==1) return

  !For 3-D case with additional layers, use the first layer as initial mesh
  !and then adjust the coordinates for additional y, z offsets. Each additional
  !layer has the same number of grid points as the first.
  do i=1,nscale_iter(3)-1
    BE_coords(:,i*BE_no+1:(i+1)*BE_no)=BE_coords(:,(i-1)*BE_no+1:i*BE_no)
    !Fix z coordinate
    BE_coords(3,i*BE_no+1:(i+1)*BE_no)=BE_coords(3,i*BE_no+1:(i+1)*BE_no)+BE_spacing(3)
    !Fix y coordinate
    if (mod(i,2)==1) then !Even layer
      BE_coords(2,i*BE_no+1:(i+1)*BE_no)=BE_coords(2,i*BE_no+1:(i+1)*BE_no)+BE_spacing(2)
    else  !Odd layer
      BE_coords(2,i*BE_no+1:(i+1)*BE_no)=BE_coords(2,i*BE_no+1:(i+1)*BE_no)-BE_spacing(2)
    end if
  end do

  !Now recenter the grid along the y coordinate
  BE_coords(2,:)=BE_coords(2,:)-(maxval(BE_coords(2,:))-BE_range(2,2))
  BE_coords(2,:)=BE_coords(2,:)+(BE_range(1,2)-minval(BE_coords(2,:)))/2

!*********************************************************************************************************
end subroutine make_hex_grid
!*********************************************************************************************************

!*********************************************************************************************************
subroutine make_sobol_grid(BE_range,skip)
!This subroutine makes a BE grid using a Sobol' low discrepancy sequence.
!*********************************************************************************************************

  use global_all; use sobol

  implicit none

  double precision, intent(in) :: BE_range(2,natoms_scale)
  double precision :: sobol_seq(natoms_scale,nBE_coords)
  integer, intent(in) :: skip

  !Generate the Sobol' sequence on the uniform hypercube
  call i8_sobol_generate(natoms_scale,nBE_coords,skip,sobol_seq)

  !Allocate memory for the BE array
  if (.not. allocated(BE_coords)) allocate(BE_coords(natoms_scale,nBE_coords))

  !Convert the Sobol' sequence into the BE array via centering and scaling
  BE_coords=(spread(BE_range(2,:),2,nBE_coords)-spread(BE_range(1,:),2,nBE_coords))*&
    sobol_seq+spread(BE_range(1,:),2,nBE_coords)

!*********************************************************************************************************
end subroutine make_sobol_grid
!*********************************************************************************************************

!*********************************************************************************************************
subroutine find_thermo_basis(nbasis_spec,basis_spec,nbasis_elem,basis_elem,&
  npert_spec,pert_spec,trans_mat)
!This subroutine is used to find the transformation matrix needed for
!perturbing species energies/entropies in a thermodynamically consistent
!manner.
!*********************************************************************************************************

  use global_all

  implicit none

  !Number of basis species, basis elements, and perturbed species.
  integer, intent(in) :: nbasis_spec, nbasis_elem, npert_spec
  !Arrays of IDs for the basis, perturbed species, and basis elements
  integer, intent(in) :: basis_spec(nbasis_spec), pert_spec(npert_spec)
  integer, intent(in) :: basis_elem(nbasis_elem)
  !Transformation matrix to generate the remaining species compositions from
  !the basis species.
  double precision, intent(out) :: trans_mat(npert_spec,nbasis_spec)
  !Elemental compositions for perturbed and basis species
  integer :: pert_spec_comp(npert_spec,nbasis_elem)
  integer :: basis_spec_comp(nbasis_spec,nbasis_elem)
  !Work arrays for inversion routine
  double precision :: work(nbasis_spec), det(2), rcond
  double precision :: basis_inv(nbasis_spec,nbasis_spec)
  integer :: ipvt(nbasis_spec)
  !Counters
  integer :: i, j, ii, jj

  !Initialize compositions for perturbed species
  ii=1
  do i=1,kmax
    if (.not. any(pert_spec==i)) cycle
    jj=1
    do j=1,nelm
      if (.not. any(basis_elem==j)) cycle
      pert_spec_comp(ii,jj)=elem_comp(j,i)
      jj=jj+1
    end do
    ii=ii+1
  end do

  !Initialize compositions for basis species
  ii=1
  do i=1,kmax
    if (.not. any(basis_spec==i)) cycle
    jj=1
    do j=1,nelm
      if (.not. any(basis_elem==j)) cycle
      basis_spec_comp(ii,jj)=elem_comp(j,i)
      jj=jj+1
    end do
    ii=ii+1
  end do

  !Perform the inversion
  basis_inv=dble(basis_spec_comp)
  call dgeco(basis_inv,nbasis_spec,nbasis_spec,ipvt,rcond,work)
  if (1.0 + rcond .eq. 1.0) then
    write(*,*) 'Singular determinant in inverting basis matrix for global &
      &sensitivity analysis'
    STOP
  end if
  call dgedi(basis_inv,nbasis_spec,nbasis_spec,ipvt,det,work,11)
  if (1.0 + (det(1) * 10**det(2)) .eq. 1.0) then ! Double-check for singular matrix
    write(*,*) 'Singular determinant in inverting basis matrix for global &
      &sensitivity analysis'
    STOP
  end if

  !Calculate the transformation matrix
  trans_mat=matmul(pert_spec_comp,basis_inv)

!*********************************************************************************************************
end subroutine find_thermo_basis
!*********************************************************************************************************

!*********************************************************************************************************
subroutine find_rxn_pert_mat(nspec,nrxns,st_local,st_factor)
!This subroutine is responsible for calculating and returning the rxn_pert_mat matrix
!*********************************************************************************************************

  implicit none

  integer, intent(in) :: nspec, nrxns
  integer, intent(in) :: st_local(nrxns,nspec)
  double precision, intent(out) :: st_factor(nspec,nrxns)
  integer :: ipvt(nspec), i
  double precision :: work(nspec), det(2), rcond
  double precision :: stTst_inv(nspec,nspec)

  !Calculate nuT*nu (where nu is the stoichiometric matrix)
  !Note that st_local=nu!
  stTst_inv=dble(matmul(transpose(st_local),st_local))

  !Perform the inversion
  call dgeco(stTst_inv,nspec,nspec,ipvt,rcond,work)
  if (1.0 + rcond .eq. 1.0) then
    write(*,*) 'Singular determinant in inverting matrix for global &
      &sensitivity analysis for regression of species properties from &
      &reaction properties.'
    STOP
  end if
  call dgedi(stTst_inv,nspec,nspec,ipvt,det,work,11)
  if (1.0 + (det(1) * 10**det(2)) .eq. 1.0) then ! Double-check for singular matrix
    write(*,*) 'Singular determinant in inverting matrix for global &
      &sensitivity analysis for regression of species properties from &
      &reaction properties.'
    STOP
  end if

  !Complete the calculation of the pre-factor matrix
  st_factor=matmul(stTst_inv,transpose(st_local))

!*********************************************************************************************************
end subroutine find_rxn_pert_mat
!*********************************************************************************************************

!*********************************************************************************************************
subroutine mvt_init(file_type)
!This subroutine is responsible for reading the MVT input files and
!initializing the distribution parameters.
!*********************************************************************************************************

  use file_handles; use globalSA; use global_all; use stats

  implicit none

  integer, intent(in) :: file_type
  integer :: nspec, file_LU, i
  double precision, allocatable :: ch_work(:,:)  !Work array for Cholesky decomp.

  if (file_type/=globalSA_H .and. file_type/=globalSA_S) return
  if (file_type==globalSA_H) then
    file_LU=imvtHinp
    open(unit=file_LU,file='INP.d/mvt_spec_H.inp',action='read',status='old')
  else
    file_LU=imvtSinp
    open(unit=file_LU,file='INP.d/mvt_spec_S.inp',action='read',status='old')
  end if

  !Read number of species
  call skip_comment(file_LU)
  read(file_LU,*) nspec
  if (nspec<=0) then
    write(*,*) 'Non-positive value specified for number of species in&
      & MVT input file. Stopping.'
    STOP
  end if

  !Allocate memory for the distribution parameters; these arrays have
  !dimensions of (nspec+1,nspec). The additional row is for the vector
  !of means, while the remaining rows are the covariance matrix.
  if (file_type==globalSA_H) then
    if (.not. allocated(Hspec_pert_dist)) allocate(Hspec_pert_dist(nspec+1,nspec))
  else
    if (.not. allocated(Sspec_pert_dist)) allocate(Sspec_pert_dist(nspec+1,nspec))
  end if

  !Read the degrees of freedom
  call skip_comment(file_LU)
  read(file_LU,*) mvt_dof(file_type)
  if (mvt_dof(file_type)<1) then
    write(*,*) 'Invalid number of degrees of freedom specified in MVT&
      & input file. Recommended number is at least 2. Stopping.'
    STOP
  end if

  !Read the mean vector and normalize as needed
  call skip_comment(file_LU)
  if (file_type==globalSA_H) then
    read(file_LU,*) Hspec_pert_dist(1,:)
    Hspec_pert_dist(1,:)=Hspec_pert_dist(1,:)/Rgas_kcal
  else
    read(file_LU,*) Sspec_pert_dist(1,:)
    Sspec_pert_dist(1,:)=Sspec_pert_dist(1,:)/Rgas_cal
  end if

  !Read the covariance matrix and normalize as needed
  call skip_comment(file_LU)
  if (file_type==globalSA_H) then
    do i=2,nspec+1
      read(file_LU,*) Hspec_pert_dist(i,:)
    end do
    Hspec_pert_dist(2:,:)=Hspec_pert_dist(2:,:)/(Rgas_kcal**2)
  else
    do i=2,nspec+1
      read(file_LU,*) Sspec_pert_dist(i,:)
    end do
    Sspec_pert_dist(2:,:)=Sspec_pert_dist(2:,:)/(Rgas_cal**2)
  end if

  close(file_LU)

  !Perform the Cholesky decomposition
  allocate(ch_work(nspec,nspec))
  if (file_type==globalSA_H) then
    call cholesky(nspec,Hspec_pert_dist(2:,:),ch_work)
    Hspec_pert_dist(2:,:)=ch_work
  else
    call cholesky(nspec,Sspec_pert_dist(2:,:),ch_work)
    Sspec_pert_dist(2:,:)=ch_work
  end if
  deallocate(ch_work)

!*********************************************************************************************************
end subroutine mvt_init
!*********************************************************************************************************

!*********************************************************************************************************
subroutine mvn_init()
!This subroutine is responsible for reading the MVN input file and
!initializing the distribution parameters.
!*********************************************************************************************************

  use file_handles; use globalSA; use global_all

  implicit none

  integer :: i, j

  open(unit=imvninp,file='INP.d/mvn_EA.inp',action='read',status='old')

  !Read number of homologous series
  call skip_comment(imvninp)
  read(imvninp,*) nmvn_EA_series

  !Allocate memory for the homologous series distributions. Each level of
  !the 3-D array corresponds to a different homologous series and is of
  !dimensions (3,2). The extra row is for the means, just as with the MVT
  !distribution array.
  if (.not. allocated(EA_mvn_pert_dist)) &
    allocate(EA_mvn_pert_dist(3,2,nmvn_EA_series))

  !Now read the distribution parameters and normalize them
  do i=1,nmvn_EA_series
    call skip_comment(imvninp)
    do j=1,3
      read(imvninp,*) EA_mvn_pert_dist(j,:,i)
    end do
  end do
  EA_mvn_pert_dist(1,:,:)=EA_mvn_pert_dist(1,:,:)/Rgas_kcal
  EA_mvn_pert_dist(2:,:,:)=EA_mvn_pert_dist(2:,:,:)/(Rgas_kcal**2)

  !Done with the file; close it
  close(imvninp)

!*********************************************************************************************************
end subroutine mvn_init
!*********************************************************************************************************

!*********************************************************************************************************
end module input
!*********************************************************************************************************
