!*********************************************************************************************************
module output
!This module has subroutines used primarily for writing output.
!*********************************************************************************************************

  implicit none

contains

!*********************************************************************************************************
subroutine open_output
!This subroutine is responsible for opening the output files.
!*********************************************************************************************************

  use global_all; use global_tube1; use update_tube2; use file_handles; use globalSA; use mt19937_module

  implicit none

  integer :: mem_stat, i, j, k, str_len, ierr
  character(len=80) :: fmt_str, moment_str
  character(len=190) :: buffer
  character(len=20) :: stat, pos

  !Open output files
  !Save the version number and other information about the problem to assist in
  !automated (programmatic) post-processing
  call write_genout

  !Conversion
  open(unit=itubeconvout, file='OUT.d/tube_conv.out', status='replace', action='write')
  write(itubeconvout,*)' Run Temp. [K]    Conversion    Rate [TOF]'

  !Steady state reaction rates (forward and net)
  if (.not. allocated(rxn_rates_ss)) allocate(rxn_rates_ss(nrxn,2,nruns,nBE_coords),STAT=mem_stat)
  open(unit=irxnratess, file='OUT.d/rxn_rate_ss.out', status='replace', action='write')

  !Transient output
  if (ltra) then
    open(unit=icovtraout, file='OUT.d/tube_cov_tra.out', status='replace', action='write')
    open(unit=igasmoletraout, file='OUT.d/tube_gasmole_tra.out', status='replace', action='write')
    open(unit=igasmasstraout, file='OUT.d/tube_gasmass_tra.out', status='replace', action='write')
    open(unit=imassbaltraout, file='OUT.d/tube_mass_bal_tra.out', status='replace', action='write')
    if (surf_chem) then
      open(unit=isdottraout, file='OUT.d/tube_sdot_tra.out', status='replace', action='write')
    end if
    if (gas_chem) then
      open(unit=igdottraout, file='OUT.d/tube_gdot_tra.out', status='replace', action='write')
    end if
  end if

  !Steady state output
  !Note: itpd==1 means UHV TPD, 0 means off, 2 means high pressure
  !      irxtr==1 means batch reactor
  if (irxtr/=1 .and. itpd/=1) then  !Don't save SS info for the /batch/ reactor or UHV TPD
    open(unit=igasmolessout, file='OUT.d/tube_gasmole_ss.out', status='replace', action='write')
    open(unit=igasmassssout, file='OUT.d/tube_gasmass_ss.out', status='replace', action='write')
    open(unit=imassbalssout, file='OUT.d/tube_mass_bal_ss.out', status='replace', action='write')
    if (surf_chem) then
      open(unit=icovssout, file='OUT.d/tube_cov_ss.out', status='replace', action='write')
      open(unit=isdotssout, file='OUT.d/tube_sdot_ss.out', status='replace', action='write')
    end if
    if (gas_chem) then
      open(unit=igdotssout, file='OUT.d/tube_gdot_ss.out', status='replace', action='write')
    end if
  end if

  !Sensitivity analysis
  if (lsenbrute) then

    !Determine whether to append to files or not
    if (SA_restart .and. isenbrute>2) then  !This only applies to a GSA
      stat='old'
      pos='append'
    else
      stat='replace'
      pos='rewind'
    end if

    !Sensitivity coefficients
    open(unit=ifiletubesen, file='OUT.d/tube_sen.out', status=stat, action='write',position=pos)


    !Raw solution vectors -- for reconstructing all derived results with post-processing
    if (SA_restart .and. isenbrute>2) then !Read in the last replicate number
      open(unit=iwSAout, file='OUT.d/tube_w_sen.out', status='old', action='readwrite', position='append')
      !Read the number of completed iterations
      backspace(iwSAout)
      read(iwSAout,*) nSAruns_old
      !Rewind to start of the data block.
      select case (isenbrute)
      case (3)  !VGSA
        do i=1,nparams_total+5  !The 5 is for 2 design points + 1 initial composition + 2 header lines
          backspace(iwSAout)
        end do
      case (4)  !DGSA
        do i=1,nparams_total+4  !The 4 is for 1 design point + 1 initial composition + 2 header lines
          backspace(iwSAout)
        end do
      case (5)  !UQ -- no extra parameter-related entries
        do i=1,5  !The 5 is for 2 design points + 1 initial composition + 2 header lines
          backspace(iwSAout)
        end do
      end select
      !Read the number of binding energies completed for this replicate
      if (iScale>0) then
        backspace(iwSAout)  !Account for binding energy header
        read(iwSAout,'(1X,A42,I5)') buffer, nBEcoords_old
      else
        nBEcoords_old=1
      end if
      !Read the number of completed experiments
      read(iwSAout,'(1X,A24,I3)') buffer, nexp_old
      close(iwSAout)
      !Check for full completion of previous set of GSA replicates
      if (nSAruns_old>=nSAruns .and. nBEcoords_old==nBE_coords .and. nexp_old==nruns) then
        write(*,*) 'Last replicate number specified in sensitivity analysis output files &
          &is at least as large as the number of perturbations specified.'
        write(*,*) 'Please increase the number of perturbations allowed. Stopping.'
        STOP
      end if
    else
      nSAruns_old=0 !No previous SA runs
      if (iScale>0) then
        nBEcoords_old=nBE_coords  !'Finished' the last round of binding energies
      else
        nBEcoords_old=1 !Always at least 1 BE
      end if
      nexp_old=1  !Always at least 1 expt
    end if
    open(unit=iwSAout, file='OUT.d/tube_w_sen.out', status=stat, action='write',position=pos)
    if (pos/='append') write(iwSAout,*) 'Complete solution vectors for each perturbation state point'
    open(unit=irSAout, file='OUT.d/tube_rate_sen.out', status=stat, action='write',position=pos)
    if (pos/='append') write(irSAout,*) 'Complete reaction rate vectors for each perturbation state point'

    !Open output files for the perturbation parameters and the PRNG state
    if (isenbrute>2) then
      open(unit=iPRNGstateout,file='OUT.d/PRNG_save.inp',status='replace',action='write')

      do i=1,max_SA_param_types
        if (.not. perturb_param(i)) cycle
        select case (i)
        case (globalSA_H)
          open(unit=iHfpertout,file='OUT.d/Hform_pert.out',status=stat,action='write',position=pos)
          if (pos/='append') write(iHfpertout,*) 'Dimensionless Perturbations to Heats of Formation'
        case (globalSA_S)
          open(unit=iSfpertout,file='OUT.d/Sform_pert.out',status=stat,action='write',position=pos)
          if (pos/='append') write(iSfpertout,*) 'Dimensionless Perturbations to Entropies of Formation'
        case (globalSA_A)
          open(unit=iApertout,file='OUT.d/Preexp_pert.out',status=stat,action='write',position=pos)
          if (pos/='append') write(iApertout,*) 'Perturbations to Reaction Pre-exponentials'
        case (globalSA_B)
          open(unit=iBpertout,file='OUT.d/Beta_pert.out',status=stat,action='write',position=pos)
          if (pos/='append') write(iBpertout,*) 'Perturbations to Reaction Beta Values'
        case (globalSA_EA)
          open(unit=iEApertout,file='OUT.d/EA_pert.out',status=stat,action='write',position=pos)
          if (pos/='append') write(iEApertout,*) 'Dimensionless Perturbations to Reaction Activation Energies'
        case (globalSA_BEP)
          open(unit=iBEPpertout,file='OUT.d/BEP_pert.out',status=stat,action='write',position=pos)
          if (pos/='append') write(iBEPpertout,*) 'Dimensionless Perturbations to BEP Coefficients'
        case (globalSA_LSR)
          open(unit=iLSRpertout,file='OUT.d/LSR_pert.out',status=stat,action='write',position=pos)
          if (pos/='append') write(iLSRpertout,*) 'Dimensionless Perturbations to LSR Coefficients'
        end select
      end do
    end if

  end if

  if (irxtr>1) then !Restart is only applicable to a flow reactor
    open(unit=ituberestartout, file='OUT.d/tube_molerestart.out', status='replace', action='write')
    write(ituberestartout,163)itube_restart+1,'itube_restart'

    !Allocate memory for the restart_save array
    if (.not. allocated(restart_save)) &
      allocate(restart_save(kgmax+ksmax+1,nruns), STAT=mem_stat)
    if (.not. allocated(restart_save)) then
      write(*,*) 'Output array for restart not allocated. Stopping.'
      STOP
    end if
  end if

  !Reaction path analysis
  if (mrpa>0) then
    open(unit=ifiletuberpa, file='OUT.d/tube_rpa.out', status='replace', action='write')
    open(unit=irpavisout, file='OUT.d/rpa_vis_output.txt', status='replace', action='write')
    if (gas_chem .and. surf_chem) open(unit=ituberateout, file='OUT.d/tube_rates.out', status='replace', action='write')
    if (mrpa==2) open(unit=ifiletrpa, file='OUT.d/trpa.out', status='replace', action='write')
  end if

  !Kinetic and thermodynamic parameters
  open(unit=iStoichout,file='OUT.d/Stoich.out',status='replace',action='write')
  open(unit=ielemcomp,file='OUT.d/elem_comp.out',status='replace',action='write')
  call write_stoich
  open(unit=ikout, file='OUT.d/Preex_out.out', status='replace', action='write')
  open(unit=ibetaout, file='OUT.d/Beta_out.out', status='replace', action='write')
  open(unit=iEARTout, file='OUT.d/Ea_over_RT.out', status='replace', action='write')
  open(unit=iEQKCout, file='OUT.d/EQKC_out.out', status='replace', action='write')
  open(unit=iDSRout, file='OUT.d/Srxn_out.out', status='replace', action='write')
  open(unit=iDHRTout, file='OUT.d/Hrxn_out.out', status='replace', action='write')
  open(unit=iDGRTout, file='OUT.d/Grxn_out.out', status='replace', action='write')
  open(unit=iHRTfout, file='OUT.d/Hform_out.out', status='replace', action='write')
  open(unit=iSRfout, file='OUT.d/Sform_out.out', status='replace', action='write')

163   format(2x,I1,14x,A16)
!*********************************************************************************************************
end subroutine open_output
!*********************************************************************************************************

!*********************************************************************************************************
subroutine close_output
!This subroutine is responsible for closing the output files.
!*********************************************************************************************************

  use file_handles

  implicit none

  integer :: kk, ierr

  !The IOSTAT clause ensures that closing any file which is already closed doesn't crash the program
  do kk=output_first,output_last
    close(UNIT=kk, IOSTAT=ierr)
  end do

!*********************************************************************************************************
end subroutine close_output
!*********************************************************************************************************

!*********************************************************************************************************
subroutine flush_SA_output
!This subroutine is responsible for flushing the output buffers for the
!sensitivity analysis output files.
!*********************************************************************************************************

  use file_handles

  implicit none

  integer :: kk, ierr

  !The IOSTAT clause ensures that flushing any nonexistent file doesn't crash the program
  do kk=SAout_first,SAout_last
    flush(UNIT=kk, IOSTAT=ierr)
  end do

!*********************************************************************************************************
end subroutine flush_SA_output
!*********************************************************************************************************

!*********************************************************************************************************
subroutine write_genout
!This subroutine is responsible for writing some general informational output
!(e.g., version number, problem sizes, options used, etc.)
!*********************************************************************************************************

  use global_all; use global_tube1; use update_tube2; use file_handles; use globalSA

  implicit none

  !General indices
  integer :: i, j
  !The following is to store how many transient and steady state points were
  !written to disk
  integer :: results(2)
  !The number of species in the inlet
  integer :: nspec_in
  !For feed condition schedule
  double precision :: TPFSV_local(4)
  !For formatting output
  character(len=80) :: fmt_str
  !For determining which perturbation options were used for GSA
  integer :: SA_opts(3)

  !Open the file
  open(unit=verout, file='OUT.d/general_info.out', status='replace', action='write')

  !Mechanism size
  write(verout,'(1X,A,I0)') 'This output was generated using the Vlachos Group reactor code, revision #', revno
  write(verout,*) ''
  write(verout,'(1X,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
    'Mechanism contains ', kgmax, ' gas, ', ksmax, ' surface, and ', kbmax, &
    ' bulk species made of ', nelm, ' elements with ', nsphases, ' surface site types'
  write(verout,'(1X,A,I0,A,I0,A)') &
    'Mechanism contains ', ngrxn, ' gas and ', nsrxn, ' surface reactions'

  !Reactor type
  write(verout,*) ''
  select case (irxtr)
  case (0)
    write(verout,'(1X,A,I0,A)') 'Reactor type ', irxtr, ' (UHV/Molecular Beam)'
  case (1)
    write(verout,'(1X,A,I0,A)') 'Reactor type ', irxtr, ' (Batch)'
  case (2)
    write(verout,'(1X,A,I0,A)') 'Reactor type ', irxtr, ' (CSTR)'
  case (3)
    write(verout,'(1X,A,I0,A)') 'Reactor type ', irxtr, ' (PFR)'
    write(verout,'(1X,I0,A)') nnodes, ' CSTRs used'
  end select

  !Reactor volume
  write(verout,*) ''
  write(verout,'(1X,A,ES12.5)') 'Reactor size is ', rlen

  !Isothermality
  write(verout,*) ''
  if (iiso==0) then
    write(verout,*) 'Isothermal system'
  else
    if (itpd>0) then
      write(verout,*) 'Temperature programmed desorption experiment'
    else
      write(verout,*) 'Temperature controlled by external heat transfer'
    end if
  end if

  !Type of output written
  results=0
  write(verout,*) ''
  if (.not. ltra) then
    results=(/0,nnodes+1/)  !SS only
    write(verout,*) 'Steady state results only saved'
  else
    if (irxtr/=1 .or. itpd/=1) then !Both transient and SS
      results=(/ncalls+1,nnodes+1/)
      write(verout,*) 'Steady state and transient results saved'
    else  !Transient only
      results=(/ncalls+1,0/)
      write(verout,*) 'Transient results only saved'
    end if
  end if
  write(verout,'(1X,A,I0,A,I0,A)') 'Results saved for ', results(1), &
    ' transient and ', results(2), ' steady state points per simulation point'

  !Schedule of T, P, F, S/V, and compositions
  nspec_in=size(act_all,1)
  write(verout,*) ''
  write(verout,'(1X,A,I0,A,I0,A)') 'Run conditions for ', nruns, &
    ' simulation points and mole fractions/coverages for ', nspec_in, &
    ' non-zero inlet species'
  write(verout,*) 'Flow rates given at STP'
  fmt_str='(1X,A12,1X,A12,1X,A12,1X,A12,'//trim(int2char(nspec_in))//'(1X,A16))'
  write(verout,fmt_str) 'Temp [K]', 'Press [atm]', 'Flow [cm3/s]', 'S/V [1/cm]', knams(iMole_spec)
  fmt_str='(1X,G12.5,1X,G12.5,1X,G12.5,1X,G12.5,'//trim(int2char(nspec_in))//'(1X,G16.6))'
  do i=1,nruns
    TPFSV_local=TPFSV(i,:)
    !Correct flow rate to be at STP -- do this while p is still in dyne/cm^2
    if (.not. lstp) TPFSV_local(3)=TPFSV_local(3)*(273.15/TPFSV_local(1))*(TPFSV_local(2)/patm)
    TPFSV_local(2)=TPFSV_local(2)/patm !Convert dynes/cm2 -> atm
    write(verout,fmt_str) TPFSV_local, act_all(:,i)
  end do

  !Sensitivity analysis stuff

  !Count the number of options used for writing to output later. This
  !array uses octal values similar to the Unix permission system.
  !It assumes there are a maximum of three options per parameter type.
  !The three options are assigned values of 1, 2, and 4. When each one
  !is applied, its associated value is added to the SA_opts element.
  !Element 1 is the species enthalpy, element 2 is the species entropy,
  !and element 3 is the activation energy. Other parameter types only
  !have one option each and don't need to be included here. Although
  !somewhat confusing, it does substantially simplify the logic later.
  SA_opts=0 !Initialize to zero (no options used)
  !Species enthalpies (element 1)
  if (lDFTH_SA .or. isenbrute==4) SA_opts(1)=SA_opts(1)+1   !DFT/DGSA
  if (lGAH_SA) SA_opts(1)=SA_opts(1)+2    !GA
  if (lLSR_SA) SA_opts(1)=SA_opts(1)+4    !LSR
  !Species entropies
  if (lDFTS_SA .or. isenbrute==4) SA_opts(2)=SA_opts(2)+1   !DFT/DGSA
  if (lGAS_SA) SA_opts(2)=SA_opts(2)+2    !GA
  !Activation energies
  if (lDFTEA_SA .or. isenbrute==4) SA_opts(3)=SA_opts(3)+1  !DFT/DGSA
  if (lBEP_SA) SA_opts(3)=SA_opts(3)+2    !BEP

  write(verout,*) ''
  if (lsenbrute) then
    select case (isenbrute)
    case (1)  !FIM-based SA
      write(verout,*) 'Sensitivity analysis type 1 (FIM-based)'
    case (2:5)  !Brute force SA and UQ
      if (isenbrute==2) then
        write(verout,*) 'Sensitivity analysis type 2 (local brute force)'
        if (rel_pert) then
          write(verout,*) 'Relative perturbations used'
        else
          write(verout,*) 'Absolute perturbations used'
        end if
      else if (isenbrute==3) then
        write(verout,*) 'Sensitivity analysis type 3 (global variance based)'
        write(verout,*) nSAruns, ' replicated points used'
      else if (isenbrute==4) then
        write(verout,*) 'Sensitivity analysis type 4 (global derivative based)'
        write(verout,*) nSAruns, ' replicated points used'
      else
        write(verout,*) 'Sensitivity analysis type 5 (uncertainty quantification)'
        write(verout,*) nSAruns, ' replicated points used'
      end if

      !Species parameters
      write(verout,'(1X,I0,A)') count(perturb_param((/globalSA_H,globalSA_S/))), &
        ' species parameter types perturbed'
      if (count(perturb_param((/globalSA_H,globalSA_S/)))>0) then
        if (fix_EA) then
          write(verout,*) 'Activation energies fixed'
        else
          write(verout,*) 'Transition state energies fixed'
        end if
      end if
      if (perturb_param(globalSA_H)) then
        write(verout,'(1X,A,I0,A)') 'Parameter type ', globalSA_H
        if (isenbrute>2) then
          select case (SA_opts(1))
          case (0)
            write(verout,*) '0 species enthalpy perturbation contributions'
          case (1)
            write(verout,*) '1 species enthalpy perturbation contribution (DFT)'
          case (2)
            write(verout,*) '1 species enthalpy perturbation contribution (GA)'
          case (3)
            write(verout,*) '2 species enthalpy perturbation contributions (DFT, GA)'
          case (4)
            write(verout,*) '1 species enthalpy perturbation contribution (LSR)'
          case (5)
            write(verout,*) '2 species enthalpy perturbation contributions (DFT, LSR)'
          case (6)
            write(verout,*) '2 species enthalpy perturbation contributions (GA, LSR)'
          case (7)
            write(verout,*) '3 species enthalpy perturbation contributions (DFT, GA, LSR)'
          end select
        end if
        if (isenbrute==2) then
          if (rel_pert) then
            write(verout,*) 'Gas and surface species enthalpy perturbations [fractional]'
            write(verout,'(1X,G12.5,1X,G12.5)') pert_local_dist(:,globalSA_H)
          else
            write(verout,*) 'Gas and surface species enthalpy perturbations [kcal/mol]'
            write(verout,'(1X,G12.5,1X,G12.5)') pert_local_dist(:,globalSA_H)*Rgas_kcal
          end if
        end if
      end if
      if (perturb_param(globalSA_S)) then
        write(verout,'(1X,A,I0,A)') 'Parameter type ', globalSA_S
        if (isenbrute>2) then
          select case (SA_opts(2))
          case (0)
            write(verout,*) '0 species entropy perturbation contributions'
          case (1)
            write(verout,*) '1 species entropy perturbation contribution (DFT)'
          case (2)
            write(verout,*) '1 species entropy perturbation contribution (GA)'
          case (3)
            write(verout,*) '2 species entropy perturbation contributions (DFT, GA)'
          end select
        end if
        if (isenbrute==2) then
          if (rel_pert) then
            write(verout,*) 'Gas and surface species entropy perturbations [fractional]'
            write(verout,'(1X,G12.5,1X,G12.5)') pert_local_dist(:,globalSA_S)
          else
            write(verout,*) 'Gas and surface species entropy perturbations [cal/mol-K]'
            write(verout,'(1X,G12.5,1X,G12.5)') pert_local_dist(:,globalSA_S)*Rgas_cal
          end if
        end if
      end if
      if (allocated(spec_pert)) then
        write(verout,'(1X,I0,A)') count(spec_pert), ' species perturbed'
        do i=1,size(spec_pert)
           if (spec_pert(i)) write(verout,*) knams(i)
        end do
      else
        write(verout,'(1X,I0,A)') 0, ' species perturbed'
      end if

      !Reaction parameters
      write(verout,'(1X,I0,A)') count(perturb_param((/globalSA_A,globalSA_B,globalSA_EA/))), &
        ' reaction parameter types perturbed'
      if (perturb_param(globalSA_A)) then
        write(verout,'(1X,A,I0,A)') 'Parameter type ', globalSA_A, &
          ': Gas and surface reaction pre-exponentials (multiplicative)'
        if (isenbrute==2) write(verout,'(1X,G12.5,1X,G12.5)') &
          exp(pert_local_dist(:,globalSA_A))
      end if
      if (perturb_param(globalSA_B)) then
        write(verout,'(1X,A,I0,A)') 'Parameter type ', globalSA_B, &
          ': Gas and surface reaction temperature exponent perturbations'
        if (isenbrute==2) write(verout,'(1X,G12.5,1X,G12.5)') &
          pert_local_dist(:,globalSA_B)
      end if
      if (perturb_param(globalSA_EA)) then
        write(verout,'(1X,A,I0,A)') 'Parameter type ', globalSA_EA
        if (isenbrute>2) then
          select case (SA_opts(2))
          case (0)
            write(verout,*) '0 activation energy perturbation contributions'
          case (1)
            write(verout,*) '1 activation energy perturbation contribution (DFT)'
          case (2)
            write(verout,*) '1 activation energy perturbation contribution (BEP)'
          case (3)
            write(verout,*) '2 activation energy perturbation contributions (DFT, BEP)'
          end select
        end if
        if (isenbrute==2) then
          if (rel_pert) then
            write(verout,*) 'Gas and surface reaction activation energy perturbations [fractional]'
            write(verout,'(1X,G12.5,1X,G12.5)') pert_local_dist(:,globalSA_EA)
          else
            write(verout,*) 'Gas and surface reaction activation energy perturbations [kcal/mol]'
            write(verout,'(1X,G12.5,1X,G12.5)') pert_local_dist(:,globalSA_EA)*Rgas_kcal
          end if
        end if
      end if
      if (allocated(rxn_pert)) then
        write(verout,'(1X,I0,A)') count(rxn_pert), ' reactions perturbed'
        do i=1,size(rxn_pert)
          if (rxn_pert(i)) write(verout,*) rxn_str(i)
        end do
      else
        write(verout,'(1X,I0,A)') 0, ' reactions perturbed'
      end if

      !Correlations
      write(verout,'(1X,I0,A)') count(perturb_param((/globalSA_BEP,globalSA_LSR/))), &
        ' correlation parameter types perturbed'
      if (perturb_param(globalSA_BEP)) then
        write(verout,'(1X,A,I0,A)') 'Parameter type ', globalSA_BEP, ': BEP coefficients'
        if (isenbrute==2) then
          write(verout,'(1X,I0,A)') nBEP, ' correlations perturbed'
          if (rel_pert) then
            write(verout,'(1X,A12,1X,A12)') 'Slope', 'Int [fractional]'
          else
            write(verout,'(1X,A12,1X,A12)') 'Slope', 'Int [kcal/mol]'
          end if
          do j=1,nBEP
            if (rel_pert) then
              write(verout,'(1X,G12.5,1X,G12.5)') BEP_coeff_pert_dist(j,1), &
                 BEP_coeff_pert_dist(j,3)
            else
              write(verout,'(1X,G12.5,1X,G12.5)') BEP_coeff_pert_dist(j,1), &
                 BEP_coeff_pert_dist(j,3)*Rgas_kcal
            end if
          end do
        end if
      end if
      if (perturb_param(globalSA_LSR)) then
        write(verout,'(1X,A,I0,A)') 'Parameter type ', globalSA_LSR, ': LSR coefficients'
        if (isenbrute==2) then
          write(verout,'(1X,I0,A)') nscale, ' correlations perturbed'
          if (rel_pert) then
            write(verout,'(1X,A12,1X,A12)') 'Slope', 'Int [fractional]'
          else
            write(verout,'(1X,A12,1X,A12)') 'Slope', 'Int [kcal/mol]'
          end if
          do j=1,nscale
            if (rel_pert) then
              write(verout,'(1X,G12.5,1X,G12.5)') LSR_coeff_pert_dist(j,1), &
                 LSR_coeff_pert_dist(j,3)
            else
              write(verout,'(1X,G12.5,1X,G12.5)') LSR_coeff_pert_dist(j,1), &
                 LSR_coeff_pert_dist(j,3)*Rgas_kcal
            end if
          end do
        end if
      end if
    end select
  else
    write(verout,*) 'Sensitivity analysis type 0 (off)'
  end if

  !Binding energy maps
  write(verout,*) ''
  if (iScale==0) then
    write(verout,*) 'Calculations performed at 0 mapped binding energies per simulation point'
  else
    write(verout,'(1X,A,I0,A)') 'Calculations performed at ', &
      nBE_coords,' mapped binding energies per simulation point'
    write(verout,'(1X,I0,A)') natoms_scale, ' atomic binding energy coordinates [kcal/mol]'
    fmt_str='(1X,A8,'//trim(int2char(natoms_scale))//'(1X,A16))'
    write(verout,fmt_str) 'BE Point', knams(nint(scale_targ(:,1)))
    fmt_str='(1X,I8,'//trim(int2char(natoms_scale))//'(1X,G16.6))'
    do i=1,nBE_coords
      write(verout,fmt_str) i, BE_coords(:,i)*Rgas_kcal
    end do
  end if

  !Finish up
  close(verout)

!*********************************************************************************************************
end subroutine write_genout
!*********************************************************************************************************

!*********************************************************************************************************
subroutine write_rpa(z,T)
!This subroutine is responsible for writing the RPA files.
!*********************************************************************************************************

  use global_all; use global_tube1; use file_handles; use update_tube2

  implicit none

  double precision, intent(in) :: z
  double precision, intent(in) :: T
  integer :: file_LU(2)
  double precision :: rpa_local(kmax,nrxn), rpa_vis(kmax,nrxn), pei(nrxn)
  double precision :: rop_rpa(nrxn), ropfwd_rpa(nrxn)
  integer :: irxns, irpa, kspec, kk
  character(len=80) :: fmt_str
  double precision  ::  PreEx(nrxn),Beta(nrxn),Ea(nrxn),EQKC(nrxn)
  double precision :: EaOverRT(nrxn)
  integer :: ISFLAG(nsrxn)

  !Convert input rates to mol/s
  rop_rpa(1:ngrxn)=rop(1:ngrxn)*z
  rop_rpa(ngrxn+1:ngrxn+nsrxn)=rop(ngrxn+1:ngrxn+nsrxn)*abyv*z
  ropfwd_rpa(1:ngrxn)=ropfwd(1:ngrxn)*z
  ropfwd_rpa(ngrxn+1:ngrxn+nsrxn)=ropfwd(ngrxn+1:ngrxn+nsrxn)*abyv*z
  if (irxtr==3 .and. verbose_rpa) then
    rop_rpa=rop_rpa*nnodes
    ropfwd_rpa=ropfwd_rpa*nnodes
  end if
  !Calculate the rpa matrices
  call rpa_calc(rop_rpa,ropfwd_rpa,rpa_local,rpa_vis,pei)

  !Extract the actual activation energies and equilibrium constants
  call ckabe(ickwrk,rckwrk,PreEx(1:ngrxn),Beta(1:ngrxn),Ea(1:ngrxn))
  call skabe(iskwrk,rskwrk,PreEx(ngrxn+1:nrxn),Beta(ngrxn+1:nrxn),Ea(ngrxn+1:nrxn),ISFLAG)
  call ckeqxp(P,T,act(1:kgmax),ickwrk,rckwrk,EQKC(1:ngrxn))
  call skeq(P,T,act,sden,iskwrk,rskwrk,EQKC(ngrxn+1:nrxn),Tref_beta)
  Ea=Ea*Rgas_kcal !Convert to kcal/mol

  !Choose the output files
  if (mrpa==1) then
    file_LU(1)=ifiletuberpa
    file_LU(2)=irpavisout
  else
    file_LU(1)=ifiletrpa
  end if

  !Write header information
  write(file_LU(1),008)' ********* Experimental Condition # ',i_tube,' *********'
  write(file_LU(1),'(A,F7.2)') 'Temperature [K]: ',T
  if (iScale>0) call write_BE_coords(file_LU(1))
  if (irxtr==3 .and. verbose_rpa) write(file_LU(1),'(A,ES12.5)') 'Position [cm3]:  ', z
  write(file_LU(1),*)'Rates (mol/s) and Partial Equilibrium Analysis:'
  write(file_LU(1),*)'  Rxn #    Net Rate       Fwd Rate       PE Ratio       EA (kcal/mol)  Equil Const'
  if (mrpa==1) then
    write(file_LU(2),'(A,F7.2,A)') 'Temperature =',T,' K; Rates in mol/s'
    if (iScale>0) call write_BE_coords(file_LU(2))
    write(file_LU(2),*)'-----------------------------------------------------'
    write(file_LU(2),'(30X,A)') 'Percentages of production & consumption'
    fmt_str='(A,'//trim(int2char(kgmax+ksmax))//'(A20))'
    write(file_LU(2),fmt_str) 'Rxn #    Net Rate      PEI    ', knams(1:kgmax+ksmax)
  end if

  !Write body information

  !Rates as an exhaustive list and rpa_vis
  do irxns=1,nrxn
    write(file_LU(1),9000)irxns,rop_rpa(irxns),ropfwd_rpa(irxns),pei(irxns),Ea(irxns),&
      EQKC(irxns),"     "//trim(rxn_str(irxns))
    if (mrpa==1) then
      fmt_str='(1X,I3,4X,ES11.4,3X,F5.2,3X,'//trim(int2char(kgmax+ksmax))//'(ES20.4))'
      write(file_LU(2),fmt_str) irxns, rop_rpa(irxns),pei(irxns),&
           rpa_vis(1:kgmax+ksmax,irxns)
    end if
  end do

  !Rates by species
  write(file_LU(1),*)'----------Reaction Path Analysis----------'
  d_07: do irpa=1,kmax-kbmax
    kspec=irpa
    write(file_LU(1),*) '------------------------------------------'
    write(file_LU(1),*) 'going to species number: ',irpa, ' and name: ',knams(irpa)
    write(file_LU(1),*) '  Rxn #   % Prod/Cons    Net Rate        Fwd Rate       PE Ratio'

    d_09: do irxns=1,nsrxn
      if (rpa_local(irpa,irxns)>0.0.or.rpa_local(irpa,irxns)<0.0) then
        write(file_LU(1),9010) irxns,rpa_local(irpa,irxns),rop_rpa(irxns),ropfwd_rpa(irxns),&
          pei(irxns), "     "//trim(rxn_str(irxns))
      endif
    end do d_09
  end do d_07

008   format(A,I2,A)
9000  format(1I8,5ES15.5,85A)
9010  format(1I8,4ES15.5,85A)

!*********************************************************************************************************
end subroutine write_rpa
!*********************************************************************************************************

!*********************************************************************************************************
subroutine write_stoich
!This subroutine writes the stoichiometric matrix to disk.
!*********************************************************************************************************

  use global_all; use file_handles

  implicit none

  integer :: i
  character(len=80) :: fmt_str
  double precision :: atomic_mass(nelm)

!-------------------------------------------------------------------------------------------------------

  !Initialize the format string
  fmt_str='('//trim(int2char(nrxn))//'I4)'

  !Write the stoichiometric matrix
  do i=1,kmax
    write(iStoichout,fmt_str) stoich_matrix(i,:)
  end do

  !Write the surface phase densities and phase names
  do i=2,nsphases+1
    write(ielemcomp,'(1X,A,1X,ES12.5)') PhaseNames(i), sden(i)
  end do
  write(ielemcomp,*) ''

  !Get the atomic masses
  call ckawt(ickwrk,rckwrk,atomic_mass)

  !Write the elemental compositions and masses
  fmt_str='('//trim(int2char(nelm+1))//'A16)'
  write(ielemcomp,fmt_str) 'Species         ', enams(:)
  fmt_str='(A16,'//trim(int2char(nelm))//'(F15.6,1X))'
  write(ielemcomp,fmt_str) 'Atomic Weights  ', atomic_mass(:)
  fmt_str='(A16,'//trim(int2char(nelm))//'(I15,1X))'
  do i=1,kmax
    write(ielemcomp,fmt_str) knams(i), elem_comp(:,i)
  end do

!*********************************************************************************************************
end subroutine write_stoich
!*********************************************************************************************************

!*********************************************************************************************************
subroutine rpa_calc(rop,ropfwd,rpa,rpa_vis,pei)
!This subroutine performs the calculations necessary to write the reaction
!path analysis.
!*********************************************************************************************************

  use global_all; use file_handles

  implicit none

  double precision, intent(in) :: rop(nrxn), ropfwd(nrxn)
  double precision, intent(out) :: rpa(kmax,nrxn), rpa_vis(kmax,nrxn)
  double precision, intent(out) :: pei(nrxn)
  double precision :: cik(nrxn)
  integer :: irpa, kspec, irxns
  double precision :: NaN

!-------------------------------------------------------------------------------------------------------

  !Calculate the partial equilibrium indices
  pei=ropfwd/(-rop+2*ropfwd)

  !Calculate the rates of production of each species as a function of the
  !reaction. SPREAD(A,b,c) will tile array A along dimension b a total of c
  !times (similar to Matlab's repmat() function).
  rpa=spread(rop,1,kmax)  !rop is length nrxn to yield an array of dimensions kmax * nrxn
  rpa=rpa*stoich_matrix  !Element-by-element multiplication

  !Calculate the % production/consumption
  do irpa=1,kmax
    cik=rpa(irpa,:)
    where (cik>=0.0 .and. sum(cik,cik>=0.0)/=0)
      rpa(irpa,:)=100*cik/sum(cik,cik>=0.0)
    elsewhere (cik<0 .and. sum(cik,cik<0.0)/=0)
      rpa(irpa,:)=-100*cik/sum(cik,cik<0.0)
    end where
  end do

  NaN=-1.0
  NaN=sqrt(NaN)
  where (stoich_matrix/=0)
    rpa_vis=rpa
  elsewhere
    rpa_vis=NaN  !This will put NaNs everywhere
  end where

!*********************************************************************************************************
end subroutine rpa_calc
!*********************************************************************************************************

!*********************************************************************************************************
subroutine compare_rates(grate,srate,abyv,ifilerates)
!*********************************************************************************************************

  use global_all

  implicit none

  integer, intent(in) :: ifilerates
  double precision, intent(in) :: grate(kgmax), srate(kmax), abyv
  double precision :: rate_tot(kgmax), contri_gas(kgmax), contri_surf(kgmax)
  integer :: ii

!-------------------------------------------------------------------------------------------------------
  forall (ii=1:kgmax) rate_tot(ii)=abs(grate(ii))+abyv*abs(srate(ii))
  contri_gas=0.0d0
  contri_surf=0.0d0
  do ii=1,kgmax
    if (rate_tot(ii)==0.0d0) then
      write(ifilerates,*)'Net rate zero for species: ',knams(ii),'skipped'
      cycle
    else
      contri_gas(ii)=100.0d0*grate(ii)/rate_tot(ii)
      contri_surf(ii)=100.0d0*abyv*srate(ii)/rate_tot(ii)
    endif
    write(ifilerates,150)knams(ii),contri_gas(ii),contri_surf(ii),grate(ii),abyv*srate(ii)
  end do
150 format(A13,2(1x,f9.2),2(1x,es10.3))

!*********************************************************************************************************
end subroutine compare_rates
!*********************************************************************************************************

!*********************************************************************************************************
subroutine write_model_gas_surf(z,gdot_save,sdot_save,write_option)
!This subroutine compares gas and surface contributions.
!*********************************************************************************************************

  use global_all; use global_tube1; use update_tube2; use file_handles

  implicit none
  double precision, intent(in) :: z, gdot_save(kgmax,nnodes), sdot_save(kgmax,nnodes)
  integer, intent(in) :: write_option
  double precision :: gdot_integral(kgmax), sdot_integral(kgmax)

  if (write_option==1) then
    write(ituberateout,'(A,I3)') 'Experimental condition #', i_tube
    write(ituberateout,*) 'Volume v [cm3] =',z,'-------------------'
    write(ituberateout,*) 'Species,     contri_gas, contri_surf, gdot,  abyv*sdot'
  end if
  if (write_option<3) call compare_rates(gdot,sdot,abyv,ituberateout)
  if (write_option==3) then
    write(ituberateout,*) '----------------- Integrated rates and contributions -----------------'
    write(ituberateout,*) 'Species,    contri_gas, contri_surf, gdot_integral, abyv*sdot_integral'
    call integrate_data(kgmax,nnodes,gdot_save,rlen,gdot_integral)
    call integrate_data(kgmax,nnodes,sdot_save,rlen,sdot_integral)
    call compare_rates(gdot_integral,sdot_integral,abyv,ituberateout)
  end if

!*********************************************************************************************************
end subroutine write_model_gas_surf
!*********************************************************************************************************

!*********************************************************************************************************
subroutine integrate_data(dim1,dim2,data_int,span_int,integral_value)
!*********************************************************************************************************

  use global_all

  implicit none

  integer, intent(in) :: dim1, dim2
  double precision, intent(in) :: data_int(dim1,dim2)
  double precision, intent(out) :: integral_value(dim1)
  double precision :: span_int, h_int
  integer :: ii

!-------------------------------------------------------------------------------------------------------
  h_int=span_int/float(dim2)
  forall (ii=1:dim1) integral_value(ii)=(data_int(ii,1)+data_int(ii,dim2))*(h_int/2.0d0)+sum(data_int(ii,2:(dim2-1))*h_int)

!*********************************************************************************************************
end subroutine integrate_data
!*********************************************************************************************************

!*********************************************************************************************************
subroutine convfwd(w)
!This subroutine calculates solution vector w, given x,y and coverages.
!*********************************************************************************************************

  use global_all; use global_tube1; use update_tube2

  implicit none

  double precision, intent(out) :: w(neqns)
  integer :: i

!-------------------------------------------------------------------------------------------------------

  x(1:kgmax)=act(1:kgmax)
  xin(1:kgmax)=actin(1:kgmax)
  call ckxty(x,ickwrk,rckwrk,y)
!  call ckytx(y,ickwrk,rckwrk,x)
  call ckxty(xin,ickwrk,rckwrk,yin)
!  call ckytx(yin,ickwrk,rckwrk,xin)
!  act(1:kgmax)=x(1:kgmax)
!  actin(1:kgmax)=xin(1:kgmax)
  w(1:kgmax)=y(1:kgmax)
  forall (i=kgmax+1:kmax-kbmax) w(i)=act(i)
  act(kmax-kbmax+1:kmax)=1.0
  actin(kmax-kbmax+1:kmax)=1.0

!*********************************************************************************************************
end subroutine convfwd
!*********************************************************************************************************

!*********************************************************************************************************
subroutine convbkd(w)
!This subroutine calculates x,y and coverages, given solution vector w.
!*********************************************************************************************************

  use global_all; use global_tube1; use update_tube2

  implicit none

  double precision, intent(in) :: w(neqns)
  integer :: i

!-------------------------------------------------------------------------------------------------------

  y(1:kgmax)=w(1:kgmax)
  xin(1:kgmax)=actin(1:kgmax)
  call ckytx(y,ickwrk,rckwrk,x)
!  call ckxty(x,ickwrk,rckwrk,y)
  call ckxty(xin,ickwrk,rckwrk,yin)
!  call ckytx(yin,ickwrk,rckwrk,xin)
  act(1:kgmax)=x(1:kgmax)
!  actin(1:kgmax)=xin(1:kgmax)
  forall (i=kgmax+1:kmax-kbmax) act(i)=w(i)
  act(kmax-kbmax+1:kmax)=1.0
  actin(kmax-kbmax+1:kmax)=1.0

!*********************************************************************************************************
end subroutine convbkd
!*********************************************************************************************************

!*********************************************************************************************************
subroutine check(neqns,w,t,jflag)
! This subroutine checks for
! -Negative mole/mass fractions or coverages
! -Mole/Mass fractions or coverages greater than 1.0
! -Normalizes such variables
! -Abnormal temperature (which could crash CHEMKIN/SURFACE-CHEMKIN)
!*********************************************************************************************************

  use global_all

  implicit none

  integer, intent(in) :: neqns
  integer, intent(out) :: jflag
  double precision :: w(neqns), t
  double precision :: sumgas, sumsurf
  integer :: i

!-------------------------------------------------------------------------------------------------------

  jflag=0
  sumgas=sum(abs(w(1:kgmax)))
  sumsurf=sum(abs(w(kgmax+1:kmax-kbmax)))
  if (abs(sumgas)<1.e-10) write(*,*)'sum of gas mass fractions close to zero'
  if (abs(sumsurf)<1.e-10) write(*,*)'sum of surface coverages close to zero'
  forall (i=1:kgmax, sumgas/=0.0) w(i)=abs(w(i))/sumgas
  forall (i=kgmax+1:kmax-kbmax, sumsurf/=0.0) w(i)=abs(w(i))/sumsurf
  sumgas=sum(abs(w(1:kgmax)))
  sumsurf=sum(abs(w(kgmax+1:kmax-kbmax)))
  if (abs(sumgas-1.0)>1.e-5.or.abs(sumsurf-1.0)>1.e-5) print*, 'wrong normalization','gas',sumgas,'surf',sumsurf
  if (t<50.0.or.t>3000.0) then
  write(*,*)'problem in temperature'
  write(*,*)'t = ',t
  write(*,*)'stopping flag activated'
  jflag=1
  return
  endif

!*********************************************************************************************************
end subroutine check
!*********************************************************************************************************

!*********************************************************************************************************
subroutine NSC_matrix_write(scale_iter)
!This subroutine writes the normalized sensitivity coefficient matrix to a file.
!The NSC matrix can be for any of the SA types.
!*********************************************************************************************************

  use global_all; use global_tube1; use update_tube2; use file_handles; use globalSA

  implicit none

  integer, intent(in) :: scale_iter
  integer :: i, j, jlo, jhi
  integer :: par_number, par_type
  character (len=81) :: buffer_long
  character (len=17) :: buffer_short
  integer :: buffer_long_max
  character(len=80) :: fmt_str

!---------------------------------------------------------------------------------------------------------

  write(ifiletubesen,008)' ********* Experimental Condition # ',i_tube,' *********'
  if (iScale>0) call write_BE_coords(ifiletubesen)

  select case (isenbrute)

    case (1)  !FIM-based approximate SA

      !Find the length of the longest reaction string
      buffer_long_max=0
      do i=1,nrxn
        buffer_long=rxn_str(i)
        buffer_long_max=max(buffer_long_max,len_trim(buffer_long))
      end do

      write(ifiletubesen,*) 'FIM-Based Sensitivity Coefficients'
      write(ifiletubesen,*) 'Reaction String'//repeat(' ',buffer_long_max-14)//'FIM Element Fwd/Rev'

      fmt_str='(1X,A,1X,ES15.5E3,1X,ES15.5E3,1X)'
      do i=1,nrxn
        buffer_long=rxn_str(i)
        write(ifiletubesen,fmt_str) buffer_long(1:buffer_long_max+1), FIM_SA_matrix(i,:)
      end do

    case (2, 3)  !LSA, GSA

      buffer_long='Species/Reaction/Correlation'
      buffer_long_max=len_trim(buffer_long)

      !Find the length of the longest string
      do i=1,max_SA_param_types
        if (.not. perturb_param(i)) cycle
        select case (i)
        case (globalSA_H, globalSA_S)
          do j=1,kgmax+ksmax
            buffer_long=knams(j)
            buffer_long_max=max(buffer_long_max,len_trim(buffer_long))
          end do
        case (globalSA_A, globalSA_B, globalSA_EA)
          do j=1,nrxn
            buffer_long=rxn_str(j)
            buffer_long_max=max(buffer_long_max,len_trim(buffer_long))
          end do
        case default
          buffer_long_max=max(buffer_long_max,11)
        end select
      end do

      !Check for sensitivity analysis type and assign loop limits
      jlo=1 !Always start at 1
      if (isenbrute==2) then  !LSA, only one type of metric
        jhi=1
      else  !GSA, three metrics
        jhi=3
      end if

      !Now we write the sensitivity coefficients/indices based on (1) the type
      !of SA and (2) which coefficients have already been written
      do j=jlo,jhi
        fmt_str='(1X,A'//trim(int2char(buffer_long_max))//',1X,A17,1X'//trim(int2char(kgmax))//'A16)'
        if (isenbrute==2) then
          write(ifiletubesen,*) 'Local Sensitivity Analysis: Semi-Normalized Sensitivity Coefficients (x*dr/dx)'
        else
          select case (j)
          case (1)
            write(ifiletubesen,*) 'Global Sensitivity Analysis: Unnormalized Elementary Effect Indices'
          case (2)
            write(ifiletubesen,*) 'Global Sensitivity Analysis: First Order Effect Indices'
          case (3)
            write(ifiletubesen,*) 'Global Sensitivity Analysis: Total Effect Indices'
          end select
        end if
        write(ifiletubesen,fmt_str) 'Species/Reaction/Correlation'//repeat(' ',buffer_long_max-27),'Type'//repeat(' ',14),(knams(i), i=1,kgmax)

        fmt_str='(1X,A'//trim(int2char(buffer_long_max))//',1X,A17,1X'//trim(int2char(kgmax))//'(ES15.5E3,1X))'

        do i=1,nparams_total
          par_type=SApert_data(i,1)%par_type
          par_number=SApert_data(i,1)%par_number
          select case (abs(par_type))
          case (globalSA_H)
            buffer_long=knams(par_number)
            buffer_short='Enthalpy'
          case (globalSA_S)
            buffer_long=knams(par_number)
            buffer_short='Entropy'
          case (globalSA_A)
            buffer_long=rxn_str(par_number)
            buffer_short='Pre-exponential'
          case (globalSA_B)
            buffer_long=rxn_str(par_number)
            buffer_short='Beta'
          case (globalSA_EA)
            buffer_long=rxn_str(par_number)
            buffer_short='Activation Energy'
          case (globalSA_BEP)
            if (par_type>0) then
              buffer_long='BEP Slope'
            else
              buffer_long='BEP Intercept'
            end if
            write(buffer_short,'(I17)') par_number
          case (globalSA_LSR)
            if (par_type>0) then
              buffer_long='LSR Slope'
            else
              cycle !Put this cycle in because we don't currently perturb the LSR intercept
              !This cycle should be removed if this ever changes
              buffer_long='LSR Intercept'
            end if
            write(buffer_short,'(I17)') par_number
          end select
          if (isenbrute==2) then
            write(ifiletubesen,fmt_str) buffer_long, buffer_short, NSC_matrix(i,:)
          else
            select case (j)
            case (1)
              write(ifiletubesen,fmt_str) buffer_long, buffer_short, elem_effect_matrix(i,:,scale_iter)
            case (2)
              write(ifiletubesen,fmt_str) buffer_long, buffer_short, first_order_SI_matrix(i,:,scale_iter)
            case (3)
              write(ifiletubesen,fmt_str) buffer_long, buffer_short, total_SI_matrix(i,:,scale_iter)
            end select
          end if
        end do
        if (isenbrute==3) write(ifiletubesen,*) repeat('-',80)
      end do

  end select

008  format(A,I2,A)

!*********************************************************************************************************
end subroutine NSC_matrix_write
!*********************************************************************************************************

!!*********************************************************************************************************
!subroutine write_matrix(A,M,N)
!!*********************************************************************************************************
!
!  implicit none
!
!  double precision, intent(in) :: A(M,N)
!  integer :: i, j
!
!!-------------------------------------------------------------------------------------------------------
!  write(*,*)'-------------------------------------------------------------------------'
!  do i = lbound(A,1),ubound(A,1)
!  write(*,*)'ROW NUMBER',i
!  write(*,*) (A(i,j),j=lbound(A,2),ubound(A,2))
!  end do
!  write(*,*)'-------------------------------------------------------------------------'
!!*********************************************************************************************************
!end subroutine write_matrix
!!*********************************************************************************************************

!*********************************************************************************************************
subroutine write_conversion(i_tube,T,X,rate)
!This subroutine writes the conversion data to disk.
!*********************************************************************************************************
  use file_handles

  implicit none

  integer, intent(in) :: i_tube
  double precision, intent(in) :: T, X, rate  !Temperature, Conversion, Rate

  write(itubeconvout,100) i_tube, T, X, rate

100 format(I3,4X,F8.2,4X,F7.2,7X,G10.3)

!*********************************************************************************************************
end subroutine write_conversion
!*********************************************************************************************************

!*********************************************************************************************************
subroutine write_thermokin(T)
!This subroutine will write the thermodynamic and kinetic parameter values
!present at the start of the simulation. It is called once for every call of
!tube_sub.
!*********************************************************************************************************

  use global_all; use global_tube1; use update_tube2; use file_handles

  implicit none

  double precision, intent(in) :: T
  double precision  ::  HORT(kmax),SOR(kmax),PreEx(nrxn),Beta(nrxn),Ea(nrxn),EQKC(nrxn)
  double precision      ::      DELS(nrxn),DELH(nrxn),DELG(nrxn)
  integer :: Rxn1, Spec1, ig
  integer :: ISFLAG(nsrxn)

  call sksor(T,iskwrk,rskwrk,SOR)
  call skhort(T,iskwrk,rskwrk,HORT)
  call ckabe(ickwrk,rckwrk,PreEx(1:ngrxn),Beta(1:ngrxn),Ea(1:ngrxn))
  call skabe(iskwrk,rskwrk,PreEx(ngrxn+1:nrxn),Beta(ngrxn+1:nrxn),Ea(ngrxn+1:nrxn),ISFLAG)
  call ckeqxp(P,T,act(1:kgmax),ickwrk,rckwrk,EQKC(1:ngrxn))
  call skeq(P,T,act,sden,iskwrk,rskwrk,EQKC(ngrxn+1:nrxn),Tref_beta)

  !Calculate reaction enthalpies, entropies, and free energies
  call rxnenth((/1,nrxn/),HORT,delh)  !Enthalpy
  call rxnenth((/1,nrxn/),SOR,dels)  !Entropy
  delg=delh-dels !Free energy

  Ea = Ea/T
  write(iDHRTout,158) 'Dimensionless Enthalpies of Reaction H/RT at T = ', T, ' K. Run #',i_tube
  if (iScale>0) call write_BE_coords(iDHRTout)
  write(iDSRout,158) 'Dimensionless Entropies of Reaction S/R at T = ', T, ' K. Run #',i_tube
  if (iScale>0) call write_BE_coords(iDSRout)
  write(iDGRTout,158) 'Dimensionless Gibbs Energies of Reaction G/RT at T = ', T, ' K. Run #',i_tube
  if (iScale>0) call write_BE_coords(iDGRTout)
  write(ikout,158) 'Pre-exponentials at T = ', T, ' K. Run #',i_tube
  if (iScale>0) call write_BE_coords(ikout)
  write(ibetaout,158) 'Beta Values at T = ', T, ' K. Run #',i_tube
  if (iScale>0) call write_BE_coords(ibetaout)
  write(iEARTout,158) 'Dimensionless Activation Energies EA/RT at T = ', T, ' K. Run #',i_tube
  if (iScale>0) call write_BE_coords(iEARTout)
  write(iEQKCout,158) 'Equilibrium Constants at T = ', T, ' K. Run #',i_tube
  if (iScale>0) call write_BE_coords(iEQKCout)
  do Rxn1=1,nrxn
    write(iDHRTout,157)DELH(Rxn1), rxn_str(Rxn1)
    write(iDSRout,157)DELS(Rxn1), rxn_str(Rxn1)
    write(iDGRTout,157)DELG(Rxn1), rxn_str(Rxn1)
    write(ikout,157)PreEx(Rxn1), rxn_str(Rxn1)
    write(ibetaout,157)Beta(Rxn1), rxn_str(Rxn1)
    write(iEARTout,157)Ea(Rxn1), rxn_str(Rxn1)
    write(iEQKCout,157)EQKC(Rxn1), rxn_str(Rxn1)
  end do

  write(iHRTfout,158) 'Dimensionless Enthalpies of Formation H/RT at T = ', T, ' K. Run #',i_tube
  if (iScale>0) call write_BE_coords(iHRTfout)
  write(iHRTfout,157) (HORT(ig), KNAMS(ig), ig=1,kmax)
  write(iSRfout,158) 'Dimensionless Entropies of Formation S/R at T = ', T, ' K. Run #',i_tube
  if (iScale>0) call write_BE_coords(iSRfout)
  write(iSRfout,157) (SOR(ig), KNAMS(ig), ig=1,kmax)

157 format (1X,ES12.5,1X,A)
158 format (1X,A,F6.2,A,I3)
!*********************************************************************************************************
end subroutine write_thermokin
!*********************************************************************************************************

!*********************************************************************************************************
subroutine write_model_output(tt,w,wtprime,header_write,node)
!This subroutine will write to the output files at each time/distance step.
!*********************************************************************************************************

  use global_all; use global_tube1; use update_tube2; use file_handles

  implicit none

  double precision, intent(in) :: tt, w(neqns)
  integer, intent(in) :: node
  logical, intent(in) :: header_write
  integer :: kk, name_len, specID, pert
  character(len=80) :: fmt_str(2)
  character(len=17) :: atom_names(natoms_scale)
  double precision :: pressure
  double precision :: wtprime(neqns), wta

  !Create some format statements; the additional E3 in the first statement
  !simply makes the field width of the exponent to be 3 characters wide so
  !that very large/small numbers are properly represented.
  fmt_str(1)='('//trim(int2char(kmax+3))//'(ES18.8E3,1X))'
  fmt_str(2)='(A10,9X,'//trim(int2char(kmax+2))//'(A19))'

  call convbkd(w)
  pressure=p/patm !Convert to atmospheres
  call ckwyp(p,w(indext),y,ickwrk,rckwrk,gdot)
  call ckqyp(p,w(indext),y,ickwrk,rckwrk,rop(1:ngrxn))
  call ckqypf(p,w(indext),y,ickwrk,rckwrk,ropfwd(1:ngrxn))
  call skrat(p,w(indext),act,sden,iskwrk,rskwrk,sdot,sitdot,Tref_beta)
  call skrop(p,w(indext),act,sden,iskwrk,rskwrk,rop(ngrxn+1:ngrxn+nsrxn),Tref_beta)
  call skropfwd(p,w(indext),act,sden,iskwrk,rskwrk,ropfwd(ngrxn+1:ngrxn+nsrxn),Tref_beta)
  call ckmmwy(y,ickwrk,rckwrk,wta)

  if (node>0) then  !Transient output
    if (header_write) then
      if (iScale>0 .and. node==1) call write_BE_coords(igasmoletraout)
      write(igasmoletraout,'(1X,2(A,I3))') 'Experimental condition #', i_tube, '  Reactor #', node
      write(igasmoletraout,fmt_str(2)) 'Time [s]  ',knams(1:kgmax), 'Temp. [K]', 'Pressure [atm]'
    end if
    write(igasmoletraout,fmt_str(1))tt,(act(kk),kk=1,kgmax),w(indext),pressure
    if (header_write) then
      if (iScale>0 .and. node==1) call write_BE_coords(igasmasstraout)
      write(igasmasstraout,'(1X,2(A,I3))') 'Experimental condition #', i_tube, '  Reactor #', node
      write(igasmasstraout,fmt_str(2)) 'Time [s]  ',knams(1:kgmax), 'Temp. [K]', 'Pressure [atm]'
    end if
    write(igasmasstraout,fmt_str(1))tt,(y(kk),kk=1,kgmax),w(indext),pressure
    if (surf_chem) then
      if (header_write) then
        if (iScale>0 .and. node==1) call write_BE_coords(icovtraout)
        write(icovtraout,'(1X,2(A,I3))') 'Experimental condition #', i_tube, '  Reactor #', node
        write(icovtraout,fmt_str(2)) 'Time [s]  ',knams(kgmax+1:kmax-kbmax),'Temp. [K]'
        if (iScale>0 .and. node==1) call write_BE_coords(isdottraout)
        write(isdottraout,'(1X,2(A,I3))') 'Experimental condition #', i_tube, '  Reactor #', node
        write(isdottraout,fmt_str(2)) 'Time [s]  ',knams(1:kmax-kbmax),'Temp. [K]'
      end if
      write(icovtraout,fmt_str(1))tt,(w(kk),kk=kgmax+1,kmax-kbmax),w(indext)
      write(isdottraout,fmt_str(1))tt,(sdot(kk),kk=1,kmax-kbmax),w(indext)
    end if
    if (gas_chem) then
      if (header_write) then
        if (iScale>0 .and. node==1) call write_BE_coords(igdottraout)
        write(igdottraout,'(1X,2(A,I3))') 'Experimental condition #', i_tube, '  Reactor #', node
        write(igdottraout,fmt_str(2)) 'Time [s]  ',knams(1:kgmax),'Temp. [K]'
      end if
      write(igdottraout,fmt_str(1))tt,(gdot(kk),kk=1,kgmax),w(indext)
    end if
    if (header_write) then !Mass balance
      write(imassbaltraout,'(1X,2(A,I3))') 'Experimental condition #', i_tube, '  Reactor #', node
      write(imassbaltraout,fmt_str(2)) 'Time [s]  ', 'Avg MM [g/mol]', 'rho [g/cm3]', 'drho/dt [g/cm3s]','dm/dt [g/s]','mdot_in [g/s]','mdot_out [g/s]','m_gen_surf [g/s]'
    end if
    write(imassbaltraout,fmt_str(1))tt,wta,w(indexrho),wtprime(indexrho),wtprime(indexrho)*rlen,fluxmass,fluxmass-wtprime(indexrho)*rlen+abyv*rlen*sum(sdot(1:kgmax)*wt),abyv*rlen*sum(sdot(1:kgmax)*wt)
  else if (node==0) then !Steady state output
    if (header_write) then
      if (iScale>0) call write_BE_coords(igasmolessout)
      write(igasmolessout,'(1X,A,I3)') 'Experimental condition #', i_tube
      write(igasmolessout,fmt_str(2)) 'Vol [cm3] ',knams(1:kgmax), 'Temp. [K]', 'Pressure [atm]'
    end if
    write(igasmolessout,fmt_str(1))tt,(act(kk),kk=1,kgmax),w(indext),pressure
    if (header_write) then
      if (iScale>0) call write_BE_coords(igasmassssout)
      write(igasmassssout,'(1X,A,I3)') 'Experimental condition #', i_tube
      write(igasmassssout,fmt_str(2)) 'Vol [cm3] ',knams(1:kgmax), 'Temp. [K]', 'Pressure [atm]'
    end if
    write(igasmassssout,fmt_str(1))tt,(y(kk),kk=1,kgmax),w(indext),pressure
    if (surf_chem) then
      if (header_write) then
        if (iScale>0) call write_BE_coords(icovssout)
        write(icovssout,'(1X,A,I3)') 'Experimental condition #', i_tube
        write(icovssout,fmt_str(2)) 'Vol [cm3] ',knams(kgmax+1:kmax-kbmax),'Temp. [K]'
        if (iScale>0) call write_BE_coords(isdotssout)
        write(isdotssout,'(1X,A,I3)') 'Experimental condition #', i_tube
        write(isdotssout,fmt_str(2)) 'Vol [cm3] ',knams(1:kmax-kbmax),'Temp. [K]'
      end if
      write(icovssout,fmt_str(1))tt,(w(kk),kk=kgmax+1,kmax-kbmax),w(indext)
      write(isdotssout,fmt_str(1))tt,(sdot(kk),kk=1,kmax-kbmax),w(indext)
    end if
    if (gas_chem) then
      if (header_write) then
        if (iScale>0) call write_BE_coords(igdotssout)
        write(igdotssout,'(1X,A,I3)') 'Experimental condition #', i_tube
        write(igdotssout,fmt_str(2)) 'Vol [cm3] ',knams(1:kgmax),'Temp. [K]'
      end if
      write(igdotssout,fmt_str(1))tt,(gdot(kk),kk=1,kgmax),w(indext)
    end if
    if (header_write) then !Mass balance
      write(imassbalssout,'(1X,2(A,I3))') 'Experimental condition #', i_tube
      write(imassbalssout,fmt_str(2)) 'Vol [cm3]  ', 'Avg MM [g/mol]', 'rho [g/cm3]', 'drho/dt [g/cm3s]','dm/dt [g/s]','mdot_in [g/s]','mdot_out [g/s]','m_gen_surf [g/s]'
    end if
    write(imassbalssout,fmt_str(1))tt,wta,w(indexrho),wtprime(indexrho),wtprime(indexrho)*rlen,fluxmass,fluxmass-wtprime(indexrho)*rlen+abyv*rlen*sum(sdot(1:kgmax)*wt),sum(sdot(1:kgmax)*wt)
  end if

!*********************************************************************************************************
end subroutine write_model_output
!*********************************************************************************************************

!*********************************************************************************************************
subroutine write_tube_restart
!This subroutine writes the compositions and mass fluxes for restarting the
!calculation in a new reactor.
!*********************************************************************************************************

  use global_all; use global_tube1; use update_tube2; use file_handles

  implicit none

  integer :: i, j, phase_no
  character(len=80) :: restart_fmt, i_tube_str
  character(len=15) :: phase
  character(len=16) :: sp_name
  character(len=35) :: ctemp

  write(i_tube_str,'(I10)') nruns  !Internal write converts integer to a string
  restart_fmt='(1X,A35,1X,'//trim(adjustl(i_tube_str))//'(ES15.7E3))'

  write(ituberestartout,*) kgmax+ksmax+1, 'Number of nonzero entries'

  phase_no=1
  do i=1,kgmax+ksmax
    if (i>site_type_lastID(phase_no)) phase_no=phase_no+1
    phase=PhaseNames(phase_no)
    sp_name=knams(i)
    ctemp="'"//trim(adjustl(sp_name))//'/'//trim(adjustl(phase))//'/'//"'"
    write(ituberestartout,restart_fmt) ctemp, restart_save(i,:)
  end do
  ctemp(1:8)='fluxmass'
  ctemp(9:35)=' '
  write(ituberestartout,restart_fmt) ctemp, restart_save(kgmax+ksmax+1,:)
  write(ituberestartout,*) 'EOF'

!*********************************************************************************************************
end subroutine write_tube_restart
!*********************************************************************************************************

!*********************************************************************************************************
subroutine write_rxn_rates_ss
!This subroutine writes net and forward reaction rates.
!*********************************************************************************************************

  use global_all; use global_tube1; use update_tube2; use file_handles

  implicit none

  integer :: i, j
  character(len=80) :: fmt_str

  !Convert to mol/s
  rxn_rates_ss(1:ngrxn,:,:,:)=rxn_rates_ss(1:ngrxn,:,:,:)*rlen
  rxn_rates_ss(ngrxn+1:ngrxn+nsrxn,:,:,:)=rxn_rates_ss(ngrxn+1:ngrxn+nsrxn,:,:,:)*rlen*abyv

  !Header
  write(irxnratess,'(1X,A)') 'Alternating forward and net reaction rates (mol/s) at steady state'

  !Write the body
  fmt_str='(1X,A80,1X,'//trim(int2char(2*nruns))//'(ES18.8E3,1X))'
  do j=1,nBE_coords
    if (iScale>0) then
      scale_targ(:,2)=BE_coords(:,j)
      BE_coord_idx=j
      call write_BE_coords(irxnratess)
    end if
    do i=1,nrxn
      write(irxnratess,fmt_str) rxn_str(i), reshape(rxn_rates_ss(i,:,:,j),(/1,2*nruns/))
    end do
  end do

!*********************************************************************************************************
end subroutine write_rxn_rates_ss
!*********************************************************************************************************

!*********************************************************************************************************
subroutine write_SA_pert(run,T)
!This subroutine writes the perturbation vectors to disk.
!*********************************************************************************************************

  use global_all; use global_tube1; use update_tube2; use file_handles; use globalSA

  implicit none

  integer, intent(in) :: run
  double precision, intent(in) :: T
  integer :: i, j, k, l, field_width
  character(len=80) :: fmt_str(max_SA_param_types), fw, wp
  character, parameter :: dpoint(2)=(/'A','B'/)
  double precision :: pert_val(max(kgmax+ksmax,nrxn,2*nBEP,2*nscale))

  !Create the formatting strings
  field_width=10+write_precision  !This allows for dynamic precision changes
  fw=int2char(field_width)
  wp=int2char(write_precision)

  !Write the output
  do i=1,max_SA_param_types
    if (.not. perturb_param(i)) cycle
    select case (i)
    case (globalSA_H)
      fmt_str(i)='(1X,A10,1X,I8,A1,1X,A8,'//trim(int2char(nspec_params))//&
        '(1X,ES'//trim(fw)//'.'//trim(wp)//'E3))'
      if (run==1) write(iHfpertout,'(1X,A,I3)') 'Experimental condition #', i_tube
    case (globalSA_S)
      fmt_str(i)='(1X,A10,1X,I8,A1,1X,A8,'//trim(int2char(nspec_params))//&
        '(1X,ES'//trim(fw)//'.'//trim(wp)//'E3))'
      if (run==1) write(iSfpertout,'(1X,A,I3)') 'Experimental condition #', i_tube
    case (globalSA_A)
      fmt_str(i)='(1X,A10,1X,I8,A1,1X,A8,'//trim(int2char(nrxn_params))//&
        '(1X,ES'//trim(fw)//'.'//trim(wp)//'E3))'
      if (run==1) write(iApertout,'(1X,A,I3)') 'Experimental condition #', i_tube
    case (globalSA_B)
      fmt_str(i)='(1X,A10,1X,I8,A1,1X,A8,'//trim(int2char(nrxn_params))//&
        '(1X,ES'//trim(fw)//'.'//trim(wp)//'E3))'
      if (run==1) write(iBpertout,'(1X,A,I3)') 'Experimental condition #', i_tube
    case (globalSA_EA)
      fmt_str(i)='(1X,A10,1X,I8,A1,1X,A8,'//trim(int2char(nrxn_params))//&
        '(1X,ES'//trim(fw)//'.'//trim(wp)//'E3))'
      if (run==1) write(iEApertout,'(1X,A,I3)') 'Experimental condition #', i_tube
    case (globalSA_BEP)
      fmt_str(i)='(1X,A10,1X,I8,A1,1X,A10,'//trim(int2char(nBEP))//&
        '(1X,ES'//trim(fw)//'.'//trim(wp)//'E3))'
      if (run==1) write(iBEPpertout,'(1X,A,I3)') 'Experimental condition #', i_tube
    case (globalSA_LSR)
      fmt_str(i)='(1X,A10,1X,I8,A1,1X,A10,'//trim(int2char(nscale*2))//&
        '(1X,ES'//trim(fw)//'.'//trim(wp)//'E3))'
      if (run==1) write(iLSRpertout,'(1X,A,I3)') 'Experimental condition #', i_tube
    end select
  end do

  do k=1,2
    j=1 !This is a counter to keep track of the block of data we are accessing
    do i=1,max_SA_param_types
      if (.not. perturb_param(i)) cycle
      select case (i)
      case (globalSA_H)
        if (lDFTH_SA .or. isenbrute==4) then
          pert_val(1:nspec_params)=SApert_data(j:j+nspec_params,k)%pert_val(H_DFT_loc)
          write(iHfpertout,fmt_str(i)) 'Replicate', run, dpoint(k), 'DFT', pert_val(1:nspec_params)/T
        end if
        if (lGAH_SA) then
          pert_val(1:nspec_params)=SApert_data(j:j+nspec_params,k)%pert_val(H_GA_loc)
          write(iHfpertout,fmt_str(i)) 'Replicate', run, dpoint(k), 'GA', pert_val(1:nspec_params)/T
        end if
        if (lLSR_SA) then
          pert_val(1:nspec_params)=SApert_data(j:j+nspec_params,k)%pert_val(H_LSR_loc)
          write(iHfpertout,fmt_str(i)) 'Replicate', run, dpoint(k), 'LSR', pert_val(1:nspec_params)/T
        end if
        j=j+nspec_params
      case (globalSA_S)
        if (lDFTS_SA .or. isenbrute==4) then
          pert_val(1:nspec_params)=SApert_data(j:j+nspec_params,k)%pert_val(S_DFT_loc)
          write(iSfpertout,fmt_str(i)) 'Replicate', run, dpoint(k), 'DFT', pert_val(1:nspec_params)
        end if
        if (lGAH_SA) then
          pert_val(1:nspec_params)=SApert_data(j:j+nspec_params,k)%pert_val(S_GA_loc)
          write(iSfpertout,fmt_str(i)) 'Replicate', run, dpoint(k), 'GA', pert_val(1:nspec_params)
        end if
        j=j+nspec_params
      case (globalSA_A)
        pert_val(1:nrxn_params)=SApert_data(j:j+nrxn_params,k)%pert_val(Preexp_loc)
        write(iApertout,fmt_str(i)) 'Replicate', run, dpoint(k), pert_val(1:nrxn_params)
        j=j+nrxn_params
      case (globalSA_B)
        pert_val(1:nrxn_params)=SApert_data(j:j+nrxn_params,k)%pert_val(Beta_loc)
        write(iBpertout,fmt_str(i)) 'Replicate', run, dpoint(k), pert_val(1:nrxn_params)
        j=j+nrxn_params
      case (globalSA_EA)
        if (lDFTEA_SA .or. isenbrute==4) then
          pert_val(1:nrxn_params)=SApert_data(j:j+nrxn_params,k)%pert_val(EA_DFT_loc)
          write(iEApertout,fmt_str(i)) 'Replicate', run, dpoint(k), 'DFT', pert_val(1:nrxn_params)/T
        end if
        if (lBEP_SA) then
          pert_val(1:nrxn_params)=SApert_data(j:j+nrxn_params,k)%pert_val(EA_BEP_loc)
          write(iEApertout,fmt_str(i)) 'Replicate', run, dpoint(k), 'BEP', pert_val(1:nrxn_params)/T
        end if
        j=j+nrxn_params
      case (globalSA_BEP)
        do l=j,j+2*nBEP
          if (SApert_data(j,k)%par_type>0) then  !Slope
            pert_val(l)=SApert_data(j,k)%pert_val(BEP_slope_loc)
          else
            pert_val(l)=SApert_data(j,k)%pert_val(BEP_int_loc)
          end if
        end do
        write(iBEPpertout,fmt_str(i)) 'Replicate', run, dpoint(k), 'Slope', pert_val(1:2*nBEP-1:2)
        write(iBEPpertout,fmt_str(i)) 'Replicate', run, dpoint(k), 'Intercept', pert_val(2:2*nBEP:2)
        j=j+2*nBEP
      case (globalSA_LSR)
        do l=j,j+2*nscale
          if (SApert_data(j,k)%par_type>0) then  !Slope
            pert_val(l)=SApert_data(j,k)%pert_val(LSR_slope_loc)
          else
            pert_val(l)=SApert_data(j,k)%pert_val(LSR_int_loc)
          end if
        end do
        write(iLSRpertout,fmt_str(i)) 'Replicate', run, dpoint(k), 'Slope', pert_val(1:2*nscale-1:2)
        !Uncomment the following line if the intercepts of the LSRs are ever perturbed
!        write(iLSRpertout,fmt_str(i)) 'Replicate', run, dpoint(k), 'Intercept', pert_val(2:2*nscale:2)
        j=j+2*nscale
      end select
    end do
  end do

!*********************************************************************************************************
end subroutine write_SA_pert
!*********************************************************************************************************

!*********************************************************************************************************
subroutine write_SA_output(rep)
!This subroutine is responsible for writing the array of solution vectors to
!disk.
!*********************************************************************************************************

  use file_handles; use globalSA; use global_all; use update_tube2; use global_tube1

  implicit none

  integer, intent(in) :: rep  !The number of the replicate
  integer :: i, j, k, field_width
  character(len=80) :: fmt_str(3)
  character, parameter :: crun(2)=(/'A','B'/)
  double precision :: pressure

  field_width=10+write_precision  !This allows for dynamic precision changes
  fmt_str(1)='(2(I15,1X),'//trim(int2char(kmax+2))//'(ES'//trim(int2char(field_width))//'.'&
             //trim(int2char(write_precision))//'E3,1X))'
  fmt_str(2)='(2(A10,6X),'//trim(int2char(kmax+2))//'(A'//trim(int2char(field_width+1))//'))'
  fmt_str(3)='(I15,1X,14X,A1,1X,'//trim(int2char(kmax+2))//'(ES'//trim(int2char(field_width))//'.'&
             //trim(int2char(write_precision))//'E3,1X))'

  !Write a header line at the start of every expt. block
  if (iScale>0) call write_BE_coords(iwSAout)
  write(iwSAout,'(1X,A,I3)') 'Experimental condition #', i_tube
  write(iwSAout,fmt_str(2)) 'Replicate ', 'Perturb   ', knams(1:kgmax+ksmax), 'Temperature [K]', 'Pressure [atm]'

  !Write the initial conditions
  call ckpy(win(indexrho),win(indext),win(1:kgmax),ickwrk,rckwrk,pressure)
  pressure=pressure/patm
  write(iwSAout, fmt_str(1)) rep, 0, win(1:indext), pressure

  !Write the solutions for each perturbation
  do i=1,nparams_total+2
    !Skip the 'B' vector for LSA & DGSA
    if ((isenbrute==2 .or. isenbrute==4) .and. i==2) cycle

    !Skip anything above the 'A' and 'B' vectors for UQ
    if (isenbrute==5 .and. i>2) exit

    !Write the replicate number and perturbation number followed by the solution vector
    call ckpy(w_SA(i,indexrho),w_SA(i,indext),w_SA(i,1:kgmax),ickwrk,rckwrk,pressure)
    pressure=pressure/patm
    if (i<3) then
      write(iwSAout, fmt_str(3)) rep, crun(i), w_SA(i,1:indext), pressure
    else
      write(iwSAout, fmt_str(1)) rep, i-2, w_SA(i,1:indext), pressure
    end if
  end do

  fmt_str(1)='(2(I15,1X),'//trim(int2char(2*nrxn))//'(ES'//trim(int2char(field_width))//'.'&
             //trim(int2char(write_precision))//'E3,1X))'
  fmt_str(2)='(2(A10,6X),A)'
  fmt_str(3)='(I15,1X,14X,A1,1X,'//trim(int2char(2*nrxn))//'(ES'//trim(int2char(field_width))//'.'&
             //trim(int2char(write_precision))//'E3,1X))'

  !Write a header line
  if (iScale>0) call write_BE_coords(irSAout)
  write(irSAout,'(1X,A,I3)') 'Experimental condition #', i_tube
  write(irSAout, fmt_str(2)) 'Replicate ', 'Perturb   ', 'Forward Rates'//repeat(' ',nrxn*(field_width+1)-13)//'Net Rates'

  !Write the solutions for each perturbation
  do i=1,nparams_total+2
    !Skip the 'B' vector for LSA & DGSA
    if ((isenbrute==2 .or. isenbrute==4) .and. i==2) cycle

    !Skip anything above the 'A' and 'B' vectors for UQ
    if (isenbrute==5 .and. i>2) exit

    !Write the replicate number and perturbation number followed by the solution vector
    if (i<3) then
      write(irSAout, fmt_str(3)) rep, crun(i), rop_SA(i,:)
    else
      write(irSAout, fmt_str(1)) rep, i-2, rop_SA(i,:)
    end if
  end do

!*********************************************************************************************************
end subroutine write_SA_output
!*********************************************************************************************************

!*********************************************************************************************************
character(len=8) function int2char(val)
!This function is a wrapper for an internal write to convert an integer into
!a string.
!*********************************************************************************************************

  implicit none

  integer, intent(in) :: val
  character(len=8) :: work

  write(work,'(I8)') val
  int2char=adjustl(work)

!*********************************************************************************************************
end function int2char
!*********************************************************************************************************

!*********************************************************************************************************
subroutine write_BE_coords(file_LU)
!This subroutine is used for writing the binding energy coordinates to the
!specified logical unit number.
!*********************************************************************************************************

  use global_all

  implicit none

  integer, intent(in) :: file_LU
  integer :: kk, name_len, specID
  character(len=80) :: fmt_str
  character(len=17) :: atom_names(natoms_scale)

  fmt_str='(1X,A,I5,A1,1X'//trim(int2char(natoms_scale*2))//'(A17,1X,F9.3,1X))'
  do kk=1,natoms_scale
    specID=nint(scale_targ(kk,1))
    name_len=len(trim(adjustl(knams(specID))))
    atom_names(kk)=trim(adjustl(knams(specID)))//':'//repeat(' ',17-name_len-1)
  end do
  write(file_LU,fmt_str) 'Atomic binding energies (kcal/mol), Set # ', BE_coord_idx, ':', &
    (atom_names(kk), scale_targ(kk,2)*Rgas_kcal, kk=1,natoms_scale)

!*********************************************************************************************************
end subroutine write_BE_coords
!*********************************************************************************************************

!*********************************************************************************************************
subroutine rxnenth(irxn,hospec,dhr)
!This subroutine calculates heats of reactions.
!*********************************************************************************************************

  use global_all

  implicit none

  integer, intent(in) :: irxn(2) !Vector of start and stop reaction numbers
  double precision, intent(in) :: hospec(kmax)
  double precision, intent(out) :: dhr(irxn(2)-irxn(1)+1)
  integer :: i, j, spec

  !Initialize dhr
  dhr=0.0

  !Loop over reactions calculating reaction energy
  do i=irxn(1),irxn(2)
    do j=2,stoich_matrix_sparse(i,1)+1
      spec=stoich_matrix_sparse(i,j)
      dhr(i)=dhr(i)+stoich_matrix(spec,i)*hospec(spec)
    end do
  end do

!*********************************************************************************************************
end subroutine rxnenth
!*********************************************************************************************************

!*********************************************************************************************************
end module output
!*********************************************************************************************************
