!This is the main program source file with subroutines which are related to
!memory allocation, initialization, etc. Make sure to update the revision
!number in the output header before you commit changes to the repository! The
!value (revno) is a parameter defined in modules.f90.
!*********************************************************************************************************
      program main
!*********************************************************************************************************
        use global_all
        implicit none
!-------------------------------------------------------------------------------------------------------
        double precision  :: glob_time_var_start, glob_time_var_stop
        character(LEN=10) :: glob_time_str(2)
        character(len=10) :: time
        character(len=13) :: pretty_time  !The time in hh:mm:ss.ms format
!-------------------------------------------------------------------------------------------------------
        CALL DATE_AND_TIME(TIME=glob_time_str(1))
        CALL CPU_TIME(glob_time_var_start)
        write(*,*)'**************************************************'
        write(*,*)'*******    Welcome to the Vlachos Group    *******'
        write(*,*)'*******       SURFACE CHEMKIN-driven       *******'
        write(*,*)'******* Microkinetic Modeling Reactor Code *******'
        write(*,'(1X,A,I0,A)')&
          '*******            Revision #',revno,'           *******'
        write(*,*)'**************************************************'
!-------------------------------------------------------------------------------------------------------
        call mech_init      !Initialize the work arrays and other mechanistic information
        call model_allocate !Allocate memory for the model-specific arrays
        call tube_driver    !Call the driver routine for solving the model
        CALL DATE_AND_TIME(TIME=glob_time_str(2))
        CALL CPU_TIME(glob_time_var_stop)
        time=glob_time_str(1)
        pretty_time=time(1:2)//':'//time(3:4)//':'//time(5:10)
        WRITE(*,*) 'Global Start time =    ', pretty_time
        time=glob_time_str(2)
        pretty_time=time(1:2)//':'//time(3:4)//':'//time(5:10)
        WRITE(*,*) 'Global Stop time  =    ', pretty_time
        WRITE(*,'(1X,A,F10.2,A)') 'Global Elapsed Time = ',&
          glob_time_var_stop-glob_time_var_start,'  s'
        stop
!*********************************************************************************************************
      end program main
!*********************************************************************************************************

!*********************************************************************************************************
      subroutine tube_driver
!This subroutine is the main driver subroutine for the code. It calls the input routines, the model
!simulation routine, and the output routines.
!*********************************************************************************************************
        use global_all; use global_tube1; use update_tube2; use input; use model
        use output; use sensitivity

        implicit none

        double precision  ::  t, z, tt, w(neqns)
        integer     ::  kk, i, scale_iter
        integer     ::  ierr
        double precision  ::  veloin, conv, HORT(ksmax), rate
        type(runs),allocatable :: DOE_data(:)     !Data read in from DOE.inp
        character(len=80) :: asterisks

!-------------------------------------------------------------------------------------------------------

        asterisks=repeat('*',80)
        call inco_tube     ! Get input parameters
        call get_tube_conditions  !Get all run-specific parameters and store them
        call open_output !Open all the output files for writing on the first iterations

        !Initialize DOE variables
        if (lDOE) then
          if (.not. allocated(DOE_data)) allocate(DOE_data(nruns-1),STAT=ierr)
          if (.not. allocated(DOE_data)) then
            write(*,*) 'Error allocating DOE memory. Stopping.'
            STOP
          end if
          call inco_DOE(DOE_data)
        endif

!-------------------------------------------------------------------------------------------------------
!.........Loop over number of model responses:
        do_01: do i_tube=1,nruns

          write(*,*) asterisks(1:50)
          write(*,008) ' ******* BEGIN Experimental Condition # ',&
            i_tube,' ********'
          write(*,*) asterisks(1:50)

          do_scale: do scale_iter=1,nBE_coords !Number of scaling relations points

            !Check for global SA or UQ and skip the nominal solution.
            if (abs(isenbrute)>1) then
              write(*,*) ' Skipping nominal solution point &
                &for SA or UQ calculation'
              exit
            end if

            !This is a run requiring the nominal point, keep going
            BE_coord_idx=scale_iter
            rckwrk=rckwrk_save  !Reset the work arrays
            rskwrk=rskwrk_save
            call update_tube_conditions

            if (iScale>0) scale_targ(:,2)=BE_coords(:,scale_iter)
            if (lBEP) call zero_BEP_EA   ! Reset activation energies for reactions using BEPs to 0
            if (lDOE) then     ! Reset the BEP and LSR arrays which may have been altered
              if (lBEP) BEP_def=BEP_def_save
              if (iScale>0) scale_slope=scale_slope_save
            end if

            T=win(indext)
            if (ltra .or. i_tube==1) then !For purely transient problem or first run, use input conditions
              w=win
            else  !Using continuation for the surface
              w(1:kgmax)=win(1:kgmax)  !keep coverage for next run, reset gas
              w(indext)=T    ! Temperature
              w(indexrho)=rhoin
            end if

            call convbkd(w)
            if (lDOE .and. i_tube > 1) call apply_DOE(DOE_data(i_tube-1))

            !Calculate the change in heat of formation using the scaling relations. This subroutine
            !is in sklib.f
            if (iScale==3) then
              write(*,*) asterisks(1:50)
              write(*,'(1X,A,I0,A,I0)') 'Binding energy coordinate #', &
                scale_iter,' of ', nBE_coords
              write(*,*) 'Atomic binding energy (kcal/mol) values are:'
              do i=1,natoms_scale
                write(*,'(2X,A,F6.2)') knams(nint(scale_targ(i,1))), &
                  scale_targ(i,2)*Rgas_kcal
              end do
              write(*,*) asterisks(1:50)
            end if

            !Apply the LSR corrections
            CALL SKSCALE(ISKWRK,RSKWRK,nscale,natoms_scale,iScale,&
              scale_bind_mode,scale_slope,scale_ref,scale_targ,&
              cov_factor,cov_factor_adjust)

            if (lStatpQ) call SKSTATPQ(ISKWRK,RSKWRK,T,StatpQ) !Apply the StatpQ correction

            call save_params(T) !Store the original parameters
            parasurf=rka; paragas=rkag                    ! Pre-exponentials

            !Calculate coverage effects and BEP estimates
            if (lcov) call SKCOVERAGE(ISKWRK,RSKWRK,cov_matrix,&
              cov_matrix2,cov_matrix3,thresh_cov,thresh_cov2,A6,DH_save,EA_save,w,omega,cov_factor)
            if (lBEP) then
              call ckBEP(ICKWRK,RCKWRK,T,nBEP,BEP_def,BEP_rxn(1:ngrxn,:),reg)
              call skBEP(ISKWRK,RSKWRK,T,nBEP,BEP_def,BEP_rxn(ngrxn+1:ngrxn+nsrxn,:),re)
            end if

!-----------------------------------------------------------------------------------------------------
            if (scale_iter==1) then
              write(*,*) asterisks(1:50)
              write(*,*)'Inlet MASS fractions:'
              do i=1,kgmax
                if (yin(i)>1.e-8) write(*,457)knams(i),yin(i)
              end do
              write(*,*)'Inlet MOLE fractions:'
              do i=1,kgmax
                if (xin(i)>1.e-8) write(*,457)knams(i),xin(i)
              end do
              write(*,*)'Operating conditions:'
              write(*,164)'pressure [atm]', p/patm
              write(*,164)'temperature [K]', t
               write(*,164)'density [gm/cm3]', rhoin
              write(*,164)'reactor volume [cm3]', rlen
              if (irxtr>1) then
                write(*,164)'velocity [cm3/s]', velo
                write(*,164)'residence time [s]', rlen/velo
                write(*,164)'mass flow [gm/s]', fluxmass
              end if
              write(*,164)'Ac/Vr [1/cm]', abyv
              if (itube_restart/=0 .and. irxtr>1) then
                fluxmass=fluxmass_inlet
                write(*,164)'mass flux restored', fluxmass
              end if
            end if
!-----------------------------------------------------------------------------------------------------

!This is the main routine where the model is actually solved.
!===================================================================================
            call tube_sub(w,conv,rate)  !=== Solve reactor model for given T and Y =
!===================================================================================

            !Write basic output common to all reactor models
            t=w(indext)
            call write_conversion(i_tube,t,conv,rate)
            if (.not. verbose_rpa .and. mrpa/=0) call write_rpa(rlen,T)

            !Save reaction rates
            rxn_rates_ss(:,1,i_tube,scale_iter)=ropfwd
            rxn_rates_ss(:,2,i_tube,scale_iter)=rop

            if (lsenbrute .and. isenbrute==1) then
              if (irxtr==3) then
                FIM_SA_matrix=FIM_SA_matrix/rlen
              else
                FIM_SA_matrix=FIM_SA_matrix/rtime
              end if
              !convert to mol/s
              FIM_SA_matrix(1:ngrxn,:)=FIM_SA_matrix(1:ngrxn,:)*rlen
              FIM_SA_matrix(ngrxn+1:ngrxn+nsrxn,2)=FIM_SA_matrix(ngrxn+1:ngrxn+nsrxn,2)*rlen*abyv
              call NSC_matrix_write(1)
            end if

!---------------------------------------------------------------------------------------------------

          end do do_scale

          !Sensitivity Analysis
          if (lsenbrute .and. isenbrute>1) then
            SA_called=.true.
            write(*,*) asterisks(1:50)
            write(*,*) 'Beginning sensitivity analysis'
            write(*,*) asterisks(1:50)
            call tube_SA
            write(*,*) asterisks(1:50)
            write(*,*) 'Sensitivity analysis complete'
            write(*,*) asterisks(1:50)
            SA_called=.false. !Reset the flag so that subsequent runs print correctly
          end if

          write(*,*) asterisks(1:50)
          write(*,008)' ********* END Experimental Condition # ',&
            i_tube,' ********'

        end do do_01
!-------------------------------------------------------------------------------------------------------

        if (irxtr>1) call write_tube_restart  !Only write it for flow reactors
        call write_rxn_rates_ss !Write steady state reaction rates
        call close_output

008     format(A,I2,A)
164     format(A25,es15.8)
457     format(A25,e15.8)
!*********************************************************************************************************
      end subroutine tube_driver
!*********************************************************************************************************

!*********************************************************************************************************
      subroutine mech_init
!       This subroutine is responsible for getting the dimensions of the work arrays
!       as well as various other useful indices.
!*********************************************************************************************************
        use global_all; use file_handles

        implicit none

        integer :: nfit, kk  !CKINDX dummy values
        integer :: nfsurf, nlsurf, nfbulk, nlbulk !SKINDX dummy values
        integer, allocatable :: kfirst(:), klast(:)  !SKPKK values
        integer :: mem_stat, mem_stat_total=0
        logical :: kerr

        !Get the sizes of the work arrays from the binary linking files
        open(unit=linc,form='unformatted',file=cklink)
        open(unit=linksk,form='unformatted',file=sklink)
        call cklen(linc,lout,leniwk,lenrwk,lencwk)    !Gas-phase dimensions
        call sklen(linksk,lout,lenisk,lenrsk,lencsk)  !Surface dimensions
        close(linc)
        close(linksk)

        call work_allocate  !Allocates the memory for the work arrays

        !Initialize the work arrays so they can be used to extract useful
        !information about the problem dimensions.
        call ckinit(leniwk,lenrwk,lencwk,linc,lout,ickwrk,rckwrk,cckwrk)
        call skinit(lenisk,lenrsk,lencsk,linksk,lout,iskwrk,rskwrk,cskwrk)
        rckwrk_save=rckwrk  !Make backup copies of the as-initialized work arrays
        rskwrk_save=rskwrk

        !Get the number of species, phases, reactions, etc.
        call ckindx(ickwrk,rckwrk,nelm,kk,ngrxn,nfit)  !Gas-phase
        call skindx(iskwrk,nelm,kgmax,ksmax,kbmax,kmax,&
            nphases,nsphases,nfsurf,nlsurf,nbphases,nfbulk,nlbulk,nsrxn)  !Surface
        ngphases=nphases-nsphases-nbphases  !Should be 1
        nrxn=ngrxn+nsrxn

        !Extract the # of temperatures used for the thermo of each species (currently 3)
        call skmxtp(iskwrk,maxtp)

        !Do some basic checking of the problem formulation and stop if errors
        if (kk/=kgmax) then
          write(*,*) 'kk = ', kk, ' kgmax = ', kgmax
          write(*,*) &
            'Error in the number of gas-phase species in the surface &
            mechanism. Stopping.'
          STOP
        end if
        if (ngphases/=1) then
          write(*,*) 'Wrong number of gaseous phases. Stopping.'
          STOP
        end if
        if (kbmax/=nbphases) then
          write(*,*) &
            'Number of bulk phases and bulk species must be equal. &
            Stopping.'
          STOP
        end if

        !Allocate memory for the SKPKK and SKSYMP return values
        if (.not. allocated(SpecTot)) &
          allocate(SpecTot(nphases),stat=mem_stat)
        if (.not. allocated(SpecTot)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(kfirst)) &
          allocate(kfirst(nphases),stat=mem_stat)
        if (.not. allocated(kfirst)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(klast)) &
          allocate(klast(nphases),stat=mem_stat)
        if (.not. allocated(klast)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(PhaseNames)) &
          allocate(PhaseNames(nphases),stat=mem_stat)
        if (.not. allocated(PhaseNames)) mem_stat_total=mem_stat_total+1

        if (mem_stat_total>0) then
          write(*,*) 'Error allocating memory for SKPKK arrays. &
            Stopping.'
          STOP
        end if

        !Get the species indices for each phase and the phase names
        call skpkk(iskwrk,SpecTot,kfirst,klast)
        call sksymp(iskwrk,cskwrk,lout,PhaseNames,kerr)

        !Allocate memory for the MultiSite arrays and initialize them
        if (.not. allocated(site_type_lastID)) &
          allocate(site_type_lastID(nsphases+1), STAT=mem_stat)
        if (.not. allocated(site_type_lastID)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(site_type_firstID)) &
          allocate(site_type_firstID(nsphases+1), STAT=mem_stat)
        if (.not. allocated(site_type_firstID)) mem_stat_total=mem_stat_total+1
        if (mem_stat_total>0) then
          write(*,*) 'Error allocating memory for &
            site type array. Stopping.'
          STOP
        end if
        do kk=1,nsphases+1
          site_type_lastID(kk)=klast(kk) !Assumes that the free site is the last species
          site_type_firstID(kk)=kfirst(kk)
        end do

        !Can now deallocate the SKPKK arrays
        deallocate(kfirst)
        deallocate(klast)

        !Now allocate and initialize various mechanism related arrays
        call mech_allocate  !Allocate space for mechanism-dependent arrays
        call sksyms (iskwrk,cskwrk,lout,knams,kerr) ! Species names
        call sksyme (iskwrk,cskwrk,lout,enams,kerr) ! Element names
        call ckwt (ickwrk,rckwrk,wt)            ! Molecular weights
        call sksden (iskwrk,rskwrk,sden)          ! Site density
        call stoich                ! Stoichiometry
        call skcov (iskwrk,kocc)                             ! Site occupancy
        call get_rxn_str !Reaction strings
        call skncf(nelm,iskwrk,elem_comp) !Elemental compositions

!*********************************************************************************************************
      end subroutine mech_init
!*********************************************************************************************************

!*********************************************************************************************************
      subroutine work_allocate
!       This subroutine handles the allocation of memory for the work arrays.
!*********************************************************************************************************
        use global_all
        implicit none

        integer               ::      mem_stat, mem_stat_total=0

        !1. Integer work arrays
        if (.not. allocated(ickwrk)) allocate(ickwrk(leniwk), STAT=mem_stat)
        if (.not. allocated(ickwrk)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(iskwrk)) allocate(iskwrk(lenisk), STAT=mem_stat)
        if (.not. allocated(iskwrk)) mem_stat_total=mem_stat_total+1

        !2. Real work arrays (both working and saved)
        if (.not. allocated(rckwrk)) allocate(rckwrk(lenrwk), STAT=mem_stat)
        if (.not. allocated(rckwrk)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(rskwrk)) allocate(rskwrk(lenrsk), STAT=mem_stat)
        if (.not. allocated(rskwrk)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(rckwrk_save)) allocate(rckwrk_save(lenrwk), STAT=mem_stat)
        if (.not. allocated(rckwrk_save)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(rskwrk_save)) allocate(rskwrk_save(lenrsk), STAT=mem_stat)
        if (.not. allocated(rskwrk_save)) mem_stat_total=mem_stat_total+1

        !3. Character work arrays
        if (.not. allocated(cckwrk)) allocate(cckwrk(lencwk), STAT=mem_stat)
        if (.not. allocated(cckwrk)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(cskwrk)) allocate(cskwrk(lencsk), STAT=mem_stat)
        if (.not. allocated(cskwrk)) mem_stat_total=mem_stat_total+1

        if (mem_stat_total>0) then
          write(*,*) 'Total number of memory allocation errors is', mem_stat_total
          write(*,*) 'Stopping due to memory allocation failure'
          STOP
        end if

!*********************************************************************************************************
      end subroutine work_allocate
!*********************************************************************************************************

!*********************************************************************************************************
      subroutine mech_allocate
!       This subroutine handles the allocation of memory for chemistry
!       related arrays.
!*********************************************************************************************************
        use global_all
        implicit none
        integer :: mem_stat, mem_stat_total=0, i

        !1. Chemistry related (stoichiometry, etc.)
        if (.not. allocated(sden)) allocate(sden(nphases), STAT=mem_stat)
        if (.not. allocated(sden)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(wt)) allocate(wt(kgmax), STAT=mem_stat)
        if (.not. allocated(wt)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(stoich_matrix)) allocate(stoich_matrix(kmax,nrxn), STAT=mem_stat)
        if (.not. allocated(stoich_matrix)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(rxn_str)) allocate(rxn_str(nrxn), STAT=mem_stat)
        if (.not. allocated(rxn_str)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(kocc)) allocate(kocc(kmax), STAT=mem_stat)
        if (.not. allocated(kocc)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(knams)) allocate(knams(kmax), STAT=mem_stat)
        if (.not. allocated(knams)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(enams)) allocate(enams(nelm), STAT=mem_stat)
        if (.not. allocated(enams)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(elem_comp)) allocate(elem_comp(nelm,kmax), STAT=mem_stat)
        if (.not. allocated(elem_comp)) mem_stat_total=mem_stat_total+1

        !2. Arrhenius & coverage parameters
        if (.not. allocated(rka)) allocate(rka(nsrxn), STAT=mem_stat)
        if (.not. allocated(rka)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(rb)) allocate(rb(nsrxn), STAT=mem_stat)
        if (.not. allocated(rb)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(re)) allocate(re(nsrxn), STAT=mem_stat)
        if (.not. allocated(re)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(rkag)) allocate(rkag(ngrxn), STAT=mem_stat)
        if (.not. allocated(rkag)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(rbg)) allocate(rbg(ngrxn), STAT=mem_stat)
        if (.not. allocated(rbg)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(reg)) allocate(reg(ngrxn), STAT=mem_stat)
        if (.not. allocated(reg)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(A6)) allocate(A6(kmax,2), STAT=mem_stat)
        if (.not. allocated(A6)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(A7)) allocate(A7(kmax,2), STAT=mem_stat)
        if (.not. allocated(A7)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(cov_matrix)) &
           allocate(cov_matrix(ksmax,ksmax), STAT=mem_stat)
        if (.not. allocated(cov_matrix)) mem_stat_total=mem_stat_total+1
		if (.not. allocated(cov_matrix2)) &
           allocate(cov_matrix2(ksmax,ksmax), STAT=mem_stat)
        if (.not. allocated(cov_matrix2)) mem_stat_total=mem_stat_total+1
		if (.not. allocated(cov_matrix3)) &
           allocate(cov_matrix3(ksmax,ksmax), STAT=mem_stat)
        if (.not. allocated(cov_matrix3)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(thresh_cov)) &
           allocate(thresh_cov(ksmax,ksmax), STAT=mem_stat)
        if (.not. allocated(thresh_cov)) mem_stat_total=mem_stat_total+1
		if (.not. allocated(thresh_cov2)) &
           allocate(thresh_cov2(ksmax,ksmax), STAT=mem_stat)
        if (.not. allocated(thresh_cov2)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(cov_factor)) allocate(cov_factor(ksmax))
        if (.not. allocated(cov_factor)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(DH_save)) &
           allocate(DH_save(nsrxn), STAT=mem_stat)
        if (.not. allocated(DH_save)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(EA_save)) &
           allocate(EA_save(nsrxn), STAT=mem_stat)
        if (.not. allocated(EA_save)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(StatpQ)) allocate(StatpQ(ksmax), STAT=mem_stat)
        if (.not. allocated(StatpQ)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(omega)) allocate(omega(nsrxn), STAT=mem_stat)
        if (.not. allocated(omega)) mem_stat_total=mem_stat_total+1

        if (mem_stat_total>0) then
          write(*,*) 'Total number of memory allocation errors is', mem_stat_total
          write(*,*) 'Stopping due to memory allocation failure'
          STOP
        end if

!*********************************************************************************************************
      end subroutine mech_allocate
!*********************************************************************************************************

!*********************************************************************************************************
      subroutine model_allocate
!This subroutine allocates memory for the arrays in the modules global_tube1 and update_tube2
!*********************************************************************************************************
        use global_all; use global_tube1; use update_tube2

        implicit none

        integer :: mem_stat, mem_stat_total=0

        neqns=kgmax+ksmax+ktt+1

        !Part 1 -- global_tube1 variables
        indext=neqns-1  !Temperature
        indexrho=neqns    !Density
        if (ngrxn==0) then
          gas_chem=.false.
        else
          gas_chem=.true.
        end if
        if (nsrxn==0) then
          surf_chem=.false.
        else
          surf_chem=.true.
        end if

        !Part 2 -- update_tube2 variables and arrays
        if (.not. allocated(actin)) allocate(actin(kmax), STAT=mem_stat)
        if (.not. allocated(actin)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(act)) allocate(act(kmax), STAT=mem_stat)
        if (.not. allocated(act)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(xin)) allocate(xin(kgmax), STAT=mem_stat)
        if (.not. allocated(xin)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(x)) allocate(x(kgmax), STAT=mem_stat)
        if (.not. allocated(x)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(yin)) allocate(yin(kgmax), STAT=mem_stat)
        if (.not. allocated(yin)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(y)) allocate(y(kgmax), STAT=mem_stat)
        if (.not. allocated(y)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(win)) allocate(win(neqns), STAT=mem_stat)
        if (.not. allocated(win)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(gdot)) allocate(gdot(kgmax), STAT=mem_stat)
        if (.not. allocated(gdot)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(sdot)) allocate(sdot(kmax), STAT=mem_stat)
        if (.not. allocated(sdot)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(sitdot)) allocate(sitdot(nphases), STAT=mem_stat)
        if (.not. allocated(sitdot)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(rop)) allocate(rop(nrxn), STAT=mem_stat)
        if (.not. allocated(rop)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(ropfwd)) allocate(ropfwd(nrxn), STAT=mem_stat)
        if (.not. allocated(ropfwd)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(rpa)) allocate(rpa(kmax-kbmax,nsrxn), STAT=mem_stat)
        if (.not. allocated(rpa)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(rpag)) allocate(rpag(kgmax,ngrxn), STAT=mem_stat)
        if (.not. allocated(rpag)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(parasurf)) allocate(parasurf(nsrxn), STAT=mem_stat)
        if (.not. allocated(parasurf)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(paragas)) allocate(paragas(ngrxn), STAT=mem_stat)
        if (.not. allocated(paragas)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(wentr)) allocate(wentr(neqns), STAT=mem_stat)
        if (.not. allocated(wentr)) mem_stat_total=mem_stat_total+1
        if (.not. allocated(wexit)) allocate(wexit(neqns), STAT=mem_stat)
        if (.not. allocated(wexit)) mem_stat_total=mem_stat_total+1

        if (mem_stat_total>0) then
          write(*,*) 'Total number of memory allocation errors is', mem_stat_total
          write(*,*) 'Stopping due to memory allocation failure'
          STOP
        end if

!*********************************************************************************************************
      end subroutine model_allocate
!*********************************************************************************************************

!*********************************************************************************************************
      subroutine stoich
!       This subroutine gets the stoichiometric coefficients.
!*********************************************************************************************************
        use global_all

        implicit none

        integer :: nuki(kgmax,ngrxn)  !Gas phase stoichiometric coefficients
        integer :: kstoic(nsrxn,kmax), nstoic(nsrxn,nphases)  !Surface stoichiometric coefficients
        integer :: nspec_list(nrxn)  !Number of species per reaction
        integer :: i, j, k

        !Get gas phase stoichiometric coefficients
        call cknu(kgmax,ickwrk,rckwrk,nuki)

        !Get surface phase stoichiometric coefficients
        call sknu(nsrxn,iskwrk,rskwrk,kstoic,nstoic)

        !Initialize the total stoichiometric matrix
        stoich_matrix=0
        stoich_matrix(1:kgmax,1:ngrxn)=nuki  !nuki is kgmax * ngrxn dimensions
        stoich_matrix(1:kmax,ngrxn+1:ngrxn+nsrxn)=transpose(kstoic) !kstoic is nsrxn * kmax dimensions

        !Count the total number of species needed
        nspec_list=count(stoich_matrix/=0,1)

        !Allocate the sparse stoichiometric matrix. This is simply a list
        !of all the species involved in the given reaction. The stoichiometric
        !coefficients can then be extracted from the full matrix. The first
        !column contains the total number of non-zero entries.
        if (.not. allocated(stoich_matrix_sparse)) &
     &    allocate(stoich_matrix_sparse(nrxn,maxval(nspec_list)+1))

        !Initialize the sparse stoichiometric matrix
        stoich_matrix_sparse=0
        stoich_matrix_sparse(:,1)=nspec_list

        !Assign the values
        do j=1,nrxn
          k=1 !This keeps track of where the current non-zero entry is
          do i=1,kmax
            if (stoich_matrix(i,j)/=0) then
              k=k+1
              stoich_matrix_sparse(j,k)=i
            end if
          end do
        end do

!*********************************************************************************************************
      end subroutine stoich
!*********************************************************************************************************

!*********************************************************************************************************
      subroutine get_rxn_str
!Assembles an array with the reaction equations in string form using CHEMKIN-inherent procedures.
!*********************************************************************************************************

        use global_all; use global_tube1; use update_tube2; use file_handles

        implicit none

        integer :: irxn
        integer :: lt
        logical :: kerr

!-------------------------------------------------------------------------------------------------------

        do irxn=1,ngrxn
          call cksymr(irxn,lout,ickwrk,rckwrk,cckwrk,lt,rxn_str(irxn),kerr)
        end do
        do irxn=1,nsrxn
          call sksymr(irxn,lout,iskwrk,rskwrk,cskwrk,lt,rxn_str(ngrxn+irxn),kerr)
        end do

!*********************************************************************************************************
      end subroutine get_rxn_str
!*********************************************************************************************************
