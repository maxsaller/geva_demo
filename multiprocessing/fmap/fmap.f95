! Simple correlation functions via mapping approach methods !
! including LSC-IVR, PBME and a number of mLSC approaches.  !
! Integration via diagonalisation or velocity verlet.       !
! Copyright: Max Saller 2020                                !
program fmap

    use variables
    implicit none
    integer :: t,ts,pctg, count
    character(len=1) :: cr = char(13)

    ! START TIMING
    t_eom = 0.d0
    t_ops = 0.d0
    t_traj = 0.d0
    call cpu_time( t_total )

    ! HEADER
    write(6,"(a50,//,23x,a4,23x,//,2x,a46,2x,//,a50,/)") repeat("#",50),&
    "FMAP","Correlation functions via the mapping approach",repeat("#",50)

		! INITIALIZE RANDOM NUMBER GENERATOR FROM SYSTEM CLOCK
    call system_clock(count)
    call RLUXGO(4,count,0,0)

    ! READ INPUT FILE
    call read_input()

    ! ALLOCATE ARRAYS
    call allocate_arrays()

    ! SET UP BATH
    call spectral_density()

    ! MONTE CARLO LOOP
    do t = 1,ntraj

        ! TIMING
        call cpu_time( start )

        ! REPORT TRAJECTORY PROGRESS
        if ( t > 1 .and. mod(t,(ntraj/100)) == 0 ) then
            pctg = floor(1.d2*t/dble(ntraj))
            write(6, "(a22,i3,a1,1x,52a)", advance="no") &
            "RUNNING TRAJECTORIES: ", pctg, "%", &
            "[", repeat("#", pctg/2), repeat(" ", 50-pctg/2), "]", cr
            flush(6)
        end if

        ! SAMPLE INITIAL CONDITIONS
        call sample_nuclear()
        call sample_electronic()

        ! CALCULATE TIME ZERO OPERATORS AND OBSERVABLE VALUES
        call cpu_time( start1 )
        call time_zero_ops()
        call accumulate_obs(1)
        !call proj_free_inputs(1)
        call cpu_time( end1 )
        t_ops = t_ops + end1 - start1

        ! TRAJECTORY LOOP
        do ts = 1,tsteps

            ! MAKE DYNAMICS STEP
            call cpu_time( start1 )
            if ( intgt == "diagonalise" .or. intgt == "diag" ) then
                call step_diag()
            else if ( intgt == "velocityverlet" .or. intgt == "vv") then
                call step_vverlet()
            else
                write(6,*) "ERROR: integrator type must be either ",&
                           "'diagonalise'/'diag' or 'velocityverlet'/'vv'!"
                stop
            end if
            call cpu_time( end1 )
            t_eom = t_eom + end1 - start1

            ! CALCULATE TIME t OPERATORS AND ACCUMULATE OBSERVABLES
            call cpu_time( start1 )
            call time_t_ops()
            call accumulate_obs(ts+1)
            !call proj_free_inputs(ts+1)
            call cpu_time( end1 )
            t_ops = t_ops + end1 - start1

        end do

        ! TIMING
        call cpu_time( end )
        t_traj = t_traj + end - start

    end do

    ! AVERAGE AND OUTPUT OBSERVABLES
    call average_obs() 

    ! DEALLOCATE ARRAYS
    call deallocate_arrays()

    ! FINISH TIMING
    call cpu_time( end )
    t_total = end - t_total
    t_traj = t_traj / dble(ntraj)
    t_ops = t_ops / dble(ntraj)
    t_eom = t_eom / dble(ntraj)
    write(6, "(/'TIMING:')")
    write(6, "('- Total time:            ',1x,f12.3,1x,'s')") t_total
    write(6, "('- Avg. trajectory time:  ',1x,f12.3,1x,'s')") t_traj
    write(6, "('- Avg. integration time: ',1x,f12.3,1x,'s')") t_eom
    write(6, "('- Avg. op/obs calc. time:',1x,f12.3,1x,'s')") t_ops

    write(6,"(/'SIMULATION COMPLETE!')")

end program fmap


! Reads the "input" file and prints input arguments to log
subroutine read_input()

    use variables
    implicit none
    integer :: iost=1
    character(len=80) :: dum

    open(11, file="input", action="read")

    read(11, '(A)') dum
    read(11, *) dum, S
    read(11, *) dum, F
    read(11, *) dum, epsilon
    read(11, *) dum, delta
    read(11, *) dum, kondo
    read(11, *) dum, omegac
    read(11, *) dum, omegamax
    read(11, *) dum, beta
    read(11, *) dum, discretize
    read(11, '(A)') dum
    read(11, *)

    read(11, '(A)') dum
    read(11, *) dum, ntraj
    read(11, *) dum, tsteps
    read(11, *) dum, dt
    read(11, *) dum, intgt
    read(11, '(A)') dum
    read(11, *)

    read(11, '(A)') dum
    read(11, *) dum, Aop
    read(11, *) dum, Bop
    read(11, *) dum, electronic
    read(11, '(A)') dum

    rewind(11)
    write(6,"('INPUT ARGUMENTS PARSED:'/)")
    do while ( .true. )
        read(11,'(A)',iostat=iost) dum
        if ( iost < 0 ) exit
        write(6,*) "   ", dum
    end do
    write(6,*)

    close(11)

end subroutine read_input


! Allocates all relevant arrays (and zeroes where required)
subroutine allocate_arrays()

    use variables
    implicit none
    double precision :: zero(1) = 0.d0

    ! LAPACK WORK ARRAY
    allocate( work(1) )
    call dsyev("V", "U", S, zero, S, zero, work, -1, info)
    lenwork = int(work(1))
    deallocate(work)
    allocate( work(lenwork) )

    ! TRAJECTORY ARRAYS
    allocate( xn(F) )
    allocate( pn(F) )
    allocate( XE(S) )
    allocate( PE(S) )

    ! BATH ARRAYS
    allocate( coeff(F) )
    allocate( omega(F) )

    ! POTENTIAL AND FORCE MATRICES
    allocate( G0(F) )
    allocate( V(S,S) )
    allocate( G(F,S,S) )

    ! POPULATION ARRAYS
    allocate( pop_0(S) )
    allocate( pop_t(S) )
    allocate( coh_0(S) )
    allocate( coh_t(S) )
    allocate( Qop_0(S) )
    allocate( Qop_t(S) )
    allocate( Cpop(tsteps+1,S,S) )
    allocate( Cimp(tsteps+1,S,S) )
    allocate( Cdot(tsteps+1,S,S) )
    allocate( Cddot(tsteps+1,S,S) )
    Cpop(:,:,:) = 0.d0
    Cimp(:,:,:) = 0.d0
    Cdot(:,:,:) = 0.d0
    Cddot(:,:,:) = 0.d0

    ! allocate( Epop(tsteps+1,S) )
    ! allocate( Eimp(tsteps+1,S) )
    ! Epop(:,:) = 0.d0
    ! Eimp(:,:) = 0.d0

    ! PROJECTION FREE INPUT ARRAYS
    ! allocate( F1(tsteps+1,S,S,S,S) )
    ! allocate( F2(tsteps+1,S,S,S,S) )
    ! allocate( K1(tsteps+1,S,S,S,S) )
    ! allocate( K3(tsteps+1,S,S,S,S) 
    ! F1(:,:,:,:,:) = 0.d0
    ! F2(:,:,:,:,:) = 0.d0
    ! K1(:,:,:,:,:) = 0.d0
    ! K3(:,:,:,:,:) = 0.d0
    ! G1_0 = 0.d0
    ! G1_t = 0.d0
    ! G2_0 = 0.d0

end subroutine allocate_arrays


! Deallocates all relevant arrays
subroutine deallocate_arrays()

    use variables
    implicit none

    ! LAPACK WORK ARRAY
    deallocate(work)

    ! TRAJECTORY ARRAYS
    deallocate( xn )
    deallocate( pn )
    deallocate( XE )
    deallocate( PE )

    ! BATH ARRAYS
    deallocate( coeff )
    deallocate( omega )

    ! POTENTIAL AND FORCE MATRICES
    deallocate( V )
    deallocate( G )
    deallocate( G0 )

    ! POPULATION ARRAYS
    deallocate( Cpop )
    deallocate( Cimp )
    deallocate( Cdot )
    deallocate( Cddot )
    deallocate( pop_0 )
    deallocate( pop_t )
    deallocate( coh_0 )
    deallocate( coh_t )
    deallocate( Qop_0 )
    deallocate( Qop_t )

    ! deallocate( Epop )
    ! deallocate( Eimp )

    ! PROJECTION FREE INPUT ARRAYS

    ! deallocate( F1 )
    ! deallocate( F2 )
    ! deallocate( K1 )
    ! deallocate( K3 )

end subroutine deallocate_arrays


! Sets bath frequencies and coupling coefficients
subroutine spectral_density()

    use variables
    implicit none
    integer :: i
    double precision :: eta, omega0

    eta = 0.5d0 * pi * kondo

    if ( discretize == "craig" ) then
        ! Ohmic bath, discretization a la Craig&Mano
        do i = 1,F
            omega(i) = -omegac * log( (i-0.5d0) / dble(F) )
            coeff(i) = omega(i) * dsqrt( (2 * eta * omegac) / (pi * F) )
        end do
    else if ( discretize == "makri") then
        ! Ohmic bath, discretization a la Ellen
        omega0 = omegac/dble(F) * ( 1 - exp(-omegamax/omegac) )
        do i = 1,F
            omega(i) = -omegac * log( 1 - i*omega0/omegac )
            coeff(i) = - dsqrt( kondo * omega0 ) * omega(i) ! Note factor of -1!
        end do
    else
        write(6,*) "ERROR: Discretization method must be either ",&
                   "'craig' or 'makri'!"
    end if

end subroutine spectral_density

! Samples nuclear positions and momenta
subroutine sample_nuclear()

    use variables
    implicit none
    integer :: i
    double precision :: tv, r1, r2, xn_stdev, pn_stdev

    do i = 1,F
        tv = tanh(0.5d0 * beta * omega(i))
        xn_stdev = dsqrt( 1.d0 / (2.d0 * omega(i) * tv) ) 
        pn_stdev = dsqrt( omega(i) / (2.d0 * tv) )
        call gauss(r1, r2)
        xn(i) = xn_stdev * r1
        pn(i) = pn_stdev * r2
    end do

end subroutine sample_nuclear


! Samples mapping variable positions and momenta
subroutine sample_electronic()

    use variables
    implicit none
    integer :: i
    double precision :: r1, r2, XE_stdev, PE_stdev

    if ( electronic == "phi" ) then
        XE_stdev = 1.d0 / dsqrt(2.d0) 
        PE_stdev = 1.d0 / dsqrt(2.d0) 
        do i = 1,S
            call gauss(r1, r2)
            XE(i) = XE_stdev * r1
            PE(i) = PE_stdev * r2
        end do
    else if ( electronic == "phi2" ) then
        XE_stdev = 1.d0 / 2.d0
        PE_stdev = 1.d0 / 2.d0 
        do i = 1,S
            call gauss(r1, r2)
            XE(i) = XE_stdev * r1
            PE(i) = PE_stdev * r2
        end do
    else
        write(6,*) "ERROR: electronic_sampling must be either 'phi' or 'phi2'!"
        stop
    end if

end subroutine sample_electronic


! Calculates time-zero operators
subroutine time_zero_ops()

    use variables
    implicit none
    integer :: i,j
    double precision :: zpe

    ! COHERENCES
    coh_0(2) = dcmplx(XE(1)*XE(2) + PE(1)*PE(2), XE(1)*PE(2) - PE(1)*XE(2)) / 2.d0
    coh_0(1) = dcmplx(XE(2)*XE(1) + PE(2)*PE(1), XE(2)*PE(1) - PE(2)*XE(1)) / 2.d0

    ! TRADITIONAL POPULATION OPERATORS
    if ( Aop == "seo" ) then
        zpe = 0.5d0
    else if ( Aop == "wigner" ) then
        zpe = 1.d0
    else
        write(6,*) "ERROR: A-operator type must be 'seo' or 'wigner'!"
        stop
    end if

    do i = 1,S
        pop_0(i) = 0.5d0 * ( XE(i)**2 + PE(i)**2 - zpe )
    end do

    ! IMPROVED POPULATION OPERATORS
    do i = 1,S
        Qop_0(i) = 0.5d0 * dble(S-1) * ( XE(i)**2 + PE(i)**2 )
        do j = 1,S
            if ( j /= i ) then
                Qop_0(i) = Qop_0(i) - 0.5d0 * ( XE(j)**2 + PE(j)**2 )
            end if
        end do
    end do

    ! PFI BATH FUNCTIONS
    ! do i = 1,F
    !     G1_0 = G1_0 - coeff(i) * xn(i)
    !     G2_0 = G2_0 + eye * coeff(i) * pn(i) * &
    !            tanh(beta * omega(i)/2.d0) / omega(i)
    ! end do

    ! CALCULATE NORMS
    ! if ( Aop == "seo" .and. Bop == "seo" ) then
    !     pop_norm = 16.d0
    ! else if ( Aop == "wigner" .and. Bop == "wigner" ) then
    !     write(6,*) "ERROR: Having both the A- and B-operator be of type ",&
    !                "'wigner' does not make sense! At least one operator ",&
    !                "must be projected onto onto the SEO subspace!"
    !     stop
    ! else
    !     pop_norm = 4.d0
    ! endif

    ! if ( electronic == "phi" ) then
    !     imp_norm = 4.d0
    ! else if ( electronic == "phi2" ) then
    !     imp_norm = 16.d0
    ! end if

end subroutine time_zero_ops


! Calculates time-t opeartors
subroutine time_t_ops()

    use variables
    implicit none
    integer :: i,j
    double precision :: zpe

    ! COHERENCES
    coh_t(2) = dcmplx(XE(1)*XE(2) + PE(1)*PE(2), XE(1)*PE(2) - PE(1)*XE(2)) / 2.d0
    coh_t(1) = dcmplx(XE(2)*XE(1) + PE(2)*PE(1), XE(2)*PE(1) - PE(2)*XE(1)) / 2.d0

    ! TRADITIONAL POPULATION OPERATORS
    if ( Bop == "seo" ) then
        zpe = 0.5d0
    else if ( Bop == "wigner" ) then
        zpe = 1.d0
    else
        write(6,*) "ERROR: B-operator type must be 'seo' or 'wigner'!"
        stop
    end if

    do i = 1,S
        pop_t(i) = 0.5d0 * ( XE(i)**2 + PE(i)**2 - zpe )
    end do

    ! IMPROVED POPULATION OPERATORS
    do i = 1,S
        Qop_t(i) = 0.5d0 * dble(S-1) * ( XE(i)**2 + PE(i)**2 )
        do j = 1,S
            if ( j /= i ) then
                Qop_t(i) = Qop_t(i) - 0.5d0 * ( XE(j)**2 + PE(j)**2 )
            end if
        end do
    end do

    ! PFI BATH FUNCTIONS
    ! do i = 1,F
    !     G1_t = G1_t - coeff(i) * xn(i)
    ! end do

end subroutine time_t_ops


! Accumulates observables
subroutine accumulate_obs(ts)

    use variables
    implicit none
    integer :: i,j,k,m
    integer, intent(in) :: ts
    double precision :: temp, np, norm

    ! CALCULATE NORMS
    if ( Aop == "seo" .and. Bop == "seo" ) then
        norm = 16.d0
    else if ( Aop == "wigner" .and. Bop == "wigner" ) then
        write(6,*) "ERROR: Having both the A- and B-operator be of type ",&
                   "'wigner' does not make sense! At least one operator ",&
                   "must be projected onto onto the SEO subspace!"
        stop
    else
        norm = 4.d0
    endif

    ! USE TIME 0 VALUES IF AT TIMESTEP 0 (FORTRAN: 0=1)
    if ( ts == 1 ) then
        coh_t(:) = coh_0(:)
        pop_t(:) = pop_0(:)
        Qop_t(:) = Qop_0(:)
    end if

    ! TRADITIONAL POPULATION OPERATORS
    do i = 1,S
        do j = 1,S
            Cpop(ts,i,j) = Cpop(ts,i,j) + norm * pop_0(i) * pop_t(j)
        end do
    end do

    ! IMPROVED POPULATION OPERATORS
    do i = 1,S
        do j = 1,S
            Cimp(ts,i,j) = Cimp(ts,i,j) + &
            ( S + norm * Qop_t(j) + norm * Qop_0(i)*Qop_t(j) ) / &
            dble(S**2)
        end do
    end do


    ! DIFFERENTIAL POPULATION ACF
    do i = 1,S
        Cdot(ts,i,1) = Cdot(ts,i,1) + norm * pop_0(i) * 2.d0*delta*AIMAG(coh_t(2))
        Cdot(ts,i,2) = Cdot(ts,i,2) + norm * pop_0(i) * 2.d0*delta*AIMAG(coh_t(1))

        Cddot(ts,i,1) = Cddot(ts,i,1) + norm * pop_0(i) * -delta * (&
                        -4.d0*epsilon*REAL(coh_t(2)) &
                        +2.d0*delta*(pop_t(1)-pop_t(2)) &
                        -4.d0*sum(coeff*xn*REAL(coh_t(2))) )

        Cddot(ts,i,2) = Cddot(ts,i,2) + norm * pop_0(i) * -delta * (&
                        +4.d0*epsilon*REAL(coh_t(1)) &
                        +2.d0*delta*(pop_t(2)-pop_t(1)) &
                        +4.d0*sum(coeff*xn*REAL(coh_t(1))) )
    end do

    ! Bath energy
    ! do i = 1, S
    !     do j = 1, F
    !         np = 0.5d0 * (pn(i)**2 + xn(i)**2 * omega(i)**2)
            
    !         Epop(ts,i) = Epop(ts,i) + norm * pop_0(i) * sum(pop_t) * np
            
    !         ! Improved: Unity Mapping
    !         Eimp(ts,i) = Eimp(ts,i) + norm * (1.d0 + Qop_0(i))/dble(S) * np

    !         ! Improved: Expand-Improve-Expand-LSC(1)
    !         ! temp = 0.d0
    !         ! do k = 1,S
    !         !     do m = 1,S
    !         !         ! Tr[rho phi e^iHt |m><m| N e^-iHt]
    !         !         temp = temp + pop_norm * pop_t(m) * np 
    !         !     end do
    !         !     ! Tr[rho 1 e^iHt Q_k N e^-iHt]
    !         !     temp = temp + pop_norm * Qop_t(k) * np
    !         !     ! Tr[rho Q_j e^iHt 1 N e^-iHt]
    !         !     temp = temp + pop_norm * Qop_0(i) * np
    !         !     ! Tr[rho Q_j e^iHt Q_k N e^-iHt]
    !         !     temp = temp + pop_norm * Qop_0(i) * Qop_t(k) * np
    !         ! end do
    !         ! Nimp(ts,i) = Nimp(ts,i) + temp/dble(S*S)
    !     end do
    ! end do

end subroutine accumulate_obs


! Averages and outputs observable arrays
subroutine average_obs()

    use variables
    implicit none
    integer :: i,j,k,a,b,c,d
    character(len=80) :: F1name, F2name, K1name, K3name
    character(len=120) :: fmt

    write(6,"(//'AVERAGING OBSERVABLES:')")

    ! TRADITIONAL POPULATION OPERATORS
    open(11, file="Cpop.out", action="write", status="unknown")
    do i = 1,tsteps+1
        Cpop(i,:,:) = Cpop(i,:,:) / dble(ntraj)
        write(11,'(F10.4,2x,4(ES13.6,2x))') &
        dble(i-1)*dt, Cpop(i,1,1), Cpop(i,1,2), Cpop(i,2,1), Cpop(i,2,2)
    end do
    close(11)
    write(6,"('- Saved population autocorrelation functions to Cpop.out')")

    ! IMPROVED POPULATION OPERATORS
    open(11, file="Cimp.out", action="write", status="unknown")
    do i = 1,tsteps+1
        Cimp(i,:,:) = Cimp(i,:,:) / dble(ntraj)
        write(11,'(F10.4,2x,4(ES13.6,2x))') &
        dble(i-1)*dt, Cimp(i,1,1), Cimp(i,1,2), Cimp(i,2,1), Cimp(i,2,2)
    end do
    close(11)
    write(6,"('- Saved improved population operator corr. fn. to Cimp.out')")

    ! DIFFERENTIAL POPULATION ACF
    open(11, file="Cdot.out", action="write", status="unknown")
    do i = 1,tsteps+1
        Cdot(i,:,:) = Cdot(i,:,:) / dble(ntraj)
        write(11,'(F10.4,2x,4(ES13.6,2x))') &
        dble(i-1)*dt, Cdot(i,1,1), Cdot(i,1,2), Cdot(i,2,1), Cdot(i,2,2)
    end do
    close(11)
    write(6,"('- Saved time derivative of population corr. fn. to Cdot.out')")

    open(11, file="Cddot.out", action="write", status="unknown")
    do i = 1,tsteps+1
        Cddot(i,:,:) = Cddot(i,:,:) / dble(ntraj)
        write(11,'(F10.4,2x,4(ES13.6,2x))') &
        dble(i-1)*dt, Cddot(i,1,1), Cddot(i,1,2), Cddot(i,2,1), Cddot(i,2,2)
    end do
    close(11)
    write(6,"('- Saved 2nd time derivative of population corr. fn. to Cddot.out')")

    ! open(11, file="Epop.out", status="unknown", action="write")
    ! write(fmt,'(a7,i3,a12)') "(f10.4,",S,"(2x,ES13.5))"
    ! Epop(:,:) = Epop(:,:)/dble(ntraj)
    ! do i = 1, tsteps+1
    !     write(11,fmt) (i-1) * dt,  Epop(i,1),  Epop(i,2)
    ! end do
    ! close(11)

    ! open(11, file="Eimp.out", status="unknown", action="write")
    ! write(fmt,'(a7,i3,a12)') "(f10.4,",S,"(2x,ES13.5))"
    ! Eimp(:,:) = Eimp(:,:)/dble(ntraj)
    ! do i = 1, tsteps+1
    !     write(11,fmt) (i-1) * dt, Eimp(i,1), Eimp(i,2)
    ! end do
    ! write(6,*) "- Wrote bath energy to Epop.out and Eimp.out"

    ! PROJECTION FREE INPUTS
    ! do a = 1,2   
    ! do b = 1,2
    ! do c = 1,2
    ! do d = 1,2
    !     write(F1name,'("F1_",4(i1),".out")') a,b,c,d
    !     write(F2name,'("F2_",4(i1),".out")') a,b,c,d
    !     write(K1name,'("K1_",4(i1),".out")') a,b,c,d
    !     write(K3name,'("K3_",4(i1),".out")') a,b,c,d
    !     open(11, file=F1name, action="write", status="unknown")
    !     open(12, file=F2name, action="write", status="unknown")
    !     open(13, file=K1name, action="write", status="unknown")
    !     open(14, file=K3name, action="write", status="unknown")
    !     do i = 1,tsteps+1
    !         write(11,'(F10.4,2x,ES13.6,2x,ES13.6)') &
    !         dble(i-1)*dt, F1(i,a,b,c,d) / dble(ntraj)
    !         write(12,'(F10.4,2x,ES13.6,2x,ES13.6)') &
    !         dble(i-1)*dt, F2(i,a,b,c,d) / dble(ntraj)
    !         write(13,'(F10.4,2x,ES13.6,2x,ES13.6)') &
    !         dble(i-1)*dt, K1(i,a,b,c,d) / dble(ntraj)
    !         write(14,'(F10.4,2x,ES13.6,2x,ES13.6)') &
    !         dble(i-1)*dt, K3(i,a,b,c,d) / dble(ntraj)
    !     end do
    !     close(11)
    !     close(12)
    !     close(13)
    !     close(14)
    ! end do
    ! end do
    ! end do
    ! end do
    ! write(6,"('- Saved projection free inputs to F1_abcd.out & F2_abcd.out')")
    ! write(6,"('- Saved memory kernels to K1_abcd.out & K3_abcd.out')")

end subroutine average_obs


! Calculate and accumulate GQME projection free intputs
subroutine proj_free_inputs(ts)

    use variables
    implicit none
    integer :: a,b,c,d
    integer, intent(in) :: ts

    ! USE TIME 0 VALUES IF AT TIMESTEP 0 (FORTRAN: 0=1)
    if ( ts == 1 ) then
        G1_t = G1_0
        pop_t(:) = pop_0(:)
        Qop_t(:) = Qop_0(:)
    end if

    do a = 1,S
        do b = 1,S
            do c = 1,S
                do d = 1,S
                    if ( a /= b ) then
                        if ( c == d ) then
                            F1(ts,a,b,c,d) = F1(ts,a,b,c,d) + &
                            (-1.d0)**(a-1) * (-1.d0)**(c-1) * 4.d0 * &
                            pop_0(c) * eye * G2_0 * (epsilon + G1_t) * coh_t(b)

                            F2(ts,a,b,c,d) = F2(ts,a,b,c,d) + &
                            (-1.d0)**(a-1) * 2.d0 * &
                            pop_0(c) * (epsilon + G1_t) * coh_t(b)

                            K1(ts,a,b,c,d) = K1(ts,a,b,c,d) + &
                            (-1.d0)**(a-1) * (-1.d0)**(c-1) * 4.d0 * &
                            pop_0(c) * eye * G2_0 * G1_t * coh_t(b)

                            K3(ts,a,b,c,d) = K3(ts,a,b,c,d) + &
                            (-1.d0)**(c-1) * 2.d0 * pop_0(c) * &
                            eye * G2_0 * coh_t(b)    
                        else
                            F1(ts,a,b,c,d) = F1(ts,a,b,c,d) + &
                            (-1.d0)**(a-1) * (-1.d0)**(c-1) * 4.d0 * coh_0(c) *&
                            (epsilon + G1_0) * (epsilon + G1_t) * coh_t(b)

                            F2(ts,a,b,c,d) = F2(ts,a,b,c,d) + &
                            (-1.d0)**(a-1) * 2.d0 * &
                            coh_0(c) * (epsilon + G1_t) * coh_t(b)

                            K1(ts,a,b,c,d) = K1(ts,a,b,c,d) + &
                            (-1.d0)**(a-1) * (-1.d0)**(c-1) * 4.d0 * &
                            coh_0(c) * G1_0 * G1_t * coh_t(b)

                            K3(ts,a,b,c,d) = K3(ts,a,b,c,d) + &
                            (-1.d0)**(c-1) * 2.d0 * coh_0(c) * &
                            G1_0 * coh_t(b)
                        end if
                    else
                        if ( c == d ) then
                            K3(ts,a,b,c,d) = K3(ts,a,b,c,d) + &
                            (-1.d0)**(c-1) * 2.d0 * pop_0(c) * &
                            eye * G2_0 * pop_t(b)
                        else 
                            K3(ts,a,b,c,d) = K3(ts,a,b,c,d) + &
                            (-1.d0)**(c-1) * 2.d0 * coh_0(c) * &
                            G1_0 * pop_t(b)
                        end if
                    end if
                end do
            end do
        end do
    end do

end subroutine

! Makes a single trajectory step using velocity verlet
subroutine step_vverlet()

    use variables
    implicit none
    integer :: i,j
    double precision :: hdt, qdt

    hdt = 0.5d0*dt
    qdt = 0.5d0*hdt

    call potential_force()

    ! HALF STEP IN NUCLEAR MOMENTA
    do i = 1,F
        pn(i) = pn(i) - hdt * G0(i)
        do j = 1,S
            pn(i) = pn(i) - qdt * (G(i,j,j)*(XE(j)**2 + PE(j)**2))
        end do
    end do

    ! HALF STEP IN MAPPING MOMENTA
    PE = PE - hdt * matmul(V, XE)

    ! FULL STEP IN NUCLEAR POSITIONS
    do i = 1,F
        xn(i) = xn(i) + dt * pn(i)
    end do

    ! FULL STEP IN MAPPING POSITIONS
    call potential_force()
    XE = XE + dt * matmul(V, PE)

    ! HALF STEP IN MAPPING MOMENTA
    call potential_force()
    PE = PE - hdt * matmul(V, XE)

    ! HALF STEP IN NUCLEAR MOMENTA
    do i = 1,F
        pn(i) = pn(i) - hdt * G0(i)
        do j = 1,S
            pn(i) = pn(i) - qdt * (G(i,j,j)*(XE(j)**2 + PE(j)**2))
        end do
    end do

end subroutine step_vverlet


! Makes a single trajectory step using diagonalisation of the potential matrix
subroutine step_diag

    use variables
    implicit none
    integer :: i,j
    double complex :: propagator(S) 
    double precision :: hdt, qdt, eval(S), evec(S,S), XEU(S), PEU(S)

    hdt = 0.5d0*dt
    qdt = 0.5d0*hdt

    call potential_force()

    ! HALF STEP IN MAPPING VARIABLES
    evec(:,:) = V(:,:)
    call dsyev("V", "U", S, evec, S, eval, work, lenwork, info)
    XEU = matmul(XE, evec)
    PEU = matmul(PE, evec)
    propagator = exp(-eye * hdt * eval) * dcmplx(XEU, PEU)
    XE = matmul(evec, real(propagator))
    PE = matmul(evec, aimag(propagator))

    ! HALF STEP IN NUCLEAR MOMENTA
    do i = 1,F
        pn(i) = pn(i) - hdt * G0(i)
        do j = 1,S
            pn(i) = pn(i) - qdt * (G(i,j,j)*(XE(j)**2 + PE(j)**2))
        end do
    end do

    ! FULL STEP IN NUCLEAR POSITIONS
    do i = 1,F
        xn(i) = xn(i) + dt * pn(i)
    end do

    ! RECALCULATE POTENTIAL AND FORCE WITH NEW NUCLEAR POSITIONS
    call potential_force()

    ! HALF STEP IN NUCLEAR MOMENTA
    do i = 1,F
        pn(i) = pn(i) - hdt * G0(i)
        do j = 1,S
            pn(i) = pn(i) - qdt * (G(i,j,j)*(XE(j)**2 + PE(j)**2))
        end do
    end do

    ! HALF STEP IN MAPPING VARIABLES
    evec(:,:) = V(:,:)
    call dsyev("V", "U", S, evec, S, eval, work, lenwork, info)
    XEU = matmul(XE, evec)
    PEU = matmul(PE, evec)
    propagator = exp(-eye * hdt * eval) * dcmplx(XEU, PEU)
    XE = matmul(evec, real(propagator))
    PE = matmul(evec, aimag(propagator))

end subroutine step_diag


! Calculates potential energy matrices
subroutine potential_force()

    use variables
    implicit none
    integer :: i,j
    double precision :: tr, trace

    ! STATE-INDEPENDENT POTENTIAL AND FORCE
    V0 = 0.d0
    do i = 1,F
        V0 = V0 + 0.5d0 * omega(i)**2 * xn(i)**2
        G0(i) = omega(i)**2 * xn(i)
    end do

    ! POTENTIAL ENERGY MATRIX AND FORCE TENSOR
    V(1,1) = epsilon
    V(1,2) = delta
    V(2,1) = delta
    V(2,2) = -epsilon
    G(:,:,:) = 0.d0
    do i = 1,F
        V(1,1) = V(1,1) + coeff(i) * xn(i)
        V(2,2) = V(2,2) - coeff(i) * xn(i)
        G(i,1,1) = coeff(i)
        G(i,2,2) = -coeff(i)
    end do

    ! SHIFT TRACE OF V and G to V0 and G0
    tr = trace(V,S)/dble(S)
    V0 = V0 + tr
    do i = 1,S
        V(i,i) = V(i,i) - tr
    end do

    do i = 1,F
        tr = trace(G(i,:,:),S)/dble(S)
        G0(i) = G0(i) + tr
        do j = 1,S
            G(i,j,j) = G(i,j,j) - tr
        end do
    end do

end subroutine potential_force


! Returns two numbers sampled from the standard normal distribution
! which are obtained via the Box-Muller transform of a RANLUX uniform call
subroutine gauss(r1, r2)

    use variables
    implicit none
    real(4) :: yield(2)
    double precision, intent(out) :: r1, r2

    call RANLUX(yield, 2)

    r1 = dsqrt(-2.d0 * log(yield(1))) * cos(2.d0 * pi * yield(2))
    r2 = dsqrt(-2.d0 * log(yield(1))) * sin(2.d0 * pi * yield(2))

end subroutine gauss


! Calculates the trace of a square matrix
double precision function trace(A,D)

    implicit none
    integer :: i
    integer, intent(in) :: D
    double precision, intent(in) :: A(D,D)
    
    trace = 0.d0

    do i = 1, D
        trace = trace + A(i,i)
    end do

end function trace
