module variables

    implicit none

    ! CONSTANTS
    double precision, parameter :: pi = 3.14159265359d0
    double complex, parameter :: eye = dcmplx(0.d0, 1.d0)

    ! PARALLEL PARAMETERS
    integer :: threads                          ! Number of max OMP threads

    ! INPUT PARAMETERS
    integer :: F                                ! Nuclear DoFs
    integer :: S                                ! Electronic states
    integer :: ntraj                            ! Number of trajectories
    integer :: tsteps                           ! Number of time-steps
    double precision :: dt                      ! Time-step duration
    double precision :: beta                    ! Inverse temperature
    double precision :: kondo                   ! Bath kondo parameter
    double precision :: delta                   ! Constant electronic coupling
    double precision :: omegac                  ! Bath cutoff frequency
    double precision :: omegamax                ! Maximum frequency for Makri
    double precision :: epsilon                 ! Electronic energy bias
    character(len=6) :: Aop                     ! Electronic sampling type
    character(len=6) :: Bop                     ! Electronic sampling type
    character(len=4) :: electronic              ! Electronic sampling type
    character(len=14):: intgt                   ! Type of integrator to use
    character(len=5) :: discretize              ! Type of discretization to use

    ! TRAJECTORY VARIABLES
    double precision, allocatable :: xn(:)      ! Nuclear position
    double precision, allocatable :: pn(:)      ! Nuclear momentum
    double precision, allocatable :: XE(:)      ! Mapping variable position
    double precision, allocatable :: PE(:)      ! Mapping variable momentum

    ! POTENTIAL AND FORCE MATRICES
    double precision :: V0                       ! State-independent potential
    double precision, allocatable :: V(:,:)      ! Potential energy matrix
    double precision, allocatable :: G0(:)       ! State-independent force
    double precision, allocatable :: G(:,:,:)    ! Force tensor

    ! BATH VARIABLES
    double precision, allocatable :: omega(:)   ! Bath frequencies
    double precision, allocatable :: coeff(:)   ! Bath coupling coefficients

    ! POPULATION OPERATORS AND OBSERVABLES
    double precision :: pop_norm                ! Normalization for trad. pops.
    double precision :: imp_norm                ! Normalization for imp. ops.
    double precision, allocatable :: pop_0(:)   ! Time-zero populations
    double precision, allocatable :: pop_t(:)   ! Time-t populations
    double complex, allocatable :: coh_0(:)     ! Time-zero coherences
    double complex, allocatable :: coh_t(:)     ! Time-t coherences
    double precision, allocatable :: Qop_0(:)   ! Time-zero improved operator
    double precision, allocatable :: Qop_t(:)   ! Time-t improved operator
    double precision, allocatable :: Cpop(:,:,:)! Population correlation fn.
    double precision, allocatable :: Cimp(:,:,:)! Improved operator corr. fn.
    double precision, allocatable :: Cdot(:,:,:)! Improved operator corr. fn.
    double precision, allocatable :: Cddot(:,:,:)! Improved operator corr. fn.

    ! LAPACK PARAMETERS
    integer :: info, lenwork                    ! Integer parameters for LAPACK
    double precision, allocatable :: work(:)    ! Work array for LAPACK

    ! TIMING PARAMETERS
    double precision :: start, end, start1, end1! Dummy variables
    double precision :: t_total                 ! Total time
    double precision :: t_traj                  ! Trajectory time
    double precision :: t_ops                   ! Operator calculation time
    double precision :: t_eom                   ! Equations of motion time

    ! PROJECTION FREE INPUTS
    double complex, allocatable :: F1(:,:,:,:,:)! F1 input
    double complex, allocatable :: F2(:,:,:,:,:)! F2 input
    double complex, allocatable :: K1(:,:,:,:,:)! K1 auxiliary kernel
    double complex, allocatable :: K3(:,:,:,:,:)! K3 auxiliary kernel
    double precision :: G1_0                    ! T.-dep. bath input fns. t=0
    double precision :: G1_t                    ! T.-dep. bath input fns. t=t
    double precision :: G2_0                    ! T.-indep. bath input fns.

    ! Photon nunber
    double precision, allocatable :: Epop(:,:)
    double precision, allocatable :: Eimp(:,:)

end module variables
