!=================================================================================================
!> CARTESIAN QUADTREE ADAPTIVE MESH REFINEMENT LIBRARY
!> AUTHOR: VAN-DAT THANG
!> E-MAIL: datthangva@gmail.com
!> E-MAIL: vandatthang@gamil.com
!> SOURCE CODE LINK: https://github.com/dattv/QTAdaptive
!=================================================================================================  
    MODULE MODULE_CFD_DATA
    
        use MODULE_PRECISION
        use MODULE_CONSTANTS
        
        type :: cfd_data
            integer(ip)                         :: NQ       ! total equations
        	real(rp), dimension(:), allocatable :: u        ! conservative variables
        	real(rp), dimension(:), allocatable :: w        ! primative variables
        	real(rp), dimension(:), allocatable :: u_old    ! conservative variables
        	real(rp), dimension(:), allocatable :: u_new    ! conservative variables
        	real(rp), dimension(:), allocatable :: res      ! residual
        	real(rp)                            :: wsn      ! wave speed
        end type cfd_data
    END MODULE MODULE_CFD_DATA