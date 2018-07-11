!=================================================================================================
!> CARTESIAN QUADTREE ADAPTIVE MESH REFINEMENT LIBRARY
!> AUTHOR: VAN-DAT THANG
!> E-MAIL: datthangva@gmail.com
!> E-MAIL: vandatthang@gmail.com
!> SOURCE CODE LINK: https://github.com/dattv/QTAdaptive
!=================================================================================================  
    MODULE MODULE_CFD_DATA
    
        use MODULE_PRECISION
        use MODULE_CONSTANTS
        
        type :: reconstruction_data
        	real(rp)    :: x_l
        	real(rp)    :: x_r
        	real(rp)    :: y_l
        	real(rp)    :: y_r
        end type reconstruction_data
        
        type :: cfd_data
            integer(ip)                         :: NQ       ! total equations
        	real(rp), dimension(:), allocatable :: u        ! conservative variables
        	real(rp), dimension(:), allocatable :: w        ! primative variables
        	real(rp), dimension(:), allocatable :: u_old    ! conservative variables
        	real(rp), dimension(:), allocatable :: u_new    ! conservative variables
        	real(rp), dimension(:), allocatable :: res      ! residual
        	real(rp)                            :: wsn      ! wave speed
            
            ! ===> reconstruction data <=====================================
            type(reconstruction_data), dimension(:), allocatable    :: recons
            
        contains
        procedure   ::    new => new_cfd_data
        procedure   :: delete => delete_cfd_data
        
        procedure   :: add_cfd_data
        generic     :: operator (+) => add_cfd_data
        
        end type cfd_data

!=================== INTERFACE =====================
        interface assignment(=)
            module procedure asign_cfd
        end interface 
        
    contains
!================================================================================================= 
    subroutine new_cfd_data(this, NQ)
    implicit none
    class(cfd_data), intent(inout)  :: this
    integer(ip), intent(in)         :: NQ
    integer(ip)                     :: temp_NQ
    
    ! body
    this%NQ = NQ
    if (.not. allocated(this%u)) then 
        allocate(this%u(NQ))
    else
        temp_NQ = sizeof(this%u)
        if (temp_NQ /= NQ) then 
            deallocate(this%u)
            allocate(this%u(NQ))
        end if
    end if
    
    if (.not. allocated(this%w)) then 
        allocate(this%w(NQ))
    else
        temp_NQ = sizeof(this%w)
        if (temp_NQ /= NQ) then 
            deallocate(this%w)
            allocate(this%w(NQ))
        end if
    end if
    
    if (.not. allocated(this%u_old)) then 
        allocate(this%u_old(NQ))
    else
        temp_NQ = sizeof(this%u_old)
        if (temp_NQ /= NQ) then 
            deallocate(this%u_old)
            allocate(this%u_old(NQ))
        end if
    end if
    
    if (.not. allocated(this%u_new)) then 
        allocate(this%u_new(NQ))
    else
        temp_NQ = sizeof(this%u_new)
        if (temp_NQ /= NQ) then 
            deallocate(this%u_new)
            allocate(this%u_new(NQ))
        end if
    end if
    
    if (.not. allocated(this%res)) then 
        allocate(this%res(NQ))
    else
        temp_NQ = sizeof(this%res)
        if (temp_NQ /= NQ) then 
            deallocate(this%res)
            allocate(this%res(NQ))
        end if
    end if
    
    if (.not. allocated(this%recons)) then 
        allocate(this%recons(NQ))
    else
        temp_NQ = sizeof(this%recons)
        if (temp_NQ /= NQ) then 
            deallocate(this%recons)
            allocate(this%recons(NQ))
        end if
    end if     
    
    this%u     = zero
    this%w     = zero
    this%u_old = zero
    this%u_new = zero
    this%res   = zero
    this%wsn   = zero
    
    this%recons%x_l = zero   
    this%recons%x_r = zero   
    this%recons%y_l = zero   
    this%recons%y_r = zero   
    return
    end subroutine new_cfd_data
!================================================================================================= 
    subroutine delete_cfd_data(this)
    implicit none
    class(cfd_data), intent(inout)  :: this
    
    ! body
    if (allocated(this%u     )) deallocate(this%u     )  ! conservative variables
    if (allocated(this%w     )) deallocate(this%w     )  ! primative variables
    if (allocated(this%u_old )) deallocate(this%u_old )  ! conservative variables
    if (allocated(this%u_new )) deallocate(this%u_new )  ! conservative variables
    if (allocated(this%res   )) deallocate(this%res   )  ! residual  
    if (allocated(this%recons)) deallocate(this%recons)
    return
    end subroutine delete_cfd_data
!=================================================================================================
    subroutine asign_cfd(out_cfd_data, in_cfd_data)
    implicit none
    type(cfd_data), intent(out)  :: out_cfd_data
    type(cfd_data), intent(in)   :: in_cfd_data
    
    ! body
    call out_cfd_data%new(in_cfd_data%NQ)
    
    out_cfd_data%NQ     = in_cfd_data%NQ       ! total equations
    out_cfd_data%u      = in_cfd_data%u        ! conservative variables
    out_cfd_data%w      = in_cfd_data%w        ! primative variables
    out_cfd_data%u_old  = in_cfd_data%u_old    ! conservative variables
    out_cfd_data%u_new  = in_cfd_data%u_new    ! conservative variables
    out_cfd_data%res    = in_cfd_data%res      ! residual
    out_cfd_data%wsn    = in_cfd_data%wsn      ! wave speed
    
    out_cfd_data%recons%x_l = in_cfd_data%recons%x_l
    out_cfd_data%recons%x_r = in_cfd_data%recons%x_r
    out_cfd_data%recons%y_l = in_cfd_data%recons%y_l
    out_cfd_data%recons%y_r = in_cfd_data%recons%y_r
    return
    end subroutine  asign_cfd
!=================================================================================================  
    function add_cfd_data(this, in_cfd_data) result(res)
    implicit none
    class(cfd_data), intent(in) :: this
    type(cfd_data), intent(in)  :: in_cfd_data
    type(cfd_data)              :: res
    
    call res%new(in_cfd_data%NQ)
    
    res%NQ     = in_cfd_data%NQ       ! total equations
    
    res%u      = this%u     + in_cfd_data%u        ! conservative variables
    res%w      = this%w     + in_cfd_data%w        ! primative variables
    res%u_old  = this%u_old + in_cfd_data%u_old    ! conservative variables
    res%u_new  = this%u_new + in_cfd_data%u_new    ! conservative variables
    res%res    = this%res   + in_cfd_data%res      ! residual
    res%wsn    = this%wsn   + in_cfd_data%wsn      ! wave speed
    
    res%recons%x_l = this%recons%x_l + in_cfd_data%recons%x_l
    res%recons%x_r = this%recons%x_r + in_cfd_data%recons%x_r
    res%recons%y_l = this%recons%y_l + in_cfd_data%recons%y_l
    res%recons%y_r = this%recons%y_r + in_cfd_data%recons%y_r
    
    return
    end function add_cfd_data
!=================================================================================================  
!=================================================================================================  
!=================================================================================================  
!=================================================================================================  
    
    END MODULE MODULE_CFD_DATA