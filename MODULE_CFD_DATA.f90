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
        
        type :: cfd_data
            integer(ip)                         :: NQ       ! total equations
        	real(rp), dimension(:), allocatable :: u        ! conservative variables
        	real(rp), dimension(:), allocatable :: w        ! primative variables
        	real(rp), dimension(:), allocatable :: u_old    ! conservative variables
        	real(rp), dimension(:), allocatable :: u_new    ! conservative variables
        	real(rp), dimension(:), allocatable :: res      ! residual
        	real(rp)                            :: wsn      ! wave speed
            
        contains
        procedure   ::    new => new_cfd_data
        procedure   :: delete => delete_cfd_data
        
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
        allocate(this%u(NQ))
    else
        temp_NQ = sizeof(this%res)
        if (temp_NQ /= NQ) then 
            deallocate(this%res)
            allocate(this%res(NQ))
        end if
    end if    
    
    return
    end subroutine new_cfd_data
!================================================================================================= 
    subroutine delete_cfd_data(this)
    implicit none
    class(cfd_data), intent(inout)  :: this
    
    ! body
    if (allocated(this%u    )) deallocate(this%u    )  ! conservative variables
    if (allocated(this%w    )) deallocate(this%w    )  ! primative variables
    if (allocated(this%u_old)) deallocate(this%u_old)  ! conservative variables
    if (allocated(this%u_new)) deallocate(this%u_new)  ! conservative variables
    if (allocated(this%res  )) deallocate(this%res  )  ! residual    
    return
    end subroutine delete_cfd_data
!=================================================================================================
    subroutine asign_cfd(out_cfd_data, in_cfd_data)
    implicit none
    type(cfd_data), intent(out)  :: out_cfd_data
    type(cfd_data), intent(in)   :: in_cfd_data
    
    call out_cfd_data%new(in_cfd_data%NQ)
    
    out_cfd_data%NQ     = in_cfd_data%NQ       ! total equations
    out_cfd_data%u      = in_cfd_data%u        ! conservative variables
    out_cfd_data%w      = in_cfd_data%w        ! primative variables
    out_cfd_data%u_old  = in_cfd_data%u_old    ! conservative variables
    out_cfd_data%u_new  = in_cfd_data%u_new    ! conservative variables
    out_cfd_data%res    = in_cfd_data%res      ! residual
    out_cfd_data%wsn    = in_cfd_data%wsn      ! wave speed
    return
    end subroutine  asign_cfd
!=================================================================================================  
!=================================================================================================  
!=================================================================================================  
!=================================================================================================  
!=================================================================================================  
    
    END MODULE MODULE_CFD_DATA