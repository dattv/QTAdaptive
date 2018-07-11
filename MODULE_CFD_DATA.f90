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
        procedure   :: new => new_cfd_data
        procedure   :: delete => delete_cfd_data
        end type cfd_data
        
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
    return
    end subroutine delete_cfd_data
!=================================================================================================  
!=================================================================================================  
!=================================================================================================  
!=================================================================================================  
!=================================================================================================  
!=================================================================================================  
    
    END MODULE MODULE_CFD_DATA