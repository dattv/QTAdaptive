MODULE MODULE_NODECOORD
    
    use MODULE_PRECISION
    use MODULE_CONSTANTS
    
    type :: nodeCoord
    	integer(ip)                         :: ndim
        real(rp), dimension(:), allocatable :: coord
        
    contains
    procedure   :: new => new_nodeCoord
    procedure   :: delete => delete_nodeCoord
    end type nodeCoord
    
!========================= INTERFACE ========================
    interface assignment(=)
    module procedure    asign_coord
    end interface
    
    contains
!==================================================================================================    
    subroutine new_nodeCoord(this, nDim, coord)
    implicit none
    class(nodeCoord), intent(inout)         :: this
    integer(ip), intent(in)                 :: nDim
    real(rp), dimension(nDim), intent(in)   :: coord
    
    this%ndim = ndim
    if (.not. allocated(this%coord))    allocate(this%coord(nDim))
    this%coord(:) = coord(:)
    
    return
    end subroutine new_nodeCoord
!==================================================================================================  
    subroutine delete_nodeCoord(this)
    class(nodeCoord), intent(inout)  :: this
    
    if (allocated(this%coord)) deallocate(this%coord)
    
    return
    end subroutine delete_nodeCoord
!================================================================================================== 
    subroutine asign_coord(node1, node2)
    implicit none
    type(nodeCoord), intent(out)    :: node1
    type(nodeCoord), intent(in)     :: node2
    
    call node1%new(node2%ndim, node2%coord)
    
    return
    end subroutine asign_coord
!==================================================================================================    
!==================================================================================================    
!==================================================================================================    
!==================================================================================================    
!==================================================================================================    
    
END MODULE MODULE_NODECOORD    