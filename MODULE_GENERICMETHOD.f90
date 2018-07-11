!=================================================================================================
!> CARTESIAN QUADTREE ADAPTIVE MESH REFINEMENT LIBRARY
!> AUTHOR: VAN-DAT THANG
!> E-MAIL: datthangva@gmail.com
!> E-MAIL: vandatthang@gmaill.com
!> SOURCE CODE LINK: https://github.com/dattv/QTAdaptive
!=================================================================================================  
MODULE MODULE_GENERICMETHOD
    
    use MODULE_QUADTREE
    
    contains
!==================================================================================================    
    subroutine loop_on_quadtree_array(first, last, tree, workFunction)
    implicit none
    integer(ip), intent(in) :: first, last
    type(quadtree), dimension(:), intent(inout), target :: tree
    
    interface   
        subroutine workFunction(work)
            use MODULE_QUADTREE
            type(quadtree), pointer, intent(inout)  :: work
        end subroutine workFunction
    end interface
    integer(ip) :: i        
    
    do i = first, last, 1
        call loop_on_quadtree_single_level(tree(i), workFunction)        
    end do
    
    return
    end subroutine loop_on_quadtree_array
!==================================================================================================
    recursive subroutine loop_on_quadtree_single_level(tree, workFunction)
    implicit none
    type(quadTree), intent(inout), target   :: tree
    
    interface   
        subroutine workFunction(work)
            use MODULE_QUADTREE
            type(quadtree), pointer, intent(inout)  :: work
        end subroutine workFunction
    end interface
    
    type(quadTree), pointer :: work
    
    if (.not. tree%is_leaf) then 
        call loop_on_quadtree_single_level(tree%north_west, workFunction)        
        call loop_on_quadtree_single_level(tree%north_east, workFunction)        
        call loop_on_quadtree_single_level(tree%south_east, workFunction)        
        call loop_on_quadtree_single_level(tree%south_west, workFunction)        
    else
        work => tree
        call workFunction(work)
    end if
    
    end subroutine loop_on_quadtree_single_level
END MODULE MODULE_GENERICMETHOD    