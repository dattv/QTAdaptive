MODULE MODULE_AMR
    use MODULE_PRECISION
    use MODULE_CONSTANTS
    use MODULE_QUADTREE
    use MODULE_GENERICMETHOD
    
    integer(ip)                         ::         LEVEL = 1
    integer(ip)                         ::     AMR_INDEX = 1
    real(rp)                            :: AMR_THRESHOLD
    contains
!==================================================================================================
    subroutine AMR_whole_domain(first, last, tree)
    implicit none
    integer(ip), intent(in)                                         :: first, last
    type(quadtree), dimension(first:last), intent(inout), target    :: tree
    
    ! >>> ADAPTIVE MESH [REFINE] REFINEMENT LEVEL 1 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   
    LEVEL = 1;  AMR_INDEX = 6; AMR_THRESHOLD = 0.3_rp
    call loop_on_quadtree_array(first, last, tree, AMR_finner_single_cell)   
    
    ! >>> ADAPTIVE MESH [REFINE] REFINEMENT LEVEL 2 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    LEVEL = 2;  AMR_INDEX = 6; AMR_THRESHOLD = 0.5_rp
    call loop_on_quadtree_array(first, last, tree, AMR_finner_single_cell)

    ! >>> ADAPTIVE MESH [REFINE] REFINEMENT LEVEL 3 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    LEVEL = 3;  AMR_INDEX = 6; AMR_THRESHOLD = 0.7_rp
    call loop_on_quadtree_array(first, last, tree, AMR_finner_single_cell)    
    
    ! >>> ADAPTIVE MESH [COARSER] REFINEMENT LEVEL 3 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    LEVEL = 3;  AMR_INDEX = 6; AMR_THRESHOLD = 0.7_rp
    call AMR_coarser_loop_on_quadtree_array(first, last, tree)
    
    ! >>> ADAPTIVE MESH [COARSER] REFINEMENT LEVEL 2 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    LEVEL = 2;  AMR_INDEX = 6; AMR_THRESHOLD = 0.5_rp
    call AMR_coarser_loop_on_quadtree_array(first, last, tree)
    
    ! >>> ADAPTIVE MESH [COARSER] REFINEMENT LEVEL 1 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    LEVEL = 1;  AMR_INDEX = 6; AMR_THRESHOLD = 0.3_rp    
    call AMR_coarser_loop_on_quadtree_array(first, last, tree)
    
    return
    end subroutine AMR_whole_domain
!==================================================================================================  
    subroutine AMR_finner_single_cell(tree)
    implicit none
    type(quadTree), pointer, intent(inout)  :: tree
    integer(ip)                             :: i
    real(rp)                                :: max_grad

    max_grad = compute_max_gradient(tree)
    
    if (max_grad >= AMR_THRESHOLD) then
        i = LEVEL
        call AMR_refiner_level(i, tree)             
    end if
    
    return
    end subroutine  AMR_finner_single_cell
!==================================================================================================
    subroutine AMR_refiner_level(level, tree)
    implicit none
    integer(ip), intent(in)         :: level
    type(quadtree), intent(inout)   :: tree
    
    if (tree%m_level == level - 1) then 
        call tree%subdivide()
    end if
    return
    end subroutine  AMR_refiner_level
!==================================================================================================
    subroutine AMR_coarser_loop_on_quadtree_array(first, last, tree)
    implicit none
    integer(ip), intent(in)                                         :: first, last
    type(quadtree), dimension(first:last), intent(inout), target    :: tree
    integer(ip)                                                     :: i
    type(quadtree), pointer                                         :: p_cell
    
    do i = first, last
        p_cell => tree(i)
        call AMR_coarser_loop_on_single_quadtree(p_cell)
    end do
    
    return
    end subroutine AMR_coarser_loop_on_quadtree_array
!==================================================================================================
    recursive subroutine AMR_coarser_loop_on_single_quadtree(tree)
    implicit none
    type(quadtree), intent(inout), target   :: tree
    type(quadtree), pointer                 :: p_cell 
    real(rp)                                :: max_grad_nw, max_grad_ne, max_grad_sw, max_grad_se
    
    if (.not. tree%is_leaf) then 
        if (associated(tree%north_west)) call AMR_coarser_loop_on_single_quadtree(tree%north_west)    
        if (associated(tree%north_east)) call AMR_coarser_loop_on_single_quadtree(tree%north_east)    
        if (associated(tree%south_east)) call AMR_coarser_loop_on_single_quadtree(tree%south_east)    
        if (associated(tree%south_west)) call AMR_coarser_loop_on_single_quadtree(tree%south_west)
       
    else
        if (.not. associated(tree%father)) return
        p_cell => tree%father
        ! ===> CHECK CELL TO REMOVE <=============================================================
        max_grad_nw = compute_max_gradient(p_cell%north_west)
        max_grad_ne = compute_max_gradient(p_cell%north_east)
        max_grad_se = compute_max_gradient(p_cell%south_east)
        max_grad_sw = compute_max_gradient(p_cell%south_west)
        
        if (max(max_grad_nw, max_grad_ne, max_grad_se, max_grad_sw) < AMR_THRESHOLD .and. p_cell%m_level == LEVEL - 1) then 
            
            ! CHANGE ADJOINT ELEMENT OF CELL NORTH WEST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (associated(p_cell%north_west%adj_north)) p_cell%north_west%adj_north%adj_south => p_cell 
            if (.not. p_cell%north_west%adj_north%is_leaf) then 
                p_cell%north_west%adj_north%south_west%adj_south => p_cell
                p_cell%north_west%adj_north%south_east%adj_south => p_cell
            end if
            
            if (associated(p_cell%north_west%adj_west))    p_cell%north_west%adj_west%adj_east => p_cell
            if (.not. p_cell%north_west%adj_west%is_leaf) then 
                p_cell%north_west%adj_west%north_east%adj_east => p_cell
                p_cell%north_west%adj_west%south_east%adj_east => p_cell
            end if
            
            call p_cell%north_west%delete()     
            
            ! CHANGE ADJOINT ELEMENT OF CELL NORTH EAST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (associated(p_cell%north_east%adj_north)) p_cell%north_east%adj_north%adj_south => p_cell
            if (.not. p_cell%north_east%adj_north%is_leaf) then 
                p_cell%north_east%adj_north%south_west%adj_south => p_cell                
                p_cell%north_east%adj_north%south_east%adj_south => p_cell                
            end if
            
            if (associated(p_cell%north_east%adj_east))    p_cell%north_east%adj_east%adj_west => p_cell
            if (.not. p_cell%north_east%adj_east%is_leaf) then
                      p_cell%north_east%adj_east%north_west%adj_west => p_cell        
                      p_cell%north_east%adj_east%south_west%adj_west => p_cell        
            end if
            
            call p_cell%north_east%delete()        
            
            ! CHANGE ADJOINT ELEMENT OF CELL SOUTH WEST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (associated(p_cell%south_west%adj_west))    p_cell%south_west%adj_west%adj_east => p_cell
            if (.not. p_cell%south_west%adj_west%is_leaf) then
                p_cell%south_west%adj_west%north_east%adj_east => p_cell      
                p_cell%south_west%adj_west%south_east%adj_east => p_cell      
            end if
            
            if (associated(p_cell%south_west%adj_south)) p_cell%south_west%adj_south%adj_north => p_cell
            if (.not. p_cell%south_west%adj_south%is_leaf) then
                      p_cell%south_west%adj_south%north_west%adj_north => p_cell   
                      p_cell%south_west%adj_south%north_east%adj_north => p_cell   
            end if
            
            call p_cell%south_west%delete() 
            
            ! CHANGE ADJOINT ELEMENT OF CELL SOUTH EAST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (associated(p_cell%south_east%adj_south)) p_cell%south_east%adj_south%adj_north => p_cell
            if (.not. p_cell%south_east%adj_south%is_leaf) then 
                p_cell%south_east%adj_south%north_east%adj_north => p_cell  
                p_cell%south_east%adj_south%north_west%adj_north => p_cell  
            end if
            
            if (associated(p_cell%south_east%adj_east))    p_cell%south_east%adj_east%adj_west => p_cell
            if (.not. p_cell%south_east%adj_east%is_leaf) then 
                p_cell%south_east%adj_east%north_west%adj_west => p_cell      
                p_cell%south_east%adj_east%south_west%adj_west => p_cell      
            end if
 
            call p_cell%south_east%delete() 
            
            deallocate(p_cell%north_west)
            deallocate(p_cell%north_east)
            deallocate(p_cell%south_east)
            deallocate(p_cell%south_west)
            p_cell%is_leaf = .true.
            continue          
        end if       
    end if
    
    return
    end subroutine AMR_coarser_loop_on_single_quadtree
!==================================================================================================
    subroutine AMR_coarser_single_cell(tree)
    implicit none
    type(quadTree), pointer, intent(inout)  :: tree
    integer(ip)                             :: i
    real(rp)                                :: max_grad

    max_grad = compute_max_gradient(tree)
    
    if (max_grad < AMR_THRESHOLD .and. tree%m_level >= LEVEL) then
        i = LEVEL
        call AMR_coarser_level(i, tree)           
    else
        
    end if
    
    continue
    
    return
    end subroutine  AMR_coarser_single_cell
!================================================================================================== 
    subroutine AMR_coarser_level(level, tree)
    implicit none
    integer(ip), intent(in)         :: level
    type(quadtree), intent(inout)   :: tree
    type(quadtree), pointer         :: temp_cell, father
    
    ! CELL NORTH
    if (associated(tree%adj_north)) tree%adj_north%adj_south => tree%father
    if (associated(tree%adj_east))  tree%adj_east%adj_west => tree%father
    if (associated(tree%adj_south)) tree%adj_south%adj_north => tree%father    
    if (associated(tree%adj_west))  tree%adj_west%adj_east => tree%father    
    
    call tree%delete()
    
    father => tree%father
    
    tree%index = UNDERFINED_VALUE
    
    if (father%north_west%index == UNDERFINED_VALUE) nullify(father%north_west)
    if (father%north_east%index == UNDERFINED_VALUE) nullify(father%north_east)
    if (father%south_west%index == UNDERFINED_VALUE) nullify(father%south_west)
    if (father%south_east%index == UNDERFINED_VALUE) nullify(father%south_east)
    
    return
    end subroutine AMR_coarser_level
!================================================================================================== 
    function compute_max_gradient(tree) result(res)
    implicit none
    type(quadtree), intent(in), target  :: tree
    real(rp)                            :: res
    
    type(quadTree), pointer             :: tree_n, tree_e, tree_s, tree_w
    type(quadtree), pointer             :: tree_nw, tree_ne, tree_sw, tree_se
    real(rp)                            :: grad_n, grad_e, grad_s, grad_w
    real(rp)                            :: grad_nw, grad_ne, grad_sw, grad_se
    
    ! COMPUTE GRADIENT OF CELL
    if (associated(tree%adj_north)) then 
        grad_n = abs(tree%w(AMR_INDEX) - tree%adj_north%w(AMR_INDEX)   )
    else
        grad_n = zero
    end if
    
    if (associated(tree%adj_west)) then 
        grad_w = abs(tree%w(AMR_INDEX) - tree%adj_west%w(AMR_INDEX)    )
    else
        grad_w = zero
    end if
    
    if (associated(tree%adj_south)) then 
        grad_s = abs(tree%w(AMR_INDEX) - tree%adj_south%w(AMR_INDEX)   )
    else
        grad_s = zero
    end if
    
    if (associated(tree%adj_east)) then 
        grad_e = abs(tree%w(AMR_INDEX) - tree%adj_east%w(AMR_INDEX)    )
    else
        grad_e = zero
    end if
    
    if (associated(tree%adj_north)) then 
        if (associated(tree%adj_north%adj_west)) then 
            grad_nw = abs(tree%w(AMR_INDEX)  - tree%adj_north%adj_west%w(AMR_INDEX))
        else
            grad_nw = zero
        end if
    else
        grad_nw = zero
    end if
    
    if (associated(tree%adj_north)) then 
        if (associated(tree%adj_north%adj_east)) then 
            grad_ne = abs(tree%w(AMR_INDEX)  - tree%adj_north%adj_east%w(AMR_INDEX))
        else
            grad_ne = zero
        end if
    else
        grad_ne = zero
    end if
    
    if (associated(tree%adj_south)) then 
        if (associated(tree%adj_south%adj_east)) then 
            grad_se = abs(tree%w(AMR_INDEX)  - tree%adj_south%adj_east%w(AMR_INDEX))
        else 
            grad_se = zero
        end if
    else
        grad_se = zero
    end if
    
    if (associated(tree%adj_south)) then 
        if (associated(tree%adj_south%adj_west)) then 
            grad_sw = abs(tree%w(AMR_INDEX)  - tree%adj_south%adj_west%w(AMR_INDEX))
        else
            grad_sw = zero
        end if
    else
        grad_sw = zero
    end if
    
    res = max( grad_n , &
               grad_w , &
               grad_s , &
               grad_e , &
               grad_nw, &
               grad_ne, &
               grad_se, &
               grad_sw  )
    return
    end function compute_max_gradient
!==================================================================================================

END MODULE MODULE_AMR    