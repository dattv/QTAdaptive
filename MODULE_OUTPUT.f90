MODULE MODULE_OUTPUT
    
    use MODULE_PRECISION
    use MODULE_CONSTANTS
    use MODULE_UTILITY
    use MODULE_QUADTREE
    
    contains
!==================================================================================================  
    subroutine output_2D(fileName, nelm, tree)
    implicit none
    
    ! variables
    character(len = *), intent(in)                          :: fileName
    integer(ip), intent(in)                                 :: nelm
    type(quadtree), dimension(1:nelm), intent(in), target   :: tree
    integer(ip) :: i, k, count_elm, count_node
    
    real(rp), dimensioN(:,:), allocatable                   :: temp_node
    integer(ip), dimension(:,:), allocatable                :: temp_elm
    real(rp), dimension(:,:), allocatable                   :: temp_cfd_data
    integer(ip)                                             :: iUnit, NQ
    character(len = 80)                                     :: first_level_fileName
    
    ! COUNT TOTAL NODE AND TOTAL ELEMENT OF GRID SYSTEM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     count_elm = 0
    count_node = 0
    do i = 1, nelm
        if (i == 100) then 
            continue
        end if
        
        call count_quadtree(tree(i), count_elm, count_node)
    end do
    
    ! GET DATA <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    NQ = tree(1)%nq
    allocate(temp_node(2,count_node))
    allocate(temp_elm(4,count_elm))
    allocate(temp_cfd_data(tree(1)%nq,count_elm))
    count_node = 1
     count_elm = 1
    do i = 1, nelm
        call fill_elm_node(tree(i), count_node, temp_node, count_elm, temp_elm, temp_cfd_data)        
    end do
    
    ! WRITE TO FILE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
    iUnit = getUnit(100)
    open(unit = iUnit, file = trim(fileName), action = "write", status = "unknown")
        print*, "START WRITING DATA TO FILE: ", TRIM(FILENAME)
        ! START WRITING <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        write(iUnit,*) 'TITLE = "GRID"'
        write(iUnit,'(100a)') 'VARIABLES = "X","Y","RHO1_A1","RHO2_A2","U","V","P","A1"'
        write(iUnit,*) 'ZONE N =',count_node-1,',E =',count_elm-1, ',DATAPACKING=BLOCK',',ZONETYPE=FEQUADRILATERAL',',  VARLOCATION=([3-8]=CELLCENTERED)'
    
        ! WRITE X COORDINATE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        do i = 1, count_node - 1
            write(iUnit, 20     ) temp_node(1, i)
        end do
        ! WRITE Y COORDINATE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        do i = 1, count_node - 1
            write(iUnit, 20     ) temp_node(2, i)
        end do        
        ! WRITE SOLUTION <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        do k = 1, NQ
            do i = 1, count_elm - 1 
                write(iUnit, 20     ) temp_cfd_data(k,i)
            end do
        end do        
        ! WRITE ELEMENT CONNECTIVIEY <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        do i = 1, count_elm - 1
            write(iUnit, 10        ) temp_elm(1, i), temp_elm(2, i), temp_elm(3, i), temp_elm(4, i)
        end do
        ! END FILE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        print*, "END WRITING TO FILE"
    close(unit = iUnit)
    
    ! >>> WRITE FIRST LEVEL
    iUnit = getUnit(100);   first_level_fileName = "FIRST_LEVEL_"//trim(fileName) 
    open(unit = iUnit, file = first_level_fileName, action = "write", status = "unknown")
        print*, "START WRITING DATA TO FILE: ", TRIM(first_level_fileName)
        ! START WRITING <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        write(iUnit,*       ) 'TITLE = "GRID"'
        write(iUnit,'(100a)') 'VARIABLES = "X","Y","RHO1_A1","RHO2_A2","U","V","P","A1"'
        write(iUnit,*       ) 'ZONE N =',nelm*4,',E =',nelm, ',DATAPACKING=BLOCK',',ZONETYPE=FEQUADRILATERAL',',  VARLOCATION=([3-8]=CELLCENTERED)'
    
        ! WRITE X COORDINATE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        do i = 1, nelm
            write(iUnit, 20     ) tree(i)%pts(1)%coord(1)
            write(iUnit, 20     ) tree(i)%pts(2)%coord(1)
            write(iUnit, 20     ) tree(i)%pts(3)%coord(1)
            write(iUnit, 20     ) tree(i)%pts(4)%coord(1)
        end do
        ! WRITE Y COORDINATE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        do i = 1, nelm
            write(iUnit, 20     ) tree(i)%pts(1)%coord(2)
            write(iUnit, 20     ) tree(i)%pts(2)%coord(2)
            write(iUnit, 20     ) tree(i)%pts(3)%coord(2)
            write(iUnit, 20     ) tree(i)%pts(4)%coord(2)
        end do        
        ! WRITE SOLUTION <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        do k = 1, NQ
            do i = 1, nelm 
                write(iUnit, 20     ) tree(i)%w(k)
            end do
        end do        
        ! WRITE ELEMENT CONNECTIVIEY <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        do i = 1, nelm
            write(iUnit, 10        ) (i-1)*4 + 1, (i-1)*4 + 2, (i-1)*4 + 3, (i-1)*4 + 4
        end do
        ! END FILE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        print*, "END WRITING TO FILE"
    close(unit = iUnit)
    continue
    deallocate(temp_node)
    deallocate(temp_elm)
    deallocate(temp_cfd_data)
    return
10  format(4(i8,x)) 
20  format(f22.13)    
    end subroutine output_2D
!================================================================================================== 
    recursive subroutine count_quadtree(tree, count_elm, count_node)
    implicit none
    type(quadTree), intent(in), target  :: tree
    integer(ip), intent(inout)          :: count_elm, count_node
    type(quadtree), pointer             :: p_tree
    
    if (.not. tree%is_leaf) then 
        p_tree => tree%north_west
        call count_quadtree(p_tree, count_elm, count_node)
        
        p_tree => tree%north_east
        call count_quadtree(p_tree, count_elm, count_node)
        
        p_tree => tree%south_west
        call count_quadtree(p_tree, count_elm, count_node)
        
        p_tree => tree%south_east
        call count_quadtree(p_tree, count_elm, count_node)        
    else
         count_elm = count_elm + 1
        count_node = count_node + 4
    end if
    
    return
    end subroutine count_quadtree
!==================================================================================================  
    recursive subroutine fill_elm_node(tree, count_node, node, count_elm, elm, cfd_data)
    implicit none
    type(quadtree),                     intent(in),     target  :: tree
    integer(ip),                        intent(inout)           :: count_node, count_elm
    real(rp), dimension(:,:),           intent(inout)           :: node
    integer(ip), dimension(:,:),        intent(inout)           :: elm
    real(rp), dimension(:,:),           intent(inout)           :: cfd_data
    
    ! variables
    type(quadtree),                                     pointer :: p_tree
    
    ! body 
    if (.not. tree%is_leaf) then 
        p_tree => tree%north_west
        call fill_elm_node(p_tree, count_node, node, count_elm, elm, cfd_data)
        
        p_tree => tree%north_east
        call fill_elm_node(p_tree, count_node, node, count_elm, elm, cfd_data)
        
        p_tree => tree%south_east
        call fill_elm_node(p_tree, count_node, node, count_elm, elm, cfd_data)
        
        p_tree => tree%south_west
        call fill_elm_node(p_tree, count_node, node, count_elm, elm, cfd_data)   
    else
        
        node(:,count_node  ) = tree%pts(1)%coord(:)
        node(:,count_node+1) = tree%pts(2)%coord(:)
        node(:,count_node+2) = tree%pts(3)%coord(:)
        node(:,count_node+3) = tree%pts(4)%coord(:)

        elm(1,count_elm) = count_node + 0
        elm(2,count_elm) = count_node + 1
        elm(3,count_elm) = count_node + 2
        elm(4,count_elm) = count_node + 3
        
        cfd_data(:,count_elm) = tree%w(:)
        count_node = count_node + 4
         count_elm = count_elm + 1
        
    end if
    
    return
    end subroutine fill_elm_node
!==================================================================================================    
!==================================================================================================    
!==================================================================================================    
!==================================================================================================    
END MODULE MODULE_OUTPUT    