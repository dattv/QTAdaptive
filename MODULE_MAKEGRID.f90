MODULE MODULE_MAKEGRID
    
    use MODULE_PRECISION
    use MODULE_CONSTANTS
    use MODULE_NODECOORD
    use MODULE_QUADTREE
    
    contains
!==================================================================================================
    subroutine make_grid_2d(nelm, x, y, nq, tree)
    implicit none
    
    ! I/O variables
    integer(ip),                                intent(in)              :: nelm
    real(rp), dimension(2),                     intent(in)              :: x, y
    integer(ip),                                intent(in)              :: nq
    type(quadtree), dimension(:), allocatable,  intent(inout), target   :: tree
    type(nodecoord), dimension(5)                                       :: points
    real(rp)                                                            :: dx
    
    ! variables
    integer(ip)                                                         :: nx, ny, i, j, k, id
    
    ! body 
    nx = int(sqrt(real(nelm, rp)))
    ny = nx
    dx = abs(x(2) - x(1))/real(nx,rp)
    
    do i = 1, 5
        call points(i)%new(2, (/zero, zero/))
    end do
    
    allocate(tree(nelm))
    do j = 1, ny
        do i = 1, nx
            id = (j-1)*nx + i
            
            points(1)%coord(1) = (i-1)*dx
            points(1)%coord(2) = (j  )*dx
            
            points(2)%coord(1) = (i  )*dx
            points(2)%coord(2) = (j  )*dx
            
            points(3)%coord(1) = (i  )*dx
            points(3)%coord(2) = (j-1)*dx
            
            points(4)%coord(1) = (i-1)*dx
            points(4)%coord(2) = (j-1)*dx
            
            points(5)%coord(1) = (i  )*dx - dx*half
            points(5)%coord(2) = (j  )*dx - dx*half
            
            call tree(id)%new(id, 0, half*dx, points, nq)
            
            if (i > 1 .and. i < nx .and. j > 1 .and. j < ny) then 
                tree(id)%adj_north => tree((j  )*nx + i     )
                tree(id)%adj_east  => tree((j-1)*nx + i + 1 )
                tree(id)%adj_south => tree((j-2)*nx + i     )
                 tree(id)%adj_west => tree((j-1)*nx + i - 1 )
            else if (i == 1 .and.  j > 1 .and. j < ny) then 
                tree(id)%adj_north => tree((j  )*nx + i     )
                tree(id)%adj_east  => tree((j-1)*nx + i + 1 )
                tree(id)%adj_south => tree((j-2)*nx + i     )
                tree(id)%adj_west => NULL() 
            else if (i == nx .and. j > 1 .and. j < ny) then 
                tree(id)%adj_north => tree((j  )*nx + i     )
                tree(id)%adj_east  => NULL()
                tree(id)%adj_south => tree((j-2)*nx + i     )
                 tree(id)%adj_west => tree((j-1)*nx + i - 1 )
            else if (i > 1 .and. i < nx .and. j == 1) then 
                tree(id)%adj_north => tree((j  )*nx + i     )
                tree(id)%adj_east  => tree((j-1)*nx + i + 1 )
                tree(id)%adj_south => NULL()
                 tree(id)%adj_west => tree((j-1)*nx + i - 1 )  
            else if (i > 1 .and. i < nx .and. j == ny) then 
                tree(id)%adj_north => NULL()
                tree(id)%adj_east  => tree((j-1)*nx + i + 1 )
                tree(id)%adj_south => tree((j-2)*nx + i     )
                 tree(id)%adj_west => tree((j-1)*nx + i - 1 ) 
            else if (i == 1 .and. j == 1) then 
                tree(id)%adj_north => tree((j  )*nx + i     )
                tree(id)%adj_east  => tree((j-1)*nx + i + 1 )
                tree(id)%adj_south => NULL()
                 tree(id)%adj_west => NULL()  
            else if (i==nx .and. j == 1) then 
                tree(id)%adj_north => tree((j )*nx + i      )
                tree(id)%adj_east  => NULL()
                tree(id)%adj_south => NULL()
                tree(id)%adj_west => tree((j-1)*nx + i - 1  )  
            else if (i == nx .and. j == ny) then 
                tree(id)%adj_north => NULL()
                tree(id)%adj_east  => NULL()
                tree(id)%adj_south => tree((j-2)*nx + i     )
                 tree(id)%adj_west => tree((j-1)*nx + i - 1 )   
            else if (i == 1 .and. j == ny) then 
                tree(id)%adj_north => NULL()
                tree(id)%adj_east  => tree((j-1)*nx + i + 1 )
                tree(id)%adj_south => tree((j-1)*nx + i     )
                 tree(id)%adj_west => NULL()                
            else
            end if
            
        end do
    end do
    
    do i = 1, 5
        call points(i)%delete()
    end do    
    return
    end subroutine make_grid_2d
    
!==================================================================================================    
!==================================================================================================    
!==================================================================================================    
!==================================================================================================    
!==================================================================================================    
!==================================================================================================    
END MODULE MODULE_MAKEGRID    