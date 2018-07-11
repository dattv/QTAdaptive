!=================================================================================================
!> CARTESIAN QUADTREE ADAPTIVE MESH REFINEMENT LIBRARY
!> AUTHOR: VAN-DAT THANG
!> E-MAIL: datthangva@gmail.com
!> E-MAIL: vandatthang@gmail.com
!> SOURCE CODE LINK: https://github.com/dattv/QTAdaptive
!================================================================================================= 
MODULE MODULE_QUADTREE  
    
    use MODULE_PRECISION
    use MODULE_CONSTANTS
    use MODULE_NODECOORD
    use MODULE_CFD_DATA
    
    type :: quadtree
    	integer(ip) :: index
        integer(ip) :: m_level
        real(rp)    :: halfDim
        
        type(nodeCoord), dimension(5)   :: pts
        
        type(quadtree), pointer :: adj_west
        type(quadtree), pointer :: adj_east
        type(quadtree), pointer :: adj_north
        type(quadtree), pointer :: adj_south
        
        type(cfd_data)  :: data
        
        reaL(rp)    :: vol
        
        type(quadtree), pointer :: father
        type(quadtree), pointer :: north_west
        type(quadtree), pointer :: north_east
        type(quadtree), pointer :: south_west
        type(quadtree), pointer :: south_east
        logical :: is_leaf
        
    contains
    procedure   ::       new => new_quadtree
    procedure   ::    delete => delete_quatree
    procedure   ::    insert => insert_point_2_quadtree
    procedure   :: subdivide => subdivide_quadtree
    end type quadtree
    
    contains
!================================================================================================== 
    subroutine new_quadtree(this, index, m_level, halfDim, points, nq)
    implicit none
    class(quadtree), intent(inout), target      :: this
    integer(ip), intent(in)                     :: index, m_level
    real(rp), intent(in)                        :: halfDim
    type(nodeCoord), dimension(5), intent(in)   :: points
    integer(ip), intent(in)                     :: nq
    integer(ip)                                 :: i
    
    ! body
      this%index = index
    this%m_level = m_level
    this%halfDim = halfDim
    
    do i = 1, 5
        this%pts(i) = points(i)
    end do
    
    call this%data%new(nq)
    
    this%is_leaf = .true.
    return
    end subroutine new_quadtree
!================================================================================================== 
    subroutine delete_quatree(this)
    implicit none
    class(quadtree), intent(inout)  :: this
    integer(ip) :: i
    
    do i = 1, 5
        call this%pts(i)%delete()
    end do

    call this%data%delete()
    
    this%is_leaf = .false.

    return
    end subroutine delete_quatree
!==================================================================================================  
    subroutine subdivide_quadtree(this)
    implicit none
    class(quadtree), intent(inout), target  :: this
    real(rp)                                :: halfDim
    type(nodeCoord), dimension(5)           :: points
    integer(ip)                             :: i
    type(quadtree), pointer                 :: child_NW, child_NE, child_SE, child_SW, temp_cell
    logical                                 :: keep_going
    ! body
    
    ! CHANGE CURRENT CELL NOT LEAF CELL <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    this%is_leaf = .false.
    
    ! DECLARE FOUR CHILDS CELL OF CURRENT CELL <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    allocate(this%north_west);  this%north_west%father => this
    allocate(this%north_east);  this%north_east%father => this
    allocate(this%south_west);  this%south_west%father => this
    allocate(this%south_east);  this%south_east%father => this
    
    ! DECLARE FIVE POINT OF EACH CHILD NODE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    do i = 1, 5
        call points(i)%new(2, (/zero, zero/))
    end do
    child_NW => this%north_west
    child_NE => this%north_east
    child_SE => this%south_east
    child_SW => this%south_west
    
    halfDim = this%halfDim*half
    
    ! MAKE CHILD CELL NORTH WEST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    points(5)%coord(1) = this%pts(5)%coord(1) - halfDim
    points(5)%coord(2) = this%pts(5)%coord(2) + halfDim
    
    points(1)%coord(1) = this%pts(1)%coord(1)
    points(1)%coord(2) = this%pts(1)%coord(2)

    points(2)%coord(1) = points(5)%coord(1) + halfDim
    points(2)%coord(2) = points(5)%coord(2) + halfDim
    
    points(3)%coord(1) = this%pts(5)%coord(1)             
    points(3)%coord(2) = this%pts(5)%coord(2)
    
    points(4)%coord(1) = this%pts(1)%coord(1)
    points(4)%coord(2) = this%pts(5)%coord(2) 
    
    call child_NW%new(0, this%m_level+1, halfDim, points, this%data%nq) 
    
    child_NW%data = this%data
    
    ! >> SET ELEMENT ADJOINT OF CELL NORTH WEST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if (associated(this%adj_north)) then 
        temp_cell => this%adj_north        
        do while(.not. temp_cell%is_leaf .and. temp_cell%m_level < child_NW%m_level)
            temp_cell => temp_cell%south_West
        end do
        
        child_NW%adj_north => temp_cell
        if (temp_cell%m_level == child_NW%m_level) temp_cell%adj_south => child_NW
        ! AT THIS POINT, I DEVELOPED THIS CODE WITH AN SIMPLE ASSUME
    end if
    
     child_NW%adj_east => child_NE
    child_NW%adj_south => child_SW
    
    if (associated(this%adj_west)) then 
        temp_cell => this%adj_west
        do while(.not. temp_cell%is_leaf .and. temp_cell%m_level < child_NW%m_level) 
            temp_cell => temp_cell%north_east
        end do
        
        child_NW%adj_west => temp_cell
        if (temp_cell%m_level == child_NW%m_level) temp_cell%adj_east => child_NW
        ! AT THIS POINT, I DEVELOPED THIS CODE WITH AN SIMPLE ASSUME            
    end if
    
    ! MAKE CHILD CELL NORTH EAST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    
    points(5)%coord(1) = this%pts(5)%coord(1) + halfDim
    points(5)%coord(2) = this%pts(5)%coord(2) + halfDim
    
    points(1)%coord(1) =  child_NW%pts(2)%coord(1)
    points(1)%coord(2) =  child_NW%pts(2)%coord(2)
    
    points(2)%coord(1) = this%pts(2)%coord(1)
    points(2)%coord(2) = this%pts(2)%coord(2)
    
    points(3)%coord(1) = points(5)%coord(1) + halfDim
    points(3)%coord(2) = points(5)%coord(2) - halfDim
    
    points(4)%coord(1) = this%pts(5)%coord(1)
    points(4)%coord(2) = this%pts(5)%coord(2)
    
    call child_NE%new(0, this%m_level+1, halfDim, points, this%data%nq)

    child_NE%data = this%data
    
    ! >> SET ELEMENT ADJOINT OF CELL NORTH EAST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if (associated(this%adj_north)) then 
        temp_cell => this%adj_north
        do while(.not. temp_cell%is_leaf .and. temp_cell%m_level < child_NE%m_level)
            temp_cell => temp_cell%south_east            
        end do
        child_NE%adj_north => temp_cell
        if (temp_Cell%m_level == child_NE%m_level) temp_cell%adj_south => child_NE
    end if
    
    if (associated(this%adj_east)) then 
        temp_cell => this%adj_east
        do while(.not. temp_cell%is_leaf .and. temp_cell%m_level < child_NE%m_level)
            temp_cell => temp_cell%north_west
        end do
        child_NE%adj_east => temp_Cell
        if (temp_Cell%m_level == child_NE%m_level) temp_cell%adj_west => child_NE
    end if
    
    child_NE%adj_south => child_SE
    
    child_NE%adj_west => child_NW
        
    ! MAKE CHILD CELL SOURTH EAST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    points(5)%coord(1) = this%pts(5)%coord(1) + halfDim
    points(5)%coord(2) = this%pts(5)%coord(2) - halfDim
    
    points(1)%coord(1) = this%pts(5)%coord(1)
    points(1)%coord(2) = this%pts(5)%coord(2)
    
    points(2)%coord(1) = child_NE%pts(3)%coord(1)
    points(2)%coord(2) = child_NE%pts(3)%coord(2)
    
    points(3)%coord(1) = this%pts(3)%coord(1)
    points(3)%coord(2) = this%pts(3)%coord(2)
    
    points(4)%coord(1) = this%pts(5)%coord(1)
    points(4)%coord(2) = points(5)%coord(2) - halfDim
    
    call child_SE%new(0, this%m_level+1, halfDim, points, this%data%nq)
    
    child_SE%data = this%data
    
    ! >> SET ELEMENT ADJOINT OF CELL SOUTH EAST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    child_SE%adj_north => child_NE
    
    if (associated(this%adj_east)) then 
        temp_cell => this%adj_east
        do while(.not. temp_cell%is_leaf .and. temp_cell%m_level < child_SE%m_level)
            temp_cell => temp_cell%south_west
        end do
        child_SE%adj_east => temp_cell
        if (temp_Cell%m_level == child_SE%m_level) temp_Cell%adj_west => child_SE
    end if
    
    if (associated(this%adj_south)) then 
        temp_cell => this%adj_south
        do while (.not. temp_Cell%is_leaf .and. temp_cell%m_level < child_NE%m_level)
            temp_cell => temp_cell%north_east
        end do
        child_SE%adj_south => temp_cell
        if (temp_cell%m_level == child_SE%m_level) temp_cell%adj_north => child_SE
    end if
    
    child_SE%adj_West => child_SW
    
    ! MAKE CHILD CELL SOUTH WES <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    points(5)%coord(1) = this%pts(5)%coord(1) - halfDim
    points(5)%coord(2) = this%pts(5)%coord(2) - halfDim
    
    points(1)%coord(1) = child_NW%pts(4)%coord(1)
    points(1)%coord(2) = child_NW%pts(4)%coord(2)
    
    points(2)%coord(1) = this%pts(5)%coord(1)
    points(2)%coord(2) = this%pts(5)%coord(2)
    
    points(3)%coord(1) = child_SE%pts(4)%coord(1)
    points(3)%coord(2) = child_SE%pts(4)%coord(2)
    
    points(4)%coord(1) = this%pts(4)%coord(1)
    points(4)%coord(2) = this%pts(4)%coord(2)
    
    call child_SW%new(0, this%m_level+1, halfDim, points, this%data%nq)

    child_SW%data = this%data
    
    ! >> SET ELEMENT ADJOINT OF CELL SOUTH WEST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
    child_SW%adj_north => child_NW
     child_SW%adj_east => child_SE
    
    if(associated(this%adj_south)) then 
        temp_Cell => this%adj_south
        do while( .not. temp_Cell%is_leaf .and. temp_Cell%m_level < child_SW%m_level)
            temp_Cell => temp_Cell%north_West
        end do
        child_SW%adj_south => temp_Cell
        if (temp_Cell%m_level == child_SW%m_level) temp_cell%adj_north => child_SW
    end if
    
    if (associated(this%adj_west)) then 
        temp_cell => this%adj_west
        do while(.not. temp_cell%is_leaf .and. temp_cell%m_level < child_SW%m_level)
            temp_cell => temp_cell%south_east
        end do
        child_SW%adj_west => temp_cell
        if (temp_cell%m_level == child_SW%m_level) temp_cell%adj_east => child_SW
    end if
    
    ! DELETE FIVE POINT 
    do i = 1, 5
        call points(i)%delete()
    end do
    
    return
    end subroutine subdivide_quadtree
!================================================================================================== 
    subroutine insert_point_2_quadtree(this, index, point)
    implicit none
    class(quadtree), intent(inout), target  :: this
    integer(ip), intent(in)                 :: index
    type(nodecoord), intent(in)             :: point
    
    call this%pts(index)%new(point%ndim, point%coord)
    
    return
    end subroutine insert_point_2_quadtree
!==================================================================================================  
    subroutine update_back_whole_domain(first, last, tree)
    implicit none
    integer(ip), intent(in)                                         :: first, last
    type(quadtree), dimension(first:last), intent(inout), target    :: tree
    integer(ip)                                                     :: i
    
    do i = first, last
        call update_back_single_cell(tree(i))     
    end do
    
    return
    end subroutine update_back_whole_domain    
!================================================================================================== 
    recursive subroutine update_back_single_cell(tree)
    implicit none
    type(quadtree), intent(inout), target   :: tree
    
    if (.not. tree%is_leaf) then 
        call update_back_single_cell(tree%north_west)
        call update_back_single_cell(tree%north_east)
        call update_back_single_cell(tree%south_east)
        call update_back_single_cell(tree%south_west)
        
            tree%data%u = (tree%north_west%data%u     + tree%north_east%data%u     + tree%south_east%data%u     + tree%south_west%data%u        )/four
        tree%data%u_old = (tree%north_west%data%u_old + tree%north_east%data%u_old + tree%south_east%data%u_old + tree%south_west%data%u_old    )/four
        tree%data%u_new = (tree%north_west%data%u_new + tree%north_east%data%u_new + tree%south_east%data%u_new + tree%south_west%data%u_new    )/four
            tree%data%w = (tree%north_west%data%w     + tree%north_east%data%w     + tree%south_east%data%w     + tree%south_west%data%w        )/four
          tree%data%res = (tree%north_west%data%res   + tree%north_east%data%res   + tree%south_east%data%res   + tree%south_west%data%res      )/four
          tree%data%wsn = (tree%north_west%data%wsn   + tree%north_east%data%wsn   + tree%south_east%data%wsn   + tree%south_west%data%wsn      )/four
    else
        return
    end if
    return
    end subroutine update_back_single_cell
!==================================================================================================    
!==================================================================================================    
!==================================================================================================    
    
END MODULE MODULE_QUADTREE    