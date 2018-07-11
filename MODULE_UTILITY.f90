!=================================================================================================
!> CARTESIAN QUADTREE ADAPTIVE MESH REFINEMENT LIBRARY
!> AUTHOR: VAN-DAT THANG
!> E-MAIL: datthangva@gmail.com
!> E-MAIL: vandatthang@gmail.com
!> SOURCE CODE LINK: https://github.com/dattv/QTAdaptive
!================================================================================================= 
MODULE MODULE_UTILITY    
    
    use MODULE_PRECISION
    use MODULE_CONSTANTS
    use MODULE_NODECOORD
    
    contains
!========================================================================================    
    function distance(point1, point2) result(dist)
    implicit none
    type(nodecoord), intent(in)    :: point1, point2
    real(rp)    :: dist, s
    integer(ip) :: nDim, i
    
    if (point1%nDim /= point2%nDim) then 
        write(*, *) "point1, point2 are not in the same dimension", point1%nDim, point2%nDim
        write(*, *) "pause"
        pause
        stop
    end if
    
    nDim = point1%nDim
    s = zero
    do i = 1, nDim
        s = s + (point2%coord(i) - point1%coord(i))**2
    end do
    
    dist = sqrt(s)
    return
    end function distance
!========================================================================================
    integer function lens(string)
    
       character(len=*) string
       
       do i=len(string),1,-1
         if( string(i:i) .ne. ' ')goto 10
       end do
       i = 0
    10 continue
    
       lens = i
        
    end function lens    
!========================================================================================
    subroutine sort(nnodes, arr, tar)
    implicit none
    integer, intent(in)    :: nnodes
    integer, dimension(nnodes) :: arr
    real(rp), dimension(nnodes)    :: tar
    integer    :: i, j, id1, id2
    
    do i = 1, nnodes
       do j = i, nnodes
           id1 = arr(i);   id2 = arr(j)
           if (tar(i) >= tar(j)) then 
               arr(i) = id2;   arr(j) = id1
           end if
       end do
    end do
    
    return
    end subroutine sort
!========================================================================================
    function getUnit(maxUnit) result(iUnit)
    implicit none
    integer, intent(in)    :: maxUnit
    integer    :: iUnit
    logical    :: opened
    
    do iUnit = 1, maxUnit
       inquire(unit = iUnit, opened = opened)
       if (.not. opened) then 
           return
       end if
    end do
    print*, "Cannot find out the adequate IO threads"
    stop
    return
    end function getUnit    
!========================================================================================      
END MODULE MODULE_UTILITY    