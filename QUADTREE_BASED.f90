!  QUADTREE_BASED.f90 
!
!  FUNCTIONS:
!  QUADTREE_BASED - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: QUADTREE_BASED
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program QUADTREE_BASED

    use MODULE_PRECISION
    use MODULE_NODECOORD
    use MODULE_CONSTANTS
    use MODULE_GENERICMETHOD
    use MODULE_MAKEGRID
    use MODULE_QUADTREE
    use MODULE_INITIALCONDITION
    use MODULE_AMR
    use MODULE_OUTPUT
    
    implicit none

    ! Variables
    integer(ip) :: nelm
    integer(ip) :: NQ
    type(quadtree), dimension(:), allocatable, target   :: tree
    type(quadtree), pointer :: p
    
    ! Body of QUADTREE_BASED
    nelm = 900; NQ = 6
    
    call make_grid_2d(nelm, (/zero, one/), (/zero, one/), NQ, tree)
    
    call initial_2D_simle_check_grid(nelm, tree)
    
    ! ===> TEST ADAPTIVE MESH REFINEMENT METHOD <================================================== 
    call AMR_whole_domain(1, nelm, tree)
    ! ===> END TEST ADAPTIVE MESH REFINEMENT METHOD <==============================================
    
    ! ===> PRINT OUT INITIAL CONDITION <===========================================================
    call output_2d("INITIAL.TEC", nelm, tree)
    ! ===> END PRINT OUT INITIAL CONDITION <=======================================================    
    
    ! ===> CHANGE CFD DATA ALITTLE TO CHECK ARM FEATURE <==========================================    
    call loop_on_quadtree_array(1, nelm, tree, initial_2d_02_simple)
    call update_back_whole_domain(1, nelm, tree)
    ! ===> END <===================================================================================
    
    ! ===> TEST AMR <==============================================================================    
    call AMR_whole_domain(1, nelm, tree)
    ! ===> END <===================================================================================
    
    ! ===> PRINT OUT NEW GRID <====================================================================
    call output_2d("AMR_MESH.TEC", nelm, tree)
    ! ===> END PRINT OUT <=========================================================================
    
    print *, 'Hello World'

    end program QUADTREE_BASED

