!=================================================================================================
!> CARTESIAN QUADTREE ADAPTIVE MESH REFINEMENT LIBRARY
!> AUTHOR: VAN-DAT THANG
!> E-MAIL: datthangva@gmail.com
!> E-MAIL: vandatthang@gmail.com
!> SOURCE CODE LINK: https://github.com/dattv/QTAdaptive
!=================================================================================================  
MODULE MODULE_PRECISION
    
    private
    
    public  :: ip, rp
    
    integer, parameter  :: sp = kind(1.d0)
    integer, parameter  :: p2 = selected_real_kind(2*precision(1._sp))
    integer, parameter  :: i5 = selected_int_kind(5)
    integer, parameter  :: i15 = selected_int_kind(15)
    integer, parameter  :: ip =i5
    integer, parameter  :: rp = p2
    
END MODULE MODULE_PRECISION    