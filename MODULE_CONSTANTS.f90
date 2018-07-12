!=================================================================================================
!> CARTESIAN QUADTREE ADAPTIVE MESH REFINEMENT LIBRARY
!> AUTHOR: VAN-DAT THANG
!> E-MAIL: datthangva@gmail.com
!> E-MAIL: vandatthang@gmail.com
!> SOURCE CODE LINK: https://github.com/dattv/QTAdaptive
!================================================================================================= 
MODULE MODULE_CONSTANTS
    
    use MODULE_PRECISION
    
    private
    
    public  :: zero     , &
                one     , &
                two     , &
                three   , &
                four    , &
                five    , &
                 six    , &
                seven   , &
                eight   , &
                nine    , &
                half    , &
               third    , &
              fourth    , & 
               fifth    , & 
               sixth    , & 
           two_third    , & 
          four_third    , & 
        three_fourth    , & 
             twelfth    , & 
    one_twentyfourth     
    
    public  :: MPI
    public  :: UNDERFINED_VALUE
    public  :: tolerance
    
    real(rp), parameter :: zero = 0._rp         , &
                            one = 1._rp         , &
                            two = 2._rp         , &
                          three = 3._rp         , &
                           four = 4._rp         , &
                           five = 5._rp         , &
                            six = 6._rp         , &
                          seven = 7._rp         , &
                          eight = 8._rp         , &
                           nine = 9._rp         , &
                           half = 0.5_rp        , &
                          third = 1.0_rp/ 3.0_rp, &
                         fourth = 1.0_rp/ 4.0_rp, &
                          fifth = 1.0_rp/ 5.0_rp, &
                          sixth = 1.0_rp/ 6.0_rp, &
                      two_third = 2.0_rp/ 3.0_rp, &
                     four_third = 4.0_rp/ 3.0_rp, &
                   three_fourth = 3.0_rp/ 4.0_rp, &
                        twelfth = 1.0_rp/12.0_rp, &
               one_twentyfourth = 1.0_rp/24.0_rp
    
    real(rp), parameter     ::              MPI = four*atan(one)
    integer(ip), parameter  :: UNDERFINED_VALUE = -9999
    
    real(rp), parameter     ::        tolerance = 1.e-8_rp
    
END MODULE MODULE_CONSTANTS    