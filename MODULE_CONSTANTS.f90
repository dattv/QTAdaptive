!=================================================================================================
!> CARTESIAN QUADTREE ADAPTIVE MESH REFINEMENT LIBRARY
!> AUTHOR: VAN-DAT THANG
!> E-MAIL: datthangva@gmail.com
!> E-MAIL: vandatthang@gamil.com
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
                half
    
    public  :: MPI
    public  :: UNDERFINED_VALUE
    
    real(rp), parameter :: zero = 0._rp, &
                            one = 1._rp, &
                            two = 2._rp, &
                          three = 3._rp, &
                           four = 4._rp, &
                           five = 5._rp, &
                            six = 6._rp, &
                          seven = 7._rp, &
                          eight = 8._rp, &
                           nine = 9._rp, &
                           half = 0.5_rp
    
    real(rp), parameter     ::              MPI = four*atan(one)
    integer(ip), parameter  :: UNDERFINED_VALUE = -9999   
    
END MODULE MODULE_CONSTANTS    