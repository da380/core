program arrayresize
  implicit none
  integer, allocatable :: array1(:)
  integer, pointer :: array2(:,:)

  ! allocate and initialize array1
  allocate(array1(6))
  array1 = (/1,2,3,4,5,6/)

  ! This starts out initialized to 2
  print *, 'array1(2) = ', array1(2)

  ! Point array2 to same data as array1. The shape of array2
  ! is passed in as an array of intergers because C_F_POINTER
  ! uses and array of intergers as a SIZE parameter.
  array2 => getArray(array1, (/2,3/))

  ! Change the value at array2(2,1) (same as array1(2))
  array2(2,1) = 5

  ! Show that data in array1(2) was modified by changing
  ! array2(2,1)
  print *, 'array(2,1) = array1(2) = ', array1(2)

contains

  function getArray(array, shape_) result(aptr)
    use iso_c_binding, only: C_LOC, C_F_POINTER
    ! Pass in the array as an array of fixed size so that there
    ! is no array descriptor associated with it. This means we
    ! can get a pointer to the location of the data using C_LOC
    integer, target :: array(1)
    integer :: shape_(:)
    integer, pointer :: aptr(:,:)

    ! Use C_LOC to get the start location of the array data, and
    ! use C_F_POINTER to turn this into a fortran pointer (aptr).
    ! Note that we need to specify the shape of the pointer using an
    ! integer array.
    call C_F_POINTER(C_LOC(array), aptr, shape_)
  end function
end program
