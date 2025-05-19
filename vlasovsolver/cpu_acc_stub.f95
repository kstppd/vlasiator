module cpu_acceleration
  use iso_c_binding
  implicit none
  interface
    integer(c_int) function printf(fmt) bind(C, name="printf")
      import :: c_char, c_int
      character(kind=c_char), intent(in) :: fmt(*)
    end function printf
  end interface

  contains

  subroutine cpu_acc_dense() bind(C, name="cpu_acc_dense")
    use iso_c_binding
    implicit none
  ! character(kind=c_char,len=*), parameter :: msg = &
           ! 'Hello from Vlasiator X Fortran!' // char(10, kind=c_char) // c_null_char
  ! integer(c_int) :: istat
  ! istat=printf(msg)
  end subroutine cpu_acc_dense

end module cpu_acceleration
