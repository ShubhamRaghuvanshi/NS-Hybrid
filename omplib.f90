module omplib
use omp_lib
use mod_serial_fluid 
public
contains
!! ---------------------------
! export OMP_NUM_THREADS=4
subroutine get_total_threads(nthreads)
integer :: nthreads
nthreads = OMP_GET_NUM_THREADS() 
endsubroutine get_total_threads
!! ---------------------------
subroutine get_thread_id(id)
integer :: id
id = omp_get_thread_num() 
endsubroutine get_thread_id 
!! ---------------------------
subroutine wall_clock_time(time)
  double precision :: time
  time=omp_get_wtime()
endsubroutine wall_clock_time
!! ---------------------------
endmodule omplib
