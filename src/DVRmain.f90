program DVR

  use omp_lib
!  use ReadInput, only : ReadInputMain

  integer ( kind = 4 ) id
  real ( kind = 8 ) wtime
! real, dimension(:,:), allocatable :: myVec

! allocate(myVec(3,2))    

! myVec=1.

! myVec(100,3)=10.

! write(*,*) myVec(1,1)
! write(*,*) myVec(100,4)

!  wtime = omp_get_wtime ( )

! write ( *, '(a,i8)' ) &
!   '  The number of processors available = ', omp_get_num_procs ( )
! write ( *, '(a,i8)' ) &
!   '  The number of threads available    = ', omp_get_max_threads ( )

!  INSIDE THE PARALLEL REGION, have each thread say hello.
!
!$omp parallel &
!$omp private ( id )
! id = omp_get_thread_num ( )

! write ( *, '(a,i8,a,i8)' ) '  HELLO from process ', id

!$omp end parallel
!
  call DVRCore()

!  Finish up by measuring the elapsed time.
!
!  wtime = omp_get_wtime ( ) - wtime

!write ( *, '(a)' ) ' '
!write ( *, '(a)' ) '  Back OUTSIDE the parallel region.'
!write ( *, '(a)' ) ' '
! write ( *, '(a,g14.6)' ) '  Elapsed wall clock time = ', wtime

end program DVR
