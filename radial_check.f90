program radial_check

  use radial_mod

  implicit none

  real(idp), allocatable     :: pot(:)
  type(para_t)               :: para
  type(grid_t)               :: grid
  real(idp), allocatable     :: eigen_vals(:)
  real(idp), allocatable     :: matrix(:,:), matrix_full(:,:)
  real(idp), allocatable     :: matrix_all(:,:), matrix_inv_all(:,:), unity(:,:)
  real(idp), allocatable     :: Tkin_cardinal(:,:)
  integer                    :: i, j, a, b, l, l_val
  real(idp),  allocatable    :: file_r(:), file_pot(:)
  integer, allocatable       :: ipiv(:)
  real(idp)                  :: full_r_max
  real(idp), allocatable     :: work(:)
  integer                    :: info

  l_val=0

  ! IMPORTANT: It is crucial that the parameters used here for the radial
  !            matrix elements are the same that will be used when obtaining
  !            the eigenfunctions of the radial Schroedinger equation!

  para%pottype       = 'file' ! 'density_file', 'file', 'analytical'
  para%pot_filename  = 'input_pot.in' 
  para%r_min         = 0.0
  para%r_max         = 300.0
  para%m             = 200
  para%nl            = 5
  para%nr            = 1001 !nr = m * nl + 1
  para%l             = l_val !Rotational quantum number
  para%mass          = 1.0

  para%mapped_grid   = .false.
  para%maptype       = 'diff'
  para%read_envelope = ''
  para%beta          = 0.015
  para%E_max         = 1d-5

  call init_grid_dim_GLL(grid, para)

  ! Note that later on the first and final grid point are removed from the grid
  ! Since we need for the matrix elements the value of r_max of the original,
  ! full grid, it will be stored here
  full_r_max = maxval(grid%r)
 
  if (para%pottype == 'analytical') then
    write(*,*) 'Using analytical potential'
    allocate(pot(size(grid%r)))
    do i = 1,size(pot)
       pot(i) = analytical_potential(grid%r(i), para%mass)
    end do
    !open(11, file="input_pot.ana", form="formatted", &
    !&    action="write")
    !write(11,*) '# File emulating a V = - 1/r hydrogen potential'
    !do i = 1, 100000
    !  write(11,*) 300.0d0 * (real(i,idp)/100000d0),                            &
    !  &           - one / (300.0d0 * (real(i,idp)/100000d0))
    !end do
    !close(11)
  elseif ((para%pottype == 'file') .or. (para%pottype == 'density_file')) then
    write(*,*) 'Using potential from file '//trim(para%pot_filename)
    call init_grid_op_file_1d(file_pot, file_r, para)
    call map_op(grid%r, pot, file_r, file_pot) !Perform splining
  else
    write(*,*) "ERROR: Invalid pottype"
  end if
  
  call init_work_cardinalbase(Tkin_cardinal, grid, para%mass)
  call redefine_ops_cardinal(pot)
  call redefine_GLL_grid_1d(grid)
  
  do j = 1, 20 
    open(11,file="explicit_dvr_primitive"//trim(int2str(j))//".dat",form="formatted",&
    &    action="write", recl=100000)
    do i = 1, size(grid%r)
      write(11,'(3ES25.17)') grid%r(i), grid%weights(i),                       &
      &        dvr_primitive_val(j, para, grid, grid%r(i))*sqrt(grid%weights(i))
    end do
    close(11)
  end do

end program radial_check
