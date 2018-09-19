program dvr_diag 

  use dvr_diag_mod

  implicit none

  real(idp), allocatable     :: pot(:)
  type(para_t)               :: para
  type(grid_t)               :: grid
  real(idp), allocatable     :: eigen_vals(:)
  real(idp), allocatable     :: matrix(:,:)
  real(idp), allocatable     :: Tkin_cardinal(:,:)
  integer                    :: i, j, l_val, ith, jth, a, b
  real(idp)                  :: dr, curr_r, curr_r_prime, integral_val,        &
  &                             integral_val_ana, integral_val_ana2, th_1,     &
  &                             th_2, dth
  real(idp),  allocatable    :: file_r(:), file_pot(:)
  logical                    :: only_bound ! Only writes out bound states
  
  l_val=0

  para%pottype       = 'analytical' ! 'density_file', 'file', 'analytical'
  para%pot_filename  = 'input_pot.in' 
  para%r_min         = 0.0
  para%r_max         = 30.0
  para%m             = 100
  para%nl            = 5
  para%nr            = 501  !nr = m * nl + 1
  para%l             = l_val !Multipole order
  para%Z             = 1
  para%mass          = 1.0

  para%mapped_grid   = .false.
  para%maptype       = 'diff'
  para%read_envelope = ''
  para%beta          = 0.015
  para%E_max         = 1d-5

  only_bound         = .true.
        
  call init_grid_dim_GLL(grid, para) 
 
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

  ! Write potential 
  open(11, file="input_pot.out", form="formatted", &
  &    action="write")
  write(11,*) '# input potential after splining and adjusting to GLL grid'
  do i = 1, size(pot(:))
    write(11,*) grid%r(i), pot(i)
  end do
  close(11)

  !! Get banded storage format of Hamiltonian matrix in the FEM-DVR basis
  call get_real_surf_matrix_cardinal(matrix, grid, pot, Tkin_cardinal)

  !! Diagonalize Hamiltonian matrix which is stored in banded format.
  !! nev specifies the first nev eigenvalues and eigenvectors to be extracted.
  !! If needed, just add more
  call diag_arpack_real_sym_matrix(matrix, formt='banded', n=size(matrix(1,:)),&
  &                   nev=nint(0.5*para%nr), which='SA', eigenvals=eigen_vals, &
  &                   rvec=.true.)

  ! Write eigenvalues.
  open(11, file="eigenvalues_GLL.dat", form="formatted", &
  &    action="write")
  write(11,*) '#_______________________________________________________________'
  write(11,*) "#              eigenvalues for hydrogen with l = 0 "
  write(11,*) '#_______________________________________________________________'
  write(11,*) "#     index    -    eigenvalue    -    eigenvector normalization"
  write(11,*) ""
  do i = 1, size(eigen_vals(:))
    write(11,'(I8,3ES25.17)') i-1, eigen_vals(i),                              &
    &                         - one / (two * (real(i+para%l, idp))**2),        &
    &                         dot_product(matrix(:,i), matrix(:,i))
  end do
  close(11)

  write(*,*) size(grid%weights)

  ! Write eigenvectors. Here, if we want to represent the eigenvectors in the
  ! physical grid instead of numerical grid defined by the normalized FEM-BASIS,
  ! we should divide by the square root of the Gaussian interpolating weights,
  ! see eg. arXiv:1611.09034 (2016).
  ! We furthermore divide by r such that we obtain the proper radial
  ! wavefunction psi(r) instead of the rescaled object u(r) = psi(r) * r

  if (only_bound) then
    open(11, file="eigenvectors_GLL.dat", form="formatted",&
    &    action="write", recl=100000)
    do i = 1, size(matrix(:,1))
      write(11,'(ES25.17, 1x)', advance = 'No') grid%r(i)
      do j = 1, size(matrix(i,:))
        if (eigen_vals(j) > zero) cycle
        write(11,'(ES25.17, 1x)', advance = 'No')                              &
        & matrix(i,j) / (sqrt(grid%weights(i)) * grid%r(i))
      end do
      write(11,*) ' '
    write(76,'(i5,f15.8)') i, grid%r(i)
    end do
    close(11)
  else
    open(11, file="eigenvectors_GLL.dat", form="formatted",&
    &    action="write", recl=100000)
    do i = 1, size(matrix(:,1))
      write(11,*) grid%r(i),                                                   &
      & (matrix(i,j) / (sqrt(grid%weights(i)) * grid%r(i)),                    &
      & j = 1, size(matrix(i,:)))
    end do
    close(11)
  end if

  do a = 1, 3
  do b = 1, 3



  dr = 0.1d0
  write(*,*) 'Doing it for a, b: ', a, b
  integral_val_ana = zero
  do i = 1, 501
    curr_r = (i-1) * dr
    do j = 1, 501
!     if (a.ne.b) write(*,*) i, j
      curr_r_prime = (j-1) * dr
      if (curr_r + curr_r_prime < 1d-16) cycle
      integral_val_ana = integral_val_ana +                                    &
      &                  radial_hydrogen_ef(curr_r, l_val+a, l_val)**2 *       &
      &                  radial_hydrogen_ef(curr_r_prime, l_val+b, l_val)**2 * &
      &                  ( min(curr_r, curr_r_prime)**(l_val) /                &
      &                    max(curr_r, curr_r_prime)**(l_val+1) ) * dr * dr *  &
      &                  curr_r * curr_r * curr_r_prime * curr_r_prime
    end do
    write(79,'(2f20.12)') curr_r, radial_hydrogen_ef(curr_r, l_val+a, l_val)
  end do

  !do i = 1, 1001
  !  curr_r = (i-1) * dr
  !  do j = 1, 1001
  !    curr_r_prime = (j-1) * dr
  !    if (curr_r + curr_r_prime < 1d-16) cycle
  !    integral_val_ana = integral_val_ana +                                    &
  !    &                  two * exp(- curr_r) * two * exp(- curr_r) *           &
  !    &    sqrt(two)**(-3) * (two - curr_r_prime) * exp(-curr_r_prime / two) * &
  !    &    sqrt(two)**(-3) * (two - curr_r_prime) * exp(-curr_r_prime / two) * &
  !    &                  ( min(curr_r, curr_r_prime)**(l_val) /                &
  !    &                    max(curr_r, curr_r_prime)**(l_val+1) ) * dr * dr *  &
  !    &                  curr_r * curr_r * curr_r_prime * curr_r_prime
  !  end do
  !end do
  integral_val = zero
  do i = 1, para%nr-2
    curr_r = grid%r(i)
    do j = 1, para%nr-2
      curr_r_prime = grid%r(j) 
      if (curr_r + curr_r_prime < 1d-16) cycle
      integral_val = integral_val +                                            &
      &              matrix(i,a) * matrix(i,a) * matrix(j,b) * matrix(j,b) *   &
      &              ( min(curr_r, curr_r_prime)**(l_val) /                    &
      &                max(curr_r, curr_r_prime)**(l_val+1) )
!     write(80,'(4i4, 4f16.8)') a, b, i, j, matrix(i,a), matrix(j,b),( min(curr_r, curr_r_prime)**(l_val) /                    &
!           &                max(curr_r, curr_r_prime)**(l_val+1)), integral_val
    end do
  end do
! do i = 1, para%nr-2
!   curr_r = grid%r(i)
!   !write(*,*) two * exp(-curr_r), matrix(i,1)/(sqrt(grid%weights(i))*grid%r(i)), &
!   !&          radial_hydrogen_ef(curr_r, 1, 0)
!   !write(*,*) sqrt(two)**(-3) * (two - curr_r) * exp(-curr_r / two),          &
!   !&          matrix(i,2)/(sqrt(grid%weights(i))*grid%r(i)),                  &
!   !&          radial_hydrogen_ef(curr_r, 2, 0)
! end do
  write(78,'(2I4,2ES25.17)') a, b, integral_val, integral_val_ana
! write(78,'(2I4,ES25.17)') a, b, integral_val
  end do
  end do
  
  ! Write transformation matrix. This is the same as the eigenvectors except
  ! without the grid as the first column 
  if (only_bound) then
    open(11, file="transformation_matrix_l"//trim(int2str(para%l))//".dat",    &
    &    form="formatted", action="write")
    do i = 1, size(matrix(:,1))
      do j = 1, size(matrix(i,:))
        if (eigen_vals(j) > zero) cycle
        write(11,'(ES25.17, 1x)', advance = 'No')                              &
        & matrix(i,j) / (sqrt(grid%weights(i)) * grid%r(i))
      end do
      write(11,*) ' '
    end do
    close(11)
  else
    open(11, file="transformation_matrix_l"//trim(int2str(para%l))//".dat",    &
    &    form="formatted", action="write")
    do i = 1, size(matrix(:,1))
      write(11,*)                                                              &
      & (matrix(i,j) / (sqrt(grid%weights(i)) * grid%r(i)),                    &
      & j = 1, size(matrix(i,:)))
    end do
    close(11)
  end if
  
end program dvr_diag 
