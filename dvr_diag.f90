program dvr_diag 

  use dvr_diag_mod

  implicit none

  real(idp), allocatable     :: pot(:)
  type(para_t)               :: para
  type(grid_t)               :: grid
  real(idp), allocatable     :: eigen_vals(:), eigen_vals_full(:)
  real(idp), allocatable     :: matrix(:,:), matrix_full(:,:)
  real(idp), allocatable     :: Tkin_cardinal(:,:)
  integer                    :: i, j, l_val, ith, jth, a, b
  real(idp)                  :: dr, curr_r, curr_r_prime, integral_val,        &
  &                             integral_val_ana, integral_val_ana2, th_1,     &
  &                             th_2, dth, rad_hyd_ef_val_a
  real(idp),  allocatable    :: file_r(:), file_pot(:), an(:)
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

  para%mapped_grid   = .true.
  para%maptype       = 'inner_outer'
  para%r_max1        = 5.0
  para%r_max2        = 30.0
  para%m1            = 50 
  para%m2            = 50
  para%read_envelope = ''
  para%beta          = 0.015
  para%E_max         = 1d-5

  para%diagtype      = 'only_outer'

  only_bound         = .false.
        
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

  !! Convert banded matrix to full matrix
  !! Watch for the para%nr-2, because the end points are not included
  call mat_banded_to_full(matrix_full, matrix, para%nr-2, 0,                   &
  &                       para%nl)
    
  open(11, file="full_matrix_l"//trim(int2str(para%l))//".dat",                &
  &    form="formatted", action="write")
  do i = 1, size(matrix_full(:,1))
    do j = 1, size(matrix_full(i,:))
      write(11,'(ES14.6, 1x)', advance = 'No')                                &
      & matrix_full(i,j)
    end do
    write(11,*) ''
  end do
  close(11)

  if (para%maptype == 'inner_outer') then
    if (para%diagtype == 'only_inner') then
      do i = (para%m1 * para%nl) + 1, para%nr-2
        do j = 1, para%nr-2
          matrix_full(i,j) = zero
        end do
      end do
      do i = 1, para%nr-2
        do j = (para%m1 * para%nl) + 1, para%nr-2
          matrix_full(i,j) = zero
        end do
      end do
    elseif (para%diagtype == 'only_outer') then
      do i = 1, (para%m1 * para%nl)
        do j = 1, para%nr-2 
          matrix_full(i,j) = zero
        end do
      end do
      do i = 1, (para%m1 * para%nl)
        do j = 1, para%nr-2 
          matrix_full(i,j) = zero
        end do
      end do
    end if
  end if
  
  open(11, file="trunc_matrix_l"//trim(int2str(para%l))//".dat",               &
  &    form="formatted", action="write")
  do i = 1, size(matrix_full(:,1))
    do j = 1, size(matrix_full(i,:))
      write(11,'(ES14.6, 1x)', advance = 'No')                                &
      & matrix_full(i,j)
    end do
    write(11,*) ''
  end do
  close(11)

  !! Diagonalize Hamiltonian matrix which is stored in banded format.
  !! nev specifies the first nev eigenvalues and eigenvectors to be extracted.
  !! If needed, just add more
  !call diag_arpack_real_sym_matrix(matrix, formt='banded', n=size(matrix(1,:)), &
  !&                   nev=para%nev, which='SA', eigenvals=eigenval_p, &
  !&                   rvec=.true.)
  
  ! Since the arpack diagonalisation becomes wonky if nev=size(matrix(1,:))
  ! we will instead perform a full BLAS diagonalisation of the symmetric
  ! matrix, matrix_full
  call diag_matrix(matrix_full, eigen_vals_full)
  
  ! Now we transfer the eigenvectors from matrix_full to matrix to keep
  ! the rest of the code intact (since it operates on matrix instead of 
  ! matrix_full). Note that we need to allocate matrix to have less columns
  ! than matrix full if para%nev < size(grid%r)
  deallocate(matrix)
  allocate(matrix(size(matrix_full(1,:)),size(matrix_full(1,:))))
  allocate(eigen_vals(size(matrix_full(1,:))))
  do i = 1, size(matrix_full(1,:))
    do j = 1, size(matrix_full(1,:))
      matrix(i,j) = matrix_full(i,j)
      eigen_vals(j) = eigen_vals_full(j)
    end do
  end do

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

  do a = 1, 15
  do b = 1, 15
  dr = 0.1d0
  integral_val_ana = zero
  do i = 1, 1001
    !write(*,*) a, b, i
    curr_r = (i-1) * dr
    rad_hyd_ef_val_a = radial_hydrogen_ef(curr_r, l_val+a, l_val, an) 
    do j = 1, 1001
      curr_r_prime = (j-1) * dr
      if (curr_r + curr_r_prime < 1d-16) cycle
      integral_val_ana = integral_val_ana +                                    &
      &              rad_hyd_ef_val_a**2 *                                     &
      &              radial_hydrogen_ef(curr_r_prime, l_val+b, l_val, an)**2 * &
      &              ( min(curr_r, curr_r_prime)**(l_val) /                    &
      &                max(curr_r, curr_r_prime)**(l_val+1) ) * dr * dr *      &
      &              curr_r * curr_r * curr_r_prime * curr_r_prime
    end do
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
    end do
  end do
  do i = 1, para%nr-2
    curr_r = grid%r(i)
    !write(*,*) two * exp(-curr_r), matrix(i,1)/(sqrt(grid%weights(i))*grid%r(i)), &
    !&          radial_hydrogen_ef(curr_r, 1, 0)
    !write(*,*) sqrt(two)**(-3) * (two - curr_r) * exp(-curr_r / two),          &
    !&          matrix(i,2)/(sqrt(grid%weights(i))*grid%r(i)),                  &
    !&          radial_hydrogen_ef(curr_r, 2, 0)
  end do
  write(*,'(2I4,2ES25.17)') a, b, integral_val, integral_val_ana
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
