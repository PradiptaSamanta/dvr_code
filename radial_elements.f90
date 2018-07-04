program radial_elements

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
  logical                    :: inversion_check
  real(idp),  allocatable    :: file_r(:), file_pot(:)
  integer, allocatable       :: ipiv(:)
  real(idp)                  :: full_r_max
  real(idp), allocatable     :: work(:)
  integer                    :: info

  l_val=100

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

  inversion_check    = .true.
        
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
  
  ! Set the potential to zero to only treat the kinetic part in the following
  pot = zero
  
  ! Get banded storage format of Hamiltonian matrix in the FEM-DVR basis
  call get_real_surf_matrix_cardinal(matrix, grid, pot, Tkin_cardinal)

  ! Remove 1/(2m) factor from kinetic operator to obtain properly
  ! scaled radial matrix elements
  matrix = (two * para%mass) * matrix

  !! Convert banded matrix to full matrix
  !! Watch for the para%nr-2, because the end points are not included anymore
  call mat_banded_to_full(matrix_full, matrix, para%nr-2, 0, para%nl)

  if (inversion_check) then
    allocate(matrix_all(size(matrix_full(:,1)),size(matrix_full(:,1))))
    do i = 1, size(matrix_all(:,1))
      do j = 1, size(matrix_all(:,1))
        if (i .le. j) then
          matrix_all(i,j) = matrix_full(i,j)
        else
          matrix_all(i,j) = matrix_full(j,i)
        end if
      end do
    end do
  end if
  !! Debug statements to print out initial matrix before inversion
  !open(11, file="full_matrix.dat", form="formatted",&
  !&    action="write", recl=100000)
  !do i = 1, size(matrix_full(:,1))
  !  write(11,*)                                                                &
  !  & (matrix_full(i,j),                                                       &
  !  & j = 1, size(matrix_full(i,:)))
  !end do
  !close(11)
  !open(11, file="full_matrix_all.dat", form="formatted",&
  !&    action="write", recl=100000)
  !do i = 1, size(matrix_all(:,1))
  !  write(11,*)                                                                &
  !  & (matrix_all(i,j),                                                       &
  !  & j = 1, size(matrix_all(i,:)))
  !end do
  !close(11)

  ! Invert radial kinetic matrix
  allocate(ipiv(size(matrix_full(1,:))))
  allocate(work(size(matrix_full(1,:))))
  call wrap_dsytrf(matrix_full, size(matrix_full(1,:)), ipiv, work)
  call wrap_dsytri(matrix_full, size(matrix_full(1,:)), ipiv, work)

  ! Now compute radial integral <ab|r_{<}^{l}/r_{>}^{l+1}|cd>, note that
  ! it is zero unless a = b and c = d so we only store the nonzero entries
  ! The order of storage is such that c is the inner loop, i.e.
  ! a = 1, c = 1
  ! a = 1, c = 2
  ! a = 1, c = 3
  ! ...
  ! a = 2, c = 1
  ! a = 2, c = 2
  ! ...
  ! and so on.
  ! TODO
  ! Consider due to symmetry of the matrix to also only store entries for a>c
  open(11, file="rad_elements_l"//trim(int2str(para%l))//".dat",               &
  &    form="formatted", action="write")
  write(11,*) '# Primitive Radial Matrix Elements for the Four-index '//       &
  &           'Integrals for l = '//trim(int2str(para%l))
  do a = 1, size(matrix_full(1,:))
    do b = 1, size(matrix_full(1,:))
      l = para%l
      write(11, '(2I8,ES25.17)') a, b,                                         &
      & ((real(2*l+1, idp) / (grid%r(a) * sqrt(grid%weights(a)) *              &
      &     grid%r(b) * sqrt(grid%weights(b)))) * matrix_full(a,b))            &
      & + ((grid%r(a)**l * grid%r(b)**l) / full_r_max**(2*l+1))
    end do
  end do
  close(11)
 
  if (inversion_check) then
    allocate(matrix_inv_all(size(matrix_full(:,1)),size(matrix_full(:,1))))
    do i = 1, size(matrix_inv_all(:,1))
      do j = 1, size(matrix_inv_all(:,1))
        if (i .le. j) then
          matrix_inv_all(i,j) = matrix_full(i,j)
        else
          matrix_inv_all(i,j) = matrix_full(j,i)
        end if
      end do
    end do
    allocate(unity(size(matrix_full(:,1)),size(matrix_full(:,1))))
    unity = matmul(matrix_inv_all, matrix_all)
    do i = 1, size(unity(:,1))
      do j = 1, size(unity(:,1))
        if (i == j) cycle
        if (abs(unity(i,j)) > 1d-10) then
          write(*,*) "WARNING: Inversion not successful with desired precision."
        end if
      end do
    end do
  end if
  !! Debug statements to print out matrix after full inversion
  !open(11, file="full_matrix_inverse_all.dat", form="formatted",&
  !&    action="write", recl=100000)
  !do i = 1, size(matrix_inv_all(:,1))
  !  write(11,*)                                                                &
  !  & (matrix_inv_all(i,j),                                                    &
  !  & j = 1, size(matrix_inv_all(i,:)))
  !end do
  !close(11)
  !open(11, file="unity.dat", form="formatted",&
  !&    action="write", recl=100000)
  !do i = 1, size(unity(:,1))
  !  write(11,*)                                                                &
  !  & (unity(i,j),                                                             &
  !  & j = 1, size(unity(i,:)))
  !end do
  !close(11)
  
end program radial_elements
