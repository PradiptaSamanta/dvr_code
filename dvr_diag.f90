program dvr_diag 

  use dvr_diag_mod

  implicit none

  real(idp), allocatable     :: pot(:)
  type(para_t)               :: para
  real(idp)                  :: mass
  type(grid_t)               :: grid
  real(idp), allocatable     :: eigen_vals(:)
  real(idp), allocatable     :: matrix(:,:)
  real(idp), allocatable     :: Tkin_cardinal(:,:)
  integer                    :: i, j

  write(*,*) 'Bonjour, tout le monde.'
  call dummy_dummy()
        
  call init_grid_dim_GLL(grid, para, .false.) 
  call init_work_cardinalbase(Tkin_cardinal, grid, mass)
  call redefine_ops_cardinal(pot)
  call redefine_GLL_grid_1d(grid)

  !! Set up potential (for now: 1/r for hydrogen, TODO: add spline module) 
  do i = 1,size(pot)
     pot(i) = - one / grid%r(i)
  end do

  !! Write out potential as it was interpolated on the Gauss-Lobatto grid
  !call write_op(gen%ham, pulses, pulse_val_i=0, grid=grid, op_type='pot',      &
  !&            op_surf=1, filename="pot_on_grid.dat")

  !! Get banded storage format of Hamiltonian matrix in the FEM-DVR basis
  call get_real_surf_matrix_cardinal(matrix, grid, pot, Tkin_cardinal)

  !! Diagonalize Hamiltonian matrix which is stored in banded format.
  !! nev specifies the first nev eigenvalues and eigenvectors to be extracted.
  !! If needed, just add more
  call diag_arpack_real_sym_matrix(matrix, formt='banded', n=size(matrix(1,:)),&
  &                                nev=90, which='SA', eigenvals=eigen_vals,   &
  &                                rvec=.true.)

  ! Write eigenvalues.
  open(11, file="eigenvalues_GLL.dat", form="formatted", &
  &    action="write")
  write(11,*) '#_______________________________________________________________'
  write(11,*) "#              eigenvalues for hydrogen with l = 0 "
  write(11,*) '#_______________________________________________________________'
  write(11,*) "#     index    -    eigenvalue    -    eigenvector normalization"
  write(11,*) ""
  do i = 1, size(eigen_vals(:))
    write(11,*) i-1, eigen_vals(i), dot_product(matrix(:,i), matrix(:,i))
  end do
  close(11)

  ! Write eigenvectors. Here, if we want to represent the eigenvectors in the
  ! physical grid instead of numerical grid defined by the normalized FEM-BASIS,
  ! we should divide by the square root of the Gaussian interpolating weights,
  ! see eg. arXiv:1611.09034 (2016).
  open(11, file="eigenvectors_GLL.dat", form="formatted",&
  &    action="write", recl=100000)
  do i = 1, size(matrix(:,1))
    write(11,*) grid%r(i),                                                     &
    & (matrix(i,j) / sqrt(grid%weights(i)), j = 1, size(matrix(i,:)))
  end do
  close(11)

end program dvr_diag 