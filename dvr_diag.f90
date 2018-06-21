program dvr_diag 

  use dvr_diag_mod

  implicit none

  real(idp), allocatable      :: pot(:)
  !type(state_t), target      :: psi, psi_temp
  !type(para_t)               :: para
  !real(idp), allocatable     :: eigen_vals(:)
  !real(idp), allocatable     :: matrix(:,:)
  !real(idp), allocatable     :: Haux(:,:), T_cardinal(:,:)
  !type(grid_t)               :: grid
  !type(dyn_generator_t)      :: gen
  !type(prop_work_t)          :: work
  !type(pulse_t), allocatable :: pulses(:)
  !integer                    :: i, j
  !character(len=runfolder_l) :: runfolder

  write(*,*) 'Bonjour, tout le monde.'
  call dummy_dummy()

  !call get_command_argument(1, runfolder)
  !call read_para(para, runfolder, 'config')
  !call init(para, grid=grid, gen=gen, pulses=pulses, quiet=.false.)
  !write(*,*) " * Spectral radius: ", gen%ham%de
  !write(*,*) " * Minimum energy: ", gen%ham%emin
  !write(*,'("")')


  !! Set up potential (for now: 1/r for hydrogen) 
  !do i = 1,size(gen%ham%ops(1)%a)
  !   gen%ham%ops(1)%a(i) = - one / grid%dim(1)%r(i)
  !end do

  !! Write out potential as it was interpolated on the Gauss-Lobatto grid
  !call write_op(gen%ham, pulses, pulse_val_i=0, grid=grid, op_type='pot',      &
  !&            op_surf=1, filename=join_path(runfolder, "pot_on_grid.dat"))

  !! Get banded storage format of Hamiltonian matrix in the FEM-DVR basis
  !call get_real_surf_matrix_cardinal(matrix, grid, gen%ham, pulses, isurf=1,   &
  !&                                  pulse_val_i=0)

  !! Diagonalize Hamiltonian matrix which is stored in banded format.
  !! nev specifies the first nev eigenvalues and eigenvectors to be extracted.
  !! If needed, just add more
  !call diag_arpack_real_sym_matrix(matrix, formt='banded', n=size(matrix(1,:)),&
  !&                                nev=90, which='SA', eigenvals=eigen_vals,   &
  !&                                rvec=.true.)

  !! Write eigenvalues.
  !open(11, file=join_path(runfolder, "eigenvalues_GLL.dat"), form="formatted", &
  !&    action="write")
  !write(11,*) '#_______________________________________________________________'
  !write(11,*) "#              eigenvalues for hydrogen with l = 0 "
  !write(11,*) '#_______________________________________________________________'
  !write(11,*) "#     index    -    eigenvalue    -    eigenvector normalization"
  !write(11,*) ""
  !do i = 1, size(eigen_vals(:))
  !  write(11,*) i-1, eigen_vals(i), dot_product(matrix(:,i), matrix(:,i))
  !end do
  !close(11)

  !! Write eigenvectors. Here, if we want to represent the eigenvectors in the
  !! physical grid instead of numerical grid defined by the normalized FEM-BASIS,
  !! we should divide by the square root of the Gaussian interpolating weights,
  !! see eg. arXiv:1611.09034 (2016).
  !open(11, file=join_path(runfolder, "eigenvectors_GLL.dat"), form="formatted",&
  !&    action="write", recl=100000)
  !do i = 1, size(matrix(:,1))
  !  write(11,*) grid%dim(1)%r(i),                                              &
  !  & (matrix(i,j) / sqrt(grid%dim(1)%weights(i)), j = 1, size(matrix(i,:)))
  !end do
  !close(11)


  !! Important, for propagation purposes, ALL vectors that need to be propagated
  !! must be written in the normalized FEM-DVR basis. This is, if for any
  !! reason an arbitrary vector (not obtained from diagonalization in the
  !! normalized FEM-DVR basis) needs to be propagated, it should be interpolated
  !! at the Gauss-Lobatto points AND multiplied by the square root of the
  !! Gaussian interpolating weights associated with each collocation point.
  !call alloc_state(psi, grid, 1)
  !call mold_state(psi_temp, psi)
  !forall (i = 1:size(psi%psi,1)) psi%psi(i,1,1,1,1) = matrix(i,10)
  !                                                !/sqrt(grid%dim(1)%weights(i))
  !! State with index 10 = highly excited but still bound
  !open(11, file=join_path(runfolder, "psi_initial_time.dat"), form="formatted",&
  !&     action="write", recl=100000)
  !do i = 1, size(matrix(:,1))
  !  psi_temp%psi(i,1,1,1,1) = psi%psi(i,1,1,1,1) / sqrt(grid%dim(1)%weights(i))
  !  write(11,*) grid%dim(1)%r(i),                                              &
  !    &         real(psi%psi(i,1,1,1,1))  / sqrt(grid%dim(1)%weights(i)),      &
  !    &         aimag(psi%psi(i,1,1,1,1)) / sqrt(grid%dim(1)%weights(i))
  !end do
  !close(11)

end program dvr_diag 
