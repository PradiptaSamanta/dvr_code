module DVRDiag

  use DVRData
  use InputData
  use dvr_diag_mod

  implicit none

  contains

  subroutine SetDVR()

    integer  :: i

    para%pottype       = 'file' ! 'density_file', 'file', 'analytical'
    para%pot_filename  = 'input_pot.in'
    para%r_min         = r_min
    para%r_max         = r_max
    para%m             = m
    para%nl            = nl
    para%nr            = m*nl + 1
    para%l             = l_max !Rotational quantum number
    para%Z             = z
    para%mass          = mass
 
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
    write(6,*) full_r_max

    if (para%pottype == 'analytical') then
      write(*,*) 'Using analytical potential'
      allocate(pot(size(grid%r)))
      do i = 1,size(pot)
         pot(i) = analytical_potential(grid%r(i), para%mass)
      end do
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
 
  end subroutine SetDVR

  subroutine DVRDiagonalization()

    integer                    :: i, j, l_val
    logical                    :: only_bound ! Only writes out bound states

    only_bound         = .true.
        
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

  end subroutine DVRDiagonalization

end module DVRDiag
