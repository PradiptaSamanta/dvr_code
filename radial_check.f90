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
  integer                    :: i, j, a, b, c, d, l, l_val
  real(idp),  allocatable    :: file_r(:), file_pot(:)
  integer, allocatable       :: ipiv(:)
  real(idp)                  :: curr_r, curr_r_prime, integral_val, dr
  real(idp), allocatable     :: work(:)
  integer                    :: info

  l_val=0

  para%pottype       = 'analytical' ! 'density_file', 'file', 'analytical'
  para%pot_filename  = 'input_pot.in' 
  para%r_min         = 0.0
  para%r_max         = 10.0
  para%m             = 200
  para%nl            = 5
  para%nr            = 1001  !nr = m * nl + 1
  para%l             = l_val !Multipole order
  para%Z             = 1
  para%mass          = 1.0

  para%mapped_grid   = .false.
  para%maptype       = 'diff'
  para%read_envelope = ''
  para%beta          = 0.015
  para%E_max         = 1d-5

  call init_grid_dim_GLL(grid, para)
 
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
  
  do j = 1, 15
    open(11,file="explicit_dvr_primitive"//trim(int2str(j))//".dat",           &
    &    form="formatted", action="write", recl=100000)
    open(12,file="full_explicit_dvr_primitive"//trim(int2str(j))//".dat",      &
    &    form="formatted", action="write", recl=100000)
    do i = 1, size(grid%r)
      write(11,'(4ES25.17)') grid%r(i), sqrt(grid%weights(i)), sqrt(grid%gllw(mod(i,grid%nl))), &
      &        dvr_primitive_val_esteban(j, para, grid, grid%r(i))
    end do
    do i = 1, 10001
      curr_r = (i-1) * 0.00001d0
      write(12,'(2ES25.17)') curr_r,                                           &
      &        dvr_primitive_val_esteban(j, para, grid, curr_r)!/sqrt(grid%gllw(mod(j,para%nl)))
    end do
    close(11)
    close(12)
  end do

  dr = 0.0001d0  
  do a = 1, 10
  do b = 1, 10
    do c = 1, 10
    do d = 1, 10
      open(13,file="numerical_integrals_l"//trim(int2str(l_val))//".dat",      &
      &    form="formatted", action="write", recl=100000)
      integral_val = zero
      do i = 1, 2001
        curr_r = (i-1) * dr
        if (abs(dvr_primitive_val_esteban(a, para, grid, curr_r)) < 1d-16) cycle
        if (abs(dvr_primitive_val_esteban(c, para, grid, curr_r)) < 1d-16) cycle
        do j = 1, 2001
          curr_r_prime = (j-1) * dr
          if (curr_r + curr_r_prime < 1d-16) cycle
          integral_val = integral_val +                                        &
          &              dvr_primitive_val_esteban(a, para, grid, curr_r) *            &
          &              dvr_primitive_val_esteban(c, para, grid, curr_r) *            &
          &              dvr_primitive_val_esteban(b, para, grid, curr_r_prime) *      &
          &              dvr_primitive_val_esteban(d, para, grid, curr_r_prime) *      &
          &              ( min(curr_r, curr_r_prime)**(l_val) /                &
          &                max(curr_r, curr_r_prime)**(l_val+1) ) * dr * dr
          !write(*,*) integral_val, dvr_primitive_val(a, para, grid, curr_r),   &
          !& dvr_primitive_val(b, para, grid, curr_r_prime), min(curr_r, curr_r_prime),&
          !& max(curr_r, curr_r_prime)
        end do
      end do
      write(13,'(4I4, ES25.17)') a, b, c, d, integral_val
      write(*,*) 'Writing integral', a, b, c, d
    end do
    end do
  end do
  end do
  close(13)
  
  do a = 1, 10
  do b = 1, 10
    do c = 1, 10
    do d = 1, 10 
      open(14,file="GLL_numerical_integrals_l"//trim(int2str(l_val))//".dat",  &
      &    form="formatted", action="write", recl=100000)
      integral_val = zero
      do i = 1, para%nr 
        curr_r = grid%r(i)
        do j = 1, para%nr 
          curr_r_prime = grid%r(j) 
          if (curr_r + curr_r_prime < 1d-16) cycle
          integral_val = integral_val +                                        &
          &              dvr_primitive_val_esteban(a, para, grid, curr_r) *            &
          &              dvr_primitive_val_esteban(c, para, grid, curr_r) *            &
          &              dvr_primitive_val_esteban(b, para, grid, curr_r_prime) *      &
          &              dvr_primitive_val_esteban(d, para, grid, curr_r_prime) *      &
          &              ( min(curr_r, curr_r_prime)**(l_val) /                &
          &                max(curr_r, curr_r_prime)**(l_val+1) ) *            &
          &              grid%weights(i) * grid%weights(j)
          !write(*,*) integral_val, dvr_primitive_val(a, para, grid, curr_r),   &
          !& dvr_primitive_val(b, para, grid, curr_r_prime), min(curr_r, curr_r_prime),&
          !& max(curr_r, curr_r_prime)
        end do
      end do
      if (abs(integral_val) < 1d-16) integral_val = zero !truncate num. zeroes
      write(14,'(4I4, ES25.17)') a, b, c, d, integral_val
      write(*,*) 'Writing GLL integral', a, b, c, d
    end do
    end do
  end do
  end do
  close(14)

end program radial_check
