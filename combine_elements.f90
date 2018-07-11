program comb_elements

  use combine_mod

  implicit none

  real(idp), allocatable     :: pot(:)
  type(para_t)               :: para
  type(grid_t)               :: grid
  real(idp), allocatable     :: rad_sp_elements(:,:)
  real(idp), allocatable     :: rad_tp_elements(:,:)
  real(idp), allocatable     :: ang_sp_elements(:,:)
  real(idp), allocatable     :: ang_tp_elements(:,:)
  real(idp), allocatable     :: combined_sp_elements(:)
  real(idp), allocatable     :: combined_tp_elements(:)
  integer                    :: i, j, a, b, c, d, l, l_ang_max, k, k_max
  integer                    :: a_rad, a_ang, b_rad, b_ang 
  integer                    :: c_rad, c_ang, d_rad, d_ang 
  integer                    :: running_index_tot
  integer                    :: running_index_rad, running_index_ang
  logical                    :: inversion_check, alternative_formula
  real(idp),  allocatable    :: file_r(:), file_pot(:)
  integer, allocatable       :: ipiv(:)
  real(idp)                  :: full_r_max
  real(idp), allocatable     :: work(:)
  integer                    :: info
  integer                    :: n_lm, n_n, n_all

  k_max = 6    ! Multipole order
  
  ! Initialise with zeroth-order multipole

  call read_ascii_table(rad_tp_elements,                                       &
  &                     'Results/twoparticle_rad_elements_l0.dat')

  call read_ascii_table(ang_tp_elements,                                       &
  &                     'Results/ang_element_l0.dat')

  n_n = nint(sqrt(real(size(rad_tp_elements(:,1)),kind=idp)))
  n_lm = nint(sqrt(sqrt(real(size(ang_tp_elements(:,1)),kind=idp))))
  n_all = n_n * n_lm 

  allocate(combined_sp_elements(n_all**2))
  allocate(combined_tp_elements(n_all**4))
  combined_sp_elements = zero
  combined_tp_elements = zero
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! Two-Index Matrix Elements !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! <a | b> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do k = 0, k_max
    if (allocated(rad_sp_elements)) deallocate(rad_sp_elements)
    call read_ascii_table(rad_sp_elements,                                     &
    &         'Results/singleparticle_rad_elements_l'//trim(int2str(k))//'.dat')
    !write(*,*) size(rad_sp_elements(1,:)), size(rad_sp_elements(:,1))
    running_index_tot = 0
    running_index_rad = 0
    do a_rad = 1, n_n
      do b_rad = 1, n_n
        running_index_rad = running_index_rad + 1
        do a_ang = 1, n_lm
          do b_ang = 1, n_lm
            running_index_tot = running_index_tot + 1
            !write(*,*) a_rad, b_rad
            !write(*,*) a_ang, b_ang
            !write(*,*) k, running_index_tot
            if (a_ang .ne. b_ang) cycle !Kronecker Delta in angular part
            ! Speherical harmonics are normalized, so multiplication with one
            ! due to angular part (skipped due to redundancy)
            combined_sp_elements(running_index_tot) =                          &
            & combined_sp_elements(running_index_tot) +                        &
            & rad_sp_elements(running_index_rad,3)
          end do
        end do
      end do
    end do
  end do
  
  open(11, file="singleparticle_combined_elements.dat",                        &
  &    form="formatted", action="write")
  write(11,*) '# Primitive Combined Matrix Elements for the Two-index '//      &
  &           'Integrals for Multipole Order up to k = '//trim(int2str(k_max))
  running_index_tot = 0
  do a = 1, n_all
    do b = 1, n_all
      running_index_tot = running_index_tot + 1
      !write(*,*) a, b, running_index_tot
      write(11, '(2I8,ES25.17)') a, b, combined_sp_elements(running_index_tot)
    end do
  end do
  close(11)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!! Four-Index Matrix Elements !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! <a b | c d> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  do k = 0, k_max
    if (allocated(rad_tp_elements)) deallocate(rad_tp_elements)
    if (allocated(ang_tp_elements)) deallocate(ang_tp_elements)
    call read_ascii_table(rad_tp_elements,                                     &
    &            'Results/twoparticle_rad_elements_l'//trim(int2str(k))//'.dat')
    call read_ascii_table(ang_tp_elements,                                     &
    &                         'Results/ang_element_l'//trim(int2str(k))//'.dat')
    running_index_tot = 0
    running_index_rad = 0
    do a_rad = 1, n_n
    do b_rad = 1, n_n
      ! Already at this stage due to Kronecker Delta in radial part
      ! \delta_ac \delta_bd 
      running_index_rad = running_index_rad + 1
      do c_rad = 1, n_n
      do d_rad = 1, n_n
      running_index_ang = 0
        do a_ang = 1, n_lm
        do b_ang = 1, n_lm
          do c_ang = 1, n_lm
          do d_ang = 1, n_lm
            running_index_ang = running_index_ang + 1
            running_index_tot = running_index_tot + 1
            if (a_rad .ne. c_rad) cycle ! Radial Kronecker Delta \delta_ac
            if (b_rad .ne. d_rad) cycle ! Radial Kronecker Delta \delta_bd
            combined_tp_elements(running_index_tot) =                          &
            & combined_tp_elements(running_index_tot) +                        &
            & rad_tp_elements(running_index_rad,3) *                           &
            & ang_tp_elements(running_index_ang,5)
          end do
          end do
        end do
        end do
      end do
      end do
    end do
    end do
  end do

  open(11, file="twoparticle_combined_elements.dat",                           &
  &    form="formatted", action="write")
  write(11,*) '# Primitive Combined Matrix Elements for the Four-index '//     &
  &           'Integrals for Multipole Order up to k = '//trim(int2str(k_max))
  running_index_tot = 0
  do a = 1, n_all
    do b = 1, n_all
      do c = 1, n_all
        do d = 1, n_all
          running_index_tot = running_index_tot + 1
          write(11, '(4I8,ES25.17)') a, b, c, d,                               &
          &                              combined_tp_elements(running_index_tot)
        end do
      end do
    end do
  end do
  close(11)
  
end program comb_elements
