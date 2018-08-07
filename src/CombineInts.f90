module CombineInts

  use constants
  use util_mod
  use DVRData, only : combined_two_e_int, para, grid

  implicit none

  contains

  subroutine CombineIntegrals() 

  real(idp), allocatable     :: rad_tp_elements(:,:)
  real(idp), allocatable     :: ang_tp_elements(:,:)
  integer                    :: i, j, a, b, c, d, l, l_ang_max, k, k_max
  integer                    :: a_rad, a_ang, b_rad, b_ang 
  integer                    :: c_rad, c_ang, d_rad, d_ang 
  integer                    :: running_index_tot
  integer                    :: running_index_rad, running_index_ang
  integer                    :: n_lm, n_n, n_all

  k_max = 2*para%l    ! Multipole order
  
  call read_ascii_table(rad_tp_elements,                                       &
  &                     'Results/twoparticle_rad_elements_l0.dat')

  call read_ascii_table(ang_tp_elements,                                       &
  &                     'Results/ang_element_l0.dat')

  n_n = nint(sqrt(real(size(rad_tp_elements(:,1)),kind=idp)))
  n_lm = nint(sqrt(sqrt(real(size(ang_tp_elements(:,1)),kind=idp))))
  n_all = n_n * n_lm 

  write(*,*) n_n, n_lm, n_all
  allocate(combined_two_e_int(n_all**4))
  combined_two_e_int = zero
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!! Four-Index Matrix Elements !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! <a b | c d> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  do k = 0, k_max
    
    write(iout,*) 'Combining two electron radial and angular integrals for l:', k
    if (allocated(rad_tp_elements)) deallocate(rad_tp_elements)
    if (allocated(ang_tp_elements)) deallocate(ang_tp_elements)

    write(*,*) 'Get here 1'
    call read_ascii_table(rad_tp_elements,                                     &
    &            'Results/twoparticle_rad_elements_l'//trim(int2str(k))//'.dat')
    call read_ascii_table(ang_tp_elements,                                     &
    &                         'Results/ang_element_l'//trim(int2str(k))//'.dat')
    write(*,*) 'Get here 2'
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
            if (a_rad .ne. c_rad) cycle ! Radial Kronecker Delta \delta_ac
            if (b_rad .ne. d_rad) cycle ! Radial Kronecker Delta \delta_bd
            running_index_ang = running_index_ang + 1
            running_index_tot = running_index_tot + 1
            write(76,*) running_index_rad, running_index_ang, running_index_tot
            flush(6)
            combined_two_e_int(running_index_tot) =                          &
            & combined_two_e_int(running_index_tot) +                        &
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
    
    write(*,*) 'Get here 3'
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
          &                              combined_two_e_int(running_index_tot)
        end do
      end do
    end do
  end do
  close(11)

  end subroutine CombineIntegrals

  subroutine CombineOrbInts()

  end subroutine CombineOrbInts
end module CombineInts
