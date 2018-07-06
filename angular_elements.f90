program angular_elements 

  use angular_mod

  implicit none

  type(sph_harm_t)        :: sph_harm
  real(idp), allocatable  :: integrals_ang(:,:,:,:,:)

  logical                 :: use_selection

  use_selection      = .true.
  sph_harm%n_l       = 4
  sph_harm%n_mp      = 6

  call allocate_int_ang(integrals_ang, sph_harm)

  call calc_int_angular(integrals_ang, sph_harm)

! if 
  call write_int_angular(integrals_ang, sph_harm, use_selection)

  write(*,*) wigner3j(1.0d0, 0.0d0, 1.0d0, -1.0d0, 0.0d0, 0.0d0 )
  write(*,*) wigner3j(1.0d0, 1.0d0, 0.0d0, -1.0d0, 0.0d0, 0.0d0 )
  write(*,*) wigner3j(1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 )
end program angular_elements
