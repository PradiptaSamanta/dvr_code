program angular_elements 

  use angular_mod

  implicit none

  type(sph_harm_t)        :: sph_harm
  real(idp), allocatable  :: integrals_ang(:,:,:,:,:)

  sph_harm%n_l       = 4
  sph_harm%n_mp      = 7

  call allocate_int_ang(integrals_ang, sph_harm)

  call calc_int_angular(integrals_ang, sph_harm)

  call write_int_angular(integrals_ang, sph_harm)

! write(*,*) wigner3j(1.0d0, 0.0d0, 1.0d0, -1.0d0, 0.0d0, 0.0d0 )
end program angular_elements
