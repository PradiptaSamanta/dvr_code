program angular_elements 

  use angular_mod

  implicit none

  type(sph_harm_t)           :: sph_harm
  real(idp), allocatable     :: integrals_ang(:,:,:,:,:)
  complex(idp), allocatable  :: integrals_realbasis(:,:,:,:,:)
  logical                    :: all_int

  sph_harm%n_l       = 2
  sph_harm%n_mp      = 3

  all_int = .true.
!  all_int = .false.

  call allocate_int_ang(integrals_ang, sph_harm)
  call allocate_int_ang_cmplx(integrals_realbasis, sph_harm)

! call calc_int_angular(integrals_ang, sph_harm)
  call calc_int_angular_main(integrals_ang, sph_harm)
  call calc_int_angular_compact_realbasis(integrals_realbasis, integrals_ang,  &
  &                                       sph_harm) 

  call write_int_angular(integrals_ang, sph_harm, all_int)
  call write_int_angular_realbasis(integrals_realbasis, sph_harm, all_int)

! write(*,*) wigner3j(1.0d0, 0.0d0, 1.0d0, -1.0d0, 0.0d0, 0.0d0 )
end program angular_elements
