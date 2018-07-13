program transf_ints

  use gen_ints_mod
  use transf_ints_2e_mod

  implicit none

  type(integrals_t) :: ao_integrals, mo_integrals

  integer  :: nbasis

  nbasis = 200

  call init_integrals(ao_integrals, nbasis)
  call init_integrals(mo_integrals, nbasis)
  
  write(*,*) "Coming Soon"

end program transf_ints
