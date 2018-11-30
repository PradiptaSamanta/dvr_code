module InputData

  use constants
  use input_mod

  implicit none

  save 

  real(dp) :: r_min, r_max, r_max1, r_max2, mass, beta, e_max, nev_fac
  integer  :: m, m1, m2, nl, nr, l_max, z, n_max, two_e_int, nfrz
  logical  :: mapped_grid, only_bound, dvr_diag, dvr_integrals, trans_integrals, orbital_ints
  character(255) :: maptype, read_envelop, pottype, pot_filename, read_envelope, diagtype


end module InputData
