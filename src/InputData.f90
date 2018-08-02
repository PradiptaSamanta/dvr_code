module InputData

  use constants
  use input_mod

  implicit none

  save 

  real(dp) :: r_min, r_max, mass, beta, e_max, nev_fac
  integer  :: m, nl, nr, l_max, z
  logical  :: mapped_grid, only_bound, dvr_diag, dvr_integrals, trans_integrals
  character(255) :: maptype, read_envelop, pottype, pot_filename, read_envelope


end module InputData
