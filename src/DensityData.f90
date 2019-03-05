module DensityData

    use constants

    implicit none

    save

    real(dp) :: r_val
    integer  :: n_orb_den(2), tot_orb
    character(len=32) :: file_1rdm, file_2rdm
    complex(idp), allocatable :: DensOrb1e(:), DensOrb2e(:,:)

end module DensityData
