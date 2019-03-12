module DensityData

    use constants

    implicit none

    save

    real(dp) :: r_val, Yield1e, Yield2e
    integer  :: n_orb_den(2), tot_orb
    character(len=32) :: file_1rdm, file_2rdm
    complex(idp), allocatable :: DensOrb1e(:), DensOrb2e(:,:)
    complex(idp), allocatable :: PrimDens1e(:), PrimDens2e(:)

end module DensityData
