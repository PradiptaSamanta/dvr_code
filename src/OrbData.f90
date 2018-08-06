module OrbData

  use constants

  implicit none

  !! This is a type to define the parameters related to the orbitals 
  type orb_t

    integer   :: n_max

  end type orb_t


  type(orb_t)            :: orb

  ! Define matrices to store one and two electron 

  integer                :: nSpatialOrbs

  integer, allocatable   :: SpatialOrbInd(:,:,:)
  real(dp), allocatable  :: OneEInt(:,:) 
  real(dp), allocatable  :: TwoEInt(:,:,:,:)

end module OrbData
