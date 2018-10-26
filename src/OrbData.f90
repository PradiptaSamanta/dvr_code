module OrbData

  use constants

  implicit none

  !! This is a type to define the parameters related to the orbitals 
  type orb_t

    integer   :: n_max
    integer   :: nSpatialOrbs

  end type orb_t


  type(orb_t)            :: orb

  ! Define matrices to store one and two electron 

  integer                :: nSpatialOrbs, nFrozen

  integer, allocatable   :: SpatialOrbInd(:,:,:)
  real(dp), allocatable  :: OneEInts(:,:) 
  real(dp), allocatable  :: TwoERadOrbInts_old(:,:,:)
  real(dp), allocatable  :: TwoERadOrbInts(:,:,:,:,:,:,:,:,:)
  real(dp), allocatable  :: TwoERadOrbInts_dr(:,:,:,:,:,:,:)
  real(dp), allocatable  :: TwoERadOrbInts_xc(:,:,:,:,:,:,:)
  real(dp), allocatable  :: TwoEInts(:,:,:,:)

  character(32)          :: file_int

end module OrbData
