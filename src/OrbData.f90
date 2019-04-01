module OrbData

  use constants

  implicit none

  !! This is a type to define the parameters related to the orbitals 
  type orb_t

    integer   :: n_max
    integer   :: n_inner, n_outer
    integer   :: nSpatialOrbs
    logical   :: shift_int
    logical   :: reduce_orb
    integer   :: break
    integer   :: n_red
    integer   :: n_shift_out_orb

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

  ! whether the inner orbitals are divided further
  logical :: break_inn
  ! break_orb(1) : Number of orbital before the break
  ! break_orb(2) : Number of orbital to remove
  ! break_orb(3) : Number of orbital after the break
  integer :: break_orb(3)

end module OrbData
