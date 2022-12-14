module DVRData
  use constants

  implicit none

  public

  !! @description: Data structure for splining data
  !! @param: n    Number of (x,f(x)) pairs of the spline
  !! @param: cm   Array of interpolation help coefficients
  !! @param: x    Array of x axis values
  !! @param: y    Array of corresponding y axis values
  type spline_t
    integer                :: n
    real(idp), allocatable :: cm(:)
    real(idp), allocatable :: x(:)
    real(idp), allocatable :: y(:)
  end type spline_t

  !! @description: Parameters 
  !! @param: mass            Mass of particle (usually electron mass = 1 amu)
  !! @param: Z               Nuclear charge
  !! @param: r_min           Minimum $r$
  !! @param: r_max           Maximum $r$
  !! @param: r_max1          Maximum $r$ for B1 region (for inner_outer mapping)
  !! @param: r_max2          Maximum $r$ for B2 region (for inner_outer mapping)
  !! @param: pottype         Type of potential {'analytical'|'file'}
  !! @param: pot_filename    Name of file with potential data
  !! @param: mapped_grid     Decides whether mapping is employed for the grid 
  !! @param: split_grid      Logical to divide the total grid into upto two sub-region
  !! @param: maptype         Type of mapping {'diff'|'int'}
  !! @param: diagtype        Region to diagonalise {'only_inner'|'only_outer'}
  !!                         (only if maptype = 'inner_outer')
  !! @param: read_envelope   If set to something other than '', read the mapping
  !!                         envelope potential from the given filename.
  !! @param: beta            $\beta$ mapping factor
  !! @param: E_max           $E_{\max}$ mapping factor
  !! @param: dr_max          Maximal size of `dr` when using mapping.
  !!                         [huge(zero)]
  !! @param: nr              Initial number of grid points
  !! @param: ng              Actual number of grid points
  !! @param: nl              When using cardinal grid, $nl + 1$ is the number of
  !!                         grid points in each element
  !! @param: m               When using cardinal grid, $m$ is the number of
  !!                         elements
  !! @param: m1              When using cardinal grid, $m$ is the number of
  !!                         elements in the B1 region (for inner_outer mapping)
  !! @param: m2              When using cardinal grid, $m$ is the number of
  !!                         elements in the B2 region (for inner_outer mapping)
  !! @param: moveable        Whether or not the grid should move
  !! @param: coord_type      Label for coordinate system {'cartesian'|'gll'|
  !!                         'spherical'}
  !! @param: spher_method    Name of the method to use for spherical
  !!                         coordinates {'dvr'|'fbr'|'fbrgll'}
  type para_t
    real(dp)                   :: mass
    real(dp)                   :: Z 
    real(dp)                   :: r_min
    real(dp)                   :: r_max
    real(dp)                   :: r_max1
    real(dp)                   :: r_max2
    character(len=pottype_l)   :: pottype
    character(len=file_l)      :: pot_filename 
    logical                    :: mapped_grid
    logical                    :: split_grid
    logical                    :: limit_ml
    character(len=maptype_l)   :: maptype
    character(len=diagtype_l)  :: diagtype
    character(len=file_l)      :: read_envelope
    real(dp)                   :: beta
    real(dp)                   :: E_max
    real(dp)                   :: dr_max
    integer                    :: nr
    integer                    :: ng
    integer                    :: nl
    integer                    :: ml_max
    integer                    :: m
    integer                    :: m1
    integer                    :: m2
    integer                    :: l 
    integer                    :: dim_l 
    integer                    :: nev
  end type para_t
  
  !! @description: Spatial grid in a specific dimension
  !! @param: r            Array of $r$-values on the grid
  !! @param: k            Array of $k$-values on the grid
  !! @param: J            Jacobian to translate between physical grid and
  !!                      working grid, when mapping is used
  !! @param: weights      Arrays of Gauss - Legendre weights for angular
  !!                      coordinates
  !! @param: dvr_matrix   Unitary transformation matrix between the DVR and FBR
  !!                      representations
  !! @param: dr           Grid spacing (after mapping)
  !! @param: dk           $k$-spacing
  !! @param: k_max        Maximum value of $k$
  !! @param: nl           When using cardinal grid, $nl + 1$ is the number of
  !!                      grid points in each element
  !! @param: m            When using cardinal grid, $m$ is the number of
  !!                      elements
  !! @param: jac          Jacobian for the cardinal grid
  !! @param: gllp         Legendre-Gauss-Lobatto points in $[-1,1]$
  !! @param: gllw         Gaussian weights in $[-1,1]$ for Gauss-Lobatto grid
  !! @param: D_primitive  For cardinal grid, Kinetic matrix representation in
  !!                      $[-1,1]$
  type grid_t
    real(dp), allocatable   :: r(:)
    real(dp), allocatable   :: J(:)
    real(dp), allocatable   :: weights(:)
    real(dp), allocatable   :: dvr_matrix(:,:)
    real(dp)                :: dr
    integer                 :: nl
    integer                 :: m
    real(dp), allocatable   :: jac(:)
    real(dp), allocatable   :: gllp(:)
    real(dp), allocatable   :: gllw(:)
    real(dp), allocatable   :: D_primitive(:,:)
  end type grid_t

  type sph_harm_t
      integer           :: n_l  ! number of L quantum number
      integer           :: n_mp ! orders of multipole included
      ! The following are some new parameters which need to be calculated and 
      ! stored in case there is a limit on the value of the m_l qunatum number
      integer, allocatable :: n_m(:) ! Array containing the number of m_l quantum for each l
      integer, allocatable :: m_init(:) ! Array containing the number of m_l quantum for each l

  end type sph_harm_t

  type(para_t)              :: para
  type(grid_t)              :: grid
  real(dp), allocatable     :: Tkin_cardinal(:,:)
  real(dp), target, allocatable     :: pot(:,:) ! Array to store the potential for all values of l
  real(dp)               :: full_r_max

  type(sph_harm_t)       :: sph_harm

  ! The eigen values and eigen vectors obtained as a solution of the radial Schroedinger equation are defined
  real(dp), target, allocatable     :: eigen_vals(:,:), eigen_vecs(:,:,:)

  ! The angular part of two electron integrals are defined in the spherical harmonics basis
  complex(idp), allocatable  :: integrals_ang(:,:,:,:,:)

  ! The radial part of the one and two electron integrals are defined for the primitive DVR basis
  real(dp), target, allocatable :: one_e_rad_int(:,:,:)
  real(dp), target, allocatable :: two_e_rad_int(:,:,:)

  real(dp), target, allocatable :: two_e_rad_int_num(:,:,:,:,:)
  real(dp), target, allocatable :: two_e_rad_int_ana_dr(:,:,:,:,:,:,:)
  real(dp), target, allocatable :: two_e_rad_int_ana_xc(:,:,:,:,:,:,:)

  real(dp), allocatable  :: combined_two_e_int(:)
  integer                :: debug
  integer                :: nFields
  logical                :: direct_2e
  logical                :: with_field
  logical                :: prim_integrals
  character(len=32), allocatable     :: FieldComp(:)


end module DVRData
