module angular_mod 

  use dvr_spline_mod
  
  implicit none

  public
  
  !! @description: Parameters 
  !! @param: mass            Mass of particle (usually electron mass = 1 amu)
  !! @param: r_min           Minimum $r$
  !! @param: r_max           Maximum $r$
  !! @param: pottype         Type of potential {'analytical'|'file'}
  !! @param: pot_filename    Name of file with potential data
  !! @param: mapped_grid     Decides whether mapping is employed for the grid 
  !! @param: maptype         Type of mapping {'diff'|'int'}
  !! @param: read_envelope   If set to something other than '', read the mapping
  !!                         envelope potential from the given filename.
  !! @param: beta            $\beta$ mapping factor
  !! @param: E_max           $E_{\max}$ mapping factor
  !! @param: dr_max          Maximal size of `dr` when using mapping.
  !!                         [huge(zero)]
  !! @param: nr              Numer of grid points
  !! @param: nl              When using cardinal grid, $nl + 1$ is the number of
  !!                         grid points in each element
  !! @param: m               When using cardinal grid, $m$ is the number of
  !!                         elements
  !! @param: moveable        Whether or not the grid should move
  !! @param: coord_type      Label for coordinate system {'cartesian'|'gll'|
  !!                         'spherical'}
  !! @param: spher_method    Name of the method to use for spherical
  !!                         coordinates {'dvr'|'fbr'|'fbrgll'}
  type para_t
    real(idp)                   :: mass
    real(idp)                   :: r_min
    real(idp)                   :: r_max
    character(len=pottype_l)    :: pottype
    character(len=file_l)       :: pot_filename 
    logical                     :: mapped_grid
    character(len=maptype_l)    :: maptype
    character(len=file_l)       :: read_envelope
    real(idp)                   :: beta
    real(idp)                   :: E_max
    real(idp)                   :: dr_max
    integer                     :: nr
    integer                     :: nl
    integer                     :: m
    integer                     :: l 
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
    real(idp), allocatable   :: r(:)
    real(idp), allocatable   :: J(:)
    real(idp), allocatable   :: weights(:)
    real(idp), allocatable   :: dvr_matrix(:,:)
    real(idp)                :: dr
    integer                  :: nl
    integer                  :: m
    real(idp), allocatable   :: jac(:)
    real(idp), allocatable   :: gllp(:)
    real(idp), allocatable   :: gllw(:)
    real(idp), allocatable   :: D_primitive(:,:)
  end type grid_t

contains
  
  !! @description: Evaluate Wigner 3-j symbol
  !! @param: J1 J1 angular momentum
  !! @param: J2 J2 angular momentum
  !! @param: J3 J3 angular momentum
  !! @param: M1 Projection of J1 angular momentum
  !! @param: M2 Projection of J2 angular momentum
  !! @param: M3 Projection of J3 angular momentum
  real(idp) function wigner3j(J1,J2,J3,M1,M2,M3)

    real(idp), intent(in) :: J1
    real(idp), intent(in) :: J2
    real(idp), intent(in) :: J3
    real(idp), intent(in) :: M1
    real(idp), intent(in) :: M2
    real(idp), intent(in) :: M3

    integer   :: i, k, intexp
    real(idp) :: C, sumk, term
    real(idp), dimension(0:99) :: fact

    ! Compute table of factorials
    fact(0) = one
    do i = 1, 99
      fact(i) = i * fact(i-1)
    end do

    ! Check for invalid input
    if (isfrac(J1+J2+J3) .or. isfrac(J1+M1) .or. isfrac(J2+M2) .or. &
    &   isfrac(J3-M3) .or. isfrac(-J1+J3-M2) .or. isfrac(-J2+J3+M1)) then
      write(*,'(A)') "=================== ERROR ========================"
      write(*,'(A)') "Message: Invalid input of Wigner 3j symbol"
      write(*,'(A)') "=================================================="
      stop
    end if

    ! Compute Clebsch-Gordan coefficient C
    if ( (J3 < abs(J1-J2)) .or.  &
    &  (J3 > (J1+J2)) .or.       &
    &  (abs(M1) > J1) .or.       &
    &  (abs(M2) > J2) .or.       &
    &  (abs(M3) > J3)) then
      C = zero
    else
      C = sqrt((J3+J3+one)/fact(nint(J1+J2+J3+one)))
      C = C * sqrt(fact(nint(J1+J2-J3))*fact(nint(J2+J3-J1))* &
      &   fact(nint(J3+J1-J2)))
      C = C * sqrt(fact(nint(J1+M1))*fact(nint(J1-M1))*fact(nint(J2+M2))* &
      &   fact(nint(J2-M2))*fact(nint(J3-M3))*fact(nint(J3+M3)))
      sumk = zero
      do k = 0, 99
        if (J1+J2-J3-K < zero) cycle
        if (J3-J1-M2+K < zero) cycle
        if (J3-J2+M1+K < zero) cycle
        if (J1-M1-K    < zero) cycle
        if (J2+M2-K    < zero) cycle
        term = fact(nint(J1+J2-J3-k))*fact(nint(J3-J1-M2+k))* &
        &      fact(nint(J3-J2+M1+k))*fact(nint(J1-M1-k))*    &
        &      fact(nint(J2+M2-k))*fact(k)
        if (mod(k,2) == 1) term = -term
        sumk = sumk + one/term
      end do
      C = C * sumk
    end if

    ! calculate 3j symbol from Clebsch-Gordan coefficient
    ! Note: Nagfor treats expressions like (-1)^n with real n as illegal
    !       (because it is) and will throw a floating invalid operation error.
    !       So in order to evaluate the wigner3j symbol we first have to
    !       convert the exponent to an integer expression.
    intexp = nint(J1-J2-M3)
    wigner3j = (-one) ** intexp / sqrt(two * J3 + one) * C

  end function wigner3j
  
  !! @description: Check if argument is fractional
  !! @param: x Argument
  logical function isfrac(x)

    real(idp), intent(in) :: x

    real(idp), parameter  :: eps = 1.0d-8

      if ((abs(x)-int(abs(x))) > eps) then
         isfrac = .true.
      else
         isfrac = .false.
      end if

  end function isfrac

  
  !! @description: Report an allocation error. Print out given the given
  !!               `module_name`, the given `routine_name`, a description of the
  !!               error that occured, and optionally an additional message
  !! @param: i              Error code (`stat` in allocate)
  !! @param: error_message  Additional error message to print
  subroutine allocerror(i, error_message)

    integer,                    intent(in) :: i
    character(len=*), optional, intent(in) :: error_message

    integer :: proc_id, error
    logical :: mpi_is_initialized

    if (i > 0) then
      write(*,'("")')
      write(*,'("================= ALLOCERROR =====================")')
      write(*,'("ERROR: Could not allocate memory")')
      stop
    end if

  end subroutine allocerror

end module angular_mod 
