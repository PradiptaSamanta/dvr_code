module dvr_diag_mod
  
  implicit none

  public

  integer, parameter :: idp = kind(1.0d0) !Value of internal double precision
  integer, parameter :: datline_l    = 255 ! Maximum length of a line in a file
  integer, parameter :: file_l       = 512 ! Maximum length of filenames/paths
  integer, parameter :: maptype_l    = 4   ! length of a maptype
  
  ! named real numbers
  real(idp), parameter :: zero  = 0.0_idp, one = 1.0_idp, two = 2.0_idp
  real(idp), parameter :: eight = 8.0_idp
  real(idp), parameter :: half  = 0.5_idp, quart = 0.25_idp
  real(idp), parameter :: three = 3.0_idp, four = 4.0_idp
  real(idp), parameter :: third = one/three
  real(idp), parameter :: Pi     = 3.14159265358979323846_idp
  real(idp), parameter :: Pihalf = 1.5707963267948965600_idp
  real(idp), parameter :: FourPi = 4.0_idp * Pi

  ! named complex numbers
  complex(idp), parameter :: czero = (0.0_idp, 0.0_idp)
  complex(idp), parameter :: cone  = (1.0_idp, 0.0_idp)
  complex(idp), parameter :: ci    = (0.0_idp, 1.0_idp)
  complex(idp), parameter :: cid   = (1.0_idp, 1.0_idp)
  
  !! @description: Parameters 
  !! @param: mass            Mass of particle (usually electron mass = 1 amu)
  !! @param: r_min           Minimum $r$
  !! @param: r_max           Maximum $r$
  !! @param: maptype         Type of mapping {'diff'|'int'}
  !!                         mapping)
  !! @param: read_envelope   If set to something other than '', read the mapping
  !!                         envelope potential from the given filename.
  !! @param: write_envelope  If set to something other than '', write the
  !!                         mapping envelope potential to the given file.
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
    !character(len=maptype_l)    :: maptype
    !character(len=file_l)       :: read_envelope
    !character(len=file_l)       :: write_envelope
    !real(idp)                   :: beta
    !real(idp)                   :: E_max
    !real(idp)                   :: dr_max
    integer                     :: nr
    integer                     :: nl
    integer                     :: m
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
  !! @param: mapped       Whether or not a mapped grid is used in this spatial
  !!                      coordinate
  !! @param: maptype      If `mapped` is `.true.`, type of mapping (`diff`,
  !!                      `int`)
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
    !real(idp), allocatable   :: k(:)
    real(idp), allocatable   :: J(:)
    real(idp), allocatable   :: weights(:)
    real(idp), allocatable   :: dvr_matrix(:,:)
    real(idp)                :: dr
    !real(idp)                :: dk
    !real(idp)                :: k_max
    logical                  :: mapped
    character(len=maptype_l) :: maptype
    integer                  :: nl
    integer                  :: m
    real(idp), allocatable   :: jac(:)
    real(idp), allocatable   :: gllp(:)
    real(idp), allocatable   :: gllw(:)
    real(idp), allocatable   :: D_primitive(:,:)
  end type grid_t

contains

  !! @description: Initializes the dimension `dim` of the given grid as a GLL
  !! (Gauss-Lobatto-Legendre) grid according to the grid parameters given in
  !! `para`. If mapped is set to true the grid will be initialized as a mapped
  !! grid with the maptype specified in `para`. Otherwise the grid will be
  !! initialized as a unmapped grid.
  !! @param: grid     Grid for which the dimension dim should be initialized
  !!                  as a cartesian grid.
  !! @param: para     Parameters defining the spatial grid.
  !! @param: mapped   Defines if the grid should be initialized as a mapped or
  !!                  constant grid.
  subroutine init_grid_dim_GLL(grid, para, mapped)
    type(grid_t), intent(inout) :: grid
    type(para_t), intent(in)    :: para
    logical,      intent(in)    :: mapped

    real(idp), allocatable      :: r_env(:), V_env(:), weight(:), jac(:),  &
    &                              X_mapped(:), lobatto(:), weights(:)
    real(idp)                   :: r_min, r_max, V_dr, min_dr, beta, E_max
    integer                     :: error, nl, m, l, j, nr_env, th, nr
    character(len=datline_l)    :: header
    character(len=file_l)       :: write_envelope, read_envelope
    !character(len=maptype_l)    :: maptype
    !character (len = maptype_l) :: basis_for_mapping

    nr    = para%nr
    nl    = para%nl
    m     = para%m
    r_min = para%r_min
    r_max = para%r_max

    if (m * nl + 1 .ne. nr) then
      write(*,*) "ERROR: For the GLL grid, nr must equal m * nl + 1"
    end if

    if (mapped) then ! Initialize dimension dim as mapped GLL grid

      write(*,*) "ERROR: Mapped grids not yet implemented."
      stop

      !beta = para%beta
      !E_max = para%E_max
      !maptype = para%maptype
      !read_envelope = para%read_envelope
      !write_envelope = para%write_envelope

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Mapping Block - disabled for now !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !! sanity checks
      !if (maptype .eq. 'diff') then
      !  basis_for_mapping = 'sin'
      !elseif(maptype .eq. 'int') then
      !  basis_for_mapping = 'exp'
      !end if

      !grid%nl = nl
      !grid%m = m
      !! Get envelope potential
      !! This goes together with a constant spatial envelop grid r_env
      !if (read_envelope /= '') then
      !  call read_ascii(r_env, V_env, read_envelope)
      !else
      !  call min_dr_1d(para, 1, min_dr)
      !  nr_env = floor((r_max - r_min) / min_dr) + 1
      !  if (nr_env < 1) then
      !    write(*,*) "********* nr_env =", nr_env
      !    call abort_error(module_name, routine_name, "Internal error in &
      !    &calculating the temporary grid, for mapping. Possibly some &
      !    &potential data files could not be read correctly.")
      !  end if
      !  call get_spatial_coord(r_env, V_dr, weight, r_min, r_max, nr_env,      &
      !  &                      'const')
      !  call envelope_pot_1d(V_env, r_env, para, 1)
      !end if
      !if (write_envelope /= '') then
      !  header = "#  r[iu]                   V[iu]"
      !  call write_ascii(r_env, V_env, write_envelope, header)
      !end if
      !! Map the grid: calculate the coordinates r(:) and the jacobian J(:)
      !select case (maptype)
      !  case('diff')
      !    grid%maptype = maptype
      !    call map_diff_1d(X_mapped, grid%J, r_min, r_max, r_env,     &
      !    &                V_env, nr, beta, para%mass, E_max,        &
      !    &                basis_for_mapping)
      !    grid%dr = one

      !    ! Calculate momentum grid
      !    call get_mom_coord(grid%k, grid%dk,                &
      !    &                  grid%k_max, grid%dr, nr,        &
      !    &                  mapped=.true.)

      !  case('int')
      !    call map_int_1d(X_mapped, grid%J, r_min, r_max, r_env,      &
      !    &               V_env, nr, beta, para%mass, E_max,         &
      !    &               basis_for_mapping)
      !    deallocate(grid%J)

      !    ! Calculate momentum grid
      !    call get_mom_coord(grid%k, grid%dk,                &
      !    &                  grid%k_max, grid%dr, nr,        &
      !    &                  mapped=.true.)

      ! case default
      !  call abort_error(module_name, routine_name, "Unable to initialize &
      !  &mapped GLL grid for dimension "//trim(nr2str(1))//"! (Reason: &
      !  &Unknown maptype: `"//trim(maptype)//"`)")
      !end select

      !! Set variables indicating that the grid is of type mapped
      !grid%mapped = .true.
      !grid%maptype = maptype
      !

     
      !! Calculate the global weights and GLL points
      !allocate(jac(m), stat = error)
      !call allocerror(error)
      !call get_lobatto_points(nl, grid%gllp, grid%gllw,      &
      !&                       grid%D_primitive)
      !allocate(lobatto(0:nl),stat=error)
      !call allocerror(error)
      !allocate(weights(0:nl),stat=error)
      !call allocerror(error)

      !do j = 0, nl
      !  lobatto(j) = grid%gllp(j)
      !  weights(j) = grid%gllw(j)
      !end do
      !allocate(grid%weights(size(X_mapped)), stat = error)
      !call allocerror(error)
      !allocate(grid%jac(m),stat=error)
      !call allocerror(error)
      !do l = 1,m
      !  th = nl*(l-1) + nl +1
      !  jac(l) = half*(X_mapped(th) - X_mapped(th-nl))
      !  grid%jac(l) = half*(X_mapped(th) - X_mapped(th-nl))
      !end do

      !do j = 0, nl
      !  grid%weights(j+1) = jac(1)*weights(j)
      !end do
      !do l=2,m
      ! do j = 0,nl
      !   if (j == 0) then
      !     grid%weights(j+nl*(l-1)+1)                                 &
      !     & = grid%weights(nl+nl*(l-2)+1) + jac(l)*weights(0)
      !    else
      !     grid%weights(j+nl*(l-1)+1) = jac(l)*weights(j)
      !   end if
      !  end do
      ! end do
      !deallocate(weights,stat=error)
      !call allocerror(error)

      !if (allocated(grid%r)) deallocate(grid%r)
      !call allocerror(error)

      !allocate(grid%r(size(X_mapped)), stat=error)
      !call allocerror(error)

      !do l = 1,m
      !  th = nl*(l-1) + nl + 1
      !  do j = 0,nl
      !    grid%r(j+nl*(l-1) + 1 )                                     &
      !    & =   half*(X_mapped(th) - X_mapped(th-nl))* lobatto(j)              &
      !    &   + half*(X_mapped(th) + X_mapped(th-nl))
      !  end do
      !end do

      !! Cleanup
      !deallocate(r_env, V_env, lobatto, X_mapped)

    else ! Initialize dimension dim as unmapped GLL grid

      grid%nl = nl
      grid%m = m
      grid%dr = one
      call get_lobatto_points(nl, grid%gllp, grid%gllw,      &
      &                       grid%D_primitive)

      allocate (grid%r(nl * m + 1), stat=error)
      call allocerror(error)
      do l = 1, m, 1
        do j = 0, nl-1, 1
          grid%r(j + nl * (l - 1) + 1) =                              &
          &        grid%gllp(j) * (r_max - r_min) / (real(2 * m,idp)) &
          &        + r_min + (real(l,idp)- half) * (r_max - r_min) / real(m,idp)
        end do
      end do

      grid%r(nl * m + 1) =                                            &
      &            grid%gllp(nl) * (r_max - r_min) / (real(2 * m,idp))&
      &            + r_min + (real(m,idp)- half) * (r_max - r_min) / real(m,idp)

      allocate (grid%jac(m), stat=error)
      call allocerror(error)
      do j = 1, m, 1
        grid%jac(j) = (r_max - r_min) / (real(2 * m,idp))
      end do

      allocate(grid%weights(nl * m + 1), stat = error)
      call allocerror(error)
      do j = 0, nl, 1
        grid%weights(j+1) = grid%gllw(j)
      end do

      do l = 2, m, 1
        do j = 0, nl, 1
          if (j == 0) then
            grid%weights(j + nl * (l - 1) + 1) =                      &
            &         grid%weights(nl + nl * (l - 2) + 1)             &
            &         + grid%gllw(0)
          else
            grid%weights(j + nl * (l - 1) + 1) = grid%gllw(j)
          end if
        end do
      end do
      grid%weights = grid%jac(1) * grid%weights

      !grid%maptype = para%maptype
      !grid%mapped = .false.

    endif

  end subroutine init_grid_dim_GLL
  
  !! @description: Initialize `nl+1` Legendre-Gauss-Lobatto  grid points in
  !!               the primitive interval $[-1,1]$
  !! @param: nl         Highest degree of Legendre Polynomials use for
  !!                    interpolation. The number of points in each element is
  !!                    `nl+1`
  !! @param: lobatto    Vector representing the `nl+1` points Gauss-Lobatto-
  !!                    Legendre in [-1,1]
  !! @param: weights    Gaussian weights associated to each of the `nl+1`
  !!                    points
  !! @param: tau        First order derivative matrix in $[-1,1]$
  subroutine get_lobatto_points(nl, lobatto, weights, tau)

    integer,                intent(in)  :: nl
    real(idp), allocatable, intent(out) :: lobatto(:)
    real(idp), allocatable, intent(out) :: weights(:)
    real(idp), allocatable, intent(out) :: tau(:,:)

    integer                             :: i, j, n_grid, n_eff, error, info
    real(idp)                           :: temp, aux_a, aux_b
    real(idp), allocatable              :: JJ(:,:), T(:,:), DCij(:,:)
    real(idp)                           :: vl(1,1), vr(1,1) ! never used
    real(idp), allocatable              :: work(:)
    real(idp), allocatable              :: wr(:), wi(:)
    integer                             :: lwork
    real(idp), allocatable              :: TT(:,:), trans_wc(:,:)

    n_grid = nl + 1
    n_eff  = n_grid - 2
    allocate(JJ(n_eff,n_eff), stat=error)
    call allocerror(error)
    allocate(wr(n_eff), stat=error)
    call allocerror(error)
    allocate(wi(n_eff), stat=error)
    call allocerror(error)

    JJ = zero

    do i = 2, n_eff - 1
      do j = 1, n_eff -1
         if (j .eq. i) then
          JJ(i,j-1) = real(i+1,idp)/real(2*i+1,idp)
          JJ(i,j+1) = real(i,idp)/real(2*i+1,idp)
         end if
      end do
    end do
    JJ(1,2) = third
    JJ(n_eff,n_eff-1) = real(1+n_eff,idp)/real(2*n_eff+1,idp)
    lwork = -1
    allocate(work(1), stat=error)
    call allocerror(error)
 
    call dgeev('N','N' , n_eff, JJ,n_eff, wr, wi, vl, n_eff, vr, n_eff, work,  &
    &           lwork, info )
    if (info /= 0) then
      write(*,*) "ERROR: Could not calculate optimal sizes for WORK array!"
    end if

    lwork = work(1) ! XXX possible change of value in conversion real -> int
    deallocate(work)
    allocate(work(lwork), stat=error)

    call allocerror(error)

    call dgeev('N','N', n_eff, JJ,n_eff, wr, wi, vl, n_eff, vr, n_eff, work,   &
    &           lwork, info )
    ! Reordering eigenvalues(Gauss-Lobatto-Legendre points) in ascending order
    deallocate(wi,work)
    do i = 1,n_eff
      do j = i,n_eff
       if (wr(i) >= wr(j)) then
        temp = WR(j)
        WR(j) = WR(i)
        WR(i) = temp
       end if
      end do
    end do
    allocate(lobatto(0:nl),stat=error)
    call allocerror(error)
    lobatto(0) = -one
    lobatto(nl) = one
    do i = 1,n_eff
    lobatto(i) = wr(i)
    end do
    deallocate(wr)
    ! Gaussian weights calculation
    ! Previous Legendre Polynomials calculation needed
    allocate(T(0:nl,0:nl),stat=error)
    call allocerror(error)

    do i = 0,n_grid-1
     do j = 0,n_grid-1
      if (i == 0) then
        T(i,j) = one
       elseif (i == 1) then
        T(i,j) = lobatto(j)
       else
        aux_a = real(2*(i-1)+1,idp)/real(i,idp)
        aux_b = real(i-1,idp)/real(i,idp)
          T(i,j) = aux_a*T(i-1,j)*lobatto(j) - aux_b*T(i-2,j)
      end if
     end do
    end do
    ! Computing Gaussian weights for legendre polynomials
    allocate(weights(0:nl), stat = error)
    call allocerror(error)
    allocate(TT(0:nl,0:nl),stat=error)
    call allocerror(error)
    TT(:,:) = transpose(T(:,:))

    do i = 0,n_grid - 1
      weights(i) = (real(n_grid,idp) - one)*(T(n_grid-1,i)*T(n_grid-1,i))
      weights(i) = weights(i) * half * real(n_grid,idp)
      weights(i) = one / weights(i)
    end do

    deallocate(T)
    allocate(DCij(0:nl,0:nl),stat=error)
    call allocerror(error)

    do i = 0,NL
     do j = 0,NL
       if ( (i == 0) .and. (j == 0)) then
         DCij(i,j) = -0.25_idp * real(nl, idp)*(real(nl,idp) + one)
        elseif ( (j == i) .and. (i == NL) ) then
         DCij(i,j) =  0.25_idp * real(nl,idp)*(real(nl, idp) + real(1,idp))
        elseif ((i > 0) .and. (i < nl) .and. (j == i)) then
         DCij(i,j) = zero
        else
         DCij(j,i) = TT(i,nl) / (TT(j,nl)*(lobatto(i)-lobatto(j)))
       end if
     end do
    end do
    deallocate(TT)

    allocate(trans_wc(0:nl,0:nl),stat=error)
    call allocerror(error)
    allocate(tau(0:nl,0:nl),stat=error)
    call allocerror(error)
    trans_wc = transpose(DCij)
    do i = 0,nl
      do j = 0,nl
       trans_wc(i,j) = weights(i)*trans_wc(i,j)
      end do
    end do
    tau = matmul(DCij,trans_wc)
    deallocate(trans_wc)

  end subroutine get_lobatto_points

  !! @description: Initialize the given cardinal-base work data structure for
  !!               the application of the kinetic operator in the cardinal
  !!               basis.
  !! @param: Tkin_carinal  Cardinal-base work data structure
  !! @param: grid          Spatial grid
  !! @param: mass          Mass involved in the kinetic operator
  subroutine init_work_cardinalbase(Tkin_cardinal, grid, mass)
    
    real(idp), allocatable,     intent(inout) :: Tkin_cardinal(:,:)
    type(grid_t),               intent(inout) :: grid
    real(idp),                  intent(in)    :: mass

    real(idp), allocatable :: temp(:,:)
    integer :: error, i, j, nl, m

    nl = grid%nl

    allocate(temp(grid%nl+1, size(grid%r)-2),stat=error)
    call allocerror(error)
    call get_kin_cardinal_banded_matrix(temp, grid, mass)
    allocate(Tkin_cardinal(nl+1, size(grid%r)-2),stat=error)
    call allocerror(error)

    do i = 1, nl+1
      do j = 1, size(temp(1,:))
        Tkin_cardinal(i,j) = temp(i,j)
      end do
    end do

  end subroutine init_work_cardinalbase
  
  !! @description: Initializing the mapped kinetic energy, written in
  !!               Legendre-Gauss-Lobatto interpolation polynomials basis
  !!               with bounded support in each element.
  !! @param: T     The kinetic energy matrix, stored in the Lapack storage
  !!               format for banded matrix.
  !! @param: grid  Coordinate grid corresponding to the mapped Legendre-
  !!               Gauss-Lobatto points.
  !! @param: mass  Mass in kinetic operator
  subroutine get_kin_cardinal_banded_matrix(T, grid,  mass)

    real(idp), allocatable, intent(inout) :: T(:,:)
    type(grid_t),           intent(in)    :: grid
    real(idp),              intent(in)    :: mass

    real(idp),    allocatable   :: tau(:,:)
    integer :: error, ku, u, v, a, b, th, l, nl,m
    integer :: i, j
    integer :: xu

    nl = grid%nl
    m  = grid%m
    ku = nl  !! Number of superdiagonals
    if (.not. allocated(T)) then
      allocate(T(nl+1, nl*m-1), stat=error)
      call allocerror(error)
    end if
    allocate(tau(0:nl,0:nl), stat=error)
    call allocerror(error)
    do i = 0, nl
      do j = 0, nl
        tau(i,j) = grid%D_primitive(i,j)
      end do
    end do
    T = real(0,idp)
    xu = 0
    xu = xu +1

    do l = 1, m
      do a = 0, nl
        do b = a, nl
           u = a + nl*(l-1)
           v = b + nl*(l-1)
           th = nl*(l-1) + 1
           if (((u < nl*m) .and. (u > 0)) .and. ((v < nl*m) .and. (v > 0))) then
             if ((l > 1) .and. (u .eq. nl*(l-1)) .and. (v .eq. nl*(l-1))) then
               T(ku+1+u-v,v) = (   tau(0,0)  / grid%jac(l)              &
               &                 + tau(nl,nl)/ grid%jac(l-1) )          &
               &               / grid%weights(th)
             else
               T(ku+1+u-v,v) =  (tau(a,b)/grid%jac(l))                  &
               &                / ( sqrt(grid%weights(th+a))            &
               &                   *sqrt(grid%weights(th+b)))
             end if
           end if
        end do
      end do
    end do
    T = T / (two*mass)

  end subroutine get_kin_cardinal_banded_matrix
  
  !! @description: Redefine the projection of operators in the modifiedy
  !!               Gaus-Lobatto-Legendre grid points to take into account
  !!               the homogenous Direchlet boundary conditions.
  !! @param: pot   Potential 
  subroutine redefine_ops_cardinal(pot)

    real(idp), allocatable, intent(inout) :: pot(:)

    real(idp), allocatable :: temp(:)
    integer :: error, i, j, k

    allocate(temp(size(pot)-2),stat=error)
    call allocerror(error)
    do j = 2, size(pot)-1
       temp(j-1) = pot(j)
    end do
    deallocate(pot)
    allocate(pot(size(temp)))
    do k = 1, size(temp)
       pot(k) = temp(k)
    end do
    temp = zero
    deallocate(temp,stat=error)

  end subroutine redefine_ops_cardinal
  
  !! @description: Initialize a simple 1D Legendre-Gauss-Lobatto grid between
  !!               `r_min` and `r_max` (exclusively) with `nl\times m-1` grid
  !!               points,\ where $nl$ and $m$ stand for the highest
  !!               polynomial interpolation in each element and the total number
  !!               of elements respectively, as a consequence of the homogeneous
  !!               Dirichlet boundary conditions at the edges of the global
  !!               grid.
  !! @param: grid  Grid variable to redefine.
  subroutine redefine_GLL_grid_1d(grid)

    type(grid_t), intent(inout)  :: grid

    real(idp), allocatable :: r(:), k(:), Jac(:), global_weights(:)
    integer                :: nr, error, j
    logical                :: J_alloc, w_alloc

    nr = size(grid%r) - 2

    allocate(r(nr), stat=error)
    call allocerror(error)
    allocate(k(nr), stat=error)
    call allocerror(error)
    allocate(Jac(nr), stat=error)
    call allocerror(error)
    allocate(global_weights(nr), stat=error)
    call allocerror(error)

    J_alloc = allocated(grid%J)
    w_alloc = allocated(grid%weights)

    do j = 1, nr
      r(j)                           = grid%r(j+1)
      if (J_alloc) Jac(j)            = grid%J(j+1)
      if (w_alloc) global_weights(j) = grid%weights(j+1)
    end do

    deallocate(grid%r)
    allocate(grid%r(nr), stat=error)
    call allocerror(error)
    if (J_alloc) then
      deallocate(grid%J)
      allocate(grid%J(nr), stat=error)
      call allocerror(error)
    end if
    if (w_alloc) then
      deallocate(grid%weights)
      allocate(grid%weights(nr), stat=error)
      call allocerror(error)
    end if

    do j = 1, nr
      grid%r(j)                    = r(j)
      if (J_alloc) grid%J(j)       = Jac(j)
      if (w_alloc) grid%weights(j) = global_weights(j)
    end do

    deallocate(r, Jac, global_weights)

  end subroutine redefine_GLL_grid_1d
  
  !! @description: Construct the explicit `nl+1` $\times$ `nr` matrix (because
  !!               of banded storage!) that is the Hamiltonian for a single
  !!               surface, consisting of the potential for that surface
  !!               (ignoring imaginary potentials, hence 'real'), and the
  !!               kinetic operator.
  !! @param: matrix       On output, explicit Hamiltonian matrix
  !! @param: grid         Coordinate Grid based on Gauss-Legendre-Lobatto points
  !! @param: Tkin_carinal Work array 
  subroutine get_real_surf_matrix_cardinal(matrix, grid, pot, Tkin_cardinal)

    real(idp), allocatable, intent(inout) :: matrix(:,:)
    type(grid_t),           intent(in)    :: grid
    real(idp), allocatable, intent(in)    :: pot(:) 
    real(idp), allocatable, intent(in)    :: Tkin_cardinal(:,:)

    integer                 :: nl, ku
    integer                 :: i, ir, error, j
    complex(idp)            :: cpulse
    real(idp)               :: rpulse

    ku = grid%nl

    if (allocated(matrix)) deallocate(matrix)
    allocate(matrix(ku+1, size(grid%r)), stat=error)
    call allocerror(error)

    ! kinetic operator
    do i = 1,size(Tkin_cardinal(:,1))
      do j = 1,size(Tkin_cardinal(1,:))
        matrix(i,j) = Tkin_cardinal(i,j)
      end do
    end do

    ku = grid%nl

    do ir = 1, size(pot)
      matrix(ku + 1,ir) = matrix(ku+1,ir) + pot(ir)
    end do

  end subroutine get_real_surf_matrix_cardinal
  
  !! @description: This subroutine returns the converged approximations to
  !! to the problem $\hat{A}_{red}z = \lambda z$ via a call to the `Arpack`
  !! routines `dsaupd` and `dseupd`.
  !! @param: matrix      Matrix that should be diagonalized, will be replaced
  !!                     by matrix of `nev` eigenvectors of dimension `n`
  !! @param: formt       Format of the input matrix to be diagonzalied. 'formt'
  !!                     *   `formt=full`: Full storage matrix format
  !!                     *   `format=banded`: Matrix stored in `Lapack` band
  !!                          storage format.
  !!                     could be one of the following:
  !! @param: n           Dimension of the eigenproblem
  !! @param: nev         Number of eigenvalues to be computed
  !! @param: which       Specifies which eigenvalues to be
  !!                     computed. `which` could be one of the following:
  !!                     *   `which=LA`: get the `nev` eigenvalues of largest
  !!                          amplitude.
  !!                     *   `which=SA': get the `nev` eigenvalues of smallest
  !!                          amplitude
  !!                     *   `which=BE`: Compute `nev` eigenvalues, half from
  !!                          each end of the spectrum.  When `nev` is odd,
  !!                          compute one more from the high end than from the
  !!                          low end.
  !! @param: eigenvals   Array of eigenvalues of the matrix containing the `nev`
  !!                     desired eigenvalues
  !! @param: rvec        Specifies wheather the associated `nev` eigenvector
  !!                     will be calculated or not.
  subroutine diag_arpack_real_sym_matrix(matrix, formt, n, nev, which,         &
  & eigenvals, rvec)

    real(idp), allocatable, intent(inout) :: matrix(:,:)
    character(len=*),       intent(in)    :: formt
    integer,                intent(in)    :: n
    integer,                intent(in)    :: nev
    character(len=2),       intent(in)    :: which
    real(idp), allocatable, intent(inout) :: eigenvals(:)
    logical,                intent(in)    :: rvec

    integer                :: ldv, mode, maxitr, lworkl
    integer                :: i, j, ido, error, iter
    real(idp)              :: sigma, tol
    character(len=1)       :: bmat
    integer                :: ncv, info
    integer                :: iparam(11)
    logical,   allocatable :: selct(:)
    real(idp), allocatable :: resid(:)
    real(idp), allocatable :: v(:,:)
    real(idp), allocatable :: workd(:)
    real(idp), allocatable :: workl(:)
    integer                :: ipntr(14)
    real(idp), allocatable :: d(:)

    ido    = 0
    info   = 0
    sigma  = zero
    tol    = zero
    ldv    = n
    ncv    = n
    lworkl = ncv**2 + 8*ncv
    maxitr = 900
    mode   = 1
    bmat   = 'I'
    iparam(1) = 1 !
    iparam(3) = maxitr
    iparam(4) = 1
    iparam(7) = mode

    allocate(eigenvals(nev),stat=error)
    call allocerror(error)
    allocate(v(ldv,ncv),stat=error)
    call allocerror(error)
    allocate(resid(n),stat=error)
    call allocerror(error)
    allocate(workd(3*n),stat=error)
    call allocerror(error)
    allocate(workl(lworkl),stat=error)
    call allocerror(error)
    allocate(selct(ncv),stat=error)
    call allocerror(error)
    allocate(d(nev),stat=error)
    call allocerror(error)

    iter = 0

    select case (formt)

      case ('full')

        do
          iter = iter + 1
          call dsaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,       &
          &           iparam, ipntr, workd, workl, lworkl, info)
          select case (abs(ido))
            case (1)
              ! compute  Y = OP * X  where
              ! IPNTR(1) is the pointer into WORKD for X,
              ! IPNTR(2) is the pointer into WORKD for Y.
              ! This is for the initialization phase to force the
              ! starting vector into the range of OP.
              call wrap_dsymv(n, one, matrix, workd(ipntr(1):ipntr(1)+n-1),    &
              &               zero, workd(ipntr(2):ipntr(2)+n-1))
              if (ido > 0) then       !no need to do this loop B=1
                ! compute Z = B * X for B = unity
                ! IPNTR(1) is the pointer into WORKD for X,
                ! IPNTR(3) is the pointer into WORKD for Z.
                do i = 1, n ! simply copy
                  workd(ipntr(3)+i-1) = workd(ipntr(1)+i-1)
                end do
              end if
            case (99)
              exit
            case default
              write(*,*) "ERROR: Unknown value for ido: ", ido
          end select
        end do

      case ('banded')

        do
          iter = iter + 1
          call dsaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,       &
          &           iparam, ipntr, workd, workl, lworkl, info)
          select case (abs(ido))
            case (1)
              ! compute  Y = OP * X  where
              ! IPNTR(1) is the pointer into WORKD for X,
              ! IPNTR(2) is the pointer into WORKD for Y.
              ! This is for the initialization phase to force the
              ! starting vector into the range of OP.
              call wrap_dsbmv(size(matrix(1,:)),size(matrix(:,1))-1, one,      &
              &               matrix, workd(ipntr(1):ipntr(1)+n-1), zero,      &
              &               workd(ipntr(2):ipntr(2)+n-1))
              if (ido > 0) then       !no need to do this loop B=1
                ! compute Z = B * X for B = unity
                ! IPNTR(1) is the pointer into WORKD for X,
                ! IPNTR(3) is the pointer into WORKD for Z.
                do i = 1, n ! simply copy
                  workd(ipntr(3)+i-1) = workd(ipntr(1)+i-1)
                end do
              end if
            case (99)
              exit
            case default
              write(*,*) "ERROR: Unknown value for ido: ", ido 
          end select
        end do

      case default

        write(*,*) "ERROR: Unknown format"

    end select ! formt

    if (allocated(matrix)) deallocate(matrix)
    allocate(matrix(n,nev), stat=error)
    call allocerror(error)

    do j = 1, nev
      do i = 1, n
       matrix(i,j) = v(i,j)
      end do
    end do

    call dseupd(rvec, 'A', selct, d, matrix, ldv, sigma, bmat, n, which, nev,  &
    &           tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl,  &
    &           info)

    if (.not. rvec) then
      do i = 1,nev
        eigenvals(i) = d(nev+1-i)
      end do
    else
      do i = 1,nev
        eigenvals(i) = d(i)
      end do
    end if

    deallocate(selct, d, v, workd)

  end subroutine diag_arpack_real_sym_matrix

  ! Symmetric matrix-vector multiplication
  !! @description: See BLAS
  !! [dsymv](http://www.netlib.org/blas/dsymv.f) documentation.
  !!
  !!  $\vec{y} = \alpha \hat{a} \vec{x} + \beta \vec{y}$
  !!
  !! @param: n      Dimension of `a`, `x`, `y`
  !! @param: alpha  Scalar for `x`
  !! @param: a      Matrix
  !! @param: x      Vector
  !! @param: beta   Scalar for `y`
  !! @param: y      Vector
  subroutine wrap_dsymv(n, alpha, a, x, beta, y)

    integer,       intent(in)    :: n
    real(idp),     intent(in)    :: alpha
    real(idp),     intent(in)    :: a(n,n)
    real(idp),     intent(in)    :: x(n)
    real(idp),     intent(in)    :: beta
    real(idp),     intent(inout) :: y(n)

    call dsymv('U', n, alpha, a, n, x, 1, beta, y, 1)

  end subroutine wrap_dsymv

  ! Symmetric banded matrix-vector multiplication
  !! @description: See BLAS
  !! [dsbmv](http://www.netlib.org/blas/dsbmv.f) documentation.
  !!
  !!  $\vec{y} = \alpha \hat{a} \vec{x} + \beta \vec{y}$
  !!
  !! @param: n      Dimension of `a`, `x`, `y`
  !! @param: ku     Number of super-diagonals
  !! @param: alpha  Scalar for `x`
  !! @param: a      Matrix (banded storage)
  !! @param: x      Vector
  !! @param: beta   Scalar for `y`
  !! @param: y      Vector
  subroutine wrap_dsbmv(n, ku, alpha, a, x, beta, y)

    integer,       intent(in)    :: n
    integer,       intent(in)    :: ku
    real(idp),     intent(in)    :: alpha
    real(idp),     intent(in)    :: a(ku+1,n)
    real(idp),     intent(in)    :: x(n)
    real(idp),     intent(in)    :: beta
    real(idp),     intent(inout) :: y(n)

    call dsbmv('U', n, ku, alpha, a, ku+1, x, 1, beta, y, 1)

  end subroutine wrap_dsbmv
  
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

end module dvr_diag_mod
