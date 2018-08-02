module dvr_init_mod

  use constants
  use DVRData, only : grid_t, para_t
  use dvr_spline_mod
  use util_mod
  
  implicit none

! public
  
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
  subroutine init_grid_dim_GLL(grid, para)
    type(grid_t), intent(inout) :: grid
    type(para_t), intent(in)    :: para

    real(idp), allocatable      :: r_env(:), V_env(:), weight(:), jac(:),  &
    &                              X_mapped(:), lobatto(:), weights(:)
    real(idp)                   :: r_min, r_max, V_dr, min_dr, beta, E_max
    integer                     :: error, nl, m, l, j, nr_env, th, nr
    character(len=datline_l)    :: header
    character(len=file_l)       :: read_envelope
    character(len=maptype_l)    :: maptype
    character(len=maptype_l)    :: basis_for_mapping

    nr    = para%nr
    nl    = para%nl
    m     = para%m
    r_min = para%r_min
    r_max = para%r_max

    if (m * nl + 1 .ne. nr) then
      write(*,*) "ERROR: For the GLL grid, nr must equal m * nl + 1"
      stop
    end if

    write(iout, *) 'Initializing the GLL grid ...'

    if (para%mapped_grid) then ! Initialize mapped GLL grid

      beta = para%beta
      E_max = para%E_max
      maptype = para%maptype
      read_envelope = para%read_envelope

      ! sanity checks
      if (maptype .eq. 'diff') then
        basis_for_mapping = 'sin'
      elseif(maptype .eq. 'int') then
        basis_for_mapping = 'exp'
        write(*,*) "ERROR: Integral mapping not implemented."
        stop
      end if

      grid%nl = nl
      grid%m = m
      ! Get envelope potential
      ! This goes together with a constant spatial envelop grid r_env
      if (read_envelope /= '') then
        call read_ascii(r_env, V_env, read_envelope)
      else
        call min_dr_1d(para, min_dr)
        nr_env = floor((r_max - r_min) / min_dr) + 1
        !write(*,*) 'Calculating Envelope, initial nr_env:', nr_env
        if (nr_env < 1) then
          write(*,*) "********* nr_env =", nr_env
          write(*,*) "ERROR: Internal error in &
          &calculating the temporary grid, for mapping. Possibly some &
          &potential data files could not be read correctly."
        end if
        call get_spatial_coord(r_env, V_dr, weight, r_min, r_max, nr_env,      &
        &                      'const')
        call envelope_pot_1d(V_env, r_env, para, 1)
      end if
      ! Map the grid: calculate the coordinates r(:) and the jacobian J(:)
      select case (maptype)
        case('diff')
          call map_diff_1d(X_mapped, grid%J, r_min, r_max, r_env,     &
          &                V_env, nr, beta, para%mass, E_max,        &
          &                basis_for_mapping)
          grid%dr = one

        case('int')
          call map_int_1d(X_mapped, grid%J, r_min, r_max, r_env,      &
          &               V_env, nr, beta, para%mass, E_max,         &
          &               basis_for_mapping)
          deallocate(grid%J)

       case default
         write(*,*) "ERROR: Cannot initialize mapped GLL grid! Unknown maptype."
         stop
      end select
     
      ! Calculate the global weights and GLL points
      allocate(jac(m), stat = error)
      call allocerror(error)
      call get_lobatto_points(nl, grid%gllp, grid%gllw,      &
      &                       grid%D_primitive)
      allocate(lobatto(0:nl),stat=error)
      call allocerror(error)
      allocate(weights(0:nl),stat=error)
      call allocerror(error)

      do j = 0, nl
        lobatto(j) = grid%gllp(j)
        weights(j) = grid%gllw(j)
      end do
      allocate(grid%weights(size(X_mapped)), stat = error)
      call allocerror(error)
      allocate(grid%jac(m),stat=error)
      call allocerror(error)
      do l = 1,m
        th = nl*(l-1) + nl +1
        jac(l) = half*(X_mapped(th) - X_mapped(th-nl))
        grid%jac(l) = half*(X_mapped(th) - X_mapped(th-nl))
      end do

      do j = 0, nl
        grid%weights(j+1) = jac(1)*weights(j)
      end do
      do l=2,m
       do j = 0,nl
         if (j == 0) then
           grid%weights(j+nl*(l-1)+1)                                 &
           & = grid%weights(nl+nl*(l-2)+1) + jac(l)*weights(0)
          else
           grid%weights(j+nl*(l-1)+1) = jac(l)*weights(j)
         end if
        end do
       end do
      deallocate(weights,stat=error)
      call allocerror(error)

      if (allocated(grid%r)) deallocate(grid%r)
      call allocerror(error)

      allocate(grid%r(size(X_mapped)), stat=error)
      call allocerror(error)

      do l = 1,m
        th = nl*(l-1) + nl + 1
        do j = 0,nl
          grid%r(j+nl*(l-1) + 1 )                                     &
          & =   half*(X_mapped(th) - X_mapped(th-nl))* lobatto(j)              &
          &   + half*(X_mapped(th) + X_mapped(th-nl))
        end do
      end do

      ! Cleanup
      deallocate(r_env, V_env, lobatto, X_mapped)

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

    endif

  end subroutine init_grid_dim_GLL
  
  !! @description: Initialize an one dimensional operator via reading data from
  !! a given file. The file must contain two columns where the first column is
  !! treated as the x-axis and the second column as the y-axis.
  !! @param: op_a       One dimensional operator vector
  !! @param: r          Spatial grid
  !! @param: para       Parameter data 
  subroutine init_grid_op_file_1d(op_a, r, para)

    real(idp),        allocatable,  intent(out) :: op_a(:)
    real(idp),        allocatable,  intent(out) :: r(:)
    type(para_t),                   intent(in)  :: para

    integer                 :: i

    call read_ascii(r, op_a, para%pot_filename)
    do i = 1, size(r)-1
      if (r(i+1) < r(i)) then
        write(*,*) "ERROR: Input file must be sorted with an increasing grid!"
        stop
      end if
    end do
    !If input file is a density, convert to a potential
    if (para%pottype == 'density_file') then
      do i = 1, size(r)
        op_a(i) = op_a(i)
        !TODO Needs to be implemented 
      end do
    end if
    !Add rotational barrier
    do i = 1, size(r)
      op_a(i) = op_a(i) +                                                      &
      &         real(para%l * (para%l + 1), idp) / (two * para%mass * r(i)**2)
    end do

  end subroutine init_grid_op_file_1d
  
  subroutine init_grid_op_file_1d_all(op_a, r, para)

    real(idp),     allocatable,  intent(inout) :: op_a(:,:)
    real(idp),     allocatable,  intent(out) :: r(:)
    type(para_t),                   intent(in)  :: para

    real(idp), allocatable  :: tmp_op_a(:)
    integer                 :: i, l, l_val
    
    call read_ascii(r, tmp_op_a, para%pot_filename)
    do i = 1, size(r)-1
      if (r(i+1) < r(i)) then
        write(*,*) "ERROR: Input file must be sorted with an increasing grid!"
        stop
      end if
    end do
    allocate(op_a(size(r), para%l+1))
    !If input file is a density, convert to a potential
    if (para%pottype == 'density_file') then
      do i = 1, size(r)
        op_a(i,:) = tmp_op_a(i)
        !TODO Needs to be implemented 
      end do
    end if
    !Add rotational barrier
    do l = 1, para%l+1
      l_val = l-1
      do i = 1, size(r)
        op_a(i,l) = tmp_op_a(i) +                                                      &
        &         real(l_val * (l_val + 1), idp) / (two * para%mass * r(i)**2)
      end do
    end do

  end subroutine init_grid_op_file_1d_all
  
  !! @description: Load a given operator with it's own number of grid points
  !! onto a new set of grid points. This is done by interpolating the given
  !! operator via splining onto the new grid points.
  !! @param: new_r      New set of grid points
  !! @param: new_op_a   New operator data
  !! @param: old_r      Old set of grid points
  !! @param: old_op_a   Old operator data
  subroutine map_op(new_r, new_op_a, old_r, old_op_a)

    real(idp),   intent(in)  :: new_r(:)
    real(idp),   intent(out) :: new_op_a(:)
    real(idp),   intent(in)  :: old_r(:)
    real(idp),   intent(in)  :: old_op_a(:)

    type(spline_t)          :: spline
    integer                 :: i, error

    call init_spline(spline, old_r, old_op_a)
!   allocate(new_op_a(size(new_r)), stat=error)
!   call allocerror(error)
    do i = 1, size(new_r)
      call spline_value(spline, new_r(i), new_op_a(i))
    end do
    call delete_spline_t(spline)

  end subroutine map_op

  !! @description: Delete a variable of type `spline_t`
  !! @param: var        Variable to delete
  !! @param: thorough  If given as `.true.`, zero all components
  subroutine delete_spline_t(var, thorough)

    type(spline_t), intent(inout) :: var
    logical, optional, intent(in) :: thorough

    integer :: i1
    integer :: error
    logical :: l_thorough

    l_thorough=.false.
    if (present(thorough)) l_thorough=thorough

    if (allocated(var%cm)) deallocate(var%cm)
    if (allocated(var%x)) deallocate(var%x)
    if (allocated(var%y)) deallocate(var%y)
    if (l_thorough) then
      var%n = 0
    end if

  end subroutine delete_spline_t
  
  !! @description: Calculate the minimal `dr` for all the potentials in the same
  !! system as the `grid`
  !! @param: para    Parameter
  !! @param: min_dr  Resulting minimal `dr`
  subroutine min_dr_1d(para, min_dr)

    type (para_t),    intent(in)  :: para
    real(idp),        intent(out) :: min_dr

    integer                 :: j, k, nr, npot
    real(idp)               :: rmax, rmin
    real(idp),  allocatable :: pot(:), r(:)

    rmin = para%r_min
    rmax = para%r_max
    nr = para%nr

    min_dr = huge(zero)

    if (para%pottype == 'analytical') then 
      if (((rmax - rmin) / real(nr-1,idp)) < min_dr) then
        min_dr = (rmax - rmin) / real(nr-1,idp)
      end if
    elseif ((para%pottype == 'file') .or. (para%pottype == 'density_file')) then
      call init_grid_op_file_1d(pot, r, para)
      ! Check for lower dr in current potential file
      do k = 1, size(r)-1
        if ((r(k+1)-r(k)) < min_dr) then
          if (r(k+1) > r(k)) then
            min_dr = r(k+1) - r(k)
          else
            write(*,*) "ERROR: dr at grid point ", k, " in file "              &
            &          //trim(para%pot_filename)//" is negative"
          end if
        end if
      end do
      if (allocated(pot)) deallocate(pot)
      if (allocated(r)) deallocate(r)
    else
      write(*,*) "ERROR: Unknown pottype."
      stop
    end if

  end subroutine min_dr_1d
  
  !! @description: Calculate an envelope potential (defined on grid r) for all
  !! potentials in para that have the same label as `para%grid(gi)`. The
  !! envelope potential is calculated as the minimum of all potentials at any
  !! given position. Potentials loaded from file are spline-interpolated at the
  !! given positions. Note that the differentiability of the envelope potential
  !! cannot be guaranteed (e.g. if there are crossing potentials).
  !! @param: venv    Envelope potential
  !! @param: r       Grid vector of the envelope potential
  !! @param: para    Parameter
  !! @param: gi      Index of grid that decides the system
  subroutine envelope_pot_1d(venv, r, para, gi)

    real(idp),  allocatable, intent(out) :: venv(:)
    real(idp),               intent(in)  :: r(:)
    type(para_t),            intent(in)  :: para
    integer,                 intent(in)  :: gi

    type(spline_t)          :: spline
    real(idp),  allocatable :: pot(:), file_pot(:), file_r(:)
    real(idp)               :: newpotval
    integer                 :: i, j, ir, nr, error, npot

    nr = size(r)
    allocate (venv(nr), stat=error)
    call allocerror(error)
    venv = huge(zero)
    if (para%pottype == 'analytical') then
      allocate(pot(size(r)))
      do j = 1, size(r)
        pot(j) = analytical_potential(r(j), para%mass)
        if (pot(j) < venv(j)) then
          venv(j) = pot(j)
        end if
      end do
      if (allocated(pot)) deallocate(pot)
    elseif ((para%pottype == 'file') .or. (para%pottype == 'density_file')) then
      call init_grid_op_file_1d(file_pot, file_r, para)
      call init_spline(spline, file_r, file_pot)
      do j=1, size(r)
        call spline_value(spline, r(j), newpotval)
        if (newpotval < venv(j)) then
          venv(j) = newpotval
        end if
      end do
      call delete_spline_t(spline)
      if (allocated(file_pot)) deallocate(file_pot)
      if (allocated(file_r)) deallocate(file_r)
    else
      write(*,*) "ERROR: Unknown pottype."
      stop
    end if

    ! Modify envelopment potential so that it is monotonically increasing.
    ! By doing this the derivative of the Jacobian becomes continuous
    do ir=nr-1, 1, -1
      venv(ir) = min(venv(ir),venv(ir+1))
    end do

  end subroutine envelope_pot_1d
  
  !! @description:
  !! Initializes a set of real space grid points on the given interval
  !! $I=\left[r_{\text{min}},\,r_{\text{max}}\right]. If the given type is
  !! `const` the spatial points will be calculated using a constant
  !! stepsize. If the type is given as `exp` the grid points will be calculated
  !! for the grid $y$ with constant stepsize $\text{d}y$ to which the real
  !! spatial grid is mapped by $y=\log(r)$. The spatial grid points are then
  !! obtained by inverting the prevous formula and stored to the $r$-aray.
  !! @param: r       Vector representing the spatial $r$-grid.
  !! @param: dr      Stepsize of the spatial grid.
  !! @param: weights Weights of coordinate transformation.
  !! @param: r_min   Minimum $r$ value.
  !! @param: r_max   Maximum $r$ value.
  !! @param: nr      Number of grid points.
  !! @param: type    Maptype of the spatial grid. Allowed values are either
  !!                 `const` or `exp`.
  subroutine get_spatial_coord(r, dr, weights, r_min, r_max, nr, type)

    real(idp), allocatable, intent(out) :: r(:)
    real(idp),              intent(out) :: dr
    real(idp), allocatable, intent(out) :: weights(:)
    real(idp),              intent(in)  :: r_min
    real(idp),              intent(in)  :: r_max
    integer,                intent(in)  :: nr
    character(len=*),       intent(in)  :: type

    integer :: ir, error

    allocate(r(nr), stat=error)
    call allocerror(error)
    allocate(weights(nr), stat=error)
    call allocerror(error)
    r       = zero
    weights = zero

    if (nr .eq. 1) then
      dr         = one
      r(1)       = zero
      weights(1) = one
    else
      select case (type)
        case('const')
          dr = (r_max - r_min) / real(nr - 1,idp)
          do ir = 1, nr, 1
            r(ir) = r_min + real(ir - 1,idp) * dr
            weights(ir) = one
          end do

        case('exp')
          dr = (log(r_max) - log(r_min)) / real(nr - 1,idp)
          do ir=1, nr
            r(ir) = exp(log(r_min) + real(ir - 1,idp) * dr)
            weights(ir) = one
          end do

        case default
          write(*,*) "ERROR: Unable to initialize spatial grid! Unknown type!"
          stop
      end select
    endif

  end subroutine get_spatial_coord
  
  !! @description: Map the given grid vector r by the differential mapping
  !!               method
  !! @param: r      Grid vector r
  !! @param: J      Jacobian
  !! @param: r_min  Minimum r value
  !! @param: r_max  Maximum r value
  !! @param: r_env  Grid vector of the envelope potential
  !! @param: V_env  Potential vector of the envelope potential
  !! @param: nr     Number of grid points
  !! @param: beta   Mapping Parameter
  !! @param: mass   Mass of the molecule
  !! @param: E_max  Maximum energy cutoff
  !! @param: base   Base to be used in later propagation
  subroutine map_diff_1d(r, J, r_min, r_max, r_env, V_env, nr, beta, mass,     &
  &                     E_max, base)

    real(idp), allocatable, intent(out)   :: r(:)
    real(idp), allocatable, intent(out)   :: J(:)
    real(idp),              intent(in)    :: r_min
    real(idp),              intent(in)    :: r_max
    real(idp),              intent(in)    :: r_env(:)
    real(idp),              intent(in)    :: V_env(:)
    integer,                intent(inout) :: nr
    real(idp),              intent(in)    :: beta
    real(idp),              intent(in)    :: mass
    real(idp),              intent(in)    :: E_max
    character(len=*),       intent(in)    :: base

    type (spline_t) :: spline
    real (idp), allocatable :: x(:)
    real (idp) :: dr, V_eff
    integer :: ir, error

    allocate (r(nr), stat=error)
    call allocerror(error)
    allocate (J(nr), stat=error)
    call allocerror(error)

    !calculate r(:)
    r(1) = r_min
    call init_spline(spline, r_env, V_env) ! prepare interpolation V_env(r)
    do ir=2, nr
      call spline_value(spline, r(ir-1), V_eff) ! V_eff = V_env(r), via spline
      dr = beta * Pi / sqrt(two * mass * abs(E_max - V_eff))
      r(ir) = r(ir-1) + dr
    end do
    call delete_spline_t(spline)

    write (*,'(A61)') "    Original grid parameters were changed by &
                       &differential mapping"
    write (*,'(A53)') "               original                 after mapping"
    write (*,'(A11,2ES25.17)') "    r_min  ", r_min,       r_min
    write (*,'(A11,2ES25.17)') "    r_max  ", r_max,       r(nr)
    write (*,'(A11,I8, A20, I8)') "    nr     ", nr , " ", nr

    !if ((abs(r_max - r(nr)) / r_max) > 0.5_idp) then
    !  write(*,*) "ERROR: r_max differs by more than 50% of the original value."
    !  stop
    !end if

    ! calculate J(:)
    allocate (x(nr), stat=error)
    call allocerror(error)
    do ir=1, nr
      x(ir) = real(ir,idp)
    end do
    call init_spline(spline, x, r) ! prepare interpolation r(x)
    do ir=1, nr
      call getderivative(spline, x(ir), J(ir)) ! J = r'(x)
    end do
    call delete_spline_t(spline)
    deallocate (x)

  end subroutine map_diff_1d
  
  !! @description: Map the given grid vector r by the integral mapping
  !!               method
  !! @param: r      Grid vector r
  !! @param: J      Jacobian
  !! @param: r_min  Minimum r value
  !! @param: r_max  Maximum r value
  !! @param: r_env  Grid vector of the envelope potential
  !! @param: V_env  Potential vector of the envelope potential
  !! @param: nr     Maximum number of grid points
  !! @param: beta   Mapping Parameter
  !! @param: mass   Mass of the molecule
  !! @param: E_max  Maximum energy cutoff
  !! @param: base   Base to be used in later propagation
  subroutine map_int_1d(r, J, r_min, r_max, r_env, V_env, nr, beta, mass,      &
  &                     E_max, base)
    real(idp), allocatable, intent(out)   :: r(:)
    real(idp), allocatable, intent(out)   :: J(:)
    real(idp),              intent(in)    :: r_min
    real(idp),              intent(in)    :: r_max
    real(idp),              intent(in)    :: r_env(:)
    real(idp),              intent(in)    :: V_env(:)
    integer,                intent(inout) :: nr
    real(idp),              intent(in)    :: beta
    real(idp),              intent(in)    :: mass
    real(idp),              intent(in)    :: E_max
    character(len=*),       intent(in)    :: base

    real(idp), allocatable :: temp_r(:)
    real(idp), allocatable :: temp_J(:)
    real(idp) :: V_eff
    type (spline_t) :: V_spline
    real(idp) :: r_prev
    real(idp) :: dummy
    real(idp) :: extension_dr
    integer :: error
    integer :: i, k
    integer :: new_nr

    ! initialize
    allocate (temp_r(nr), stat=error)
    call allocerror(error)
    allocate (temp_J(nr), stat=error)
    call allocerror(error)

    call init_spline(V_spline, r_env, V_env) ! prepare interpolation V_env(r)

    call spline_value(V_spline, r_max, V_eff) ! V_eff = V_env(r_max), via spline
    temp_r(nr) = r_max
    temp_J(nr) = Pi/sqrt(two * mass * abs(E_max - V_eff))
    dummy = intmap_eq(r_min, r_min + beta, V_spline, E_max, beta, mass) ! init

    ! calculate r and J
    r_prev = r_max

    i = nr
    do while (r_prev >= r_min)
      i = i - 1
      if (i .le. 0) then
         write(*,*) "ERROR: nr is too small"
         stop
      end if
      r_prev = temp_r(i+1) - beta
      call find_bracket(r_prev, temp_r(i+1), V_spline, beta)
      r_prev = find_root(r_prev, temp_r(i+1), V_spline)
      temp_r(i) = r_prev
      call spline_value(V_spline, r_prev, V_eff)
      temp_J(i) =  Pi/sqrt(two * mass * abs(E_max - V_eff))
    end do
    !Once r becomes negative, continue linearly with last known spacing
    extension_dr = temp_r(i+1) - temp_r(i)
    do while (i > 1)
      i = i - 1
      r_prev    = temp_r(i+1) - extension_dr
      temp_r(i) = r_prev
      call spline_value(V_spline, r_prev, V_eff)
      temp_J(i) =  Pi/sqrt(two * mass * abs(E_max - V_eff))
    end do

    new_nr = nr

    allocate (r(new_nr), stat=error)
    call allocerror(error)
    allocate (J(new_nr), stat=error)
    call allocerror(error)

    do k = i, new_nr
      r(k) = temp_r(k)
      J(k) = temp_J(k)
    end do
    write (*,'(A61)') "    Original grid parameters were changed by &
                       &integral mapping"
    write (*,'(A53)') "               original                 after mapping"
    write (*,'(A11,2ES25.17)') "    r_min  ", r_min,       r_prev
    write (*,'(A11,2ES25.17)') "    r_max  ", r_max,       r(new_nr)
    write (*,'(A11,I8, A20, I8)') "    nr     ", nr , " ", new_nr

    !if ((abs(r_max - r(new_nr)) / r_max) > 0.5_idp) then
    !  write(*,*) "ERROR: r_max differs by more than 50% of the original value."
    !  stop
    !end if
    !if ((abs(real(nr - new_nr,idp)) / real(nr, idp)) > 0.5_idp) then
    !  write(*,*) "ERROR: nr differs by more than 50% of the original value."
    !  stop
    !end if

    nr = new_nr

  end subroutine map_int_1d
  
  !! @description: You need to call this at least once with the optional
  !!               parameters `E_max`, `beta`, and `mass`, for initialization,
  !!               or all hell will break loose
  !! @param: rl        Lower limit of integral
  !! @param: ru        Upper limit of integral
  !! @param: V_spline  `V(r)` in paper
  !! @param: E_max     `E_add` in paper (initialize once)
  !! @param: beta      Mapping parameter (initialize once)
  !! @param: mass      Mass (initialize once)
  real(idp) function intmap_eq(rl, ru, V_spline, E_max, beta, mass)

    real(idp),           intent(in) :: rl
    real(idp),           intent(in) :: ru
    type (spline_t),     intent(in) :: V_spline
    real(idp), optional, intent(in) :: E_max
    real(idp), optional, intent(in) :: beta
    real(idp), optional, intent(in) :: mass

    real(idp), save :: E_max_s
    real(idp), save :: beta_s
    real(idp), save :: mass_s

    real(idp)       :: V_eff, r, rstart
    real(idp)       :: fsum
    real(idp)       :: mapping_l, mapping_u
    real(idp)       :: qromb, dqromb
    real(idp)       :: del
    real(idp), allocatable       :: h(:), s(:)
    real(idp), allocatable       :: c(:), d(:), den(:)

    real(idp), parameter :: eps = 1.0e-6_idp
    integer, parameter :: jmax=60, jmaxp=jmax+1, k=5, km=k-1
    integer             :: i, j, it, ns, m, error

    if (present(E_max)) then
      E_max_s = E_max
      intmap_eq = zero
    end if
    if (present(beta)) then
      beta_s  = beta
      intmap_eq = zero
    end if
    if (present(mass)) then
      mass_s  = mass
      intmap_eq = zero
      return
    end if

    allocate (h(jmaxp), stat=error)
    call allocerror(error)
    allocate (s(jmaxp), stat=error)
    call allocerror(error)
    allocate (c(k), stat=error)
    call allocerror(error)
    allocate (d(k), stat=error)
    call allocerror(error)
    allocate (den(k), stat=error)
    call allocerror(error)

    if (ru < rl) then
      write(*,*) "ERROR: rl should be smaller then ru"
      stop
    end if
    if ( rl == ru) then
      intmap_eq = - beta_s
      return
    end if
    call spline_value(V_spline, rl, V_eff)
    mapping_l = sqrt( two * mass_s * (E_max_s-V_eff) ) / pi
    call spline_value(V_spline, ru, V_eff)
    mapping_u = sqrt( two * mass_s * (E_max_s-V_eff) ) / pi
    h(1) = one
    do j = 1, jmax
      if (j == 1) then
        s(j) = half * (ru-rl) * (mapping_l+ mapping_u)
      else
        it = 2**(j-2)
        del = (ru-rl) / real(it,idp) ! This is the spacing of the points to be
                                     ! added.
        rstart = rl+half*del
        r = rstart
        fsum = zero
        do i = 1, it
          call spline_value(V_spline, r, V_eff)
          fsum = fsum + sqrt( two * mass_s * (E_max_s-V_eff) ) / pi
          r = rstart +  del * real(i,idp)
        end do
        s(j) = half*(s(j) + del*fsum) ! This replaces s by its refined value.
      endif
     if (j >= K) then
       c = s(j-km:j) ; d = s(j-km:j)
       ns = minloc(abs(h(j-km:j)), dim=1)
       qromb = s(ns)
       ns = ns-1
       do m = 1, km
         den(1:k-m)=h(1:k-m)-h(1+m:k)
         if (any(den(1:k-m) == 0.0)) then
           write(*,*) "ERROR: calculation failure"
           stop
         end if
         den(1:k-m) = (c(2:k-m+1)-d(1:k-m)) / den(1:k-m)
         d(1:k-m) = h(1+m:k)*den(1:k-m)
         c(1:k-m) = h(1:k-m)*den(1:k-m)
         if (2*ns < k-m) then
           dqromb = c(ns+1)
         else
           dqromb = d(ns)
           ns = ns-1
         end if
         qromb = qromb + dqromb
       end do
       if (abs(dqromb) <= eps*abs(qromb)) exit
     end if
       s(j+1) = s(j)
       h(j+1) = 0.25d0*h(j)
    end do
    if (j == JMAX) then
      write(*,*) "ERROR: too many step in qromb"
    end if
    intmap_eq = qromb - beta_s

    deallocate(h); deallocate(s)
    deallocate(c); deallocate(d)
    deallocate(den)

    return

  end function intmap_eq
  
  !! @description: move down `xvar` so that `intmap_eq(x, xfix, V_spline)` has a
  !!               root in the interval $x = [xvar, xfix]$
  !! @param: xvar     Interval point that moves
  !! @param: xfix     Interval point that stays fixed
  !! @param: V_spline Spline for envelope potential
  !! @param: beta     beta
  subroutine find_bracket(xvar, xfix, V_spline, beta)

    real(idp),       intent(inout) :: xvar
    real(idp),       intent(in)    :: xfix
    type(spline_t),  intent(in)    :: V_spline
    real(idp),       intent(in)    :: beta

    integer, parameter :: max_tries = 100
    real(idp), parameter :: factor = 1.6_idp

    integer :: j
    real(idp) :: f1, f2

    f2 = -beta ! = intmap_eq(xfix, xfix, V_spline)
    f1 = intmap_eq(xvar, xfix, V_spline)
    do j = 1, max_tries
     if ( f1 * f2 < zero ) return
     xvar = xvar - factor*(xfix-xvar)
     if (xvar < zero) exit
     f1 = intmap_eq(xvar, xfix, V_spline)
    end do
    !Two possibilities to have not found success in the loop: Either max_tries
    !was exceeded, then an error will be returned, or the resulting xvar
    !is negative, in this case this value will be returned and treated in the
    !mapping routine itself
    if (xvar > zero) then
      write(*,*) "ERROR: Can't find bracket"
      stop
    end if

  end subroutine find_bracket
  
  !! @description: Assuming that there exists a root `r` so that
  !!               `intmap_eq(r, xfix, V_spline)` is zero and `r` is in the
  !!               interval $[r1:r2]$, find `r`.
  !!               This is implemented as Brent's alogrithm
  !! @param: r1       Lower limit of bracketed range
  !! @param: r2       Upper limit of bracketed range
  !! @param: V_spline Spline for envelope potential
  real(idp) function find_root(r1, r2, V_spline)

    real(idp),       intent(in)    :: r1
    real(idp),       intent(inout) :: r2
    type(spline_t),  intent(in)    :: V_spline

    real(idp), parameter :: TOL = 1.0d-8 ! get within TOL of the root
    real, parameter :: EPS = 1.0e-15_idp
    integer, parameter :: iter_max = 100

    integer :: iter
    real(idp) :: a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm, ur

    find_root = zero
    a = r1
    b = r2
    ur = r2
    fa = intmap_eq(a, ur, V_spline)
    fb = intmap_eq(b, ur, V_spline)

    if ( (fa > zero .and. fb > 0.) .or. (fa < zero .and. fb < zero) ) then 
      write(*,*) "ERROR: root must be bracketed"
      stop
    end if

    c = b
    fc = fb
    do iter = 1, iter_max
      if ( (fb > zero .and. fc > zero) .or. (fb < zero .and. fc < zero) )then
        c = a
        fc = fa
        d = b-a
        e = d
      endif
      if (abs(fc) < abs(fb)) then
        a = b
        b = c
        c = a
        fa = fb
        fb = fc
        fc = fa
      endif
      tol1 = two * EPS * abs(b) + half * TOL
      xm = half * (c - b)
      if (abs(xm) <= tol1 .or. fb == zero) then
        find_root = b
        return
      endif
      if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
        s = fb / fa
        if (a == c) then
          p = two * xm * s
          q = one - s
        else
          q = fa / fc
          r = fb / fc
          p = s * (two*xm*q*(q-r) - (b-a)*(r-one))
          q = (q-one)*(r-one)*(s-one)
        endif
        if (p > zero) q = -q
        p = abs(p)
        if( two*p < min( three*xm*q-abs(tol1*q), abs(e*q) ) ) then
          e = d
          d = p / q
        else
          d = xm
          e = d
        endif
      else
         d = xm
         e = d
      endif
      a = b
      fa = fb
      if(abs(d) > tol1) then
         b = b + d
      else
         b = b + sign(tol1, xm)
      endif
      fb = intmap_eq(b, ur, V_spline)
    end do
    write(*,*) "ERROR: root finding exceeding maximum iterations"

  end function find_root
  
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
      stop
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

    deallocate(temp)
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

    deallocate(tau)
  end subroutine get_kin_cardinal_banded_matrix
  
  !! @description: Redefine the projection of operators in the modifiedy
  !!               Gaus-Lobatto-Legendre grid points to take into account
  !!               the homogenous Direchlet boundary conditions.
  !! @param: pot   Potential 
  subroutine redefine_ops_cardinal(para, pot)

    type(para_t),            intent(in)   :: para
    real(idp), allocatable, intent(inout) :: pot(:,:)

    real(idp), allocatable :: temp(:,:)
    integer :: error, i, j, k

    allocate(temp(size(pot(:,1))-2,para%l+1),stat=error)
    call allocerror(error)
    do j = 2, size(pot(:,1))-1
       temp(j-1,:) = pot(j,:)
    end do
    deallocate(pot)

    allocate(pot(size(temp(:,1)), para%l+1))

    do k = 1, size(temp(:,1))
       pot(k,:) = temp(k,:)
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

    deallocate(r, k, Jac, global_weights)

  end subroutine redefine_GLL_grid_1d
  
  !! @description: Returns potential values of an analytical potential
  !! @param: r_val         Value of radial coordinate 
  real(idp) function analytical_potential(r_val, mass)
  
    real(idp), intent(in) :: r_val
    real(idp), intent(in) :: mass
  
    integer :: l
  
    l = 0 !Roational quantum number
  
    !Regularised -1/r potential with rotational barrier
    !Negative values for r_val are treated with a flat, positive constant
    !(Only relevant for integral mapping)
    if (r_val > zero) then
      analytical_potential = - r_val / (r_val**2 + 1d-16) +                    &
      &                      real(l*(l+1),idp) / (two * mass * r_val**2 + 1d-16)
    else
      analytical_potential = 1d16
    end if
  
  end function analytical_potential

end module dvr_init_mod
