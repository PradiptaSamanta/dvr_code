module dvr_diag_mod

  use dvr_spline_mod
  
  implicit none

  public
  
  !! @description: Parameters 
  !! @param: mass            Mass of particle (usually electron mass = 1 amu)
  !! @param: r_min           Minimum $r$
  !! @param: r_max           Maximum $r$
  !! @param: pottype         Type of potential {'analytical'|'file'}
  !! @param: maptype         Type of mapping {'diff'|'int'}
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
    character(len=pottype_l)    :: pottype
    character(len=file_l)       :: pot_filename 
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
      stop
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
  
  !! @description: Initialize an one dimensional operator via reading data from
  !! a given file. The file must contain two columns where the first column is
  !! treated as the x-axis and the second column as the y-axis.
  !! @param: op_a       One dimensional operator vector
  !! @param: r          Spatial grid
  !! @param: filename   Filename of the operator data
  subroutine init_grid_op_file_1d(op_a, r, filename)

    real(idp),        allocatable,  intent(out) :: op_a(:)
    real(idp),        allocatable,  intent(out) :: r(:)
    character(len=*),               intent(in)  :: filename

    integer                 :: i

    call read_ascii(r, op_a, filename)
    do i = 1, size(r)-1
      if (r(i+1) < r(i)) then
        write(*,*) "ERROR: Input file must be sorted with an increasing grid!"
        stop
      end if
    end do

  end subroutine init_grid_op_file_1d
  
  !! @description: Read two columns from an ascii file.
  !!               If there are not enough columns in the file, fill the
  !!               remaining ones with zero.
  !! @param: col1     Unallocated array for the first column
  !! @param: col2     Unallocated array for the second column
  !! @param: filename File name
  subroutine read_ascii(col1, col2, filename)

    real(idp), allocatable, intent(inout) :: col1(:), col2(:)
    character(len=*),       intent(in)    :: filename

    integer :: n_rows, n_cols
    integer :: error
    integer :: shape_table(2)
    real(idp), allocatable :: table(:,:)

    call read_ascii_table(table, filename)

    shape_table = shape(table)
    n_rows = shape_table(1)
    n_cols = shape_table(2)

    if (allocated(col1)) then
      if (size(col1) /= n_rows) deallocate(col1)
    end if
    if (.not. allocated(col1)) then
      allocate(col1(n_rows), stat=error)
      call allocerror(error)
    end if
    if (allocated(col2)) then
      if (size(col2) /= n_rows) deallocate(col2)
    end if
    if (.not. allocated(col2)) then
      allocate(col2(n_rows), stat=error)
      call allocerror(error)
    end if

    select case (n_cols)
    case (1)
      col1(:) = table(:,1)
      col2 = zero
    case (2:)
      col1(:) = table(:,1)
      col2(:) = table(:,2)
    end select

  end subroutine read_ascii
  
  !! @description: Read a table from an ascii file.
  !! @param: table       Unallocated table of reals
  !! @param: filename    Name of file from which data is read
  subroutine read_ascii_table(table, filename)

    real(idp), allocatable, intent(out) :: table(:,:)
    character(len=*),       intent(in)  :: filename

    integer :: row, line_nr
    integer :: error
    integer :: n_rows, n_cols
    integer, allocatable :: non_data_lines(:)
    integer :: cur_non_data_line
    real(idp), allocatable :: row_data(:)
    integer :: j ! pointer to where cur_non_data_line is in non_data_lines
    character(len=datline_l) :: line
    integer :: u

    call file_shape(filename, n_rows, non_data_lines=non_data_lines)
    if (n_rows == 0) then
      write(*,*) "ERROR: The file contains no data rows"
      stop
    end if
    n_cols = 2 

    allocate(table(n_rows, n_cols), stat=error)
    call allocerror(error)
    allocate(row_data(n_cols), stat=error)
    call allocerror(error)

    open(newunit=u, file=filename, action='READ', iostat=error)

    row = 1
    line_nr = 1
    j = 1
    do
      if (j > size(non_data_lines)) then
        cur_non_data_line = 0
      else
        cur_non_data_line = non_data_lines(j)
      end if
      if (line_nr == cur_non_data_line) then
        read(u, '(A)', iostat=error) line
        j = j + 1
        if (error < 0) exit ! End of File
      else
        ! we must use row_data as a temporary store. If we read directly
        ! from the file into table(row,:), we would get a segmentation fault
        ! when the EOF is reached.
        read(u, *, iostat=error) row_data(:)
        if (error < 0) exit ! End of File
        table(row,:) = row_data(:)
        row = row + 1
      end if
      line_nr = line_nr + 1
    end do

    close(u, iostat=error)
    deallocate(non_data_lines, row_data)

  end subroutine read_ascii_table
  
  !! @description: Return .true. if a line is a comment. Comments are lines
  !!               that start with '#', or the given `comment_char`. Leading
  !!               whitespace in the line is ignored.
  !! @param: line          String read from a file
  !! @param: comment_char  If given, character that indicates a comment.
  !!                       Defaults to `'#'`
  logical function is_comment(line, comment_char)

    character(len=*),           intent(in) :: line
    character(len=1), optional, intent(in) :: comment_char

    character(len=1) :: first_char

    is_comment = .false.
    first_char = adjustl(line)
    if (present(comment_char)) then
      if (first_char == comment_char) is_comment = .true.
    else
      if (first_char == "#") is_comment = .true.
    end if

  end function is_comment
  
  !! @description: Calculate shape of a data file, i.e. the number
  !!               of rows and the number of columns in the file.
  !!               Columns are separated by spaces. Empty lines or
  !!               lines starting with '#' are not counted.
  !! @param: filename        File name
  !! @param: number_of_rows  Number of rows in the file
  !! @param: non_data_lines  If given, a sorted combination of comment_lines and
  !!                         blank_lines
  subroutine file_shape(filename, number_of_rows, non_data_lines)

    character(len=*),               intent(in)    :: filename
    integer,                        intent(out)   :: number_of_rows
    integer, allocatable, optional, intent(inout) :: non_data_lines(:)

    logical :: ex
    integer :: error, i, j, k
    integer :: u
    integer :: c_i, b_i ! index for comments, blanks
    integer :: line_nr
    integer :: columns
    character(len=datline_l) :: line
    integer, pointer :: comments(:)
    integer, pointer :: blanks(:)
    integer, pointer :: temp(:)

    allocate(comments(5), stat=error)
    call allocerror(error)
    comments = 0
    allocate(blanks(5), stat=error)
    call allocerror(error)
    blanks = 0
    nullify(temp)

    inquire(file=filename, exist=ex)
    if (.not. ex) then
      write(*,*) "ERROR: Cannot find file: "//trim(filename)
      stop
    end if

    open(newunit=u, file=filename, action='READ', iostat=error)
    number_of_rows = 0
    c_i = 0
    b_i = 0
    line_nr = 0
    loop_over_lines: do
      read(u, '(A)', iostat=error) line
      line_nr = line_nr + 1
      if (error < 0) then
        exit loop_over_lines
      end if
      if (is_comment(line)) then
        c_i = c_i + 1
        if (c_i > size(comments)) then
          ! resize
          allocate(temp(2*size(comments)), stat=error)
          call allocerror(error)
          temp = 0
          temp(1:size(comments)) = comments(:)
          deallocate(comments)
          comments => temp
          nullify(temp)
        end if
        comments(c_i) = line_nr
      elseif (line == '') then
        b_i = b_i + 1
        if (b_i > size(blanks)) then
          ! resize
          allocate(temp(2*size(blanks)), stat=error)
          call allocerror(error)
          temp = 0
          temp(1:size(blanks)) = blanks(:)
          deallocate(blanks)
          blanks => temp
          nullify(temp)
        end if
        blanks(b_i) = line_nr
      else
        number_of_rows = number_of_rows + 1
      end if
    end do loop_over_lines
    close(u, iostat=error)

    if (present(non_data_lines)) then
      allocate(non_data_lines(b_i+c_i), stat=error)
      call allocerror(error)
      j = 1
      k = 1
      do i = 1, size(non_data_lines)
        if (((comments(j) < blanks(k)) .or. (k > b_i)) .and. (j <= c_i)) then
          non_data_lines(i) = comments(j)
          j = j + 1
        else
          non_data_lines(i) = blanks(k)
          k = k + 1
        end if
      end do
    end if

    if (associated(comments)) deallocate(comments)
    if (associated(blanks)) deallocate(blanks)
    if (associated(temp)) deallocate(temp)

  end subroutine file_shape
  
  !! @description: Load a given operator with it's own number of grid points
  !! onto a new set of grid points. This is done by interpolating the given
  !! operator via splining onto the new grid points.
  !! @param: new_r      New set of grid points
  !! @param: new_op_a   New operator data
  !! @param: old_r      Old set of grid points
  !! @param: old_op_a   Old operator data
  subroutine map_op(new_r, new_op_a, old_r, old_op_a)

    real(idp),              intent(in)  :: new_r(:)
    real(idp), allocatable, intent(out) :: new_op_a(:)
    real(idp),              intent(in)  :: old_r(:)
    real(idp),              intent(in)  :: old_op_a(:)

    type(spline_t)          :: spline
    integer                 :: i, error

    call init_spline(spline, old_r, old_op_a)
    allocate(new_op_a(size(new_r)), stat=error)
    call allocerror(error)
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
              stop
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
              stop 
          end select
        end do

      case default

        write(*,*) "ERROR: Unknown format"
        stop

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
