module dvr_diag_mod

  use dvr_spline_mod
  
  implicit none

  public
  
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
  !! @param: maptype         Type of mapping {'diff'|'int'}
  !! @param: diagtype        Region to diagonalise {'only_inner'|'only_outer'}
  !!                         (only if maptype = 'inner_outer')
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
    real(idp)                   :: mass
    real(idp)                   :: Z 
    real(idp)                   :: r_min
    real(idp)                   :: r_max
    real(idp)                   :: r_max1
    real(idp)                   :: r_max2
    character(len=pottype_l)    :: pottype
    character(len=file_l)       :: pot_filename 
    logical                     :: mapped_grid
    character(len=maptype_l)    :: maptype
    character(len=diagtype_l)   :: diagtype
    character(len=file_l)       :: read_envelope
    real(idp)                   :: beta
    real(idp)                   :: E_max
    real(idp)                   :: dr_max
    integer                     :: nr
    integer                     :: nl
    integer                     :: m
    integer                     :: m1
    integer                     :: m2
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
    real(idp)                   :: r_min, r_max, V_dr, min_dr, beta, E_max, &
    &                              r_max1, r_max2
    integer                     :: error, nl, m, l, j, nr_env, th, nr, m1, m2
    character(len=datline_l)    :: header
    character(len=file_l)       :: read_envelope
    character(len=maptype_l)    :: maptype
    character(len=maptype_l)    :: basis_for_mapping

    nr    = para%nr
    if (para%maptype == 'inner_outer') then
      nl    = para%nl
      m     = para%m1 + para%m2
      m1    = para%m1
      m2    = para%m2
      r_min = para%r_min
      r_max = para%r_max2
      r_max1= para%r_max1
      r_max2= para%r_max2
    else
      nl    = para%nl
      m     = para%m
      r_min = para%r_min
      r_max = para%r_max
    end if

    if (m * nl + 1 .ne. nr) then
      write(*,*) "ERROR: For the GLL grid, nr must equal m * nl + 1"
      stop
    end if

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
     
      ! Calculate the global weights and GLL points
      allocate(jac(m), stat = error)
      call allocerror(error)
      call get_lobatto_points(nl, grid%gllp, grid%gllw,      &
      &                       grid%D_primitive)
      allocate(lobatto(0:nl),stat=error)
      call allocerror(error)
      allocate(weights(0:nl),stat=error)
      call allocerror(error)

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
        
        case('inner_outer')
          allocate(X_mapped(nr), stat=error)
          call allocerror(error)
          !We do not need a Jacobian so we simply do not allocate it
          do l = 1, m1, 1
            do j = 0, nl-1, 1
              X_mapped(j + nl * (l - 1) + 1) =                                 &
              &   grid%gllp(j) * (r_max1 - r_min) / (real(2 * m1,idp))         &
              &   + r_min + (real(l,idp)- half) * (r_max1 - r_min)             &
              &                                                   / real(m1,idp)
            end do
          end do
          do l = m1+1, m1+m2, 1
            do j = 0, nl-1, 1
              X_mapped(j + nl * (l - 1) + 1) =                                 &
              &   grid%gllp(j) * (r_max2 - r_max1) / (real(2 * m2,idp))        &
              &   + r_max1 + (real(l-m1,idp)- half) * (r_max2 - r_max1) /      &
              &                                                     real(m2,idp)
            end do
          end do

          X_mapped(nl * m + 1) =                                               &
          &   grid%gllp(nl) * (r_max2 - r_max1) / (real(2 * m2,idp))           &
          &   + r_max1 + (real(m2,idp)- half) * (r_max2 - r_max1) /            &
          &                                                         real(m2,idp)

          write(*,*) X_mapped

       case default
         write(*,*) "ERROR: Cannot initialize mapped GLL grid! Unknown maptype."
         stop
      end select

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

    call file_shape(filename, n_rows, n_cols, non_data_lines=non_data_lines)
    if (n_rows == 0) then
      write(*,*) "ERROR: The file contains no data rows"
      stop
    end if
    if (n_cols == 0) then
      write(*,*) "ERROR: The file contains no data columns"
      stop
    end if

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
  !! @param: number_of_rows  Number of columns in the file
  !! @param: non_data_lines  If given, a sorted combination of comment_lines and
  !!                         blank_lines
  subroutine file_shape(filename, number_of_rows, number_of_cols,              &
  &                     non_data_lines)

    character(len=*),               intent(in)    :: filename
    integer,                        intent(out)   :: number_of_rows
    integer,              optional, intent(out)   :: number_of_cols
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
    
    if (present(number_of_cols)) then
      call file_columns(filename, columns)
      number_of_cols = columns
    end if

    if (associated(comments)) deallocate(comments)
    if (associated(blanks)) deallocate(blanks)
    if (associated(temp)) deallocate(temp)

  end subroutine file_shape
  
  !! @description: Calculate the number of columns in a data file
  !! @param: filename        File name
  !! @param: number_of_cols  Number of columns in the file
  subroutine file_columns(filename, number_of_cols)

    character(len=*),  intent(in)    :: filename
    integer, optional, intent(out)   :: number_of_cols

    integer :: ic ! character counter
    integer :: u
    character(len=1) :: c
    integer :: error
    integer :: noc ! number_of_cols
    logical :: is_comment
    logical :: is_blank
    logical :: searching_for_col

    open(newunit=u, file=filename, access='direct', action='READ',             &
    &    form='formatted', recl=1)
    ! Opening the file with direct formatted access means that we can read it
    ! character-by-character
    ic = 0
    searching_for_col = .true.
    is_blank = .true.
    is_comment = .false.
    noc = 0
    do
      ic = ic + 1
      if (ic == huge(ic)) then
        ! overflow is imminent
        write(*,*) "ERROR: Could not find columns by reading a sane amount "// &
        &          "of characters"
      end if
      read(u, '(A1)', rec=ic, iostat=error) c
      if (error /= 0) then
        exit ! end of line; we're done
      end if
      ! ASCII 32: space, 9: tab, 10: line-feed, 13: carriage-return
      if (ichar(c)/=32 .and. ichar(c)/=9 .and. ichar(c)/=10 .and.              &
      & ichar(c)/=13) then
        is_blank = .false.
        if (c == '#') then
          ! here, we assume that if the charcter '#' occurs anywhere in a line,
          ! the entire line is a comment.
          is_comment = .true.
        else
          if (searching_for_col) then
            ! We found a new column!
            noc = noc + 1
            searching_for_col = .false.
          end if
        end if
      elseif (ichar(c)==32 .or. ichar(c)==9) then
        ! we found a space, i.e. the delimiter between columns
        searching_for_col = .true.
      elseif (ichar(c)==9 .or. ichar(c)==10) then
        ! end of line
        if (.not. is_comment .and. .not. is_blank) then
          ! if the line was not a comment, we're done. number_of_cols (noc)
          ! should hold the correct value for the line we just finished
          exit
        else
          ! if the line was not a data line, we just reset everything and
          ! continue
          noc = 0
          is_comment = .false.
          is_blank = .true.
          searching_for_col = .true.
        end if
      end if
    end do

    if (present(number_of_cols)) number_of_cols = noc

    close(u, iostat=error)

  end subroutine file_columns
  
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
  subroutine get_real_surf_matrix_cardinal(matrix, grid, pot, Tkin_cardinal,   &
  &                                        two_index_variant, para)

    real(idp), allocatable, intent(inout) :: matrix(:,:)
    type(grid_t),           intent(in)    :: grid
    real(idp), allocatable, intent(in)    :: pot(:) 
    real(idp), allocatable, intent(in)    :: Tkin_cardinal(:,:)
    logical, optional,      intent(in)    :: two_index_variant
    type(para_t), optional, intent(in)    :: para

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
    
    if (present(two_index_variant)) then
      if (two_index_variant) then
        if (.not.(present(two_index_variant))) then
          write(*,*) 'get_real_surf_matrix_cardinal needs para if '//          &
          &          'the flag two_index_variant is set to .true.'
          stop
        end if
        do i = 1,size(Tkin_cardinal(:,1))
          do j = 1,size(Tkin_cardinal(1,:))
            matrix(i,j) = two * para%mass * matrix(i,j)
          end do
        end do
        do ir = 1, size(pot)
          matrix(ku + 1,ir) = matrix(ku+1,ir) + two * para%mass * pot(ir)
        end do
        return
      end if
    end if

    do ir = 1, size(pot)
      matrix(ku + 1,ir) = matrix(ku+1,ir) + pot(ir)
    end do

  end subroutine get_real_surf_matrix_cardinal
  
  !! @description: Convert banded storage matrix to full matrix, assuming the
  !! layout indicated by the combination of `mode`, `kl`, and `ku`.
  !!
  !! * `mode='g'` (general matrix) or `mode='t'` (triangular matrix):
  !!   convert directly, assuming `full(i,j)` is stored in
  !!   `packed(ku+1+i-j,j)` for `max(1,j-ku) <= i <= min(m, j+kl)
  !! * `mode='h'` (Hermitian matrix) or `mode='s'` (symmetric matrix):
  !!   If `kl=0`, assume information in `packed` describes the upper triangle
  !!   of the matrix. In this case `full(i,j)` is stored in
  !!   `packed(kd+1+i-j,j)` for `max(1,j-ku) <= i <= j` and `full(j,i)` is
  !!   either the complex conjugate of `full(i,j)` (`'h'`) or identical to
  !!   `full(i,j)` (`'s'`).
  !!   If `ku=0`, assume information in `packed` describes the lower triangle of
  !!   the matrix. `full(i,j)` is stored in `packed(1+i-j,j)` for
  !!   `j <= i <= min(n, j+kd)`. Again `full(j,i)` is again completed as
  !!   Hermitian or symmetric
  !!
  !! @param: full              Storage for full matrix. If already allocated,
  !!                           must be of size $n \times n$.
  !! @param: banded            Banded data, for input.
  !! @param: n                 Dimension of full matrix.
  !! @param: kl                Number of sub-diagonals. Must be 0 for Hermitian
  !!                           matrices
  !! @param: ku                Number of super-diagonals
  subroutine mat_banded_to_full(full, banded, n, kl, ku)

    real(idp), allocatable,    intent(inout) :: full(:,:)
    real(idp),                 intent(in)    :: banded(1+kl+ku,n)
    integer,                   intent(in)    :: n
    integer,                   intent(in)    :: kl
    integer,                   intent(in)    :: ku

    integer :: i, j, error, kd

    if (allocated(full)) then
      if ((lbound(full, 1) /= 1) .or.  (ubound(full, 1) /= n) .or.             &
      &   (lbound(full, 2) /= 1) .or.  (ubound(full, 2) /= n)) then
        write(*,*) 'ERROR: Full matrix is allocated to the wrong size.'
      end if
    else
      allocate(full(n,n), stat=error)
      call allocerror(error)
    end if

    full = czero

    ! See formula in http://www.netlib.org/lapack/lug/node124.html
    do j = 1, n
      do i = max(1, j-ku), min(n, j+kl)
        full(i,j) = full(i,j) + banded(ku+1+i-j,j)
      end do
    end do

  end subroutine mat_banded_to_full
  
  
  !! @description: Diagonalize the given real symmetric matrix via a call to the
  !!               Lapack routine `dsyevd`. The calculated eigenvectors are
  !!               saved in the columns of the matrix.
  !! @param: eigen_vecs  Matrix that should be diagonalized, will be replaced
  !!                     by matrix of eigenvectors
  !! @param: eigen_vals  Array of eigenvalues of the matrix
  subroutine diag_matrix(eigen_vecs, eigen_vals)

    real(idp),              intent(inout) :: eigen_vecs(:,:)
    real(idp), allocatable, intent(inout) :: eigen_vals(:)

    integer :: nn, lwork, liwork, error
    integer ,  allocatable :: iwork(:)
    real(idp), allocatable :: work(:)

    nn = size(eigen_vecs(:,1))

    if (allocated(eigen_vals)) deallocate(eigen_vals)
    allocate(eigen_vals(nn), stat=error)
    call allocerror(error)
    allocate (work(1), stat=error)
    call allocerror(error)
    allocate (iwork(1), stat=error)
    call allocerror(error)

    ! Perform workspace query: dsyevd only calculates the optimal sizes of the
    ! WORK and IWORK arrays, returns these values as the first entries of the
    ! WORK and IWORK arrays, cf. http://www.netlib.org/lapack/double/dsyevd.f
    lwork = -1; liwork = -1 ! indicate workspace query should be done
    call dsyevd('v', 'u', nn, eigen_vecs, nn, eigen_vals,                      &
    &           work, lwork, iwork, liwork, error)
    if (error /= 0) then
      write(*,*) 'ERROR: ' //                                                  &
      & "Could not calculate optimal sizes for WORK and IWORK arrays!"
    end if

    ! Now we can re-allocate WORK and IWORK with the optimal size, obtained from
    ! the first call to dsyevd
    lwork = work(1)
    liwork = iwork(1)
    deallocate(work, iwork)
    allocate(work(lwork), stat=error)
    call allocerror(error)
    allocate(iwork(liwork), stat=error)
    call allocerror(error)

    ! The second call to dsyevd performs the actual diagonalization
    call dsyevd('v', 'u', nn, eigen_vecs, nn, eigen_vals,                      &
    &           work, lwork, iwork, liwork, error)

    ! Check for erroneous exit status of DSYEVD
    if (error .lt. 0) then
      write(*,*) 'ERROR: ' //                                                  &
      &'A variable passed to the DSYEVD routine had an                         &
      &illegal entry!'
    elseif (error .gt. 0) then
      write(*,*) 'ERROR: '//                                                   &
      &'Routine DSYEVD failed to compute an eigenvalue!'
    end if

    deallocate(work, iwork)

  end subroutine diag_matrix

  
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
  
  !! @description: Diagonalize the given complex Hermitian matrix via a call to
  !!               the Lapack routine `zheevd`. The calculated eigenvectors are
  !!               saved in the columns of the matrix.
  !! @param: eigen_vecs  Matrix that should be diagonalized, will be replaced
  !!                     by matrix of (complex) eigenvectors
  !! @param: eigen_vals  Array of (real) eigenvalues of the matrix
  subroutine diag_hermitian_matrix(eigen_vecs, eigen_vals)

    complex(idp),              intent(inout) :: eigen_vecs(:,:)
    real(idp),    allocatable, intent(inout) :: eigen_vals(:)

    integer                   :: nn, lwork, lrwork, liwork, error
    integer ,     allocatable :: iwork(:)
    real(idp),    allocatable :: rwork(:)
    complex(idp), allocatable :: work(:)

    nn = size(eigen_vecs(:,1))

    if (allocated(eigen_vals)) deallocate(eigen_vals)
    allocate(eigen_vals(nn), stat=error)
    call allocerror(error)
    allocate (work(1), stat=error)
    call allocerror(error)
    allocate (iwork(1), stat=error)
    call allocerror(error)
    allocate (rwork(1), stat=error)
    call allocerror(error)

    ! Perform workspace query: zheevd only calculates the optimal sizes of the
    ! WORK, RWORK and IWORK arrays, returns these values as the first entries of
    ! the WORK, RWORK and IWORK arrays,
    ! cf. http://www.netlib.org/lapack/complex16/zheevd.f
    lwork = -1; liwork = -1; lrwork = -1 ! indicate workspace query
    call zheevd('v', 'u', nn, eigen_vecs, nn, eigen_vals,                      &
    &           work, lwork, rwork, lrwork, iwork, liwork, error)
    if (error /= 0) then
      write(*,*) "Could not calculate sizes for WORK, IWORK, RWORK arrays!"
    end if

    ! Now we can re-allocate WORK, IWORK, RWORK with the optimal size, obtained
    ! from the first call to zheevd
    lwork = work(1)
    liwork = iwork(1)
    lrwork = rwork(1)
    deallocate(work, iwork, rwork)
    allocate(work(lwork), stat=error)
    call allocerror(error)
    allocate(iwork(liwork), stat=error)
    call allocerror(error)
    allocate(rwork(lrwork), stat=error)
    call allocerror( error)

    ! The second call to zheevd performs the actual diagonalization
    call zheevd('v', 'u', nn, eigen_vecs, nn, eigen_vals,                      &
    &           work, lwork, rwork, lrwork, iwork, liwork, error)
    if (error /= 0) then
      write(*,*) "An argument of zheevd had an illegal entry!"
    end if

    deallocate(work, iwork, rwork)

  end subroutine diag_hermitian_matrix

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

  character(len=converted_l) function int2str(i, format)

    integer,                    intent(in) :: i
    character(len=*), optional, intent(in) :: format

    if (present(format)) then
      write(int2str, format) i
    else
      write(int2str, '(I25)') i
    end if
    int2str = adjustl(int2str)

  end function int2str
  
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
  
  !! @description: Returns value of radial hydrogen eigenfunction at given point 
  !! @param: r_val Value of radial coordinate 
  !! @param: n     Principal quantum number
  !! @param: l     Orbital angular momentum quantum number
  !! @param: an    Coefficient array (should be unallocated on first call) 
  real(idp) function radial_hydrogen_ef(r_val, n, l, an)
  
    real(idp), intent(in) :: r_val
    integer, intent(in)   :: n
    integer, intent(in)   :: l
    real(idp), allocatable :: an(:)
  
    integer :: i, k
    real(idp) :: val, dr, curr_r, drho, curr_rho
    
    integer, save :: prev_n = -1, prev_l = -1
    real(idp), save :: norm = - one

    if ((n .ne. prev_n) .or. (l .ne. prev_l)) then
      prev_n = n
      prev_l = l
      norm = zero
      if (allocated(an)) deallocate(an)
      allocate(an(0:n-l-1))
      an(0) = one
      do k = 0, n-l-2
        an(k+1) = an(k) * (real(k + l + 1 - n,kind=idp) /                      &
        &                  real((k + 1) * (k + 2*l + 2),kind=idp))
      end do
      do i = 1, 10001
        val = zero
        dr = 0.01_idp
        drho = (two / real(n, kind=idp)) * dr 
        curr_r = (i-1) * dr
        curr_rho = (i-1) * drho
        do k = 0, n-l-1
          val = val + an(k) * curr_rho**k *exp(-curr_rho/two)
        end do
        val = curr_rho**l * val
        norm = norm + (val**2 * dr * curr_r**2)
      end do
    end if

    curr_rho = (two / real(n, kind=idp)) * r_val
    val = zero
    do k = 0, n-l-1
      val = val + an(k) * curr_rho**k *exp(-curr_rho/two)
    end do
    val = curr_rho**l * val

    radial_hydrogen_ef = val / sqrt(norm)
  
  end function radial_hydrogen_ef 

end module dvr_diag_mod
