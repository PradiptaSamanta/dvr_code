module util_mod

  use constants
  use DVRData, only : datline_l, converted_l

  contains


  function get_free_unit() result(free_unit)

    ! Returns:
    !    The first free file unit above 10 and less than or equal to
    !    the paramater max_unit (currently set to 200).
    !
    !    If max_unit is exceeded, the function returns -1

    integer, parameter :: max_unit = 100
    integer :: free_unit
    integer :: i
    logical :: t_open, t_exist

    free_unit = -1
    do i = 10, max_unit
        inquire(unit=i, opened=t_open, exist=t_exist)
        if (.not.t_open .and. t_exist) then
            free_unit = i
            exit
        end if
    end do
    if (i == max_unit+1) call stop_all('get_free_unit','Cannot find a free unit below max_unit.')

  end function get_free_unit

  subroutine stop_all (sub_name, error_msg)

    ! Stop calculation due to an error. Exit with code 999?
    !
    ! In: sub_name    - Calling routine
    !     error_msg   - Error message

    implicit none

!   interface
!       subroutine print_backtrace_neci () bind(c)
!       end subroutine
!   end interface

    character(*), intent(in) :: sub_name, error_msg

    ! It seems that giving STOP a string is far more portable.
    ! MPI_Abort requires an integer though.
    character(3), parameter :: error_str='999'

    write (*,'(/a7)') 'ERROR.'
    write (*,'(a35,a)') ' The code DVR stops in subroutine: ',adjustl(sub_name)
    write (*,'(X,a,26X,a)') 'Reason: ',adjustl(error_msg)
    write (*,'(a11)') 'EXITING...'

    ! Also push this to the stderr unit, so it hopefully ends up somewhere
    ! more useful.
    write (7,'(/a7)') 'ERROR.'
    write (7,'(a35,a)') ' The code DVR stops in subroutine: ',adjustl(sub_name)
    write (7,'(a9,26X,a)') 'Reason: ',adjustl(error_msg)
    write (7,'(a11)') 'EXITING...'

!   call print_backtrace_neci()

!    stop error_str
    stop 1

end subroutine stop_all

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
      write(*,*) 'Debugging:', i
      write(*,'("")')
      write(*,'("================= ALLOCERROR =====================")')
      write(*,'("ERROR: Could not allocate memory")')
      stop
    end if

  end subroutine allocerror
  
end module util_mod
