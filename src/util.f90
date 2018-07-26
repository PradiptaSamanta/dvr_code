module util_mod

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
    write (*,'(a27,a)') ' The code DVR stops in subroutine: ',adjustl(sub_name)
    write (*,'(a9,18X,a)') 'Reason: ',adjustl(error_msg)
    write (*,'(a11)') 'EXITING...'

    ! Also push this to the stderr unit, so it hopefully ends up somewhere
    ! more useful.
    write (7,'(/a7)') 'ERROR.'
    write (7,'(a27,a)') ' The code DVR stops in subroutine: ',adjustl(sub_name)
    write (7,'(a9,18X,a)') 'Reason: ',adjustl(error_msg)
    write (7,'(a11)') 'EXITING...'

!   call print_backtrace_neci()

!    stop error_str
    stop 1

end subroutine stop_all

end module util_mod
