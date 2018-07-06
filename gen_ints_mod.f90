module gen_ints_mod

  use dvr_spline_mod
  implicit none

  public

  type integrals_t

    integer :: n_b ! Total number of basis functions used in the integrals

    real(idp), allocatable   :: one_e_ints(:,:)
    real(idp), allocatable   :: two_e_ints(:,:,:,:)

  end type integrals_t
    
  contains

  subroutine init_integrals(integrals, nbasis)

    type(integrals_t), intent(inout)  :: integrals
    integer                           :: nbasis

    integer   :: error
    integrals%n_b = nbasis

    allocate(integrals%one_e_ints(nbasis, nbasis), stat=error)
    call allocerror(error)

    allocate(integrals%two_e_ints(nbasis, nbasis, nbasis, nbasis), stat=error)
    call allocerror(error)


  end subroutine init_integrals






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

end module gen_ints_mod
