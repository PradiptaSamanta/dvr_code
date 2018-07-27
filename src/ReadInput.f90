module ReadInput

  use util_mod, only : get_free_unit, stop_all
  use input_mod
  use InputData

  implicit none

  contains

  subroutine ReadInputMain(filename)

!   use DVRData

    character(64), intent(in) :: filename
    
    integer :: ir, ios
    logical :: tEof
    character(len=100)  :: w

    ir = get_free_unit()

    open (ir, File=filename, status='OLD', form= "FORMATTED", iostat = ios)

    write(6,*) 'Reading the input file : ', trim(filename)

    call SetDVRInpDefaults()

    do 
    call read_line(tEof, ir)
    flush(6)
    if(tEof) exit
    call readu(w)
    write(6,*) trim(w)
    select case(w)
    case("DVR")
      ! Read here the data defining DVR for the calculation
      call DVRInput(ir)
    case ("END")
      exit
    case default
      call stop_all('ReadInputMain', 'Keyword not recognized')
    end select
    end do

  end subroutine  ReadInputMain

  subroutine SetDVRInpDefaults()

    ! Setting up the defauls for all the parameters which can further
    ! be modified after reading the input
    pottype       = 'file' ! 'density_file', 'file', 'analytical'
    pot_filename  = 'input_pot.in'
    r_min         = 0.0
    r_max         = 300.0
    m             = 200
    nl            = 5
    nr            = m * nl + 1
    l_max             = 3 !Rotational quantum number
    mass          = 1.0

    mapped_grid   = .false.
    maptype       = 'diff'
    read_envelope = ''
    beta          = 0.015
    E_max         = 1d-5
    only_bound    = .true.
    
    dvr_diag = .false.
    dvr_integrals = .false.
    trans_integrals = .false.

  end subroutine SetDVRInpDefaults

  subroutine DVRInput(ir)

    integer, intent(in) :: ir
    logical :: eof
    character (len=100)  :: w

    do
      call read_line(eof, ir)
      if (eof) then 
        exit
      end if  
      call readu(w)

      select case(w)
      
      ! Minimum and maximum limits for the radial distance
      case("RLIMITS")
        call getf(r_min,1.0d0)
        write(6,*) r_min
        call getf(r_max,1.0d0)
        write(6,*) r_max

      ! Number of Finite Elements grids
      case("N-GRIDS")
        call geti(m)
        write(6,*) m
      ! Number of point in the Gauss-Lobato quadrature
      case("N-GL")
        call geti(nl)
        write(6,*) nl
      ! Maximum number of the angular momentum to be involved
      case("ROT-QN-MAX")
        call geti(l_max)
        write(6,*) l_max
      ! Mass of the Nucleus
      case("NUCLEAR-CHARGE")
        call geti(z)
      case("MASS")
        call getf(mass)
        write(6,*) mass
      case("DVR-DIAG")
        dvr_diag = .true.
        write(6,*) dvr_diag
      case("DVR-INTEGRALS")
        dvr_integrals = .true.
      case("ENDDVR")
        exit
      case default
        call stop_all('DVRInput', 'Keyword not recognized')
      end select

    end do

  end subroutine DVRInput
end module ReadInput
