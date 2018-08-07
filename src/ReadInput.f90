module ReadInput

  use constants
  use util_mod, only : get_free_unit, stop_all
  use input_mod
  use InputData
  use DVRData, only: debug

  implicit none

  contains

  subroutine ReadInputMain(filename)

!   use DVRData

    character(64), intent(in) :: filename
    
    integer :: ir, ios, check
    logical :: tEof
    character(len=100)  :: w

    ir = get_free_unit()

    open (ir, File=filename, status='OLD', form= "FORMATTED", iostat = ios)

    write(iout, *) 'Reading the input from the file : ', trim(filename)

    call SetDVRInpDefaults()

    check = 0
    do 
      call read_line(tEof, ir)
      flush(6)
      if(tEof) exit
      call readu(w)
      if (w(1:1).eq.'!'.or.w(1:1).eq.'#') cycle
      select case(w)
      case("DVR")
        ! Read here the data defining DVR for the calculation
        call DVRInput(ir)
      case("ORBITAL")
        ! Read here the data related to orbitals
        orbital_ints = .true.
        call OrbitalInput(ir)
      case ("")
        check = check + 1
        if (check.gt.50) call report('input is empty', .true.)
      case ("END")
        exit
      case default
        call report('Keyword not recognized',.true.)
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
    l_max         = 3 !Rotational quantum number
    mass          = 1.0
    nev_fac       = 0.5d0

    mapped_grid   = .false.
    maptype       = 'diff'
    read_envelope = ''
    beta          = 0.015
    E_max         = 1d-5
    only_bound    = .true.
    
    debug = 1
    dvr_diag = .false.
    dvr_integrals = .false.
    trans_integrals = .false.

    n_max = 10

  end subroutine SetDVRInpDefaults

  subroutine DVRInput(ir)

    integer, intent(in) :: ir
    integer             :: check
    logical :: eof
    character (len=100)  :: w

    check = 0
    do
      call read_line(eof, ir)
      if (eof) then 
        exit
      end if  
      call readu(w)
      
      if (w(1:1).eq.'!'.or.w(1:1).eq.'#') cycle

      select case(w)
      
      ! Minimum and maximum limits for the radial distance
      case("RLIMITS")
        call getf(r_min,1.0d0)
        call getf(r_max,1.0d0)

      ! Number of Finite Elements grids
      case("N-GRIDS")
        call geti(m)
      ! Number of point in the Gauss-Lobato quadrature
      case("N-GL")
        call geti(nl)
      ! Maximum number of the angular momentum to be involved
      case("ROT-QN-MAX")
        call geti(l_max)
      ! The fraction of the total number of grid points to be 
      ! obtained as eigenvalue after the diagonalization
      case("NUM-EIG-VALUE-FAC")
        call getf(nev_fac)
      ! Mass of the Nucleus
      case("NUCLEAR-CHARGE")
        call geti(z)
      case("MASS")
        call getf(mass)
      case("DVR-DIAG")
        dvr_diag = .true.
      case("DVR-INTEGRALS")
        dvr_integrals = .true.
      case("DEBUG")
        call geti(debug)
      case("")
        check = check + 1
        if (check.gt.50) call report('many empty lines in the input')
      case("ENDDVR")
        exit
      case default
        call stop_all('DVRInput', 'Keyword not recognized')
      end select

    end do

  end subroutine DVRInput


  subroutine OrbitalInput(ir)

    integer, intent(in) :: ir
    logical :: eof
    character (len=100)  :: w

    do
      call read_line(eof, ir)
      if (eof) then 
        exit
      end if  
      call readu(w)

      if (w(1:1).eq.'!'.or.w(1:1).eq.'#') cycle

      select case(w)
      
      ! Number of Finite Elements grids
      case("NUM-N-QN")
        call geti(n_max)
      case("ENDORBITAL")
        exit
      case default
        call stop_all('DVRInput', 'Keyword not recognized')
      end select

    end do

  end subroutine OrbitalInput

end module ReadInput
