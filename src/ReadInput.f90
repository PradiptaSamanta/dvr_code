module ReadInput

  use constants
  use util_mod, only : get_free_unit, stop_all
  use input_mod
  use InputData
  use DVRData, only: debug, direct_2e, with_field, nFields, FieldComp, prim_integrals
  use RHFData, only: tRHF, n_rhf, maxit, DenTol

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
      case("RHF")
        ! Read here the data defining DVR for the calculation
        call RHFInput(ir)
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
    direct_2e = .false.
    prim_integrals = .false.

    n_max = 10
    two_e_int = 1
    nfrz = 0
    with_field = .false.
    nFields = 0

    tRHF = .false.
    n_rhf = 0
    maxit = 10

  end subroutine SetDVRInpDefaults

  subroutine DVRInput(ir)

    integer, intent(in) :: ir
    integer :: int_dummy
    integer             :: check
    logical :: eof
    character (len=100)  :: w
    character(len=32) :: Comp(3)

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
        r_max2 = r_max
      ! Grid point which separates B1 from B2 region 
      case("RMAX1")
        call getf(r_max1,1.0d0)
      ! Number of Finite Elements grids
      case("N-GRIDS")
        call geti(m)
      ! Number of Finite Elements in B1 region
      case("N-GRIDS1")
        call geti(m1)
      ! Number of Finite Elements in B2 region
      case("N-GRIDS2")
        call geti(m2)
        m = m1 + m2
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
      case("DIRECT-INTEGRALS")
        direct_2e = .true.
      case("DEBUG")
        call geti(debug)
      case("MAPPED-GRID")
        mapped_grid = .true.
        call getf(beta)
      case("MAP_INNER_OUTER")
        mapped_grid = .true.
        maptype = 'inner_outer'
      case("DIAGTYPE")
        call geti(int_dummy)
        if (int_dummy == 1) then
          diagtype = 'only_inner'
        elseif (int_dummy == 2) then
          diagtype = 'only_outer'
        else
          call stop_all('DVRInput', 'Invalid Integer for diagtype')
        end if
      case("FIELD")
        with_field = .true.
        if (nitems==1) then
          call stop_all('DVRInput', 'Please specify a component for the field')
        end if
        do while (item.lt.nitems)
          nFields = nFields + 1
          if (nFields.gt.3) then
            call stop_all('DVRInput','More than three fields are not allowed')
          end if
          call readu(Comp(nFields))
        end do
        allocate(FieldComp(nFields))
        FieldComp(1:nFields) = Comp(1:nFields)
      case("")
        check = check + 1
        if (check.gt.50) call report('many empty lines in the input')
      case("ENDDVR")
        exit
      case default
        stop
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
      case("TWO-E-INT")
        call geti(two_e_int)
      case("FROZEN")
        call geti(nfrz)
      case("PRIMITIVES")
        prim_integrals = .true.
      case("ENDORBITAL")
        exit
      case default
        call stop_all('OrbitalInput', 'Keyword not recognized')
      end select

    end do

  end subroutine OrbitalInput

  subroutine RHFInput(ir)

    integer, intent(in) :: ir
    logical :: eof
    character (len=100)  :: w

    tRHF = .true.

    do
      call read_line(eof, ir)
      if (eof) then 
        exit
      end if  
      call readu(w)

      if (w(1:1).eq.'!'.or.w(1:1).eq.'#') cycle

      select case(w)
      
      ! Number of n quntum to be stored from the Hartree-Fock solution
      case("NUM-N-QN")
        call geti(n_rhf)
      case("MAXITER")
        call geti(maxit)
      case("DENSITY-TOL")
        call getf(DenTol)
      case("ENDRHF")
        exit
      case default
        call stop_all('RHFInput', 'Keyword not recognized')
      end select

    end do

  end subroutine RHFInput

end module ReadInput
