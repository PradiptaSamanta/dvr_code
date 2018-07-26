module DVRUtils

  use constants
  use input_mod

  implicit none

  save 

  real(dp) :: r_min, r_max, mass, beta, e_max
  integer  :: m, nl, nr, l
  logical  :: mapped_grid, only_bound
  character(255) :: maptype, read_envelop, pottype, pot_filename, read_envelope


  contains

  subroutine SetDVRInpDefaults()

    pottype       = 'file' ! 'density_file', 'file', 'analytical'
    pot_filename  = 'input_pot.in'
    r_min         = 0.0
    r_max         = 300.0
    m             = 200
    nl            = 5
    nr            = m * nl + 1
    l             = 3 !Rotational quantum number
    mass          = 1.0

    mapped_grid   = .false.
    maptype       = 'diff'
    read_envelope = ''
    beta          = 0.015
    E_max         = 1d-5
    only_bound    = .true.

  end subroutine SetDVRInpDefaults

  subroutine DVRInput()

    logical :: eof
    character (len=100)  :: w

    do
      call read_line(eof)
      if (eof) then 
        exit
      end if  
      call readu(w)

      select case(w)

      end select

    end do

  end subroutine DVRInput

end module DVRUtils


