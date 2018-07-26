module ReadInput

  use util_mod, only : get_free_unit, stop_all
  use input_mod
  use DVRUtils

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
      call DVRInput()
    case ("END")
      exit
    case default
      call stop_all('ReadInputMain', 'Keyword not recognized')
    end select
    end do

  end subroutine  ReadInputMain

end module ReadInput
