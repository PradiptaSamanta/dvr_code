module dvr_diag_mod
  
  implicit none

  private :: module_name

  public

  character(len=200) :: module_name = 'dvr_diag_mod'

contains

  subroutine dummy_dummy()

    write(*,*) module_name

  end subroutine dummy_dummy

end module dvr_diag_mod
