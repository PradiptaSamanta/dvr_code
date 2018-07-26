subroutine DVRCore()

    use ReadInput, only : ReadInputMain
!   use dvr_spline_mod
!   use dvr_diag_mod
!   use dvr_diag


    character(64) :: filename

    filename = 'dvr.inp'

    call ReadInputMain(filename)

end subroutine DVRCore
