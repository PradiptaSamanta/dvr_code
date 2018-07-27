subroutine DVRCore()

  use ReadInput, only : ReadInputMain
  use InputData, only : dvr_diag, dvr_integrals
  use DVRDiag, only : SetDVR, DVRDiagonalization
  use DVRIntRad, only : GetRadialElements
  use DVRIntAng, only : GetAngularElements

  character(64) :: filename

  filename = 'dvr.inp'

  ! Reading the input from an input file
  call ReadInputMain(filename)

  ! Initializing the FE-DVR basis
  call SetDVR()

  if (dvr_diag) call DVRDiagonalization()

  if (dvr_integrals) then
    call GetRadialElements()
    call GetAngularElements()
  endif

end subroutine DVRCore
