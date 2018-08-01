subroutine DVRCore()

  use constants
  use ReadInput, only : ReadInputMain
  use InputData, only : dvr_diag, dvr_integrals
  use DVRDiag, only : SetDVR, DVRDiagonalization
  use DVRIntRad, only : GetRadialElements
  use DVRIntAng, only : GetAngularElements
  use DVRData, only: iout

  implicit none

  character(64) :: file_in, file_out
  real(dp)  :: start, finish

  call cpu_time(start)

  file_in = 'dvr.inp'
  file_out = 'dvr.out'

  ! open the output file
  open (iout, File=file_out, status='OLD', form= "FORMATTED")

  ! Reading the input from an input file
  call ReadInputMain(file_in)

  ! Initializing the FE-DVR basis
  call SetDVR()

  if (dvr_diag) call DVRDiagonalization()

  if (dvr_integrals) then
    call GetRadialElements()
    call GetAngularElements()
  endif


  write(iout, *) 'The calculation is done. '
  
  call cpu_time(finish)
  write(iout,*) 'Time = ', finish-start, 'seconds.'
  close(iout)

end subroutine DVRCore
