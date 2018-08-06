subroutine DVRCore()

  use constants
  use ReadInput, only : ReadInputMain
  use InputData, only : dvr_diag, dvr_integrals, orbital_ints
  use DVRDiag, only : SetDVR, DVRDiagonalization
  use DVRIntRad, only : GetRadialElements
  use DVRIntAng, only : GetAngularElements
  use OrbInts, only  : GetOrbInts

  implicit none

  character(64) :: file_in, file_out
  real(dp)  :: start, finish

  call cpu_time(start)

  file_in = 'dvr.inp'
  file_out = 'dvr.out'

  ! open the output file
  open (iout, File=file_out, status='UNKNOWN', form= "FORMATTED")

  ! Reading the input from an input file
  call ReadInputMain(file_in)

  ! Initializing the FE-DVR basis
  call SetDVR()

  if (dvr_diag) call DVRDiagonalization()

  if (dvr_integrals) then
    call GetRadialElements()
    call GetAngularElements()
  endif

  if (orbital_ints) then
    call GetOrbInts()
  end if

  call DeallocMatrices()

  write(iout, *) 'The calculation is done. '
  
  call cpu_time(finish)
  write(iout,*) 'Time = ', finish-start, 'seconds.'
  close(iout)

end subroutine DVRCore

subroutine DeallocMatrices()
   
  use constants, only : iout
  use DVRData, only : eigen_vals, eigen_vecs, one_e_rad_int, &
    &                 two_e_rad_int, integrals_ang
  use OrbData, only : SpatialOrbInd, OneEInt, TwoEInt
  
  if (allocated(eigen_vals)) deallocate(eigen_vals)
  if (allocated(eigen_vecs)) deallocate(eigen_vecs)
  if (allocated(one_e_rad_int)) deallocate(one_e_rad_int)
  if (allocated(two_e_rad_int)) deallocate(two_e_rad_int)
  if (allocated(integrals_ang)) deallocate(integrals_ang)

  if (allocated(SpatialOrbInd)) deallocate(SpatialOrbInd)
  if (allocated(OneEInt)) deallocate(OneEInt)
  if (allocated(TwoEInt)) deallocate(TwoEInt)
  
  write(iout, *) 'Deallocated all the main matrices.'

end subroutine DeallocMatrices

