module Density

    use constants
    use util_mod

    implicit none

    contains

    subroutine CalcDensity()

      use DensityData
    
      call ReadOneRDM(DensOrb1e, file_1rdm, tot_orb)

    end subroutine CalcDensity

    subroutine ReadOneRDM(Dens1e, file_1, tot_orb)

      complex(idp), allocatable, intent(inout) :: Dens1e(:)
      character(len=32), intent(in) :: file_1
      integer, intent(in) :: tot_orb

      integer :: n_elements, i, j, error, ij
      real(dp) :: val
      character(len=32) :: filename
      logical :: file_exists

      n_elements = tot_orb*(tot_orb + 1)/2

      allocate(Dens1e(n_elements))

      Dens1e = czero

      ! First read the real part of the one body reduced density matrix
      filename = trim(file_1)//"1"

      write(iout, *)  'Reading the real part of the 1-body RDM from the file ', trim(filename)
      inquire(file=trim(filename), exist=file_exists)

      if (file_exists) then
        open(12, file=trim(filename), form="formatted",&
        &    action="read", recl=100000)
      else
        call stop_all('ReadOneRDM', 'File for OneRDM is not present for reading')
      end if

      do
        read(12, *, iostat=error) i, j, val
        if  (error < 0) exit 
        if (i.le.j) then
          ij = i*(i-1)/2 + j 
          Dens1e(ij) = cmplx(val, zero)
          write(81, '(2i5, f25.17)') i, j, val
        else
          cycle
        endif
      end do

      close(12)

      ! Read the imaginary part of the one body reduced density matrix
      filename = trim(file_1)//"2"

      write(iout, *)  'Reading the imaginary part of the 1-body RDM from the file ', trim(filename)
      inquire(file=trim(filename), exist=file_exists)

      if (file_exists) then
        open(12, file=trim(filename), form="formatted",&
        &    action="read", recl=100000)
      else
        call stop_all('ReadOneRDM', 'File for OneRDM is not present for reading')
      end if

      do
        read(12, *, iostat=error) i, j, val
        if  (error < 0) exit 
        if (i.le.j) then
          ij = i*(i-1)/2 + j 
          Dens1e(ij) = cmplx(real(Dens1e(ij)), val)
          write(81, '(2i5, f25.17)') i, j, val
        else
          cycle
        endif
      end do

      close(12)

    end subroutine ReadOneRDM

end module Density
