module DVRrhf

  use constants
  use RHFData
  use RHFMod
  use DVRData, only: para, grid, debug
  use dvr_diag_mod, only: allocerror
  use util_mod, only: stop_all, int2str

  implicit none

  contains
  
  subroutine DoRHF(ng, EigVec, OneInts, TwoInts)

    use DVRData, only:  integrals_ang, sph_harm

    integer, intent(in)   :: ng
    real(dp), allocatable, intent(inout) :: EigVec(:,:,:)
    real(dp), allocatable, intent(in) :: OneInts(:,:,:)
    real(dp), allocatable, intent(in) :: TwoInts(:,:,:)
    real(dp), allocatable :: MOCoeffs(:,:), Den(:,:), DenOld(:,:), F(:,:)
    real(dp), allocatable :: hcore(:,:)
    real(dp), allocatable :: TwoEInts(:,:,:,:)
    real(dp), allocatable :: Vred(:,:)
    real(dp), allocatable :: OrbEn(:)
    real(dp), allocatable :: OrbEn_l(:,:)
    real(dp), allocatable :: overlap(:,:)
    integer, allocatable :: OrbInd(:,:,:), n_m(:)
    integer :: n_l, n_ml, n_nqn, error, indx, l_val, nTotOrbs
    integer :: i, j, k, iter, l, m, occ
    integer :: vopt
    real(dp) :: Energy, Del, start, finish, tol, val
    logical :: RemoveL

    vopt = 1
    tol = DenTol
    RemoveL = .false.

    write(iout,*) '***** Doing RHF Calculation using FE-DVR basis *****'

    n_l = para%l + 1
    ! Total number of spherical harmonics
    n_ml = n_l**2
    ! Total number of n_quantum number, which equal to the total number of radial grids
    n_nqn = ng

    !!!!!
    ! Here the total number product basis is calculated and linked to the 
    ! individual sets of (n,l,m) quantum numbers


    allocate(OrbInd(ng,para%l+1,2*para%l+1), stat=error)
    call allocerror(error)
    allocate(OrbEn_l(ng,para%l+1), stat=error)
    call allocerror(error)

    OrbInd = 0
    nTotOrbs = 0

    indx = 0
    indx = 0
    
    allocate(n_m(n_l))
    n_m = sph_harm%n_m

    if (RemoveL) then
      do i = 1, ng
        do j = 1, min(n_l,i)
          !do k = 1, 2*j-1
          do k = 1, n_m(j)
            indx = indx + 1
            OrbInd(i,j,k) = indx
          end do
        end do
      end do
      nTotOrbs = indx
    else
      do i = 1, ng
        do j = 1, n_l
          !do k = 1, 2*j-1
          do k = 1, n_m(j)
            indx = indx + 1
            OrbInd(i,j,k) = indx
          end do
        end do
      end do
      nTotOrbs = indx
    end if

    write(iout,'(a,i6)')  ' Total number of spatial orbitals involved in the RHF calculations: ', nTotOrbs
!   if (nTotOrbs.ne.(n_nqn*n_ml)) call stop_all('DVRrhf','Total number of orbitals does not match')

    !!!!!

    ! allocate all the necessary matrices
    ! Den will store the density whereas DenOld will store those from the previous iteration
    allocate(Den(nTotOrbs,nTotOrbs), stat=error)
    call allocerror(error)
    allocate(DenOld(nTotOrbs,nTotOrbs), stat=error)
    call allocerror(error)
    ! F stores the Fock matrix
    allocate(F(nTotOrbs,nTotOrbs), stat=error)
    call allocerror(error)
    ! hcore stores the one-body part of the Fock matrix
    allocate(hcore(nTotOrbs,nTotOrbs), stat=error)
    call allocerror(error)
    allocate(Vred(nTotOrbs,nTotOrbs), stat=error)
    call allocerror(error)
    allocate(OrbEn(nTotOrbs), stat=error)
    call allocerror(error)
    allocate(MOCoeffs(nTotOrbs,nTotOrbs), stat=error)
    call allocerror(error)



    ! First expand the matrix of eigen vectors that we already get from solving 
    ! the radial Schroedinger into a product basis of R_{n,l}*Y{l,m}
    call ExpandBasis(EigVec, MOCoeffs, OrbInd, ng, n_l, n_nqn, n_m, &
    &                RemoveL)

!   do i = 1, ng
!     do l = 1, n_l
!       do m = 1, 2*l-1
!         indx = OrbInd(i,l,m)
!         if (abs(MOCoeffs(1,indx)).gt.1e-12) then
!           write(82,*) (l-1)**2 + m, grid%r(i), MOCoeffs(1,indx)
!         end if
!       end do
!     end do
!   end do

!   do i = 1, nTotOrbs
!     do j = 1, nTotOrbs
!       if (abs(MOCoeffs(i,j)).gt.1e-12) then
!         write(83, '(2i4, f15.8)') i, j, MOCoeffs(i,j)
!       end if
!     end do
!   end do

!   write(83, *) '##### Done'

    ! Calculate the initial density which then be used to calculate the initial energy
    call GetDensity(Den, DenOld, MOCoeffs,  nTotOrbs)

    ! Calculate the one-body part of the Fock matrix from the one-electron integrals stores in one_rad_int
    call CalcHCore(OneInts, hcore, OrbInd, ng, n_m, RemoveL)

!   if (vopt.eq.1) then

    call CalcVred(TwoInts, integrals_ang, Den, Vred, OrbInd, n_nqn, n_l, n_m)

!   else 
!     call Calc2ePrimOrbInts(TwoEInts, OrbInd, n_l, nTotOrbs)
!     call CalcVred_2(TwoEInts, Den, Vred, OrbInd, nTotOrbs, n_l)
!   end if

!   Vred = 0.0d0

    call GetFock(hcore, Vred, F, nTotOrbs)

    call CalcEnergy(Den, hcore, F, nTotOrbs, Energy)

    write(iout, *) 'Energy at the beginning: ', Energy

    ! Iterations to converge RHF energy

    write(iout,*) 'iter      error          Energy     time'

    do iter = 1, maxit

      call cpu_time(start)

      ! First diagonalise the Fock matrix
      call DiagFock(F , MOCoeffs, ng, n_m, OrbEn)

!     do i = 1, nTotOrbs
!       do j = 1, nTotOrbs
!         if (abs(MOCoeffs(i,j)).gt.1e-12) then
!           write(83, '(2i4, f15.8)') i, j, MOCoeffs(i,j)
!         end if
!       end do
!     end do

!     write(83, *) '##### Done'

      ! Calculate the new density from updated MOCoeffs
      call GetDensity(Den, DenOld, MOCoeffs, nTotOrbs)

      ! Check whether the density is converged or not
      call CalcDelDen(Den, DenOld, nTotOrbs, Del)

      ! Calculate the new Fock matrix
      if (vopt.eq.1) then
        call CalcVred(TwoInts, integrals_ang, Den, Vred, OrbInd, n_nqn, n_l, n_m)
      else
        call CalcVred_2(TwoEInts, Den, Vred, nTotOrbs)
      end if

      call GetFock(hcore, Vred, F, nTotOrbs)

      ! Calculate the energy at this iteration
      call CalcEnergy(Den, hcore, F, nTotOrbs, Energy)

      if (Del.lt.tol) then 

        call cpu_time(finish)

        write(iout,'(i4,2f15.8,f10.5)') iter, Del, Energy, finish-start

        write(iout,*) 'RHF calculation is converged in:', iter, 'iterations...'

        write(iout,'(a,f20.12,a)') ' Final RHF energy:', Energy, ' Hartree.'

        do i = 1, min(20,nTotorbs)

          if (2*i.le.para%Z) then
            occ = 2
          else
            occ = 0
          end if

          write(iout, '(a,i3,a,f20.12,a,i2)') ' MO #', i, 'energy= ',  OrbEn(i), ' occ= ', occ

        end do

        call ContractBasis(EigVec, MOCoeffs, OrbEn_l, OrbEn, OrbInd, ng, n_l, n_nqn)

        if (debug.gt.4) then
          do l = 1, n_l
            l_val = l - 1
            open(11, file="rhfenergies_GLL"//trim(int2str(l_val))//".dat", form="formatted",&
            &    action="write", recl=100000)
            do i = 1, ng
              write(11,*) i, OrbEn_l(i,l)
            end do
            close(11)
          end do
          do l = 1, n_l
            l_val = l - 1
            open(11, file="rhfvectors_GLL"//trim(int2str(l_val))//".dat", form="formatted",&
            &    action="write", recl=100000)
            do i = 1, size(EigVec(:,1,l))
              write(11,*) grid%r(i), (EigVec(i,j,l)/(sqrt(grid%weights(i)) * grid%r(i)), &
              & j = 1, size(EigVec(i,:,l)))
              !write(11,*) grid%r(i), (EigVec(i,j,l), &
              !& j = 1, size(matrix(i,:)))
            end do
            close(11)
          end do
        end if

        if (debug.gt.5) then

          allocate(overlap(para%nev,para%nev))
          do l = 1, n_l
            l_val = l - 1
            open(11, file="rhfoverlap"//trim(int2str(l_val))//".dat", form="formatted", &
            &    action="write")
  
            overlap = 0.0d0
            do i = 1, para%nev
              do j = 1, i
                val = 0.0d0
                do k = 1, size(grid%r)
                  val = val + EigVec(k,i,l)*EigVec(k,j,l)
                end do
                overlap(i,j) = val
                overlap(j,i) = val
                if (abs(val).gt.1e-12) &
                & write(11,'(2i5,f20.12)') i, j, val
              end do
            end do
            close(11)
          end do
          deallocate(overlap)
        end if
       
        exit

      else

        call cpu_time(finish)

        write(iout,'(i4,2f15.8,f10.5)') iter, Del, Energy, finish-start

        if (iter.eq.maxit) then
           write(iout,*) 'RHF calculation has not converged after ', maxit, 'iterations..'
        end if

      end if

    end do

    do j = 1, 10
      do i = 1, ng
        if (RemoveL) then
          l_val = min(i,n_l)
        else
          l_val = n_l
        end if
        do l = 1, l_val
          do m = 1, 2*l-1
            indx = OrbInd(i,l,m)
!           if (abs(MOCoeffs(j,indx)).gt.1e-12) then
!           write(81,'(3i4,f15.8)') i, l, m, MOCoeffs(j,indx)
!           end if
          end do
        end do
      end do
      write(81,*) '#### j=', j
    end do

    write(iout,*) '***** RHF is done...'

    deallocate(Den, DenOld, F, hcore, Vred, OrbInd, OrbEn, MOCoeffs, n_m)

  end subroutine DoRHF
end module DVRrhf
