module DVRrhf

  use constants
  use RHFData
  use DVRData, only: para
  use dvr_diag_mod, only: allocerror
  use util_mod, only: stop_all

  implicit none

  contains
  
  subroutine DoRHF()

    use DVRData, only:  eigen_vecs, one_e_rad_int, two_e_rad_int, integrals_ang

    real(dp), allocatable :: MOCoeffs(:,:), Den(:,:), DenOld(:,:), F(:,:)
    real(dp), allocatable :: hcore(:,:)
    real(dp), allocatable :: Vred(:,:)
    integer, allocatable :: OrbInd(:,:,:)
    integer :: n_l, n_ml, n_nqn, error, indx, l_val, nTotOrbs
    integer :: i, j, k, maxit, iter
    real(dp) :: Energy, Del, start, finish

    write(iout,*) '***** Doing RHF Calculation using FE-DVR basis *****'

    n_l = para%l + 1
    ! Total number of spherical harmonics
    n_ml = n_l**2
    ! Total number of n_quantum number, which equal to the total number of radial grids
    n_nqn = para%ng

    maxit = 10

    !!!!!
    ! Here the total number product basis is calculated and linked to the 
    ! individual sets of (n,l,m) quantum numbers
    OrbInd = 0
    nTotOrbs = 0

    indx = 0
    l_val = para%l+1
!   write(*,*) 'limit for l_val:', l_val

    allocate(OrbInd(para%ng,para%l+1,2*para%l+1), stat=error)
    call allocerror(error)

    do i = 1, para%ng
      do j = 1, l_val
        do k = 1, 2*j-1
          indx = indx + 1
          OrbInd(i,j,k) = indx
!         write(*, *) i, j, k, OrbInd(i, j, k)
        end do
      end do
    end do
      
    nTotOrbs = indx
    if (nTotOrbs.ne.(n_nqn*n_ml)) call stop_all('DVRrhf','Total number of orbitals does not match')

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


    ! First expand the matrix of eigen vectors that we already get from solving 
    ! the radial Schroedinger into a product basis of R_{n,l}*Y{l,m}
    call ExpandBasis(eigen_vecs, MOCoeffs, OrbInd, n_l, n_nqn, nTotOrbs)

    ! Calculate the initial density which then be used to calculate the initial energy
    call GetDensity(Den, DenOld, MOCoeffs, nTotOrbs)

    ! Calculate the one-body part of the Fock matrix from the one-electron integrals stores in one_rad_int
    call CalcHCore(one_e_rad_int, hcore, OrbInd)

    call CalcVred(two_e_rad_int, integrals_ang, Den, Vred, OrbInd, nTotOrbs, n_nqn, n_l)

    call GetFock(hcore, Vred, F, nTotOrbs)

    call CalcEnergy(Den, hcore, F, nTotOrbs, Energy)
    
    write(iout, *) 'Energy at the beginning: ', Energy

    ! Iterations to converge RHF energy

    write(iout,*) 'iter      error          Energy     time'
    do iter = 1, maxit

      call cpu_time(start)

      ! First diagonalise the Fock matrix
      call DiagFock(F , MOCoeffs, nTotOrbs)

      ! Calculate the new density from updated MOCoeffs
      call GetDensity(Den, DenOld, MOCoeffs, nTotOrbs)

      ! Check whether the density is converged or not
      call CalcDelDen(Den, DenOld, nTotOrbs, Del)

      ! Calculate the new Fock matrix
      call CalcVred(two_e_rad_int, integrals_ang, Den, Vred, OrbInd, nTotOrbs, n_nqn, n_l)
      call GetFock(hcore, Vred, F, nTotOrbs)

      ! Calculate the energy at this iteration
      call CalcEnergy(Den, hcore, F, nTotOrbs, Energy)

      call cpu_time(finish)

      write(iout,'(i4,2f15.8,f10.5)') iter, Del, Energy, finish-start

    end do

    write(iout,*) '***** RHF is done...'

    deallocate(Den, DenOld, F, hcore, Vred, OrbInd)

  end subroutine DoRHF

  subroutine ExpandBasis(EigenVecs, MOCoeffs, OrbInd, n_l, n_nqn, n_COs)

    real(dp), allocatable, intent(in) :: EigenVecs(:,:,:)
    real(dp), allocatable, intent(inout) :: MOCoeffs(:,:)
    integer, allocatable, intent(in) :: OrbInd(:,:,:)
    integer, intent(in) :: n_l, n_nqn, n_COs

    integer :: i, n, l, m, indx_1, indx_2, error
    real(dp) :: val

    allocate(MOCoeffs(n_COs,n_COs), stat=error)

    call allocerror(error)

    indx_1 = 0
    indx_2 = 0
    do n = 1, n_nqn
      do i = 1, para%ng
        do l = 1, n_l
          val = EigenVecs(i,n,l)
          do m = 1, 2*l-1
            indx_1 = OrbInd(n,l,m)
            indx_2 = OrbInd(i,l,m)
            MOCoeffs(indx_1, indx_2) = val
!           write(77,'(6i4,f15.8)') n, l, m, i, indx_1, indx_2, MOCoeffs(indx_1,indx_2)
          end do
        end do
      end do
    end do

  end subroutine ExpandBasis

  subroutine GetDensity(Den, DenOld, MOCoeffs, n_COs)

    real(dp), allocatable, intent(inout) :: Den(:,:), DenOld(:,:)
    real(dp), allocatable, intent(in) :: MOCoeffs(:,:)
    integer, intent(in) :: n_COs

    integer :: i, j, nocc, iocc, nelec
    real(dp) :: val


    nelec =  para%Z
    
    if(mod(nelec,2).eq.0) then 
      nocc = nelec/2
    else
      call stop_all('GetDensity','RHF can not be done yet for open shels atoms')
    end if

    ! Before computing the new density copy the old density in the matrix Den
    ! into DenOld
    DenOld = Den

    Den = 0.0d0
    do i = 1, n_COs
      do j = 1, i
        val = 0.0d0
        do iocc = 1, nocc
          val = val + 2.0d0*MOCoeffs(iocc, i)*MOCoeffs(iocc, j)
        end do
        Den(i,j) = val
        Den(j,i) = val
!       write(78,*)  i, j, Den(i,j)
      end do
    end do

  end subroutine GetDensity

  subroutine CalcHCore(OneEInts, hcore, OrbInd)

    real(dp), allocatable, intent(in) :: OneEInts(:,:,:)
    real(dp), allocatable, intent(inout) :: hcore(:,:)
    integer, allocatable, intent(in) :: OrbInd(:,:,:)
    integer :: i, j, l, m, ind_1, ind_2
    real(dp) :: val

    do i = 1, para%ng
      do j = 1, i
        do l = 1, para%l+1 
          do m = 1, 2*l-1
            ind_1 = OrbInd(i, l, m)
            ind_2 = OrbInd(j, l, m)
            val = OneEInts(i,j,l)
            hcore(ind_1, ind_2) = val
            hcore(ind_2, ind_1) = val
!           write(76,*) ind_1, ind_2, hcore(ind_1, ind_2)
          end do
        end do
      end do
    end do

  end subroutine CalcHCore

  subroutine CalcVred(TwoERadInts, AngInts, Den, Vred, OrbInd, nTotOrbs, n_nqn, n_l)

    real(dp), allocatable, intent(in) :: TwoERadInts(:,:,:), Den(:,:)
    complex(idp), allocatable, intent(in) :: AngInts(:,:,:,:,:)
    real(dp), allocatable, intent(inout) :: Vred(:,:)
    integer, allocatable, intent(in) :: OrbInd(:,:,:)
    integer, intent(in) :: n_nqn, n_l, nTotOrbs

    integer :: k1, k2, l1, l2, l3, l4, m1, m2, m3, m4, lm1, lm2, lm3, lm4, nl2
    integer :: n_m1, n_m2, n_m3, n_m4, klm_1, klm_2, klm_3, klm_4, l, error
    complex(idp) :: int_value, val
    complex(idp), allocatable :: AngEls(:,:,:,:,:), inter_int(:,:,:,:)
    real(dp) :: start, finish

    nl2 = n_l**2

    call cpu_time(start)

    allocate(AngEls(2*para%l+1, nl2, nl2, nl2, nl2), stat=error)
    call allocerror(error)

    do l4 = 1, n_l
      n_m4 = 2*l4 - 1
      do m4 = 1, n_m4
        lm4 = (l4-1)**2 + m4
        do l3 = 1, n_l
          n_m3 = 2*l3 - 1
          do m3 = 1, n_m3
            lm3 = (l3-1)**2 + m3
            do l2 = 1, n_l
              n_m2 = 2*l2 - 1
              do m2 = 1, n_m2
                lm2 = (l2-1)**2 + m2
                do l1 = 1, n_l
                  n_m1 = 2*l1 - 1
                  do m1 = 1, n_m1
                    lm1 = (l1-1)**2 + m1

                    do l = 1, 2*para%l + 1
                        AngEls(l, lm1, lm2, lm3, lm4) =  AngInts(l, lm1, lm2, lm3, lm4) - &
                          & 0.5d0*AngInts(l, lm1, lm2, lm4, lm3)
                    end do

                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

    Vred = zero

    do lm4 = 1, nl2
      do lm3 = 1, nl2
        do lm2 = 1, nl2
          do lm1 = 1, nl2
            do l = 1, 2*para%l + 1
               val = AngEls(l, lm1, lm2, lm3, lm4)
               if (abs(val).gt.1e-12) then
                 do k2 = 1, n_nqn
                   klm_2 = (k2 -1)*nl2 + lm2 
                   klm_4 = (k2 -1)*nl2 + lm4 
                   do k1 = 1, n_nqn
                     klm_1 = (k1 -1)*nl2 + lm1
                     klm_3 = (k1 -1)*nl2 + lm3 
                     Vred(klm_1, klm_2) = Vred(klm_1, klm_2) + &
                     &  val * Den(klm_3, klm_4) * TwoERadInts(k1,k2,l)
                   end do
                 end do
               end if
            end do
          end do
        end do
      end do
    end do

!   do k1 = 1, n_nqn
!     do l1 = 1, n_l
!       n_m1 = 2*l1 - 1
!       do m1 = 1, n_m1
!         lm1 = (l1-1)**2 + m1
!         klm_1 = OrbInd(k1,l1,m1)

!         do k2 = 1, n_nqn
!           do l2 = 1, n_l
!             n_m2 = 2*l2 - 1
!             do m2 = 1, n_m2
!               lm2 = (l2-1)**2 + m2
!               klm_2 = OrbInd(k2,l2,m2)

!!              int_value = 0.0d0
!!              do l3 = 1, n_l
!!                n_m3 = 2*l3 - 1
!!                do m3 = 1, n_m3
!!                  lm3 = (l3-1)**2 + m3
!!                  klm_3 = OrbInd(k1,l3,m3)


!!                  do l4 = 1, n_l
!!                    n_m4 = 2*l4 - 1
!!                    do m4 = 1, n_m4
!!                      lm4 = (l4-1)**2 + m4
!!                      klm_4 = OrbInd(k2,l4,m4)

!!                      do l = 1, 2*para%l+1
!!                        int_value = int_value + (AngInts(l, lm1, lm2, lm3, lm4) - &
!!                        & 0.5d0*AngInts(l, lm1, lm2, lm4, lm3))* &
!!                        &  TwoERadInts(k1,k2,l)
!!                      end do

!!                    end do
!!                  end do
!! 
!!                end do
!!              end do

!!              Vred(klm_1, klm_2) = int_value
!               write(79,'(2i4,f15.8)') klm_1, klm_2, real(Vred(klm_1, klm_2))

!             end do
!           end do
!         end do

!       end do
!     end do
!   end do

    deallocate(AngEls)
    call cpu_time(finish)

!   write(iout, *) 'The reduced V is calculated in ', (finish-start), 'seconds...'

  end subroutine CalcVred

  subroutine GetFock(hcore, V, F, n_COs)

    real(dp), allocatable, intent(in)    :: hcore(:,:)
    real(dp), allocatable, intent(in)    :: V(:,:)
    real(dp), allocatable, intent(inout) :: F(:,:)
    integer,               intent(in)    :: n_COs

    integer  :: i, j 
    real(dp) :: start, finish, val

    call cpu_time(start)

    F = zero

    do j = 1, n_COs
      do i = 1, j
        val = hcore(i,j) + V(i,j)
        F(i,j) = val
        F(j,i) = val
!       write(79,*) i, j, F(i,j)
      end do
    end do

    call cpu_time(finish)

!   write(iout, *) 'Fock matrix is calculated in ', (finish-start), 'seconds...'

  end subroutine GetFock

  subroutine CalcEnergy(Den, hcore, F, nTotOrbs, En)

    real(dp), allocatable, intent(in)  :: Den(:,:), hcore(:,:), F(:,:)
    integer,               intent(in)  :: nTotOrbs
    real(dp),              intent(out) :: En

    integer  :: i, j
    real(dp) :: start, finish

    call cpu_time(start)

    En = 0.0d0
    do j = 1, nTotOrbs
      do i = 1, nTotOrbs
        En = En + 0.5d0*Den(i,j)*(hcore(i,j) + F(i,j))
      end do
    end do

    call cpu_time(finish)

!   write(iout, *) 'Energy is calculated in ', (finish-start), 'seconds...'

  end subroutine CalcEnergy

  subroutine DiagFock(F, MOCoeff, n_COs)

    use dvr_diag_mod, only : diag_matrix

    real(dp), allocatable, intent(inout) :: F(:,:), MOCoeff(:,:)
    integer, intent(in) :: n_COs

    real(dp), allocatable :: evals(:) 
    integer :: error, i, j

    allocate(evals(n_COs), stat=error)

    call diag_matrix(F, evals)
    
!   do i = 1, 20
!     do j = 1, 20
!       write(81, '(2i4, f15.8)') i, j, F(i,j)
!     end do
!   end do

!   do i = 1, n_COs
!     write(80,'(i4,19f15.8)') i, evals(i), (F(i,j),j=1,18)
!   end do

    MOCoeff = F

!   do i = 1, 20
!     do j = 1, 20
!       write(82, '(2i4, f15.8)') i, j, MOCoeff(i,j)
!     end do
!   end do

  end subroutine DiagFock

  subroutine CalcDelDen(Den, DenOld, nTotOrbs, Del)
    
    real(dp), allocatable, intent(in)  :: Den(:,:), DenOld(:,:)
    integer, intent(in) :: nTotOrbs
    real(dp), intent(out) :: Del

    integer :: i, j

    Del = 0.0d0
    do i = 1, nTotOrbs
      do j = 1, nTotOrbs
        Del = Del + (Den(i,j) - DenOld(i,j))**2
      end do
    end do

    Del = sqrt(Del/4.0d0)

  end subroutine CalcDelDen

end module DVRrhf
