module RHFMod

  use constants
  use RHFData
  use DVRData, only: para, grid, debug
  use dvr_diag_mod, only: allocerror
  use util_mod, only: stop_all, int2str

  implicit none

  contains

  subroutine ExpandBasis(EigenVecs, MOCoeffs, OrbInd, ng, n_l, n_nqn, n_m, &
    &                RemoveL)

    real(dp), allocatable, intent(in) :: EigenVecs(:,:,:)
    real(dp), allocatable, intent(inout) :: MOCoeffs(:,:)
    integer, allocatable, intent(in) :: OrbInd(:,:,:)
    integer, intent(in) :: ng, n_l, n_nqn
    integer, allocatable, intent(in) :: n_m(:)
    logical, intent(in) :: RemoveL

    integer :: i, n, l, m, indx_1, indx_2, l_val
    real(dp) :: val

    !allocate(MOCoeffs(nTotOrbsNew,nTotOrbsOld), stat=error)
    indx_1 = 0
    indx_2 = 0
    do n = 1, n_nqn
      do i = 1, ng
        if (RemoveL) then
          l_val =  min(i,n,n_l)
        else 
          l_val = n_l
        end if
        do l = 1, l_val
          !val = EigenVecs(i,n,l)/(sqrt(grid%weights(i)) * grid%r(i))
          val = EigenVecs(i,n,l)
          !do m = 1, 2*l-1
          do m = 1, n_m(l)
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

!   do iocc = 1, nocc
!     norm = 0.0d0
!     do i = 1, n_COs
!       norm = norm + MOCoeffs(iocc, i)*MOCoeffs(iocc, i)
!     end do
!     write(iout,*) 'Norm: ', iocc, norm
!   end do

    ! Before computing the new density copy the old density in the matrix Den
    ! into DenOld
    DenOld = Den

    Den = 0.0d0
    do i = 1, n_COs
      do j = 1, n_COs
        val = 0.0d0
        do iocc = 1, nocc
          val = val + 2.0d0*MOCoeffs(iocc, i)*MOCoeffs(iocc, j)
        end do
        Den(i,j) = val
!       Den(j,i) = val
!       if (abs(val).gt.1e-12) then
!         write(78,*)  i, j, Den(i,j)
!       end if
      end do
    end do
!   write(78,*)  '####### Done'

  end subroutine GetDensity

  subroutine CalcHCore(OneEInts, hcore, OrbInd, ng, n_m, RemoveL)

    real(dp), allocatable, intent(in) :: OneEInts(:,:,:)
    real(dp), allocatable, intent(inout) :: hcore(:,:)
    integer, intent(in) :: ng
    integer, allocatable, intent(in) :: n_m(:)
    logical, intent(in) :: RemoveL
    integer, allocatable, intent(in) :: OrbInd(:,:,:)
    integer :: i, j, l, m, ind_1, ind_2, l_prime
    real(dp) :: val
    real(dp), allocatable :: delta(:,:)

    allocate(delta(ng, ng))

    delta = 0.0d0
    do i = 1, ng
      delta(i,i) = 1.0d0
    end do

    do i = 1, ng
      do j = 1, i
        if (RemoveL) then
          l_prime = min(i,j,para%l+1)
        else
          l_prime =para%l+1
        end if
        do l = 1, l_prime
!         l_val = real(l * (l - 1), idp) / (two * para%mass * grid%r(i)**2) * delta(i,j)
          !do m = 1, 2*l-1
          do m = 1, n_m(l)
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

  subroutine CalcVred(TwoERadInts, AngInts, Den, Vred, OrbInd, n_nqn, n_l, n_m)

    real(dp), allocatable, intent(in) :: TwoERadInts(:,:,:), Den(:,:)
    complex(idp), allocatable, intent(in) :: AngInts(:,:,:,:,:)
    real(dp), allocatable, intent(inout) :: Vred(:,:)
    integer, allocatable, intent(in) :: OrbInd(:,:,:)
    integer, intent(in) :: n_nqn, n_l
    integer, allocatable, intent(in) :: n_m(:)

    integer :: k1, k3, l1, l2, l3, l4, m1, m2, m3, m4, lm1, lm2, lm3, lm4, nl2
    integer :: n_m1, n_m2, n_m3, n_m4, klm_1, klm_2, klm_3, klm_4, l, error
    complex(idp) :: val
    complex(idp), allocatable :: AngEls(:,:,:,:,:)
    real(dp) :: start, finish
    integer, allocatable :: ml_interm(:)

    !nl2 = n_l**2
    nl2 = para%dim_l

    call cpu_time(start)

    allocate(AngEls(2*para%l+1, nl2, nl2, nl2, nl2), stat=error)
    call allocerror(error)

    allocate(ml_interm(n_l), stat=error)
    call allocerror(error)

    ml_interm(1) = 0
    do l = 2, n_l
      ml_interm(l) = sum(n_m(1:l-1))
    end do

    do l4 = 1, n_l
      n_m4 = 2*l4 - 1
      do m4 = 1, n_m(l4)
        !lm4 = (l4-1)**2 + m4
        lm4 = ml_interm(l4) + m4
        do l3 = 1, n_l
          n_m3 = 2*l3 - 1
          do m3 = 1, n_m(l3)
            lm3 = ml_interm(l3) + m3
            do l2 = 1, n_l
              n_m2 = 2*l2 - 1
              do m2 = 1, n_m(l2)
                lm2 = ml_interm(l2) + m2
                do l1 = 1, n_l
                  n_m1 = 2*l1 - 1
                  do m1 = 1, n_m(l1)
                    lm1 = ml_interm(l1) + m1

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

!   do lm4 = 1, nl2
!     do lm3 = 1, nl2
!       do lm2 = 1, nl2
!         do lm1 = 1, nl2
!           do l = 1, 2*para%l + 1
!              val = AngEls(l, lm1, lm2, lm3, lm4)
!              if (abs(val).gt.1e-12) then
!                do k3 = 1, n_nqn
!                  klm_3 = (k3 -1)*nl2 + lm3 
!                  klm_4 = (k3 -1)*nl2 + lm4 
!                  if (abs(Den(klm_3, klm_4)).gt.1e-12) then
!                    do k1 = 1, n_nqn
!                      klm_1 = (k1 -1)*nl2 + lm1
!                      klm_2 = (k1 -1)*nl2 + lm2 
!                      Vred(klm_1, klm_2) = Vred(klm_1, klm_2) + &
!                      &  val * Den(klm_3, klm_4) * TwoERadInts(k1,k3,l)
!                    end do
!                  end if
!                end do
!              end if
!           end do
!         end do
!       end do
!     end do
!   end do

    do l4 = 1, para%l + 1
    do m4 = 1, n_m(l4)
      lm4 = ml_interm(l4) + m4
      do l3 = 1, para%l + 1
      do m3 = 1, n_m(l3)
        lm3 = ml_interm(l3) + m3
        do l2 = 1, para%l + 1
        do m2 = 1, n_m(l2)
          lm2 = ml_interm(l2) + m2
          do l1 = 1, para%l + 1
          do m1 = 1, n_m(l1)
            lm1 = ml_interm(l1) + m1
            do l = 1, 2*para%l + 1
               val = AngEls(l, lm1, lm2, lm3, lm4)
               if (abs(val).gt.1e-12) then
                 do k3 = 1, n_nqn
                   klm_3 = OrbInd(k3,l3,m3)
                   klm_4 = OrbInd(k3,l4,m4)
                   if (klm_3.eq.0.or.klm_4.eq.0) cycle
                   if (abs(Den(klm_3, klm_4)).gt.1e-12) then
                     !write(80,'(2i4,f15.10)') klm_3, klm_4, Den(klm_3, klm_4)
                     do k1 = 1, n_nqn
                       klm_1 = OrbInd(k1,l1,m1)
                       klm_2 = OrbInd(k1,l2,m2)
                       if (klm_1.eq.0.or.klm_2.eq.0) cycle
                       Vred(klm_1, klm_2) = Vred(klm_1, klm_2) + &
                       &  val * Den(klm_3, klm_4) * TwoERadInts(k1,k3,l)
                     end do
                   end if
                 end do
               end if
            end do
          end do
          end do
        end do
        end do
      end do
      end do
    end do
    end do

!   do k1 = 1, ng
!     do l1 = 1, n_l
!       n_m1 = 2*l1 - 1
!       do m1 = 1, n_m1
!         lm1 = (l1-1)**2 + m1
!         klm_1 = OrbInd(k1,l1,m1)

!           do l2 = 1, n_l
!             n_m2 = 2*l2 - 1
!             do m2 = 1, n_m2
!               lm2 = (l2-1)**2 + m2
!               klm_2 = OrbInd(k1,l2,m2)

!               int_value = 0.0d0
!               do k3 = 1, ng
!               do l3 = 1, n_l
!                 n_m3 = 2*l3 - 1
!                 do m3 = 1, n_m3
!                   lm3 = (l3-1)**2 + m3
!                   klm_3 = OrbInd(k3,l3,m3)


!                   do l4 = 1, n_l
!                     n_m4 = 2*l4 - 1
!                     do m4 = 1, n_m4
!                       lm4 = (l4-1)**2 + m4
!                       klm_4 = OrbInd(k3,l4,m4)

!                       if (abs(Den(klm_3,klm_4)).gt.1e-12) then
!                         do l = 1, 2*para%l+1
!                           int_value = int_value + Den(klm_3,klm_4)*(AngInts(l, lm1, lm3, lm2, lm4) - &
!                           & 0.5d0*AngInts(l, lm1, lm3, lm4, lm2)) * TwoERadInts(k1,k3,l)
!                         end do
!                       end if

!                     end do
!                   end do
!  
!                 end do
!               end do
!               end do

!               Vred(klm_1, klm_2) = int_value

!               if (abs(int_value).gt.1e-12) then
!               write(79,'(2i4,f15.8)') klm_1, klm_2, real(int_value)
!               end if

!             end do
!           end do

!       end do
!     end do
!   end do

    deallocate(AngEls, ml_interm)
    call cpu_time(finish)

    if (debug.gt.6) &
    &   write(iout, *) 'The reduced V is calculated in ', (finish-start), 'seconds...'

  end subroutine CalcVred

  subroutine CalcVred_2(TwoEInts, Den, Vred, nTotOrbs)

    real(dp), allocatable, intent(in) :: TwoEInts(:,:,:,:), Den(:,:)
    real(dp), allocatable, intent(inout) :: Vred(:,:)
    integer, intent(in) ::  nTotOrbs

    integer :: klm_1, klm_2, klm_3, klm_4
    complex(idp) :: int_value
    real(dp) :: start, finish

    call cpu_time(start)

    Vred = 0.0d0
    do klm_1 = 1, nTotOrbs
      do klm_2 = 1, nTotOrbs
        int_value = 0.0d0
        do klm_3 = 1, nTotOrbs
          do klm_4 = 1, nTotOrbs

            if (abs(Den(klm_3,klm_4)).gt.1e-12) then
              int_value = int_value + Den(klm_3,klm_4)*(TwoEInts(klm_1,klm_3,klm_2, klm_4) - 0.5d0*TwoEInts(klm_1,klm_3,klm_4,klm_2))
          end if

          end do
        end do
        Vred(klm_1, klm_2) = int_value
        if (abs(int_value).gt.1e-12) then
        write(79,'(2i4,f15.8)') klm_1, klm_2, real(int_value)
        end if
      end do
    end do
    write(79,*) 'Done####'

    call cpu_time(finish)

    if (debug.gt.6) &
    & write(iout, *) 'The reduced V is calculated in ', (finish-start), 'seconds...'

  end subroutine CalcVred_2

  subroutine GetFock(hcore, V, F, n_COs)

    real(dp), allocatable, intent(in)    :: hcore(:,:)
    real(dp), allocatable, intent(in)    :: V(:,:)
    real(dp), allocatable, intent(inout) :: F(:,:)
    integer,               intent(in)    :: n_COs

    integer  :: i, j, len_1, len_2, n_l
    real(dp) :: start, finish, val

    call cpu_time(start)

    F = zero

    if (para%split_grid) then

      if (para%diagtype == 'only_inner'.or.para%diagtype == 'only_outer') &
      &  call stop_all('GetFock', 'RHF is not implemented yet for only_inner or only_outer')

      n_l = para%l + 1
      len_1 = para%m1*para%nl*n_l*n_l
      len_2 = (para%m2*para%nl-1)*n_l*n_l

      do j = 1, len_1
        do i = 1, j
          val = hcore(i,j) + V(i,j)
          F(i,j) = val
          F(j,i) = val
!         if (abs(F(i,j)).gt.1e-12) write(84,*) i, j, F(i,j)
        end do
      end do
      do j = len_1+1, len_1+len_2
        do i = len_1+1, j
          val = hcore(i,j) + V(i,j)
          F(i,j) = val
          F(j,i) = val
!         if (abs(F(i,j)).gt.1e-12) write(84,*) i, j, F(i,j)
        end do
      end do

    else
      do j = 1, n_COs
        do i = 1, j
          val = hcore(i,j) + V(i,j)
          F(i,j) = val
          F(j,i) = val
!         if (abs(F(i,j)).gt.1e-12) write(84,*) i, j, F(i,j)
        end do
      end do
    end if
!   write(84, *) 'Done####'

    call cpu_time(finish)

    if (debug.gt.6) &
    & write(iout, *) 'Fock matrix is calculated in ', (finish-start), 'seconds...'

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

    if (debug.gt.6) &
    &  write(iout, *) 'Energy is calculated in ', (finish-start), 'seconds...'

  end subroutine CalcEnergy

  subroutine DiagFock(F, MOCoeff, ng, n_m, OrbEn)

    use dvr_diag_mod, only : diag_matrix

    real(dp), allocatable, intent(inout) :: F(:,:), MOCoeff(:,:)
    real(dp), allocatable, intent(inout) :: OrbEn(:)
    integer, intent(in) :: ng
    integer, allocatable, intent(in) :: n_m(:)

    real(dp), allocatable :: F_r(:,:), En(:)
    integer, allocatable :: OrbInd(:,:,:)
    integer, allocatable :: NewOrbInd(:,:,:)
    integer :: error, i, j, indx
    integer :: ind_1, ind_2, l, m, len_1, len_2
    real(dp) :: val, start, finish


    call cpu_time(start)

!   allocate(F_r(n_COs,n_COs), stat=error)
    allocate(F_r(ng,ng), stat=error)
    allocate(En(ng), stat=error)
    allocate(NewOrbInd(ng, para%l+1, 2*para%l+1), stat=error)
    allocate(OrbInd(ng, para%l+1, 2*para%l+1), stat=error)

    indx = 0
    do i = 1, ng
      do l = 1, para%l+1
        do m = 1, n_m(l)
          indx = indx + 1
          OrbInd(i, l, m) = indx
        end do
      end do
    end do

    indx = 0
    do l = 1, para%l+1
      do m = 1, n_m(l)
        do i = 1, ng
          indx = indx + 1
          NewOrbInd(i, l, m) = indx
        end do
      end do
    end do

!   F_r = 0.0d0
!   do i = 1, ng
!     do j = 1, i
!       do l = 1, para%l+1 
!         do m = 1, 2*l-1
!           ind_1 = OrbInd(i, l, m)
!           ind_2 = OrbInd(j, l, m)
!           ind_3 = NewOrbInd(i, l, m)
!           ind_4 = NewOrbInd(j, l, m)
!           val = F(ind_1, ind_2)
!           F_r(ind_3, ind_4) = val
!           F_r(ind_4, ind_3) = val
!         end do
!       end do
!     end do
!   end do


    OrbEn = 0.0d0

    if (para%split_grid) then

      len_1 = para%m1*para%nl
      len_2 = para%m2*para%nl - 1
      if (allocated(F_r)) deallocate(F_r)
      allocate(F_r(len_1, len_1))
      if (allocated(En)) deallocate(En)
      allocate(En(len_1))

      do l = 1, para%l+1
        do m = 1, n_m(l)
          F_r = 0.0d0
          do i = 1, len_1
            do j = 1, i
              ind_1 = OrbInd(i, l, m)
              ind_2 = OrbInd(j, l, m)
              val = F(ind_1, ind_2)
              F_r(i,j) = val
              F_r(j,i) = val
            end do
          end do
          call diag_matrix(F_r, En)
          do i = 1, len_1
            ind_1 = OrbInd(i, l, m)
            do j = 1, len_1
              ind_2 = OrbInd(j, l, m)
              val = F(j,i)
              F(ind_1,ind_2) = F_r(i,j)
              F(ind_2,ind_1) = F_r(j,i)
            end do
            OrbEn(ind_1) = En(i)
          end do
        end do
      end do

      deallocate(F_r)
      allocate(F_r(len_2, len_2))
      deallocate(En)
      allocate(En(len_2))

      do l = 1, para%l+1
        do m = 1, n_m(l)
          F_r = 0.0d0
          do i = 1, len_2
            do j = 1, i
              ind_1 = OrbInd(i+len_1, l, m)
              ind_2 = OrbInd(j+len_1, l, m)
              val = F(ind_1, ind_2)
              F_r(i,j) = val
              F_r(j,i) = val
            end do
          end do
          call diag_matrix(F_r, En)
          do i = 1, len_2
            ind_1 = OrbInd(i+len_1, l, m)
            do j = 1, len_2
              ind_2 = OrbInd(j+len_1, l, m)
              val = F(j,i)
              F(ind_1,ind_2) = F_r(i,j)
              F(ind_2,ind_1) = F_r(j,i)
            end do
            OrbEn(ind_1) = En(i)
          end do
        end do
      end do

    else

      do l = 1, para%l+1
        do m = 1, n_m(l)
          F_r = 0.0d0
          do i = 1, ng
            do j = 1, i
              ind_1 = OrbInd(i, l, m)
              ind_2 = OrbInd(j, l, m)
              val = F(ind_1, ind_2)
              F_r(i,j) = val
              F_r(j,i) = val
            end do
          end do
          call diag_matrix(F_r, En)
          do i = 1, ng
            ind_1 = OrbInd(i, l, m)
            do j = 1, ng
              ind_2 = OrbInd(j, l, m)
              val = F(j,i)
              F(ind_1,ind_2) = F_r(i,j)
              F(ind_2,ind_1) = F_r(j,i)
            end do
            OrbEn(ind_1) = En(i)
          end do
        end do
      end do

    end if

!   call diag_matrix(F_r, OrbEn)

!   F = 0.0d0
!   do i =1, n_COs
!     do j = 1, ng
!       do l = 1, para%l+1 
!         do m = 1, 2*l-1
!           ind_1 = OrbInd(j, l, m)
!           ind_2 = NewOrbInd(j, l, m)
!           val = F_r(ind_2, i)
!           F(ind_1, i) = val
!         end do
!       end do
!     end do
!   end do


!   call diag_matrix(F, OrbEn)

!   do i = 1, n_COs
!     write(80,'(i4,f15.8)') i, OrbEn(i)
!   end do
!   write(80,*) 'Done####'

    MOCoeff = transpose(F)

    call cpu_time(finish)
    if (debug.gt.6) &
    & write(iout, *) 'Fock matrix is diagonalized in: ', (finish - start), 'secs...'

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
!       Del = Del + (DenOld(i,j))**2
!       val = Den(i,j) - DenOld(i,j)
!       if (abs(val).gt.1e-12) write(83,*) i, j, val
      end do
    end do

    Del = sqrt(Del/4.0d0)

  end subroutine CalcDelDen

  subroutine Calc2ePrimOrbInts(TwoEInts, OrbInd, ng, n_l, n_orb)

    use DVRData, only : two_e_rad_int, para, integrals_ang

    real(dp), allocatable, intent(inout) :: TwoEInts(:,:,:,:)
    integer, allocatable, intent(in) :: OrbInd(:,:,:)
    integer, intent(in) :: ng, n_l, n_orb
    integer  :: n_mp, l1, l2, l3, l4, m1, m2, m3, m4, k1, k2, l
    integer  :: klm_1, klm_2, klm_3, klm_4, error
    integer  :: n_m1, n_m2, n_m3, n_m4, lm1, lm2, lm3, lm4
    real(dp) :: start, finish, int_value, tol
    integer :: i, j, k, ij, kl, ijkl

    call cpu_time(start)

    tol = 1e-12
    n_mp = 2*n_l - 1

    allocate(TwoEInts(n_orb, n_orb, n_orb, n_orb), stat= error)
    call allocerror(error)

    TwoEints = zero

    do l1 = 1, n_l
      n_m1 = 2*l1 - 1
      do l2 = 1,  n_l
        n_m2 = 2*l2 - 1
        do l3 = 1,  n_l
          n_m3 = 2*l3 - 1
          do l4 = 1,  n_l
            n_m4 = 2*l4 - 1
            do m1 = 1, n_m1
              lm1 = (l1-1)**2 + m1
              do m2 = 1, n_m2
                lm2 = (l2-1)**2 + m2
                do m3 = 1, n_m3
                  lm3 = (l3-1)**2 + m3
                  do m4 = 1, n_m4
                    lm4 = (l4-1)**2 + m4
                    do k1 = 1, ng
                      do k2 = 1, ng
                        klm_1 = OrbInd(k1,l1,m1)
                        klm_2 = OrbInd(k2,l2,m2)
                        klm_3 = OrbInd(k1,l3,m3)
                        klm_4 = OrbInd(k2,l4,m4)
                      
                        int_value = 0.0d0
                        do l = 1, 2*para%l + 1
                          int_value = int_value + (integrals_ang(l, lm1, lm2, lm3, lm4)* &
                          &  two_e_rad_int(k1,k2,l))
                        end do

                        TwoEInts(klm_1, klm_2, klm_3, klm_4) = int_value

                      end do
                     end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

    ij  = 0
    ijkl = 0
    do i = 1, n_orb
      do j = 1, i
        kl = 0
        do k = 1, i
          do l = 1, k
            if (ij.ge.kl) then
              if (abs(TwoEInts(i,k,j,l)).gt.tol) &
              & write(96, 1005) TwoEInts(i,k,j,l), i, j, k, l
              ijkl = ijkl + 1
            end if
            kl = kl + 1
          end do
        end do
        ij = ij + 1
      end do
    end do

1005 format(f20.16,x,5i5)

    call cpu_time(finish)

  end subroutine Calc2ePrimOrbInts

  subroutine ContractBasis(EigenVecs, MOCoeffs, OrbEn_l, OrbEn, OrbInd, ng, n_l, n_nqn)

    real(dp), allocatable, intent(inout) :: EigenVecs(:,:,:), OrbEn_l(:,:)
    real(dp), allocatable, intent(in) :: MOCoeffs(:,:), OrbEn(:)
    integer, allocatable, intent(in) :: OrbInd(:,:,:)
    integer, intent(in) :: ng, n_l, n_nqn

    integer :: i, n, l, indx_1, indx_2
    real(dp) :: val

    EigenVecs = 0.0d0
    indx_1 = 0
    indx_2 = 0
    do l = 1, n_l
      !do n = 1, n_nqn - (l - 1)
      do n = 1, n_nqn
        indx_1 = OrbInd(n,l,1)
        do i = 1, ng
          !indx_1 = OrbInd(n+l-1,l,1)
          indx_2 = OrbInd(i,l,1)
          val = MOCoeffs(indx_1, indx_2) 
          EigenVecs(i,n,l) = val 
!         do m = 1, 2*l-1
!         end do
        end do
        OrbEn_l(n,l) = OrbEn(indx_1)
      end do
    end do

  end subroutine ContractBasis

end module RHFMod
