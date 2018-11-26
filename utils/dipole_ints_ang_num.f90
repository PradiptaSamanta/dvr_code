   ! The angular part of one-electron integrals for dipole operators which involves only the spherical harmonics are calculated here numerically using Cartesian grid points. 
   ! The grid points (x,y,z) are chosen from a total arrays of 3D grid points for which x^2 + y^2 + z^2 = r^2, where r is a fixed radius of a sphere around the nuclues.
   ! The values of spherical harmonics for different l and m are calculated and stored in arrays. Products of these spherical harmonics are then used to calculate the integrals
   program NumericalInts
   implicit none
   integer i,j,k,l,m,ngrid,n,mm,ic,igrid, n_l, nao, orb, j_val
   integer :: m_val, l_val
   double precision check,r,r12,rm,rn,pi,check6d,rval
   double precision, allocatable ::  wt(:),x(:),y(:),z(:), coords(:,:)
   double precision ::  wt_val, x_val, y_val, z_val, val_1, val_2, start, finish, fac_1,fac_2
   double precision, allocatable :: angwt(:),angx(:),angy(:),angz(:)
   double precision, allocatable :: coords2(:,:)
   complex*16, allocatable :: zpsi(:,:),zpsi_alt(:,:),zoverlap(:,:)
   complex*16, allocatable :: dipx(:,:), dipy(:,:), dipz(:,:)
   complex*16 :: zi,zcheck1,zcheck2,zcheck, zval, zval_1
   complex*16, allocatable :: zall(:,:,:,:,:)
   complex*16 :: transf_l_0, transf_l_1(3,3), transf_l_2(5,6)
   complex*16  :: get_val
   character*80 grid
   character*32 :: file_name_syntx 
   character*8 :: int2str
   logical :: real_func
   

   ! First define some common constants

   zi = dcmplx(0.0d0,1.0d0)
   pi = 3.141592653589d0
    
   ! Integrals can be obtained for spherical harmonics in both of their complex and real forms.
   ! The real spherical harmonics are obtained by combining the original descriptions of spherical
   ! harmonics which has the complex form
   real_func = .true.


   ! The common initial part of the names of files where the final integrals are written out
   !file_name_syntx = 'angular_element'

   ! Total number of 'l' quantum number to be used to define the spherical harmonics
   n_l = 2
   ! Total number of spherical harmonics
   nao = n_l**2
   write(*,*) 'nao:', nao

   allocate(zall(nao,nao,nao,nao,n_l))

   ! Sph. Harm. for l=0 does not depend on x, y and z
   transf_l_0 = 1.d0/2.d0/sqrt(pi)

   ! Initialising the transformation matrices for l=1,2
   transf_l_1 = dcmplx(0.0d0,0.0d0)

   transf_l_2 = dcmplx(0.0d0,0.0d0)

   ! transf_l_1 is contracted with the column vector (x,z,y) to produce the 3 spher. harm. for l=1
   ! transf_l_2 is contracted with the column vector (xy,yz,zx,x^2,y^2,z^2) to produce the 5 spher. harm. for l=2
   if (real_func) then
     val_1 = 1.0d0/2.d0*sqrt(3.d0/pi)
     ! Elements of transf_l_1 to obtain real sph. harm. 
     transf_l_1(3,1) = val_1
     transf_l_1(2,2) = val_1
     transf_l_1(1,3) = val_1

     ! Elements of transf_l_2 to obtain real sph. harm. 
     val_1 = 1.0d0/2.d0*sqrt(15.d0/pi)
     val_2 = 1.0d0/4.d0*sqrt(5.d0/pi)
     transf_l_2(1,1) = val_1
     transf_l_2(2,2) = val_1
     transf_l_2(3,4) = -1.0d0*val_2
     transf_l_2(3,5) = -1.0d0*val_2
     transf_l_2(3,6) =  2.0d0*val_2
     transf_l_2(4,3) = val_1
     transf_l_2(5,4) = 0.5d0*val_1
     transf_l_2(5,5) = -0.5d0*val_1
   else 
     ! Elements of transf_l_1 to obtain complex sph. harm. 
     transf_l_1(1,1) = 1.0d0/2.d0*sqrt(3.d0/(2.d0*pi))
     transf_l_1(1,3) = -1.0d0*zi/2.d0*sqrt(3.d0/(2.d0*pi))
     transf_l_1(2,2) = 1.0d0/2.d0*sqrt(3.d0/(pi))
     transf_l_1(3,1) = -1.0d0/2.d0*sqrt(3.d0/(2.d0*pi))
     transf_l_1(3,3) = -1.0d0*zi/2.d0*sqrt(3.d0/(2.d0*pi))

     ! Elements of transf_l_2 to obtain complex sph. harm. 
     val_1 = 1.0d0/4.d0*sqrt(15.d0/2.0d0/pi)
     transf_l_2(1,1) = dconjg(zi)*2.0d0*val_1
     transf_l_2(1,4) = val_1
     transf_l_2(1,5) = -1.0d0*val_1
     transf_l_2(5,1) = zi*2.0d0*val_1
     transf_l_2(5,4) = val_1
     transf_l_2(5,5) = -1.0d0*val_1
     zval_1 = 1.0d0/2.d0*sqrt(15.d0/2.0d0/pi)
     transf_l_2(2,2) = dconjg(zi)*zval_1
     transf_l_2(2,3) = -1.0d0*zval_1
     zval_1 = 1.0d0/2.d0*sqrt(15.d0/2.0d0/pi)
     transf_l_2(4,2) = dconjg(zi)*zval_1
     transf_l_2(4,3) =  1.0d0*zval_1
     val_1 = 1.0d0/4.d0*sqrt(5.d0/pi)
     transf_l_2(3,4) = -1.0d0*val_1
     transf_l_2(3,5) = -1.0d0*val_1
     transf_l_2(3,6) = 2.0d0*val_1
   end if

   zall = dcmplx(0.d0,0.d0)
   ! The file to read all the 3d grid points. This file is generated from a calculation using the
   ! quantum chemistry package 'PySCF'. A sample input to generate this file can be found in the current directory
   write(grid,'(a,i0)') '3dgrid'
   open(81,file=trim(grid))
   ! grids will be chosen for a sphere with the radius rval
   rval = 4.9114443141353972d0

   ! read the total number of 3d grid points
   read(81,*) ngrid
   
   allocate(angwt(ngrid), angx(ngrid), angy(ngrid), angz(ngrid))

   angwt = 0.0d0
   angx = 0.0d0
   angy = 0.0d0
   angz = 0.0d0
   
   ic = 0
   do m = 1,ngrid
     read(81,*) x_val,y_val,z_val,wt_val
     r = dsqrt(x_val**2+y_val**2+z_val**2)
!    if(dabs(r-1.391916542938d0/100).lt.1.d-8) then
     if(dabs(r-4.9114443141353972).lt.1.d-5) then
       ic = ic + 1
!      write(6,*) m,r,wt(m)
       angwt(ic) = wt_val
       angx(ic)  = x_val
       angy(ic)  = y_val
       angz(ic)  = z_val
     endif
   enddo
   
   write(6,*) ic
   
   ngrid = ic
   
   ! allocate the matrices which store weights, x, y and z coordinates
   allocate(wt(ngrid), x(ngrid), y(ngrid), z(ngrid), coords(3,ngrid), coords2(6,ngrid))
   
   wt(1:ic) = angwt(1:ic)/rval**2
   x(1:ic) = angx(1:ic)/(rval)
   y(1:ic) = angy(1:ic)/(rval)
   z(1:ic) = angz(1:ic)/(rval)
   ! coords contains the grid points for the column vector (x,z,y)
   coords(1,1:ic) = angx(1:ic)/(rval)
   coords(2,1:ic) = angz(1:ic)/(rval)
   coords(3,1:ic) = angy(1:ic)/(rval)
   
   ! coords2 contains the grid points for the column vector (xy,yz,zx,x^2,y^2,z^2) 
   do i = 1, ngrid
     coords2(1,i) = x(i)*y(i)
     coords2(2,i) = y(i)*z(i)
     coords2(3,i) = z(i)*x(i)
     coords2(4,i) = x(i)*x(i)
     coords2(5,i) = y(i)*y(i)
     coords2(6,i) = z(i)*z(i)
   end do

! 3d check

   ! Checking first if all the spherical harmonics are orthonormal or not
   check  = 0.d0
   zcheck1 = 0.d0
   zcheck2 = 0.d0
   
   ! zpsi  contains the value of sph. harmonics in each of the grid points
   ! zoverlap contains the value of overlap between two sph. harmonics
   allocate(zpsi(ngrid,nao), zoverlap(nao,nao))
   ! allocate matrices those store the dipole integrals
   allocate(dipx(nao,nao), dipy(nao,nao), dipz(nao,nao))

   if (real_func) allocate(zpsi_alt(ngrid,nao))

   zoverlap = 0.d0
   
   zpsi = dcmplx(0.d0, 0.d0)
   ! Get zpsi for l=0
   zpsi(:,1) = 1.d0/2.d0/sqrt(pi)

   do m = 1,ngrid

     i = 2
     ! Get zpsi for l=1
     do j = 1, 2*i-1
       orb = (i-1)**2 + j
       do k = 1,3
         zpsi(m,orb) = zpsi(m,orb) + transf_l_1(j,k)*coords(k,m)
       end do
     end do

     ! Get zpsi for l=2
     if (n_l.gt.2) then
       i = 3
       do j = 1, 2*i-1
         orb = (i-1)**2 + j
         do k = 1,6
           zpsi(m,orb) = zpsi(m,orb) + transf_l_2(j,k)*coords2(k,m)
         end do
       end do
     end if

     ! Calculating the overlap between different harmonics
     do i = 1,nao
      do j = 1,nao
        zoverlap(i,j) = zoverlap(i,j) + wt(m)*dconjg(zpsi(m,i))*zpsi(m,j)
      enddo
     enddo
   enddo
   
   !write(*,*) 'norm:', zoverlap(1,1)
   ! divide the weight to finally obtain harmonics which normalize to one
   wt = wt/(zoverlap(1,1))
   
   zoverlap = 0.d0
   do m = 1,ngrid
     do i = 1,nao
      do j = 1,nao
        zoverlap(i,j) = zoverlap(i,j) + wt(m)*dconjg(zpsi(m,i))*zpsi(m,j)
      enddo
     enddo
   enddo
   
   
   do i = 1,nao
    do j = 1,nao
      if (abs(zoverlap(i,j)).gt.1e-12) write(6,*) i,j,zoverlap(i,j)
    enddo
   enddo
!  write(6,*) 'Check3d',check,pi,ngrid,zcheck2,real(zcheck1)

! 6d check

   zcheck  = 0.d0
   zcheck1 = 0.d0
!$omp parallel do reduction(+:check6d)  private(m,n,rm,rn,r12) &
!$omp  shared (x,y,z,wt,ngrid) & 
!$omp  default (none) 

   call cpu_time(start)

   ! Calculating the dipole integrals here where the dipole field is also represented using the spherical harmonics
   val_1 = 1.0d0/2.d0*sqrt(3.d0/pi)
   dipz = 0.0d0
   dipy = 0.0d0
   dipx = 0.0d0
   fac_1 = sqrt(4.0d0*pi/3.0d0)
   fac_2 = sqrt(2.0d0*pi/3.0d0)
!  fac_1 = 1.0d0
!  fac_2 = 1.0d0
   do m = 1, ngrid
     r = sqrt(x(m)**2 + y(m)**2 + z(m)**2)
     do i =  1, nao
       do j =  1, nao
          dipz(i,j) = dipz(i,j) + fac_1*wt(m)*dconjg(zpsi(m,i))*r*zpsi(m,3)*zpsi(m,j)
          dipx(i,j) = dipx(i,j) + fac_1*wt(m)*dconjg(zpsi(m,i))*r*zpsi(m,4)*zpsi(m,j)
          dipy(i,j) = dipy(i,j) + fac_1*wt(m)*dconjg(zpsi(m,i))*r*zpsi(m,2)*zpsi(m,j)
       end do
     end do
   end do

   write(*,*) 'DIPX:'
   do i = 1,nao
     do j = 1, nao
       if (abs(dipx(i,j)).gt.1e-12) write(6,*) i,j,dipx(i,j)
     enddo
   enddo

   write(*,*) 'DIPY:'
   do i = 1,nao
     do j = 1, nao
       if (abs(dipy(i,j)).gt.1e-12) write(6,*) i,j,dipy(i,j)
     enddo
   enddo

   write(*,*) 'DIPZ:'
   do i = 1,nao
     do j = 1, nao
       if (abs(dipz(i,j)).gt.1e-12) write(6,*) i,j,dipz(i,j)
     enddo
   enddo

   call cpu_time(finish)
!$OMP END PARALLEL DO
   write(6,*) 'Integrals are calculated in', (finish-start), 'sec'
  
   deallocate(wt,x,y,z,angwt,angx,angy,angz,zpsi,zoverlap, coords)
   
   end program NumericalInts

   pure function get_val(x,y,z,l,m) result(zval)

   implicit none

   double precision, intent(in) :: x,y,z
   integer, intent(in) :: l, m
   complex*16 zval, zi
   double precision :: pi

   pi = 3.141592653589d0

   zi = dcmplx(0.0,1.0)

   zval = dcmplx(0.0d0,0.0d0)

   if (l.eq.1) then
     zval = 1.d0/2.d0/sqrt(pi)
   elseif (l.eq.2) then
     if (m.eq.-1) then
       zval = (x - zi*y)/2.d0*sqrt(3/(2*pi))
     elseif(m.eq.0) then
       zval = z/2.d0*sqrt(3/(pi))
     else
       zval = -(x+zi*y)/2.d0*sqrt(3/(2*pi))
     end if
   end if
     
   end function get_val

!      do i = 2, n_l
!        do j = 1, 2*i-1
!          j_val = -1*i + j 
!          orb = (i-1)**2 + j 
!          zpsi(m,orb) =  get_val(x(m),y(m),z(m),i,j_val)
!        end do
!      end do

!      zpsi(m,1) = 1.d0/2.d0/sqrt(pi)
!      zpsi(m,2) = (x(m) - zi*y(m))/2.d0*sqrt(3/(2*pi))
!      zpsi(m,3) = z(m)/2.d0*sqrt(3/(pi))
!      zpsi(m,4) =  -(x(m)+zi*y(m))/2.d0*sqrt(3/(2*pi))

  character(len=8) function int2str(i, format)

    integer,                    intent(in) :: i
    character(len=*), optional, intent(in) :: format

    if (present(format)) then
      write(int2str, format) i
    else
      write(int2str, '(I7)') i
    end if
    int2str = adjustl(int2str)

  end function int2str
