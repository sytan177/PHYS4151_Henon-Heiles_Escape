program orbit_integrate
implicit none

integer, parameter :: numofpar = 41966
integer :: neqn
parameter (neqn=4)

intrinsic dble
real*8 :: ti, tf, y_init(4*numofpar), particle_list(4*numofpar)
real*8 :: yi(neqn), yf(neqn)
real*8, parameter :: t0 = 0.d0
real*8, parameter :: tmax = 10000.d0
real*8, parameter :: dt = 1.0d-2
real*8 :: x_i, y_i, dx_i, r_s, dy_i, dis

character(len=50) :: init_file
character(len=50) :: exit_file, final_file
integer :: i, j, k, ind
external hhdrv, hhdrv_i, rk4n

open (unit = 15, file = '/Users/shuyut/Downloads/PHYS4151_Project/Initial_coords/extra_energy/hh_initial_300_0.45.txt')
    do i = 1, 4*numofpar
    read(15, *, end=99) particle_list(i)
    particle_list(i) = dble(particle_list(i))
    y_init(i) = particle_list(i)
    end do
 99 close(15)
print *, particle_list(1)
print *, particle_list(9999)
print *, particle_list(1235)


! Initial values of particles
r_s = 2.0d0
exit_file = "escape_list_0.45_10000s.txt"
final_file = "final_state_0.45_10000s.txt"

open(24, file = exit_file)
open(7, file = final_file)

i = 1
do while (i <= numofpar)
    ind = 4*(i-1)
    yi(1) = y_init(ind+1)
    yi(2) = y_init(ind+2)
    yi(3) = y_init(ind+3)
    yi(4) = y_init(ind+4)
    ti = t0
    do while (ti <= tmax)
        tf = ti+dt
        call rk4n(hhdrv_i,ti,tf,yi,yf,neqn)
        dis = yf(1)**2 + yf(2)**2
        if (dis > 100.0d0 .and. yf(2) > 1.0d0) then
            write(24,*)  i, y_init(ind+1), y_init(ind+2), y_init(ind+4), 1, tf
            print *, i
            go to 101
        else if (dis > 100.0d0 .and. yf(1) < - 1.0d0 .and. yf(2) < - 1.0d0) then
            write(24,*)  i, y_init(ind+1), y_init(ind+2), y_init(ind+4), 2, tf
            go to 101
        else if (dis > 100.0d0 .and. yf(1) > sqrt(3.0d0)/2.0d0 .and. yf(2) < - 1.0d0) then
            write(24,*)  i, y_init(ind+1), y_init(ind+2), y_init(ind+4), 3, tf
            go to 101
        end if
        ti = tf
        do j = 1, 4
            yi(j) = yf(j)
        end do
    end do
    do k = 1, neqn
    write(7,*) yf(k)
    end do
101 i = i + 1
end do

end program orbit_integrate


subroutine hhdrv_i(t, y, dy, n)
! y: x = y(1); y = y(2); dx = y(3); dy = y(4)
! dy: dx = dy(1); dy = dy(2); dx^2 = dy(3); dy^2 = dy(4)
implicit none
integer :: n
real*8 :: t, y(n), dy(n)
real*8 :: fx, fy

! Derivative subroutine for the hh potential
fx = -y(1) - 2*y(1)*y(2)
fy = -y(1)**2.0d0 - y(2) + y(2)**2
dy(1) = y(3)
dy(2) = y(4)
dy(3) = fx
dy(4) = fy

end subroutine hhdrv_i


subroutine hhdrv(t,y,dy,n)
! y: x = y(1); y = y(2); dx = y(3); dy = y(4)
! dy: dx = dy(1); dy = dy(2); dx^2 = dy(3); dy^2 = dy(4)
implicit none
integer, parameter :: numofpar = 18709
integer :: n
real*8 :: t, y(n), dy(n)
real*8 :: fx, fy
integer :: i, ind

! Derivative subroutine for the hh potential
do i = 1, numofpar
ind = 4*(i-1)
!hhpot = 0.5d0*(y(1+ind)**2.0d0+y(2+ind)**2.0d0)+y(1+ind)**2.0d0*y(2+ind)-1.0d0/3.0d0*y(2+ind)**3.0d0
fx = -y(1+ind) - 2*y(1+ind)*y(2+ind)
fy = -y(1+ind)**2.0d0 - y(2+ind) + y(2+ind)**2
dy(1+ind) = y(3+ind)
dy(2+ind) = y(4+ind)
dy(3+ind) = fx
dy(4+ind) = fy
end do
return
end subroutine hhdrv

subroutine rk4n(fcn,ti,tf,yi,yf,n)
!===========================================================
! Solution for a system of n first-order ODEs
! Method:  Runge-Kutta 4th-order
! Comment: can be easily used for n/2 second order ODEs
! Alex G. February 2010
!-----------------------------------------------------------
! call ...
! fcn(t,x,dx,n)- functions dx/dt   (supplied by a user)
! input ...
! ti    - initial time
! tf    - solution time
! yi()  - initial values
! n     - number of first order equations
! output ...
! yf()  - solutions
!===========================================================
implicit none
integer :: n
real*8 :: ti, tf
real*8 :: yi(n), yf(n)

integer :: j
real*8 :: h, t
real*8 :: y(n), dy(n)
real*8 :: k1(n),k2(n),k3(n),k4(n)

h = tf-ti
t = ti

!* evaluate k1
call fcn(t, yi, dy, n)
do j=1,n
   k1(j) = h*dy(j)
   y(j)  = yi(j) + k1(j)/2.0d0
end do

!* evaluate k2
call fcn(t+h/2.0d0,y,dy,n)
do j=1,n
   k2(j) = h*dy(j)
   y(j)  = yi(j) + k2(j)/2.d0
end do

!* evaluate k3
call fcn(t+h/2.0d0,y,dy,n)
do j=1,n
   k3(j) = h*dy(j)
   y(j)  = yi(j) + k3(j)
end do

!* evaluate k4 and the result
call fcn(t+h,y,dy,n)
do j=1,n
   k4(j) = h*dy(j)
   yf(j) = yi(j) + k1(j)/6.0d0+k2(j)/3.0d0+k3(j)/3.0d0+k4(j)/6.0d0
end do

end subroutine rk4n

