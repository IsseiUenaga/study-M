
program isseiuenaga
implicit none
integer :: ite1, ite2, N1, N2
real(8) ::  tmin, tmax, dt

open(11,file="result-i_v.dat",status="replace")
open(12,file="test.dat",status="replace")

call set_parameter(N1,N2,tmin,tmax)
call time_interval(N1,tmin,tmax,dt)
call calculation(N1,N2,tmin,dt)

close(11)
close(12)

end program isseiuenaga



subroutine set_parameter(N1,N2,tmin,tmax)
implicit none
integer,intent(out) :: N1, N2
real(8),intent(out) :: tmin, tmax
N1 = 100
N2 = 100
tmin = 0.0d0
tmax = 1.0d0
end subroutine set_parameter


subroutine time_interval(N,tmin,tmax,dt)
implicit none
integer,intent(in) :: N
real(8),intent(in) :: tmin, tmax
real(8),intent(out) :: dt

dt = (tmax -tmin)/dble(N)

end subroutine time_interval


subroutine calculation(N1,N2,tmin,dt)
implicit none
integer :: ite1, ite2
integer,intent(in) :: N1, N2
real(8),intent(in) :: dt, tmin
real(8) :: dtp, t, tmax, told
real(8) :: Fold, Fnew, f

do ite1 = 0,N1

   Fold = 0.0d0
   t = dble(ite1)*dt
   tmax = t

   call time_interval(N2,tmin,tmax,dtp)
   call calculation_step_doubling(dtp,t,told,tmax,Fold)

   write(11,*) told, Fold

end do

end subroutine calculation


subroutine Runge_Kutta_method_of_4th_order(t,dt,Fold,Fnew)
implicit none
real(8) :: RK(1:4), f
real(8),intent(in) :: t, dt, Fold
real(8),intent(out) :: Fnew

RK(1) = f(t)
RK(2) = f(t+dt/2.0d0)
RK(3) = f(t+dt/2.0d0)
RK(4) = f(t+dt)

Fnew = Fold +dt*(RK(1) +2.0d0*RK(2) +2.0d0*RK(3) +RK(4))/6.0d0

end subroutine Runge_Kutta_method_of_4th_order


subroutine step_doubling(dt1,t,Fold,Fnew1,Fnew2)
implicit none
real(8) :: dt1, dt2, t
real(8) :: Fold, Fnew
real(8),intent(out) :: Fnew1, Fnew2

call Runge_Kutta_method_of_4th_order(t,dt1,Fold,Fnew1)

dt2 = dt1/2.0d0
call Runge_Kutta_method_of_4th_order(t,dt2,Fold,Fnew)
t = t+dt2
Fold = Fnew
call Runge_Kutta_method_of_4th_order(t,dt2,Fold,Fnew2)

end subroutine step_doubling


subroutine calculation_step_doubling(dt,t,told,tmax,Fold)
implicit none
integer,parameter :: N = 10**9
integer :: ite
real(8) :: delta, delta0, delta1, error
real(8) :: dt1, t, told
real(8) :: Fold, Fnew1, Fnew2
real(8),intent(in) :: dt, tmax

delta0 = 1.0d-5

told = t

do ite = 1,N
   t = told
   if(dt <= 1.0d-50) exit
   dt1 = dt
   write(12,*) "iteration :", ite, "dt =", dt
   do
     call step_doubling(dt1,t,Fold,Fnew1,Fnew2)
     error = abs(Fnew2 -Fnew1)
     delta1 = error
     if(delta1/dt1 > delta0) then
        delta = (delta0/(2.0d0*(delta1/dt1)))**(1.0d0/4.0d0)
        if(delta <= 0.1d0) then
           dt1 = 0.1d0*dt1
        elseif(delta >= 4.0d0) then
               dt1 = 4.0d0*dt1
        elseif(0.1d0 < delta .and. delta < 4.0d0) then
               dt1 = delta*dt1
        end if
        if(dt1 > dt) then
           dt1 = dt
        end if
     elseif(delta1/dt1 < delta0) then
            exit
     end if
     write(12,*) "dt1 =", dt1
   end do

   if(t >= tmax) exit
   told = t +dt1
   Fold = Fnew1

end do

end subroutine calculation_step_doubling



function f(t)
real(8) :: t
real(8) :: f
f = t
end function f

