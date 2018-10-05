
program isseiuenaga
implicit none
integer,parameter :: N = 500
real(8) :: tmin, tmax, t_vec(0:N)

open(11,file="result-i_v.dat",status="replace")
open(12,file="test.dat",status="replace")

call set_parameter(tmin,tmax)
call time_interval(N,tmin,tmax,t_vec)
call calculation(N,t_vec)

close(11)
close(12)

end program isseiuenaga



subroutine set_parameter(tmin,tmax)
implicit none
real(8),intent(out) :: tmin, tmax
tmin = 0.0d0
tmax = 1.0d0
end subroutine set_parameter


subroutine time_interval(N,tmin,tmax,t_vec)
implicit none
integer :: i
integer,intent(in) :: N
real(8) :: dt, f3
real(8),intent(in) :: tmin, tmax
real(8),intent(out) :: t_vec(0:N)
dt = (tmax -tmin)/dble(N)
do i = 0,N
   t_vec(i) = f3(dble(i)*dt)
end do
end subroutine time_interval


subroutine calculation(N,t_vec)
implicit none
integer :: ite1, ite2
integer,intent(in) :: N
real(8) :: t, dt
real(8) :: Fold, Fnew, Fexact, f1, f2
real(8) :: t_vec(0:N), f_vec(0:N)

do ite1 = 0,N

   Fold = 0.0d0
   t = t_vec(ite1)
   f_vec(ite1) = f1(t)

   do ite2 = 0,ite1
      if(ite2 == 0) cycle
      dt = t_vec(ite2) -t_vec(ite2-1)
      Fnew = Fold +dt*f_vec(ite2)
      Fold = Fnew
   end do

   Fexact = f2(t)

   write(11,*) t, Fnew, Fexact

end do

end subroutine calculation



function f1(t)
real(8) :: t
real(8) :: f1
f1 = t
end function f1


function f2(t)
real(8) :: t
real(8) :: f2
f2 = t**2/2.0d0
end function f2


function f3(t)
real(8) :: t
real(8) :: f3
f3 = t**2
end function f3


