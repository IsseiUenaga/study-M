
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
real(8) :: Fold, Fnew, Fexact, f1, f2

do ite1 = 0,N1

   Fold = 0.0d0
   t = dble(ite1)*dt
   tmax = t

   call time_interval(N2,tmin,tmax,dtp)
   do ite2 = 0,N2
      t = dble(ite2)*dtp
      Fnew = Fold +dtp*f1(t)
      Fold = Fnew
   end do

   Fexact = f2(t)

   told = tmax
   write(11,*) told, Fold, Fexact

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
