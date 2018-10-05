
program isseiuenaga
implicit none
integer :: i
integer,parameter :: N = 100
real(8) :: xs, alpha, eta, mu_0, f_0, pi, dt, f
real(8) :: f_p_p(0:N), f_p_m(0:N), t(0:N), w(0:N), delta_p(0:N)

open(11,file="result-i_w_1.dat",status="replace")

dt = 5.0d0/dble(N)
pi = 4.0d0*atan(1.0d0)
rs = 1.0d0
alpha = 0.1d0
eta = 7.46d0*10d-5
mu_0 = 4.0d0*pi*10d-7
f_0 = 1.0d0

w(0) = 0.5d0
t(0) = 0.0d0
f = f_0*exp(-alpha*xs)

do i = 0,N-1
   t(i+1) = dble(i+1)*dt
   f_p_m(i) = -alpha*f_0*exp(alpha*w(i))
   f_p_p(i) = -alpha*f_0*exp(-alpha*w(i))
   delta_p(i) = (f_p_p(i) -f_p_m(i))/f
   w(i+1) = w(i) +eta*delta_p(i)/(2.0d0*mu_0)*dt
end do

do i = 0,N
   write(11,*) t(i), w(i)
end do

close(11)

end
