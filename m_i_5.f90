
program isseiuenaga
implicit none
integer :: i
integer,parameter :: NN = 10000, m = 2
real(8) :: dt, t, W(0:NN), domega_theta(0:NN), domega_z(0:NN)
real(8) :: r_p(0:NN), r_m(0:NN), cap_delta_p(0:NN)
real(8) :: T_theta_EM(0:NN), T_z_EM(0:NN)
real(8) :: T_theta_VS(0:NN), T_z_VS(0:NN)

open(11,file="result-m_i_5.dat",status="replace")

dt = 10.d0/dble(NN)
W(0) = 0.001d0
domega_theta(0) = 0.01d0
domega_z(0) = 0.01d0

do i = 0,NN
   r_p(i) = 1.0d0 +W(i)/2.0d0
   r_m(i) = 1.0d0 -W(i)/2.0d0
   T_theta_EM(i) = -r_p(i)**(2*m)*W(i)**4/(1.0d0 +(1.0d0 -r_p(i)**(2*m))**2)
   T_z_EM(i) = -T_theta_EM(i)
   T_theta_VS(i) = -domega_theta(i)*(r_p(i)**3 +r_m(i)**3)
   T_theta_VS(i) = -domega_z(i)/(log(1.0d0/r_m(i)) +1.0d0)
   cap_delta_p(i) = (1.0d0 -5.0d0*W(i)) -r_p(i)**(2*m)*(1-r_p(i)**(2*m))/(1.0d0 +(1.0d0 -r_p(i)**(2*m))**2)
   W(i+1) = W(i) +cap_delta_p(i)*dt
   domega_theta(i+1) = domega_theta(i) +(T_theta_VS(i) +T_theta_EM(i))*dt
   domega_z(i+1) = domega_z(i) +(T_z_VS(i) +T_z_EM(i))*dt
end do
do i = 0,NN
   t = dble(i)*dt
   write(11,*) t,W(i),domega_theta(i),domega_z(i)
end do

close(11)

end
