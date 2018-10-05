
program isseiuenaga
implicit none
integer :: i
integer,parameter :: m = 2, n = 1, NN = 5000000
real(8) :: dt, t, I_1, y
real(8) :: Re, Mm, RK(1:4,1:3)
real(8) :: rho, kappa_X, f_1, omega
real(8) :: cap_delta_mode_p, r_m_s, r_p_s
real(8) :: epsilon_a, epsilon_r_w, epsilon_delta_D, epsilon_tau_w, epsilon_R0
real(8) :: W_s, domega_theta_s, domega_z_s
real(8) :: W_s_vec(0:NN), domega_theta_s_vec(0:NN), domega_z_s_vec(0:NN)

open(11,file="result-m_i_6.dat",status="replace")
open(12,file="result1-m_i_6.dat",status="replace")
open(13,file="result2-m_i_6.dat",status="replace")
open(14,file="result3-m_i_6.dat",status="replace")
open(15,file="result4-m_i_6.dat",status="replace")

dt = 1000.0d0/dble(NN)
I_1 = 0.8227d0
!plasma density
rho = 1.0d0
!angular frequency of rotation of the tearing mode
epsilon_a = 2.0d0
epsilon_r_w = 2.5d0
epsilon_R0 = 20.0d0
epsilon_delta_D = 0.001d0
epsilon_tau_w = 1.0d0
!the ratio of momentum diffusion velocity from the bulk fulid into the boundary layer to the wall/limeter
kappa_X = 0.5d0
!Reynolds number
Re = 1000.0d0
!magnetic Mach numebr
Mm = 0.001d0

W_s_vec(0) = 0.001d0
W_s = W_s_vec(0)
domega_theta_s_vec(0) = 0.01d0
domega_theta_s = domega_theta_s_vec(0)
domega_z_s_vec(0) = 0.01d0
domega_z_s = domega_z_s_vec(0)



do i = 1,NN
call Runge_Kutta_method_of_4th_order(W_s,domega_theta_s,domega_z_s,rho,kappa_X,Re,Mm, &
                                     epsilon_delta_D,epsilon_tau_w,epsilon_r_w,epsilon_a,epsilon_R0,RK,I_1,dt,m,n)
   W_s_vec(i) = W_s +(RK(1,1) +2.0d0*RK(2,1) +2.0d0*RK(3,1) +RK(4,1))*dt/6.0d0
   domega_theta_s_vec(i) = domega_theta_s +(RK(1,2) +2.0d0*RK(2,2) +2.0d0*RK(3,2) +RK(4,2))*dt/6.0d0
   domega_z_s_vec(i) = domega_z_s +(RK(1,3) +2.0d0*RK(2,3) +2.0d0*RK(3,3) +RK(4,3))*dt/6.0d0
   W_s = W_s_vec(i)
   domega_theta_s = domega_theta_s_vec(i)
   domega_z_s = domega_z_s_vec(i)
   omega = m*domega_theta_s -n*domega_z_s
   t = dble(i)*dt
   write(15,*) t, domega_theta_s, domega_z_s, omega
end do

do i = 0,NN
   t = dble(i)*dt
   write(11,*) t, W_s_vec(i), domega_theta_s_vec(i), domega_z_s_vec(i)
   if(mod(i, 1000) == 0) then
      write(*,*) t, W_s_vec(i), W_vec(i), domega_theta_s_vec(i), domega_z_s_vec(i)
   end if
end do

close(11)
close(12)
close(13)
close(14)
close(15)

end program isseiuenaga


subroutine Runge_Kutta_method_of_4th_order(W_s,domega_theta_s,domega_z_s,rho,kappa_X,Re,Mm, &
                                           epsilon_delta_D,epsilon_tau_w,epsilon_r_w,epsilon_a,epsilon_R0,RK,I_1,dt,m,n)
implicit none
integer :: m, n
real(8) :: RK(1:4,1:3), dt
real(8) :: W_s, domega_theta_s, domega_z_s
real(8) :: rho, kappa_X
real(8) :: Re, Mm, I_1, cap_delta_mode_p
real(8) :: epsilon_delta_D, epsilon_tau_w, epsilon_r_w, epsilon_a, epsilon_R0
RK(1,1) = f_1(domega_theta_s,domega_z_s,W_s, &
              epsilon_tau_w,epsilon_r_w,I_1,m,n)
RK(1,2) = f_2(domega_theta_s,domega_z_s,W_s, &
              rho,epsilon_delta_D,epsilon_tau_w,epsilon_r_w,Re,Mm,m,n)
RK(1,3) = f_3(domega_theta_s,domega_z_s,W_s, &
              rho,epsilon_a,epsilon_tau_w,epsilon_r_w,epsilon_R0,kappa_X,Re,Mm,m,n)

RK(2,1) = f_1(domega_theta_s+dt*RK(1,2)/2.0d0,domega_z_s+dt*RK(1,3)/2.0d0,W_s+dt*RK(1,1)/2.0d0, &
              epsilon_tau_w,epsilon_r_w,I_1,m,n)
RK(2,2) = f_2(domega_theta_s+dt*RK(1,2)/2.0d0,domega_z_s+dt*RK(1,3)/2.0d0,W_s+dt*RK(1,1)/2.0d0, &
              rho,epsilon_delta_D,epsilon_tau_w,epsilon_r_w,Re,Mm,m,n)
RK(2,3) = f_3(domega_theta_s+dt*RK(1,2)/2.0d0,domega_z_s+dt*RK(1,3)/2.0d0,W_s+dt*RK(1,1)/2.0d0, &
              rho,epsilon_a,epsilon_tau_w,epsilon_r_w,epsilon_R0,kappa_X,Re,Mm,m,n)

RK(3,1) = f_1(domega_theta_s+dt*RK(2,2)/2.0d0,domega_z_s+dt*RK(2,3)/2.0d0,W_s+dt*RK(2,1)/2.0d0, &
              epsilon_tau_w,epsilon_r_w,I_1,m,n)
RK(3,2) = f_2(domega_theta_s+dt*RK(2,2)/2.0d0,domega_z_s+dt*RK(2,3)/2.0d0,W_s+dt*RK(2,1)/2.0d0, &
              rho,epsilon_delta_D,epsilon_tau_w,epsilon_r_w,Re,Mm,m,n)
RK(3,3) = f_3(domega_theta_s+dt*RK(2,2)/2.0d0,domega_z_s+dt*RK(2,3)/2.0d0,W_s+dt*RK(2,1)/2.0d0, &
              rho,epsilon_a,epsilon_tau_w,epsilon_r_w,epsilon_R0,kappa_X,Re,Mm,m,n)

RK(4,1) = f_1(domega_theta_s+dt*RK(3,2),domega_z_s+dt*RK(3,3),W_s+dt*RK(3,1), &
              epsilon_tau_w,epsilon_r_w,I_1,m,n)
RK(4,2) = f_2(domega_theta_s+dt*RK(3,2),domega_z_s+dt*RK(3,3),W_s+dt*RK(3,1), &
              rho,epsilon_delta_D,epsilon_tau_w,epsilon_r_w,Re,Mm,m,n)
RK(4,3) = f_3(domega_theta_s+dt*RK(3,2),domega_z_s+dt*RK(3,3),W_s+dt*RK(3,1), &
              rho,epsilon_a,epsilon_tau_w,epsilon_r_w,epsilon_R0,kappa_X,Re,Mm,m,n)

write(12,*) RK(1,1), RK(2,1), RK(3,1), RK(4,1)
write(13,*) RK(1,2), RK(2,2), RK(3,2), RK(4,2)
write(14,*) RK(1,3), RK(2,3), RK(3,3), RK(4,3)

contains

function f_1(domega_theta_s,domega_z_s,W_s,epsilon_tau_w,epsilon_r_w,I_1,m,n)
integer :: m, n
real(8) :: omega, epsilon_tau_w, epsilon_r_w
real(8) :: f_1, r_p_s, W0, cap_delta_mode_0_p
real(8) :: cap_delta_mode_p, cap_delta_wall_p, I_1
real(8) :: domega_theta_s, domega_z_s, W_s
r_p_s = 1.0d0 +W_s/2.0d0
W0 = 0.2d0
cap_delta_mode_0_p = 2.0d0
omega = m*domega_theta_s -n*domega_z_s
cap_delta_wall_p = -2.0d0*dble(m)*(omega*epsilon_tau_w)**2*(r_p_s/epsilon_r_w)**(2*m)* &
                   (1.0d0 -(r_p_s/epsilon_r_w)**(2*m))/(1.0d0 + &
                   (omega*epsilon_tau_w)**(2*m)*(1.0d0 -(r_p_s/epsilon_r_w)**(2*m))**2)
cap_delta_mode_p = cap_delta_mode_0_p*(1.0d0 -W_s/W0)
f_1 = (cap_delta_mode_p +cap_delta_wall_p)/I_1
end function f_1

function f_2(domega_theta_s,domega_z_s,W_s,rho,epsilon_delta_D,epsilon_tau_w,epsilon_r_w,Re,Mm,m,n)
integer :: m, n
real(8) :: domega_theta_s, domega_z_s, W_s
real(8) :: Re, Mm, epsilon_delta_D, omega, rho
real(8) :: epsilon_tau_w, epsilon_r_w, f_2
real(8) :: T_theta_VS_s, T_theta_EM_s, m_p
real(8) :: r_p_s, r_m_s
r_m_s = 1.0d0 -W_s/2.0d0
r_p_s = 1.0d0 +W_s/2.0d0
omega = m*domega_theta_s -n*domega_z_s
T_theta_VS_s = -domega_theta_s*(r_m_s**3 +r_p_s**3)/(Re*epsilon_delta_D)
T_theta_EM_s = -m**2*((omega*epsilon_tau_w)*(r_p_s/epsilon_r_w)**(2*m))/ &
              (1.0d0 +(omega*epsilon_tau_w)**2*(1 -(r_p_s/epsilon_r_w)**(2*m))**2) &
               *W_s**4/Mm**2
m_p = rho*(r_p_s**4 -r_m_s**4)/4.0d0
f_2 = (T_theta_VS_s +T_theta_EM_s)/m_p
end function f_2

function f_3(domega_theta_s,domega_z_s,W_s,rho,epsilon_a,epsilon_tau_w, &
             epsilon_r_w,epsilon_R0,kappa_X,Re,Mm,m,n)
integer :: m, n
real(8) :: domega_theta_s, domega_z_s, W_s
real(8) :: epsilon_a, epsilon_tau_w, epsilon_r_w, epsilon_R0
real(8) :: omega, rho, kappa_X, Re, Mm, f_3
real(8) :: T_z_VS_s, T_z_EM_s, m_t
real(8) :: r_p_s, r_m_s
r_m_s = 1.0d0 -W_s/2.0d0
r_p_s = 1.0d0 +W_s/2.0d0
omega = m*domega_theta_s -n*domega_z_s
T_z_VS_s = -domega_z_s/((log(epsilon_a/r_m_s) +kappa_X)*re)
T_z_EM_s = m*n*((omega*epsilon_tau_w)*(r_p_s/epsilon_r_w)**(2*m))/ &
          (1.0d0 +(omega*epsilon_tau_w)**2*(1 -(r_p_s/epsilon_r_w)**(2*m))**2)* &
           W_s**4/(Mm**2*epsilon_R0**2)
m_t = rho*(r_p_s**2 -r_m_s**2)/2.0d0
f_3 = (T_z_VS_s +T_z_EM_s)/m_t
end function f_3

end subroutine Runge_Kutta_method_of_4th_order
