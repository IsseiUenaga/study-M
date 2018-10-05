
program isseiuenaga
implicit none
integer :: i
integer,parameter :: m = 2, n = 1, NN = 10000
real(8) :: dt, epsilon_a, epsilon_delta_D, epsilon_r_w, epsilon_tau_w, kappa_X, Re, Mm, gamma
real(8) :: W_s(0:NN), r_p_s, r_m_s, alpha, beta, T_theta_EM_s, T_z_EM_s, cap_delta_mode_p, cap_delta_wall_p
real(8) :: domega_theta_s(0:NN), domega_z_s(0:NN), t, cap_delta_p_s, omega, rho, m_p, m_t, I_1

open(11,file="result-m_i_4.dat",status="replace")

dt = 10.0d0/dble(NN)
I_1 = 0.8227d0
omega = 1.0d0
rho = 1.0d0
gamma = 5.0d0
epsilon_a = 2.0d0
epsilon_delta_D = 1.0d0
epsilon_r_w = 2.5d0
epsilon_tau_w = 1.0d0
kappa_X = 0.5d0
Re = 10.0d0
Mm = 0.1d0

W_s(0) = 0.001d0
domega_theta_s(0) = 0.01d0
domega_z_s(0) = 0.01d0
r_p_s = 1.0d0 +W_s(0)/2.0d0
r_m_s = 1.0d0 -W_s(0)/2.0d0
m_p = rho*(r_p_s**4 -r_m_s**4)/4.0d0
m_t = rho*(r_p_s**2 -r_m_s**2)/2.0d0
alpha = r_p_s**3 +r_m_s**3
beta = 1.0d0/(log(epsilon_a/r_m_s) +kappa_X)
T_theta_EM_s = -m**3*n*W_s(0)*((omega*epsilon_tau_w)*(r_p_s/epsilon_r_w)**(2*m))/ &
               (1.0d0 +(omega*epsilon_tau_w)**2*(1.0d0 -(r_p_s/epsilon_r_w)**(2*m))**2)
T_z_EM_s = -dble(n)*T_theta_EM_s/dble(m)
cap_delta_mode_p = 2.0d0
cap_delta_wall_p = -2.0d0*dble(m)*((omega*epsilon_tau_w)**2*(r_p_s/epsilon_r_w)**(2*m)* &
                   (1.0d0 -(r_p_s/epsilon_r_w)**(2*m)))/(1.0d0 + &
                   (omega*epsilon_tau_w)**(2*m)*(1.0d0 -(r_p_s/epsilon_r_w)**(2*m))**2)
cap_delta_p_s = cap_delta_mode_p +cap_delta_wall_p

!write(*,*) "r_p_s =", r_p_s, "r_m_s =", r_m_s
!write(*,*) "m_p =", m_p, "m_t =", m_t
!write(*,*) "alpha =", alpha, "beta =", beta
!write(*,*) "T_theta_EM_s =", T_theta_EM_s, "T_z_EM_s =", T_z_EM_s
!write(*,*) "cap_delta_p_s =", cap_delta_p_s

do i = 0,NN
   W_s(i+1) = W_s(i) +(cap_delta_p_s -gamma*W_s(i))*dt/I_1
   domega_theta_s(i+1) = domega_theta_s(i) +(-alpha*domega_theta_s(i)/(Re*epsilon_delta_D) + &
                         T_theta_EM_s/Mm**2)*dt/m_p
   domega_z_s(i+1) = domega_z_s(i) +(-beta*domega_z_s(i)/Re +T_z_EM_s/Mm**2)*dt/m_t
end do
do i = 0,NN
   t = dble(i)*dt
   write(11,*) t, W_s(i), domega_theta_s(i), domega_z_s(i)
end do

close(11)
end

