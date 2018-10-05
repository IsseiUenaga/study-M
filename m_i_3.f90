
program isseiuenaga
implicit none
integer :: i
integer,parameter :: m = 2, n = 1, NN = 10
real(8) :: dt, pi, I_1
real(8) :: rho, omega, t
real(8) :: kappa_X, Mm, Re
real(8) :: epsilon_tau_w, epsilon_tau_X
real(8) :: epsilon_delta_D, epsilon_delta_X, epsilon_a, epsilon_r_w
real(8) :: cap_delta_mode_p, cap_delta_wall_p(0:NN), cap_delta_p_s(0:NN)
real(8) :: r_p_s(0:NN), r_m_s(0:NN), m_p(0:NN), m_t(0:NN)
real(8) :: T_theta_EM_s(0:NN), T_z_EM_s(0:NN)
real(8) :: T_theta_VS_s(0:NN), T_z_VS_s(0:NN)
real(8) :: W_s(0:NN), domega_theta_s(0:NN), domega_z_s(0:NN)

open(11,file="result-m_i_3.dat",status="replace")

dt = 1.0d0/dble(NN)
pi = 4.0d0*atan(1.0d0)
I_1 = 0.8227d0
!plasma density
rho = 1.0d0
!angular frequency of rotation of the tearing mode
omega = 1.0d0
!standard stability index for a free boundary tearing mode
cap_delta_mode_p = 1.0d0
epsilon_a = 2.0d0
epsilon_r_w = 1.5d0
epsilon_delta_D = 10d-3
epsilon_delta_X = 10d-3
epsilon_tau_w = 10d-1
epsilon_tau_X = 10d-1
!Reynolds number
Re = 100.0d0
!magnetic Mach numebr
Mm = 0.001d0
!the ratio of momentum diffusion velocity from the bulk fulid into the boundary layer to the wall/limeter
kappa_X = epsilon_tau_X/(epsilon_delta_X*epsilon_a*Re)

W_s(0) = 0.01d0
domega_theta_s(0) = 0.01d0
domega_z_s(0) = 0.01d0

do i = 0,NN-1
   r_p_s(i) = 1.0d0 +W_s(i)/2.0d0
   r_m_s(i) = 1.0d0 -W_s(i)/2.0d0
   m_p(i) = rho*(r_p_s(i)**4 -r_m_s(i)**4)/4.0d0
   m_t(i) = rho*(r_p_s(i)**2 -r_m_s(i)**2)/2.0d0
!index of the stabilizing effect of the conducting wall on the rotating mode
   cap_delta_wall_p(i) = -2.0d0*dble(m)*((omega*epsilon_tau_w)**2*(r_p_s(i)/epsilon_r_w)**(2*m)* &
                        (1.0d0 -(r_p_s(i)/epsilon_r_w)**(2*m)))/(1.0d0 + &
                         (omega*epsilon_tau_w)**(2*m)*(1.0d0 -(r_p_s(i)/epsilon_r_w)**(2*m))**2)
!tearing mode index
   cap_delta_p_s(i) = cap_delta_mode_p + cap_delta_wall_p(i)
!the component of the total poloidal electromagnetic torque
   T_theta_EM_s (i) = -dble(m**3*n)*((omega*epsilon_tau_w)*(r_p_s(i)/epsilon_r_w)**(2*m))/(1.0d0 + &
                      (omega*epsilon_tau_w)**2*(1 -(r_p_s(i)/epsilon_r_w)**(2*m))**2)*W_s(i)**2
!the component of the total toroidal electromagnetic torque
   T_z_EM_s (i) = -dble(n)/dble(m)*T_theta_EM_s(i)
!the component of the total poloidal viscous torque
   T_theta_VS_s(i) = -domega_theta_s(i)*(r_m_s(i)**3 +r_p_s(i)**3)
!the component of the total toroidal viscous torque
   T_z_VS_s(i) = -domega_z_s(i)/(log(epsilon_a/r_m_s(i)) +kappa_X)
!island width
   W_s(i+1) = W_s(i) +cap_delta_p_s(i)*dt/I_1
!the flux surface averaged change(poloidal)
   domega_theta_s(i+1) = domega_theta_s(i) +(-T_theta_VS_s(i)/(Re*epsilon_delta_D) -T_theta_EM_s(i)/Mm**2)*dt/m_p(i)
!the flux surface averaged change(toroidal)
   domega_z_s(i+1) = domega_z_s(i) +(-T_z_VS_s(i)/Re +T_z_EM_s(i)/Mm**2)*dt/m_t(i)
end do

do i = 0,NN
   t = dble(i)*dt
   write(11,*) t, W_s(i), domega_theta_s(i), domega_z_s(i)
end do

close(11)

end

