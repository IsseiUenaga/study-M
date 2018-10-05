
program isseiuenaga
implicit none
integer :: i
integer,parameter :: m = , n = , NN =
real(8) :: dt, pi, I_1
real(8) :: a, R_0, rho, omega, mu_0, mu_perp
real(8) :: eta_para, sigma_w_zz, kappa_X
real(8) :: tau_D, tau_R, tau_w, tau_X
real(8) :: delta_D, delta_w, delta_X
real(8) :: cap_delta_mode_p(0:NN), cap_delta_wall_p(0:NN), cap_delta_p_s(0:NN)
real(8) :: r_s, r_w, r_p_s(0:NN), r_m_s(0:NN)
real(8) :: B_z, B_theta, B_theta_p
real(8) :: q_s, q_p_s, s_s
real(8) :: T_theta_EM_s(0:NN), T_z_EM_s(0:NN)
real(8) :: T_theta_VS_s(0:NN), T_z_VS_s(0:NN)
real(8) :: W_s(0:NN), domega_theta_s(0:NN), domega_z_s(0:NN)

dt = 1.0d0/dble(NN)
pi = 4.0d0*atan(1.0d0)
I_1 = 0.8227d0
a =  !minor radius
R_0 =  !simulated major radius
rho =  !plasma density
omega =  !angular frequency of rotation of the tearing mode
mu_0 = 4.0d0*pi**2*10d-7 !permeability of vacuum magnetic constant
mu_perp =  !perpendicular viscosity
eta_para =  !electric resistivity along magnetic field lines
sigma_w_zz =  !conductivity tensor of zz component
tau_D =  !poloidal flow dampimg time-scale
tau_X =  !average momentum exchange time between the plasma in the boundary layer and the wall/limiter
tau_R = mu_0*r_s**2/eta_para !resistive time-scale
delta_X =  !thickness of thin boundary layer
delta_w =  !thickness of conducting wall
delta_D = sqrt(mu_perp*tau_D/rho)  !localization length scale
tau_w = mu_0*sigma_w_zz*delta_w*r_w/(2.0d0*dble(m)) !wall time constant
r_s = !rational surface
r_w =  !radius of conducting wall
B_z =  !toroidal magnetic field strength
B_theta =  !poloidal magnetic field
B_theta_p =  !derivative of B_theta
q_s = r_s*B_z/(R_0*B_theta_s)  !safe factor
q_p_s = B_z/(R_0*B_theta_s)*(1.0d0 -r_s*B_theta_p_s/B_theta_s)  !derivative of q
s_s = r_s*q_p_s/q_s  !local magnetic shear
kappa_X = tau_X*mu_perp/(delta_X*a*rho)  !the ratio of momentum
                                         !diffusion velocity from the bulk fulid into
                                         !the boundary layer to the wall/limeter

do i = 0,NN-1
   r_p_s(i) = r_s +W_s(i)
   r_m_s(i) = r_s -W_s(i)
   norm_cap_psi(i) = B_z*W_s(i)**2*norm(s_s)/(16.0d0*R_0*q_s)  !reconected magnetic flux
   cap_delta_mode_p(i) = dble(m)*((r_p_s(i)/r_m_s(i))**m -1.0d0)  !standard stability index for a free boundary tearing mode
   cap_delta_wall_p(i) = -2.0d0*dble(m)*((omega*tau_w)**2*(r_p_s(i)/r_w)**(2*m)* &
                    (1.0d0 -(r_p_s(i)/r_w)**(2*m)))/(1.0d0 +(omega*tau_w)**2* &
                    (1.0d0 -(r_p_s(i)/r_w)**(2*m))**2)  !index of the stabilizing effect of the conducting wall on the rotating mode
   cap_delta_p_s(i) = delta_mode_p(i) + delta_wall_p(i)  !tearing mode index
   T_theta_EM_s (i) = -4.0d0*pi**2*R_0*dble(m)**2/mu_0*(omega*tau_w)*(r_p_s(i)/r_w)**(2*m)/ &
                     (1.0d0 +(omega*tau_w)**2*(1 -(r_p_s(i)/r_w)**(2*m))**2)*norm_cap_psi(i)**2  !the component of the total poloidal electromagnetic torque
   T_z_EM_s(i) = -dble(n)/dble(m)*T_theta_EM_s(i)  !the component of the total toroidal electromagnetic torque
   T_theta_VS_s(i) = -4.0d0*pi**2*R_0*mu_perp*domega_theta_s(i)/delta_D*(r_m_s(i)**3 +r_p_s(i)**3)  !the component of the total poloidal poloidal viscous torque
   T_z_VS_s(i) = -4.0d0*pi**2*R_0**3*mu_perp/(log(a/r_p_s(i)) +kappa_X)  !the component of the total poloidal toroidal viscous torque
   W_s(i+1) = W_s(i) +r_s*dt*delta_p_s(i)/(I_1*tau_R)  !island width
   domega_theta_s(i+1) = domega_theta_s(i) +(T_theta_VS_s(i) +T_theta_EM_s(i))*dt/ &
                        (rho*pi**2*R_0*(r_p_s(i)**4) -r_m_s(i)**4))  !the flux surface averaged change(poloidal)
   domega_z_s(i+1) = domega_z_s(i) +(T_z_VS_s(i) -(dble(n)/dble(m))*T_theta_EM_s(i))/ &
                     (2.0d0*rho*pi**2*R_0**3*(r_p_s**2 -r_m_s**2))  !the flux surface averaged change(toroidal)
end do

end
