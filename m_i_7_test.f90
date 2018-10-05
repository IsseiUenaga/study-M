
program isseiuenaga
implicit none
integer :: i, j
integer,parameter :: m = 2, n = 1, NN = 10**2
real(8) :: dt, dt0, t, t_0, t_N, I_1
real(8) :: delta_0, delta_1, S, error(1:3)
real(8) :: Re, Mm, rho, kappa_X
real(8) :: epsilon_a, epsilon_r_w, epsilon_delta_D, epsilon_tau_w, epsilon_R0
real(8) :: cap_delta_mode_0_p, W_0
real(8) :: W_s, domega_theta_s, domega_z_s
real(8) :: W_s_vec(0:NN), domega_theta_s_vec(0:NN), domega_z_s_vec(0:NN)
real(8) :: W_s_vec_new(0:NN), domega_theta_s_vec_new(0:NN), domega_z_s_vec_new(0:NN)
real(8) :: subf(1:14), RK(1:4,1:3)

open(11,file="result-m_i_7-main.dat",status="replace")
open(12,file="result-m_i_7-sub.dat",status="replace")
open(13,file="result-m_i_7-test1.dat",status="replace")
open(14,file="result-m_i_7-test2.dat",status="replace")

S = 0.9d0
t_0 = 0.0d0
t_N = 1.0d0
dt = (t_N -t_0)/dble(NN)
delta_0 = 10d-3
I_1 = 0.8227d0
rho = 1.0d0
W_0 = 0.2d0
cap_delta_mode_0_p = 2.0d0
epsilon_a = 2.0d0
epsilon_r_w = 2.5d0
epsilon_R0 = 20.0d0
epsilon_delta_D = 0.001d0
epsilon_tau_w = 1.0d0
kappa_X = 0.5d0
Re = 1000.0d0
Mm = 0.001d0

i = 0
W_s_vec(0) = 0.001d0
W_s = W_s_vec(0)
domega_theta_s_vec(0) = 0.01d0
domega_theta_s = domega_theta_s_vec(0)
domega_z_s_vec(0) = 0.01d0
domega_z_s = domega_z_s_vec(0)

call subfunctions(m,n,i,dt,W_s,domega_theta_s, domega_z_s,W_0, cap_delta_mode_0_p, &
                  epsilon_r_w,epsilon_tau_w,epsilon_delta_D,epsilon_a,kappa_X,epsilon_R0, &
                  Re,Mm,rho)
do i = 1,NN
   dt0 = dt
   if(mod(i,5) == 0) then
!      write(*,*) ""
!      write(*,*) "step : ", i, "island width : ", W_s
!      write(*,*) ""
   end if
   write(14,*) i, dt0
   call calculate(W_s,domega_theta_s,domega_z_s,W_s_vec,domega_theta_s_vec,domega_z_s_vec, &
                  W_s_vec_new,domega_theta_s_vec_new,domega_z_s_vec_new,rho,kappa_X,Re,Mm, &
                  epsilon_delta_D,epsilon_tau_w,epsilon_r_w,epsilon_a, &
                  epsilon_R0,RK,I_1,dt0,cap_delta_mode_0_p,W_0,m,n,i,NN)
   do
     write(14,*) i, dt0
     error(1:3) = (/abs(W_s_vec_new(i) -W_s_vec(i)), abs(domega_theta_s_vec_new(i) -domega_theta_s_vec(i)), &
                   abs(domega_z_s_vec_new(i) -domega_z_s_vec(i))/)
     delta_1 = sqrt(dot_product(error,error))
     if(delta_1 < delta_0*dt0/(t_N-t_0)) then
        write(13,*) i, "delta_1 =", delta_1, "delta_0*dt0/(t_N-t_0) =", delta_0*dt0/(t_N-t_0) &
                     , "delta_1 < delta_0*dt0/(t_N-t_0)"
        exit
     elseif(delta_1 > delta_0*dt0/(t_N-t_0)) then
            write(13,*) i, "delta_1 =", delta_1, "delta_0*dt0/(t_N-t_0) =", delta_0*dt0/(t_N-t_0) &
                         , "delta_1 > delta_0*dt0/(t_N-t_0)"
            if(delta_1 > delta_0) then
               dt0 = S*dt0*abs(delta_0/delta_1)**(0.2)
               write(13,*) i, "delta_0= ", delta_0, "delta1>delta0"
               call calculate(W_s,domega_theta_s,domega_z_s,W_s_vec,domega_theta_s_vec,domega_z_s_vec, &
                              W_s_vec_new,domega_theta_s_vec_new,domega_z_s_vec_new,rho,kappa_X,Re,Mm, &
                              epsilon_delta_D,epsilon_tau_w,epsilon_r_w,epsilon_a, &
                              epsilon_R0,RK,I_1,dt0,cap_delta_mode_0_p,W_0,m,n,i,NN)
            elseif(delta_1 < delta_0) then
                   dt0 = S*dt0*abs(delta_0/delta_1)**(0.25)
                   write(13,*) i, "delta_0 =", delta_0, "delta1<delta0"
                   call calculate(W_s,domega_theta_s,domega_z_s,W_s_vec,domega_theta_s_vec,domega_z_s_vec, &
                                  W_s_vec_new,domega_theta_s_vec_new,domega_z_s_vec_new,rho,kappa_X,Re,Mm, &
                                  epsilon_delta_D,epsilon_tau_w,epsilon_r_w,epsilon_a, &
                                  epsilon_R0,RK,I_1,dt0,cap_delta_mode_0_p,W_0,m,n,i,NN)
            end if
     end if
   write(13,*) ""
   end do
   write(13,*) ""
   W_s = W_s_vec(i)
   domega_theta_s = domega_theta_s_vec(i)
   domega_z_s = domega_z_s_vec(i)
   call subfunctions(m,m,i,dt,W_s,domega_theta_s, domega_z_s,W_0, cap_delta_mode_0_p, &
                     epsilon_r_w,epsilon_tau_w,epsilon_delta_D,epsilon_a,kappa_X,epsilon_R0, &
                     Re,Mm,rho)
end do

do i = 0,NN
   t = dble(i)*dt
   write(11,*) t, W_s_vec(i), domega_theta_s_vec(i), domega_z_s_vec(i)
end do

close(11)
close(12)
close(13)
close(14)

end program isseiuenaga


subroutine Runge_Kutta_method_of_4th_order(W_s,domega_theta_s,domega_z_s,rho,kappa_X,Re,Mm, &
                                           epsilon_delta_D,epsilon_tau_w,epsilon_r_w,epsilon_a, &
                                           epsilon_R0,RK,I_1,dt,cap_delta_mode_0_p,W_0,m,n)
implicit none
integer,intent(in) :: m, n
real(8),intent(in) :: W_s, domega_theta_s, domega_z_s
real(8),intent(in) :: rho, kappa_X, Re, Mm, I_1, dt, cap_delta_mode_0_p, W_0
real(8),intent(in) :: epsilon_delta_D,epsilon_tau_w,epsilon_r_w,epsilon_a,epsilon_R0
real(8) :: r_p_s, r_m_s, omega
real(8) :: cap_delta_mode_p, cap_delta_wall_p
real(8) :: T_theta_VS_s, T_theta_EM_s
real(8) :: T_z_VS_s, T_z_EM_s
real(8) :: m_p, m_t
real(8) :: f_1, f_2, f_3
real(8) :: r_p_s_vec(1:4), r_m_s_vec(1:4), omega_vec(1:4)
real(8) :: cap_delta_mode_p_vec(1:4), cap_delta_wall_p_vec(1:4)
real(8) :: T_theta_VS_s_vec(1:4), T_theta_EM_s_vec(1:4)
real(8) :: T_z_VS_s_vec(1:4), T_z_EM_s_vec(1:4)
real(8) :: m_p_vec(1:4), m_t_vec(1:4)
real(8),intent(out) :: RK(1:4,1:3)

r_p_s_vec(1) = r_p_s(W_s)
r_m_s_vec(1) = r_m_s(W_s)
omega_vec(1) = omega(domega_theta_s,domega_z_s,m,n)
cap_delta_mode_p_vec(1) = cap_delta_mode_p(W_s,W_0,cap_delta_mode_0_p)
cap_delta_wall_p_vec(1) = cap_delta_wall_p(omega_vec(1),r_p_s_vec(1),epsilon_r_w,epsilon_tau_w,m)
T_theta_VS_s_vec(1) = T_theta_VS_s(domega_theta_s,r_p_s_vec(1),r_m_s_vec(1),Re,epsilon_delta_D)
T_theta_EM_s_vec(1) = T_theta_EM_s(W_s,omega_vec(1),r_p_s_vec(1),epsilon_r_w,epsilon_tau_w,Mm,m)
T_z_VS_s_vec(1) = T_z_VS_s(domega_z_s,r_m_s_vec(1),epsilon_a,kappa_X,re)
T_z_EM_s_vec(1) = T_z_EM_s(omega_vec(1),r_p_s_vec(1),W_s,Mm,epsilon_r_w,epsilon_tau_w,epsilon_R0,m,n)
m_p_vec(1) = m_p(r_p_s_vec(1),r_m_s_vec(1),rho)
m_t_vec(1) = m_t(r_p_s_vec(1),r_m_s_vec(1),rho)
RK(1,1) = f_1(cap_delta_mode_p_vec(1),cap_delta_wall_p_vec(1),I_1)
RK(1,2) = f_2(T_theta_VS_s_vec(1),T_theta_EM_s_vec(1),m_p_vec(1))
RK(1,3) = f_3(T_z_VS_s_vec(1),T_z_EM_s_vec(1),m_t_vec(1))

r_p_s_vec(2) = r_p_s(W_s+dt*RK(1,1)/2.0d0)
r_m_s_vec(2) = r_m_s(W_s+dt*RK(1,1)/2.0d0)
omega_vec(2) = omega(domega_theta_s+dt*RK(1,2)/2.0d0,domega_z_s+dt*RK(1,3)/2.0d0,m,n)
cap_delta_mode_p_vec(2) = cap_delta_mode_p(W_s+dt*RK(1,1)/2.0d0,W_0,cap_delta_mode_0_p)
cap_delta_wall_p_vec(2) = cap_delta_wall_p(omega_vec(2),r_p_s_vec(2),epsilon_r_w,epsilon_tau_w,m)
T_theta_VS_s_vec(2) = T_theta_VS_s(domega_theta_s+dt*RK(1,2)/2.0d0,r_p_s_vec(2),r_m_s_vec(2),Re,epsilon_delta_D)
T_theta_EM_s_vec(2) = T_theta_EM_s(W_s+dt*RK(1,1)/2.0d0,omega_vec(2),r_p_s_vec(2),epsilon_r_w,epsilon_tau_w,Mm,m)
T_z_VS_s_vec(2) = T_z_VS_s(domega_z_s+dt*RK(1,3)/2.0d0,r_m_s_vec(2),epsilon_a,kappa_X,re)
T_z_EM_s_vec(2) = T_z_EM_s(omega_vec(2),r_p_s_vec(2),W_s+dt*RK(1,1)/2.0d0,Mm,epsilon_r_w,epsilon_tau_w,epsilon_R0,m,n)
m_p_vec(2) = m_p(r_p_s_vec(2),r_m_s_vec(2),rho)
m_t_vec(2) = m_t(r_p_s_vec(2),r_m_s_vec(2),rho)
RK(2,1) = f_1(cap_delta_mode_p_vec(2),cap_delta_wall_p_vec(2),I_1)
RK(2,2) = f_2(T_theta_VS_s_vec(2),T_theta_EM_s_vec(2),m_p_vec(2))
RK(2,3) = f_3(T_z_VS_s_vec(2),T_z_EM_s_vec(2),m_t_vec(2))

r_p_s_vec(3) = r_p_s(W_s+dt*RK(2,1)/2.0d0)
r_m_s_vec(3) = r_m_s(W_s+dt*RK(2,1)/2.0d0)
omega_vec(3) = omega(domega_theta_s+dt*RK(2,2)/2.0d0,domega_z_s+dt*RK(2,3)/2.0d0,m,n)
cap_delta_mode_p_vec(3) = cap_delta_mode_p(W_s+dt*RK(2,1)/2.0d0,W_0,cap_delta_mode_0_p)
cap_delta_wall_p_vec(3) = cap_delta_wall_p(omega_vec(3),r_p_s_vec(3),epsilon_r_w,epsilon_tau_w,m)
T_theta_VS_s_vec(3) = T_theta_VS_s(domega_theta_s+dt*RK(2,2)/2.0d0,r_p_s_vec(3),r_m_s_vec(3),Re,epsilon_delta_D)
T_theta_EM_s_vec(3) = T_theta_EM_s(W_s+dt*RK(2,1)/2.0d0,omega_vec(3),r_p_s_vec(3),epsilon_r_w,epsilon_tau_w,Mm,m)
T_z_VS_s_vec(3) = T_z_VS_s(domega_z_s+dt*RK(2,3)/2.0d0,r_m_s_vec(3),epsilon_a,kappa_X,re)
T_z_EM_s_vec(3) = T_z_EM_s(omega_vec(3),r_p_s_vec(3),W_s+dt*RK(2,1)/2.0d0,Mm,epsilon_r_w,epsilon_tau_w,epsilon_R0,m,n)
m_p_vec(3) = m_p(r_p_s_vec(3),r_m_s_vec(3),rho)
m_t_vec(3) = m_t(r_p_s_vec(3),r_m_s_vec(3),rho)
RK(3,1) = f_1(cap_delta_mode_p_vec(3),cap_delta_wall_p_vec(3),I_1)
RK(3,2) = f_2(T_theta_VS_s_vec(3),T_theta_EM_s_vec(3),m_p_vec(3))
RK(3,3) = f_3(T_z_VS_s_vec(3),T_z_EM_s_vec(3),m_t_vec(3))

r_p_s_vec(4) = r_p_s(W_s+dt*RK(3,1))
r_m_s_vec(4) = r_m_s(W_s+dt*RK(3,1))
omega_vec(4) = omega(domega_theta_s+dt*RK(3,2),domega_z_s+dt*RK(3,3),m,n)
cap_delta_mode_p_vec(4) = cap_delta_mode_p(W_s+dt*RK(3,1),W_0,cap_delta_mode_0_p)
cap_delta_wall_p_vec(4) = cap_delta_wall_p(omega_vec(4),r_p_s_vec(4),epsilon_r_w,epsilon_tau_w,m)
T_theta_VS_s_vec(4) = T_theta_VS_s(domega_theta_s+dt*RK(3,2),r_p_s_vec(4),r_m_s_vec(4),Re,epsilon_delta_D)
T_theta_EM_s_vec(4) = T_theta_EM_s(W_s+dt*RK(3,1),omega_vec(4),r_p_s_vec(4),epsilon_r_w,epsilon_tau_w,Mm,m)
T_z_VS_s_vec(4) = T_z_VS_s(domega_z_s+dt*RK(3,3),r_m_s_vec(4),epsilon_a,kappa_X,re)
T_z_EM_s_vec(4) = T_z_EM_s(omega_vec(4),r_p_s_vec(4),W_s,Mm,epsilon_r_w,epsilon_tau_w,epsilon_R0,m,n)
m_p_vec(4) = m_p(r_p_s_vec(4),r_m_s_vec(4),rho)
m_t_vec(4) = m_t(r_p_s_vec(4),r_m_s_vec(4),rho)
RK(4,1) = f_1(cap_delta_mode_p_vec(4),T_theta_EM_s_vec(4),I_1)
RK(4,2) = f_2(T_theta_VS_s_vec(4),T_theta_EM_s_vec(4),m_p_vec(4))
RK(4,3) = f_3(T_z_VS_s_vec(4),T_z_EM_s_vec(4),m_t_vec(4))

end subroutine Runge_Kutta_method_of_4th_order


subroutine calculate(W_s,domega_theta_s,domega_z_s,W_s_vec,domega_theta_s_vec,domega_z_s_vec, &
                     W_s_vec_new,domega_theta_s_vec_new,domega_z_s_vec_new,rho,kappa_X,Re,Mm, &
                     epsilon_delta_D,epsilon_tau_w,epsilon_r_w,epsilon_a, &
                     epsilon_R0,RK,I_1,dt0,cap_delta_mode_0_p,W_0,m,n,i,NN)
implicit none
integer,intent(in) :: m, n, i, NN
real(8) :: W_s, domega_theta_s, domega_z_s, RK(1:4,1:3)
real(8) :: dt0, dt1
real(8),intent(in) :: rho, kappa_X, Re, Mm
real(8),intent(in) :: epsilon_delta_D, epsilon_tau_w, epsilon_r_w, epsilon_a, epsilon_R0
real(8),intent(in) :: I_1, cap_delta_mode_0_p, W_0
real(8) :: W_s_vec(0:NN), domega_theta_s_vec(0:NN), domega_z_s_vec(0:NN)
real(8) :: W_s_vec_new(0:NN), domega_theta_s_vec_new(0:NN), domega_z_s_vec_new(0:NN)

call Runge_Kutta_method_of_4th_order(W_s,domega_theta_s,domega_z_s,rho,kappa_X,Re,Mm, &
                                     epsilon_delta_D,epsilon_tau_w,epsilon_r_w,epsilon_a, &
                                     epsilon_R0,RK,I_1,dt0,cap_delta_mode_0_p,W_0,m,n)
W_s_vec(i) = W_s +dt0*(RK(1,1) +2.0d0*RK(2,1) +2.0d0*RK(3,1) +RK(4,1))/6.0d0
domega_theta_s_vec(i) = domega_theta_s +dt0*(RK(1,2) +2.0d0*RK(2,2) +2.0d0*RK(3,2) +RK(4,2))/6.0d0
domega_z_s_vec(i) = domega_z_s +dt0*(RK(1,3) +2.0d0*RK(2,3) +2.0d0*RK(3,3) +RK(4,3))/6.0d0

dt1 = dt0/2.0d0
call Runge_Kutta_method_of_4th_order(W_s,domega_theta_s,domega_z_s,rho,kappa_X,Re,Mm, &
                                     epsilon_delta_D,epsilon_tau_w,epsilon_r_w,epsilon_a, &
                                     epsilon_R0,RK,I_1,dt1,cap_delta_mode_0_p,W_0,m,n)
W_s_vec_new(i) = W_s +dt1*(RK(1,1) +2.0d0*RK(2,1) +2.0d0*RK(3,1) +RK(4,1))/6.0d0
domega_theta_s_vec_new(i) = domega_theta_s +dt1*(RK(1,2) +2.0d0*RK(2,2) +2.0d0*RK(3,2) +RK(4,2))/6.0d0
domega_z_s_vec_new(i) = domega_z_s +dt1*(RK(1,3) +2.0d0*RK(2,3) +2.0d0*RK(3,3) +RK(4,3))/6.0d0
W_s = W_s_vec_new(i)
domega_theta_s = domega_theta_s_vec_new(i)
domega_z_s = domega_z_s_vec_new(i)
call Runge_Kutta_method_of_4th_order(W_s,domega_theta_s,domega_z_s,rho,kappa_X,Re,Mm, &
                                     epsilon_delta_D,epsilon_tau_w,epsilon_r_w,epsilon_a, &
                                     epsilon_R0,RK,I_1,dt1,cap_delta_mode_0_p,W_0,m,n)
W_s_vec_new(i) = W_s +dt1*(RK(1,1) +2.0d0*RK(2,1) +2.0d0*RK(3,1) +RK(4,1))/6.0d0
domega_theta_s_vec_new(i) = domega_theta_s +dt1*(RK(1,2) +2.0d0*RK(2,2) +2.0d0*RK(3,2) +RK(4,2))/6.0d0
domega_z_s_vec_new(i) = domega_z_s +dt1*(RK(1,3) +2.0d0*RK(2,3) +2.0d0*RK(3,3) +RK(4,3))/6.0d0

end subroutine calculate


subroutine subfunctions(m,n,i,dt,W_s,domega_theta_s,domega_z_s,W_0,cap_delta_mode_0_p, &
epsilon_r_w,epsilon_tau_w,epsilon_delta_D,epsilon_a,kappa_X,epsilon_R0, &
Re,Mm,rho)
implicit none
integer :: m, n, i
real(8) :: t, dt, Re, Mm, rho
real(8) :: W_s, domega_theta_s, domega_z_s
real(8) :: W_0, cap_delta_mode_0_p
real(8) :: epsilon_r_w, epsilon_tau_w, epsilon_delta_D, epsilon_a,kappa_X, epsilon_R0
real(8) :: r_p_s, r_m_s, omega
real(8) :: cap_delta_mode_p, cap_delta_wall_p
real(8) :: T_theta_VS_s, T_theta_EM_s
real(8) :: T_z_VS_s, T_z_EM_s
real(8) :: m_p, m_t
real(8) :: subf(1:11)

t = dble(i)*dt
subf(1) = r_p_s(W_s)
subf(2) = r_m_s(W_s)
subf(3) = omega(domega_theta_s,domega_z_s,m,n)
subf(4) = cap_delta_mode_p(W_s,W_0,cap_delta_mode_0_p)
subf(5) = cap_delta_wall_p(subf(3),subf(1),epsilon_r_w,epsilon_tau_w,m)
subf(6) = T_theta_VS_s(domega_theta_s,subf(1),subf(2),Re,epsilon_delta_D)
subf(7) = T_theta_EM_s(W_s,subf(3),subf(1),epsilon_r_w,epsilon_tau_w,Mm,m)
subf(8) = T_z_VS_s(domega_z_s,subf(2),epsilon_a,kappa_X,re)
subf(9) = T_z_EM_s(subf(3),subf(1),W_s,Mm,epsilon_r_w,epsilon_tau_w,epsilon_R0,m,n)
subf(10) = m_p(subf(1),subf(2),rho)
subf(11) = m_t(subf(1),subf(2),rho)

write(12,*) t, subf(1:11)

end subroutine subfunctions


function r_p_s(W_s)
real(8),intent(in) :: W_s
real(8) :: r_p_s
r_p_s = 1.0d0 +W_s/2.0d0
end function r_p_s

function r_m_s(W_s)
real(8),intent(in) :: W_s
real(8) :: r_m_s
r_m_s = 1.0d0 -W_s/2.0d0
end function r_m_s

function omega(domega_theta_s,domega_z_s,m,n)
integer,intent(in) :: m, n
real(8),intent(in) :: domega_theta_s, domega_z_s
real(8) :: omega
omega = m*domega_theta_s -n*domega_z_s
end function omega

function cap_delta_mode_p(W_s,W_0,cap_delta_mode_0_p)
real(8),intent(in) :: W_s
real(8),intent(in) :: cap_delta_mode_0_p, W_0
real(8) :: cap_delta_mode_p
cap_delta_mode_p = cap_delta_mode_0_p*(1.0d0 -W_s/W_0)
end function cap_delta_mode_p

function cap_delta_wall_p(omega,r_p_s,epsilon_r_w,epsilon_tau_w,m)
integer,intent(in) :: m
real(8),intent(in) :: omega, r_p_s
real(8),intent(in) :: epsilon_r_w, epsilon_tau_w
real(8) :: cap_delta_wall_p
cap_delta_wall_p = -2.0d0*dble(m)*(omega*epsilon_tau_w)**2*(r_p_s/epsilon_r_w)**(2*m)* &
(1.0d0 -(r_p_s/epsilon_r_w)**(2*m))/(1.0d0 + &
(omega*epsilon_tau_w)**(2*m)*(1.0d0 -(r_p_s/epsilon_r_w)**(2*m))**2)
end function cap_delta_wall_p

function T_theta_VS_s(domega_theta_s,r_p_s,r_m_s,Re,epsilon_delta_D)
real(8),intent(in) :: domega_theta_s, r_p_s, r_m_s
real(8),intent(in) :: Re, epsilon_delta_D
real(8) :: T_theta_VS_s
T_theta_VS_s = -domega_theta_s*(r_m_s**3 +r_p_s**3)/(Re*epsilon_delta_D)
end function T_theta_VS_s

function T_theta_EM_s(W_s,omega,r_p_s,epsilon_r_w,epsilon_tau_w,Mm,m)
integer,intent(in) :: m
real(8),intent(in) :: omega, r_p_s, W_s
real(8),intent(in) :: Mm, epsilon_r_w, epsilon_tau_w
real(8) :: T_theta_EM_s
T_theta_EM_s = -m**2*((omega*epsilon_tau_w)*(r_p_s/epsilon_r_w)**(2*m))/ &
(1.0d0 +(omega*epsilon_tau_w)**2*(1 -(r_p_s/epsilon_r_w)**(2*m))**2)* &
W_s**4/Mm**2
end function T_theta_EM_s

function T_z_VS_s(domega_z_s,r_m_s,epsilon_a,kappa_X,re)
real(8),intent(in) :: domega_z_s, r_m_s
real(8),intent(in) :: epsilon_a, kappa_X, re
real(8) :: T_z_VS_s
T_z_VS_s = -domega_z_s/((log(epsilon_a/r_m_s) +kappa_X)*re)
end function T_z_VS_s

function T_z_EM_s(omega,r_p_s,W_s,Mm,epsilon_r_w,epsilon_tau_w,epsilon_R0,m,n)
integer,intent(in) :: m, n
real(8),intent(in) :: omega, r_p_s, W_s
real(8),intent(in) :: Mm, epsilon_r_w, epsilon_tau_w, epsilon_R0
real(8) :: T_z_EM_s
T_z_EM_s = m*n*((omega*epsilon_tau_w)*(r_p_s/epsilon_r_w)**(2*m))/ &
(1.0d0 +(omega*epsilon_tau_w)**2*(1 -(r_p_s/epsilon_r_w)**(2*m))**2)* &
W_s**4/(Mm**2*epsilon_R0**2)
end function T_z_EM_s

function m_p(r_p_s,r_m_s,rho)
real(8),intent(in) :: r_p_s, r_m_s
real(8),intent(in) :: rho
real(8) :: m_p
m_p = rho*(r_p_s**4 -r_m_s**4)/4.0d0
end function m_p

function m_t(r_p_s,r_m_s,rho)
real(8),intent(in) :: r_p_s, r_m_s
real(8),intent(in) :: rho
real(8) :: m_t
m_t = rho*(r_p_s**2 -r_m_s**2)/2.0d0
end function m_t

function f_1(cap_delta_mode_p,cap_delta_wall_p,I_1)
real(8),intent(in) :: cap_delta_mode_p, cap_delta_wall_p
real(8),intent(in) :: I_1
real(8) :: f_1
f_1 = (cap_delta_mode_p +cap_delta_wall_p)/I_1
end function f_1

function f_2(T_theta_VS_s,T_theta_EM_s,m_p)
real(8),intent(in) :: T_theta_VS_s, T_theta_EM_s, m_p
real(8) :: f_2
f_2 = (T_theta_VS_s +T_theta_EM_s)/m_p
end function f_2

function f_3(T_z_VS_s,T_z_EM_s,m_t)
real(8),intent(in) :: T_z_VS_s, T_z_EM_s, m_t
real(8) :: f_3
f_3 = (T_z_VS_s +T_z_EM_s)/m_t
end function f_3
