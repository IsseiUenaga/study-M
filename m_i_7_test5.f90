
program isseiuenaga
implicit none
integer(8) :: i, j
integer(8),parameter :: m = 2, n = 1, NN = 10**6, ite = 10**9
real(8) :: dt, dt1, dt2, t, I_1, beta
real(8) :: S, t0, tN, tnet
real(8) :: Re, Mm, rho, kappa_X, omega_0_vec(0:50)
real(8) :: epsilon_a, epsilon_r_w, epsilon_delta_D, epsilon_tau_w, epsilon_R0
real(8) :: cap_delta_mode_0_p, W_0, omega_0, domega_0
real(8) :: delta0, delta1, delta
real(8) :: W_s0, domega_theta_s0, domega_z_s0
real(8) :: W_s, domega_theta_s, domega_z_s
real(8) :: W_s1, domega_theta_s1, domega_z_s1
real(8) :: W_s2, domega_theta_s2, domega_z_s2
real(8),allocatable :: W_s_vec(:,:), domega_theta_s_vec(:,:)
real(8),allocatable :: domega_z_s_vec(:,:), omega_vec(:,:)
real(8) :: subf(1:11), RK(1:4,1:3), error(1:3)
character(len=50) :: tmp1, tmp2, filename1, filename2, filename3

allocate(W_s_vec(0:ite,0:50))
allocate(domega_theta_s_vec(0:ite,0:50))
allocate(domega_z_s_vec(0:ite,0:50))
allocate(omega_vec(0:ite,0:50))

call set_parameter(m,NN,tN,dt,domega_0,delta0,W_0,beta,Re,Mm,I_1,rho,kappa_X, &
                   epsilon_a,epsilon_r_w,epsilon_R0,epsilon_delta_D,epsilon_tau_w, &
                   cap_delta_mode_0_p)

write(tmp1,'(1F4.2)') beta
write(tmp2,'(1F4.2)') W_0

filename1 = "result-m_i_7-main_beta="//trim(adjustl(tmp1))//"_W0="//trim(adjustl(tmp2))//".dat"
filename2 = "result-m_i_7-sub_beta="//trim(adjustl(tmp1))//"_W0="//trim(adjustl(tmp2))//".dat"
filename3 = "result-m_i_7-steady_beta="//trim(adjustl(tmp1))//"_W0="//trim(adjustl(tmp2))//".dat"

open(11,file=filename1,status="replace")
open(12,file=filename2,status="replace")
open(13,file=filename3,status="replace")

do j = 0,50

   call initial_calculation(m,n,j,ite,t,dt,omega_0,domega_0,W_s_vec,domega_theta_s_vec,domega_z_s_vec, &
                            omega_0_vec,W_s0,domega_theta_s0,domega_z_s0,W_0,epsilon_r_w,epsilon_tau_w, &
                            epsilon_delta_D,epsilon_a,kappa_X,epsilon_R0,Re,Mm,rho,subf,cap_delta_mode_0_p)

   call calculation(m,n,ite,j,dt,tN,delta0,W_s0,domega_theta_s0,domega_z_s0,rho,kappa_X,Re,Mm,I_1,W_0,omega_0, &
                    epsilon_delta_D,epsilon_tau_w,epsilon_r_w,epsilon_a,epsilon_R0,cap_delta_mode_0_p)

end do

deallocate(W_s_vec)
deallocate(domega_theta_s_vec)
deallocate(domega_z_s_vec)
deallocate(omega_vec)

close(11)
close(12)
close(13)
close(14)

end program isseiuenaga


subroutine set_parameter(m,NN,tN,dt,domega_0,delta0,W_0,beta,Re,Mm,I_1,rho,kappa_X, &
                         epsilon_a,epsilon_r_w,epsilon_R0,epsilon_delta_D,epsilon_tau_w, &
                         cap_delta_mode_0_p)
implicit none
integer(8),intent(in) :: m, NN
real(8) :: t0, tN, tnet, dt
real(8) :: domega_0, delta0
real(8) :: W_0, beta
real(8) :: Re, Mm
real(8) :: I_1, rho, kappa_X
real(8) :: epsilon_a, epsilon_r_w
real(8) :: epsilon_R0, epsilon_delta_D, epsilon_tau_w
real(8) :: cap_delta_mode_0_p

t0 = 0.0d0
tN = 700.0d0
tnet = tN -t0
dt = tnet/dble(NN)

domega_0 = 10.0d0/50.0d0

delta0 = 10d-4

write(*,*) "beta ="
read(*,*) beta
write(*,*) "W0 ="
read(*,*) W_0
write(*,*) ""
write(*,*) ""

Re = 1000.0d0
Mm = 0.001d0

I_1 = 0.8227d0
rho = 1.0d0
epsilon_a = 2.0d0
epsilon_r_w = 2.5d0
epsilon_R0 = 20.0d0
epsilon_delta_D = 0.001d0
epsilon_tau_w = 1.0d0
kappa_X = 0.5d0

cap_delta_mode_0_p = 2.0d0*m/((1.0d0 -beta)*(epsilon_r_w**(2*m) -1.0d0))

end subroutine set_parameter

subroutine initial_conditions(j,ite,t,W_s0,domega_theta_s0,domega_z_s0,omega_0,domega_0, &
                              W_s_vec,domega_theta_s_vec,domega_z_s_vec,omega_0_vec)
implicit none
integer(8),intent(in) :: j, ite
real(8) :: t, W_s0, domega_theta_s0, domega_z_s0, omega_0, domega_0
real(8) :: W_s_vec(0:ite,0:50), domega_theta_s_vec(0:ite,0:50)
real(8) :: domega_z_s_vec(0:ite,0:50), omega_0_vec(0:50)

t = 0.0d0
W_s0 = 0.001d0
domega_theta_s0 = 0.01d0
domega_z_s0 = 0.01d0
omega_0 = dble(j)*domega_0

W_s_vec(0,j) = W_s0
domega_theta_s_vec(0,j) = domega_theta_s0
domega_z_s_vec(0,j) = domega_z_s0
omega_0_vec(j) = omega_0

end subroutine initial_conditions

subroutine initial_calculation(m,n,j,ite,t,dt,omega_0,domega_0,W_s_vec,domega_theta_s_vec,domega_z_s_vec, &
                               omega_0_vec,W_s0,domega_theta_s0,domega_z_s0,W_0,epsilon_r_w,epsilon_tau_w, &
                               epsilon_delta_D,epsilon_a,kappa_X,epsilon_R0,Re,Mm,rho,subf,cap_delta_mode_0_p)
implicit none
integer(8) :: m, n, j, ite
real(8) :: t, dt
real(8) :: W_s0, domega_theta_s0, domega_z_s0
real(8) :: rho, kappa_X, Re, Mm, W_0, omega_0, domega_0
real(8) :: epsilon_delta_D, epsilon_tau_w, epsilon_r_w, epsilon_a, epsilon_R0
real(8) :: subf(1:11), cap_delta_mode_0_p
real(8) :: W_s_vec(0:ite,0:50), domega_theta_s_vec(0:ite,0:50)
real(8) :: domega_z_s_vec(0:ite,0:50), omega_vec(0:ite,0:50), omega_0_vec(0:50)

call initial_conditions(j,ite,t,W_s0,domega_theta_s0,domega_z_s0,omega_0,domega_0, &
                        W_s_vec,domega_theta_s_vec,domega_z_s_vec,omega_0_vec)

write(*,*) "STEP :", j, "omega_0 =", omega_0
write(*,*) ""

write(11,'(1X,4F15.7)') t, W_s0, domega_theta_s0, domega_z_s0

call subfunctions(m,n,dt,W_s0,domega_theta_s0,domega_z_s0,W_0,cap_delta_mode_0_p,omega_0, &
                  epsilon_r_w,epsilon_tau_w,epsilon_delta_D,epsilon_a,kappa_X,epsilon_R0, &
                  Re,Mm,rho,subf,t)

write(12,'(1X,12F15.7)') t, subf(1:11)

write(*,*) "TIME :", t, "island width :", W_s0, "island frequency :", subf(3)

end subroutine initial_calculation

subroutine calculation(m,n,ite,j,dt,tN,delta0,W_s0,domega_theta_s0,domega_z_s0,rho,kappa_X,Re,Mm,I_1,W_0,omega_0, &
                        epsilon_delta_D,epsilon_tau_w,epsilon_r_w,epsilon_a,epsilon_R0,cap_delta_mode_0_p)
implicit none
integer(8) :: m, n, ite, i, j
real(8) :: dt, dt1, dt2, t, tN
real(8) :: delta, delta0, delta1, delta2, error(1:3)
real(8) :: W_s0, domega_theta_s0, domega_z_s0
real(8) :: W_s, domega_theta_s, domega_z_s
real(8) :: W_s1, domega_theta_s1, domega_z_s1
real(8) :: W_s2, domega_theta_s2, domega_z_s2
real(8) :: rho, kappa_X, Re, Mm, I_1, W_0, omega_0, omega_result(0:50)
real(8) :: epsilon_delta_D,epsilon_tau_w,epsilon_r_w,epsilon_a,epsilon_R0
real(8) :: cap_delta_mode_0_p, subf(1:11), omega_0_vec(0:50)
real(8) :: W_s_vec(0:ite,0:50), domega_theta_s_vec(0:ite,0:50)
real(8) :: domega_z_s_vec(0:ite,0:50), omega_vec(0:ite,0:50)

do i = 1,ite
   dt1 = dt
   do
     call step_doubling(W_s0,domega_theta_s0,domega_z_s0,W_s,domega_theta_s,domega_z_s, &
                        W_s1,domega_theta_s1,domega_z_s1,W_s2,domega_theta_s2,domega_z_s2, &
                        rho,kappa_X,Re,Mm,epsilon_delta_D,epsilon_tau_w,epsilon_r_w,epsilon_a, &
                        epsilon_R0,I_1,cap_delta_mode_0_p,W_0,omega_0,m,n,i,dt1,dt2)

     error(1:3) = (/abs(W_s2 -W_s1), abs(domega_theta_s2 -domega_theta_s1), abs(domega_z_s2 -domega_z_s1)/)
     delta1 = sqrt(dot_product(error,error))

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
   end do

   t = t +dt1
   if(t >= tN) exit

   call subfunctions(m,n,dt1,W_s1,domega_theta_s1,domega_z_s1,W_0,cap_delta_mode_0_p,omega_0, &
                     epsilon_r_w,epsilon_tau_w,epsilon_delta_D,epsilon_a,kappa_X,epsilon_R0, &
                     Re,Mm,rho,subf,t)

   W_s_vec(i,j) = W_s1
   domega_theta_s_vec(i,j) = domega_theta_s1
   domega_z_s_vec(i,j) = domega_z_s1
   omega_vec(i,j) = subf(3)

   if(abs((omega_vec(i,j)-omega_vec(i-1,j))/omega_vec(i,j))>1.0d-10) then
      write(11,'(1X,4E20.10e2)') t, W_s1, domega_theta_s1, domega_z_s1
      write(12,'(1X,12E20.10e2)') t, subf(1:11)
   end if

   if((tN-10.0d0<t .and. t<tN) .or. abs((omega_vec(i,j)-omega_vec(i-1,j))/omega_vec(i,j))<1.0d-10) then
        write(13,*) omega_0, omega_vec(i,j)
        exit
   end if

   W_s0 = W_s1
   domega_theta_s0 = domega_theta_s1
   domega_z_s0 = domega_z_s1

   if(mod(i,10000) == 0 .and. abs((omega_vec(i,j)-omega_vec(i-1,j))/omega_vec(i,j))>1.0d-10) then
      write(*,*) "TIME :", t, "island width :", W_s, "island frequency :", subf(3)
   end if

end do

write(11,*) ""
write(12,*) ""
write(*,*) ""
write(*,*) ""

end subroutine calculation

subroutine Runge_Kutta_method_of_4th_order(W_s0,domega_theta_s0,domega_z_s0,W_s,domega_theta_s,domega_z_s, &
                                           rho,kappa_X,Re,Mm,epsilon_delta_D,epsilon_tau_w,epsilon_r_w, &
                                           epsilon_a,epsilon_R0,I_1,dt,cap_delta_mode_0_p,W_0,omega_0,m,n,i)
implicit none
integer(8) :: m, n, i
real(8) :: W_s0, domega_theta_s0, domega_z_s0
real(8) :: W_s, domega_theta_s, domega_z_s
real(8) :: rho, kappa_X, Re, Mm, I_1, dt, cap_delta_mode_0_p, W_0, omega_0
real(8) :: epsilon_delta_D,epsilon_tau_w,epsilon_r_w,epsilon_a,epsilon_R0
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
real(8) :: RK(1:4,1:3)

r_p_s_vec(1) = r_p_s(W_s0)
r_m_s_vec(1) = r_m_s(W_s0)
omega_vec(1) = omega(domega_theta_s0,domega_z_s0,omega_0,m,n)
cap_delta_mode_p_vec(1) = cap_delta_mode_p(W_s0,W_0,cap_delta_mode_0_p)
cap_delta_wall_p_vec(1) = cap_delta_wall_p(omega_vec(1),r_p_s_vec(1),epsilon_r_w,epsilon_tau_w,m)
T_theta_VS_s_vec(1) = T_theta_VS_s(domega_theta_s0,r_p_s_vec(1),r_m_s_vec(1),Re,epsilon_delta_D)
T_theta_EM_s_vec(1) = T_theta_EM_s(W_s0,omega_vec(1),r_p_s_vec(1),epsilon_r_w,epsilon_tau_w,Mm,m)
T_z_VS_s_vec(1) = T_z_VS_s(domega_z_s0,r_m_s_vec(1),epsilon_a,kappa_X,re)
T_z_EM_s_vec(1) = T_z_EM_s(omega_vec(1),r_p_s_vec(1),W_s0,Mm,epsilon_r_w,epsilon_tau_w,epsilon_R0,m,n)
m_p_vec(1) = m_p(r_p_s_vec(1),r_m_s_vec(1),rho)
m_t_vec(1) = m_t(r_p_s_vec(1),r_m_s_vec(1),rho)
RK(1,1) = f_1(cap_delta_mode_p_vec(1),cap_delta_wall_p_vec(1),I_1)
RK(1,2) = f_2(T_theta_VS_s_vec(1),T_theta_EM_s_vec(1),m_p_vec(1))
RK(1,3) = f_3(T_z_VS_s_vec(1),T_z_EM_s_vec(1),m_t_vec(1))

r_p_s_vec(2) = r_p_s(W_s0+dt*RK(1,1)/2.0d0)
r_m_s_vec(2) = r_m_s(W_s0+dt*RK(1,1)/2.0d0)
omega_vec(2) = omega(domega_theta_s0+dt*RK(1,2)/2.0d0,domega_z_s0+dt*RK(1,3)/2.0d0,omega_0,m,n)
cap_delta_mode_p_vec(2) = cap_delta_mode_p(W_s0+dt*RK(1,1)/2.0d0,W_0,cap_delta_mode_0_p)
cap_delta_wall_p_vec(2) = cap_delta_wall_p(omega_vec(2),r_p_s_vec(2),epsilon_r_w,epsilon_tau_w,m)
T_theta_VS_s_vec(2) = T_theta_VS_s(domega_theta_s0+dt*RK(1,2)/2.0d0,r_p_s_vec(2),r_m_s_vec(2),Re,epsilon_delta_D)
T_theta_EM_s_vec(2) = T_theta_EM_s(W_s0+dt*RK(1,1)/2.0d0,omega_vec(2),r_p_s_vec(2),epsilon_r_w,epsilon_tau_w,Mm,m)
T_z_VS_s_vec(2) = T_z_VS_s(domega_z_s0+dt*RK(1,3)/2.0d0,r_m_s_vec(2),epsilon_a,kappa_X,re)
T_z_EM_s_vec(2) = T_z_EM_s(omega_vec(2),r_p_s_vec(2),W_s0+dt*RK(1,1)/2.0d0,Mm,epsilon_r_w,epsilon_tau_w,epsilon_R0,m,n)
m_p_vec(2) = m_p(r_p_s_vec(2),r_m_s_vec(2),rho)
m_t_vec(2) = m_t(r_p_s_vec(2),r_m_s_vec(2),rho)
RK(2,1) = f_1(cap_delta_mode_p_vec(2),cap_delta_wall_p_vec(2),I_1)
RK(2,2) = f_2(T_theta_VS_s_vec(2),T_theta_EM_s_vec(2),m_p_vec(2))
RK(2,3) = f_3(T_z_VS_s_vec(2),T_z_EM_s_vec(2),m_t_vec(2))

r_p_s_vec(3) = r_p_s(W_s0+dt*RK(2,1)/2.0d0)
r_m_s_vec(3) = r_m_s(W_s0+dt*RK(2,1)/2.0d0)
omega_vec(3) = omega(domega_theta_s0+dt*RK(2,2)/2.0d0,domega_z_s0+dt*RK(2,3)/2.0d0,omega_0,m,n)
cap_delta_mode_p_vec(3) = cap_delta_mode_p(W_s0+dt*RK(2,1)/2.0d0,W_0,cap_delta_mode_0_p)
cap_delta_wall_p_vec(3) = cap_delta_wall_p(omega_vec(3),r_p_s_vec(3),epsilon_r_w,epsilon_tau_w,m)
T_theta_VS_s_vec(3) = T_theta_VS_s(domega_theta_s0+dt*RK(2,2)/2.0d0,r_p_s_vec(3),r_m_s_vec(3),Re,epsilon_delta_D)
T_theta_EM_s_vec(3) = T_theta_EM_s(W_s0+dt*RK(2,1)/2.0d0,omega_vec(3),r_p_s_vec(3),epsilon_r_w,epsilon_tau_w,Mm,m)
T_z_VS_s_vec(3) = T_z_VS_s(domega_z_s0+dt*RK(2,3)/2.0d0,r_m_s_vec(3),epsilon_a,kappa_X,re)
T_z_EM_s_vec(3) = T_z_EM_s(omega_vec(3),r_p_s_vec(3),W_s0+dt*RK(2,1)/2.0d0,Mm,epsilon_r_w,epsilon_tau_w,epsilon_R0,m,n)
m_p_vec(3) = m_p(r_p_s_vec(3),r_m_s_vec(3),rho)
m_t_vec(3) = m_t(r_p_s_vec(3),r_m_s_vec(3),rho)
RK(3,1) = f_1(cap_delta_mode_p_vec(3),cap_delta_wall_p_vec(3),I_1)
RK(3,2) = f_2(T_theta_VS_s_vec(3),T_theta_EM_s_vec(3),m_p_vec(3))
RK(3,3) = f_3(T_z_VS_s_vec(3),T_z_EM_s_vec(3),m_t_vec(3))

r_p_s_vec(4) = r_p_s(W_s0+dt*RK(3,1))
r_m_s_vec(4) = r_m_s(W_s0+dt*RK(3,1))
omega_vec(4) = omega(domega_theta_s0+dt*RK(3,2),domega_z_s0+dt*RK(3,3),omega_0,m,n)
cap_delta_mode_p_vec(4) = cap_delta_mode_p(W_s0+dt*RK(3,1),W_0,cap_delta_mode_0_p)
cap_delta_wall_p_vec(4) = cap_delta_wall_p(omega_vec(4),r_p_s_vec(4),epsilon_r_w,epsilon_tau_w,m)
T_theta_VS_s_vec(4) = T_theta_VS_s(domega_theta_s0+dt*RK(3,2),r_p_s_vec(4),r_m_s_vec(4),Re,epsilon_delta_D)
T_theta_EM_s_vec(4) = T_theta_EM_s(W_s0+dt*RK(3,1),omega_vec(4),r_p_s_vec(4),epsilon_r_w,epsilon_tau_w,Mm,m)
T_z_VS_s_vec(4) = T_z_VS_s(domega_z_s0+dt*RK(3,3),r_m_s_vec(4),epsilon_a,kappa_X,re)
T_z_EM_s_vec(4) = T_z_EM_s(omega_vec(4),r_p_s_vec(4),W_s0+dt*RK(3,1),Mm,epsilon_r_w,epsilon_tau_w,epsilon_R0,m,n)
m_p_vec(4) = m_p(r_p_s_vec(4),r_m_s_vec(4),rho)
m_t_vec(4) = m_t(r_p_s_vec(4),r_m_s_vec(4),rho)
RK(4,1) = f_1(cap_delta_mode_p_vec(4),T_theta_EM_s_vec(4),I_1)
RK(4,2) = f_2(T_theta_VS_s_vec(4),T_theta_EM_s_vec(4),m_p_vec(4))
RK(4,3) = f_3(T_z_VS_s_vec(4),T_z_EM_s_vec(4),m_t_vec(4))

W_s = W_s0 +dt*(RK(1,1) +2.0d0*RK(2,1) +2.0d0*RK(3,1) +RK(4,1))/6.0d0
domega_theta_s = domega_theta_s0 +dt*(RK(1,2) +2.0d0*RK(2,2) +2.0d0*RK(3,2) +RK(4,2))/6.0d0
domega_z_s = domega_z_s0 +dt*(RK(1,3) +2.0d0*RK(2,3) +2.0d0*RK(3,3) +RK(4,3))/6.0d0

end subroutine Runge_Kutta_method_of_4th_order

subroutine step_doubling(W_s0,domega_theta_s0,domega_z_s0,W_s,domega_theta_s,domega_z_s, &
                         W_s1,domega_theta_s1,domega_z_s1,W_s2,domega_theta_s2,domega_z_s2, &
                         rho,kappa_X,Re,Mm,epsilon_delta_D,epsilon_tau_w,epsilon_r_w,epsilon_a, &
                         epsilon_R0,I_1,cap_delta_mode_0_p,W_0,omega_0,m,n,i,dt1,dt2)
implicit none
integer(8) :: m, n, i
real(8) :: W_s0, domega_theta_s0, domega_z_s0
real(8) :: W_s, domega_theta_s, domega_z_s
real(8) :: W_s1, domega_theta_s1, domega_z_s1
real(8) :: W_s2, domega_theta_s2, domega_z_s2
real(8) :: rho, kappa_X, Re, Mm, I_1
real(8) :: W_0, omega_0, cap_delta_mode_0_p
real(8) :: epsilon_delta_D, epsilon_tau_w, epsilon_r_w, epsilon_a, epsilon_R0
real(8) :: dt1, dt2

call Runge_Kutta_method_of_4th_order(W_s0,domega_theta_s0,domega_z_s0,W_s,domega_theta_s,domega_z_s, &
rho,kappa_X,Re,Mm,epsilon_delta_D,epsilon_tau_w,epsilon_r_w, &
epsilon_a,epsilon_R0,I_1,dt1,cap_delta_mode_0_p,W_0,omega_0,m,n,i)
W_s1 = W_s
domega_theta_s1 = domega_theta_s
domega_z_s1 = domega_z_s

dt2 = dt1/2.0d0
call Runge_Kutta_method_of_4th_order(W_s0,domega_theta_s0,domega_z_s0,W_s,domega_theta_s,domega_z_s, &
rho,kappa_X,Re,Mm,epsilon_delta_D,epsilon_tau_w,epsilon_r_w, &
epsilon_a,epsilon_R0,I_1,dt2,cap_delta_mode_0_p,W_0,omega_0,m,n,i)
W_s0 = W_s
domega_theta_s0 = domega_theta_s
domega_z_s0 = domega_z_s
call Runge_Kutta_method_of_4th_order(W_s0,domega_theta_s0,domega_z_s0,W_s,domega_theta_s,domega_z_s, &
rho,kappa_X,Re,Mm,epsilon_delta_D,epsilon_tau_w,epsilon_r_w, &
epsilon_a,epsilon_R0,I_1,dt2,cap_delta_mode_0_p,W_0,omega_0,m,n,i)
W_s2 = W_s
domega_theta_s2 = domega_theta_s
domega_z_s2 = domega_z_s

end subroutine step_doubling

subroutine subfunctions(m,n,dt,W_s,domega_theta_s,domega_z_s,W_0,cap_delta_mode_0_p,omega_0, &
                        epsilon_r_w,epsilon_tau_w,epsilon_delta_D,epsilon_a,kappa_X,epsilon_R0, &
                        Re,Mm,rho,subf,t)
implicit none
integer(8),intent(in) :: m, n
real(8) :: t
real(8),intent(in) :: dt, Re, Mm, rho
real(8),intent(in) :: W_s, domega_theta_s, domega_z_s
real(8),intent(in) :: W_0, cap_delta_mode_0_p, omega_0
real(8),intent(in) :: epsilon_r_w, epsilon_tau_w, epsilon_delta_D, epsilon_a,kappa_X, epsilon_R0
real(8) :: r_p_s, r_m_s, omega
real(8) :: cap_delta_mode_p, cap_delta_wall_p
real(8) :: T_theta_VS_s, T_theta_EM_s
real(8) :: T_z_VS_s, T_z_EM_s
real(8) :: m_p, m_t
real(8) :: subf(1:11)

subf(1) = r_p_s(W_s)
subf(2) = r_m_s(W_s)
subf(3) = omega(domega_theta_s,domega_z_s,omega_0,m,n)
subf(4) = cap_delta_mode_p(W_s,W_0,cap_delta_mode_0_p)
subf(5) = cap_delta_wall_p(subf(3),subf(1),epsilon_r_w,epsilon_tau_w,m)
subf(6) = T_theta_VS_s(domega_theta_s,subf(1),subf(2),Re,epsilon_delta_D)
subf(7) = T_theta_EM_s(W_s,subf(3),subf(1),epsilon_r_w,epsilon_tau_w,Mm,m)
subf(8) = T_z_VS_s(domega_z_s,subf(2),epsilon_a,kappa_X,re)
subf(9) = T_z_EM_s(subf(3),subf(1),W_s,Mm,epsilon_r_w,epsilon_tau_w,epsilon_R0,m,n)
subf(10) = m_p(subf(1),subf(2),rho)
subf(11) = m_t(subf(1),subf(2),rho)

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

function omega(domega_theta_s,domega_z_s,omega_0,m,n)
integer(8),intent(in) :: m, n
real(8),intent(in) :: domega_theta_s, domega_z_s
real(8),intent(in) :: omega_0
real(8) :: omega
omega = m*domega_theta_s -n*domega_z_s +omega_0
end function omega

function cap_delta_mode_p(W_s,W_0,cap_delta_mode_0_p)
real(8),intent(in) :: W_s
real(8),intent(in) :: cap_delta_mode_0_p, W_0
real(8) :: cap_delta_mode_p
cap_delta_mode_p = cap_delta_mode_0_p*(1.0d0 -W_s/W_0)
end function cap_delta_mode_p

function cap_delta_wall_p(omega,r_p_s,epsilon_r_w,epsilon_tau_w,m)
integer(8),intent(in) :: m
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
integer(8),intent(in) :: m
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
integer(8),intent(in) :: m, n
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
