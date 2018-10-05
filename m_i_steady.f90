
program isseiuenaga
implicit none
integer :: i, j, k
integer,parameter :: n = 1000, m = 10000
real(8) :: omega0_min, omega0_max, domega0
real(8) :: omega_min, omega_max, domega
real(8) :: epsilon, Gamma, W0, beta
real(8) :: omega
real(8) :: omega0_vec(0:n), omega_old_vec(0:n), omega_mat(0:n,0:n)

open(11,file="result-m_i_steady-beta=1_Gamma=0_W0=16.dat",status="replace")
!open(12,file="result-m_i_steady-test.dat",status="replace")

omega0_min = 0.0d0
omega0_max = 10.0d0
domega0 = (omega0_max -omega0_min)/dble(n)
omega_min = 0.0d0
omega_max = 10.0d0
domega = (omega_max -omega_min)/dble(n)
epsilon = 0.001d0
beta = 1.0d0
Gamma = 0.0d0
W0 = 2.0d0**4

do i = 0,n
   omega0_vec(i) = dble(i)*domega0
   omega_old_vec(i) = dble(i)*domega
end do

omega_mat(0,0:n) = 0.0d0

do i = 1,n
   if(mod(i,n/10) == 0) then
      write(*,*) "step :", i, "omega0 =", omega0_vec(i)
   end if
   do j = 1,n
      omega = omega_old_vec(j)
      do k = 1,m
         call Newtons_method(omega_mat,omega,omega0_vec,Gamma,W0,beta,i,j,n)
!         write(12,*) omega_mat(i,j), abs(omega_mat(i,j)/omega -1.0d0), epsilon
         if(abs((omega_mat(i,j)-omega)/omega_mat(i,j)) < epsilon) then
            exit
         elseif(k == m) then
                omega_mat(i,j) = 0.0d0
         end if
         omega = omega_mat(i,j)
      end do
!      write(12,*) ""
   end do
end do

write(*,*) ""
write(*,*) "Now writing ..."
write(*,*) ""

do i = 1,n
   do j = 1,n
      write(11,*) omega0_vec(i), omega_old_vec(j), omega_mat(i,j)
   end do
   if(mod(i,n/10) == 0) then
      write(*,*) i*100/n, "%"
   end if
end do

close(11)
!close(12)

end program isseiuenaga


subroutine Newtons_method(omega_mat,omega,omega0_vec,Gamma,W0,beta,i,j,n)
implicit none
integer,intent(in) :: i, j, n
real(8) :: omega_mat(0:n,0:n)
real(8) :: omega, omega0
real(8),intent(in) :: omega0_vec(0:n)
real(8),intent(in) :: Gamma, W0, beta
real(8) :: p, dp
omega0 = omega0_vec(i)
omega_mat(i,j) = omega -p(omega,omega0,Gamma,W0,beta)/dp(omega,omega0,Gamma,W0,beta)
end subroutine Newtons_method


function p(omega,omega0,Gamma,W0,beta)
real(8),intent(in) :: omega, omega0
real(8),intent(in) :: Gamma, W0, beta
real(8) :: p
p = (omega +Gamma -omega0)*(1.0d0 +3.0d0*omega**2)**5 + &
     8.0d0*W0**4*omega*(1.0d0 +3.0d0*beta*omega**2)**4
end function p

function dp(omega,omega0,Gamma,W0,beta)
real(8),intent(in) :: omega, omega0
real(8),intent(in) :: Gamma, W0, beta
real(8) :: dp
dp = (33.0d0*omega**2 +30.0d0*(Gamma -omega0)*omega +1.0d0)*(1.0d0 +3.0d0*omega**2)**4 + &
      8.0d0*W0**4*(27.0d0*beta*omega**2 +1.0d0)*(1.0d0 +3.0d0*beta*omega**2)**3
end function dp
