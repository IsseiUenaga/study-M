
program isseiuenaga
implicit none
integer,parameter :: n = 1000
integer :: i, j
real(8) :: da, a(0:n), amin, amax, x(0:n,0:n)

open(11,file="result-l_m.dat",status="replace")

amin = 0.0d0
amax = 4.0d0
da = (amax-amin)/dble(n)
x(0:n,0) = 0.6d0

do i = 0,n
   a(i) = dble(i)*da
   if(mod(i,n/10) == 0) then
      write(*,*) "step :", i, "a =", a(i)
   end if
   do j = 0,n-1
      x(i,j+1) = a(i)*x(i,j)*(1.0d0-x(i,j))
   end do
end do

write(*,*) "Now writing ..."
do i = 0,n
   do j = 950,n-1
      write(11,*) a(i), x(i,j), x(i,j+1)
   end do
   write(11,*) ""
   if(mod(i,n/10) == 0) then
      write(*,*) i*100/n, "%"
   end if
end do

close(11)

end program isseiuenaga
