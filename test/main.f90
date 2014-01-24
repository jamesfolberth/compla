program main
   use compla

   implicit none

   integer (kind=4), parameter :: Nr=3,Nc=5
   integer (kind=4) :: i,j
   real (kind=8), allocatable :: A(:,:)

   integer (kind=4) :: clock

   call tic(clock)

   allocate(A(Nr,Nc))
   row: do i=1,size(A,1)
      col: do j=1,size(A,2)
         A(i,j) = 1. / (real(i)+real(j))
      end do col
   end do row

   A = -A

   call print_array(A)

   deallocate(A)

   call toc(clock)

end program main
