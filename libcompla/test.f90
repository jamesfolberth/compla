program test
   use compla
   implicit none
   
   !integer (kind=4) :: t,i,j
   integer (kind=4), parameter :: Nr=3,Nc=5
   real (kind=8), allocatable :: A(:,:)

   ! Cholesky Decomposition
   ! Test chol(), chol_blk(), ...
   ! {{{
   if (allocated(A)) deallocate(A)
   allocate(A(5,5))
   A(:,1) = (/ 27,10,5,12,4 /)
   A(:,2) = (/ 10,29,11,15,14 /)
   A(:,3) = (/ 5,11,18,14,10 /)
   A(:,4) = (/ 12,15,14,27,11 /)
   A(:,5) = (/ 4,14,10,11,20 /)

   !A(:,1) = (/ 1,2,3,4,5 /)
   !A(:,2) = (/ 2,8,14,20,26 /)
   !A(:,3) = (/ 3,14,34,54,74 /)
   !A(:,4) = (/ 4,20,54,104,154 /)
   !A(:,5) = (/ 5,26,74,154,259 /)

   !allocate(A(3,3))
   !A(:,1) = (/ 9,6,9 /)
   !A(:,2) = (/ 6,5,8 /)
   !A(:,3) = (/ 9,8,17 /)

   print *, "Before chol():"
   call print_array(A)
   print *,
   call chol(A)
   print *, "After chol():"
   call print_array(A)

   print *,
   print *, "R'*R:"
   A = matmul(transpose(A),A)
   call print_array(A)
   
   deallocate(A)
   ! }}}


   !print *,"test: tic"
   !call tic(t)

   !allocate(A(Nr,Nc))
   !
   !row: do i=1,size(A,1)
   !   col: do j=1,size(A,2)
   !      A(i,j) = 1./real(i+j)
   !   end do col
   !end do row

   !print *,"test: print_array"
   !call print_array(A)

   !deallocate(A)

   !print *,"test: toc"
   !call toc(t)

end program test
