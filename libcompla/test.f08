program test
   use compla
   implicit none

   call init_random_seed()

   ! Cholesky decomposition
   !call test_chol(100)
  
   ! Forward and back solve with Cholesky decomp
   !call test_fb_solve_chol(100)

   ! Forward and back solve by blocks with Cholesky decomp (not by blocks)
   !call test_fb_solve_blk_chol(100)

   ! LU decomposition with partial pivoting
   !call test_lu(100)

   ! forward and back solve (by blocks) with LU decomp (not by blocks)
   call test_fb_solve_blk_lu(1000)
  
   ! {{{
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
   ! }}}

   contains

      subroutine test_chol(N)
         ! {{{
         integer (kind=4), intent(in) :: N

         real (kind=8), allocatable :: A(:,:), wrk(:,:)

         if (allocated(A)) deallocate(A)
         A = rand_spd_mat(N)

         !allocate(A(5,5))
         !A(:,1) = (/ 27,10,5,12,4 /)
         !A(:,2) = (/ 10,29,11,15,14 /)
         !A(:,3) = (/ 5,11,18,14,10 /)
         !A(:,4) = (/ 12,15,14,27,11 /)
         !A(:,5) = (/ 4,14,10,11,20 /)
      
         !A(:,1) = (/ 1,2,3,4,5 /)
         !A(:,2) = (/ 2,8,14,20,26 /)
         !A(:,3) = (/ 3,14,34,54,74 /)
         !A(:,4) = (/ 4,20,54,104,154 /)
         !A(:,5) = (/ 5,26,74,154,259 /)
      
         !allocate(A(3,3))
         !A(:,1) = (/ 9,6,9 /)
         !A(:,2) = (/ 6,5,8 /)
         !A(:,3) = (/ 9,8,17 /)
      
         !print *, "Before chol():"
         !call print_array(A)
         !print *,
         !call chol(A)
         !print *, "After chol():"
         !call print_array(A)
     
         wrk = A
         call chol(wrk)
         wrk = A - matmul(transpose(wrk),wrk)

         print *,
         print *, "Testing chol:"
         print *, "Number of rows: ",N
         print *, "Fro norm of A-R'*R: ", norm_f(wrk)
         
         deallocate(A)
      end subroutine test_chol
      ! }}}

      subroutine test_fb_solve_chol(N)
         ! {{{
         integer (kind=4), intent(in) :: N
         real (kind=8), allocatable :: A(:,:), wrk(:,:)
         real (kind=8), allocatable :: b(:,:), x(:,:)

         if (allocated(A)) deallocate(A)
         if (allocated(b)) deallocate(b)
         if (allocated(x)) deallocate(x)
         if (allocated(wrk)) deallocate(wrk)
         A = rand_spd_mat(N)
         b = rand_mat(N,1)
         x = b

         !allocate(A(5,5))
         !allocate(b(5,1))
         !allocate(x(5,1))
         !allocate(wrk(size(A,1),size(A,2)))
         !A(:,1) = (/ 27,10,5,12,4 /)
         !A(:,2) = (/ 10,29,11,15,14 /)
         !A(:,3) = (/ 5,11,18,14,10 /)
         !A(:,4) = (/ 12,15,14,27,11 /)
         !A(:,5) = (/ 4,14,10,11,20 /)
         !
         !b(:,1) = (/ 1,2,3,4,5 /)

         !A(:,1) = (/ 1,2,3,4,5 /)
         !A(:,2) = (/ 2,8,14,20,26 /)
         !A(:,3) = (/ 3,14,34,54,74 /)
         !A(:,4) = (/ 4,20,54,104,154 /)
         !A(:,5) = (/ 5,26,74,154,259 /)
        
         ! manual test of for_solve(), back_solve()
         ! {{{
         !if (allocated(R)) deallocate(R)
         !if (allocated(Rt)) deallocate(R)
         !allocate(R(size(A,1),size(A,2)))
         !allocate(Rt(size(A,1),size(A,2)))
         !R = A
         !call chol(R)
         !Rt = transpose(R)
      
         !call for_solve(Rt,b)
         !call print_array(b)
         !
         !print *,
         !call back_solve(R,b)
         !call print_array(b)
         !deallocate(A,R,Rt) 
         ! }}}
      
         call fb_solve_chol(A,b,x)
         !call print_array(x)
      
         x = matmul(A,x)-b
       
         print *,
         print *, "Testing fb_solve_chol:"
         print *, "Number of rows: ",N
         print *, "2-norm of residual vector = ", norm_f(x)
         !call print_array(x)

         deallocate(A,b,x)

      end subroutine test_fb_solve_chol
      ! }}}

      subroutine test_fb_solve_blk_chol(N)
         ! {{{
         integer (kind=4), intent(in) :: N
         real (kind=8), allocatable :: A(:,:)
         real (kind=8), allocatable :: b(:,:), x(:,:)

         !real (kind=8), allocatable :: R(:,:), Rt(:,:)

         A = rand_spd_mat(N)
         b = rand_mat(N,1)
         x = b ! just allocating x

         ! manual test of for_solve_blk(), back_solve_blk()
         ! {{{
         !if (allocated(R)) deallocate(R)
         !if (allocated(Rt)) deallocate(R)
         !allocate(R(size(A,1),size(A,2)))
         !allocate(Rt(size(A,1),size(A,2)))
         !R = A
         !call chol(R)
         !Rt = transpose(R)
     
         !!x=b
         !!call for_solve_blk(Rt,x)
         !!x = matmul(Rt,x)-b

         !!print *,
         !!print *, "Testing for_solve_blk:"
         !!print *, "Number of rows: ",N
         !!print *, "2-norm of residual vector = ", norm_f(x)
         !!!call print_array(x)

         !x=b
         !call back_solve_blk(R,x)
         !x = matmul(R,x)-b

         !print *,
         !print *, "Testing back_solve_blk:"
         !print *, "Number of rows: ",N
         !print *, "2-norm of residual vector = ", norm_f(x)
         !!call print_array(x)

         !deallocate(A,R,Rt,b,x) 
         ! }}}
      
         call fb_solve_blk_chol(A,b,x)
         !call print_array(x)
      
         x = matmul(A,x)-b
       
         print *,
         print *, "Testing fb_solve_blk_chol:"
         print *, "Number of rows: ",N
         print *, "2-norm of residual vector = ", norm_f(x)
         !call print_array(x)

         deallocate(A,b,x)

      end subroutine test_fb_solve_blk_chol
      ! }}}

      subroutine test_lu(N)
      ! {{{
         integer (kind=4), intent(in) :: N

         real (kind=8), allocatable :: A(:,:)
         integer (kind=4), allocatable :: p(:)
         real (kind=8), allocatable :: wrk(:,:), L(:,:), U(:,:)
         integer (kind=4) :: i

         ! Manual test matrix
         ! {{{
         !allocate(A(5,5))
         !A(:,1) = (/ 1,1,10,1,3 /)
         !A(:,2) = (/ 1,4,3,0,6 /)
         !A(:,3) = (/ 7,7,4,4,8 /)
         !A(:,4) = (/ 1,7,7,3,1 /)
         !A(:,5) = (/ 2,5,3,2,5 /)
         !allocate(p(5))
         !p(:) = (/ (i,i=1,5) /)

         !allocate(A(3,3))
         !A(:,1) = (/ 2,4,8 /)
         !A(:,2) = (/ 1,3,7 /)
         !A(:,3) = (/ 1,3,9 /)
         !allocate(p(3))
         !p(:) = (/ (i,i=1,3) /)
         ! }}}

         A = rand_mat(N,N)
         allocate(p(N))
         p(:) = (/ (i,i=1,N) /)

         wrk = A

         call lu(wrk,p)
         call apply_perm_vector(wrk,p,1)

         allocate(L(size(A,1),size(A,2)),U(size(A,1),size(A,2)))
         call form_LU(wrk,L,U)

         wrk = matmul(L,U)
         call apply_perm_vector(A,p,1)
         wrk = A - wrk ! P*A - L*U

         print *,
         print *, "Testing lu:"
         print *, "Number of rows: ",N
         print *, "Fro norm of P*A-L*U: ", norm_f(wrk)
 

      end subroutine
      ! }}}

      subroutine test_fb_solve_blk_lu(N)
         ! {{{
         integer (kind=4), intent(in) :: N
         real (kind=8), allocatable :: A(:,:)
         real (kind=8), allocatable :: b(:,:), x(:,:)

         !real (kind=8), allocatable :: R(:,:), Rt(:,:)

         A = rand_mat(N,N)
         b = rand_mat(N,1)
         x = b ! just allocating x

         
         call fb_solve_blk_lu(A,b,x)
         !call print_array(x)
      
         x = matmul(A,x)-b
       
         print *,
         print *, "Testing fb_solve_blk_lu:"
         print *, "Number of rows: ",N
         print *, "2-norm of residual vector = ", norm_f(x)
         !call print_array(x)

         deallocate(A,b,x)

      end subroutine test_fb_solve_blk_lu
      ! }}}


end program test
