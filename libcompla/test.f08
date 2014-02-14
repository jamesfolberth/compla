program test
   use compla
   implicit none

   call init_random_seed()

   ! Cholesky decomposition
   !call test_chol(100)

   ! Time row/col oriented Cholesky decomp
   !call time_chol(1000)
  
   ! Forward and back solve with Cholesky decomp
   !call test_fb_solve_chol(100)

   ! Forward and back solve by blocks with Cholesky decomp (not by blocks)
   !call test_fb_solve_blk_chol(100)

   ! LU decomposition with partial pivoting
   !call test_lu(100)

   ! LU decomposition WITHOUT partial pivoting
   !call test_lu_nopp(100)

   ! forward and back solve (not by blocks) with LU decomp (not by blocks)
   !call test_fb_solve_lu(100)
   
   ! forward and back solve (by blocks) with LU decomp (not by blocks)
   !call test_fb_solve_blk_lu(100)

   ! Test matrix condition number
   !call test_condest_lu()

   ! Test QR decomp by reflectors
   call test_qr(100)
 
  
   ! {{{
   !print *,"test: tic"
   !call tic(t)

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

         !! manual test
         !! {{{
         !!allocate(A(5,5))
         !!A(:,1) = (/ 27,10,5,12,4 /)
         !!A(:,2) = (/ 10,29,11,15,14 /)
         !!A(:,3) = (/ 5,11,18,14,10 /)
         !!A(:,4) = (/ 12,15,14,27,11 /)
         !!A(:,5) = (/ 4,14,10,11,20 /)
     
         !allocate(A(5,5))
         !A(:,1) = (/ 1,2,3,4,5 /)
         !A(:,2) = (/ 2,8,14,20,26 /)
         !A(:,3) = (/ 3,14,34,54,74 /)
         !A(:,4) = (/ 4,20,54,104,154 /)
         !A(:,5) = (/ 5,26,74,154,259 /)
      
         !!allocate(A(3,3))
         !!A(:,1) = (/ 9,6,9 /)
         !!A(:,2) = (/ 6,5,8 /)
         !!A(:,3) = (/ 9,8,17 /)
     
         !wrk = A
         !print *, "Before chol():"
         !call print_array(wrk)
         !print *,
         !!call chol(A)
         !call chol_row(wrk)
         !print *, "After chol():"
         !call print_array(wrk)

         !print *, norm_f(A-matmul(transpose(wrk),wrk))
         ! }}}
     
         wrk = A
         call chol(wrk)
         wrk = A - matmul(transpose(wrk),wrk)
         print *,
         print *, "Testing chol:"
         print *, "Number of rows: ",N
         print *, "1 norm of A-R'*R: ", norm_p(wrk,1)
         
         deallocate(A)
      end subroutine test_chol
      ! }}}


      subroutine time_chol(N)
         ! {{{
         integer (kind=4), intent(in) :: N

         real (kind=8), allocatable :: A(:,:), wrk(:,:)
         real (kind=8) :: t_0, t_1

         if (allocated(A)) deallocate(A)
         A = rand_spd_mat(N)

         wrk = A
         call cpu_time(t_0)
         call chol(wrk)
         call cpu_time(t_1)
         wrk = A - matmul(transpose(wrk),wrk)
         print *,
         print *, "Timing chol: ",t_1-t_0," CPU seconds"
         print *, "Number of rows: ",N
         print *, "1 norm of A-R'*R: ", norm_p(wrk,1)

         wrk = A
         call cpu_time(t_0)
         call chol_row(wrk)
         call cpu_time(t_1)
         wrk = A - matmul(transpose(wrk),wrk)
         print *,
         print *, "Timing chol_row: ",t_1-t_0," CPU seconds"
         print *, "Number of rows: ",N
         print *, "1 norm of A-R'*R: ", norm_p(wrk,1)

         deallocate(A,wrk)
      end subroutine time_chol
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

         ! manual test of for_solve(), back_solve()
         ! {{{
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
         print *, "2-norm of residual vector = ", norm_p(x,2)
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
         print *, "2-norm of residual vector = ", norm_p(x,2)
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

         A = 2d0*rand_mat(N,N)-1d0
         allocate(p(N))
         p(:) = (/ (i,i=1,N) /)

         wrk = A

         call lu(wrk,p)
         call apply_perm_vector(wrk,p,0)

         allocate(L(size(A,1),size(A,2)),U(size(A,1),size(A,2)))
         call form_LU(wrk,L,U)

         wrk = matmul(L,U)
         call apply_perm_vector(wrk,p,1)
         wrk = A - wrk ! P*A - L*U

         print *,
         print *, "Testing lu:"
         print *, "Number of rows: ",N
         print *, "1 norm of P*A-L*U: ", norm_p(wrk,1)

      end subroutine
      ! }}}


      subroutine test_lu_nopp(N)
      ! {{{
         integer (kind=4), intent(in) :: N

         real (kind=8), allocatable :: A(:,:)
         real (kind=8), allocatable :: wrk(:,:), L(:,:), U(:,:)

         A = 2d0*rand_mat(N,N)-1d0
         wrk = A

         call lu_nopp(wrk)

         allocate(L(size(A,1),size(A,2)),U(size(A,1),size(A,2)))
         call form_LU(wrk,L,U)

         wrk = matmul(L,U)
         wrk = A - wrk ! A - L*U

         print *,
         print *, "Testing lu_nopp:"
         print *, "Number of rows: ",N
         print *, "1 norm of A-L*U: ", norm_p(wrk,1)

      end subroutine test_lu_nopp
      ! }}}


      subroutine test_fb_solve_lu(N)
         ! {{{
         integer (kind=4), intent(in) :: N
         real (kind=8), allocatable :: A(:,:),wrk(:,:)
         real (kind=8), allocatable :: b(:,:), x(:,:)
         integer (kind=4), allocatable :: p(:)

         A = rand_mat(N,N)
         wrk = A
         b = rand_mat(N,1)
         x = b ! just allocating x

         call fb_solve_lu(wrk,b,x,p)
         !call print_array(x)
      
         x = matmul(A,x)-b
       
         print *,
         print *, "Testing fb_solve_lu:"
         print *, "Number of rows: ",N
         print *, "2norm of residual vector = ", norm_p(x,2)
         print *, "1 norm condition number = ", condest_lu(A,wrk,p)
         print *, "1 norm relative error = ", condest_lu(A,wrk,p)*norm_p(x,1)/norm_p(b,1)
 
         !call print_array(x)

         deallocate(A,b,x)

      end subroutine test_fb_solve_lu
      ! }}}


      subroutine test_fb_solve_blk_lu(N)
         ! {{{
         integer (kind=4), intent(in) :: N
         real (kind=8), allocatable :: A(:,:), wrk(:,:)
         real (kind=8), allocatable :: b(:,:), x(:,:)
         integer (kind=4), allocatable :: p(:)

         !real (kind=8), allocatable :: R(:,:), Rt(:,:)

         A = 2d0*rand_mat(N,N)-1d0
         wrk = A
         b = rand_mat(N,1)
         x = b ! just allocating x
         
         call fb_solve_blk_lu(wrk,b,x,p)
         !call print_array(x)
      
         x = matmul(A,x)-b
       
         print *,
         print *, "Testing fb_solve_blk_lu:"
         print *, "Number of rows: ",N
         print *, "2 norm of residual vector = ", norm_p(x,2)
         print *, "1 norm condition number = ", condest_lu(A,wrk,p)
         print *, "1 norm relative error = ", condest_lu(A,wrk,p)*norm_p(x,1)/norm_p(b,1)
         !call print_array(x)

         deallocate(A,b,x)

      end subroutine test_fb_solve_blk_lu
      ! }}}


      subroutine test_condest_lu()
         ! {{{
         !integer (kind=4), intent(in) :: N

         real (kind=8), allocatable :: A(:,:),wrk(:,:)
         integer (kind=4), allocatable :: p(:)
         !integer (kind=4) :: i

         ! test matrix
         ! K(A) = 1197.
         allocate(A(2,2))
         A(:,1) = (/ 1d0, -0.99d0 /)
         A(:,2) = (/ -2d0, 1.99d0 /)
         allocate(p(2))
         p = (/ 1,2 /)
         wrk = A
         call lu(wrk,p)
         call apply_perm_vector(wrk,p,0)

         !A = 2d0*rand_mat(N,N)-1d0
         !wrk = A
         !allocate(p(N))
         !p = (/ (i,i=1,N) /)
         !call lu(wrk,p)

         print *,
         print *, "Testing condest_lu:"
         print *, "Condition number should be 1197"
         print *, condest_lu(A,wrk,p)

      end subroutine test_condest_lu
      ! }}}


      subroutine test_qr(N)
      ! {{{
         integer (kind=4), intent(in) :: N

         real (kind=8), allocatable :: A(:,:), Q(:,:), R(:,:)
         real (kind=8), allocatable :: wrk(:,:), x(:,:), b(:)

         ! Manual test matrix
         ! {{{
         !allocate(A(5,5),x(5,1),b(5))
         !allocate(Q(5,5),R(5,5))
         !A(:,1) = (/ 1,1,10,1,3 /)
         !A(:,2) = (/ 1,4,3,0,6 /)
         !A(:,3) = (/ 7,7,4,4,8 /)
         !A(:,4) = (/ 1,7,7,3,1 /)
         !A(:,5) = (/ 2,5,3,2,5 /)
         !x(:,1) = (/ 1,2,3,4,5 /)
         !x = matmul(A,x)
         !b(:) = x(:,1)

         ! octave's qr
         !  -10.583005   -5.008029   -7.748272   -7.937254   -5.102520
         !    0.094491   -6.076154   -9.248685   -2.674389   -5.833695
         !    0.944911   -0.359921    6.958888    2.552910    2.228844
         !    0.094491   -0.085365   -0.508625   -5.685973   -1.267749
         !    0.283473    0.731371    0.399906   -0.910104    0.597787

         ! octave's doesn't match mine with QR overwritten on A, but when I/octave form Q,R, they're the same
         ! }}}

         A = 2d0*rand_mat(N,N)-1d0
         x = 2d0*rand_mat(N,1)-1d0
         x = matmul(A,x)
         b = x(:,1)
         wrk = A

         !call print_array(A)
         wrk = A
         call qr(wrk,b)
         x(:,1) = b(:)
         !call print_array(wrk)

         call back_solve_blk(wrk,x)
         !call print_array(x)

         allocate(Q(N,N),R(N,N))
         call form_qr(wrk,Q,R)
         wrk = A-matmul(Q,R)

         print *,
         print *, "Testing qr:"
         print *, "Number of rows (SQUARE): ",N
         print *, "1 norm of A-Q*R: ", norm_p(wrk,1)



         !call print_array(Q)
         !call print_array(R)
         !call print_array(matmul(Q,R))

   
      end subroutine test_qr
      ! }}}


end program test

! vim: set ts=3 sw=3 sts=3 et :
! vim: foldmarker={{{,}}}
