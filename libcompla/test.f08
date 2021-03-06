program test
   use compla
   implicit none

   call init_random_seed()
   !call openblas_set_num_threads(1)

   ! Cholesky decomposition
   !call test_chol(100)

   ! Time row/col oriented Cholesky decomp
   !call time_chol(2000)
   !call time_chol(200)
  
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
   !call test_qr(1500,1325)
   !call test_qr(150,132)

   ! Time dgemm from whichever BLAS implementation you linked with
   !call time_dgemm(5000)
   !call time_dgemm(1000)

   ! Test saving data to HDF5 file
   call test_save_stuff(10)
 
  
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
     
         print *,
         print *, "Testing chol:"
         print *, "Number of rows: ",N
 
         wrk = A
         call chol(wrk)
         wrk = A - matmul(transpose(wrk),wrk)
         print *, "1 norm of A-R'*R: ", norm_p(wrk,1)
        
         print *,
         print *, "Testing chol_blas:"
         print *, "Number of rows: ",N
 
         wrk = A
         call chol_blas(wrk)
         wrk = A - matmul(transpose(wrk),wrk)
         print *, "1 norm of A-R'*R: ", norm_p(wrk,1)

         print *,
         print *, "Testing chol_row:"
         print *, "Number of rows: ",N

         wrk = A
         call chol_row(wrk)
         wrk = A - matmul(transpose(wrk),wrk)
         print *, "1 norm of A-R'*R: ", norm_p(wrk,1)
         
         deallocate(A)
      end subroutine test_chol
      ! }}}


      subroutine time_chol(N)
         ! {{{
         integer (kind=4), intent(in) :: N

         real (kind=8), allocatable :: A(:,:), wrk(:,:), wrk2(:,:)
         real (kind=8) :: t_0, t_1

         print *,
         print *, "Timing chol: ",t_1-t_0," CPU seconds"
         print *, "Number of rows: ",N
 
         if (allocated(A)) deallocate(A)
         A = rand_spd_mat(N)

         wrk = A
         call cpu_time(t_0)
         call chol(wrk)
         call cpu_time(t_1)

         !wrk = A - matmul(transpose(wrk),wrk)
         wrk2 = A
         call dgemm('T','N',N,N,N,-1d0,wrk,N,wrk,N,1d0,wrk2,N)
         print *, "1 norm of A-R'*R: ", norm_p(wrk2,1)


         print *,
         print *, "Timing chol_blas: ",t_1-t_0," CPU seconds"
         print *, "Number of rows: ",N
 
         wrk = A
         call cpu_time(t_0)
         call chol_blas(wrk)
         call cpu_time(t_1)
         
         wrk2 = A
         call dgemm('T','N',N,N,N,-1d0,wrk,N,wrk,N,1d0,wrk2,N)
         print *, "1 norm of A-R'*R: ", norm_p(wrk2,1)


         print *,
         print *, "Timing chol_row: ",t_1-t_0," CPU seconds"
         print *, "Number of rows: ",N
 
         wrk = A
         call cpu_time(t_0)
         call chol_row(wrk)
         call cpu_time(t_1)
         wrk2 = A
         call dgemm('T','N',N,N,N,-1d0,wrk,N,wrk,N,1d0,wrk2,N)
         print *, "1 norm of A-R'*R: ", norm_p(wrk2,1)

         deallocate(A,wrk)
      end subroutine time_chol
      ! }}}


      subroutine test_fb_solve_chol(N)
         ! {{{
         integer (kind=4), intent(in) :: N
         real (kind=8), allocatable :: A(:,:), wrk(:,:)
         real (kind=8), allocatable :: b(:,:), x(:,:)

         print *,
         print *, "Testing fb_solve_chol:"
         print *, "Number of rows: ",N
 
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

         print *,
         print *, "Testing fb_solve_blk_chol:"
         print *, "Number of rows: ",N
 
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

         print *,
         print *, "Testing lu:"
         print *, "Number of rows: ",N
 
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

         print *, "1 norm of P*A-L*U: ", norm_p(wrk,1)

      end subroutine
      ! }}}


      subroutine test_lu_nopp(N)
      ! {{{
         integer (kind=4), intent(in) :: N

         real (kind=8), allocatable :: A(:,:)
         real (kind=8), allocatable :: wrk(:,:), L(:,:), U(:,:)

         print *,
         print *, "Testing lu_nopp:"
         print *, "Number of rows: ",N
 
         A = 2d0*rand_mat(N,N)-1d0
         wrk = A

         call lu_nopp(wrk)

         allocate(L(size(A,1),size(A,2)),U(size(A,1),size(A,2)))
         call form_LU(wrk,L,U)

         wrk = matmul(L,U)
         wrk = A - wrk ! A - L*U

         print *, "1 norm of A-L*U: ", norm_p(wrk,1)

      end subroutine test_lu_nopp
      ! }}}


      subroutine test_fb_solve_lu(N)
         ! {{{
         integer (kind=4), intent(in) :: N
         real (kind=8), allocatable :: A(:,:),wrk(:,:)
         real (kind=8), allocatable :: b(:,:), x(:,:)
         integer (kind=4), allocatable :: p(:)

         print *,
         print *, "Testing fb_solve_lu:"
         print *, "Number of rows: ",N
 
         A = rand_mat(N,N)
         wrk = A
         b = rand_mat(N,1)
         x = b ! just allocating x

         call fb_solve_lu(wrk,b,x,p)
         !call print_array(x)
      
         x = matmul(A,x)-b
       
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

         print *,
         print *, "Testing fb_solve_blk_lu:"
         print *, "Number of rows: ",N
 
         A = 2d0*rand_mat(N,N)-1d0
         wrk = A
         b = rand_mat(N,1)
         x = b ! just allocating x
         
         call fb_solve_blk_lu(wrk,b,x,p)
         !call print_array(x)
      
         x = matmul(A,x)-b
       
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

         print *,
         print *, "Testing condest_lu:"
 
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

         print *, "Condition number should be 1197"
         print *, condest_lu(A,wrk,p)

      end subroutine test_condest_lu
      ! }}}

      
      subroutine test_qr(Nr,Nc) 
         ! {{{ 
         integer (kind=4), intent(in) :: Nr,Nc
         real (kind=8), allocatable :: A(:,:), Q(:,:), R(:,:) 
         real (kind=8), allocatable :: wrk(:,:), wrk2(:,:)

         real (kind=8) :: t_0, t_1

         print *,
         print *, "Testing qr:"
         print *, "Size: ",Nr,"x ",Nc
 

         ! Manual test matrix
         ! {{{
         !allocate(A(5,5),x(5,1),b(5))
         !allocate(A(6,5))
         !allocate(Q(5,5),R(5,5))
         !A(:,1) = (/ 1,1,10,1,3 /)
         !A(:,2) = (/ 1,4,3,0,6 /)
         !A(:,3) = (/ 7,7,4,4,8 /)
         !A(:,4) = (/ 1,7,7,3,1 /)
         !A(:,5) = (/ 2,5,3,2,5 /)

         !A(:,1) = (/ 1,1,10,1,3,4 /)
         !A(:,2) = (/ 1,4,3,0,6,2 /)
         !A(:,3) = (/ 7,7,4,4,8,1 /)
         !A(:,4) = (/ 1,7,7,3,1,5 /)
         !A(:,5) = (/ 2,5,3,2,5,1 /)
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

         allocate(Q(Nr,Nr),R(Nr,Nc))
         A = 2d0*rand_mat(Nr,Nc)-1d0
         !x = 2d0*rand_mat(M,1)-1d0
         !x = matmul(A,x)
         !b = x(:,1)
         !b2 = b
         !wrk = A

         wrk = A
         call cpu_time(t_0)
         !call qr(wrk,b2)
         call qr(wrk)
         call cpu_time(t_1)
         !call print_array(wrk)
         !x(:,1) = b2(:)

         call form_qr(wrk,Q,R)
         !call print_array(Q)
         !call print_array(R)
         !wrk = A-matmul(Q,R)
         wrk2 = A
         call dgemm('N','N',Nr,Nc,Nr,-1d0,Q,Nr,R,Nr,1d0,wrk2,Nr)
         print *, "qr time: ",t_1-t_0," CPU seconds"
         print *, "1 norm of A-Q*R: ", norm_p(wrk2,1)

         !b2 = b
         !wrk = A

         print *,
         print *, "Testing qr_blas:"
         print *, "Size: ",Nr,"x ",Nc
 
         wrk = A
         call cpu_time(t_0)
         !call qr_blas(wrk,b)
         call qr_blas(wrk)
         call cpu_time(t_1)
         !call print_array(wrk)
         !x(:,1) = b2(:)

         !!call back_solve_blk(wrk,x)
         !!call print_array(x)

         call form_qr(wrk,Q,R)
         !!wrk = A-matmul(Q,R)
         wrk2 = A
         call dgemm('N','N',Nr,Nc,Nr,-1d0,Q,Nr,R,Nr,1d0,wrk2,Nr)
         print *, "qr_blas time: ",t_1-t_0," CPU seconds"
         print *, "1 norm of A-Q*R: ", norm_p(wrk2,1)


         !call print_array(Q)
         !call print_array(R)
         !call print_array(matmul(Q,R))

   
      end subroutine test_qr
      ! }}}


      subroutine time_dgemm(N)
      ! {{{
         integer (kind=4), intent(in) :: N

         real (kind=8), allocatable :: A(:,:)
         real (kind=8), allocatable :: wrk(:,:)

         real (kind=8) :: t_0, t_1

         print *, 
         print *, "Timing dgemm:"
         print *, "Number of rows (SQUARE): ",N
 
         A = 2d0*rand_mat(N,N)-1d0
         wrk = A

         call cpu_time(t_0)
         call dgemm('N','N',N,N,N,1d0,A,N,A,N,1d0,wrk,N)
         call cpu_time(t_1)

         ! This is for ATLAS LAPACK v. openBLAS LAPACK v. etc
         print *, "dgemm time: ",t_1-t_0," CPU seconds"
   
      end subroutine time_dgemm
      ! }}}


      subroutine save_stuff(savefile,vec,array)
      ! {{{
         character (len=*) :: savefile
         real (kind=8) :: vec(:),array(:,:)

         integer (kind=intk) :: h5error, file_id
   
         call h5open_f(h5error)
         call h5fcreate_f(savefile, H5F_ACC_TRUNC_F, file_id, h5error)

         call write_dset(file_id, vec, "vec")
         call write_dset(file_id, array, "array")

         call h5fclose_f(file_id,h5error)
         call h5close_f(h5error)

      end subroutine save_stuff
      ! }}}

      subroutine test_save_stuff(N)
      ! {{{
         integer (kind=4), intent(in) :: N

         real (kind=8), allocatable :: vec(:),array(:,:)

         print *,
         print *, "Testing HDF5:"
 
         array = rand_mat(N,N)
         vec = array(:,1)

         call save_stuff("save_stuff.h5",vec,array)

         print *, "Executing: octave --quiet --eval ""load(''save_stuff.h5'');array,transpose(vec)"""
         ! /usr/bin/sh has the following quotes (and more)
         !  ' - strong quote: print with '' in Fortran
         !  " - weak quote: print with "" in Fortran
         call execute_command_line("''octave --quiet --eval ""load('save_stuff.h5');array,transpose(vec)""")

      end subroutine test_save_stuff
      ! }}}


end program test

! vim: set ts=3 sw=3 sts=3 et :
! vim: foldmarker={{{,}}}
