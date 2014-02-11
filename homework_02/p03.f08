program p03
   use compla
   implicit none

   call init_random_seed()

   print *, "Comparing lu,lu_nopp:"
   call compare_pivot(40)
   print *,
   call compare_pivot(80)
   print *,
   call compare_pivot(160)
   print *,

   print *,
   print *, "Comparing lu,lu_nopp (small a_{1,1})"
   call compare_pivot_small_a11(40)
   print *,
   call compare_pivot_small_a11(80)
   print *,
   call compare_pivot_small_a11(160)

   contains

      subroutine compare_pivot(N)
      ! {{{
         integer (kind=4), intent(in) :: N

         real (kind=8), allocatable :: A(:,:)
         integer (kind=4), allocatable :: p(:)
         real (kind=8), allocatable :: wrk(:,:), L(:,:), U(:,:)
         integer (kind=4) :: i

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

         print *, "Comparing lu,lu_nopp:"
         print *, "Number of rows: ",N
         print *, "Inf norm of P*A-L*U (pp): ", norm_p(wrk,0)

         wrk = A
         call lu_nopp(wrk)
         call form_LU(wrk,L,U)
         wrk = matmul(L,U)
         wrk = A-wrk ! A - L*U

         print *, "Inf norm of L (nopp): ", norm_p(L,0)
         print *, "Inf norm of U (nopp): ", norm_p(U,0)
         print *, "Inf norm of A-L*U (nopp): ", norm_p(wrk,0)

         deallocate(A,wrk,L,U)
      end subroutine compare_pivot
      ! }}}

      subroutine compare_pivot_small_a11(N)
      ! {{{
         integer (kind=4), intent(in) :: N

         real (kind=8), allocatable :: A(:,:)
         integer (kind=4), allocatable :: p(:)
         real (kind=8), allocatable :: wrk(:,:), L(:,:), U(:,:)
         integer (kind=4) :: i

         A = 2d0*rand_mat(N,N)-1d0
         A(1,1)=0.0001d0*A(1,1)
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

         print *, "Number of rows: ",N
         print *, "Inf norm of L (pp): ", norm_p(L,0)
         print *, "Inf norm of U (pp): ", norm_p(U,0)
         print *, "Inf norm of P*A-L*U (pp): ", norm_p(wrk,0)

         wrk = A
         call lu_nopp(wrk)
         call form_LU(wrk,L,U)
         wrk = matmul(L,U)
         wrk = A-wrk ! A - L*U

         print *, "Inf norm of L (nopp): ", norm_p(L,0)
         print *, "Inf norm of U (nopp): ", norm_p(U,0)
         print *, "Inf norm of A-L*U (nopp): ", norm_p(wrk,0)

         deallocate(A,wrk,L,U)
      end subroutine compare_pivot_small_a11
      ! }}}

end program p03

! vim: set ts=3 sw=3 sts=3 et :
! vim: foldmarker={{{,}}}
