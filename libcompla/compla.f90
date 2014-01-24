! libcompla.f90
! Computational Linear Algebra library
! James Folberth - Spring 2014

! wrap lib in module so it interfaces properly
module compla

   real (kind=8), parameter :: ZERO_TOL = 10**(-10)

   contains

   !!!!!!!!!!!!!!!!!!!
   ! Cholesky Decomp !
   !!!!!!!!!!!!!!!!!!!
   ! {{{ 
   ! These both use the outer product formulation of the Cholesky decomposition
   ! Input matrix must by symmetric, positive definite.
   ! TODO optional input for transpose output or not

   ! Outer product formulation
   subroutine chol(A)
   real (kind=8), allocatable :: A(:,:)

   integer (kind=4) :: i,j,k,Nc
   real (kind=8), allocatable :: b(:,:)

   ! check that A is square
   if (size(A,1) /= size(A,2)) then
      print *, "error: compla.f90: chol: input matrix is not square"
      stop
   end if

   ! parameter to store size(A,1)
   Nc = size(A,1)
   allocate(b(Nc,1))

   ! main Cholesky loop
   recurse: do i=1,Nc-1
      ! check a_{i,i} > 0
      if (A(i,i) <= 0d0) then
         print *, "error: compla.f90: chol: input matrix is not positive definite"
         stop
      end if

      A(i,i) = sqrt(A(i,i))
      b(i+1:Nc,1) = A(i+1:Nc,i) / A(i,i)
      A(i+1:Nc,i) = b(i+1:Nc,1)
     
      ! Outer product A~ <- A^ - s*s'
      row: do j=i+1,Nc
         col: do k=i+1,Nc
            A(j,k) = A(j,k)-b(j,1)*b(k,1)
         end do col
      end do row

   end do recurse

   if (A(Nc,Nc) <= 0d0) then
      print *, "error: compla.f90: chol: input matrix is not positive definite"
      stop
   end if
  
   deallocate(b)
   A(Nc,Nc) = sqrt(A(Nc,Nc))

   ! TODO smarter way?
   A = transpose(A)
   stupid: do i=1,Nc
      dumb: do j=1,i-1
         A(i,j) = 0d0
      end do dumb
   end do stupid
         
   end subroutine chol 

   !subroutine chol_blk()

   !function sym_outer_prod()

   !subroutine check_sym()

   ! }}}


   !!!!!!!!!!!!!!!!!!
   ! For/Back Solve !
   !!!!!!!!!!!!!!!!!!
   ! {{{
   subroutine back_solve(U,b)
      real (kind=8), allocatable :: U(:,:), b(:,:)
      integer (kind=4) :: i,j,Nc

      ! call check_square(U)
      ! call check_upper_tri(U)

      Nc=size(U,1)

      row: do j=Nc,1,-1
         ! zero on diagonal
         if (abs(U(j,j)) <= ZERO_TOL) then
            print *, "error: compla.f90: back_solve: input matrix is singular to tolerance"
            stop
         end if

         b(j,1) = b(j,1)/U(j,j)
         col: do i=j-1,1,-1
            b(i,1)=b(i,1) - U(i,j)*b(j,1)
         end do col
      end do row

   end subroutine back_solve

   subroutine for_solve(L,b)
      real (kind=8), allocatable :: L(:,:), b(:,:)
      integer (kind=4) :: i,j,Nc

      ! call check_square(L)
      ! call check_lower_tri(L)

      Nc = size(L,1)

      row: do j=1,Nc
         ! zero on diagonal means singular matrix
         if (abs(L(j,j)) <= ZERO_TOL) then
            print *, "error: compla.f90: for_solve: input matrix is singular to tolerance"
            stop
         end if

         b(j,1) = b(j,1)/L(j,j)
         col: do i=j+1,Nc
            b(i,1)=b(i,1) - L(i,j)*b(j,1) 
         end do col
      end do row

   end subroutine for_solve
   
   subroutine fb_solve_chol(A,b,x)
      real (kind=8), allocatable :: A(:,:), b(:,:), x(:,:)
      real (kind=8), allocatable :: wrk(:,:)
      integer (kind=4) :: Nr,Nc

      Nr = size(A,1)
      Nc = size(A,2)
   
      if (Nr /= Nc) then
         print *, "error: compla.f90: fb_solve_chol: input matrix is not square"
         stop
      end if

      allocate(wrk(Nc,Nc))
      x = b
      wrk = A
      call chol(wrk) ! stores R in wrk
      wrk = transpose(wrk)
      call for_solve(wrk, x) !A*x = R'*R*x=b -> R'*y=b
      wrk = transpose(wrk)
      call back_solve(wrk,x) ! Rx=y

   end subroutine fb_solve_chol

   ! call back_solve_block()
   ! call forward_solve_block()
   ! call fbsolve_block()

   ! }}}


   !!!!!!!!!!!!!!!!!
   ! Misc Routines !
   !!!!!!!!!!!!!!!!!
   ! {{{
   ! behaves like Octave's tic()
   subroutine tic(t)
      integer, intent(out) :: t
      call system_clock(t)
   end subroutine tic
   
   ! returns time in seconds from now to time described by t
   real function toc_return(t)
      integer, intent(in) :: t
      integer :: now, clock_rate
   
      call system_clock(now, clock_rate)
      toc_return = real(now - t)/real(clock_rate)
   end function toc_return
   
   ! prints time in seconds from now to time described by t
   subroutine toc(t)
      integer, intent(in) :: t
      real (kind=8) :: time
      integer :: now, clock_rate
   
      call system_clock(now, clock_rate)
      time = real(now - t)/real(clock_rate)
      print *,"Elapsed time is ",time," seconds."
   end subroutine toc
   
   ! print a rank 2 allocatable array
   subroutine print_array(A)
      real (kind=8), allocatable, intent(in) :: A(:,:)
      character (len=30) :: rowfmt
   
      write(rowfmt, "(A,I4,A)") "(",size(A,2),"(1X,SS,10Es13.4))"
      row_print: do i=1,size(A,1)
         write(*, fmt=rowfmt) (A(i,j), j=1,size(A,2))
      end do row_print
   end subroutine print_array

   ! }}}

end module compla
