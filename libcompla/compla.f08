! libcompla.f90
! Computational Linear Algebra library
! James Folberth - Spring 2014

! wrap lib in module so it interfaces properly
module compla

   real (kind=8), parameter :: ZERO_TOL = 10**(-13)

   contains

   !!!!!!!!!!!!!!!!!!!
   ! Cholesky Decomp !
   !!!!!!!!!!!!!!!!!!!
   ! {{{ 
   ! This uses the outer product formulation of the Cholesky decomposition
   ! Input matrix must by symmetric, positive definite.
   ! TODO optional input for transpose output or not

   ! Outer product formulation
   subroutine chol(A)
   real (kind=8) :: A(:,:)

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
   ! TODO

   !subroutine check_sym()

   ! }}}


   !!!!!!!!!!!!!
   ! LU Decomp !
   !!!!!!!!!!!!!
   ! {{{
   subroutine lu(A,p)
      real (kind=8) :: A(:,:), p(:)
      
      integer (kind=4) :: i,m,k,Nc,Nr
      real (kind=8) :: col_max

      Nc = size(A,1)
      Nr = size(A,2)

      if (Nc /= Nr ) then
         print *, "error: compla.f90: lu: input matrix is not square"
         stop
      end if

      row: do k=1,Nc-1
         m = maxloc(abs(A(k:Nc,k)))
         col_max = A(m,k)

      end do row
         



   end subroutine lu


   !subroutine apply_perm_vector(L,p,trans)
   !   real (kind=8) :: L(:,:)
   !   integer (kind=4) :: p(:), trans
   !  
   !   ! ``compute'' P*L
   !   if (trans == 0 ) then
   !
   !   ! `` compute P'*L
   !   else
   !
   !   end if
   !
   !end subroutine apply_perm_vector





   ! }}}


   !!!!!!!!!!!!!!!!!!
   ! For/Back Solve !
   !!!!!!!!!!!!!!!!!!
   ! {{{
   subroutine back_solve(U,b)
      real (kind=8) :: U(:,:), b(:,:)
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

   subroutine back_solve_blk(U,b)
      real (kind=8) :: U(:,:), b(:,:)

      integer (kind=4), parameter :: blk_size=4
      integer (kind=4) :: i,j,s,Nc
      integer (kind=4) :: il, ih, jl, jh

      ! call check_square(U)
      ! call check_upper_tri(U)

      Nc = size(U,1)
      s = Nc / blk_size

      ! Column oriented backward substitution
      blk: do j=1,s
         if (j>1) then
            col_blk: do i=j-1,s-1
               il = Nc-blk_size*(i+1)+1
               ih = Nc-blk_size*i
               jl = Nc-blk_size*(j-1)+1
               jh = Nc-blk_size*(j-2)
               !print *, "subtract block: ", il,ih,jl,jh
               
               b(il:ih,1) = b(il:ih,1) - matmul(U(il:ih,jl:jh),b(jl:jh,1))

            end do col_blk

            ! top block (not necc. blk_size by blk_size)
            il = 1
            ih = Nc-blk_size*s
            jl = Nc-blk_size*(j-1)+1
            jh = Nc-blk_size*(j-2)
            !print *, "subtract top block: ", il,ih,jl,jh

            b(il:ih,1) = b(il:ih,1) - matmul(U(il:ih,jl:jh),b(jl:jh,1))

         end if

         ! call back_solve on the diagonal blocks
         jl = Nc-blk_size*j+1
         jh = Nc-blk_size*(j-1)
         !print *, "back_solve: ",jl,jh

         call back_solve(U(jl:jh,jl:jh),b(jl:jh,:))

      end do blk

      ! subtract final top block
      if (s>0) then
         il = 1
         ih = Nc-blk_size*s
         jl = Nc-blk_size*s+1
         jh = Nc-blk_size*(s-1)
         !print *, "subtract final top block: ", il,ih,jl,jh
         
         b(il:ih,1) = b(il:ih,1) - matmul(U(il:ih,jl:jh),b(jl:jh,1))
      end if

      ! Finish with regular back solve
      row_fin: do j=Nc-blk_size*s,1,-1
         ! zero on diagonal
         if (abs(U(j,j)) <= ZERO_TOL) then
            print *, "error: compla.f90: back_solve: input matrix is singular to tolerance"
            stop
         end if

         b(j,1) = b(j,1)/U(j,j)
         col_fin: do i=j-1,1,-1
            b(i,1)=b(i,1) - U(i,j)*b(j,1)
         end do col_fin
      end do row_fin

   end subroutine back_solve_blk

   
   subroutine for_solve(L,b)
      real (kind=8) :: L(:,:), b(:,:)
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

   subroutine for_solve_blk(L,b)
      real (kind=8) :: L(:,:), b(:,:)

      integer (kind=4), parameter :: blk_size=8
      integer (kind=4) :: i,j,s,Nc
      integer (kind=4) :: il,ih,jl,jh

      ! call check_square(L)
      ! call check_lower_tri(L)

      Nc = size(L,1)
      s = Nc / blk_size

      ! Do forward subs in square blocks as far as possible
      ! column oriented
      blk: do j=1,s
         !print *, "j=",j
         if (j>1) then
            col_blk: do i=j,s
               il = blk_size*(i-1)+1
               ih = blk_size*i
               jl = blk_size*(j-2)+1
               jh = blk_size*(j-1)
               ! print *, il,ih,jl,jh
               
               b(il:ih,1) = b(il:ih,1) - matmul(L(il:ih,jl:jh),b(jl:jh,1))

            end do col_blk

            ! subtract bottom block (not necc. of size blk_size by blk_size)
            il = blk_size*s+1
            ih = Nc
            jl = blk_size*(j-2)+1
            jh = blk_size*(j-1)
            ! print *, il,ih,jl,ih

            b(il:ih,1) = b(il:ih,1) - matmul(L(il:ih,jl:jh),b(jl:jh,1))

         end if
 
         !print *, "calling for_solve"
         jl = blk_size*(j-1)+1
         jh = blk_size*j

         call for_solve(L(jl:jh,jl:jh), b(jl:jh,:))

      end do blk

      ! subtract final bottom block
      if (s>0) then
         il = blk_size*s+1
         ih = Nc
         jl = blk_size*(s-1)+1
         jh = blk_size*s

         b(il:ih,1) = b(il:ih,1) - matmul(L(il:ih,jl:jh),b(jl:jh,1))
      end if



      ! Finish up with regular forward subs
      !print *, blk_size*s
      row_fin: do j=blk_size*s+1,Nc
         ! zero on diagonal means singular matrix
         if (abs(L(j,j)) <= ZERO_TOL) then
            print *, "error: compla.f90: for_solve_blk: input matrix is singular to tolerance"
            stop
         end if

         b(j,1) = b(j,1)/L(j,j)
         col_fin: do i=j+1,Nc
            b(i,1)=b(i,1) - L(i,j)*b(j,1) 
         end do col_fin
      end do row_fin

   end subroutine for_solve_blk

   
   subroutine fb_solve_chol(A,b,x)
      real (kind=8) :: A(:,:), b(:,:), x(:,:)
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

   subroutine fb_solve_blk_chol(A,b,x)
      real (kind=8) :: A(:,:), b(:,:), x(:,:)
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
      call for_solve_blk(wrk, x) !A*x = R'*R*x=b -> R'*y=b
      wrk = transpose(wrk)
      call back_solve_blk(wrk,x) ! Rx=y
   
   end subroutine fb_solve_blk_chol

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
      real (kind=8), intent(in) :: A(:,:)
      character (len=30) :: rowfmt
   
      write(rowfmt, "(A,I4,A)") "(",size(A,2),"(1X,SS,10Es13.4))"
      row_print: do i=1,size(A,1)
         write(*, fmt=rowfmt) (A(i,j), j=1,size(A,2))
      end do row_print
   end subroutine print_array

   ! Frobenius matrix norm/vector 2-norm
   function norm_f(A)
      real (kind=8), intent(in) :: A(:,:)
      real (kind=8) :: norm_f

      integer (kind=4) :: i,j,Nr,Nc

      Nr=size(A,1)
      Nc=size(A,2)

      norm_f = 0d0

      row: do j=1,Nc
         col: do i=1,Nr
            norm_f = norm_f + A(i,j)*A(i,j)
         end do col
      end do row

      norm_f = sqrt(norm_f)
   
   end function norm_f


   function rand_spd_mat(N)
      integer (kind=4), intent(in) :: N
      real (kind=8), allocatable :: rand_spd_mat(:,:)
      
      real (kind=8) :: rand, m
      integer (kind=4) :: i,j

      allocate(rand_spd_mat(N,N))
      m=0d0

      ! Populate with random numbers [0,10]
      row: do j=1,N
         col: do i=1,N
            call random_number(rand)
            rand = 10d0*rand
            rand_spd_mat(i,j) = rand

            if (rand > m) m=rand

         end do col
      end do row

      ! Make symmetric
      rand_spd_mat = rand_spd_mat + transpose(rand_spd_mat)

      ! Make diagonally dominant (implies positive definite)
      diag: do j=1,N
         rand_spd_mat(j,j) = rand_spd_mat(j,j) + 10d0*N
      end do diag

   end function rand_spd_mat


   function rand_mat(Nr,Nc)
      integer (kind=4), intent(in) :: Nr, Nc
      real (kind=8), allocatable :: rand_mat(:,:)

      real (kind=8) :: rand
      integer (kind=4) :: i,j

      allocate(rand_mat(Nr,Nc))

      ! Populate with random numbers [0,1]
      row: do j=1,Nc
         col: do i=1,Nr
            call random_number(rand)
            rand_mat(i,j) = rand
         end do col
      end do row
   end function rand_mat


   subroutine init_random_seed()
     ! Stolen from http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
     ! {{{
     implicit none
     integer, allocatable :: seed(:)
     integer :: i, n, un, istat, dt(8), pid, t(2), s
     integer(8) :: count, tms
   
     call random_seed(size = n)
     allocate(seed(n))
     ! First try if the OS provides a random number generator
     open(newunit=un, file="/dev/urandom", access="stream", &
          form="unformatted", action="read", status="old", iostat=istat)
     if (istat == 0) then
        read(un) seed
        close(un)
     else
        ! Fallback to XOR:ing the current time and pid. The PID is
        ! useful in case one launches multiple instances of the same
        ! program in parallel.
        call system_clock(count)
        if (count /= 0) then
           t = transfer(count, t)
        else
           call date_and_time(values=dt)
           tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                + dt(5) * 60 * 60 * 1000 &
                + dt(6) * 60 * 1000 + dt(7) * 1000 &
                + dt(8)
           t = transfer(tms, t)
        end if
        s = ieor(t(1), t(2))
        pid = getpid() + 1099279 ! Add a prime
        s = ieor(s, pid)
        if (n >= 3) then
           seed(1) = t(1) + 36269
           seed(2) = t(2) + 72551
           seed(3) = pid
           if (n > 3) then
              seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
           end if
        else
           seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
        end if
     end if
     call random_seed(put=seed)
     ! }}}
   end subroutine init_random_seed


   ! }}}

end module compla
