program p03
   use compla
   implicit none

   integer (kind=4) :: n,i
   real (kind=8), allocatable :: A(:,:), Lambda(:)

   n = 6;
   allocate(A(n,n))
   allocate(Lambda(n))

   ! build A
   A = 0d0
   do i=1,n-1
      A(i+1,i) = 1.0d0
      A(i,i) = 2.0d0
      A(i,i+1) = 1.0d0
   end do
   A(n,n) = 2.0d0

   A(4,3) = 0; A(3,4) = 0
   !call print_array(A)
   
   ! compute evals of A
   !Lambda = 0d0
   !call francis1(A,Lambda)
   !call print_vector(Lambda)

   ! compute evals of A
   Lambda = 0d0
   call francis2_driver(A,Lambda)
   call print_vector(Lambda)
   
   deallocate(A,Lambda)

   contains

      subroutine francis1(A,Lambda)
         ! {{{
         real (kind = 8), intent(inout) :: A(:,:), Lambda(:)
         logical :: PRINT_MSGS = .true.

         integer (kind=4) :: i,iter,n,ldA
         real (kind=8) :: shift,cs,sn,rotg_r,rotg_z

         ldA = size(A,1)
         n = size(A,1) ! note that n will decrease as we deflate from the
                       ! from the bottom right on up, but ldA (leading dim A)
                       ! will be fixed (needed for dgem*)

         iterate: do iter=1,max(2000,5*n)
            
            ! compute shift
            shift = shift_wilk(A,n)
            
            ! create bulge
            ! http://www.netlib.org/lapack/explore-html/de/d13/drotg_8f_source.html
            rotg_r = A(1,1) - shift
            rotg_z = A(2,1)
            call drotg(rotg_r,rotg_z,cs,sn) ! drotg changes the
            ! input da,db -> rotg_r,rotg_z.  these are the vec norm and
            ! something else. see the source at netlib.org (note I'm using
            ! ATLAS, not netlib LAPACK)
            call drot(ldA,A(1,1),ldA,A(2,1),ldA,cs,sn)
            call drot(ldA,A(1,1),1,A(1,2),1,cs,sn)

            ! chase bulge
            chase: do i=1,n-2
               rotg_r = A(i+1,i)
               rotg_z = A(i+2,i)
               call drotg(rotg_r,rotg_z,cs,sn)
               ! note that drotg overwrites the first two entries
               A(i+1,i) = rotg_r
               A(i+2,i) = 0.0d0
               call drot(n-i,A(i+1,i+1),ldA,A(i+2,i+1),ldA,cs,sn)
               call drot(ldA,A(1,i+1),1,A(1,i+2),1,cs,sn)
            end do chase

            ! attempt to deflate lower right entry
            if ( abs(A(n,n-1)) < epsilon(1.0d0)* &
                  (abs(A(n-1,n-1))+abs(A(n,n))) ) then
               ! we can deflate, so trim off bottom row and right col
               n = n-1
               if ( n == 1 ) then
                  if ( PRINT_MSGS ) then
                     print *, "number of iterations = ",iter
                  end if

                  !do i=1,ldA
                  !   Lambda(i) = A(i,i)
                  !end do
                  Lambda = diag(A)
                  return
               end if
            end if

         end do iterate

         stop('p02.f08: francis1: we have done too many iterations')

         ! }}}
      end subroutine francis1

      subroutine francis2_driver(A,Lambda)
         ! {{{
         real (kind=8), intent(inout) :: A(:,:), Lambda(:)

         integer (kind=4), allocatable :: deflates(:)
         integer (kind=4) :: def_len,n1,n2

         n1 = 1
         n2 = size(A,1)
         def_len = 0
         allocate(deflates(size(A,1)))
         deflates = 0

         call francis2(A,deflates,def_len,n1,n2)

         deallocate(deflates)
         ! }}}
      end subroutine francis2_driver

      recursive subroutine francis2(A,deflates,def_len,n1,n2)
         ! {{{
         real (kind=8), intent(inout) :: A(:,:)
         integer (kind=4), intent(inout) :: deflates(:),def_len
         integer (kind=4), intent(in) :: n1,n2

         logical :: PRINT_MSGS = .true.

         real (kind=8), allocatable :: maind(:),subd(:)
         real (kind=8) :: htr,dscr,det,root1,root2

         integer (kind=4), allocatable :: definds(:)
         integer (kind=4) :: iter,i,definds_len

         !! down to a 2x2 case; compute evals explicitly
         if ( n2 - n1 == 1 ) then
            htr = 0.5d0*(A(n1,n1)+A(n2,n2)) ! half trace
            dscr = hypot(0.5d0*(A(n1,n1)-A(n2,n2)),A(n2,n1)) ! descriminant
   
            if ( htr < ZERO_TOL ) then
               dscr = -dscr ! avoids cancellation
            end if
   
            root1 = htr+dscr ! quadradic formula, in terms of tr and det
            if ( abs(root1) < ZERO_TOL ) then
               root2 = 0d0
            else
               det = A(n1,n1)*A(n2,n2)-A(n2,n1)**2
               root2 = det/root1 ! det = root1*root2
            end if

            A(n1,n1) = root1; A(n2,n2) = root2
            A(n1,n2) = 0d0; A(n2,n1) = 0d0

            deflates(def_len+1) = n1
            deflates(def_len+2) = n2

            return

         end if

         !! try to deflate (initially)
         allocate(maind(size(A,1)), subd(size(A,1)-1)) ! store whole (sub)diag
         allocate(definds(n2-n1))
         maind = diag(A)
         subd = diag(A,-1)

         ! check for deflation condition
         definds = 0
         definds_len = 0
         do i=n1,n2-1
            ! if okay to deflate and haven't deflated
            if ( subd(i) < epsilon(1d0)*(maind(i)+maind(i+1)) &
               .and. count(i==deflates)==0 ) then
               subd(i) = 1 ! found one, overwrite subdiag vec 
               definds_len = definds_len + 1
               definds(definds_len) = i
            else
               subd(i) = 0
            end if
         end do

         !call print_vector(dble(definds(1:definds_len)))
         !stop('stuff')

         ! single deflatioin
         if ( definds_len == 1 ) then
            if ( PRINT_MSGS ) then 
               print *, 'deflate single (pre-iter): ',definds(1)+1
            end if

            def_len = def_len + 1
            deflates(def_len) = definds(1)+1
            call francis2(A,deflates,def_len,n1,definds(1))
            call francis2(A,deflates,def_len,definds(1)+1,n2)
            return

         ! multiple deflations (take them one at a time)
         elseif ( definds_len > 1 ) then
            stop('not implemented')

         end if

         ! deflate the bottom (separate check)
         ! in case the francis2 calls above already deflate this,
         ! don't try to do it again
         if ( abs(A(n2,n2-1)) < epsilon(1d0)*&
            (abs(A(n2-1,n2-1))+abs(A(n2,n2)))&
            .and. subd(n2-1) == 1 ) then
           
            ! trim off the bottom row and rightmost col
            if ( PRINT_MSGS ) print *, "deflate end (pre-iter): ",n2

            def_len = def_len + 1
            deflates(def_len) = n2
            call francis2(A,deflates,def_len,n1,n2-1)
            return
         end if

         
         deallocate(maind,subd)
         deallocate(definds)

         ! }}}
      end subroutine francis2

      ! Wilkinson shift
      function shift_wilk(A,n)
         ! {{{
         real (kind=8), intent(in) :: A(:,:)
         integer (kind=4), intent(in) :: n
         real (kind=8) :: shift_wilk

         real (kind=8) :: htr,dscr,det,root1,root2
         
         htr = 0.5d0*(A(n-1,n-1)+A(n,n)) ! half trace
         dscr = hypot(0.5d0*(A(n-1,n-1)-A(n,n)),A(n,n-1)) ! descriminant

         if ( htr < ZERO_TOL ) then
            dscr = -dscr ! avoids cancellation
         end if

         root1 = htr+dscr ! quadradic formula, in terms of tr and det
         if ( abs(root1) < ZERO_TOL ) then
            root2 = 0d0
         else
            det = A(n-1,n-1)*A(n,n)-A(n,n-1)**2
            root2 = det/root1 ! det = root1*root2
         end if

         if ( abs(A(n,n)-root1) < abs(A(n,n)-root2) ) then
            shift_wilk = root1
         else
            shift_wilk = root2
         end if
         ! }}}
      end function shift_wilk

end program p03

! vim: set ts=3 sw=3 sts=3 et :
! vim: foldmarker={{{,}}}
