! p02.f08

program p02

   use compla ! this lib has a bunch of stuff I wrote throughout the semester
      ! stuff from libcompla: intk,dblk,diag,write_dset,print_vector/array
   use HDF5 ! to save data and read into GNU Octave
      ! stuff like h5* is from HDF5

   ! BLAS calls
   ! drotg,drot
   
   implicit none

   integer (kind=intk) :: n,i
   real (kind=dblk), allocatable :: A(:,:), Lambda(:)

   n = 6
   !n = 10
   n = 50
   allocate(A(n,n))
   allocate(Lambda(n))

   ! build A
   A = 0_dblk
   do i=1,n-1
      A(i+1,i) = 1.0_dblk
      A(i,i) = 2.0_dblk
      A(i,i+1) = 1.0_dblk
   end do
   A(n,n) = 2.0_dblk

   A(4,3) = 0; A(3,4) = 0
   A(6,5) = 0; A(5,6) = 0
   A(7,6) = 0; A(6,7) = 0 
   A(34,33) = 0; A(33,34) = 0
   A(35,34) = 0; A(34,35) = 0
   A(36,35) = 0; A(35,36) = 0

   ! compute evals of A (Francis 1)
   !Lambda = 0d0
   !call francis1(A,Lambda)
   !call print_vector(Lambda)

   ! compute evals of A (Francis 2)
   Lambda = 0_dblk
   call francis2_driver(A,Lambda)
   call print_vector(Lambda)

   call save_stuff("data.h5",A,Lambda)
   print *, "data saved to ","data.h5"

   deallocate(A,Lambda)

   contains

      ! Francis algo (deg 1, sym-tridiag, various shifts, somewhat deflated)
      subroutine francis1(A,Lambda)
         ! {{{
         real (kind=dblk), intent(inout) :: A(:,:), Lambda(:)
         logical :: PRINT_MSGS = .true.

         integer (kind=intk) :: i,iter,n,ldA
         real (kind=dblk) :: shift,cs,sn,rotg_r,rotg_z

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
               A(i+2,i) = 0.0_dblk
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

      ! Francis algo driver (def 1, sym-tridiag, various shifts, deflated)
      subroutine francis2_driver(A,Lambda)
         ! {{{
         real (kind=dblk), intent(inout) :: A(:,:), Lambda(:)

         integer (kind=intk), allocatable :: deflates(:)
         integer (kind=intk) :: def_len,n1,n2

         n1 = 1
         n2 = size(A,1)
         def_len = 0
         allocate(deflates(size(A,1)))
         deflates = 0

         call francis2(A,deflates,def_len,n1,n2)

         Lambda = diag(A)

         deallocate(deflates)
         ! }}}
      end subroutine francis2_driver

      ! Francis algo (def 1, sym-tridiag, various shifts, deflated) 
      recursive subroutine francis2(A,deflates,def_len,n1in,n2in)
         ! {{{
         real (kind=dblk), intent(inout) :: A(:,:)
         integer (kind=intk), intent(inout) :: deflates(:),def_len
         integer (kind=intk), intent(in) :: n1in,n2in

         logical :: PRINT_MSGS = .true.

         ! Deflate vars
         real (kind=dblk), allocatable :: maind(:),subd(:)
         real (kind=dblk) :: htr,dscr,det,root1,root2

         integer (kind=intk), allocatable :: definds(:)
         integer (kind=intk) :: iter,i,definds_len,n1,n2

         ! Francis iter vars
         integer (kind=intk) :: ldA
         real (kind=dblk) :: shift,cs,sn,rotg_r,rotg_z

         ldA = size(A,1)
         n1 = n1in
         n2 = n2in

         !! down to a 2x2 case; compute evals explicitly
         if ( n2 - n1 == 1 ) then
            htr = 0.5_dblk*(A(n1,n1)+A(n2,n2)) ! half trace
            dscr = hypot(0.5_dblk*(A(n1,n1)-A(n2,n2)),A(n2,n1)) ! descriminant
   
            if ( htr < ZERO_TOL ) then
               dscr = -dscr ! avoids cancellation
            end if
   
            root1 = htr+dscr ! quadradic formula, in terms of tr and det
            if ( abs(root1) < ZERO_TOL ) then
               root2 = 0_dblk
            else
               det = A(n1,n1)*A(n2,n2)-A(n2,n1)**2
               root2 = det/root1 ! det = root1*root2
            end if

            A(n1,n1) = root1; A(n2,n2) = root2
            A(n1,n2) = 0_dblk; A(n2,n1) = 0_dblk

            deflates(def_len+1) = n1
            deflates(def_len+2) = n2

            return

         end if

         !! try to deflate (initially)
         allocate(maind(size(A,1)), subd(size(A,1)-1)) ! store whole (sub)diag
         allocate(definds(n2-n1))
         maind = abs(diag(A))
         subd = abs(diag(A,-1))

         ! check for deflation condition
         definds = 0
         definds_len = 0
         do i=n1,n2-2
            ! if okay to deflate and haven't deflated
            if ( subd(i) < epsilon(1d0)*(maind(i)+maind(i+1)) &
               .and. count(i==deflates)==0 ) then
              
               ! if adjacent deflations, don't deflate top one yet (we'll 
               !get it (when we deflate the end)
               if ( definds_len > 0 ) then
                  if ( definds(definds_len) == i-1 ) then
                     definds(definds_len) = i
               
                  ! non-adjacent deflations
                  else
                     definds_len = definds_len + 1
                     definds(definds_len) = i
                  end if
               else
                  definds_len = definds_len + 1
                  definds(definds_len) = i
               end if
            end if
         end do

         !call print_vector(dble(definds(1:definds_len)))
         !stop('stuff')

         ! single deflation
         if ( definds_len == 1 ) then
            if ( PRINT_MSGS ) then 
               print *, 'deflate single (pre-iter): ',definds(1)+1
            end if

            def_len = def_len + 1
            deflates(def_len) = definds(1)+1
            call francis2(A,deflates,def_len,n1,definds(1))
            call francis2(A,deflates,def_len,definds(1)+1,n2)
            return ! done with n1:n2

         ! multiple deflations (take them one at a time)
         elseif ( definds_len > 1 ) then
            if ( PRINT_MSGS ) then 
               print *, 'deflate multiple (pre-iter, start): ',definds(1)+1
            end if

            def_len = def_len + 1
            deflates(def_len) = definds(1)+1
            call francis2(A,deflates,def_len,n1,definds(1))

            do i=2,definds_len
               if ( Print_MSGS ) then
                  print *, 'deflate multiple (pre-iter, mid): ',definds(i)+1
               end if

               def_len = def_len + 1
               deflates(def_len) = definds(i)+1
               call francis2(A,deflates,def_len,definds(i-1)+1,definds(i))
            end do
            
            if ( PRINT_MSGS ) then
               print *, 'deflate multiple (pre-iter, end): ',&
                  definds(definds_len)+1
            end if
            call francis2(A,deflates,def_len,definds(definds_len)+1,n2)
               
            return ! done with n1:n2
         end if

         ! deflate the bottom (separate check)
         ! in case the francis2 calls above already deflate this,
         ! don't try to do it again
         if ( abs(A(n2,n2-1)) < epsilon(1d0)*&
            (abs(A(n2-1,n2-1))+abs(A(n2,n2)))&
            .and. count(n2-1==definds)==1 ) then
           
            ! trim off the bottom row and rightmost col
            if ( PRINT_MSGS ) print *, "deflate end (pre-iter): ",n2

            def_len = def_len + 1
            deflates(def_len) = n2
            call francis2(A,deflates,def_len,n1,n2-1)
            return ! done with n1:n2
         end if


         !! Francis iterations
         iterate: do iter=1,max(2000,5*(n2-n1))

            ! compute shift
            !print *, "pre-shift",n1,n2
            shift = shift_wilk(A,n2)
            
            ! create bulge
            ! http://www.netlib.org/lapack/explore-html/de/d13/drotg_8f_source.html
            rotg_r = A(n1,n1) - shift
            rotg_z = A(n1+1,n1)
            call drotg(rotg_r,rotg_z,cs,sn) ! drotg changes the
            ! input da,db -> rotg_r,rotg_z.  these are the vec norm and
            ! something else. see the source at netlib.org (note I'm using
            ! ATLAS, not netlib LAPACK)
            call drot(ldA,A(n1,1),ldA,A(n1+1,1),ldA,cs,sn)
            call drot(ldA,A(1,n1),1,A(1,n1+1),1,cs,sn)

            ! chase bulge
            chase: do i=n1,n2-2
               rotg_r = A(i+1,i)
               rotg_z = A(i+2,i)
               call drotg(rotg_r,rotg_z,cs,sn)
               ! note that drotg overwrites the first two entries
               A(i+1,i) = rotg_r
               A(i+2,i) = 0.0_dblk
               call drot(n2-i,A(i+1,i+1),ldA,A(i+2,i+1),ldA,cs,sn)
               call drot(ldA,A(1,i+1),1,A(1,i+2),1,cs,sn)
            end do chase


            !! try to deflate (mid-iteration)
            maind = abs(diag(A))
            subd = abs(diag(A,-1))
   
            ! check for deflation condition
            definds = 0
            definds_len = 0
            do i=n1,n2-2
               ! if okay to deflate and haven't deflated
               if ( subd(i) < epsilon(1d0)*(maind(i)+maind(i+1)) &
                  .and. count(i==deflates)==0 ) then
                 
                  ! if adjacent deflations, don't deflate top one yet (we'll 
                  !get it (when we deflate the end)
                  if ( definds_len > 0 ) then
                     if ( definds(definds_len) == i-1 ) then
                        definds(definds_len) = i
                  
                     ! non-adjacent deflations
                     else
                        definds_len = definds_len + 1
                        definds(definds_len) = i
                     end if
                  else
                     definds_len = definds_len + 1
                     definds(definds_len) = i
                  end if
               end if
            end do
  
            ! single deflation
            if ( definds_len == 1 ) then
               if ( PRINT_MSGS ) then 
                  print *, 'deflate single (mid-iter): ',definds(1)+1
               end if
   
               def_len = def_len + 1
               deflates(def_len) = definds(1)+1
               call francis2(A,deflates,def_len,n1,definds(1))
               call francis2(A,deflates,def_len,definds(1)+1,n2)
               return ! done with n1:n2
   
            ! multiple deflations (take them one at a time)
            elseif ( definds_len > 1 ) then
               if ( PRINT_MSGS ) then 
                  print *, 'deflate multiple (pre-iter, start): ',definds(1)+1
               end if
   
               def_len = def_len + 1
               deflates(def_len) = definds(1)+1
               call francis2(A,deflates,def_len,n1,definds(1))
   
               do i=2,definds_len
                  if ( Print_MSGS ) then
                     print *, 'deflate multiple (pre-iter, mid): ',definds(i)+1
                  end if
   
                  def_len = def_len + 1
                  deflates(def_len) = definds(i)+1
                  call francis2(A,deflates,def_len,definds(i-1)+1,definds(i))
               end do
               
               if ( PRINT_MSGS ) then
                  print *, 'deflate multiple (pre-iter, end): ',&
                     definds(definds_len)+1
               end if
               call francis2(A,deflates,def_len,definds(definds_len)+1,n2)
                  
               return ! done with n1:n2
            end if

            ! deflate the bottom
            if ( abs(A(n2,n2-1)) < epsilon(1d0)*&
               (abs(A(n2-1,n2-1))+abs(A(n2,n2))) ) then
              
               ! trim off the bottom row and rightmost col
               if ( PRINT_MSGS ) print *, "deflate end (mid-iter): ",n2

               def_len = def_len + 1
               deflates(def_len) = n2
               n2  = n2 - 1
              
               ! down to 2x2 case
               if ( n2 - n1 == 1 ) then
                  !print *, n1,n2
                  call francis2(A,deflates,def_len,n1,n2)
                  return
               end if
               !call francis2(A,deflates,def_len,n1,n2-1)
               !return ! done with n1:n2
            end if

         end do iterate

         
         deallocate(maind,subd)
         deallocate(definds)

         ! }}}
      end subroutine francis2

      ! Wilkinson shift
      function shift_wilk(A,n)
         ! {{{
         real (kind=dblk), intent(in) :: A(:,:)
         integer (kind=intk), intent(in) :: n
         real (kind=dblk) :: shift_wilk

         real (kind=dblk) :: htr,dscr,det,root1,root2
         
         htr = 0.5_dblk*(A(n-1,n-1)+A(n,n)) ! half trace
         dscr = hypot(0.5_dblk*(A(n-1,n-1)-A(n,n)),A(n,n-1)) ! descriminant

         if ( htr < ZERO_TOL ) then
            dscr = -dscr ! avoids cancellation
         end if

         root1 = htr+dscr ! quadradic formula, in terms of tr and det
         if ( abs(root1) < ZERO_TOL ) then
            root2 = 0_dblk
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

      ! Rayleigh shift
      function shift_ray(A,n)
         ! {{{
         real (kind=dblk), intent(in) :: A(:,:)
         real (kind=dblk) :: shift_ray
         integer (kind=intk), intent(in) :: n

         shift_ray = A(n,n)
         ! }}}
      end function shift_ray

      ! save a A,Lambda in the top dir of HDF5 file
      subroutine save_stuff(savefile,A,Lambda)
         ! {{{
         character (len=*), intent(in) :: savefile
         real (kind=dblk), intent(in) :: A(:,:),Lambda(:)

         integer (kind=intk) :: h5error, file_id

         call h5open_f(h5error)
         call h5fcreate_f(savefile, H5F_ACC_TRUNC_F, file_id, h5error)
         
         call write_dset(file_id, A, "A")
         call write_dset(file_id, Lambda, "Lambda")

         call h5fclose_f(file_id,h5error)
         call h5close_f(h5error)
         ! }}}
      end subroutine save_stuff

end program p02

! vim: set ts=3 sw=3 sts=3 et :
! vim: foldmarker={{{,}}}
