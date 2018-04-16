program x_matrix

use mpi
use w90_io
use w90_parameters
use w90_constants, only : dp, eps6
use w90_kmesh
use w90_get_oper, only: get_SS_R, SS_R
use w90_comms, only : on_root,num_nodes, comms_setup, comms_end, comms_bcast
use w90_postw90_common, only : v_matrix,pw90common_wanint_data_dist,pw90common_wanint_param_dist, &
    pw90common_wanint_setup,pw90common_wanint_get_kpoint_file,nrpts


implicit none

integer :: nwann
integer :: i,r,n,m
integer :: u
integer :: nkp,len_seedname
logical :: have_gamma
real(dp) :: t1,t2

call comms_setup

if (on_root) then
    t1 = mpi_wtime()
end if

if (on_root) then
     call io_get_seedname
     len_seedname = len(seedname)
  end if
call comms_bcast(len_seedname,1)
call comms_bcast(seedname,len_seedname)

if (on_root) then 
   call param_read
   !call param_postw90_write
   if(.not.effective_model) then
          ! Check if the q-mesh includes the gamma point
          !
          have_gamma=.false.
          do nkp=1,num_kpts
             if (all(abs(kpt_latt(:,nkp))<eps6)) have_gamma=.true.       
          end do
          if(.not. have_gamma) write(stdout,'(1x,a)')&
               'Ab-initio does not include Gamma. Interpolation may be incorrect!!!'
          !
          ! Need nntot, wb, and bk to evaluate WF matrix elements of
          ! the position operator in reciprocal space. Also need
          ! nnlist to compute the additional matrix elements entering
          ! the orbital magnetization
          !
          call kmesh_get
       endif
endif

call pw90common_wanint_param_dist
if(on_root) call param_read_chkpt()
call  pw90common_wanint_data_dist

! Read list of k-points in irreducible BZ and their weights
!
! Should this be done on root node only?
!
if(wanint_kpoint_file) call pw90common_wanint_get_kpoint_file

! Setup a number of common variables for all interpolation tasks
!
call  pw90common_wanint_setup

call get_SS_R

nwann = size(SS_R,1)

if (on_root) then
    open(newunit=u, file='smat_wann.dat')

    do i=1,3
        do r=1,nrpts
            do n=1,nwann
                do m=1,nwann
                    write(u,*) i,r,n,m,SS_R(n,m,r,i)
                end do
            end do
        end do
    end do
end if

if (on_root) then
    t2 = mpi_wtime()
    write(*,*) 'Time to get the S matrix:', t2-t1,'s'
end if

call comms_end

end program
