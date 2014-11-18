!==========================
! TEST BINARY TREE MODULE !
!========================== 
!
!> @author Yder MASSON
!> @date October 7, 2014
!> @brief This is to test the binary_tree module
!> @details This is to check the binary_tree module is working as expected: 
!> 1) generate random values
!> 2) construct the binary tree
!> 3) recursively pick and update the root node 
!> and check that the values are well sorted
!
!=======================
program test_binary_tree
!=======================
  !
  use module_binary_tree
  !
  implicit none 
  !
  integer, parameter :: nval = 10 !< number of values to be stored
  integer :: i !> looping index
  integer :: initmem !> initial dimension to allocate tree
  real :: values(nval) !< array of random values
  real :: sorted_values(nval) !< array of sorted values
  logical :: check_successful !< .true. if pass test
  !
  ! generate random values
  !
  do i = 1,nval
     values(i) = rand(0)
  enddo
  !
  ! allocate memory for tree (this is optional)
  ! we use initmem<nval to check that
  ! dynamic reallocation is working fine
  ! 
  initmem = 2
  !
  call allocate_binary_tree(initmem)
  !
  ! construct tree: add all values one after another
  !
  do i = 1,nval
     call add_branch(values,i)
  enddo 
  !
  ! Undo tree, i.e. recursively pick the root node
  !
  do i = 1,treedim  ! we should have treedim=nval now
     sorted_values(i) = values(tree(1)) !> pick the value of the root node
     call update_tree_root(values)
  enddo
  !
  check_successful = .true.
  !
  do i = 1,nval
     print*,i,values(i),sorted_values(i)
     if(i<nval.and.sorted_values(i)>sorted_values(i+1))check_successful=.false.
  enddo
  !
  if(check_successful)then
     print*,'binary tree module passed the test successfully !!!'
  else
     print*,'binary tree module failed the test :('
  endif
  !
  call deallocate_binary_tree()
  !
  stop
  !
!===========================
end program test_binary_tree
!===========================
