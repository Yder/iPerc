!==========================!
! TEST DISJOINT SET MODULE !
!==========================!
!
!> @author Yder MASSON
!> @date October 7, 2014
!> @brief This is to test the disjoint set module
!> @details 1) create a few label or classes
!> 2) check all newly created classes are canonical classes
!> 3) union some label or classes
!> 4) check that the union have been performed correctly
!
!========================
program test_disjoint_set
!========================
  !
  use module_disjoint_set
  !
  implicit none
  !
  integer, parameter :: nclass=10 !< number of class to create
  integer :: i !< looping index
  integer :: c1 !< class label
  integer :: c2 !< class label
  !
  ! allocate some memory to the data structure (this is optional)
  !
  call allocate_disjoint_set(2) ! we use 2<nclass to check that reallocation works fine
  ! 
  ! initialize largest_label (no class have yet been created)
  !
  largest_label = 0
  !
  ! create new classes
  !
  print*,
  !
  do i = 1,nclass
     print '(a,i3.0)','we just created a new class with label: ',create_set(largest_label)
  enddo
  !
  print*,
  !
  ! find the canonical class (label) of the newly created class
  ! 
  do i = 1,nclass
     c1 = i
     print '(a,i3.0,a,i3.0)','the canonical label of class ',c1,' is ',find(c1)
  enddo
  !
  print*,
  !
  ! proceed to some unions
  !
  c1 = 1 ; c2 = 3
  !
  print '(a,i3.0,a,i3.0,a,i3.0)',"let's union class ",c1,' and ',c2,&
       ', they are now part of one and the same class:',union(c1,c2)
  !
  c1 = 4 ; c2 = 7
  !
  print '(a,i3.0,a,i3.0,a,i3.0)',"let's union class ",c1,' and ',c2,&
       ', they are now part of one and the same class:',union(c1,c2)
  !
  c1 = 1 ; c2 = 4
  !
  print '(a,i3.0,a,i3.0,a,i3.0)',"let's union class ",c1,' and ',c2,&
       ', they are now part of one and the same class:',union(c1,c2)
  !
  c1 = 9 ; c2 = 10
  !
  print '(a,i3.0,a,i3.0,a,i3.0)',"let's union class ",c1,' and ',c2,&
       ', they are now part of one and the same class:',union(c1,c2)
  !
  print*,
  !
  ! check the new equivlent classes
  !
  do i = 1,nclass
     c1 = i
     print '(a,i3.0,a,i3.0)','the canonical label of class ',c1,' is ',find(c1)
  enddo
  print*,
  !
  ! clean data structure
  !
  call deallocate_disjoint_set()
  !
  stop
  !
!
!============================
end program test_disjoint_set
!============================
