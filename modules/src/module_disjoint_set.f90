!=============================!
!                             !
!                             !
!     MODULE DISJOINT SET     !  
!                             !
!                             !
!=============================!
!
!===========================================================================
!
!                            iPerc Version 1.0    
!                            -----------------
!
!    Main author: Yder MASSON, yder.masson@cal.berkeley.edu                      
!
!    This file is part of iPerc.
!
!    iPerc is a fortran library for modeling invasion percolation
!    Copyright (C) 2014  Yder MASSON
!
!    iPerc is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!===========================================================================
!
!> @author Yder Masson
!
!> @date October 6, 2014
!
!> @brief Define disjoint set data structure 
!! and the associated \b Union, \b Find and \b Create_set procedures.
!
!> @details Path compression and union by rank are implemented for efficiency. 
!! See subroutines's descriptions for more details.
!! 
!
!> @see http://en.wikipedia.org/wiki/Disjoint-set_data_structure
!
!=========================
module module_disjoint_set
!=========================
  !
  implicit none
  !> Array containing the label trees, 
  !! each label \b i points toward its parent \b labels(i)
  !! the root or canonical labels satisfies \b labels(i)=i
  integer, dimension(:), allocatable :: labels
  !> Array containing the label trees's ranks. 
  !! It is used for union by rank (also called weighted union)
  integer, dimension(:), allocatable :: ranks
  !> Number of labels currently used
  integer :: largest_label
  !
contains
  !
  !> @author Yder Masson
  !> @date October 6, 2014
  !> @brief Create a new label or class with label \b largest_label
  !! @details This function allocates more memory 
  !! to arrays \b labels and \b ranks if needed
  !> @see http://en.wikipedia.org/wiki/Disjoint-set_data_structure
  !
  !-----------------------------------------
  integer function create_set(largest_label)
  !-----------------------------------------
    !
    implicit none
    !
    !> Number of labels currently used
    integer :: largest_label
    !
    largest_label = largest_label+1
    call allocate_disjoint_set(largest_label)
    labels(largest_label) = largest_label
    ranks(largest_label) = 1
    create_set = largest_label
    !
  !----------------------
  end function create_set
  !----------------------
  !
  !> @author Yder Masson
  !> @date October 6, 2014
  !> @brief This function returns the canonical label 
  !! (tree's root) or class associated with \b label
  !! @details Path compression is implemented for more efficiency.
  !> @see http://en.wikipedia.org/wiki/Disjoint-set_data_structure
  !
  !---------------------------
  integer function find(label)
  !---------------------------
    !
    implicit none
    !
    integer :: label !< label for wich we search root or canonical label or class
    integer :: root_label !< output root label, find=root_label
    !> Used for path comperssion, i.e. first find the root label 
    !! then make all encountered tree nodes point toward the root label
    integer :: save_label 
    !
    root_label = label
    !
    do while(labels(root_label)/=root_label)
       root_label = labels(root_label)
    enddo
    !
    do while (labels(label)/=label)
       save_label = labels(label)
       labels(label) = root_label
       label = save_label
    end do
    !
    find = root_label
    !
  !----------------
  end function find
  !----------------
  !
  !> @author Yder Masson
  !> @date October 6, 2014
  !> @brief Set \b label1 and \b label2 to equivalence classes 
  !!(i.e. after the call \b label1 and \b label2 will point
  !! to the same canonical root label). 
  !! @details This function returns canonical label or class of the union.
  !! The union by rank alogorithm (also called weighted union) is implemented.
  !> @see http://en.wikipedia.org/wiki/Disjoint-set_data_structure
  !
  !------------------------------------
  integer function union(label1,label2)
  !------------------------------------ 
    !
    implicit none
    !
    integer :: label1 !< input label 1
    integer :: label2 !< input label 2
    !
    label1 = find(label1)
    label2 = find(label2)
    !
    if(ranks(label1)==ranks(label2))then
       ranks(label1) = ranks(label1)+1
       labels(label2) = label1
       union = label1
    elseif(ranks(label1)>ranks(label2))then
       labels(label2) = label1
       union = label1
    else
       labels(label1) = label2
       union = label2
    endif
    !
  !------------------
  end function  union
  !------------------
  !
  !> @author Yder Masson
  !> @date October 6, 2014
  !> @brief Allocate memory for the arrays \b labels and \b ranks.
  !! @details When needed, use the intrinsic function 
  !! \b move_alloc to increase the size of the arrays on the fly. 
  !! To avoid reallocation, allocate enough memory at first call.
  !! @see https://gcc.gnu.org/onlinedocs/gfortran/MOVE_005fALLOC.html
  !
  !-----------------------------------------
  subroutine allocate_disjoint_set(new_size)
  !-----------------------------------------
    !
    implicit none
    !
    !> Desired dimension for : \b ranks and \b labels arrays.
    integer :: new_size
    !> Temporary array for use with fortran intrinsic function \b move_alloc.
    integer, dimension(:), allocatable :: buffer
    !
    ! make sure arrays are allocated
    !
    if(.not.allocated(labels))allocate(labels(new_size))
    if(.not.allocated(ranks))allocate(ranks(new_size))
    !
    ! allocate more memory if needed
    ! 
    if(new_size>size(labels))then
       allocate(buffer(2*size(labels)))
       buffer(1:size(labels)) = labels(:)
       call move_alloc(buffer,labels)
    endif
    ! 
    if(new_size>size(ranks))then
       allocate(buffer(2*size(ranks)))
       buffer(1:size(ranks)) = ranks(:)
       call move_alloc(buffer,ranks)
    endif
    !
    return
    !
  !-----------------------------------
  end subroutine allocate_disjoint_set
  !-----------------------------------
  !
  !> @author Yder Masson
  !> @date October 6, 2014
  !> Deallocate the \b labels and \b ranks arrays.
  !
  !-----------------------------------
  subroutine deallocate_disjoint_set()
  !-----------------------------------
    !
    implicit none
    !
    if(allocated(labels))deallocate(labels)
    if(allocated(ranks))deallocate(ranks)
    !
    return
    !
  !-------------------------------------
  end subroutine deallocate_disjoint_set
  !-------------------------------------
  !
!=============================
end module module_disjoint_set
!=============================
