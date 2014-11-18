!============================!
!                            !
!                            !
!     MODULE BINARY TREE     !
!                            !
!                            !
!============================!
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
!> @author Yder MASSON
!! @date October 7, 2014
!! @brief Defines the binary tree data structure 
!! and the associated \b add_branch and \b update_root procedures
!! @see \b A \b fast \b algorithm \b for \b invasion \b percolation
!! Y Masson, SR Pride - Transport in porous media, 2014 - Springer
!! @see http://link.springer.com/article/10.1007/s11242-014-0277-8
!! @see https://sites.google.com/site/ydermasson/
!
!========================
module module_binary_tree
!========================
  !
  implicit none
  !
  !> Tree dimension (i.e. the number of sites or nodes in the tree, not the memory size)
  integer :: treedim
  !> Tree array, each site with index \b i 
  !! points toward his parent site with index \b tree(i)
  integer, dimension(:), allocatable :: tree
  !
contains
  !
  !> @author Yder Masson
  !> @date October 6, 2014
  !> @brief Add a new site to the tree
  !> @see Y Masson, SR Pride - Transport in porous media, 2014 - Springer
  !> @see http://link.springer.com/article/10.1007/s11242-014-0277-8
  !> @see https://sites.google.com/site/ydermasson/
  !
  !--------------------------------
  subroutine add_branch(values,ind)
  !--------------------------------
    !
    ! add a new site to the tree
    !
    implicit none
    !> Internal index to keep track of our position in the tree
    integer inode
    !> Index of the site to be added to the tree
    integer :: ind
    !> Array containing sites's values
    real :: values(:)
    !
    treedim = treedim+1
    !
    ! check if tree is allocated
    !
    call allocate_binary_tree(treedim)
    !
    ! add new site to tree
    !
    tree(treedim) = ind
    !
    ! put new site at the right place in the tree
    !
    inode = treedim
    !
    if(treedim==1)return
    !
    do 
       !
       select case(values(tree(inode))<=values(tree(inode/2)))
       case(.true.)
          call swap(tree(inode),tree(inode/2))
          inode=inode/2
          if(inode==1)exit
       case(.false.)
          exit
       end select
       !
    enddo
    !
    return
    !
  !------------------------
  end subroutine add_branch
  !------------------------
  !
  !> @author Yder Masson
  !> @date October 6, 2014
  !> @brief Update the root of the binary tree
  !> @see Y Masson, SR Pride - Transport in porous media, 2014 - Springer
  !> @see http://link.springer.com/article/10.1007/s11242-014-0277-8
  !> @see https://sites.google.com/site/ydermasson/
  !
  !----------------------------------
  subroutine update_tree_root(values)
  !----------------------------------
    !
    implicit none
    !
    !> Index to keep track of position in the tree
    integer :: inode
    !> array containing the values atached to sites (tree nodes)
    real :: values(:)
    ! 
    inode = 1
    if(treedim.le.1)then
       treedim = 0
       return
    endif
    !
    do
       !
       if(2*inode+1<treedim)then
          !
          select case(values(tree(2*inode))<=values(tree(2*inode+1)))
          case(.true.)
             tree(inode)=tree(inode*2)
             inode = 2*inode
          case(.false.)
             tree(inode)=tree(inode*2+1)
             inode = 2*inode+1
          end select
          !
       elseif(2*inode+1.eq.treedim)then
          !
          if(values(tree(2*inode))<=values(tree(2*inode+1)))then
             tree(inode)=tree(2*inode)
             tree(2*inode)=tree(2*inode+1)
          else
             tree(inode)=tree(2*inode+1)
          endif
          treedim = treedim-1
          exit
          !
       elseif(2*inode.eq.treedim)then
          !
          tree(inode)=tree(inode*2)
          treedim = treedim-1
          exit
          !
       else
          if(inode.eq.1)exit
          tree(inode)=tree(treedim)
          treedim = treedim-1
          do
             if(values(tree(inode))<values(tree(inode/2)))then
                call swap(tree(inode),tree(inode/2))
                inode=inode/2
                if(inode==1)exit
             else
                exit
             endif
          end do
          exit
          !
       endif
       !
    enddo
    !
    return
    !
  !------------------------------
  end subroutine update_tree_root
  !------------------------------
  !
  !> @author Yder Masson
  !> @date October 6, 2014
  !> @brief exchange \b a and \b b values.
  !-------------------
  subroutine swap(a,b)
  !-------------------
    !
    ! swap two values a becomes b and b becomes a
    !
    implicit none
    !
    integer :: a !< input value 1
    integer :: b !< input value 2
    integer :: save_b !< tmp value
    !
    save_b=b
    b=a
    a=save_b
    !
    return
    !
  !------------------
  end subroutine swap
  !------------------
  !
  !> @author Yder Masson
  !> @date October 6, 2014
  !> @brief Allocate memory for the \b tree array.
  !> @details If the array is too small, use the intrinsic function \b move_alloc
  !> to increase the size of the \b tree array.
  !> To prevent reallocation, allocate enough memory in the first call.
  !
  !----------------------------------------
  subroutine allocate_binary_tree(new_size)
  !----------------------------------------
    !
    implicit none  
    !
    !> Desired dimension for the \b tree array.
    integer :: new_size
    !> Temporary array for use with fortran intrinsic function \b move_alloc.
    integer, dimension(:), allocatable :: buffer
    !
    ! make sure array is not yet allocated before allocating
    !
    if(.not.allocated(tree))allocate(tree(new_size))
    !
    ! allocate more memory if needed
    ! 
    if(new_size>size(tree))then
       allocate(buffer(2*size(tree)))
       buffer(1:size(tree)) = tree(:)
       call move_alloc(buffer,tree)
    endif
    !
    return
    !
  !----------------------------------
  end subroutine allocate_binary_tree
  !----------------------------------
  !
  !> @author Yder Masson
  !> @date October 6, 2014
  !> @brief Deallocate array \b tree
  !
  !----------------------------------
  subroutine deallocate_binary_tree()
  !----------------------------------
    !
    implicit none 
    !
    if(allocated(tree))deallocate(tree)
    !
    return
    !
  !------------------------------------
  end subroutine deallocate_binary_tree
  !------------------------------------
  !
!============================
end module module_binary_tree
!============================
