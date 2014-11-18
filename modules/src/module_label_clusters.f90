!===============================!
!                               !
!                               !
!     MODULE LABEL CLUSTERS     !
!                               !
!                               !
!===============================!
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
!> @date October 23, 2014
!> @brief This module contains functions for cluster labeling, 
!> these are used to identify trapped regions.
!
!===========================
module module_label_clusters
!===========================
  !
  use module_disjoint_set 
  !
contains
  !
  !> @author Yder MASSON
  !> @date October 23, 2014
  !> @brief Label clusters on a 3D cubic lattice
  !-----------------------------------------------------------------------
  subroutine label_clusters_cubic(mat,nx,ny,nz,period_x,period_y,period_z)
  !-----------------------------------------------------------------------
    !
    implicit none
    !
    integer :: nx !< Grid dimension in the x direction (\b input)
    integer :: ny !< Grid dimension in the y direction (\b input)
    integer :: nz !< Grid dimension in the z direction (\b input)
    !> Matrix of clusters labels (\b input /\b output)
    !> \n At input time, <b>mat(i,j,k)=-1</b> at sites belonging to cluster to be labeled 
    !> (sites filled with the defending fluid)
    !> and <b>mat(i,j,k)= 0</b> otherwise (sites filled with the invading fluid)
    !> \n At output time, <b>mat(i,j,k)=L</b> where \b L is the label 
    !> of the cluster to which the site with index (i,j,k) belongs to
    integer :: mat(nx,ny,nz) 
    logical :: period_x !< Flag for periodic boundaries in the x direction (\b input)
    logical :: period_y !< Flag for periodic boundaries in the y direction (\b input)
    logical :: period_z !< Flag for periodic boundaries in the z direction (\b input)
    !
    ! internal variables
    ! 
    integer :: i, j, k, n, lc(3)
    !
    ! init largest label
    !
    largest_label = 0
    !
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             if(mat(i,j,k)<0)then
                n = 0
                if(i/=1.and.mat(i-1,j,k)/=0)then
                   n = n+1
                   lc(n) = mat(i-1,j,k)
                endif
                if(j/=1.and.mat(i,j-1,k)/=0)then
                   n = n+1
                   lc(n) = mat(i,j-1,k)
                endif
                if(k/=1.and.mat(i,j,k-1)/=0)then
                   n = n+1
                   lc(n) = mat(i,j,k-1)
                endif
                mat(i,j,k) = get_label_mat(n,lc)
             endif
          enddo
       enddo
    enddo
    !
    ! apply periodic boundaries
    !
    if(period_x)then
       do k = 1,nz
          do j = 1,ny
             if(mat(1,j,k)/=0.and.mat(nx,j,k)/=0)then
                mat(1,j,k)=union(mat(1,j,k),mat(nx,j,k))
             endif
          enddo
       enddo
    endif
    !
    if(period_y)then
       do k = 1,nz
          do i = 1,nx
             if(mat(i,1,k)/=0.and.mat(i,ny,k)/=0)then
                mat(i,1,k)=union(mat(i,1,k),mat(i,ny,k))
             endif
          enddo
       enddo
    endif
    !
    if(period_z)then
       do j = 1,ny
          do i = 1,nx
             if(mat(i,j,1)/=0.and.mat(i,j,nz)/=0)then
                mat(i,j,1)=union(mat(i,j,1),mat(i,j,nz))
             endif
          enddo
       enddo
    endif
    !
    ! find canonical labels
    !
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             if(mat(i,j,k)/=0)mat(i,j,k) = find(mat(i,j,k))
          enddo
       enddo
    enddo
    !
    return
    !
  !----------------------------------
  end subroutine label_clusters_cubic
  !----------------------------------
  !
  !> @author Yder MASSON
  !> @date October 23, 2014
  !> @brief  Label clusters on arbitrary lattices
  !--------------------------------------------------------------------
  subroutine label_clusters_arbitrary(mat,n_sites,offsets,connectivity)
  !--------------------------------------------------------------------
    !
    implicit none
    !
    integer :: n_sites !< Number of sites invaded (\b output)
    !> Matrix of clusters labels (\b input /\b output)
    !> \n At input time, <b>mat(i)=-1</b> at sites belonging to cluster to be labeled 
    !> (sites filled with the defending fluid)
    !> and <b>mat(i)= 0</b> otherwise (sites filled with the invading fluid)
    !> \n At output time, <b>mat(i)=L</b> where \b L is the label 
    !> of the cluster to which the site with index \b i belongs to
    integer :: mat(n_sites)
    !> Array containing the offsets of the data stored 
    !> in the connectivity array (\b Input)
    integer :: offsets(:)
    !> Array containing thelattice connectivity (\b input)
    !> \n For a given site \b i, with offset <b>j = offsets(i)</b> :
    !> \n <b>n=connectivity(j)</b> is the number of sites neighboring site \b i.
    !> \n <b>connectivity(j+1, j+2, ... ,j+n)</b> 
    !> contains the indices of the sites that are neighboring site \b i.
    integer :: connectivity(:)
    !
    ! intenal variables
    !
    integer :: n,i,i_site,i_start,i_end,n_neighbors
    integer, allocatable :: lc(:)
    !
    largest_label = 0
    !
    do i_site = 1,n_sites
       !
       if(mat(i_site)<0)then
          !
          i = offsets(i_site)
          n_neighbors = connectivity(i)
          i_start = i+1
          i_end   = i+n_neighbors
          !
          if(.not.allocated(lc))allocate(lc(n_neighbors))
          if(n_neighbors>size(lc))then
             deallocate(lc)
             allocate(lc(n_neighbors))
          endif
          !
          n = 0
          !
          do i = i_start,i_end
             if(mat(connectivity(i))>0)then
                n = n+1
                lc(n) = mat(connectivity(i))
             endif
          enddo
          !
          mat(i_site) = get_label_mat(n,lc)
          !
       endif
    enddo
    !
    if(allocated(lc))deallocate(lc)
    !
    ! find canonical label
    !
    do i = 1,n_sites
       if(mat(i)/=0)mat(i) = find(mat(i))
    enddo
    !
    return
    !
  !--------------------------------------
  end subroutine label_clusters_arbitrary
  !--------------------------------------
  !
  !> @author Yder MASSON
  !> @date October 23, 2014
  !> @brief Get label based on the neighbor's label
  !> \n If the site has neighbors, the site's label is the union of all neighbors's labels
  !> \n If the site has no neighbor, a new label is created
  !-----------------------------------
  integer function get_label_mat(n,lc)
  !-----------------------------------
    !
    implicit none
    !
    integer :: n !< Number of neighbors
    integer :: lc(:) !< Array containing the indices of the neighbors
    !
    ! internal variables
    !
    integer :: j
    !
    if(n==0)then
       get_label_mat = create_set(largest_label)
    elseif(n==1)then
       get_label_mat = find(lc(1))
    else
       do j = 2,n
          get_label_mat = union(lc(1),lc(j))
       enddo
    endif
    !
    return
    !
  !-------------------------
  end function get_label_mat
  !-------------------------
  !
!===============================
end module module_label_clusters
!===============================
