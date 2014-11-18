!=========================!
!                         !
!                         !
!     MODULE TRAPPING     !
!                         !
!                         !
!=========================!
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
!> @brief This module contains functions to identify trapped sites 
!
!=====================
module module_trapping
!=====================
  !
  use module_cubic_indices
  use module_disjoint_set
  use module_label_clusters
  use module_invasion_percolation_constants
  !
contains
  !
  !> @author Yder MASSON
  !> @date October 23, 2014
  !> @brief Ths function does the union  of the labels 
  !> of all sites connected to an exit site and return 
  !> the corresponding label that is the free cluster's label
  !-----------------------------------------------------------
  integer function get_label_free_clusters(n_sites,states,mat)
  !-----------------------------------------------------------
    !
    integer :: states(n_sites) !< Array containing the sites's states (\b input)
    integer :: mat(n_sites) !< Array containing the clusters's labels (\b input)
    integer :: label_free_clusters !< free cluster's label
    !
    ! create a new set
    !
    label_free_clusters = create_set(largest_label)
    !
    ! union all exit sites's clusters
    ! (fluid is free to escape at all sites
    ! connected to exit sites)
    !
    do i = 1,n_sites
       if(states(i)==exit_site)then
          label_free_clusters = union(mat(i),label_free_clusters)
       endif
    enddo
    !
    get_label_free_clusters = label_free_clusters
    !
    return
    !
  !-----------------------------------
  end function get_label_free_clusters
  !-----------------------------------
  !
  !> @author Yder MASSON
  !> @date October 23, 2014
  !> @brief Identify trapped sites on cubic lattices
  !> @details It uses the fast a posteriori identification of trapped node 
  !> if the <b>undo_invasion=.true.</b>. 
  !> If <b>undo_invasion=.false.</b> cluster labeling is performed at each invasion steps.
  !-------------------------------------------------------
  subroutine find_trapped_sites_cubic(nx,                &
                                      ny,                &
                                      nz,                &
                                      states,            &
                                      period_x,          &
                                      period_y,          &
                                      period_z,          &
                                      n_sites_invaded,   &
                                      invasion_list,     &
                                      undo_invasion      )
  !-------------------------------------------------------
    !
    implicit none
    !
    integer :: nx !< Grid dimension in the x direction (\b input)
    integer :: ny !< Grid dimension in the y direction (\b input)
    integer :: nz !< Grid dimension in the z direction (\b input)
    !
    integer :: states(nx,ny,nz) !< Array containing the sites's states (\b input/\b output)
    !! \n The state of trapped sites is updated to <b>state(i,j,k)=trapped</b>
    !
    logical :: period_x !< Flag for periodic boundaries in the x direction (\b input)
    logical :: period_y !< Flag for periodic boundaries in the y direction (\b input)
    logical :: period_z !< Flag for periodic boundaries in the z direction (\b input)
    !
    logical :: undo_invasion !< Flag: If <b>undo_invasion==.true.</b>
    !! Then use fast a posteriori method for trapping (\b input)
    integer :: n_sites_invaded !< Number of invaded sites (\b input)
    integer :: invasion_list(:) !< List of invaded sites 
    !! sorted in chronological order (\b input)
    !
    ! internal variables
    ! 
    integer :: label_free_clusters !< Free cluster's label
    integer, allocatable, dimension(:,:,:) :: mat !< bufer for cluster labeling
    integer :: i,j,k,n,in1,in2,jn1,jn2,kn1,kn2,is,js,ks,i_site,lc(6)
    integer :: n_sites_invaded_with_trapping
    !
    ! init mat
    !
    allocate(mat(nx,ny,nz))
    !
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             if(states(i,j,k)==invaded)then
                mat(i,j,k) = 0
             else
                mat(i,j,k) =-1
             endif
          enddo
       enddo
    enddo
    !
    ! label cluster
    !
    call label_clusters_cubic(mat,nx,ny,nz,period_x,period_y,period_z)
    !
    ! find free cluster's label
    !
    label_free_clusters = get_label_free_clusters(nx*ny*nz,states,mat)
    !
    ! update trapped sites's states
    !
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             if(mat(i,j,k)/=0.and.mat(i,j,k)/=label_free_clusters)then
                states(i,j,k) = trapped
             endif
          enddo
       enddo
    enddo
    !
    if(undo_invasion)then
       !
       do i_site = n_sites_invaded,1,-1
          !
          ! get site's 3d indices
          !
          call ind2ijk(nx,ny,invasion_list(i_site),is,js,ks)
          !   
          in1 = is-1
          in2 = is+1
          jn1 = js-1
          jn2 = js+1
          kn1 = ks-1
          kn2 = ks+1
          !
          ! account for periodic boundaries
          !
          if(in1<1 .and.period_x) in1 = in1+nx
          if(in2>nx.and.period_x) in2 = in2-nx
          if(jn1<1 .and.period_y) jn1 = jn1+ny
          if(jn2>ny.and.period_y) jn2 = jn2-ny
          if(kn1<1 .and.period_z) kn1 = kn1+nz
          if(kn2>nz.and.period_z) kn2 = kn2-nz
          !
          !  scan neighboring sites
          !
          n = 0
          !
          if(in1>=1 .and.mat(in1,js,ks)/=0)then 
             n = n + 1
             lc(n) = mat(in1,js,ks) 
          endif
          if(in2<=nx.and.mat(in2,js,ks)/=0)then
             n = n + 1
             lc(n) = mat(in2,js,ks)
          endif
          if(jn1>=1 .and.mat(is,jn1,ks)/=0)then
             n = n + 1
             lc(n) = mat(is,jn1,ks) 
          endif
          if(jn2<=ny.and.mat(is,jn2,ks)/=0)then 
             n = n + 1
             lc(n) = mat(is,jn2,ks) 
          endif
          if(kn1>=1 .and.mat(is,js,kn1)/=0)then 
             n = n + 1
             lc(n) = mat(is,js,kn1) 
          endif
          if(kn2<=nz.and.mat(is,js,kn2)/=0)then 
             n = n + 1
             lc(n) = mat(is,js,kn2) 
          endif
          !
          mat(is,js,ks) = get_label_mat(n,lc)
          !
          if(find(mat(is,js,ks))/=find(label_free_clusters))then
             states(is,js,ks) = trapped
          endif
          !
       enddo
       !
       ! remove trapped sites from invasion_list
       !
       n_sites_invaded_with_trapping = 0
       !
       do i_site = 1,n_sites_invaded
          !
          ! get site's 3d indices
          !
          call ind2ijk(nx,ny,invasion_list(i_site),is,js,ks)
          if(states(is,js,ks)/=trapped)then
             n_sites_invaded_with_trapping = n_sites_invaded_with_trapping+1
             invasion_list(n_sites_invaded_with_trapping) = invasion_list(i_site)
          endif
       enddo
       !
       n_sites_invaded = n_sites_invaded_with_trapping
       !
    endif
    !
    deallocate(mat)
    call deallocate_disjoint_set
    !
    return
    !
  !--------------------------------------
  end subroutine find_trapped_sites_cubic
  !--------------------------------------
  !
  !> @author Yder MASSON
  !> @date October 23, 2014
  !> @brief find trapped sites in arbitrary lattices
  !> @details It uses the fast a posteriori identification of trapped node 
  !> if the <b>undo_invasion=.true.</b>. 
  !> If <b>undo_invasion=.false.</b> cluster labeling is performed at each invasion steps.
  !
  !-----------------------------------------------------------
  subroutine find_trapped_sites_arbitrary(n_sites,           &
                                          states,            &
                                          offsets,           &
                                          connectivity,      &
                                          n_sites_invaded,   &
                                          invasion_list,     &
                                          undo_invasion      )
  !-----------------------------------------------------------
    !
    implicit none
    !
    integer :: n_sites ! Total number of sites in the lattice
    integer :: states(:) !< Array containing the sites's states (\b input/\b output)
    !! \n The state of trapped sites is updated to <b>state(i)=trapped</b>
    logical :: undo_invasion !< Flag: If <b>undo_invasion==.true.</b>
    !! Then use fast a posteriori method for trapping (\b input)
    integer :: n_sites_invaded ! Number of sites invaded (\b input)
    integer :: invasion_list(:) !< List of invaded sites 
    !! sorted in chronological order (\b input)
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
    ! internal variables
    !
    integer, allocatable, dimension(:) :: mat
    integer :: i,n,i_site,j,i_start,i_end,n_neighbors
    integer, dimension(:), allocatable :: lc
    integer :: label_free_clusters
    integer :: n_sites_invaded_with_trapping
    !
    ! label clusters for final state
    !
    allocate(mat(n_sites))
    ! init mat
    do i = 1,n_sites
       if(states(i)==invaded)then
          mat(i) = 0
       else
          mat(i) =-1
       endif
    enddo
    ! label cluster
    call label_clusters_arbitrary(mat,n_sites,offsets,connectivity)
    ! find free cluster's label
    label_free_clusters = get_label_free_clusters(n_sites,states,mat)
    ! update trapped sites's states
    do i = 1,n_sites
       if(mat(i)/=0.and.mat(i)/=label_free_clusters)then
          states(i) = trapped
       endif
    enddo
    !
    if(undo_invasion)then
       !
       ! undo invasion
       !
       do j = n_sites_invaded,1,-1
          !
          i_site = invasion_list(j) 
          i = offsets(i_site)
          n_neighbors = connectivity(i)
          i_start = i+1
          i_end  = i+n_neighbors
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
          if(find(mat(i_site))/=find(label_free_clusters))then
             states(i_site) = trapped
          endif
       enddo
       !
       ! remove trapped sites from invasion_list
       !
       n_sites_invaded_with_trapping = 0
       !
       do i_site = 1,n_sites_invaded
          if(states(invasion_list(i_site))/=trapped)then
             n_sites_invaded_with_trapping = n_sites_invaded_with_trapping+1
             invasion_list(n_sites_invaded_with_trapping) = invasion_list(i_site)
          endif
       enddo
       !
       n_sites_invaded = n_sites_invaded_with_trapping
       !
    endif
    !
    deallocate(mat)
    call deallocate_disjoint_set
    !
    return
    !
  !------------------------------------------
  end subroutine find_trapped_sites_arbitrary
  !------------------------------------------
  !
  !> @author Yder MASSON
  !> @date November 14, 2014
  !> @brief This functions returns the times 
  !> at which the sites got trapped.
  !> This function is for cubic lattices. 
  !
  !-------------------------------------------------------
  subroutine get_trapping_times_cubic(nx, ny, nz,        &
                                      period_x,          &
                                      period_y,          &
                                      period_z,          &
                                      states,            &
                                      n_sites_invaded,   &
                                      invasion_list,     &
                                      trapping_times     )
  !-------------------------------------------------------
    !
    implicit none
    !
    integer :: nx ! Grid dimension in the x direction (\b Input)
    integer :: ny ! Grid dimension in the x direction (\b Input)
    integer :: nz ! Grid dimension in the x direction (\b Input)
    !
    logical :: period_x !< logical flag for periodic boundaries in the x direction
    logical :: period_y !< logical flag for periodic boundaries in the y direction
    logical :: period_z !< logical flag for periodic boundaries in the z direction
    !
    integer :: states(nx,ny,nz) !< Array containing the sites's states (\b input/\b output)
    !! \n The state of trapped sites is updated to <b>state(i)=trapped</b>
    integer :: n_sites_invaded ! Number of sites invaded (\b input)
    !> List of invaded sites sorted in chronological order (\b input)
    integer :: invasion_list(:) 
    !> Array containing the times at which the sites got trapped. (\b Output)
    !> \n <b>trapping_times(i) = -1</b> at sites that have not been trapped
    integer :: trapping_times(nx,ny,nz)
    !
    ! internal variables
    ! 
    integer, dimension(:), allocatable :: labels_times
    integer, dimension(:,:,:), allocatable :: mat
    integer :: i,j,k,i_site,in1,in2,jn1,jn2,kn1,kn2,is,js,ks
    !
    !
    ! label clusters 
    !
    allocate(mat(nx,ny,nz))
    ! init mat
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             if(states(i,j,k)==trapped)then
                mat(i,j,k) =-1
             else
                mat(i,j,k) = 0
             endif
          enddo
       enddo
    enddo
    ! label cluster
    call label_clusters_cubic(mat,nx,ny,nz,period_x,period_y,period_z)
    !
    allocate(labels_times(largest_label))
    !
    labels_times(:) = -1
    !       
    do i_site = 1,n_sites_invaded
       !
       ! get site's 3d indices
       !
       call ind2ijk(nx,ny,invasion_list(i_site),is,js,ks)
       !   
       in1 = is-1
       in2 = is+1
       jn1 = js-1
       jn2 = js+1
       kn1 = ks-1
       kn2 = ks+1
       !
       ! account for periodic boundaries
       !
       if(in1<1 .and.period_x) in1 = in1+nx
       if(in2>nx.and.period_x) in2 = in2-nx
       if(jn1<1 .and.period_y) jn1 = jn1+ny
       if(jn2>ny.and.period_y) jn2 = jn2-ny
       if(kn1<1 .and.period_z) kn1 = kn1+nz
       if(kn2>nz.and.period_z) kn2 = kn2-nz
       !
       !  scan neighboring sites
       !
       if(in1>=1 .and.states(in1,js,ks)==trapped)then 
          labels_times(mat(in1,js,ks)) = i_site
       endif
       if(in2<=nx.and.states(in2,js,ks)==trapped)then
          labels_times(mat(in2,js,ks)) = i_site
       endif
       if(jn1>=1 .and.states(is,jn1,ks)==trapped)then
          labels_times(mat(is,jn1,ks)) = i_site 
       endif
       if(jn2<=ny.and.states(is,jn2,ks)==trapped)then 
          labels_times(mat(is,jn2,ks)) = i_site 
       endif
       if(kn1>=1 .and.states(is,js,kn1)==trapped)then 
          labels_times(mat(is,js,kn1)) = i_site 
       endif
       if(kn2<=nz.and.states(is,js,kn2)==trapped)then 
          labels_times(mat(is,js,kn2)) = i_site 
       endif
       !
    enddo
    !
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             if(states(i,j,k)==trapped)then
                trapping_times(i,j,k) = labels_times(mat(i,j,k))
             else
                trapping_times(i,j,k) = -1
             endif
          enddo
       enddo
    enddo
    !
    deallocate(labels_times)
    deallocate(mat)
    call deallocate_disjoint_set
    !
    return
    !
  !--------------------------------------
  end subroutine get_trapping_times_cubic
  !--------------------------------------
  !
  !> @author Yder MASSON
  !> @date November 11, 2014
  !> @brief This functions returns the times 
  !> at which the sites got trapped.
  !> This function is for arbitrary lattices. 
  !
  !-----------------------------------------------------------
  subroutine get_trapping_times_arbitrary(n_sites,           &
                                          states,            &
                                          offsets,           &
                                          connectivity,      &
                                          n_sites_invaded,   &
                                          invasion_list,     &
                                          trapping_times     )
  !-----------------------------------------------------------
    !
    implicit none
    !
    integer :: n_sites ! Total number of sites in the lattice
    integer :: states(:) !< Array containing the sites's states (\b input/\b output)
    !! \n The state of trapped sites is updated to <b>state(i)=trapped</b>
    integer :: n_sites_invaded ! Number of sites invaded (\b input)
    integer :: invasion_list(:) !< List of invaded sites 
    !! sorted in chronological order (\b input)
    !> Array containing the offsets of the data stored 
    !> in the connectivity array (\b Input)
    integer :: offsets(:)
    !> Array containing thelattice connectivity (\b input)
    !> \n For a given site \b i, with offset <b>j = offsets(i)</b> :
    !> \n <b>n=connectivity(j)</b> is the number of sites neighboring site \b i.
    !> \n <b>connectivity(j+1, j+2, ... ,j+n)</b> 
    !> contains the indices of the sites that are neighboring site \b i.
    integer :: connectivity(:)
    !> Array containing the times at which the sites got trapped. (\b Output)
    !> \n <b>trapping_times(i) = -1</b> at sites that have not been trapped
    integer :: trapping_times(:)
    !
    ! internal variables
    ! 
    integer, dimension(:), allocatable :: labels_times
    integer, allocatable, dimension(:) :: mat
    integer :: offset
    integer :: i,i_site,j,n_neighbors
    !
    !
    ! label clusters 
    !
    allocate(mat(n_sites))
    ! init mat
    do i = 1,n_sites
       if(states(i)==trapped)then
          mat(i) =-1
       else
          mat(i) = 0
       endif
    enddo
    ! label cluster
    call label_clusters_arbitrary(mat,n_sites,offsets,connectivity)
    !
    allocate(labels_times(largest_label))
    !
    labels_times(:) = -1
    !
    do j = 1,n_sites_invaded
       !
       i_site = invasion_list(j) 
       offset = offsets(i_site)
       n_neighbors = connectivity(offset)
       !
       do i = offset+1,offset+n_neighbors 
          if(states(connectivity(i))==trapped)then
             labels_times(labels(mat(connectivity(i)))) = j
          endif
       enddo
    enddo
    !
    do i_site = 1,n_sites
       if(states(i_site)==trapped)then
          trapping_times(i_site) = labels_times(mat(i_site))
       else
          trapping_times(i_site) = -1
       endif
    enddo
    !
    deallocate(labels_times)
    deallocate(mat)
    call deallocate_disjoint_set
    !
    return
    !
  !------------------------------------------
  end subroutine get_trapping_times_arbitrary
  !------------------------------------------
  !
!=========================
end module module_trapping
!=========================
