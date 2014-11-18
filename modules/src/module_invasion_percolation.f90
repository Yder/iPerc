!=====================================!
!                                     !
!                                     !
!     MODULE INVASION PERCOLATION     !
!                                     ! 
!                                     !
!=====================================!
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
!=================================
module module_invasion_percolation
!=================================
  !
  use module_trapping
  use module_gravity
  use module_cubic_indices
  use module_invasion_percolation_constants
  use module_disjoint_set
  use module_binary_tree
  use module_write_output_files
  use module_random_media
  !
contains
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! functions for cubic lattices !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !> @author Yder MASSON
  !> @date October 21, 2014
  !> @brief A simple but inefficient implementation of invasion percolation
  !> on a 3D cubic lattice
  !> @details
  !
  !--------------------------------------------------------
  subroutine invade_cubic_lattice_simple(nx,              &
                                         ny,              &
                                         nz,              &
                                         dx,              &
                                         dy,              &
                                         dz,              &
                                         period_x,        &
                                         period_y,        &
                                         period_z,        &
                                         values,          &
                                         states,          &
                                         n_sites_invaded, &
                                         invasion_list,   &
                                         gravity,         &
                                         trapping,        &
                                         sigma,           &
                                         theta_c,         &
                                         delta_rho,       &
                                         gx,              &
                                         gy,              &
                                         gz               )
  !--------------------------------------------------------
    !
    implicit none
    !
    integer :: nx !< Grid dimension in the x direction (\b input)
    integer :: ny !< Grid dimension in the y direction (\b input)
    integer :: nz !< Grid dimension in the z direction (\b input)
    !
    real :: dx !<  Grid spacing in the x direction \f$ \Delta x\f$ (\b input)
    real :: dy !<  Grid spacing in the x direction \f$ \Delta y\f$ (\b input)
    real :: dz !<  Grid spacing in the x direction \f$ \Delta z\f$ (\b input)
    !
    logical :: period_x !< Flag for periodic boundaries in the x direction (\b input)
    logical :: period_y !< Flag for periodic boundaries in the y direction (\b input)
    logical :: period_z !< Flag for periodic boundaries in the z direction (\b input)
    !
    !> Array containing the sites's sizes \f$ a_i \f$ (\b input/\b output)
    !> \n When <b>gravity==.true.</b> this array is modified 
    !> and contains the invasion potentials \f$ P_i \f$ at output time
    real :: values(nx,ny,nz) 
    !
    !> Array containing the sites's states (\b input/\b output)
    !> \n Set <b>state(i)=neighboring</b> at sites \b i 
    !> where you would like to inject the fluid
    !> \n Set <b>state(i)=exit_site</b> at sites \b i 
    !> where the defending fluid can escape 
    !>(the simulation will stop when the invading fluid percolates, 
    !> i.e. it reaches one of these sites)
    !> \n Set <b>state(i)=sealed</b> to prevent sites \b i from being invaded
    !> \n At output time, we have <b>state(i)=trapped</b> at trapped sites
    !> \n At output time, we have <b>state(i)=invaded</b> at invaded sites
    integer :: states(nx,ny,nz)
    !
    integer :: n_sites_invaded !< Number of sites invaded (\b output)
    !
    !> Array containing the list of sites invaded 
    !> sorted in chronological order (\b output)
    integer :: invasion_list(:)
    !
    logical :: gravity !< Flag for gravity (\b input)
    !! \n Set <b>gravity=.true.</b> to account for gravity
    !! \n Set <b>gravity=.false.</b> to ignore gravity
    !
    logical :: trapping !< Flag for trapping (\b input)
    !! \n Set <b>trapping=.true.</b> to account for trapping
    !! \n Set <b>trapping=.false.</b> to ignore trapping
    !
    real :: sigma !< Surface tension \f$ \sigma \f$ (\b input)
    real :: theta_c !< Equilibrium contact angle \f$ \theta_c \f$ (\b input)
    real :: delta_rho !< Fluid density contrast \f$ \Delta \rho \f$ (\b input)
    !
    real :: gx !< Acceleration of gravity \f$ g \f$ in the x direction (\b input)
    real :: gy !< Acceleration of gravity \f$ g \f$ in the y direction (\b input)
    real :: gz !< Acceleration of gravity \f$ g \f$ in the z direction (\b input)
    !
    ! internal variables :
    !
    integer :: i_site,j_site,k_site
    integer :: i_minus_one,i_plus_one  
    integer :: j_minus_one,j_plus_one  
    integer :: k_minus_one,k_plus_one  
    real    :: min_value 
    integer :: i,j,k 
    logical, parameter :: undo_invasion = .false.
    !
    ! Add gravity if desired
    !
    if(gravity)then
       !
       call add_gravity_cubic_lattice(values,    &
                                      nx,ny,nz,  &
                                      dx,dy,dz,  &
                                      sigma,     &
                                      theta_c,   &
                                      delta_rho, &
                                      gx,gy,gz)
       !
    endif
    !
    ! Set number of invaded sites to zero
    !
    n_sites_invaded = 0
    !
    ! Invasion percolation loop starts here
    !
    do
       !
       ! Find next site to invade
       !
       ! Set min_value to huge number
       !
       min_value = huge(min_value)
       !
       ! Find the index of the neighboring site
       !
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                if(states(i,j,k)==neighboring.and.values(i,j,k)<min_value)then
                   i_site = i
                   j_site = j
                   k_site = k
                   min_value = values(i_site,j_site,k_site)
                endif
             enddo
          enddo
       enddo
       !
       ! Invade site 
       !
       n_sites_invaded = n_sites_invaded+1
       states(i_site,j_site,k_site) = invaded
       invasion_list(n_sites_invaded) = ijk2ind(nx,ny,i_site,j_site,k_site)
       !
       ! Update neighboring sites's states
       !
       !
       ! Neighboring sites indices 
       !
       i_minus_one = i_site - 1
       i_plus_one  = i_site + 1
       j_minus_one = j_site - 1
       j_plus_one  = j_site + 1
       k_minus_one = k_site - 1
       k_plus_one  = k_site + 1
       !
       ! Account for periodic boundaries
       !
       if(i_minus_one<1.and.period_x) i_minus_one = i_minus_one + nx
       if(i_plus_one>nx.and.period_x) i_plus_one  = i_plus_one  - nx
       if(j_minus_one<1.and.period_y) j_minus_one = j_minus_one + ny
       if(j_plus_one>ny.and.period_y) j_plus_one  = j_plus_one  - ny
       if(k_minus_one<1.and.period_z) k_minus_one = k_minus_one + nz
       if(k_plus_one>nz.and.period_z) k_plus_one  = k_plus_one  - nz
       !
       ! Exit when cluster percolates
       !
       if(i_minus_one>=1 .and.states(i_minus_one,j_site,k_site)==exit_site) exit
       if(i_plus_one<=nx.and.states(i_plus_one,j_site,k_site)==exit_site) exit
       if(j_minus_one>=1 .and.states(i_site,j_minus_one,k_site)==exit_site) exit
       if(j_plus_one<=ny.and.states(i_site,j_plus_one,k_site)==exit_site) exit
       if(k_minus_one>=1 .and.states(i_site,j_site,k_minus_one)==exit_site) exit
       if(k_plus_one<=nz.and.states(i_site,j_site,k_plus_one)==exit_site) exit
       !
       ! Update neighboring sites's states
       !
       if(i_minus_one>=1 .and. states(i_minus_one,j_site,k_site)==not_invaded)&
            states(i_minus_one,j_site,k_site)=neighboring
       if(i_plus_one<=nx .and. states(i_plus_one,j_site,k_site)==not_invaded)&
            states(i_plus_one,j_site,k_site)=neighboring
       if(j_minus_one>=1 .and. states(i_site,j_minus_one,k_site)==not_invaded)&
            states(i_site,j_minus_one,k_site)=neighboring
       if(j_plus_one<=ny .and. states(i_site,j_plus_one,k_site)==not_invaded)&
            states(i_site,j_plus_one,k_site)=neighboring
       if(k_minus_one>=1 .and. states(i_site,j_site,k_minus_one)==not_invaded)&
            states(i_site,j_site,k_minus_one)=neighboring
       if(k_plus_one<=nz .and. states(i_site,j_site,k_plus_one)==not_invaded)&
            states(i_site,j_site,k_plus_one)=neighboring 
       !
       ! Account for trapping 
       !
       if(trapping)then
          call find_trapped_sites_cubic(nx,                &
                                        ny,                &
                                        nz,                &
                                        states,            &
                                        period_x,          &
                                        period_y,          &
                                        period_z,          &
                                        n_sites_invaded,   &
                                        invasion_list,     &
                                        undo_invasion      )
       endif
       !
    end do ! Invasion percolation loop ends here
    !
    return
    !
  !-----------------------------------------
  end subroutine invade_cubic_lattice_simple
  !-----------------------------------------
  !
  !> @author Yder MASSON
  !> @date October 21, 2014
  !> @brief An efficient implementation of invasion percolation 
  !> on a 3D cubic lattice
  !> @details
  !
  !------------------------------------------------------
  subroutine invade_cubic_lattice_fast(nx,              &
                                       ny,              &
                                       nz,              &
                                       dx,              &
                                       dy,              &
                                       dz,              &
                                       period_x,        &
                                       period_y,        &
                                       period_z,        &
                                       values,          &
                                       states,          &
                                       n_sites_invaded, &
                                       invasion_list,   &
                                       gravity,         &
                                       trapping,        &
                                       sigma,           &
                                       theta_c,         &
                                       delta_rho,       &
                                       gx,              &
                                       gy,              &
                                       gz               )
  !------------------------------------------------------
    !
    implicit none
    !
    integer :: nx !< Grid dimension in the x direction (\b input)
    integer :: ny !< Grid dimension in the y direction (\b input)
    integer :: nz !< Grid dimension in the z direction (\b input)
    !
    real :: dx !<  Grid spacing in the x direction \f$ \Delta x\f$ (\b input)
    real :: dy !<  Grid spacing in the x direction \f$ \Delta y\f$ (\b input)
    real :: dz !<  Grid spacing in the x direction \f$ \Delta z\f$ (\b input)
    !
    logical :: period_x !< Flag for periodic boundaries in the x direction (\b input)
    logical :: period_y !< Flag for periodic boundaries in the y direction (\b input)
    logical :: period_z !< Flag for periodic boundaries in the z direction (\b input)
    !
    !> Array containing the sites's sizes \f$ a_i \f$ (\b input/\b output)
    !> \n When <b>gravity==.true.</b> this array is modified 
    !> and contains the invasion potentials \f$ P_i \f$ at output time
    real :: values(nx*ny*nz) ! declared as 1D array for use with the binary tree module
    !
    !> Array containing the sites's states (\b input/\b output)
    !> \n Set <b>state(i)=neighboring</b> at sites \b i 
    !> where you would like to inject the fluid
    !> \n Set <b>state(i)=exit_site</b> at sites \b i 
    !> where the defending fluid can escape 
    !>(the simulation will stop when the invading fluid percolates, 
    !> i.e. it reaches one of these sites)
    !> \n Set <b>state(i)=sealed</b> to prevent sites \b i from being invaded
    !> \n At output time, we have <b>state(i)=trapped</b> at trapped sites
    !> \n At output time, we have <b>state(i)=invaded</b> at invaded sites
    integer :: states(nx,ny,nz)
    !
    integer :: n_sites_invaded !< Number of sites invaded (\b output)
    !
    !> Array containing the list of sites invaded 
    !> sorted in chronological order (\b output)
    integer :: invasion_list(:)
    !
    logical :: gravity !< Flag for gravity (\b input)
    !! \n Set <b>gravity=.true.</b> to account for gravity
    !! \n Set <b>gravity=.false.</b> to ignore gravity
    !
    logical :: trapping !< Flag for trapping (\b input)
    !! \n Set <b>trapping=.true.</b> to account for trapping
    !! \n Set <b>trapping=.false.</b> to ignore trapping
    !
    real :: sigma !< Surface tension \f$ \sigma \f$ (\b input)
    real :: theta_c !< Equilibrium contact angle \f$ \theta_c \f$ (\b input)
    real :: delta_rho !< Fluid density contrast \f$ \Delta \rho \f$ (\b input)
    !
    real :: gx !< Acceleration of gravity \f$ g \f$ in the x direction (\b input)
    real :: gy !< Acceleration of gravity \f$ g \f$ in the y direction (\b input)
    real :: gz !< Acceleration of gravity \f$ g \f$ in the z direction (\b input)
    !
    ! internal variables
    !
    integer :: i_site,j_site,k_site 
    integer :: i_minus_one,i_plus_one 
    integer :: j_minus_one,j_plus_one  
    integer :: k_minus_one,k_plus_one  
    integer :: i,j,k 
    logical, parameter :: undo_invasion = .true.
    !
    ! Add gravity if desired
    !
    if(gravity)then
       !
       call add_gravity_cubic_lattice(values,    &
                                      nx,ny,nz,  &
                                      dx,dy,dz,  &
                                      sigma,     &
                                      theta_c,   &
                                      delta_rho, &
                                      gx,gy,gz)
       !
    endif
    !
    ! Initialize binary tree
    ! 
    allocate(tree(1))
    treedim = 0
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             if(states(i,j,k)==neighboring)then
                call add_branch(values,ijk2ind(nx,ny,i,j,k))
             endif
          enddo
       enddo
    enddo
    !
    ! Set number of invaded sites to zero
    !
    n_sites_invaded = 0
    !
    ! Invasion percolation loop starts here
    !
    do
       !
       ! Find next site to invade
       !
       i = tree(1)
       !
       call ind2ijk(nx,ny,i,i_site,j_site,k_site)
       call update_tree_root(values)
       !
       ! Invade site
       !
       n_sites_invaded = n_sites_invaded+1
       states(i_site,j_site,k_site) = invaded
       invasion_list(n_sites_invaded) = ijk2ind(nx,ny,i_site,j_site,k_site)
       !
       ! Update neighboring sites's states 
       !
       ! Neighboring sites indices 
       !
       i_minus_one = i_site - 1
       i_plus_one  = i_site + 1
       j_minus_one = j_site - 1
       j_plus_one  = j_site + 1
       k_minus_one = k_site - 1
       k_plus_one  = k_site + 1
       !
       ! Account for periodic boundaries
       !
       if(i_minus_one<1.and.period_x) i_minus_one = i_minus_one + nx
       if(i_plus_one>nx.and.period_x) i_plus_one  = i_plus_one  - nx
       if(j_minus_one<1.and.period_y) j_minus_one = j_minus_one + ny
       if(j_plus_one>ny.and.period_y) j_plus_one  = j_plus_one  - ny
       if(k_minus_one<1.and.period_z) k_minus_one = k_minus_one + nz
       if(k_plus_one>nz.and.period_z) k_plus_one  = k_plus_one  - nz
       !
       ! Exit when cluster percolates
       !
       if(i_minus_one>=1 .and.states(i_minus_one,j_site,k_site)==exit_site) exit
       if(i_plus_one<=nx.and.states(i_plus_one,j_site,k_site)==exit_site) exit
       if(j_minus_one>=1 .and.states(i_site,j_minus_one,k_site)==exit_site) exit
       if(j_plus_one<=ny.and.states(i_site,j_plus_one,k_site)==exit_site) exit
       if(k_minus_one>=1 .and.states(i_site,j_site,k_minus_one)==exit_site) exit
       if(k_plus_one<=nz.and.states(i_site,j_site,k_plus_one)==exit_site) exit
       !
       ! Update neighboring sites's states
       !
       if(i_minus_one>=1 .and. states(i_minus_one,j_site,k_site)==not_invaded)then
          states(i_minus_one,j_site,k_site)=neighboring
          call add_branch(values,ijk2ind(nx,ny,i_minus_one,j_site,k_site))
       endif
       if(i_plus_one<=nx .and. states(i_plus_one,j_site,k_site)==not_invaded)then
          states(i_plus_one,j_site,k_site)=neighboring                    
          call add_branch(values,ijk2ind(nx,ny,i_plus_one,j_site,k_site))
       endif
       if(j_minus_one>=1 .and. states(i_site,j_minus_one,k_site)==not_invaded)then
          states(i_site,j_minus_one,k_site)=neighboring                    
          call add_branch(values,ijk2ind(nx,ny,i_site,j_minus_one,k_site))
       endif
       if(j_plus_one<=ny .and. states(i_site,j_plus_one,k_site)==not_invaded)then
          states(i_site,j_plus_one,k_site)=neighboring                    
          call add_branch(values,ijk2ind(nx,ny,i_site,j_plus_one,k_site))
       endif
       if(k_minus_one>=1 .and. states(i_site,j_site,k_minus_one)==not_invaded)then
          states(i_site,j_site,k_minus_one)=neighboring                    
          call add_branch(values,ijk2ind(nx,ny,i_site,j_site,k_minus_one))
       endif
       if(k_plus_one<=nz .and. states(i_site,j_site,k_plus_one)==not_invaded)then
          states(i_site,j_site,k_plus_one)=neighboring                    
          call add_branch(values,ijk2ind(nx,ny,i_site,j_site,k_plus_one))
       endif
       !
    end do ! invasion percolation loop ends here
    !
    ! Account for trapping a posteriori
    !
    if(trapping)then
       call find_trapped_sites_cubic(nx,                &
                                     ny,                &
                                     nz,                &
                                     states,            &
                                     period_x,          &
                                     period_y,          &
                                     period_z,          &
                                     n_sites_invaded,   &
                                     invasion_list,     &
                                     undo_invasion      )
    endif
    !
    ! Deallocate tree
    !
    call deallocate_binary_tree
    !
    return
    !
  !---------------------------------------
  end subroutine invade_cubic_lattice_fast
  !---------------------------------------
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! functions for arbitrary lattices !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !> @author Yder MASSON
  !> @date October 21, 2014
  !> @brief A simple but inefficient implementation of invasion percolation 
  !> on arbitrary lattices
  !> @details
  !
  !------------------------------------------------------------
  subroutine invade_arbitrary_lattice_simple(n_sites,         &
                                             x,               &
                                             y,               &
                                             z,               &
                                             offsets,         &
                                             connectivity,    &
                                             values,          &
                                             states,          &
                                             n_sites_invaded, &
                                             invasion_list,   &
                                             gravity,         &
                                             trapping,        &
                                             sigma,           &
                                             theta_c,         &
                                             delta_rho,       &
                                             gx,              &
                                             gy,              &
                                             gz               )
  !------------------------------------------------------------
    !
    implicit none
    !
    integer :: n_sites !< Total number of sites in the lattice (\b Input)
    !
    real,dimension(:) :: x !< Array containing the x coordinates of the sites (\b Input)
    real,dimension(:) :: y !< Array containing the y coordinates of the sites (\b Input)
    real,dimension(:) :: z !< Array containing the z coordinates of the sites (\b Input)
    !
    !> Array containing the offsets of the data stored 
    !> in the connectivity array (\b Input)
    integer :: offsets(:)
    !
    !> Array containing thelattice connectivity (\b input)
    !> \n For a given site \b i, with offset <b>j = offsets(i)</b> :
    !> \n <b>n=connectivity(j)</b> is the number of sites neighboring site \b i.
    !> \n <b>connectivity(j+1, j+2, ... ,j+n)</b> 
    !> contains the indices of the sites that are neighboring site \b i.
    integer :: connectivity(:)
    !
    !> Array containing the sites's sizes \f$ a_i \f$ (\b input/\b output)
    !> \n When <b>gravity==.true.</b> this array is modified 
    !> and contains the invasion potentials \f$ P_i \f$ at output time
    !
    real, dimension(:) :: values 
    !
    !> Array containing the sites's states (\b input/\b output)
    !> \n Set <b>state(i)=neighboring</b> at sites \b i 
    !> where you would like to inject the fluid
    !> \n Set <b>state(i)=exit_site</b> at sites \b i 
    !> where the defending fluid can escape 
    !>(the simulation will stop when the invading fluid percolates, 
    !> i.e. it reaches one of these sites)
    !> \n Set <b>state(i)=sealed</b> to prevent sites \b i from being invaded
    !> \n At output time, we have <b>state(i)=trapped</b> at trapped sites
    !> \n At output time, we have <b>state(i)=invaded</b> at invaded sites
    integer, dimension(:) :: states
    !
    integer :: n_sites_invaded !< Number of sites invaded (\b output)
    !
    !> Array containing the list of sites invaded 
    !> sorted in chronological order (\b output)
    integer :: invasion_list(:)
    !
    logical :: gravity !< Flag for gravity (\b input)
    !! \n Set <b>gravity=.true.</b> to account for gravity
    !! \n Set <b>gravity=.false.</b> to ignore gravity
    !
    logical :: trapping !< Flag for trapping (\b input)
    !! \n Set <b>trapping=.true.</b> to account for trapping
    !! \n Set <b>trapping=.false.</b> to ignore trapping
    !
    real :: sigma !< Surface tension \f$ \sigma \f$ (\b input)
    real :: theta_c !< Equilibrium contact angle \f$ \theta_c \f$ (\b input)
    real :: delta_rho !< Fluid density contrast \f$ \Delta \rho \f$ (\b input)
    !
    real :: gx !< Acceleration of gravity \f$ g \f$ in the x direction (\b input)
    real :: gy !< Acceleration of gravity \f$ g \f$ in the y direction (\b input)
    real :: gz !< Acceleration of gravity \f$ g \f$ in the z direction (\b input)
    !
    ! internal variables
    !
    integer :: i_site
    integer :: i, i_start,i_end, i_neighbor, n_neighbors
    real    :: min_value
    logical, parameter :: undo_invasion = .false.
    !
    ! Add gravity if desired
    !
    if(gravity)then
       !
       call add_gravity_arbitrary_lattice(values,    &
                                          n_sites,   &
                                          x,         &
                                          y,         &
                                          z,         &
                                          sigma,     &
                                          theta_c,   &
                                          delta_rho, &
                                          gx,        &
                                          gy,        &
                                          gz         )
       !
    endif
    !
    ! Set number of invaded sites to zero
    !
    n_sites_invaded = 0
    !
    ! Invasion percolation loop starts here
    !
    do
       !
       ! Find next site to invade
       !
       ! Set min_value to huge number
       !
       min_value = huge(min_value)
       !
       ! Find the index of the neighboring site
       !
       do i = 1,n_sites
          if(states(i)==neighboring.and.values(i)<min_value)then
             i_site = i
             min_value = values(i_site)
          endif
       enddo
       !
       ! Invade site 
       !
       n_sites_invaded = n_sites_invaded+1
       states(i_site) = invaded
       invasion_list(n_sites_invaded) = i_site
       !
       ! Update neighboring sites's states 
       !
       i_neighbor = offsets(i_site)
       n_neighbors = connectivity(i_neighbor)
       i_start = i_neighbor + 1
       i_end   = i_neighbor + n_neighbors
       !
       ! For all sites neighboring the site invaded (i_site)
       !
       do i = i_start,i_end
          !
          ! Exit when cluster percolates
          !
          if(states(connectivity(i))==exit_site)then
             return
          endif
          !
          ! Update neighbor site's state
          !
          if(states(connectivity(i))==not_invaded)then
             states(connectivity(i))=neighboring
          endif
          !
       enddo
       !
       ! Account for trapping 
       !
       if(trapping)then
          call find_trapped_sites_arbitrary(n_sites,           &
                                            states,            &
                                            offsets,           &
                                            connectivity,      &
                                            n_sites_invaded,   &
                                            invasion_list,     &
                                            undo_invasion      )
       endif
       !
    end do ! Invasion percolation loop ends here
    !
    return
    !
  !---------------------------------------------
  end subroutine invade_arbitrary_lattice_simple
  !---------------------------------------------
  !
  !> @author Yder MASSON
  !> @date October 21, 2014
  !> @brief An efficient implementation of invasion percolation 
  !> on arbitrary lattices.
  !> @details
  !
  !----------------------------------------------------------
  subroutine invade_arbitrary_lattice_fast(n_sites,         &
                                           x,               &
                                           y,               &
                                           z,               &
                                           offsets,         &
                                           connectivity,    &
                                           values,          &
                                           states,          &
                                           n_sites_invaded, &
                                           invasion_list,   &
                                           gravity,         &
                                           trapping,        &
                                           sigma,           &
                                           theta_c,         &
                                           delta_rho,       &
                                           gx,              &
                                           gy,              &
                                           gz               )
  !----------------------------------------------------------
    !
    implicit none
    !
    integer :: n_sites !< Total number of sites in the lattice (\b Input)
    !
    real,dimension(:) :: x !< Array containing the x coordinates of the sites (\b Input)
    real,dimension(:) :: y !< Array containing the y coordinates of the sites (\b Input)
    real,dimension(:) :: z !< Array containing the z coordinates of the sites (\b Input)
    !
    !> Array containing the offsets of the data stored 
    !> in the connectivity array (\b Input)
    integer :: offsets(:)
    !
    !> Array containing thelattice connectivity (\b input)
    !> \n For a given site \b i, with offset <b>j = offsets(i)</b> :
    !> \n <b>n=connectivity(j)</b> is the number of sites neighboring site \b i.
    !> \n <b>connectivity(j+1, j+2, ... ,j+n)</b> 
    !> contains the indices of the sites that are neighboring site \b i.
    integer :: connectivity(:)
    !
    !> Array containing the sites's sizes \f$ a_i \f$ (\b input/\b output)
    !> \n When <b>gravity==.true.</b> this array is modified 
    !> and contains the invasion potentials \f$ P_i \f$ at output time
    !
    real, dimension(:) :: values 
    !
    !> Array containing the sites's states (\b input/\b output)
    !> \n Set <b>state(i)=neighboring</b> at sites \b i 
    !> where you would like to inject the fluid
    !> \n Set <b>state(i)=exit_site</b> at sites \b i 
    !> where the defending fluid can escape 
    !>(the simulation will stop when the invading fluid percolates, 
    !> i.e. it reaches one of these sites)
    !> \n Set <b>state(i)=sealed</b> to prevent sites \b i from being invaded
    !> \n At output time, we have <b>state(i)=trapped</b> at trapped sites
    !> \n At output time, we have <b>state(i)=invaded</b> at invaded sites
    integer, dimension(:) :: states
    !
    integer :: n_sites_invaded !< Number of sites invaded (\b output)
    !
    !> Array containing the list of sites invaded 
    !> sorted in chronological order (\b output)
    integer :: invasion_list(:)
    !
    logical :: gravity !< Flag for gravity (\b input)
    !! \n Set <b>gravity=.true.</b> to account for gravity
    !! \n Set <b>gravity=.false.</b> to ignore gravity
    !
    logical :: trapping !< Flag for trapping (\b input)
    !! \n Set <b>trapping=.true.</b> to account for trapping
    !! \n Set <b>trapping=.false.</b> to ignore trapping
    !
    real :: sigma !< Surface tension \f$ \sigma \f$ (\b input)
    real :: theta_c !< Equilibrium contact angle \f$ \theta_c \f$ (\b input)
    real :: delta_rho !< Fluid density contrast \f$ \Delta \rho \f$ (\b input)
    !
    real :: gx !< Acceleration of gravity \f$ g \f$ in the x direction (\b input)
    real :: gy !< Acceleration of gravity \f$ g \f$ in the y direction (\b input)
    real :: gz !< Acceleration of gravity \f$ g \f$ in the z direction (\b input)
    !
    ! Internal variables
    !
    integer :: i_site
    integer :: i, i_start,i_end, i_neighbor, n_neighbors
    logical, parameter :: undo_invasion = .true.
    logical :: percolation
    !
    ! Add gravity if desired
    !
    if(gravity)then
       !
       call add_gravity_arbitrary_lattice(values,    &
                                          n_sites,   &
                                          x,         &
                                          y,         &
                                          z,         &
                                          sigma,     &
                                          theta_c,   &
                                          delta_rho, &
                                          gx,        &
                                          gy,        &
                                          gz         )
       !
    endif
    !
    ! Initialize binary tree
    ! 
    allocate(tree(1))
    treedim = 0
    do i = 1,n_sites
       if(states(i)==neighboring)call add_branch(values,i)
    enddo
    !
    ! Set number of invaded sites to zero
    !
    n_sites_invaded = 0
    !
    ! Init percolation flag
    !
    percolation = .false.
    !
    ! Invasion percolation loop starts here
    !
    do
       !
       ! Find next site to invade
       !
       i_site = tree(1)
       call update_tree_root(values)
       !
       ! Invade site
       !
       n_sites_invaded = n_sites_invaded+1
       states(i_site) = invaded
       invasion_list(n_sites_invaded) = i_site
       !
       ! update neighboring sites's states
       !
       i_neighbor = offsets(i_site)
       n_neighbors = connectivity(i_neighbor)
       i_start = i_neighbor + 1
       i_end   = i_neighbor + n_neighbors
       !
       ! For all sites neighboring the site invaded (i_site)
       !
       do i = i_start,i_end
          !
          ! Exit if cluster percolates
          !
          if(states(connectivity(i))==exit_site)then 
             !
             ! Account for trapping a posteriori
             !
             if(trapping)then
                call find_trapped_sites_arbitrary(n_sites,           &
                                                  states,            &
                                                  offsets,           &
                                                  connectivity,      &
                                                  n_sites_invaded,   &
                                                  invasion_list,     &
                                                  undo_invasion      )
             endif
             !
             return
             ! 
          endif
          !
          ! Update neighbor site's state
          !
          if(states(connectivity(i))==not_invaded)then
             states(connectivity(i))=neighboring
             call add_branch(values,connectivity(i))
          endif
          !
       enddo
       !
    end do ! Invasion percolation loop ends here
    !
    ! Deallocate tree
    !
    call deallocate_binary_tree
    !
    return
    !
  !-------------------------------------------
  end subroutine invade_arbitrary_lattice_fast
  !-------------------------------------------
  !
!=====================================
end module module_invasion_percolation
!=====================================

