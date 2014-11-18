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
!> @date November 18, 2014
!> @brief A simple example showing how to model invason percolation on arbitrary lattices.
!> \n This example permits to model invasion percolation on 2D hexagonal lattics and 3D HCP lattices
!---------------------------------------
program example_IP_on_arbitrary_lattices
!---------------------------------------
!
! 1) First load the invasion percolation module:
! """"""""""""""""""""""""""""""""""""""""""""""
!
use module_invasion_percolation
!
implicit none
!
! 2) Declare the following variables that are needed 
! """"""""""""""""""""""""""""""""""""""""""""""""""
! for modeling IP on a cubic lattice :
! """"""""""""""""""""""""""""""""""""
!
! Lattice variables:
!
integer :: n_sites
real, allocatable :: x(:),y(:),z(:)
integer, allocatable :: offsets(:)
integer, allocatable ::connectivity(:)
!
! Sites attributes:
!
real, allocatable :: values(:) 
integer, allocatable :: states(:)  
integer :: n_sites_invaded
integer, allocatable :: invasion_list(:)
!
! Trapping:
!
logical :: trapping ! Flag to account for trapping
! 
! Fluid variables: 
!
real :: sigma ! surface tension
real :: theta_c ! contact angle
real :: delta_rho ! density contrast
!
! Gravity variables:
!
logical :: gravity ! Flag to account for gravity
real :: gx ! acceleration of gravity (x component)
real :: gy ! acceleration of gravity (y component)
real :: gz ! acceleration of gravity (z component)
!
! These are dummy variables to setup the simulation: 
!
integer :: nx
integer :: ny
integer :: nz
!
integer :: i,j,k,ii,jj,kk,n,i_site,i_neighbor,offset
integer, allocatable :: neighbors(:)
!
real :: pi,dh,d
integer :: select_state
character(len=100) :: file_name
integer, parameter :: file_unit = 12
!
! get pi value
!
pi = 4.*atan(1.)
!
! 3) Setup the simulation:
! """"""""""""""""""""""""
!
! Let's first build and HCP lattice (hexagonal close packing)
! we need to construct the connectivity, offsets, and coordinates arrays
! that defines the lattice
!
! grid dimensions: 
! (these are only used to build the
! connectivity, offsets, and coordinates arrays
! they are not used for the simulation)
!
! To run a 2D simulation on an hexagonal 
! lattice, use ny = 1 here 
nx = 10*3
ny = 15*3
nz = 20*3
!
! Get the total number of sites 
! (needed for arbitrary lattices)
!
n_sites = nx*ny*nz
!
! Set the grid spacing:
! (it is only used to build the
! connectivity, offsets, and coordinates arrays
! is is not used for the simulation)
dh = 1.0
!
! allocate memory :
!
allocate(x(n_sites))
allocate(y(n_sites))
allocate(z(n_sites))
allocate(offsets(n_sites))
allocate(connectivity(n_sites*(12+1))) ! Each site have at most 12 neighbors on an HCP lattice
allocate(values(n_sites))
allocate(states(n_sites))
allocate(invasion_list(n_sites))
allocate(neighbors(12)) ! dummy array not used in the simulation
!
! setup sites's states
!
i_site = 0
!
do k = 1,nz
   do j = 1,ny
      do i = 1,nx
         !
         i_site = i_site+1
         !
         x(i_site) = (2*i+mod(j,2)+mod(k,2))*dh/2.
         y(i_site) = 2.*sqrt(6.)/3.*j*dh/2.
         z(i_site) = sqrt(3.)*(k+mod(j,2)/3.)*dh/2.
         !
         ! init state
         !
         states(i_site) = not_invaded
         !
         ! setup injection region (bottom wall)
         !
         if(k==1)states(i_site) = neighboring
         !
         ! setup exit_region (top wall)
         !
         if(k==nz)states(i_site) = exit_site
         !
      enddo
   enddo
enddo
!
! find sites's neighbors
!
i_site = 0
offset = 1
!
do k = 1,nz
   do j = 1,ny
      do i = 1,nx
         !
         i_site = i_site+1
         !
         n = 0
         !
         do kk = max(1,k-1),min(nz,k+1)
            do jj = max(1,j-1),min(ny,j+1)
               do ii = max(1,i-1),min(nx,i+1)
                  !
                  i_neighbor = ijk2ind(nx,ny,ii,jj,kk)
                  !
                  d = ( (x(i_site)-x(i_neighbor))**2 &
                       +(y(i_site)-y(i_neighbor))**2 &
                       +(z(i_site)-z(i_neighbor))**2 )
                  !
                  if(d<dh*1.01.and.i_neighbor/=i_site)then
                     n = n+1
                     if(n>12)then
                        print*,'error : n should be less than 12 for hpc lattices'
                        stop
                     endif
                     neighbors(n) = i_neighbor
                  endif
               enddo
            enddo
         enddo
         !
         ! fill the offsets and connectivity arrays 
         ! needed for the simulations
         !
         offsets(i_site) = offset
         connectivity(offset) = n
         connectivity(offset+1:offset+n) = neighbors(1:n)
         offset = offset+n+1
         !
      enddo
   enddo
enddo
!
! setup fluid properties
!
sigma = 0.0728
theta_c = pi
delta_rho  = 0.2
!
! setup gravity
!
gx = 0.00
gy = 0.00
gz =-9.81
!
gravity = .false.
!
!setup trapping
!
trapping = .true.
!
! set values
!
do i_site = 1,n_sites
   values(i_site) = rand()
enddo
!
! 4) Run the invasion percolation simulation:
! """""""""""""""""""""""""""""""""""""""""""
!
call invade_arbitrary_lattice_fast(n_sites,         &
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
!
! 5) save result to file:
! """""""""""""""""""""""
!
! 5.a) save sites based on their state:
! '''''''''''''''''''''''''''''''''''''
!
! Uncomment one of the following:
!
! select_state = trapped ! uncomment to save trapped sites only
 select_state = invaded ! uncomment to save invaded sites only
! select_state = exit_site ! uncomment to save exit_sites sites only
! select_state = not_invaded ! uncomment to save invaded_sites sites only
!
! Uncomment one of the following to choose the desired file format
!
!file_name = 'example_selected_sites_cubic_lattice.vtp' ! if you want a vtk file
file_name = 'example_selected_sites_cubic_lattice.cvs' ! if you want a cvs file
!file_name = 'example_selected_sites_cubic_lattice.dat' ! if you want an ascii file
!file_name = 'example_selected_sites_cubic_lattice.bin' ! if you want a binary file
! 
call save_selected_sites_arbitrary_lattice(select_state,    &
                                           n_sites,         &
                                           states,          &
                                           n_sites_invaded, &
                                           invasion_list,   &
                                           x,y,z,           &
                                           file_name,       &
                                           file_unit        )
!
! 5.b) save complete info:
! ''''''''''''''''''''''''
!
! Uncomment one of the following to choose the desired file format
!
!file_name = 'example_save_cubic_lattice.vtp' ! if you want a vtk file
file_name = 'example_save_sites_cubic_lattice.cvs' ! if you want a cvs file
!file_name = 'example_save_sites_cubic_lattice.dat' ! if you want an ascii file
!file_name = 'example_save_sites_cubic_lattice.bin' ! if you want a binary file
!
call save_arbitrary_lattice(n_sites,         &
                            states,          &
                            values,          &
                            n_sites_invaded, &
                            invasion_list,   &
                            offsets,         &
                            connectivity,    &
                            x,y,z,           &
                            trim(file_name), &
                            file_unit        )  
!
! 6) Clean up:
! """"""""""""
!
deallocate(values,states,invasion_list)
!
! Done !
!
stop
!
!-------------------------------------------
end program example_IP_on_arbitrary_lattices
!-------------------------------------------
