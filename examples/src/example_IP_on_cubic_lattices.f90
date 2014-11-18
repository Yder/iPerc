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
!> @date November 17, 2014
!> @brief A simple example showing how to setup and 
!> model invasion percolation on a cubic lattice
!
!-----------------------------------
program example_IP_on_cubic_lattices
!-----------------------------------
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
integer :: nx ! Grid dimension, i.e. number of sites in the x direction
integer :: ny ! Grid dimension, i.e. number of sites in the y direction
integer :: nz ! Grid dimension, i.e. number of sites in the z direction
!
real :: dx ! Grid spacing, i.e. distance between two neighbor sites in the x direction
real :: dy ! Grid spacing, i.e. distance between two neighbor sites in the y direction
real :: dz ! Grid spacing, i.e. distance between two neighbor sites in the z direction
!
logical :: period_x ! Flag for using periodic boundaries in the x direction
logical :: period_y ! Flag for using periodic boundaries in the y direction
logical :: period_z ! Flag for using periodic boundaries in the z direction
!
! Sites attributes:
!
real, allocatable :: values(:,:,:) ! Array containing the radius of the sites
integer, allocatable :: states(:,:,:) ! Array containing the state fo the sites
integer :: n_sites_invaded ! Nuber of invaded sites
integer, allocatable :: invasion_list(:) ! Indices of the invaded sites
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
integer :: i, j, k, seed, select_state
real :: pi
character(len=100) :: file_name
integer, parameter :: file_unit = 12
!
! get pi value:
! 
pi = 4.*atan(1.)
!
! 3) Setup the simulation:
! """"""""""""""""""""""""
!
! 3.a) Define the grid dimensions:
! ''''''''''''''''''''''''''''''''
!
! for 2D simulations, set one of the dimension to one, e.g. ny=1
!
nx = 100
ny = 100
nz = 200
!
! 3.b) Define the grid spacing:
! '''''''''''''''''''''''''''''
!
dx = 1.0
dy = 1.0
dz = 1.0
!
! 3.c) Setup the domain boundaries:
! ''''''''''''''''''''''''''''''''' 
!
! use period_x = .true. for periodic bondaries in the x direction
! use period_x = .false. for sealed walls in the x direction
!
period_x = .false.
period_y = .false.
period_z = .false.
!
! 3.d) Allocate memory for states, values and invasion_list arrays:
! '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!
allocate(values(nx,ny,nz))        ! values is a 3D array
allocate(states(nx,ny,nz))        ! states is a 3D array
allocate(invasion_list(nx*ny*nz)) ! invasion_list is a 1D array
!
! 3.e) Setup sites's radius:
! ''''''''''''''''''''''''''
! 
seed = 22
call srand(seed) ! init intrinsic random number generator
!
do k = 1,nz
   do j = 1,ny
      do i = 1,nx
         ! assign a random value to all sites
         values(i,j,k) = rand(0)*dx
      enddo
   enddo
enddo
!
! 3.f) Initialize sites's states:
! '''''''''''''''''''''''''''''''
!
states(:,:,:) = not_invaded ! set all sites to uninvaded state
!
! 3.g) setup injection region (bottom wall here):
! '''''''''''''''''''''''''''''''''''''''''''''''
!
! set neighboring state to all site 
! where you want to inject fluid
states(:,:,1) = neighboring 
!
! 3.h) setup exit region (top wall here):
! '''''''''''''''''''''''''''''''''''''''
!
! Set exit_site state to all sites 
! where the fluid is free to escape
states(:,:,nz) = exit_site
!
! 3.i) setup fluid properties
! '''''''''''''''''''''''''''
!
sigma = 0.0728 ! Surface tension
theta_c = pi ! Contact angle
delta_rho  = 0.2 ! Density contrast
!
! 3.j) Setup gravity
! ''''''''''''''''''
!
! acceleration of gravity vector
gx = 0.00
gy = 0.00
gz =-9.81
!
! 3.k) Choose weather or not you want
! '''''''''''''''''''''''''''''''''''
! to account for gravity:
! '''''''''''''''''''''''
!
gravity = .false. ! set gravity=.true. to acount for gravity 
!
! 3.k) Choose weather or not you want
! '''''''''''''''''''''''''''''''''''
! to account for trapping:
! ''''''''''''''''''''''''
!
trapping = .true. ! set trapping=.true. to account for trapping
!
! 4) Run the invasion percolation simulation:
! """""""""""""""""""""""""""""""""""""""""""
!
call invade_cubic_lattice_fast(nx,              &
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
call save_selected_sites_cubic_lattice(invaded,          &
                                       states,           &
                                       n_sites_invaded,  &
                                       invasion_list,    &
                                       nx,ny,nz,         &
                                       dx,dy,dz,         &
                                       trim(file_name),  &
                                       file_unit         )  
!
! 5.b) save complete info:
! ''''''''''''''''''''''''
!
! Uncomment one of the following to choose the desired file format
!
!file_name = 'example_save_cubic_lattice.vti' ! if you want a vtk file
file_name = 'example_save_sites_cubic_lattice.cvs' ! if you want a cvs file
!file_name = 'example_save_sites_cubic_lattice.dat' ! if you want an ascii file
!file_name = 'example_save_sites_cubic_lattice.bin' ! if you want a binary file
!
call save_cubic_lattice(states,          &
                        values,          &
                        n_sites_invaded, &
                        invasion_list,   &
                        nx,ny,nz,        &
                        dx,dy,dz,        &
                        period_x,        &
                        period_y,        &
                        period_z,        &
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
!---------------------------------------
end program example_IP_on_cubic_lattices
!---------------------------------------
