!==============================!
!                              !
!                              !
!     MODULE INTERPOLATION     !
!                              !
!                              !
!==============================!
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
!==========================
module module_interpolation
!==========================
!
contains
!
! @author Yder MASSON
! @date November 12, 2014
! @brief Preforms 3D trilinear interpolation 
!
!-------------------------------------------------------------------------------------
real function trilinear_iterpolation(mat,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
!-------------------------------------------------------------------------------------
!
real, dimension(0:nx-1,0:ny-1,0:nz-1) :: mat !< 3D matrix of real values 
!
integer :: nx !< grid dimension in the x direction
integer :: ny !< grid dimension in the y direction
integer :: nz !< grid dimension in the z direction
!
real :: xmin !< Lower bound of grid extent in the x direction
real :: xmax !< Upper bound of grid extent in the x direction 
real :: ymin !< Lower bound of grid extent in the y direction 
real :: ymax !< Upper bound of grid extent in the y direction
real :: zmin !< Lower bound of grid extent in the z direction
real :: zmax !< Upper bound of grid extent in the z direction
!
real :: x !< x coordinate of the interpolation point 
real :: y !< y coordinate of the interpolation point 
real :: z !< z coordinate of the interpolation point
!
! Internal variables
!
real :: dx, delta, v1, v2
integer :: i
!
dx = (xmax-xmin)/real(nx-1)
!
i = int((x-xmin)/dx)
if(i==nx)i=i-1
!
delta = x-i*dx
!
v1 = bilinear_interpolation(mat(i  ,:,:),ny,nz,ymin,ymax,zmin,zmax,y,z)
v2 = bilinear_interpolation(mat(i+1,:,:),ny,nz,ymin,ymax,zmin,zmax,y,z)
!
trilinear_interpolation = (1.E0-delta)*v1+delta*v2 
!
!----------------------------------
end function trilinear_iterpolation
!----------------------------------
!
! @author Yder MASSON
! @date November 12, 2014
! @brief Preforms 2D bilinear interpolation 
!
!----------------------------------------------------------------------
real function bilinear_interpolation(mat,nx,ny,xmin,xmax,ymin,ymax,x,y)
!----------------------------------------------------------------------
!
implicit none
!
real, dimension(0:nx-1,0:ny-1) :: mat !< 2D matrix of real values 
!
integer :: nx !< Number of grid points in the x direction
integer :: ny !< Number of grid points in the y direction
!
real :: xmin !< Lower bound of grid extent in the x direction
real :: xmax !< Upper bound of grid extent in the x direction
real :: ymin !< Lower bound of grid extent in the y direction
real :: ymax !< Upper bound of grid extent in the y direction
!
real :: x !< x coordinate of the interpolation point 
real :: y !< y coordinate of the interpolation point 
!
! Internal variables
!
real :: dx, delta, v1, v2
integer :: i
!
dx = (xmax-xmin)/real(nx-1)
!
i = int((x-xmin)/dx)
if(i==nx)i=i-1
!
delta = x-i*dx
!
v1 = linear_interpolation(mat(i  ,:),ny,ymin,ymax,y)
v2 = linear_interpolation(mat(i+1,:),ny,ymin,ymax,y)
!
bilinear_interpolation = (1.E0-delta)*v1+delta*v2 
! 
!----------------------------------
end function bilinear_interpolation
!----------------------------------
!
! @author Yder MASSON
! @date November 12, 2014
! @brief Preforms 1D linear interpolation 
!
!-----------------------------------------------------
real function linear_interpolation(mat,nx,xmin,xmax,x)
!-----------------------------------------------------
!
real, dimension(0:nx-1) :: mat !< 1D vector of real values 
!
integer :: nx !< Number of grid points 
!
real :: xmin !< Lower bound of grid extent 
real :: xmax !< Upper bound of grid extent
!
real :: x !< x coordinate of the interpolation point 
!
! Internal variables
!
real :: dx, delta, v1, v2
integer :: i
!
dx = (xmax-xmin)/real(nx-1)
!
i = int((x-xmin)/dx)
if(i==nx)i=i-1
!
delta = x-i*dx
!
v1 = mat(i)
v2 = mat(i+1)
!
linear_interpolation = (1.E0-delta)*v1+delta*v2
! 
!--------------------------------
end function linear_interpolation
!--------------------------------
!
!==============================
end module module_interpolation
!==============================
