!========================!
!                        !
!                        !
!     MODULE GRAVITY     !
!                        !
!                        !
!========================!
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
!> @date October, 7 2014
!> @brief This module contains the functions to account for gravity
!> @details The functions in this module compute the invasion potential 
!> \f$ P_i = \frac{2 \sigma \mbox{cos} \theta_c}{a_i}-\Delta \rho g (L-z_i)\f$
!> @warning Be careful when using periodic boundaries, gravity gets periodic as well...
!
!====================
module module_gravity
!====================
  !
contains
  !
  !> @author Yder MASSON
  !> @date October, 7 2014
  !> @brief Setup sites's invasion potential for cubic lattices
  !> @details The invasion potential is defined as: 
  !> \f$ P_i = \frac{2 \sigma \mbox{cos} \theta_c}{a_i}-\Delta \rho g (L-z_i)\f$
  !> @warning Be careful when using periodic boundaries, gravity gets periodic as well...
  !
  !------------------------------------------------
  subroutine add_gravity_cubic_lattice(values,    &
                                       nx,        &
                                       ny,        &
                                       nz,        &
                                       dx,        &
                                       dy,        &
                                       dz,        &
                                       sigma,     &
                                       theta_c,   &
                                       delta_rho, &
                                       gx,        &
                                       gy,        &
                                       gz         )
  !------------------------------------------------
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
    !> Output array containig the sites's invasion potentials (\b input/\b output)
    !> \n At input time, this array contains the pores's sizes \f$ a_i \f$
    !> \n At output time, this array contains the invasion potential 
    !> \f$ P_i = \frac{2 \sigma \mbox{cos} \theta_c}{a_i}-\Delta \rho g (L-z_i)\f$
    real :: values(nx,ny,nz) 
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
    real :: l_max !< Height of the system to be invaded \f$ L \f$
    !
    integer :: i !< Looping index
    integer :: j !< Looping index
    integer :: k !< Looping index
    integer :: n !< counter
    !
    real :: norm_g !< norm of the gravity vector with components gx,gy,gz.
    real :: l !< height of the site taken in the direction of unit gravity vector
    !
    real :: x !< site's coordinate x
    real :: y !< site's coordinate y
    real :: z !< site's coordinate z
    !
    ! Computing norm_g is useless, it is only used for clarity
    !
    norm_g = sqrt(gx**2+gy**2+gz**2)
    !
    ! First find l_max
    !
    n = 0
    !
    do k = 1,nz
       z = k*dz
       do j = 1,ny
          y = j*dy
          do i = 1,nx
             x = i*dx
             n = n+1
             l = (x*gx+y*gy+z*gz)/norm_g
             if(n==1)then
                l_max = l
             else
                l_max = max(l_max,l)
             endif
          enddo
       enddo
    enddo
    !
    ! Compute site's invasion potential
    !
    do k = 1,nz
       z = k*dz
       do j = 1,ny
          y = j*dy
          do i = 1,nx
             x = i*dx
             l = (x*gx+y*gy+z*gz)/norm_g
             values(i,j,k) = 2*sigma*cos(theta_c)/values(i,j,k)&
                            -delta_rho*norm_g*(l_max-l)
          enddo
       enddo
    enddo
    ! 
    return
    !
  !---------------------------------------
  end subroutine add_gravity_cubic_lattice
  !---------------------------------------
  !
  !> @author Yder MASSON
  !> @date October, 7 2014
  !> @brief Setup sites's invasion potential for arbitrary lattices
  !> @details The invasion potential is defined as:
  !> \f$ P_i = \frac{2 \sigma \mbox{cos} \theta_c}{a_i}-\Delta \rho g (L-z_i)\f$
  !> @warning Be careful when using periodic boundaries, gravity gets periodic as well...
  !
  !----------------------------------------------------
  subroutine add_gravity_arbitrary_lattice(values,    &
                                           n,         &
                                           x,         &
                                           y,         &
                                           z,         &
                                           sigma,     &
                                           theta_c,   &
                                           delta_rho, &
                                           gx,        &
                                           gy,        &
                                           gz         )
  !----------------------------------------------------
    !
    implicit none
    !
    integer :: n !< Total number of sites in the lattice (\b input) 
    !
    real :: x(:) !< Array containing the x coordinates od the sites (\b input)
    real :: y(:) !< Array containing the y coordinates od the sites (\b input)
    real :: z(:) !< Array containing the z coordinates od the sites (\b input)
    !
    !
    !> Output array containig the sites's invasion potentials (\b input/\b output)
    !> \n At input time, this array contains the pores's sizes \f$ a_i \f$
    !> \n At output time, this array contains the invasion potential 
    !> \f$ P_i = \frac{2 \sigma \mbox{cos} \theta_c}{a_i}-\Delta \rho g (L-z_i)\f$
    real :: values(:) 
    !
    real :: sigma !< Surface tension \f$ \sigma \f$ (\b input)
    real :: theta_c !< Equilibrium contact angle \f$ \theta_c \f$ (\b input)
    real :: delta_rho !< Fluid density contrast \f$ \Delta \rho \f$ (\b input)
    !
    real :: gx !< Acceleration of gravity \f$ g \f$ in the x direction (\b input)
    real :: gy !< Acceleration of gravity \f$ g \f$ in the y direction (\b input)
    real :: gz !< Acceleration of gravity \f$ g \f$ in the z direction (\b input)
    !
    real :: l_max !< Height of the system to be invaded \f$ L \f$
    integer :: i !< Looping index
    real :: norm_g !< norm of the gravity vector with components gx,gy,gz.
    real :: l !< height of the site taken in the direction of unit gravity vector
    !
    ! Computing norm_g is useless, it is only used for clarity
    !
    norm_g = sqrt(gx**2+gy**2+gz**2)
    !
    ! First find l_max
    !
    do i = 1,n
       l = (x(i)*gx+y(i)*gy+z(i)*gz)/norm_g
       if(i==1)then
          l_max = l
       else
          l_max = max(l_max,l)
       endif
    enddo
    !
    ! Compute site's invasion potential
    !
    do i = 1,n
       l = (x(i)*gx+y(i)*gy+z(i)*gz)/norm_g
       values(i) = 2*sigma*cos(theta_c)/values(i)&
                  -delta_rho*norm_g*(l_max-l)
    enddo
    ! 
    return
    !
  !-------------------------------------------
  end subroutine add_gravity_arbitrary_lattice
  !-------------------------------------------
  !
!========================
end module module_gravity
!========================
