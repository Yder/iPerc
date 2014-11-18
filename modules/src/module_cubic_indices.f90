!==============================!
!                              !
!                              !
!     MODULE CUBIC INDICES     ! 
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
!> @author Yder MASSON
!> @date October 23, 2014
!> @brief This module contains functions to transform 3D indices to 1D index and vice versa
!
!==========================
module module_cubic_indices
!==========================
  !
contains
  !
  !> Transform 3D indices (i,j,k) to 1D index (ijk2ind)
  !------------------------------------
  integer function ijk2ind(nx,ny,i,j,k)
  !------------------------------------
    !
    implicit none
    !
    integer :: nx !< Grid dimension in the x direction (\b input)
    integer :: ny !< Grid dimension in the y direction (\b input)
    integer :: i  !< i index, \f$ 1<i<\mbox{nx} \f$ (\b input)
    integer :: j  !< j index, \f$ 1<j<\mbox{ny} \f$ (\b input)
    integer :: k  !< k index, \f$ 1<k<\mbox{nz} \f$ (\b input)
    !
    ijk2ind = (k-1)*nx*ny+(j-1)*nx+i
    !
  !-------------------
  end function ijk2ind
  !-------------------
  !
  !> Transform  1D index (ijk2ind) to 3D indices (i,j,k)
  !----------------------------------
  subroutine ind2ijk(nx,ny,ind,i,j,k)
  !----------------------------------
    !
    implicit none
    !
    integer :: ind !< 1D index, \f$ 1<\mbox{ind}<\mbox{nx}\times\mbox{ny}\times\mbox{nz}\f$
    integer :: nx !< Grid dimension in the x direction (\b input)
    integer :: ny !< Grid dimension in the y direction (\b input)
    integer :: i  !< i index, \f$ 1<\mbox{i}<\mbox{nx} \f$ (\b input)
    integer :: j  !< j index, \f$ 1<\mbox{j}<\mbox{ny} \f$ (\b input)
    integer :: k  !< k index, \f$ 1<\mbox{k}<\mbox{nz} \f$ (\b input)
    !
    k = (ind-1)/nx/ny+1
    j = (ind-(k-1)*nx*ny-1)/nx+1
    i = ind-(k-1)*nx*ny-(j-1)*nx
    !
  !---------------------
  end subroutine ind2ijk
  !---------------------
  !
!==============================
end module module_cubic_indices
!==============================

