!=============================!
!                             !
!                             !
!     MODULE RANDOM MEDIA     !
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
!> @author Yder MASSON
!> @date November 9, 2014
!> @brief This module contains function for generating random media
!> having a determined correlation function
!>
!> @see KLIMEŠ, Lulěk. Correlation functions of random media.
!> Pure and applied geophysics, 2002, vol. 159, no 7-8, p. 1811-1831.
!
!=========================
module module_random_media
!=========================
  !
contains
  !
  !> @author Yder MASSON
  !> @date November 8, 2014
  !> @brief Generates a 3D matrix of random real numbers 
  !> with zero mean and unity standard deviation 
  !> and having a given correlation function
  !> @details The correlation function currently available are : \n
  !> \n \n General     : 
  !> \f$ \hat{f}(k) = \kappa \left[ a^{-2}+k^2\right]^{-\frac{d}{4}-\frac{H}{2}} 
  !> \mbox{exp}\left(\frac{a_G^2k^2}{8}\right)\f$
  !> \n \n Gaussian    : 
  !> \f$ \hat{f}(k) = \kappa \mbox{exp}\left(\frac{a_G^2k^2}{8}\right)\f$
  !> \n \n Von Karman  : 
  !> \f$ \hat{f}(k) = \kappa \left[ a^{-2}+k^2\right]^{-\frac{d}{4}-\frac{H}{2}}\f$
  !> \n \n Exponential : 
  !> \f$ \hat{f}(k) = \kappa \left[ a^{-2}+k^2\right]^{-\frac{d+1}{4}}\f$
  !> \n \n Self-affine : 
  !> \f$ \hat{f}(k) = \kappa k^{-\frac{d}{2}-H}\f$
  !> \n \n Kummer      : 
  !>\f$ \hat{f}(k) = \kappa k^{-\frac{d}{2}-H} \mbox{exp}\left(\frac{a_G^2k^2}{8}\right)\f$
  !> \n \n White noise : 
  !> \f$ \hat{f}(k) = \kappa \f$
  !> \n \n in the above expression \f$ d = 3\f$ in 3D
  !>
  !> @see KLIMEŠ, Lulěk. Correlation functions of random media. 
  !> Pure and applied geophysics, 2002, vol. 159, no 7-8, p. 1811-1831.
  !
  !----------------------------------------------------
  subroutine gen_random_media_3D (mat,                & 
                                  nx, ny, nz,         &
                                  dx, dy, dz,         &
                                  H,                  &
                                  ak,                 &
                                  ag,                 &
                                  correlation_function)
  !----------------------------------------------------
    !
    implicit none
    !
    real, dimension(nx,ny,nz) :: mat !< 3D matrix containing the random values (\b Output)
    !
    integer :: nx !< Grid dimension in the x direction (\b Input)
    integer :: ny !< Grid dimension in the y directionn (\b Input)
    integer :: nz !< Grid dimension in the z directionn (\b Input)
    !
    real :: dx !< Grid spacing in the x directionn (\b Input)
    real :: dy !< Grid spacing in the y directionn (\b Input)
    real :: dz !< Grid spacing in the z directionn (\b Input)
    !
    real :: H !< Hurst exponent \f$ H \f$ n (\b Input)
    real :: ak !< Von karman correlation length \f$ a \f$ n (\b Input)
    real :: ag !< Gaussian correlation length \f$ a_G \f$ n (\b Input)
    !
    !> Desired correlation function (character string \b Input) : \n
    !> \n use <b>correlation_function='general'</b> 
    !> for the most general correlation function
    !> \n use <b>correlation_function='gaussian'</b> 
    !> for a gaussian correlation function
    !> \n use <b>correlation_function='von_karman'</b> 
    !> for a Von Karman correlation function
    !> \n use <b>correlation_function='self_affine'</b> 
    !> for a self-affine correlation function
    !> \n use <b>correlation_function='kummer'</b> 
    !> for a Kummer correlation function
    !> \n use <b>correlation_function='white_noise'</b> 
    !>for a white noise (you may not need this function for that)
    !
    character(len=*) :: correlation_function 
    !
    ! internal variables
    !
    integer, parameter :: ndim = 3
    integer :: nn(ndim),i ,j, k, n
    real :: kx, ky, kz, dkx, dky, dkz, k2, pi, filter
    real, dimension(:), allocatable :: fft_mat
    !      
    ! define PI
    !
    PI = 4.0*atan(1.0)
    !
    ! get fft_mat dimension (nn(i) must be a power of 2)
    !
    nn(1) = 2**ceiling(log10(dble(nx))/log10(dble(2)))
    nn(2) = 2**ceiling(log10(dble(ny))/log10(dble(2)))
    nn(3) = 2**ceiling(log10(dble(nz))/log10(dble(2)))
    !
    ! allocate memory for fft_mat
    !
    allocate(fft_mat(2*nn(1)*nn(2)*nn(3)))
    !
    ! determine wavenumber spacings
    !
    dkx = 2 * pi / real(nn(1)) / real(dx)
    dky = 2 * pi / real(nn(2)) / real(dy)
    dkz = 2 * pi / real(nn(3)) / real(dz)
    !
    ! fill work array with random numbers having zero mean
    !
    do n = 1,nn(1)*nn(2)*nn(3)*2,2
       fft_mat(n) = rand()-0.5
       fft_mat(n+1) = 0.E0
    enddo
    !
    ! compute fft
    !
    call fourn (fft_mat, nn, ndim, 1)
    !
    ! force zero mean
    !
    fft_mat(1) = 0.E0 ! real part
    fft_mat(2) = 0.E0 ! imaginary part
    !
    n = 1
    !
    do k = 1, nn(3)
       !
       kz = k-1
       if(k>nn(3)/2+1)kz = kz-nn(3)
       kz = kz*dkz
       !
       do j = 1, nn(2)
          !
          ky = j-1
          if(j>nn(2)/2+1)ky = ky-nn(2)
          ky = ky*dky
          !
          do i = 1, nn(1)
             !
             kx = i-1
             if(i>nn(1)/2+1)kx = kx-nn(1)
             kx = kx*dkx
             !
             k2 = kx**2 + ky**2 + kz**2
             !
             ! get filter value
             !
             select case(correlation_function)
             case('general') 
                filter = (ak**(-2)+k2)**(-real(ndim)/4.-H/2.)*exp(-(ag**2*k2/8.))
             case ('gaussian') 
                filter = exp(-(ag**2*k2/8.))
             case ('von_karman') 
                filter = (ak**(-2)+k2)**(-real(ndim)/4.-H/2.)
             case ('exponential') 
                filter = (ak**(-2)+k2)**(-real(ndim+1)/4.)
             case ('self_affine') 
                filter = k2**(-real(ndim)/2.-H)
             case ('kummer') 
                filter = k2**(-real(ndim)/2.-H)*exp(-(ag**2*k2/8.))
             case ('white_noise') 
                filter = 1.
             case default
                print*,'Error: Wrong correlation_function in gen_random_media_3D.'
             end select
             !
             if(k2==0)filter = 0
             !
             fft_mat(n)   = filter * fft_mat(n)   ! real part
             fft_mat(n+1) = filter * fft_mat(n+1) ! imaginary part
             !
             n = n+2
             !
          enddo
       enddo
    enddo
    !
    ! compute inverse fft
    !
    call fourn (fft_mat, nn, ndim, -1)
    !
    n = 1
    !
    do k = 1,nn(3)
       do j = 1,nn(2)
          do i = 1,nn(1)
             if(i<=nx.and.j<=ny.and.k<=nz)mat(i,j,k) = fft_mat(n)
             n = n+2
          enddo
       enddo
    enddo
    !
    ! free memory
    !
    deallocate(fft_mat)
    !
    ! substract mean
    !
    mat = mat - sum(mat)/size(mat)
    !
    ! normalize standard deviation
    !
    mat = mat/sqrt(sum(mat**2)/size(mat))
    !
    !-----
    return
    !-----
    !
  !---------------------------------
  end subroutine gen_random_media_3D
  !---------------------------------
  !
  !> @author Yder MASSON
  !> @date November 8, 2014
  !> @brief Generates a 3D matrix of random real numbers 
  !> with zero mean and unity standard deviation 
  !> and having a given correlation function
  !> @details The correlation function currently available are : \n
  !> \n \n General     : 
  !> \f$ \hat{f}(k) = \kappa \left[ a^{-2}+k^2\right]^{-\frac{d}{4}-\frac{H}{2}} 
  !> \mbox{exp}\left(\frac{a_G^2k^2}{8}\right)\f$
  !> \n \n Gaussian    : 
  !> \f$ \hat{f}(k) = \kappa \mbox{exp}\left(\frac{a_G^2k^2}{8}\right)\f$
  !> \n \n Von Karman  : 
  !> \f$ \hat{f}(k) = \kappa \left[ a^{-2}+k^2\right]^{-\frac{d}{4}-\frac{H}{2}}\f$
  !> \n \n Exponential : 
  !> \f$ \hat{f}(k) = \kappa \left[ a^{-2}+k^2\right]^{-\frac{d+1}{4}}\f$
  !> \n \n Self-affine : 
  !> \f$ \hat{f}(k) = \kappa k^{-\frac{d}{2}-H}\f$
  !> \n \n Kummer      : 
  !>\f$ \hat{f}(k) = \kappa k^{-\frac{d}{2}-H} \mbox{exp}\left(\frac{a_G^2k^2}{8}\right)\f$
  !> \n \n White noise : 
  !> \f$ \hat{f}(k) = \kappa \f$
  !> \n \n in the above expression \f$ d = 2\f$ in 2D
  !>
  !> @see KLIMEŠ, Lulěk. Correlation functions of random media. 
  !> Pure and applied geophysics, 2002, vol. 159, no 7-8, p. 1811-1831.
  !
  !----------------------------------------------------
  subroutine gen_random_media_2D (mat,                & 
                                  nx, ny,             &
                                  dx, dy,             &
                                  H,                  &
                                  ak,                 &
                                  ag,                 &
                                  correlation_function)
  !----------------------------------------------------
    !
    implicit none
    !
    real, dimension(nx,ny) :: mat !< 2D matrix containing the random values (\b Output)
    !
    integer :: nx !< Grid dimension in the x direction (\b Input)
    integer :: ny !< Grid dimension in the y directionn (\b Input)
    !
    real :: dx !< Grid spacing in the x directionn (\b Input)
    real :: dy !< Grid spacing in the y directionn (\b Input)
    !
    real :: H !< Hurst exponent \f$ H \f$ n (\b Input)
    real :: ak !< Von karman correlation length \f$ a \f$ n (\b Input)
    real :: ag !< Gaussian correlation length \f$ a_G \f$ n (\b Input)
    !
    !> Desired correlation function (character string \b Input) : \n
    !> \n use <b>correlation_function='general'</b> 
    !> for the most general correlation function
    !> \n use <b>correlation_function='gaussian'</b> 
    !> for a gaussian correlation function
    !> \n use <b>correlation_function='von_karman'</b> 
    !> for a Von Karman correlation function
    !> \n use <b>correlation_function='self_affine'</b> 
    !> for a self-affine correlation function
    !> \n use <b>correlation_function='kummer'</b> 
    !> for a Kummer correlation function
    !> \n use <b>correlation_function='white_noise'</b> 
    !>for a white noise (you may not need this function for that)
    !
    character(len=*) :: correlation_function 
    !
    ! internal variables
    !
    integer, parameter :: ndim = 2
    integer :: nn(ndim),i ,j, n
    real :: kx, ky, dkx, dky, k2, pi, filter
    real, dimension(:), allocatable :: fft_mat
    !      
    ! define PI
    !
    PI = 4.0*atan(1.0)
    !
    ! get fft_mat dimension (nn(i) must be a power of 2)
    !
    nn(1) = 2**ceiling(log10(dble(nx))/log10(dble(2)))
    nn(2) = 2**ceiling(log10(dble(ny))/log10(dble(2)))
    !
    ! allocate memory for fft_mat
    !
    allocate(fft_mat(2*nn(1)*nn(2)))
    !
    ! determine wavenumber spacings
    !
    dkx = 2 * pi / real(nn(1)) / real(dx)
    dky = 2 * pi / real(nn(2)) / real(dy)
    !
    ! fill work array with random numbers having zero mean
    !
    do n = 1,nn(1)*nn(2)*2,2
       fft_mat(n) = rand()-0.5
       fft_mat(n+1) = 0.E0
    enddo
    !
    ! compute fft
    !
    call fourn (fft_mat, nn, ndim, 1)
    !
    ! force zero mean
    !
    fft_mat(1) = 0.E0 ! real part
    fft_mat(2) = 0.E0 ! imaginary part
    !
    n = 1
    !
    do j = 1, nn(2)
       !
       ky = j-1
       if(j>nn(2)/2+1)ky = ky-nn(2)
       ky = ky*dky
       !
       do i = 1, nn(1)
          !
          kx = i-1
          if(i>nn(1)/2+1)kx = kx-nn(1)
          kx = kx*dkx
          !
          k2 = kx**2 + ky**2
          !
          ! get filter value
          !
          select case(correlation_function)
          case('general') 
             filter = (ak**(-2)+k2)**(-real(ndim)/4.-H/2.)*exp(-(ag**2*k2/8.))
          case ('gaussian') 
             filter = exp(-(ag**2*k2/8.))
          case ('von_karman') 
             filter = (ak**(-2)+k2)**(-real(ndim)/4.-H/2.)
          case ('exponential') 
             filter = (ak**(-2)+k2)**(-real(ndim+1)/4.)
          case ('self_affine') 
             filter = k2**(-real(ndim)/2.-H)
          case ('kummer') 
             filter = k2**(-real(ndim)/2.-H)*exp(-(ag**2*k2/8.))
          case ('white_noise') 
             filter = 1.
          case default
             print*,'Error: Wrong correlation_function in gen_random_media_3D.'
          end select
          !
          if(k2==0)filter = 0
          !
          fft_mat(n)   = filter * fft_mat(n)   ! real part
          fft_mat(n+1) = filter * fft_mat(n+1) ! imaginary part
          !
          n = n+2
          !
       enddo
    enddo
    !
    ! compute inverse fft
    !
    call fourn (fft_mat, nn, ndim, -1)
    !
    n = 1
    !
    do j = 1,nn(2)
       do i = 1,nn(1)
          if(i<=nx.and.j<=ny)mat(i,j) = fft_mat(n)
          n = n+2
       enddo
    enddo
    !
    ! free memory
    !
    deallocate(fft_mat)
    !
    ! substract mean
    !
    mat = mat - sum(mat)/size(mat)
    !
    ! normalize standard deviation
    !
    mat = mat/sqrt(sum(mat**2)/size(mat))
    !
    !-----
    return
    !-----
    !
  !---------------------------------
  end subroutine gen_random_media_2D
  !---------------------------------
  !
  ! @brief Fast Fourier Transform
  !-----------------------------------
  subroutine fourn(data,nn,ndim,isign)
  !-----------------------------------
    !
    implicit none
    !
    integer :: isign,ndim,nn(ndim)
    real :: data(*)
    integer :: i1,i2,i2rev,i3,i3rev,ibit,idim
    integer :: ifp1,ifp2,ip1,ip2,ip3,k1,k2,n,nprev,nrem,ntot
    real :: tempi,tempr
    double precision :: theta,wi,wpi,wpr,wr,wtemp
    !
    ntot=1
    do idim=1,ndim
       ntot=ntot*nn(idim)
    enddo
    nprev=1
    do idim=1,ndim
       n=nn(idim)
       nrem=ntot/(n*nprev)
       ip1=2*nprev
       ip2=ip1*n
       ip3=ip2*nrem
       i2rev=1
       do i2=1,ip2,ip1
          if(i2.lt.i2rev)then
             do i1=i2,i2+ip1-2,2
                do i3=i1,ip3,ip2
                   i3rev=i2rev+i3-i2
                   tempr=data(i3)
                   tempi=data(i3+1)
                   data(i3)=data(i3rev)
                   data(i3+1)=data(i3rev+1)
                   data(i3rev)=tempr
                   data(i3rev+1)=tempi
                enddo
             enddo
          endif
          ibit=ip2/2
          do while((ibit.ge.ip1).and.(i2rev.gt.ibit)) 
             i2rev=i2rev-ibit
             ibit=ibit/2
          enddo
          i2rev=i2rev+ibit
       enddo
       ifp1=ip1
       do while(ifp1.lt.ip2)
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*sin(0.5d0*theta)**2
          wpi=sin(theta)
          wr=1.d0
          wi=0.d0
          do i3=1,ifp1,ip1
             do i1=i3,i3+ip1-2,2
                do i2=i1,ip3,ifp2
                   k1=i2
                   k2=k1+ifp1
                   tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                   tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                   data(k2)=data(k1)-tempr
                   data(k2+1)=data(k1+1)-tempi
                   data(k1)=data(k1)+tempr
                   data(k1+1)=data(k1+1)+tempi
                enddo
             enddo
             wtemp=wr
             wr=wr*wpr-wi*wpi+wr
             wi=wi*wpr+wtemp*wpi+wi
          enddo
          ifp1=ifp2
       enddo
       nprev=n*nprev
    enddo
    !
    return
    !
  !-------------------
  end subroutine fourn
  !-------------------
  !
!=============================
end module module_random_media
!=============================
