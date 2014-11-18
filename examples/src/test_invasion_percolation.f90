!> @author Yder MASSON
!> @date November 1, 2014
!> @brief This programs makes sure all the invasion functions 
!> are producing the same output for a given lattice
!================================
program test_invasion_percolation
!================================
!
use module_invasion_percolation
!
implicit none
!
integer :: nx
integer :: ny
integer :: nz
integer :: n_sites
!
real :: dx 
real :: dy 
real :: dz 
!
real, allocatable :: x(:), y(:), z(:)
!
logical :: period_x 
logical :: period_y
logical :: period_z 
!
real, allocatable :: values(:,:,:) 
real, allocatable :: values_cubic(:,:,:) 
real, allocatable :: values_cubic_simple(:,:,:) 
real, allocatable :: values_arbitrary(:) 
real, allocatable :: values_arbitrary_simple(:) 
!
integer, allocatable :: states(:,:,:) 
integer, allocatable :: states_cubic(:,:,:) 
integer, allocatable :: states_cubic_simple(:,:,:) 
integer, allocatable :: states_arbitrary(:) 
integer, allocatable :: states_arbitrary_simple(:) 
!
integer, allocatable :: offsets(:) 
integer, allocatable :: connectivity(:) 
!
integer :: n_sites_invaded_cubic
integer :: n_sites_invaded_cubic_simple
integer :: n_sites_invaded_arbitrary
integer :: n_sites_invaded_arbitrary_simple
!
integer, allocatable :: invasion_list_cubic(:)
integer, allocatable :: invasion_list_cubic_simple(:)
integer, allocatable :: invasion_list_arbitrary(:)
integer, allocatable :: invasion_list_arbitrary_simple(:)
!
logical :: gravity 
logical :: trapping 
!
real :: sigma
real :: theta_c 
real :: delta_rho 
real :: pi
!
real :: gx
real :: gy 
real :: gz 
!
integer :: i,j,k,l, n, offset, i_site, j_site, k_site
integer :: i_plus_one, i_minus_one
integer :: j_plus_one, j_minus_one
integer :: k_plus_one, k_minus_one
!
pi = 4.*atan(1.)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! setting up an Invasion Percolation simulation !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! setup grid dimensions :
!
nx = 20
ny = 13
nz = 45
!
n_sites = nx*ny*nz
!
dx = 1.0
dy = 1.0
dz = 1.0
!
allocate(x(n_sites))
allocate(y(n_sites))
allocate(z(n_sites))
!
l = 0
!
do k = 1,nz
   do j = 1,ny
      do i = 1,nx
         l = l+1
         x(l) = (i-1)*dx 
         y(l) = (j-1)*dy 
         z(l) = (k-1)*dz 
      enddo
   enddo
enddo
!
! setup domain boundaries : 
!
! use period_x = .true. for periodic bondaries in the x direction
! use period_x = .false. for sealed walls in the x direction
!
period_x = .false.
period_y = .false.
period_z = .false.
!
! allocate memory :
!
allocate(values(nx,ny,nz))
allocate(states(nx,ny,nz))
!
allocate(values_cubic(nx,ny,nz))
allocate(states_cubic(nx,ny,nz))
allocate(invasion_list_cubic(nx*ny*nz))
!
allocate(values_cubic_simple(nx,ny,nz))
allocate(states_cubic_simple(nx,ny,nz))
allocate(invasion_list_cubic_simple(nx*ny*nz))
!
allocate(values_arbitrary(nx*ny*nz))
allocate(states_arbitrary(nx*ny*nz))
allocate(invasion_list_arbitrary(nx*ny*nz))
!
allocate(values_arbitrary_simple(nx*ny*nz))
allocate(states_arbitrary_simple(nx*ny*nz))
allocate(invasion_list_arbitrary_simple(nx*ny*nz))
!
allocate(offsets(n_sites)) 
allocate(connectivity(n_sites+n_sites*6)) 
!
! setup arbitrary lattice structure :
!
offset = 0
!
do k_site = 1,nz
   do j_site = 1,ny
      do i_site = 1,nx
         !
         offset = offset+1
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
         ! Init number of neighboring sites
         !
         n = 0
         !
         if(i_minus_one>=1.and.nx/=1)then
            n=n+1
            offset=offset+1
            connectivity(offset) = ijk2ind(nx,ny,i_minus_one,j_site,k_site)
         endif
         if(i_plus_one<=nx.and.nx/=1)then
            n=n+1
            offset=offset+1
            connectivity(offset) = ijk2ind(nx,ny,i_plus_one,j_site,k_site)
         endif
         if(j_minus_one>=1.and.ny/=1)then
            n=n+1
            offset=offset+1
            connectivity(offset) = ijk2ind(nx,ny,i_site,j_minus_one,k_site)
         endif
         if(j_plus_one<=ny.and.ny/=1)then
            n=n+1
            offset=offset+1
            connectivity(offset) = ijk2ind(nx,ny,i_site,j_plus_one,k_site)
         endif
         if(k_minus_one>=1.and.nz/=1)then
            n=n+1
            offset=offset+1
            connectivity(offset) = ijk2ind(nx,ny,i_site,j_site,k_minus_one)
         endif
         if(k_plus_one<=nz.and.nz/=1)then
            n=n+1
            offset=offset+1
            connectivity(offset) = ijk2ind(nx,ny,i_site,j_site,k_plus_one)
         endif
         !
         ! store site's number of neighbors
         !
         connectivity(offset-n) = n
         !
         ! store offsets
         !
         offsets(ijk2ind(nx,ny,i_site,j_site,k_site)) = offset-n
         !
      enddo
   enddo
enddo
!
! setup sites's percolation potentials
!
do k = 1,nz
   do j = 1,ny
      do i = 1,nx
         values(i,j,k) = rand(0)*dx
      enddo
   enddo
enddo
!
! initialize sites's states :
!
states(:,:,:) = not_invaded
!
! setup injection region (bottom wall here) :
!
states(:,:,1) = neighboring
!
! setup exit region (top wall here) :
!
states(:,:,nz) = exit_site
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
gravity = .true.
!
! setup trapping
!
trapping = .false.
!
! same setup for all simulations
!
l = 0
!
do k = 1,nz
   do j = 1,ny
      do i = 1,nx
         !
         l = l+1
         !
         values_cubic(i,j,k) = values(i,j,k)
         values_cubic_simple(i,j,k) = values(i,j,k)
         values_arbitrary(l) = values(i,j,k)
         values_arbitrary_simple(l) = values(i,j,k)
         !
         states_cubic(i,j,k) = states(i,j,k)
         states_cubic_simple(i,j,k) = states(i,j,k)
         states_arbitrary(l) = states(i,j,k)
         states_arbitrary_simple(l) = states(i,j,k)
         !
      enddo
   enddo
enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Now run the four simulations !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
call invade_cubic_lattice_simple(nx,                          &
                                 ny,                          &
                                 nz,                          &
                                 dx,                          &
                                 dy,                          &
                                 dz,                          &
                                 period_x,                    &
                                 period_y,                    &
                                 period_z,                    &
                                 values_cubic_simple,         &
                                 states_cubic_simple,         &
                                 n_sites_invaded_cubic_simple,&
                                 invasion_list_cubic_simple,  &
                                 gravity,                     &
                                 trapping,                    &
                                 sigma,                       &
                                 theta_c,                     &
                                 delta_rho,                   &
                                 gx,                          &
                                 gy,                          &
                                 gz                           )
!
write(*,*)'simulation 1 finished'
! 
call invade_cubic_lattice_fast(nx,                    &
                               ny,                    &
                               nz,                    &
                               dx,                    &
                               dy,                    &
                               dz,                    &
                               period_x,              &
                               period_y,              &
                               period_z,              &
                               values_cubic,          &
                               states_cubic,          &
                               n_sites_invaded_cubic, &
                               invasion_list_cubic,   &
                               gravity,               &
                               trapping,              &
                               sigma,                 &
                               theta_c,               &
                               delta_rho,             &
                               gx,                    &
                               gy,                    &
                               gz                     )
!
write(*,*)'simulation 2 finished'
!
call invade_arbitrary_lattice_simple(n_sites,                          &
                                     x,                                &
                                     y,                                &
                                     z,                                &
                                     offsets,                          &
                                     connectivity,                     &
                                     values_arbitrary_simple,          &
                                     states_arbitrary_simple,          &
                                     n_sites_invaded_arbitrary_simple, &
                                     invasion_list_arbitrary_simple,   &
                                     gravity,                          &
                                     trapping,                         &
                                     sigma,                            &
                                     theta_c,                          &
                                     delta_rho,                        &
                                     gx,                               &
                                     gy,                               &
                                     gz                                )
!
write(*,*)'simulation 3 finished'
!
call invade_arbitrary_lattice_fast(n_sites,                   &
                                   x,                         &
                                   y,                         &
                                   z,                         &
                                   offsets,                   &
                                   connectivity,              &
                                   values_arbitrary,          &
                                   states_arbitrary,          &
                                   n_sites_invaded_arbitrary, &
                                   invasion_list_arbitrary,   &
                                   gravity,                   &
                                   trapping,                  &
                                   sigma,                     &
                                   theta_c,                   &
                                   delta_rho,                 &
                                   gx,                        &
                                   gy,                        &
                                   gz                         )
!
write(*,*)'simulation 4 finished'
!
! make sure all simulation gives the same result
!
write(*,*)'number of sites invaded cubic simple: ',n_sites_invaded_cubic_simple
write(*,*)'number of sites invaded cubic: ',n_sites_invaded_cubic
write(*,*)'number of sites invaded arbitrary simple: ',n_sites_invaded_arbitrary_simple
write(*,*)'number of sites invaded arbitrary: ',n_sites_invaded_arbitrary
stop
!
!====================================
end program test_invasion_percolation
!====================================
