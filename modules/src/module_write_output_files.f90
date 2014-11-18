!===================================!
!                                   !
!                                   !
!     MODULE WRITE OUTPUT FILES     !
!                                   !
!                                   !
!===================================!
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
!> @brief This modulecontains functions 
!> for writing the simulations results
!> into data file, e,g., 
!> for post-processing and visualization.
!
!===============================
module module_write_output_files
!===============================
!
use module_trapping
use module_cubic_indices
use module_invasion_percolation_constants
!
contains
!
!> @author Yder MASSON
!> @date October 28, 2014
!> @brief write arbitrary lattice info to VTK file (.vtp) for viewing with e.g. Paraview
!> @details The states and values data are attached to points (i.e. sites),
!> you can view these using a glyph filter, for example. 
!> \n The bonds linking sites can be visualized by plotting the wireframe.
!
!---------------------------------------------------
subroutine save_arbitrary_lattice(n_sites,         &
                                  states,          &
                                  values,          &
                                  n_sites_invaded, &
                                  invasion_list,   &
                                  offsets,         &
                                  connectivity,    &
                                  x,y,z,           &
                                  file_name,       &
                                  unit_vtk         )
!---------------------------------------------------
!
implicit none
!
integer, intent(in) :: n_sites !< Total number of sites in the lattice 
!> Sites's states as defiend in invasion percolation module 
integer, dimension(:), intent(in) :: states
!> Sites's values as defiend in invasion percolation module 
real,    dimension(:), intent(in) :: values 
integer, intent(in) :: n_sites_invaded !< number of sites invaded
integer, dimension(:), intent(in) :: invasion_list !< list of invaded sites
!> Array containing the offsets of the data stored 
!> in the connectivity array 
integer, intent(in) :: offsets(:)
!> Array containing thelattice connectivity 
!> \n For a given site \b i, with offset <b>j = offsets(i)</b> :
!> \n <b>n=connectivity(j)</b> is the number of sites neighboring site \b i.
!> \n <b>connectivity(j+1, j+2, ... ,j+n)</b> 
!> contains the indices of the sites that are neighboring site \b i.
integer, intent(in) :: connectivity(:)
!
real,dimension(:), intent(in) :: x !< Array containing the x coordinates of the sites 
real,dimension(:), intent(in) :: y !< Array containing the y coordinates of the sites
real,dimension(:), intent(in) :: z !< Array containing the z coordinates of the sites 
!
character(len=*), intent(in) :: file_name !< Output file name, must have the.vtu extension 
integer, intent(in) :: unit_vtk !< logical unit for output file 
!
! internal variables
!
integer :: e_io
integer :: i, j, k, l, n, n_cells
integer :: dim_arrays, offset
character(1), parameter:: end_rec = char(10)
character*20 :: cnc, cnp, coffset
integer, dimension(:,:), allocatable :: connectivity_vtk
integer, dimension(n_sites) :: invasion_times
integer, dimension(n_sites) :: trapping_times
character(len=4) :: extension
!
! get trapping times
!

call get_trapping_times_arbitrary(n_sites,           &
                                  states,            &
                                  offsets,           &
                                  connectivity,      &
                                  n_sites_invaded,   &
                                  invasion_list,     &
                                  trapping_times     )
!
! setup invasion times array
!
invasion_times(:) = -1 ! default for undinvaded sites
do i = 1,n_sites_invaded
   invasion_times(invasion_list(i)) = i
enddo
!
! get files extension
!
l = len_trim(file_name)
extension(1:4) = file_name(l-3:l)
!
! save data according to file extension
!
select case (extension)

case('.vtp')
!
! find number of cells
!
n_cells = 0
!
do i = 1,n_sites
   j = offsets(i)
   n_cells = n_cells+connectivity(j)
enddo
!
n_cells = n_cells/2
!
! store connectivity info
!
! ... allocate memory
!
allocate(connectivity_vtk(n_cells,2))
!
n_cells = 0
!
do i = 1,n_sites
   k = offsets(i)
   n = connectivity(k)
   do j = k+1,k+n
    if(i<connectivity(j))then
      n_cells = n_cells+1
      connectivity_vtk(n_cells,1) = i
      connectivity_vtk(n_cells,2) = connectivity(j)
    endif
   enddo
enddo
!
! convert integer values to strings
!
write(cnc,*)n_cells
write(cnp,*)n_sites
!
open(unit       = unit_vtk,            &
     file       = trim(file_name),     &
     form       = 'unformatted',       &
     access     = 'stream',            &
     action     = 'write',             &
     status     = 'replace',           &
     convert    = 'little_endian',     &
     iostat     = e_io                 )
!
! init offset for appended data
!
offset = 0
!
! write header
!
write(unit=Unit_VTK,iostat=E_IO) &
     '<?xml version = "1.0"?>'//end_rec
!
write(unit=Unit_VTK,iostat=E_IO)   &
     '<VTKFile'                    &
     //' type = "PolyData"'        &
     //' version="0.1"'            &
     //' byte_order="LittleEndian"'&
     //'>'//end_rec
!
! start unstructured grid infos
!
write(unit=unit_vtk,iostat=e_io) &
     '<PolyData>'//end_rec
!
! start piece infos
!
write(unit=unit_vtk,iostat=e_io)   &
     '<Piece'                      &
     //' NumberOfPoints='          &
     //'"'//trim(adjustl(cnp))//'"'&
     //' NumberOfVerts='           &
     //'"'//trim(adjustl(cnp))//'"'&
     //' NumberOfLines='           &
     //'"'//trim(adjustl(cnc))//'"'&
     //'>'//end_rec
!
! start point data infos 
!
write(unit=unit_vtk,iostat=e_io) &
     '<PointData>'//end_rec
!
! ... states data infos
!
write(coffset,'(I10)')offset
!
write(unit=unit_vtk,iostat=e_io)       &
     '<DataArray'                      &
     //' type="Int32" '                &
     //' NumberOfComponents="1"'       &
     //' Name="states"'                &
     //' format="appended"'            &
     //' offset='                      &
     //'"'//trim(adjustl(coffset))//'"'&
     //'/>'//end_rec
!
offset = offset                        &
        +int(sizeof(states(1:n_sites)) &
        +sizeof(i))
!
! ... values data infos
!
write(coffset,'(I10)')offset
!
write(unit=unit_vtk,iostat=e_io)       &
     '<DataArray'                      &
     //' type="Float32" '              &
     //' NumberOfComponents="1"'       &
     //' Name="values"'                &
     //' format="appended"'            &
     //' offset='                      &
     //'"'//trim(adjustl(coffset))//'"'&
     //'/>'//end_rec
!
offset = offset                        &
        +int(sizeof(values(1:n_sites)) &
        +sizeof(i))           
!
! ... invasion time data infos
!
write(coffset,'(I10)')offset
!
write(unit=unit_vtk,iostat=e_io)       &
     '<DataArray'                      &
     //' type="Int32" '                &
     //' NumberOfComponents="1"'       &
     //' Name="invasion times"'        &
     //' format="appended"'            &
     //' offset='                      &
     //'"'//trim(adjustl(coffset))//'"'&
     //'/>'//end_rec
!
offset = offset                                &
        +int(sizeof(invasion_times(1:n_sites)) &
        +sizeof(i))

!
! ... trapping time data infos
!
write(coffset,'(I10)')offset
!
write(unit=unit_vtk,iostat=e_io)       &
     '<DataArray'                      &
     //' type="Int32" '                &
     //' NumberOfComponents="1"'       &
     //' Name="traping times"'         &
     //' format="appended"'            &
     //' offset='                      &
     //'"'//trim(adjustl(coffset))//'"'&
     //'/>'//end_rec
!
offset = offset                                &
        +int(sizeof(trapping_times(1:n_sites)) &
        +sizeof(i))
!
! end point data infos
!
write(unit=unit_vtk,iostat=e_io)  &
     '</PointData>'//end_rec
!
! start points infos
!
write(unit=unit_vtk,iostat=e_io)  &
     '<Points>'//end_rec
!
! ... points coordinates infos
!
write(coffset,'(I10)')offset
!
write(unit=unit_vtk,iostat=e_io)       &
     '<DataArray'                      &
     //' type="Float32" '              &
     //' NumberOfComponents="3"'       &
     //' Name="Points"'                &
     //' format="appended"'            &
     //' offset='                      &
     //'"'//trim(adjustl(coffset))//'"'&
     //'/>'//end_rec
!
offset = offset+int(sizeof(x(1:n_sites))&
               +    sizeof(y(1:n_sites))&
               +    sizeof(z(1:n_sites))&
               +    sizeof(i))
!
! end points info
!
write(unit=unit_vtk,iostat=e_io)  &
     '</Points>'//end_rec
!
! start lines infos
!
write(unit=unit_vtk,iostat=e_io)  &
     '<Lines>'//end_rec
!
write(coffset,'(I10)')offset
!
write(unit=unit_vtk,iostat=e_io)       &
     '<DataArray'                      &
     //' type="Int32" '                &
     //' Name="connectivity"'          &
     //' format="appended"'            &
     //' offset='                      &
     //'"'//trim(adjustl(coffset))//'"'&
     //'/>'//end_rec
!
offset = offset                            &
        +int(sizeof(connectivity_vtk(:,:)) &
        +sizeof(i))
!
write(coffset,'(I10)')offset
!
write(unit=unit_vtk,iostat=e_io)       &
     '<DataArray'                      &
     //' type="Int32" '                &
     //' Name="offsets"'               &
     //' format="appended"'            &
     //' offset='                      &
     //'"'//trim(adjustl(coffset))//'"'&
     //'/>'//end_rec
!
offset = offset                            &
        +int(sizeof(connectivity_vtk(:,1)) & 
        +sizeof(i))
!
! end lines infos
!
write(unit=unit_vtk,iostat=e_io) &
     '</Lines>'//end_rec
!
! start verts infos
!
write(unit=unit_vtk,iostat=e_io) &
     '<Verts>'//end_rec
!
write(coffset,'(I10)')offset
!
write(unit=unit_vtk,iostat=e_io)       &
     '<DataArray'                      &
     //' type="Int32" '                &
     //' Name="connectivity"'          &
     //' format="appended"'            &
     //' offset='                      &
     //'"'//trim(adjustl(coffset))//'"'&
     //'/>'//end_rec
!
offset = offset                &
        +int(sizeof(states(:)) &
        +sizeof(i))
!
write(coffset,'(I10)')offset
!
write(unit=unit_vtk,iostat=e_io)       &
     '<DataArray'                      &
     //' type="Int32" '                &
     //' Name="offsets"'               &
     //' format="appended"'            &
     //' offset='                      &
     //'"'//trim(adjustl(coffset))//'"'&
     //'/>'//end_rec
!
offset = offset                &
        +int(sizeof(states(:)) &
        +sizeof(i))
!
! end verts infos
!
write(unit=unit_vtk,iostat=e_io)  &
     '</Verts>'//end_rec
!
! start piece infos
!
write(unit=unit_vtk,iostat=e_io)  &
     '</Piece>'//end_rec
!
! start unstructured grid infos
!
write(unit=unit_vtk,iostat=e_io)  &
     '</PolyData>'//end_rec
!
! appended data infos
!
write(unit=unit_vtk,iostat=e_io)  &
     '<AppendedData encoding="raw">'//end_rec
write(unit=unit_vtk,iostat=e_io)'_'
!
! append data 
!
! ... sites's states
!
dim_arrays = int(sizeof(states(1:n_sites)))
write(unit=unit_vtk,iostat=e_io)dim_arrays
do i = 1,n_sites
   write(unit=unit_vtk,iostat=e_io)states(i)
enddo
!
! ... sites's values
!
dim_arrays = int(sizeof(values(1:n_sites)))
write(unit=unit_vtk,iostat=e_io)dim_arrays
do i = 1,n_sites
   write(unit=unit_vtk,iostat=e_io)values(i)
enddo
!
! ... sites's invasion times
!
dim_arrays = int(sizeof(invasion_times(1:n_sites)))
write(unit=unit_vtk,iostat=e_io)dim_arrays
do i = 1,n_sites
   write(unit=unit_vtk,iostat=e_io)invasion_times(i)
enddo
!
! ... sites's trapping times
!
dim_arrays = int(sizeof(trapping_times(1:n_sites)))
write(unit=unit_vtk,iostat=e_io)dim_arrays
do i = 1,n_sites
   write(unit=unit_vtk,iostat=e_io)trapping_times(i)
enddo
!
! ... sites's coordinates
!
dim_arrays =  int( sizeof(x(1:n_sites)) &
                 + sizeof(y(1:n_sites)) &
                 + sizeof(z(1:n_sites)))
write(unit=unit_vtk,iostat=e_io)dim_arrays
do i = 1,n_sites
   write(unit=unit_vtk,iostat=e_io)x(i),y(i),z(i)
enddo
!
! ... lines
!
dim_arrays = int(sizeof(connectivity_vtk(:,:)))
write(unit=unit_vtk,iostat=e_io)dim_arrays
do i = 1,n_cells
 write(unit=unit_vtk,iostat=e_io)connectivity_vtk(i,1)-1
 write(unit=unit_vtk,iostat=e_io)connectivity_vtk(i,2)-1
enddo
!
dim_arrays = int(sizeof(connectivity_vtk(:,1)))
write(unit=unit_vtk,iostat=e_io)dim_arrays
do i = 1,n_cells
   write(unit=unit_vtk,iostat=e_io)2*i
enddo
!
! ... vertex
!
dim_arrays = int(sizeof(states(:)))
write(unit=unit_vtk,iostat=e_io)dim_arrays
do i = 1,n_sites
   write(unit=unit_vtk,iostat=e_io)i
enddo
!
dim_arrays = int(sizeof(states(:)))
write(unit=unit_vtk,iostat=e_io)dim_arrays
do i = 1,n_sites
   write(unit=unit_vtk,iostat=e_io)i
enddo
!
write(unit=unit_vtk,iostat=e_io)end_rec
!
! end appended data infos
!
write(unit=unit_vtk,iostat=e_io)  &
     '</AppendedData>'//end_rec
!
! write footer
!
write(unit=unit_vtk,iostat=e_io)  &
     '</VTKFile>'//end_rec
!
! close file
!
close(unit=unit_vtk)
!
deallocate(connectivity_vtk)
!
case('.csv')
   !
   open(unit_vtk,file=file_name,status='unknown')
   !
   write(unit_vtk,*)           &
        'x coord'       ,', ', &
        'y coord'       ,', ', &
        'z coord'       ,', ', &
        'states'        ,', ', &
        'values'        ,', ', &
        'invasion times',', ', &
        'trapping times'
   !
   do i = 1,n_sites
      write(unit_vtk,'(3(e12.4,a),i2,a,e12.4,a,i12,a,i12)')&
           x(i)             ,',', &
           y(i)             ,',', &
           z(i)             ,',', &
           states(i)        ,',', &
           values(i)        ,',', &
           invasion_times(i),',', &
           trapping_times(i)
   enddo
   !
   close(unit_vtk)
   !
case('.dat')
   !
   do i = 1,n_sites
      write(unit_vtk,'(3e12.4,i2,e12.4,i12,i12)') &
           x(i)             , &
           y(i)             , &
           z(i)             , &
           states(i)        , &
           values(i)        , &
           invasion_times(i), &
           trapping_times(i)
   enddo
   !
   close(unit_vtk)
   !
case('.bin')
   !
   open(unit       = unit_vtk,            &
        file       = file_name,           &
        form       = 'unformatted',       &
        access     = 'stream',            &
        action     = 'write',             &
        status     = 'replace',           &
        iostat     = e_io                 )
   !
write(unit_vtk)n_sites
!
do i = 1,n_sites
   write(unit_vtk)         &
        x(i)             , &
        y(i)             , &
        z(i)             , &
        states(i)        , &
        values(i)        , &
        invasion_times(i), &
        trapping_times(i)
enddo
!
close(unit_vtk)
!
case default
   print*, 'WARNING: Wrong file extension in subroutine save_arbitrary_lattice !'
   read(*,*)
end select
!
return
!
!------------------------------------
end subroutine save_arbitrary_lattice
!------------------------------------
!
!> @author Yder MASSON
!> @date October 30, 2014
!> @brief write cubic lattice info to VTK file (.vti) 
!> for viewing with e.g. Paraview
!> @details The states and values data 
!> are attached to cells (i.e. cells represent sites)
!
!-----------------------------------------------
subroutine save_cubic_lattice(states,          &
                              values,          &
                              n_sites_invaded, &
                              invasion_list,   &
                              nx,ny,nz,        &
                              dx,dy,dz,        &
                              period_x,        &
                              period_y,        &
                              period_z,        &
                              file_name,       &
                              unit_vtk         )
!-----------------------------------------------
!
implicit none
!
!> values array as defined in the invasion percolation module
real, dimension(nx,ny,nz), intent(in) :: values
!> states array as defined in the invasion percolation module 
integer, dimension(nx,ny,nz), intent(in) :: states
!
integer, intent(in) :: n_sites_invaded !< number of sites invaded
integer, dimension(:), intent(in) :: invasion_list !< list of invaded sites
!
integer, intent(in) :: nx !< grid dimension in the x direction 
integer, intent(in) :: ny !< grid dimension in the y direction
integer, intent(in) :: nz !< grid dimension in the z direction
!
real, intent(in) :: dx !< grid spacing in the x direction
real, intent(in) :: dy !< grid spacing in the y direction
real, intent(in) :: dz !< grid spacing in the z direction 
!
logical, intent(in) :: period_x !< flag for periodic boundaries in the x direction
logical, intent(in) :: period_y !< flag for periodic boundaries in the y direction
logical, intent(in) :: period_z !< flag for periodic boundaries in the z direction
!
character(len=*), intent(in) :: file_name !< Output file name, must have the.vtu extension
integer, intent(in) :: unit_vtk !< logical unit for output file
!
! internal variables
!
integer :: e_io
integer :: i, j, k, l, n
integer :: x1, x2, y1, y2, z1, z2
integer :: dim_arrays, offset
character(1), parameter:: end_rec = char(10)
character*20 :: cx1, cx2, cy1, cy2, cz1, cz2, cdx, cdy, cdz, coffset
integer, dimension(nx,ny,nz) :: invasion_times
integer, dimension(nx,ny,nz) :: trapping_times
character(len=100) :: file_name_vtk
character(len=4) :: extension
!
! get trapping times
!
call get_trapping_times_cubic(nx, ny, nz,        &
                              period_x,          &
                              period_y,          &
                              period_z,          &
                              states,            &
                              n_sites_invaded,   &
                              invasion_list,     &
                              trapping_times     )
!
! fill invasion times array
!
invasion_times(:,:,:) = -1 ! set negative time to undinvaded sites
do n = 1,n_sites_invaded
   call ind2ijk(nx,ny,invasion_list(n),i,j,k)
   invasion_times(i,j,k) = n
enddo
!
! get files extension
!
l = len_trim(file_name)
extension(1:4) = file_name(l-3:l)
!
select case(extension)
   
case('.vti')
   !
   x1 = 1
   x2 = nx
   y1 = 1
   y2 = ny
   z1 = 1
   z2 = nz
   !
   ! convert integers to strings
   !
   write(cx1,'(i10)')x1
   write(cx2,'(i10)')x2+1
   write(cy1,'(i10)')y1
   write(cy2,'(i10)')y2+1
   write(cz1,'(i10)')z1
   write(cz2,'(i10)')z2+1
   !
   write(cdx,*)dx
   write(cdy,*)dy
   write(cdz,*)dz
!
   ! open output file (add extension if not present)
!
file_name_vtk = trim(adjustl(file_name))
l = len_trim(file_name_vtk)
if(file_name_vtk(l-3:l) /= '.vti')then
   file_name_vtk(l+1:l+4) = '.vti'
endif
!
open(unit       = unit_vtk,            &
     file       = trim(file_name_vtk), &
     form       = 'unformatted',       &
     access     = 'stream',            &
     action     = 'write',             &
     status     = 'replace',           &
     convert    = 'little_endian',     &
     iostat     = e_io)
!
offset = 0
!
! write header
!
write(unit=Unit_VTK,iostat=E_IO)  &
     '<?xml version = "1.0"?>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)  &
     '<VTKFile'                  &
     //' type = "ImageData"'&
     //' version="0.1"'            &
     //' byte_order="LittleEndian"'&
     //'>'//end_rec
!
! ... grid info
!
write(unit=unit_vtk,iostat=e_io) &
     '<ImageData WholeExtent="'  &
     //trim(adjustl(cx1))//' '   &
     //trim(adjustl(cx2))//' '   &
     //trim(adjustl(cy1))//' '   &
     //trim(adjustl(cy2))//' '   &
     //trim(adjustl(cz1))//' '   &
     //trim(adjustl(cz2))//      &
     '" Origin="0 0 0"'//        &
     ' Spacing="'                &
     //trim(adjustl(cdx))//' '   &
     //trim(adjustl(cdy))//' '   &
     //trim(adjustl(cdz))//      &
     '">'//end_rec
!
write(unit=unit_vtk,iostat=e_io)  &
     '<Piece Extent="'            &
     //trim(adjustl(cx1))//' '    &
     //trim(adjustl(cx2))//' '    &
     //trim(adjustl(cy1))//' '    &
     //trim(adjustl(cy2))//' '    &
     //trim(adjustl(cz1))//' '    &
     //trim(adjustl(cz2))//'">'//end_rec
!
! ... cells states infos
!
write(unit=unit_vtk,iostat=e_io) &
     '<CellData>'//end_rec
!
write(coffset,'(I10)')offset
!
write(unit=unit_vtk,iostat=e_io)             &
     '<DataArray type="Int32" Name="states"',&
     ' format="appended" offset="'           &
     //trim(adjustl(coffset))                &
     //'" />'//end_rec
!
offset = offset             &
        +int(sizeof(states) &
        +sizeof(i))
!
write(coffset,'(I10)')offset
!
! ... cells values infos
!
write(unit=unit_vtk,iostat=e_io)               &
     '<DataArray type="Float32" Name="values"',&
     ' format="appended" offset="'             &
     //trim(adjustl(coffset))                  &
     //'" />'//end_rec
!
offset = offset             &
        +int(sizeof(values) &
        +sizeof(i))
!
write(coffset,'(I10)')offset
!
! ... cells invasion times infos
!
write(unit=unit_vtk,iostat=e_io)                     &
     '<DataArray type="Int32" Name="invasion times"',&
     ' format="appended" offset="'                   &
     //trim(adjustl(coffset))                        &
     //'" />'//end_rec
!
offset = offset                     &
        +int(sizeof(invasion_times) &
        +sizeof(i))
!
write(coffset,'(I10)')offset
!
! ... cells trapping times infos
!
write(unit=unit_vtk,iostat=e_io) &
     '<DataArray type="Int32" Name="trapping times"',&
     ' format="appended" offset="'                   &
     //trim(adjustl(coffset))                        &
     //'" />'//end_rec
!
offset = offset                      &
        +int(sizeof(trapping_times)  &
        +sizeof(i))
!
write(coffset,'(I10)')offset
!
write(unit=unit_vtk,iostat=e_io)  &
     '</CellData>'//end_rec
write(unit=unit_vtk,iostat=e_io)  &
     '</Piece>'//end_rec
write(unit=unit_vtk,iostat=e_io) &
     '</ImageData>'//end_rec
write(unit=unit_vtk,iostat=e_io) &
     '<AppendedData encoding="raw">'//end_rec
write(unit=unit_vtk,iostat=e_io)'_'
!
! write data 
!
! ... states
!
dim_arrays = int(sizeof(states))
!
write(unit=unit_vtk,iostat=e_io)dim_arrays
!
do k = z1,z2
   do j = y1,y2
      do i = x1,x2
         write(unit=unit_vtk,iostat=e_io)states(i,j,k)
      enddo
   enddo
enddo
!
! ... values
!
dim_arrays = int(sizeof(values))
write(unit=unit_vtk,iostat=e_io)dim_arrays
do k = z1,z2
   do j = y1,y2
      do i = x1,x2
         write(unit=unit_vtk,iostat=e_io)values(i,j,k)
      enddo
   enddo
enddo
!
! ... invasion times
!
dim_arrays = int(sizeof(invasion_times))
write(unit=unit_vtk,iostat=e_io)dim_arrays
do k = z1,z2
   do j = y1,y2
      do i = x1,x2
         write(unit=unit_vtk,iostat=e_io)invasion_times(i,j,k)
      enddo
   enddo
enddo
!
! ... trapping times
!
dim_arrays = int(sizeof(invasion_times))
write(unit=unit_vtk,iostat=e_io)dim_arrays
do k = z1,z2
   do j = y1,y2
      do i = x1,x2
         write(unit=unit_vtk,iostat=e_io)trapping_times(i,j,k)
      enddo
   enddo
enddo
!
write(unit=unit_vtk,iostat=e_io)end_rec
!
! write footer
!
write(unit=unit_vtk,iostat=e_io)  &
     '</AppendedData>'//end_rec
write(unit=unit_vtk,iostat=e_io)  &
     '</VTKFile>'//end_rec
!
! close file
!
close(unit=unit_vtk)
!
case('.csv')
   !
   open(unit_vtk,file=file_name,status='unknown')
   !
   write(unit_vtk,*)'x coord'       ,', ', &
                    'y coord'       ,', ', &
                    'z coord'       ,', ', &
                    'states'        ,', ', &
                    'values'        ,', ', &
                    'invasion times',', ', &
                    'trapping times'
   !
   do k = 1,nz
      do j = 1,ny
         do i = 1,nx
            write(unit_vtk,'(3(e12.4,a),i2,a,e12.4,a,i12,a,i12)')&
                 i*dx                 ,',', &
                 j*dy                 ,',', &
                 k*dz                 ,',', &
                 states(i,j,k)        ,',', &
                 values(i,j,k)        ,',', &
                 invasion_times(i,j,k),',', &
                 trapping_times(i,j,k)
         enddo
      enddo
   enddo
   !
   close(unit_vtk)
   !
!
case('.dat')
   !
   open(unit_vtk,file=trim(adjustl(file_name)),status='unknown')
   !
   do k = 1,nz
      do j = 1,ny
         do i = 1,nx
            write(unit_vtk,'(3e12.4,i2,e12.4,i12,i12)')&
                 i*dx                 , &
                 j*dy                 , &
                 k*dz                 , &
                 states(i,j,k)        , &
                 values(i,j,k)        , &
                 invasion_times(i,j,k), &
                 trapping_times(i,j,k)
         enddo
      enddo
   enddo
   !
   close(unit_vtk)
   !
case('.bin')
   !
   open(unit       = unit_vtk,            &
        file       = file_name,           &
        form       = 'unformatted',       &
        access     = 'stream',            &
        action     = 'write',             &
        status     = 'replace',           &
        iostat     = e_io                 )
   !
   write(unit_vtk)nx*ny*nz
   !
   do k = 1,nz
      do j = 1,ny
         do i = 1,nx
            write(unit_vtk)             &
                 real(i*dx)           , &
                 real(j*dy)           , &
                 real(k*dz)           , &
                 states(i,j,k)        , &
                 values(i,j,k)        , &
                 invasion_times(i,j,k), &
                 trapping_times(i,j,k)
         enddo
      enddo
   enddo
   !
   close(unit_vtk)
   !
case default
   print*, 'WARNING: Wrong file extension in subroutine save_arbitrary_lattice !'
   read(*,*)
end select
!
return
!
!--------------------------------
end subroutine save_cubic_lattice
!--------------------------------
!
!
!> @author Yder MASSON
!> @date October 30, 2014
!> @bief Produces a minimalistic 3D rendering of IP clusters inside a terninal 
!---------------------------------------
subroutine funny_3D(mat,nx,ny,nz,matval)
!---------------------------------------
!
implicit none
!
integer :: nx !< grid dimension in the x direction (\b Input)
integer :: ny !< grid dimension in the y direction (\b Input)
integer :: nz !< grid dimension in the z direction (\b Input)
integer :: mat(nx,ny,nz) !< input matrix (\b Input)
!> Value to render, i.e. only the cells where 
!> <b>mat(i,j,k) = matval</b> will be showed (\b Input)
integer :: matval 
!
! internal variables
!
integer :: i,j,k,dx,dy,dz
character(len=1), dimension(:,:), allocatable :: cmat
!
dx = nx-1
dy = ny-1
dz = nz-1
!
if(nx==1)dx=nx
if(ny==1)dy=ny
if(nz==1)dz=nz
!
allocate(cmat(2*nx+ny,nz+ny))
!
cmat(:,:) = ' '
!
do k = 1,nz
   do j = 1,ny
      do i = 1,nx
         if(mat(i,j,k)==matval)then
            cmat(2*i+j,k+j) = '*'
         endif
      enddo
   enddo
enddo
!
do k = 2,nz-1
   do j = 1,ny,dy
      do i = 1,nx,dx    
         cmat(2*i+j,k+j) = '|'
      enddo
   enddo
enddo
!
do k = 1,nz,dz
   do j = 2,ny-1
      do i = 1,nx,dx    
         cmat(2*i+j,k+j) = '/'
      enddo
   enddo
enddo
!
do k = 1,nz,dz
   do j = 1,ny,dy
      do i = 2,nx-1    
         cmat(2*i+j,k+j) = '-'
      enddo
   enddo
enddo
!
do k = 1,nz,dz
   do j = 1,ny,dy
      do i = 1,nx,dx    
         cmat(2*i+j,k+j) = '.'
      enddo
   enddo
enddo
!
print*,
do j = nz+ny+3,1,-1
   print*,(cmat(i,j),i=1,2*nx+ny)
enddo
print*,
!
deallocate(cmat)
!
return
!
!----------------------
end subroutine funny_3D
!----------------------

!------------------------------------------------------------------
subroutine save_selected_sites_arbitrary_lattice(select_state,    &
                                                 n_sites,         &
                                                 states,          &
                                                 n_sites_invaded, &
                                                 invasion_list,   &
                                                 x,y,z,           &
                                                 file_name,       &
                                                 file_unit        )
!------------------------------------------------------------------
!
implicit none
!
integer, intent(in) :: select_state !< Selection flag, 
!! all the sites where state(i)=select_state 
!! will be written to file 
!
integer, intent(in) :: n_sites !< Total number of sites in the lattice 
!> Sites's states as defiend in invasion percolation module
integer, dimension(:), intent(in) :: states
integer, intent(in) :: n_sites_invaded !< number of sites invaded
integer, dimension(:), intent(in) :: invasion_list !< list of invaded sites
!
real,dimension(:), intent(in) :: x !< Array containing the x coordinates of the sites 
real,dimension(:), intent(in) :: y !< Array containing the y coordinates of the sites 
real,dimension(:), intent(in) :: z !< Array containing the z coordinates of the sites 
!
character(len=*), intent(in) :: file_name !< Output file name, must have the.vtu extension 
integer, intent(in) :: file_unit !< logical unit for output file 
!
! internal variables
!
character(len=4) :: extension
integer :: n, l, n_points
real, dimension(:), allocatable :: x_points, y_points, z_points
!
! get files extension
!
l = len_trim(file_name)
extension(1:4) = file_name(l-3:l)
!
! find number of points to write
!
n_points = 0
!
if(select_state==invaded)then
   n_points = n_sites_invaded
else
   do n = 1,n_sites
      if(states(n)==select_state)then
         n_points = n_points+1
      endif
   enddo
endif
!
! Store points coordinates
!
allocate(x_points(n_points))
allocate(y_points(n_points))
allocate(z_points(n_points))
!
if(select_state==invaded)then
   do n = 1,n_sites_invaded
      x_points(n) = x(invasion_list(n))
      y_points(n) = y(invasion_list(n))
      z_points(n) = z(invasion_list(n))
   enddo
else
   !
   n_points = 0
   !
   do n = 1,n_sites
      if(states(n)==select_state)then
         n_points = n_points+1
         x_points(n_points) = x(n)
         y_points(n_points) = x(n)
         z_points(n_points) = x(n)
      endif
   enddo
endif
!
! write file
!
select case(extension)
   
case('.dat')
   call write_points_to_ascii(n_points,  &
                              x_points,  &
                              y_points,  &
                              z_points,  &
                              file_name, &
                              file_unit  ) 
case('.csv')
   call write_points_to_csv(n_points,  &
                            x_points,  &
                            y_points,  &
                            z_points,  &
                            file_name, &
                            file_unit  ) 
case('.vtp')
   call write_points_to_vtp(n_points,  &
                            x_points,  &
                            y_points,  &
                            z_points,  &
                            file_name, &
                            file_unit  ) 
case('.bin')
   call write_points_to_binary(n_points,  &
                               x_points,  &
                               y_points,  &
                               z_points,  &
                               file_name, &
                               file_unit  ) 
case default
   print*, 'WARNING: Wrong file extension in subroutine save_selected_sites_cubic_lattice !'
   read(*,*)
end select
!
deallocate(x_points)
deallocate(y_points)
deallocate(z_points)
!
!---------------------------------------------------
end subroutine save_selected_sites_arbitrary_lattice
!---------------------------------------------------

!--------------------------------------------------------------
subroutine save_selected_sites_cubic_lattice(select_state,    &
                                             states,          &
                                             n_sites_invaded, &
                                             invasion_list,   &
                                             nx,ny,nz,        &
                                             dx,dy,dz,        &
                                             file_name,       &
                                             file_unit        )      
!--------------------------------------------------------------
!
implicit none
!
integer, intent(in) :: select_state !< Selection flag, 
!! all the sites where state(i)=select_state 
!! will be written to file
!
!> states array as defined in the invasion percolation module 
integer, dimension(nx,ny,nz), intent(in) :: states
!
integer, intent(in) :: n_sites_invaded !< number of sites invaded
integer, dimension(:), intent(in) :: invasion_list !< list of invaded sites
!
integer, intent(in) :: nx !< grid dimension in the x direction 
integer, intent(in) :: ny !< grid dimension in the y direction 
integer, intent(in) :: nz !< grid dimension in the z direction 
!
real, intent(in) :: dx !< grid spacing in the x direction
real, intent(in) :: dy !< grid spacing in the y direction 
real, intent(in) :: dz !< grid spacing in the z direction 
!
character(len=*), intent(in) :: file_name !< Output file name, must have the.vtu extension 
integer, intent(in) :: file_unit !< logical unit for output file
!
! internal variables
!
character(len=4) :: extension
integer :: n, i, j, k, l, n_points
real, dimension(:), allocatable :: x_points, y_points, z_points
!
! get files extension
!
l = len_trim(file_name)
extension(1:4) = file_name(l-3:l)
!
! find number of points to write
!
n_points = 0
!
if(select_state==invaded)then
   n_points = n_sites_invaded
else
 do k = 1,nz
  do j = 1,ny
   do i = 1,nx
    if(states(i,j,k)==select_state)then
       n_points = n_points+1
    endif
   enddo
  enddo
 enddo
endif
!
! Store points coordinates
!
allocate(x_points(n_points))
allocate(y_points(n_points))
allocate(z_points(n_points))
!
if(select_state==invaded)then
   do n = 1,n_sites_invaded
      call ind2ijk(nx,ny,invasion_list(n),i,j,k)
      x_points(n) = (i-1)*dx
      y_points(n) = (j-1)*dy
      z_points(n) = (k-1)*dz
   enddo
else
   !
   n_points = 0
   !
   do k = 1,nz
    do j = 1,ny
     do i = 1,nx
      if(states(i,j,k)==select_state)then
        n_points = n_points+1
        x_points(n_points) = (i-1)*dx
        y_points(n_points) = (j-1)*dy
        z_points(n_points) = (k-1)*dz
      endif
     enddo
    enddo
   enddo
endif
!
! write file
!
select case(extension)
   
case('.dat')
   call write_points_to_ascii(n_points,  &
                              x_points,  &
                              y_points,  &
                              z_points,  &
                              file_name, &
                              file_unit  ) 
case('.csv')
   call write_points_to_csv(n_points,  &
                            x_points,  &
                            y_points,  &
                            z_points,  &
                            file_name, &
                            file_unit  ) 
case('.vti')
   call write_points_to_vtp(n_points,  &
                            x_points,  &
                            y_points,  &
                            z_points,  &
                            file_name, &
                            file_unit  ) 
case('.bin')
   call write_points_to_binary(n_points,  &
                               x_points,  &
                               y_points,  &
                               z_points,  &
                               file_name, &
                               file_unit  ) 
case default
   print*, 'WARNING: Wrong file extension in subroutine save_selected_sites_cubic_lattice !'
   read(*,*)
end select
!
deallocate(x_points)
deallocate(y_points)
deallocate(z_points)
!
return
!
!-----------------------------------------------
end subroutine save_selected_sites_cubic_lattice
!-----------------------------------------------
!
!> @author Yder MASSON
!> @date November 5, 2014
!> @brief Write the list of points sites to a VTK file 
!> (i.e. a polydata xml file with extension .vtp)
!> for viewing with e.g. Paraview.
!
!------------------------------------------
subroutine write_points_to_vtp(n_points,  &
                               x_points,  &
                               y_points,  &
                               z_points,  &
                               file_name, &
                               file_unit  )   
!------------------------------------------
!
implicit none
!
integer, intent(in) :: n_points !< number of points
!
real, dimension(n_points), intent(in) :: x_points !< x coordinates of the points
real, dimension(n_points), intent(in) :: y_points !< y coordinates of the points
real, dimension(n_points), intent(in) :: z_points !< z coordinates of the points
!
character (len=*), intent(in) :: file_name !< Output file name
!
integer, intent(in) :: file_unit !< Logical file unit
!
! internal variables
!
integer :: e_io
integer :: n
integer :: dim_arrays, offset
character(1), parameter:: end_rec = char(10)
character*20 :: cnp, coffset
!
! convert integer values to string
!
write(cnp,*)n_points
!
! open file
!
open(unit       = file_unit,                    &
     file       = trim(adjustl(file_name)),     &
     form       = 'unformatted',                &
     access     = 'stream',                     &
     action     = 'write',                      &
     status     = 'replace',                    &
     convert    = 'little_endian',              &
     iostat     = e_io                          )
!
! init offset for appended data
!
offset = 0
!
! write header
!
write(unit=file_unit,iostat=E_IO)  &
     '<?xml version = "1.0"?>'//end_rec
write(unit=file_unit,iostat=E_IO)  &
     '<VTKFile'                    &
     //' type = "PolyData"'        &
     //' version="0.1"'            &
     //' byte_order="LittleEndian"'&
     //'>'//end_rec
!
! start unstructured grid infos
!
write(unit=file_unit,iostat=e_io)  &
     '<PolyData>'//end_rec
!
! start piece infos
!
write(unit=file_unit,iostat=e_io)  &
     '<Piece'                      &
     //' NumberOfPoints='          &
     //'"'//trim(adjustl(cnp))//'"'&
     //'>'//end_rec
!
write(unit=file_unit,iostat=e_io) &
     '<Points>'//end_rec
!
! ... points coordinates infos
!
write(coffset,'(I10)')offset
!
write(unit=file_unit,iostat=e_io)      &
     '<DataArray'                      &
     //' type="Float32" '              &
     //' NumberOfComponents="3"'       &
     //' Name="Points"'                &
     //' format="appended"'            &
     //' offset='                      &
     //'"'//trim(adjustl(coffset))//'"'&
     //'/>'//end_rec
!
! end points info
!
write(unit=file_unit,iostat=e_io) &
     '</Points>'//end_rec
!
! end piece infos
!
write(unit=file_unit,iostat=e_io) &
     '</Piece>'//end_rec
!
! end unstructured grid infos
!
write(unit=file_unit,iostat=e_io) &
     '</PolyData>'//end_rec
!
! appended data infos
!
write(unit=file_unit,iostat=e_io) &
     '<AppendedData encoding="raw">'//end_rec
write(unit=file_unit,iostat=e_io)'_'
!
! append data 
!
dim_arrays =  int( sizeof(x_points) &
                 + sizeof(y_points) &
                 + sizeof(z_points) )
!
write(unit=file_unit,iostat=e_io)dim_arrays
!
do n = 1,n_points
   !
   write(unit=file_unit,iostat=e_io) &
        x_points(n),                 &
        y_points(n),                 &
        z_points(n)
   !
enddo
!
write(unit=file_unit,iostat=e_io)end_rec
!
! end appended data infos
!
write(unit=file_unit,iostat=e_io) &
     '</AppendedData>'//end_rec
!
! write footer
!
write(unit=file_unit,iostat=e_io) &
     '</VTKFile>'//end_rec
!
! close file
!
close(unit=file_unit)
!
return
!
!---------------------------------
end subroutine write_points_to_vtp
!---------------------------------
!
!> @author Yder MASSON
!> @date November 5, 2014
!> @brief Write the list of points sites to a .csv file 
!> (i.e. an ascii file with comma separated values)
!
!------------------------------------------
subroutine write_points_to_csv(n_points,  &
                               x_points,  &
                               y_points,  &
                               z_points,  &
                               file_name, &
                               file_unit  )   
!------------------------------------------
!
implicit none
!
integer, intent(in) :: n_points !< number of points
!
real, dimension(n_points), intent(in) :: x_points !< x coordinates of the points
real, dimension(n_points), intent(in) :: y_points !< y coordinates of the points
real, dimension(n_points), intent(in) :: z_points !< z coordinates of the points
!
character (len=*), intent(in) :: file_name !< Output file name
!
integer, intent(in) :: file_unit !< Logical file unit
!
! internal variables
!
integer :: e_io
integer :: n
character*1, parameter :: comma = ','
!
! open file
!
open(unit       = file_unit,                    &
     file       = trim(adjustl(file_name)),     &
     status     = 'unknown',                    &
     iostat     = e_io                          )
!
write(file_unit,'(a)')'x coord, y coord, z coord'
!
do n = 1,n_points
   !
   write(file_unit,'(e12.4,a,e12.4,a,e12.4)') &
        x_points(n),comma,                    &
        y_points(n),comma,                    &
        z_points(n)
   !
enddo
!
close(unit=file_unit)
!
return
!
!---------------------------------
end subroutine write_points_to_csv
!---------------------------------
!
!> @author Yder MASSON
!> @date November 5, 2014
!> @brief Write the list of points sites to an ascii file 
!
!--------------------------------------------
subroutine write_points_to_ascii(n_points,  &
                                 x_points,  &
                                 y_points,  &
                                 z_points,  &
                                 file_name, &
                                 file_unit  )   
!--------------------------------------------
!
implicit none
!
integer, intent(in) :: n_points !< number of points
!
real, dimension(n_points), intent(in) :: x_points !< x coordinates of the points
real, dimension(n_points), intent(in) :: y_points !< y coordinates of the points
real, dimension(n_points), intent(in) :: z_points !< z coordinates of the points
!
character (len=*), intent(in) :: file_name !< Output file name
!
integer, intent(in) :: file_unit !< Logical file unit
!
! internal variables
!
integer :: e_io
integer :: n
character*1, parameter :: comma = ','
!
! open file
!
open(unit       = file_unit,                &
     file       = trim(adjustl(file_name)), &
     status     = 'unknown',                &
     iostat     = e_io                      )
!
do n = 1,n_points
   !
   write(file_unit,*) &
        x_points(n),  &
        y_points(n),  &
        z_points(n)
   !
enddo
!
close(unit=file_unit)
!
return
!
!-----------------------------------
end subroutine write_points_to_ascii
!-----------------------------------
!
!> @author Yder MASSON
!> @date November 5, 2014
!> @brief Write the list of points sites to binary stream file 
!
!---------------------------------------------
subroutine write_points_to_binary(n_points,  &
                                  x_points,  &
                                  y_points,  &
                                  z_points,  &
                                  file_name, &
                                  file_unit  )   
!---------------------------------------------
!
implicit none
!
integer, intent(in) :: n_points !< number of points
!
real, dimension(n_points), intent(in) :: x_points !< x coordinates of the points
real, dimension(n_points), intent(in) :: y_points !< y coordinates of the points
real, dimension(n_points), intent(in) :: z_points !< z coordinates of the points
!
character (len=*), intent(in) :: file_name !< Output file name
!
integer, intent(in) :: file_unit !< Logical file unit
!
! internal variables
!
integer :: e_io
integer :: n
character*1, parameter :: comma = ','
!
! open file
!
open(unit       = file_unit,                &
     file       = trim(adjustl(file_name)), &
     form       = 'unformatted',            &
     access     = 'stream',                 &
     action     = 'write',                  &
     status     = 'replace',                &
     iostat     = e_io                      )
!
write(unit=file_unit)n_points
!
do n = 1,n_points
   !
   write(unit=file_unit) &
        x_points(n),     &
        y_points(n),     &
        z_points(n)
   !
enddo
!
close(unit=file_unit)
!
return
!
!------------------------------------
end subroutine write_points_to_binary
!------------------------------------
!
!===================================
end module module_write_output_files
!===================================
