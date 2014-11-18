!===============================================!
!                                               !
!                                               !
!     MODULE INVASION PERCOLATION CONSTANTS     !
!                                               !
!                                               !
!===============================================!
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
!> @author Yder Masson
!> @date October 22, 2014
!> @brief This module is used to define and share some constants 
!
!===========================================
module module_invasion_percolation_constants
!===========================================

INTEGER, PARAMETER :: NOT_INVADED = 0 !< State flag
INTEGER, PARAMETER :: NEIGHBORING = 1 !< State flag
INTEGER, PARAMETER :: TRAPPED     = 2 !< State flag
INTEGER, PARAMETER :: EXIT_SITE   = 3 !< State flag
INTEGER, PARAMETER :: INVADED     = 4 !< State flag

!===============================================
end module module_invasion_percolation_constants
!===============================================
!
!=======================
!                      !
!                      !
!     iPerc MANUAL     !
!                      !
!                      !
!======================!
!
!> @mainpage iPerc Manual
!> @author Yder Masson
!> @date October 22, 2014
!> @section intro_sec Introduction
!> \b iPerc is a software suite for modeling invasion percolation as introduced by Wilkinson and Willemsen in 1983 <em>(WILKINSON, David et WILLEMSEN, Jorge F. Invasion percolation: a new form of percolation theory. Journal of Physics A: Mathematical and General, 1983, vol. 16, no 14, p. 3365.)</em>. The code is written in Fortran 2003 and implement fast algorithms for simulating invasion percolation on arbitrary lattices. Both gravity and trapping can be modeled. This software explicitely model site percolation but it can also be used to model bond percolation. Some additional tools for generating random media and for visualization are also part of this package.  
!>
!> @section ref_sec References
!>
!> If you use this code for your own research, please cite the following articles written by the developers of the package:
!>
!>  [1] <em>MASSON, Yder et PRIDE, Steven R. A fast algorithm for invasion percolation. Transport in porous media, 2014, vol. 102, no 2, p. 301-312.</em>
!> \n \n [2] <em>MASSON, Yder et PRIDE, Steven R. A fast algorithm for invasion percolation Part II: Efficient a posteriory treatment of trapping. (To be published).</em>
!>
!> @section license_sec License
!>
!>    iPerc is a fortran library for modeling invasion percolation
!>    Copyright (C) 2014  Yder MASSON
!>
!>    \n iPerc is free software: you can redistribute it and/or modify
!>    it under the terms of the GNU General Public License as published by
!>    the Free Software Foundation, either version 3 of the License, or
!>    (at your option) any later version.
!>
!>    \n This program is distributed in the hope that it will be useful,
!>    but WITHOUT ANY WARRANTY; without even the implied warranty of
!>    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!>    GNU General Public License for more details.
!>
!>    \n You should have received a copy of the GNU General Public License
!>    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!>
!> @section getting_started_sec Getting started
!> @subsection prerequisites_sec Prerequisites
!> * You should have received a copy of the source code: <b>iPerc.1.0.tar.gz</b>
!> \n \n * You must have the GNU Fortran compiler \b gfortran installed.
!> \n \n * For a better visual experience, you may want to install the ParaView visualisation software or any 3D plotter using the Visual toolkit (e.g. Check Mayavi if you are a Python addict). <b>(This is optional)</b>
!> @note You can of course compile \b iPerc using another Fortran compiler. In this case, you have to edit the Makefile located in the  \b iPerc/ directory. Replace \b gfortran with your compiler (i.e. FC = your_fortran_compiler) and make sure to change the compiling options accordingly (i.e. FCFLAG= your_compiler_options, FCFLAG+= your_compiler_options).
!> @see http://www.paraview.org/
!> @see http://mayavi.sourceforge.net/
!> @see http://en.wikipedia.org/wiki/Gfortran
!>
!> \subsection opening_the_box_sec Opening the box...
!> Unzip the archive:
!> \code tar -zxvf iPerc.1.0.tar.gz \endcode
!> Move to the main directory:
!> \code cd iPerc/ \endcode
!> Compile the source code:
!> \code make \endcode
!> This will compile the iPerc library as well as the examples in the <b>iPerc/examples/src/</b> directory and the projects in the <b>iPerc/my_project/src/</b>.
!> Then, you can try to run the examples in the <b>iPerc/examples/bin/</b> directory>, for example, in the iPerc directory \b iPerc/, type:
!> \code ./examples/bin/name_of_the_example_you_want_to_run.exe \endcode
!> Once you have successfully run the examples, you can move on to the next section and start building your own project !
!> @section project Building, compiling and running new projects
!> @subsection new_project Create your project
!> Any new projet should be placed in the <b>iPerc/my_project/src/</b> directory and have the <b>.f90</b> extension (no upper case please, i.e., no <b>.F90</b> extension). For example, create the file:
!> \code
!> iPerc/my_project/src/the_name_of_your_project.f90
!> \endcode
!> or, find one example that does something close to what you want, rename it and place it in the <b>iPerc/my_project/src</b> directory. In the <b>iPerc/</b> directory type:
!> \code
!> cp ./examples/src/the_example_you_like.f90 ./my_project/src/the_name_of_your_project.f90
!> \endcode
!> When your project has been created, you can edit it using your favorite text editor (e.g. emacs, vim, etc...). Any iPerc project should have the following basic structure:
!> \code
!> !========!
!> ! HEADER !
!> !========!
!>
!> ! Lines starting with a exclamation mark are comment lines,
!> ! these do not need to be present in your code.
!> ! Your code project starts with program
!> ! followed by the name of your project: 
!>
!> program name_of_your_project
!>
!> ! In order to use iPerc, the following statment must be present 
!> ! at the begining of your code project:
!>
!> use module_invasion_percolation
!>
!> ! Please put the following line in any Fortran code your write !
!> ! This will save you a lot of trouble ;)
!>
!> implicit none
!> 
!> !===============!
!> ! DERCLARATIONS !
!> !===============!
!>
!> ! Delare your variables here, e.g.:
!>
!> integer :: some, integers
!> real :: a, few, reals
!> logical :: etc
!>
!> !==============!
!> ! INSTRUCTIONS !
!> !==============!
!>
!> ! This is where you write your code
!> ! See the following sections for more details
!> 
!> print*, 'This is my first iPerc project !'
!>
!> !========!
!> ! FOOTER !
!> !========!
!>
!> ! Your code project ends with end program
!> ! followed by the name of your project: 
!>
!> end program name_of_your_project
!> \endcode
!>
!> @subsection compile_project Compiling and running your project
!> To compile your projects, move to the \b iPerc/ directory and type:
!> \code make my_project \endcode
!> Or, if you want to recompile the whole iPerc library, type:
!> \code make \endcode
!> This will create the files: 
!> \code iPerc/my_project/bin/the_name_of_your_projects.exe \endcode
!> To run your project, type: 
!> \code ./my_project/bin/the_name_of_your_project.exe \endcode
!> As an exercise, you can create a project containing the code in the previous section and run it. The following message will print on sreen:
!> \code This is my frst iPerc project ! \endcode
!> @section modeling_sec Modeling invasion percolation
!>
!> \code
!> integer, dimension(:,:), allocatable i
!> call find_trapped_sites_arbitrary(n_sites,           &
!>                                   states,            &
!>                                   offsets,           &
!>                                   connectivity,      &
!>                                   n_sites_invaded,   &
!>                                   invasion_list,     &
!>                                   undo_invasion      )
!>
!> \endcode
!!  @code{.f90}
!!  class Python:
!!     pass
!!  @endcode
!>
!> lkjdl jk lkj

