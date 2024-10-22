########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2010 Sergio Galindo                                    #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# any later version.                                                   #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>  #
########################################################################


cmake_minimum_required(VERSION 3.18)
project(MECHSYSSANDBOX)

# Include the UseMECHSYS.cmake module
include(${CMAKE_CURRENT_SOURCE_DIR}/UseMECHSYS.cmake)

# Set the CUDAPROGS list
set(CUDAPROGS
   wreten
 # single
# pureTest
  #tclbm05
 
)






# Check if MechSys is found
if(MECHSYS_FOUND)
    # Include the FindDEPS.cmake module
    include(${MECHSYS_SOURCE_DIR}/Modules/FindDEPS.cmake)
    
    # Check if all dependencies are satisfied
    if(NOT DEPS_OK)
        message("Not all dependencies were satisfied")
    endif()
    
    # Include MechSys include directories
    include_directories(${MECHSYS_INCLUDE_DIRS})
else()
    message("This project requires MechSys")
endif()



# Iterate over each CUDA program in CUDAPROGS
foreach(var ${CUDAPROGS})
    # Add an executable target for the CUDA program
    add_executable(${var} "${var}.cu")
    
    # Link the CUDA libraries to the target
    target_link_libraries(${var} ${CUDA_LIBRARIES})
    
    # Set compile and link flags for the target
    set_target_properties(${var} PROPERTIES
        COMPILE_FLAGS "${FLAGS} -DUSE_CUDA"
        LINK_FLAGS "${LFLAGS}"
        CUDA_ARCHITECTURES "61"
    )
    
    # Include the directory "/home/tq/mechsys"
   # include_directories(/home/tq/mechsys)
    

# Include the directory "~/mechsys"
include_directories("$ENV{HOME}/mechsys")    

    
endforeach()
