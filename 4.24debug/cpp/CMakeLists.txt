########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             #
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



PROJECT                (tes)
CMAKE_MINIMUM_REQUIRED (VERSION 2.8)
### DEPENDENCIES ########################################################################################

# essential dependencies
INCLUDE (${tes_SOURCE_DIR}/FindMECHSYS.cmake)
#message("idiota " ${MECHSYS_SOURCE_DIR})
if(MECHSYS_FOUND)
    INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindDEPS.cmake)
    if(NOT DEPS_OK)
#        MESSAGE("Not all dependencies were satisfied")
    endif(NOT DEPS_OK)
    INCLUDE_DIRECTORIES (${MECHSYS_INCLUDE_DIRS})
else(MECHSYS_FOUND)
    MESSAGE("This project requires MechSys")
endif(MECHSYS_FOUND)


SET(PROGS
   # colloid
   # hyratfix
   # test_write
   # graflow
   # bubble
    #soil3D
   # soilana
   # streamline
   # alexgflow
   # alexDEM
   # multicomp
   # btransport
   # fracking
    wreten)

FOREACH(var ${PROGS})
    ADD_EXECUTABLE        (${var} "${var}.cpp")
    TARGET_LINK_LIBRARIES (${var} ${LIBS})
    SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}")
ENDFOREACH(var)
