# /******************************************************************************
#  *   Copyright (C) 2024  Pablo Catalan                                        *
#  *                                                                            *
#  *  This program is free software: you can redistribute it and/or modify      *
#  *  it under the terms of the GNU General Public License as published by      *
#  *  the Free Software Foundation, either version 3 of the License, or         *
#  *  (at your option) any later version.                                       *
#  *                                                                            *
#  *  This program is distributed in the hope that it will be useful,           *
#  *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
#  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
#  *  GNU General Public License for more details.                              *
#  *                                                                            *
#  *  You should have received a copy of the GNU General Public License         *
#  *  along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
#  ******************************************************************************/

############## Operating-System and Processor Architecture Detection #########
# (many thanks to Juan Antonio García-Martín for this part of the code)
# This section attempts to auto-detect the OS. This is used only for convenience
# and allows the same Makefile to be used on multiple operating systems
# without modification.
# You can circumvent auto-detection by setting these environment variables:
#   TOY_MAKE_OS   -- the operating system name (e.g. Linux, Windows, Mac)
#   TOY_MAKE_ARCH -- the target architecture (e.g. x86 or x86_64)
# (These can be set as environment variables or on the MAKE command-line. 
#    e.g. `make all TOY_MAKE_OS=Linux`)
    ifneq (${TOY_MAKE_OS},) 
      #if TOY_MAKE_OS is NOT blank, use it as the OS
      OPSYSTEM=$(TOY_MAKE_OS)
    else ifeq (${OPSYSTEM},)
      # IF both TOY_MAKE_OS and OPSYSTEM are blank, use the `uname` command 
      #   (if available) to determine the OS.
      #   Replace 'UNKNOWN' with default OS if desired. 
      OPSYSTEM_RAW:=$(shell uname -s 2>/dev/null || echo UNKNOWN) 
      # Perform some replacements to normalize the output of uname on various systems.
      # OS_REPLACEMENTS= CYGWIN%=Windows MSYS%=Windows Darwin=Mac GNU%=Linux
      OPSYSTEM := $(OPSYSTEM_RAW:CYGWIN%=Windows)
      OPSYSTEM := $(OPSYSTEM:MSYS%=Windows)
      OPSYSTEM := $(OPSYSTEM:Darwin=Mac)
      OPSYSTEM := $(OPSYSTEM:GNU%=Linux)
      $(if $(DEBUG), $(info Make: Operating System: $(OPSYSTEM)))
      export OPSYSTEM #make it available for recursive calls to make, so auto-detection is performed only once.
    endif


CXX     	= g++
CXXFLAGS        = --std=c++17 -fPIC -O4 -DNDEBUG -DDSFMT_MEXP=19937 -DHAVE_SSE2
    ifeq (${OPSYSTEM},Linux)
      ############ LINUX ##########################################################
      LDFLAGS         = -lstdc++ -lgcc -lm
    else ifeq (${OPSYSTEM},Mac)
      ############ MAC ############################################################
      LDFLAGS         = -lstdc++ -lgcc -lm
    else ifeq (,${OPSYSTEM})
      ############ NO OS ########################################################
	  $(error No Operating system defined!!!)
    else ifeq (${OPSYSTEM},Windows)
      ############ WINDOWS ########################################################
      LDFLAGS         = -lstdc++ -lgcc -lm
   else
	  ############ UNKNOWN OS ###################################################
      $(error Unknown Operating system defined: $(OPSYSTEM))
    endif

main: main.o FossilRecord.o toylife/helper_functions.o toylife/toy_plugin.o func.o dSFMT.o Population.o Prokaryote.o toyAgent.o
	$(CXX) $(CXXFLAGS) main.o FossilRecord.o toylife/helper_functions.o toylife/toy_plugin.o func.o dSFMT.o Population.o Prokaryote.o toyAgent.o $(LDFLAGS) -o main

main.o: main.cc
	$(CXX) $(CXXFLAGS) -c main.cc -g -o main.o

classes.o: FossilRecord.cc
	$(CXX) $(CXXFLAGS) -c FossilRecord.cc -g -o FossilRecord.o

toylife/helper_functions.o: toylife/helper_functions.cpp
	$(CXX) $(CXXFLAGS) -c toylife/helper_functions.cpp -g -o toylife/helper_functions.o

toylife/toy_plugin.o: toylife/toy_plugin.cpp
	$(CXX) $(CXXFLAGS) -c toylife/toy_plugin.cpp -g -o toylife/toy_plugin.o

func.o: func.cpp
	$(CXX) $(CXXFLAGS) -c func.cpp -g -o func.o

dSFMT.o: dSFMT.cc
	$(CXX) $(CXXFLAGS) -c dSFMT.cc -g -o dSFMT.o

Population.o: Population.cc
	$(CXX) $(CXXFLAGS) -c Population.cc -g -o Population.o

Prokaryote.o: Prokaryote.cc
	$(CXX) $(CXXFLAGS) -c Prokaryote.cc -g -o Prokaryote.o

toyAgent.o: toyAgent.cc
	$(CXX) $(CXXFLAGS) -c toyAgent.cc -g -o toyAgent.o

  # Test for memory leaks using Valgrind
memcheck: main
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose ./main

clean:
	rm -rf *.o
	rm -rf toylife/*.o
	rm -rf main
