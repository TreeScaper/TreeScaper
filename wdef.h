
//##########################################################################
//# This software is part of the Treescaper i
//# -- Version 0.1
//# Copyright (C) 2010 Wen Huang
//#
//# This program is free software; you can redistribute it and/or
//# modify it under the terms of the GNU General Public License
//# as published by the Free Software Foundation; either version 2
//# of the License, or (at your option) any later version.
//#
//# This program is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//# http://www.gnu.org/copyleft/gpl.html
//##########################################################################

// define the parameter for differnt platform

#ifndef WDEF_H
#define WDEF_H

//#define COMMAND_LINE_VERSION

#ifndef COMMAND_LINE_VERSION

    #ifdef _WIN64 // The following code is compiled only when this library is compiled in Windows (64-bit only)
        #define _WINDOWS
    #elif _WIN32 // The following code is compiled only when this library is compiled in Windows (both 32-bit and 64-bit only)
       //define something for Windows (32-bit and 64-bit, this part is common)
        #define _WINDOWS
    #elif __APPLE__ // The following code is compiled only when this library is compiled in MAC
        #define _MAC
    #elif __linux// The following code is compiled only when this library is compiled in Linux system
        #define _LINUX
    #elif __unix // all unices not caught above
        // Unix
        #define _LINUX
    #elif __posix
        // POSIX
    #endif // end of checking platforms

#endif


#endif // WDEF_H
