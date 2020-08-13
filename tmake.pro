HEADERS   = bc2d.h  common.h  efread.h  getpot.h  stopwatch.h  t_mesh2d_aux_funcs.h  t_mesh2d.h

SOURCES   = efread.cpp  front_from_file.cpp  t_mesh2d_dump.cpp  t_mesh2d_fist.cpp  t_mesh2d_gen.cpp  t_mesh2d_stream.cpp

TARGET    = mesh2d_V2

INCLUDEPATH = .

DEFINES =

DEPENDPATH =

CONFIG = release warn_on

LIBS =

TMAKE_CFLAGS_RELEASE	= -O3 -g0
TMAKE_CFLAGS_DEBUG	  = -g0 -O3

TMAKE_CXXFLAGS_RELEASE	= -O3 -g0
TMAKE_CXXFLAGS_DEBUG	  = -g0 -O3

TMAKE_MOC		=
TMAKE_UIC		=
