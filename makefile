#############################################################################

# Makefile for building mesh2d
# Generated by tmake at 01:30, 2007/04/20
#     Project: tmake
#    Template: app
#############################################################################

####### Compiler, tools and options

QTDIR	=	/usr
CC	=	gcc
CXX	=	g++
CFLAGS	=	-pipe -Wall -W -O0 -g3
CXXFLAGS=	-pipe -Wall -W -O0 -g3
INCPATH	=	-I.
LINK	=	g++
LFLAGS	=
LIBS	=	$(SUBLIBS)
MOC	=
UIC	=

TAR	=	tar -cf
GZIP	=	gzip -9f

####### Files

HEADERS =	bc2d.h \
		common.h \
		efread.h \
		getpot.h \
		stopwatch.h \
		t_mesh2d_aux_funcs.h \
		t_mesh2d.h
SOURCES =	efread.cpp \
		front_from_file.cpp \
		t_mesh2d_dump.cpp \
		t_mesh2d_fist.cpp \
		t_mesh2d_gen.cpp \
		t_mesh2d_stream.cpp
OBJECTS =	efread.o \
		front_from_file.o \
		t_mesh2d_dump.o \
		t_mesh2d_fist.o \
		t_mesh2d_gen.o \
		t_mesh2d_stream.o
INTERFACES =
UICDECLS =
UICIMPLS =
SRCMOC	=
OBJMOC	=
DIST	=
TARGET	=	mesh2d_V2
INTERFACE_DECL_PATH = .

####### Implicit rules

.SUFFIXES: .cpp .cxx .cc .C .c

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o $@ $<

####### Build rules


all: $(TARGET)

$(TARGET): $(UICDECLS) $(OBJECTS) $(OBJMOC)
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJMOC) $(LIBS)

moc: $(SRCMOC)

tmake:
	tmake tmake.pro

dist:
	$(TAR) tmake.tar tmake.pro $(SOURCES) $(HEADERS) $(INTERFACES) $(DIST)
	$(GZIP) tmake.tar

clean:
	-rm -f $(OBJECTS) $(OBJMOC) $(SRCMOC) $(UICIMPLS) $(UICDECLS) $(TARGET)
	-rm -f *~ core

####### Sub-libraries


###### Combined headers


####### Compile

efread.o: efread.cpp \
		efread.h

front_from_file.o: front_from_file.cpp \
		stopwatch.h \
		t_mesh2d.h \
		efread.h \
		common.h \
		t_mesh2d_aux_funcs.h \
		bc2d.h \
		getpot.h

t_mesh2d_dump.o: t_mesh2d_dump.cpp \
		t_mesh2d.h \
		efread.h \
		common.h \
		t_mesh2d_aux_funcs.h

t_mesh2d_fist.o: t_mesh2d_fist.cpp \
		t_mesh2d.h \
		efread.h \
		common.h \
		t_mesh2d_aux_funcs.h

t_mesh2d_gen.o: t_mesh2d_gen.cpp \
		t_mesh2d.h \
		efread.h \
		common.h \
		t_mesh2d_aux_funcs.h

t_mesh2d_stream.o: t_mesh2d_stream.cpp \
		t_mesh2d.h \
		efread.h \
		common.h \
		t_mesh2d_aux_funcs.h

