# ============================================================================
#  Makefile HBTeg                             Chris Plumberg, October 31, 2018
# ============================================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install	make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  

CC := g++
CFLAGS= -O3 -fopenmp -std=c++11 -g

RM		=	rm -f
O               =       .o
LDFLAGS         =       $(CFLAGS)
SYSTEMFILES     =       $(SRCGNU)
SRCDIR		=	src

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	SV.e
endif

SRC		=	main.cpp \
			$(SRCDIR)/Arsenal.cpp \
			$(SRCDIR)/ParameterReader.cpp \
			$(SRCDIR)/gauss_quadrature.cpp \
			$(SRCDIR)/Compute_radii_qmode_3D.cpp \
			$(SRCDIR)/SourceVariances.cpp

INC		= 	main.h \
			$(SRCDIR)/EventRecord.h \
			$(SRCDIR)/ParticleRecord.h \
			$(SRCDIR)/Stopwatch.h \
			$(SRCDIR)/Arsenal.h \
			$(SRCDIR)/ParameterReader.h \
			$(SRCDIR)/gauss_quadrature.h \
			$(SRCDIR)/ensemble.h \
			$(SRCDIR)/matrix.h \
			$(SRCDIR)/random_events.h \
			$(SRCDIR)/read_in_data.h \
			$(SRCDIR)/estimate_error.h \
			$(SRCDIR)/SourceVariances.h

# -------------------------------------------------

OBJDIR		=	.
SRCFILES 	= 	$(SRC) $(INC) GNUmakefile_HBTeg
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC))))
TARGET		=	$(MAIN)
INSTPATH	=	..
#INSTPATH	=	$(HOME)/local/bin

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

%.cpp:
	if [ -f $@ ] ; then touch $@ ; else false ; fi

# -------------------------------------------------

.PHONY:		all mkobjdir clean distclean install

all:		mkobjdir $(TARGET)

help:
		@grep '^##' GNUmakefile_HBTeg

mkobjdir:	
		-@mkdir -p $(OBJDIR)

$(TARGET):	$(OBJECTS)	
		$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS) 
#		strip $(TARGET)

clean:		
		-rm $(OBJECTS)

distclean:	
		-rm $(TARGET)
		-rm $(OBJECTS)
#		-rm $(MAIN)

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
$(SRCDIR)/Arsenal.cpp: $(SRCDIR)/Arsenal.h
$(SRCDIR)/ParameterReader.cpp: $(SRCDIR)/ParameterReader.h $(SRCDIR)/Arsenal.h
$(SRCDIR)/HBT_event_generator.cpp: $(SRCDIR)/HBT_event_generator.h $(SRCDIR)/Arsenal.h \
									$(SRCDIR)/Stopwatch.h $(SRCDIR)/EventRecord.h \
									$(SRCDIR)/ParticleRecord.h $(SRCDIR)/ensemble.h
./main.cpp: $(SRCDIR)/Stopwatch.h $(SRCDIR)/ParameterReader.h main.h $(SRCDIR)/EventRecord.h \
			$(SRCDIR)/ParticleRecord.h $(SRCDIR)/ensemble.h

