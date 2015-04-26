############################################################################
#
#       Project:          GenPSD makefile for gnu
#        
#       Authors:          Chris Knight October 25th 2014
#
############################################################################

TARGET      = GenPSD
SRCS      = main.cpp GenPSD.cpp integrate.cpp
OBJS      = $(SRCS:.cpp=.o)

BOOST_LIBS	= -lboost_timer -lboost_system
  
CC      = g++
COMPILE = $(CC) -std=c++11 -c $(CFLAGS)
CFLAGS      = -g -w
LINK      = $(CC) -o


$(TARGET):$(OBJS)	
	$(LINK) $(TARGET) $(OBJS) $(INCLUDES) $(OGL_LIBS) $(BOOST_LIBS)

.cpp.o:
	$(COMPILE) $(INCLUDES) $(OGL_LIBS) $(BOOST_LIBS) $< -o $@ 

all:$(TARGET)

clean:
	rm -f *.o *.x
