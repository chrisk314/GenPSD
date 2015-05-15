############################################################################
#
#       Project:          GenPSD makefile for gnu
#        
#       Authors:          Chris Knight October 25th 2014
#
############################################################################

TARGET      = GenPSD
SRCS	:= $(wildcard src/*.cpp)
OBJS	:= $(addprefix obj/,$(notdir $(SRCS:.cpp=.o)))
BOOST_LIBS	= -lboost_timer -lboost_system
  
CC      = g++
COMPILE = $(CC) $(CFLAGS) -c
CFLAGS      = -g -w -std=c++11
LINK      = $(CC) -o


$(TARGET):$(OBJS)	
	$(LINK) $(TARGET) $(OBJS) $(INCLUDES) $(BOOST_LIBS)

obj/%.o: src/%.cpp
	$(COMPILE) $(INCLUDES) $(BOOST_LIBS) -o $@ $<

all:$(TARGET)

clean:
	rm -f obj/*.o $(TARGET)
