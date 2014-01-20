ifndef $(build)
	build = Release
endif

ifeq ($(build),Debug)
	THIS_CFLAGS = $(CFLAGS) -ggdb -Wall $(PG) $(PROGRESS) 
else
	THIS_CFLAGS = $(CFLAGS) $(CPU_FLAGS) $(PG) $(PROGRESS) -O3 -ffast-math -Wall
	#CFLAGS = -O2 -Wall
endif

INCPATH = -I../glm -Iinclude
LIBS =-lglm$(build) $(GLLIBS) -lpthread $(PG) $(MPI_LIBS)
THIS_LIBPATH = $(LIBPATH) -L../glm
TARGET = ../ray$(build)
OBJ := $(addprefix $(build)/, $(subst src/,, $(patsubst %.cpp, %.o, $(wildcard src/*.cpp))))

all: makedir $(OBJ)
	@rm -f $(TARGET)
	$(CC) $(THIS_CFLAGS) $(THIS_LIBPATH) -o $(TARGET) $(OBJ) $(LIBS) $(FRAMEWORK)

$(build)/%.o:src/%.cpp
	$(CC) -c $(THIS_CFLAGS) $(INCPATH) -o $@ $<

makedir:
	$(shell if test ! -d $(build); then mkdir $(build);fi)
clean:
	@rm -f $(OBJ)
	@rm -f $(TARGET)
