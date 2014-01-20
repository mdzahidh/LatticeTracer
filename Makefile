#I dont know how to automatically figure out the CPU model in Mac computers
#so to build in Mac please specify the cpu model manually here.
CC = g++
CPU_FLAGS="-march=native"

PG =
PROGRESS = -DPRINT_PROGRESS 

MPI = 1
MPI_TRACE = 0
MPI_HEADER_DIR = -I/opt/local/include/mpich2
MPI_LIB_DIR = -L/opt/local/lib
MPI_LIBS = -lmpich

#ifneq ($(shell cat /proc/cpuinfo | grep -i "Core(TM)2" | wc -l),0)
#	CPU_FLAGS=-march=core2
#else ifneq ($(shell cat /proc/cpuinfo | grep -i "Opteron" | wc -l),0)
#	CPU_FLAGS=-march=opteron
#endif

ifeq ($(MPI),1)
	CC=mpic++
	CFLAGS+=${MPI_HEADER_DIR}
	LIBPATH="$(MPI_LIB_DIR)"
else
	CC=g++
	MPI_LIBS=#I dont know how to automatically figure out the CPU model in Mac computers
#so to build in Mac please specify the cpu model manually here.
CC = g++
CPU_FLAGS="-march=native"

PG =
PROGRESS = -DPRINT_PROGRESS 

MPI = 1
MPI_TRACE = 1
MPI_HEADER_DIR = -I/opt/local/include/mpich2
MPI_LIB_DIR = -L/opt/local/lib
MPI_LIBS = -lmpich

#ifneq ($(shell cat /proc/cpuinfo | grep -i "Core(TM)2" | wc -l),0)
#	CPU_FLAGS=-march=core2
#else ifneq ($(shell cat /proc/cpuinfo | grep -i "Opteron" | wc -l),0)
#	CPU_FLAGS=-march=opteron
#endif

ifeq ($(MPI),1)
	CC=mpic++
	CFLAGS+=${MPI_HEADER_DIR}
	LIBPATH="$(MPI_LIB_DIR)"
else
	CC=g++
	MPI_LIBS=
endif

ifeq ($(MPI_TRACE),1)	
	MPI_LIBS+=-lmpe -llmpe
endif



PLATFORM = $(shell uname -s)
ifeq ($(PLATFORM),Darwin)
	FRAMEWORK = -framework GLUT
	FRAMEWORK += -framework OpenGL
else
	FRAMEWORK =
	GLLIBS = -lGL -lGLU -lglut
endif


ifndef $(build)
	build = Release
endif


all:
	make -C glm CC="$(CC)" CFLAGS="$(CFLAGS)" LIBPATH="$(LIBPATH)" MPI_LIBS="$(MPI_LIBS)" CPU_FLAGS="$(CPU_FLAGS)" FRAMEWORK="$(FRAMEWORK)" GLLIBS="$(GLLIBS)" 
	make -C RayTracer CC="$(CC)" CFLAGS="$(CFLAGS)" LIBPATH="$(LIBPATH)" MPI_LIBS="$(MPI_LIBS)" CPU_FLAGS="$(CPU_FLAGS)" FRAMEWORK="$(FRAMEWORK)" GLLIBS="$(GLLIBS)" PG="$(PG)" PROGRESS="$(PROGRESS)"

clean:
	make -C RayTracer clean
	
clean_all:
	make -C glm clean
	make -C RayTracer clean



endif

ifeq ($(MPI_TRACE),1)	
	MPI_LIBS+=-lmpe -llmpe
endif



PLATFORM = $(shell uname -s)
ifeq ($(PLATFORM),Darwin)
	FRAMEWORK = -framework GLUT
	FRAMEWORK += -framework OpenGL
else
	FRAMEWORK =
	GLLIBS = -lGL -lGLU -lglut
endif


ifndef $(build)
	build = Release
endif


all:
	make -C glm CC="$(CC)" CFLAGS="$(CFLAGS)" LIBPATH="$(LIBPATH)" MPI_LIBS="$(MPI_LIBS)" CPU_FLAGS="$(CPU_FLAGS)" FRAMEWORK="$(FRAMEWORK)" GLLIBS="$(GLLIBS)" 
	make -C RayTracer CC="$(CC)" CFLAGS="$(CFLAGS)" LIBPATH="$(LIBPATH)" MPI_LIBS="$(MPI_LIBS)" CPU_FLAGS="$(CPU_FLAGS)" FRAMEWORK="$(FRAMEWORK)" GLLIBS="$(GLLIBS)" PG="$(PG)" PROGRESS="$(PROGRESS)"

clean:
	make -C RayTracer clean
	
clean_all:
	make -C glm clean
	make -C RayTracer clean


