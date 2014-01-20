ifndef $(build)
	build = Release
endif

ifeq ($(build),Debug)
	CFLAGS = -w -ggdb
else 
	CFLAGS = -O6 $(CPU_FLAGS) -funroll-all-loops -w -ffast-math -fomit-frame-pointer -ftree-vectorize  
endif

AR = ar
CC = gcc
INCPATH =
LIBPATH =
TARGET = libglm$(build).a
OBJ := $(addprefix $(build)/, $(patsubst %.c, %.o, $(wildcard *.c)))

all: makedir $(OBJ)
	$(AR) -rs $(TARGET) $(OBJ)

$(build)/%.o:%.c
	$(CC) -c $(CFLAGS) -o $@ $< 

makedir:
	$(shell if test ! -d $(build); then mkdir $(build);fi)

clean:
	@rm -f $(TARGET)
	@rm -f $(OBJ)
