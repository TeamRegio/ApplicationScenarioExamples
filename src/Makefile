#sources
SRC =  	PASTAA.cpp \
	PSCM_to_PSEM.cpp \
	TRAP_normalized.cpp \

platform=$(shell uname -s)

ifeq ($(platform), Darwin)

CC = clang++
#CFLAGS =  -std=c++11  -Xpreprocessor -fopenmp
#CFLAGS =  -std=c++11  -stdlib=libc++ 
LDFLAGS = -stdlib=libc++

all: $(SRC) Makefile
	$(CC) PASTAA.cpp  $(LDFLAGS) -o PASTAA
	$(CC) PSCM_to_PSEM.cpp  $(LDFLAGS) -o PSCM_to_PSEM
	$(CC) TRAP_normalized.cpp  $(LDFLAGS)  -o TRAP_normalized

else

CC = g++ 
#CFLAGS = -std=c++11 -fopenmp 
LDFLAGS = -static

all: $(SRC) Makefile
	$(CC) PASTAA.cpp  $(LDFLAGS) -o PASTAA
	$(CC) PSCM_to_PSEM.cpp  $(LDFLAGS) -o PSCM_to_PSEM
	$(CC) TRAP_normalized.cpp  $(LDFLAGS)  -o TRAP_normalized


endif
	
clean:
	-/bin/rm -f core *.o *~

