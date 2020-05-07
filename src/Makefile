#sources
SRC =  	PASTAA.cpp \
	PSCM_to_PSEM.cpp \
	TRAP.cpp \

platform=$(shell uname -s)

ifeq ($(platform), Darwin)

CC = clang++
CFLAGS =  -std=c++11  -Xpreprocessor -fopenmp
#CFLAGS =  -std=c++11  -stdlib=libc++ 
LDFLAGS = -stdlib=libc++

all: $(SRC) Makefile
	$(CC) PASTAA.cpp ${CFLAGS}  $(LDFLAGS) -o PASTAA
	$(CC) PSCM_to_PSEM.cpp ${CFLAGS} $(LDFLAGS) -o PSCM_to_PSEM
	$(CC) TRAP.cpp ${CFLAGS} $(LDFLAGS)  -o TRAP

else

#CC = clang++ 
CC = g++ 
CFLAGS = -std=c++11 -fopenmp 
LDFLAGS = -static
#WFLAG = -Wall #for warnings

all: $(SRC) Makefile
	$(CC) PASTAA.cpp -g ${CFLAGS}  $(LDFLAGS)  -o PASTAA
	$(CC) PSCM_to_PSEM.cpp ${CFLAGS} $(LDFLAGS) -o PSCM_to_PSEM
	$(CC) TRAP.cpp ${CFLAGS} $(LDFLAGS)  -o TRAP


endif
	
clean:
	-/bin/rm -f core *.o *~

