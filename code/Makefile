N=50

CC=gcc
#CC=icc

CFLAGS= -c -g  -mavx -Wall -Wshadow  #-Wconversion
#CFLAGS= -c -O2 -mavx  -Wall -Wshadow -finline-functions -funroll-loops

####################################################################################

DEFINES= -DN=$(N) #-DDEBUG 


TARGET=stressed_mf_gille_N$(N) 

LARGS= -lm -lrt

OBJS = stress_mf_gille_main.o dsvdfit.o


$(TARGET): $(OBJS)
	$(CC) $(OBJS) $(LFLAGS) $(LARGS) -g -o $(TARGET)

.c.o:
	$(CC) -c $< $(DEFINES) $(CFLAGS)

clean:
	/bin/rm -f $(OBJS) *~

tar:
	tar -czvf stress_gille.tar.gz stress_mf_gille*.[hc]  Makefile

tarall:
	tar -czvf  stress_gille.tar.gz *.c *.h *.tex *.mak *.ps input.* betas.* Makefile *.sh

#		*Individual File Dependencies*
stress_mf_gille.o: stress_mf_gille_main.c stress_gille.h dsvdfit.o   Makefile 


