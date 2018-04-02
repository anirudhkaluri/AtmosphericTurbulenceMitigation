CFLAGS=-std=c99 -Wall -Wextra -Wno-unused -pedantic -O3 -Wno-unknown-pragmas

ifdef OPENMP
	CFLAGS=-std=c99 -Wall -Wextra -Wno-unused -Wno-uninitialized -pedantic -fopenmp -O3
endif

main: main_MaoGilles.c
	$(CC) $(CFLAGS) -DNDEBUG -D_GNU_SOURCE -c iio.c
	$(CC) $(CFLAGS) -o MaoGilles main_MaoGilles.c iio.o -lm -lpng -ljpeg -ltiff
	$(CC) $(CFLAGS) -o ShiftMaoGilles main_ShiftMaoGilles.c iio.o -lm -lpng -ljpeg -ltiff
	rm iio.o

clean:
	rm -rf *.o *~ MaoGilles ShiftMaoGilles srcdoc/

srcdoc:
	doxygen doxygen.conf
