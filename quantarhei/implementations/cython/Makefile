CFLAGS=`python3-config --cflags`
LDFLAGS=`python3-config --ldflags`

all: loopit.o 
	@echo ${CFLAGS}
	@echo ${LDFLAGS}


loopit.o: loopit.c
	gcc -c loopit.c ${CFLAGS}

loopit.so: loopit.o
	gcc loopit.o -o loopit.so -shared ${LDFLAGS} 


