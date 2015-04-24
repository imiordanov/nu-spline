
CC=g++

all: nuspline

nuspline: 
	$(CC) driver.cpp -o nuspline

clean: 
	rm -f nuspline spline points knots a.out