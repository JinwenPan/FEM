CC = g++
CFLAGS  = -std=c++17 -pedantic
TARGET = waveguide

all: main.cpp 
	$(CC) $(CFLAGS) -g main.cpp -o $(TARGET)