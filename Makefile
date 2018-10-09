CC = g++

CFLAGS = -O3 -std=c++11

TARGET = mW

all: $(TARGET)

$(TARGET): $(TARGET).cpp
		$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).cpp

clean:
		$(RM) $(TARGET)
