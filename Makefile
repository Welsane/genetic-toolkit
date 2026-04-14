# ============================================================
#  Genetic Analyzer Toolkit — Makefile
#  Targets:
#    make          → build (default)
#    make clean    → remove compiled output
#    make run      → build + run
# ============================================================

CC      = gcc
CFLAGS  = -Wall -Wextra -std=c11 -pedantic

# Output binary name (add .exe on Windows automatically via uname check)
UNAME := $(shell uname 2>/dev/null || echo Windows)
ifeq ($(UNAME), Windows)
    TARGET = gat.exe
else
    TARGET = gat
endif

SRCS    = main.c dna_utils.c file_handler.c
OBJS    = $(SRCS:.c=.o)
HEADERS = dna_utils.h file_handler.h

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^
	@echo ""
	@echo "  Build successful → $(TARGET)"
	@echo "  Run with: ./$(TARGET)"
	@echo ""

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	@echo "  Cleaning build artifacts..."
	rm -f $(OBJS) $(TARGET)

run: $(TARGET)
	./$(TARGET)
