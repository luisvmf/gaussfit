all: $(OBJS)
	gcc -o gaussfit gaussfit.c -lm
clean:
	rm gaussfit
