all:
	gcc -pthread -o out ./main.c

clean:
	rm ./out
