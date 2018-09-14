CC=gcc

all: page_rank_compressed_parallel page_rank_compressed_serial gauss_seidel test_ranks

page_rank_compressed_parallel: page_rank_compressed_parallel.c
	$(CC) -O3 -o page_rank_compressed_parallel page_rank_compressed_parallel.c -fopenmp

page_rank_compressed_serial: page_rank_compressed_serial.c
	$(CC) -O3 -o page_rank_compressed_serial page_rank_compressed_serial.c -fopenmp

gauss_seidel: gauss_seidel.c
	$(CC) -O3 -o gauss_seidel gauss_seidel.c

test_ranks: test_ranks.c
	$(CC) -O3 -o test_ranks test_ranks.c

clean:
	rm page_rank_compressed_parallel page_rank_compressed_serial gauss_seidel test_ranks