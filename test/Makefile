M := 10
N := 20
target := a.out

all: subjects data

test: all
	@echo "\n SOCCS_SCE_CBLAS: reading data \n"
	@./$(target) `cat test.in`
	@echo "\n ORIGINAL_CBLAS: reading data \n"
	@./origin_$(target) `cat test.in`

subjects: $(target) origin_$(target)

data: gen
	@./$< $(M) $(N)

gen: gen.c
	gcc $^ -o $@

$(target): main.o
	gcc $^ -o $@ -L./lib/soc -lcblas

origin_$(target): main.o
	gcc $^ -o $@ -L./lib/origin -lcblas -lrefblas -lgfortran

main.o: main.c
	gcc $^ -I./include -c -o $@ 

.PHONY: clean

clean:
	rm -rf *.o *.out test.in gen
