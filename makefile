all:	variable.o	aco.o	values.o	ant.o	main.o propagacao.o	arv.o
	g++ -o exe arv.o propagacao.o variable.o ant.o aco.o values.o main.o -lm
arv.o:	arv.h	arv.cpp
	g++ -c arv.cpp
propagacao.o: propagacao.h propagacao.cpp
	g++ -c propagacao.cpp
ant.o:	ant.cpp	ant.h
	g++ -c ant.cpp
aco.o:	aco.cpp aco.h ant.h values.h variable.h parameters.h propagacao.h
	g++ -c aco.cpp
values.o:	values.cpp	values.h
	g++ -c values.cpp
main.o: main.cpp propagacao.h variable.h aco.h
	g++ -c main.cpp
go:
	./exe
clean:
	rm *.o && rm exe