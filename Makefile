bin/main: src/main.cpp src/gnuplot-iostream.h | bin
	g++ -std=c++11 -o $@ $< -lboost_iostreams -lboost_system -lboost_filesystem

bin:
	mkdir -p $@

doc:
	mkdir doc
	cd doc && doxygen ..

clean:
	$(RM) -r bin doc
