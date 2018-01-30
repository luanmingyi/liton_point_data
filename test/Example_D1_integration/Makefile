CXX = g++
proname = Example_D1_integration

%.o_debug:%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $(@) $(<)

out: CXXFLAGS = -std=c++11 -O2
out: CPPFLAGS = -DNDEBUG
out: $(proname).o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(@) $(<)

.PHONY: run
run: out
	./out
	cat example.txt

out_debug: CXXFLAGS = -std=c++11
out_debug: CPPFLAGS = -D_DEBUG
out_debug: $(proname).o_debug
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(@) $(<)

.PHONY: run_debug
run_debug: out_debug
	./out_debug
	cat *.txt

.PHONY: clean
clean:
	-rm *.o
	-rm *.o_debug
	-rm out
	-rm out_debug

.PHONY: clean_run
clean_run:
	-rm *.txt

.PHONY: clean_all
clean_all: clean clean_run
