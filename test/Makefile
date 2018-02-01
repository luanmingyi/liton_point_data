subfolder = Test_D0_PointData_C Test_D1_SizeT Test_D1_RangeT Test_D1_PointData_C Test_D1_PointData_H Test_D1_tecfile Example_D1_integration Example_D1_heat Test_D2_SizeT Test_D2_RangeT Test_D2_PointData_C Test_D2_PointData_H Test_D2_tecfile

out_flag = $(foreach n,$(subfolder),$(n).out_flag)
%.out_flag:
	echo $(@) > $(@)
	cd $(basename $(@)) && $(MAKE) out
out:$(out_flag)
	rm *.out_flag

run_flag = $(foreach n,$(subfolder),$(n).run_flag)
%.run_flag:
	echo $(@) > $(@)
	cd $(basename $(@)) && $(MAKE) run
run:$(run_flag)
	rm *.run_flag

out_debug_flag = $(foreach n,$(subfolder),$(n).out_debug_flag)
%.out_debug_flag:
	echo $(@) > $(@)
	cd $(basename $(@)) && $(MAKE) out_debug
out_debug:$(out_debug_flag)
	rm *.out_debug_flag

run_debug_flag = $(foreach n,$(subfolder),$(n).run_debug_flag)
%.run_debug_flag:
	echo $(@) > $(@)
	cd $(basename $(@)) && $(MAKE) run_debug
run_debug:$(run_debug_flag)
	rm *.run_debug_flag

clean_flag = $(foreach n,$(subfolder),$(n).clean_flag)
%.clean_flag:
	echo $(@) > $(@)
	cd $(basename $(@)) && $(MAKE) clean
clean:$(clean_flag)
	rm *.clean_flag

clean_run_flag = $(foreach n,$(subfolder),$(n).clean_run_flag)
%.clean_run_flag:
	echo $(@) > $(@)
	cd $(basename $(@)) && $(MAKE) clean_run
clean_run:$(clean_run_flag)
	rm *.clean_run_flag

clean_all_flag = $(foreach n,$(subfolder),$(n).clean_all_flag)
%.clean_all_flag:
	echo $(@) > $(@)
	cd $(basename $(@)) && $(MAKE) clean_all
clean_all:$(clean_all_flag)
	rm *.clean_all_flag
