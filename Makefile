all: mindcurrent

mindcurrent: main.cpp CellSyn.cpp CellSyn.h currents.cpp currents.h io.cpp io.h network.h network.cpp
	g++ -g -O3 -lm -Wall -fopenmp  main.cpp CellSyn.cpp io.cpp currents.cpp network.cpp -o mindcurrent
	 #icpc -O3 -Wall -fopenmp main.cpp CellSyn.cpp io.cpp currents.cpp network.cpp -o mindcurrent

#sanity check whether current version of code changes output (compared to previously stored test/.files)
check: mindcurrent
	./mindcurrent test/params.txt test connection_info2
	cd test; for f in *; do diff -u $$f .$$f; done

check-prepare: mindcurrent
	./mindcurrent test/params.txt test connection_info2
	cd test; for f in *; do cp -f  $$f .$$f; done

run: mindcurrent
	./mindcurrent params.txt out_w17_printIN connection_info2 sixclusters_nml.csv c1half_rand.csv undercut_nml_c1all.csv

runPaperDelays: mindcurrent
	./mindcurrent params.txt out /bazhlab/edelanois/tch/spindleNetwork/conn_s7_bqr_30-Aug-2017/connectionsFiles/connection_info2-delays

runPaperNoDelays: mindcurrent
	./mindcurrent params.txt out /bazhlab/edelanois/tch/spindleNetwork/conn_s7_bqr_30-Aug-2017/connectionsFiles/connection_info2-noDelays

doxy:
	doxygen ./docs/Doxyfile

clean: 
	-rm mindcurrent 

# constructs random network based on network.cfg file parameters
network:
	g++ -O2 generate_network.cpp -o generate_network
	./generate_network > connection_info2

# don't use this one
network_mri:
	g++ -O2 -g generate_network.cpp -o generate_network
	./generate_network /space/seh8/1/halgdev/projects/bqrosen/SS_SIM/out/conn/conn_s7_bqr_30-Aug-2017/DistanceMatrix_30-Aug-2017_LH.txt /space/seh8/1/halgdev/projects/bqrosen/SS_SIM/out/conn/conn_s7_bqr_30-Aug-2017/Map_642_To_10242_30-Aug-2017_LH.txt > connection_info2

# dont use
network_mri_2hems:
	g++ -O2 -g generate_network.cpp -o generate_network
	./generate_network /space/seh8/1/halgdev/projects/bqrosen/SS_SIM/out/conn/conn_s7_bqr_30-Aug-2017/DistanceMatrix_30-Aug-2017_LH.txt /space/seh8/1/halgdev/projects/bqrosen/SS_SIM/out/conn/conn_s7_bqr_30-Aug-2017/Map_642_To_10242_30-Aug-2017_LH.txt /space/seh8/1/halgdev/projects/bqrosen/SS_SIM/out/conn/conn_s7_bqr_30-Aug-2017/DistanceMatrix_30-Aug-2017_RH.txt /space/seh8/1/halgdev/projects/bqrosen/SS_SIM/out/conn/conn_s7_bqr_30-Aug-2017/Map_642_To_10242_30-Aug-2017_RH.txt /space/seh8/1/halgdev/projects/bqrosen/SS_SIM/out/conn/conn_s7_bqr_30-Aug-2017/Map_LH_RH_10242_30-Aug-2017.txt /space/seh8/1/halgdev/projects/bqrosen/SS_SIM/out/conn/conn_s7_bqr_30-Aug-2017/Map_LH_RH_642_30-Aug-2017.txt > connection_info2

# USE THIS ONE #
network_mri_2hems_local:
	g++ -O2 -g generate_network.cpp -o generate_network
	./generate_network /bazhlab/edelanois/tch/spindleNetwork/conn_s7_bqr_30-Aug-2017/DistanceMatrix_30-Aug-2017_LH.txt /bazhlab/edelanois/tch/spindleNetwork/conn_s7_bqr_30-Aug-2017/Map_642_To_10242_30-Aug-2017_LH.txt /bazhlab/edelanois/tch/spindleNetwork/conn_s7_bqr_30-Aug-2017/DistanceMatrix_30-Aug-2017_RH.txt /bazhlab/edelanois/tch/spindleNetwork/conn_s7_bqr_30-Aug-2017/Map_642_To_10242_30-Aug-2017_RH.txt /bazhlab/edelanois/tch/spindleNetwork/conn_s7_bqr_30-Aug-2017/Map_LH_RH_10242_30-Aug-2017.txt /bazhlab/edelanois/tch/spindleNetwork/conn_s7_bqr_30-Aug-2017/Map_LH_RH_642_30-Aug-2017.txt > connection_info2

# New multilayer connectivity, only LH
network_mri_LH_multilayer:
	g++ -O2 -g generate_network.cpp -o generate_network
	./generate_network /bazhlab/ysokolov/meg/newConn/conn_data/v7/dist_inMeters_prob_LH_onesOnDiag.txt /bazhlab/edelanois/tch/spindleNetwork/conn_s7_bqr_30-Aug-2017/Map_642_To_10242_30-Aug-2017_LH.txt /bazhlab/ysokolov/meg/newConn/conn_data/weights_LH.txt /bazhlab/edelanois/tch/spindleNetwork/conn_s7_bqr_30-Aug-2017/Map_642_To_10242_30-Aug-2017_RH.txt /bazhlab/edelanois/tch/spindleNetwork/conn_s7_bqr_30-Aug-2017/Map_LH_RH_10242_30-Aug-2017.txt /bazhlab/edelanois/tch/spindleNetwork/conn_s7_bqr_30-Aug-2017/Map_LH_RH_642_30-Aug-2017.txt /bazhlab/ysokolov/meg/newConn/conn_data/weights_INPY_both_hemi.txt /bazhlab/ysokolov/meg/newConn/conn_data/thalCort_dist_inMeters_LH.txt /bazhlab/ysokolov/meg/newConn/conn_data/dist_TC_LH_zeroDiag.txt > connection_info2


# New multilayer connectivity, both hems
network_mri_2hems_multilayer:
	g++ -O2 -g generate_network.cpp -o generate_network
	./generate_network /bazhlab/ysokolov/meg/newConn/conn_data/v7/dist_inMeters_prob_onesOnDiag.txt /bazhlab/edelanois/tch/spindleNetwork/conn_s7_bqr_30-Aug-2017/Map_642_To_10242_30-Aug-2017_LH.txt /bazhlab/ysokolov/meg/newConn/conn_data/weight_factor.data /bazhlab/edelanois/tch/spindleNetwork/conn_s7_bqr_30-Aug-2017/Map_642_To_10242_30-Aug-2017_RH.txt /bazhlab/edelanois/tch/spindleNetwork/conn_s7_bqr_30-Aug-2017/Map_LH_RH_10242_30-Aug-2017.txt /bazhlab/edelanois/tch/spindleNetwork/conn_s7_bqr_30-Aug-2017/Map_LH_RH_642_30-Aug-2017.txt /bazhlab/ysokolov/meg/newConn/conn_data/weights_INPY_both_hemi.txt > connection_info2

network_mri_2hems_one:
	g++ -O2 -g generate_network.cpp -o generate_network
	./generate_network /bazhlab/edelanois/tch/spindleNetwork/conn_rsfMRI_26-Mar-2019/DistanceMatrix_26-Mar-2019_LH.txt /bazhlab/edelanois/tch/spindleNetwork/conn_rsfMRI_26-Mar-2019/Map_642_To_10242_26-Mar-2019_LH.txt /bazhlab/edelanois/tch/spindleNetwork/conn_rsfMRI_26-Mar-2019/DistanceMatrix_26-Mar-2019_RH.txt /bazhlab/edelanois/tch/spindleNetwork/conn_rsfMRI_26-Mar-2019/Map_642_To_10242_26-Mar-2019_RH.txt /bazhlab/edelanois/tch/spindleNetwork/conn_rsfMRI_26-Mar-2019/Map_LH_RH_10242_14-Jul-2017.txt /bazhlab/edelanois/tch/spindleNetwork/conn_rsfMRI_26-Mar-2019/Map_LH_RH_642_14-Jul-2017.txt > connection_info2

network_mri_2hems_delay:
	g++ -O2 -g generate_network.cpp -o generate_network
	./generate_network /bazhlab/edelanois/tch/spindleNetwork/conn_rsfMRI_26-Mar-2019/DistanceMatrix_mm_14-Jul-2017_LH.txt /bazhlab/edelanois/tch/spindleNetwork/conn_rsfMRI_26-Mar-2019/Map_642_To_10242_26-Mar-2019_LH.txt /bazhlab/edelanois/tch/spindleNetwork/conn_rsfMRI_26-Mar-2019/DistanceMatrix_mm_14-Jul-2017_RH.txt /bazhlab/edelanois/tch/spindleNetwork/conn_rsfMRI_26-Mar-2019/Map_642_To_10242_26-Mar-2019_RH.txt /bazhlab/edelanois/tch/spindleNetwork/conn_rsfMRI_26-Mar-2019/Map_LH_RH_10242_14-Jul-2017.txt /bazhlab/edelanois/tch/spindleNetwork/conn_rsfMRI_26-Mar-2019/Map_LH_RH_642_14-Jul-2017.txt > connection_info2

