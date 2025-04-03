
#include <stdio.h>
#include <string>
#include "CellSyn.h" //Pair

using namespace std;

//prints all cells voltages
void print_occ(Pair *cell_sizes);
void print_connectivity(Pair *cell_sizes);
void print_used_mem();
//prints one column/row of cells voltages in a weird format
void print_freq(double **cx_base_v_SOMA, double **cx5b_base_v_SOMA,Pair *cell_sizes, double const t, int Num_Types);
void close_files(FILE **field_file, int num_field_layers);
void open_files(string output_location,FILE **field_file, int num_field_layers);


//read input from either input file or commandline
//TODO in long term we should have just single parameter class instead of passing zillion args.
void load_input_params(
  int argc,char *argv[],

  double& tmax,
  double&  t3D,
  double&  ttime,
  int&  num_mp_threads,
  int&  print_c_sten,
  int&  fre_print_cs,
  int&  LFP_local_field_effect,
  double&  LFP_lfp_scale,
  int&  LFP_num_field_layers,
  int&  homeo_boost,
  double&  homeo_amp_boost,
  double&  homeo_con_boost,
  double&  homeo_fre_boost,
  double&  homeo_target_f,
  int&  homeo_fre_window,
  int&  homeo_start,
  int&  homeo_num_regions,
  int&  homeo_num_clusters,
  int&  undercut,
  int&  undercut_start,
  double& stim_on,
  double& stim_start,
  double& stim_stop,
  double& stim_stren

 );

