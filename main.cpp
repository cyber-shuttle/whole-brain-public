///Created by M.Bazhenov, 1997-2009
//All rights reserved
//Most recently edited by Peter Lonjers(just search the net if you need to find me)
//I added all the MPI stuff but don't know much about the bio
#include <omp.h>
#include "CellSyn.h"    //header for all classes describing currents and basic cell types
#include "io.h"
#include "network.h"

int num_mp_threads; //number of openMP threads 

int print_c_sten; // if 1 prints all the connection strengths for cx every once in awile note this only works with openmp
int fre_print_cs; // how often to print connection strengths(every how many milliseconds)
int total_rec; // used for debugging number of non normal messages received in a time step
int total_rec_norm; // used for debugging number of normal(basic spiking messages) messages received in a time step

string output_location; // keeps track of what folder to put outputs in

Homeostasis homeo;

//boolean for if we are doing the very first step and send functions
int first_round = 1;
//main track of milliseconds of simulation so far. ii tracks steps of simulation locally
double t = 0;

//hou long we run the simulation
double tmax;
//time at which we print information from all cells
double t3D; 
//smaller time interval to print information from some cells
double ttime;
//how long main process runs for set in input
double run_time = 0; 

//TODO encapsulare these into single "network" object
// number of types of cells
int Num_Types; 
//gives the dimensions of the arrays of cells of each type
Pair *cell_sizes;
//total number of cells
int num_cells;
//total number of cells of each type
int *cell_numbers;
//cells info is an array of all cells containing the color of each cell and the numbers of different types of connections it has.
Cell_Info ****cells_info = NULL;
//array of all the pointers to cells belonging to this process
CellSyn **cells =NULL;

int *mapns, *hhns;

int mapi,hhi;



int awake_end = -1;
int stage2_end = -1; // 200000;
int stage3_end = 200000; //stage2_end; // -1;

double fac_gkl; 
double fac_gkl_TC;  
double fac_gkl_RE;  
      
double fac_gkv_cx;  
double fac_gkca_cx; 
double fac_gkm_cx;
      
double fac_gh_TC; 
      
double fac_AMPA_D2; 
double fac_AMPA_TC;
double fac_GABA_D2;
double fac_GABA_TC;

double fac_pL_cx_map;
double fac_pD_cx_map;


//run for openMP
//all == whether both RK and MAP should be called at the same time
//the order is important in order to ensure signaling and receiveing spikes between different steps

void runMP(int all){
    int iter,i;
#pragma omp parallel  private(i,iter)  num_threads(num_mp_threads)
  {
    // alternative schedular dynamic,group_sizes[my_color]/(num_mp_threads*2)

	// Run step for maps
	if (all){
#pragma omp for schedule(guided)  
	  for(i=0; i<mapi; i++){ 
		cells[mapns[i]]->step();
	  }
	}
	
	for(iter=0; iter<4; iter++){ 
#pragma omp for schedule(guided)    
	  for(i=0; i<hhi; i++){ 
		cells[hhns[i]]->step();
	  }
	}
    
#pragma omp for schedule(guided)   
	for(i=0; i<hhi; i++){ 
	  cells[hhns[i]]->reset_y();
	}

    // Signal and/or reset spikes here
    if (all) {
#pragma omp for schedule(guided)   
        for(i=0; i<mapi; i++) { 
            //Reset flip for maps which should be consumed in this time step (used for HH->Map connection)
            cells[mapns[i]]->reset_flip_maps();
        }

#pragma omp for schedule(guided)   
        for(i=0; i<hhi; i++) {
            cells[hhns[i]]->reset_flip_maps();
        }
	}

#pragma omp for schedule(guided)    
    for(i=0; i<mapi; i++) { 
        cells[mapns[i]]->signal_spike();
    }

#pragma omp for schedule(guided)    
    for(i=0; i<hhi; i++) {
        cells[hhns[i]]->signal_spike();
	}


  }
}

//gets index for a particular cell in a group(cells on the same process)
//ATM can be used as mapping from layer&2D location into global cells array
int get_cell_index(int type, int m, int n){
  // printf("proc_index- %d \n",cells_info[type][m][n]->proc_index);
  return cells_info[type][m][n]->cell_index;
}


//all arguments are inputs
//this functions decides what cells a process should simulate and starts them
void initialize_cells(CellSyn** &cells){
 //mapping between cell number and location, we might want to delete the whole concept later
 //loop rewritten from metis coloring routines
 int total=0;
 for(int x=0; x<Num_Types; x++)
  for(int y=0; y<cell_sizes[x].x; y++)
   for(int z=0; z<cell_sizes[x].y; z++)
		cells_info[x][y][z]->cell_index = total++;

 int type = 0;
 int m = 0;
 int n = 0;
 int total_set = 0;
  
 cells = new CellSyn*[num_cells];
 for(m = 0; m<num_cells; m++){
  cells[m]= NULL;
 }

 for(type = 0; type<Num_Types; type++){ 
	printf("Initializing cell %d for %d, %d\n",type, cell_sizes[type].x, cell_sizes[type].y);  
  for(m=0; m<cell_sizes[type].x; m++){
   for(n=0; n<cell_sizes[type].y; n++){
		cells[get_cell_index(type,m,n)]=CellSyn::initialize(type,m,n,cells_info,output_location);
	  total_set++;
   }
  }
 }

 //sanity checks
 if(total_set != num_cells){
  printf("process not given enough cells. This should not happen\n");
  printf("number given: %d\n",total_set);
  printf("number needed: %d\n",num_cells);
  exit(1);
 }

 for(m=0;m<num_cells;m++){
  if(cells[m]==NULL){
   printf("graph coloring probably failed\n");
   exit(1);
  }
 }
 //all checks passed
 return;
}


  // Scale connection strength by number of inputs
void scale_synapses(CellSyn** &cells, int num_cells){

    printf("Scaling Synapses ... \n");

    // Stupid algo to figure out size of each neuron's every type of input and scale by that number
    for(int i=0; i<num_cells; i++) {

        // // -------------- Scaling Method 1 ----------------------
        // // Uncomment the following for loop for connection specific scaling for all neurons in each group - Used before
        // for(int k=0; k < cells[i]->num_syns; k++) {
        //     if (cells[i]->syns[k]->type==E_GAP) {
        //         if ((max_connect_gap[cells[i]->type][cells[i]->syns[k]->from_type])==0) { 
        //             printf("Error during scaling GAP synapse for %d from %d =0 \n",cells[i]->type, cells[i]->syns[k]->from_type);
        //             exit(1);
        //         }
        //         cells[i]->syns[k]->strength= cells[i]->syns[k]->strength/max_connect_gap[cells[i]->type][cells[i]->syns[k]->from_type];
        //         cells[i]->syns[k]->mini_s = cells[i]->syns[k]->mini_s/max_connect_gap[cells[i]->type][cells[i]->syns[k]->from_type];
        //     }
        //     else{
        //         if ((max_connect[cells[i]->type][cells[i]->syns[k]->from_type])==0) {
        //             printf("Error during scaling synapse for %d from %d =0 \n",cells[i]->type, cells[i]->syns[k]->from_type);
        //             exit(1);
        //         }
        //         cells[i]->syns[k]->strength= cells[i]->syns[k]->strength/max_connect[cells[i]->type][cells[i]->syns[k]->from_type];
        //         cells[i]->syns[k]->mini_s = cells[i]->syns[k]->mini_s/max_connect[cells[i]->type][cells[i]->syns[k]->from_type];          
        //     }
        // }

        // -------------- Scaling Method 2 ----------------------
        // Uncomment the following for loop for neuron + connection specific scaling 
        int **syn_numbers=new int*[Num_Types];

        for(int k=0;k<Num_Types;k++)
            syn_numbers[k]= new int[ENUM_Syn_Types];

        for(int k=0;k<Num_Types;k++)
            for(int j=0;j<ENUM_Syn_Types;j++)
                syn_numbers[k][j]=0;

        for(int k=0; k < cells[i]->num_syns; k++){
            syn_numbers[cells[i]->syns[k]->from_type][cells[i]->syns[k]->type]=syn_numbers[cells[i]->syns[k]->from_type][cells[i]->syns[k]->type]+1;
        }
    
        for(int k=0; k < cells[i]->num_syns; k++){
            if (syn_numbers[cells[i]->syns[k]->from_type][cells[i]->syns[k]->type]>0) {
             // if (cells[i]->syns[k]->from_type == E_CX && cells[i]->syns[k]->type == E_AMPAMap_D1){
             //  printf("E_CX->E_CX ->%d Syns=%d", i, syn_numbers[cells[i]->syns[k]->from_type][cells[i]->syns[k]->type]);}
                cells[i]->syns[k]->strength= cells[i]->syns[k]->strength/syn_numbers[cells[i]->syns[k]->from_type][cells[i]->syns[k]->type];
                cells[i]->syns[k]->mini_s = cells[i]->syns[k]->mini_s/syn_numbers[cells[i]->syns[k]->from_type][cells[i]->syns[k]->type];
            }
        }

        for(int k=0;k<Num_Types;k++)
            delete[] syn_numbers[k];

        // -------------- Scaling Method 3 ----------------------
        // To be implemented later ...
        // // // Uncomment the following for loop for neuron + connection specific scaling -- ignores from-type
        // for(k=0; k < cells[i]->num_syns; ++k){
        //     if (syn_numbers[cells[i]->syns[k]->type]>0) {
        //         cells[i]->syns[k]->strength= cells[i]->syns[k]->strength/syn_numbers[cells[i]->syns[k]->type];
        //         cells[i]->syns[k]->mini_s = cells[i]->syns[k]->mini_s/syn_numbers[cells[i]->syns[k]->type];
        //     }
        // }
  
    }
  return;
}



void send_IStim1(int i, int j,double cx_stim, double cxa_stim){

  int group_index_CXa = get_cell_index(E_CXa,i,j);
  int group_index_CX = get_cell_index(E_CX,i,j);
  
  ((CX*)cells[group_index_CX]->base_cell)->I_Stim1 = cx_stim;
  ((CX*)cells[group_index_CXa]->base_cell)->I_Stim1 = cxa_stim;
}

//receive for getting things to the root printing
double print_receive(int m , int n, enum Cell_Type type){

  double message = 0;
  int group_index = get_cell_index(type,m,n);

  message = cells[group_index]->base_cell->get_v_soma();
  return message;
}

//receive for getting things to the root printing
double print_receive_var(int m , int n, enum Cell_Type type, int index){

  double message = 0;
  int group_index = get_cell_index(type,m,n);

  if (index==1000)
    message=cells[group_index]->e_current;
  else if (index==1001)
    message=cells[group_index]->i_current;
  else
    message = cells[group_index]->y_var;

  return message;
}

//receive for getting how many times a cell has spiked during the fre_window
int fre_receive(int m , int n, enum Cell_Type type){

  int message = 0;
  int group_index = get_cell_index(type,m,n);

  message = cells[group_index]->num_spikes;
  return message;
}

//receives dedritic voltage from connecting cell and returns it
//all args inputs
//m n and type are the index of the sending cell, my_type is the index of the receiving cell
double receive_dend(int my_type,int m, int n, enum Cell_Type type){

  int group_index = get_cell_index(type,m,n);

  return cells[group_index]->base_cell->get_v_dend();  //cell is in group so no need to receive a message to get its value
}

//receives voltage from connecting cell and returns it
//all args inputs
//m n and type are the index of the sending cell, my_type is the index of the receiving cell
int receive_spike(int m, int n, enum Cell_Type type){

  int group_index = get_cell_index(type,m,n);

  return cells[group_index]->get_flip(); //cell is in group so no need to receive a message to get its value

}

int receive_spike_for_maps(int m, int n, enum Cell_Type type){

  int group_index = get_cell_index(type,m,n);

  return cells[group_index]->get_flip_for_maps(); //cell is in group so no need to receive a message to get its value

}


//receive for getting things to the root printing
void print_receive_gsynapse(int m , int n, enum Cell_Type type, enum Syn_Type stype, FILE *fp){

 int group_index = get_cell_index(type,m,n);

 // Run through all synapses
 for(int i =0; i < cells[group_index]->num_syns; i++){
  if (cells[group_index]->syns[i]->type == stype){

   // fprintf(fp,"%d %lf ", m, ((AMPA_D3*)(cells[group_index]->syns[i]))->g_AMPA0); //D3
   // fprintf(fp,"%d %d %d %lf ", m,n,i, ((AMPAmapD1*)(cells[group_index]->syns[i]))->d); 
   fprintf(fp,"%d %d %d %lf ", m,n,i, ((GABAAmapD1*)(cells[group_index]->syns[i]))->I); 
  }
 }
 //message = cells[group_index]->base_cell->get_v_soma();
}


//receive for getting things to the root printing
void print_receive_gsynapse_index(int m , int n, enum Cell_Type from_type, enum Cell_Type type, enum Syn_Type stype, FILE *fp){

  int group_index = get_cell_index(type,m,n);

  // Run through all synapses
  for(int i =0; i < cells[group_index]->num_syns; i++){
    if ((cells[group_index]->syns[i]->type == stype) && ((cells[group_index]->syns[i]->from_type == from_type))){
      int from=cells[cells[group_index]->syns[i]->from_cell]->m;
      // print_receive_gsynapse_index((i, j, E_CX, E_CX, E_AMPA_D3, fp);

      fprintf(fp,"%d %d ", from ,m);
      // fprintf(fp," %lf ", (cells[group_index]->syns[i])->strength);
    }
  }

  //message = cells[group_index]->base_cell->get_v_soma();
}



//receive for getting things to the root printing
void print_receive_gsynapse_multLayer(int m , int n, enum Cell_Type from_type, enum Cell_Type type, enum Syn_Type stype, FILE *fp){

 int group_index = get_cell_index(type,m,n);

 // Run through all synapses
 for(int i =0; i < cells[group_index]->num_syns; i++){
  if ((cells[group_index]->syns[i]->type == stype) && (cells[group_index]->syns[i]->from_type == from_type)){

   // fprintf(fp,"%d %lf ", m, ((AMPA_D3*)(cells[group_index]->syns[i]))->g_AMPA0); //D3
   fprintf(fp,"%lf ", ((AMPAmapD1*)(cells[group_index]->syns[i]))->g); 
   // fprintf(fp,"%d %d %d %lf ", m,n,i, ((GABAAmapD1*)(cells[group_index]->syns[i]))->I); 
  }
 }
 //message = cells[group_index]->base_cell->get_v_soma();
}




extern FILE *f28; //TODO get rid of this

void load_input_data(int argc, char *argv[]){ //TODO move global variables into local ones perhaps move it to io.h/c
    //checks inputs
  if (argc < 4) {
    puts("Bad command: Should be");
    puts("input_file_name output_directory connections_file_name"); 
    exit(1);
  }
  //TODO check to make sure output location exists
  output_location = argv[2];
  output_location = output_location + "/";

  //open up the file of connections
  FILE *connections_file;
  if (!(connections_file=fopen(argv[3],"r"))) {
    printf("%s file for connections does not exist or something\n",argv[3]);
    exit(1); 
  }

  set_cell_sizes(connections_file,Num_Types,cell_sizes,cell_numbers); //scans connections file for how many cells of each type we have
  create_connections_from_file(&cells_info,cell_sizes,connections_file,Num_Types);
    
  //finds total number of cells
  num_cells = 0;
  for(int iter = 0; iter<Num_Types; iter++){
    num_cells = num_cells+cell_numbers[iter];
  }
  printf("number cells: %d\n",num_cells);

 } //load_input_data

//+++++++++++++++++++ MAIN PROGRAM +++++++++++++++++++++++++++++++++++++++++++
int main(int argc,char **argv){

  LocalFieldPotential LFP;

  load_input_data(argc,argv);
  load_input_params(argc,argv,tmax,t3D,ttime,num_mp_threads,print_c_sten,fre_print_cs,LFP.local_field_effect,LFP.lfp_scale,LFP.num_field_layers,homeo.boost,homeo.amp_boost,homeo.con_boost,homeo.fre_boost,homeo.target_f,homeo.fre_window,homeo.num_regions);

  // seed random number generator with current second of execution
  // srand(time(NULL)); 
  srand(2002); 

  homeo.allocate();

  printf("Call Initialize cells ...\n");

  initialize_cells(cells);

  printf("Scaling synapses ...\n");

  scale_synapses(cells, num_cells);

  LFP.init();
  open_files(output_location,LFP.field_file,LFP.num_field_layers);  
  LFP.allocate_state_save(cell_sizes);

  //------Main Loop--------
  printf("timer started");
  // time_t start_time = time(NULL);
  printf("\n starting main loop: time= %lf: end_time= %lf\n", t,tmax);
  // time_t current_time;
  int print_count = 0;
  int ii = 0;
  int i=0;

  mapi=0; hhi=0;
  mapns=new int[num_cells];
  hhns=new int[num_cells];
  for(i=0; i<num_cells; i++){
	if (cells[i]->ismap){
	  mapns[mapi]=i;
	  mapi=mapi+1;
	}
	else{
	  hhns[hhi]=i;
	  hhi=hhi+1;
	} 
  }

  print_used_mem();
  time_t start_time = time(NULL);
  printf("root set up tau:%lf time:%lf tmax:%lf\n",TAU,t,tmax);



  double s2_scale=1.2;
  double s3_scale=1.2; // 2.0;

  double gkl_awake_fix     = 0.8333; // 0.19; //
  double gkl_s2            = gkl_awake_fix*s2_scale; //
  double gkl_s3            = gkl_awake_fix*s3_scale; //

  double gkl_TC_awake_fix  =  0.8; //correct1_func(awake_ach_fac ,0.0);
  double gkl_TC_s2         = gkl_TC_awake_fix*s2_scale; //
  double gkl_TC_s3         = gkl_TC_awake_fix*s3_scale; //

  double gkl_RE_awake_fix  =  0.9; //correct1_func(awake_ach_fac ,0.0);
  double gkl_RE_s2         = gkl_RE_awake_fix*((2-s2_scale/2)-0.5); //
  double gkl_RE_s3         = gkl_RE_awake_fix*((2-s3_scale/2)-0.5); //

  double awake_AMPAd2_fix  = 0.8772; // 0.19; //wake_ach_fac;
  double s2_AMPAd2         =awake_AMPAd2_fix*s2_scale*0.95; //
  double s3_AMPAd2         = 0.365; //0.1; // 0.01;// 0.001; //0.5; // 1.14; // awake_AMPAd2_fix*s3_scale*1.28; //

  double gh_TC_awake       =-8.0; //
  double gh_TC_s2          = 0.0; // -2.0; 
  double gh_TC_s3          =  -1.0; 

  double awake_GABAd2_fix  = 0.8333; // 0.22; //0.6; //awake_gaba_fac;
  double s2_GABAd2         =awake_GABAd2_fix*s2_scale; //
  double s3_GABAd2         =awake_GABAd2_fix*s3_scale; //

  double awake_GABA_TC_fix =  0.6; //awake_gaba_fac;
  double s2_GABA_TC        =awake_GABA_TC_fix*s2_scale; //
  double s3_GABA_TC        =awake_GABA_TC_fix*s3_scale; //

  double awake_AMPA_TC     = 1.0; //0.5; 
  double s2_AMPA_TC        = 1.0; //0.5; 
  double s3_AMPA_TC        = 0.1; 

  double gk_cx_slow_awake  = 1.0; 
  double gk_cx_slow_s2     = 1.0; 
  double gk_cx_slow_s3     = 1.0; 

  double gk_cx_spike_awake = 1.0; 
  double gk_cx_spike_s2    = 1.0; 
  double gk_cx_spike_s3    = 1.0; 

  double pL_cx_map_awake   = 1.0;
  double pL_cx_map_s2      = 0.15;
  double pL_cx_map_s3      = 0.5;

  double pD_cx_map_awake   = 1.0;
  double pD_cx_map_s2      = 5.5;
  double pD_cx_map_s3      = 4.0;



  print_connectivity(cell_sizes);





  while( t < tmax){ 
    //printf("total_rec:%d root\n",total_rec);
    total_rec = 0;


  if ((t<=awake_end)){

      // Fix all values for awake state
      fac_AMPA_D2 = awake_AMPAd2_fix*0.7;
      fac_AMPA_TC = awake_AMPA_TC;
      fac_GABA_D2 = awake_GABAd2_fix;
      fac_GABA_TC = awake_GABA_TC_fix;
      fac_gkl_RE  = gkl_RE_awake_fix;
      fac_gkl_TC  = gkl_TC_awake_fix*0.5;
      fac_gkl     = gkl_awake_fix*0.7;
      fac_gh_TC   = gh_TC_awake*3;
      fac_gkca_cx = gk_cx_slow_awake;
      fac_gkm_cx  = gk_cx_slow_awake;
      fac_gkv_cx  = gk_cx_spike_awake;

      fac_pL_cx_map = pL_cx_map_awake;
      fac_pD_cx_map = pD_cx_map_awake;

    }else if ((t>awake_end&&t<=stage2_end)){

      // Fix all values for S2
      fac_AMPA_D2 = s2_AMPAd2;
      fac_AMPA_TC = s2_AMPA_TC;
      fac_GABA_D2 = s2_GABAd2;
      fac_GABA_TC = s2_GABA_TC;
      fac_gkl_RE  = gkl_RE_s2;
      fac_gkl_TC  = gkl_TC_s2;
      fac_gkl     = gkl_s2;
      fac_gh_TC   = gh_TC_s2;
      fac_gkca_cx = gk_cx_slow_s2;
      fac_gkm_cx  = gk_cx_slow_s2;
      fac_gkv_cx  = gk_cx_spike_s2;

      fac_pL_cx_map = pL_cx_map_s2;
      fac_pD_cx_map = pD_cx_map_s2;

    }else if ((t>stage2_end&&t<=stage3_end)){
      // Fix all values for S3
      fac_AMPA_D2 = s3_AMPAd2;
      fac_AMPA_TC = s3_AMPA_TC;
      fac_GABA_D2 = s3_GABAd2;
      fac_GABA_TC = s3_GABA_TC;
      fac_gkl_RE  = gkl_RE_s3;
      fac_gkl_TC  = gkl_TC_s3;
      fac_gkl     = gkl_s3;
      fac_gh_TC   = gh_TC_s3;
      fac_gkca_cx = gk_cx_slow_s3;
      fac_gkm_cx  = gk_cx_slow_s3;
      fac_gkv_cx  = gk_cx_spike_s3;

      fac_pL_cx_map = pL_cx_map_s3;
      fac_pD_cx_map = pD_cx_map_s3;
    }
  else { //awake after sleep
      fac_AMPA_D2 = awake_AMPAd2_fix*0.7;
      fac_AMPA_TC = awake_AMPA_TC;
      fac_GABA_D2 = awake_GABAd2_fix;
      fac_GABA_TC = awake_GABA_TC_fix;
      fac_gkl_RE  = gkl_RE_awake_fix;
      fac_gkl_TC  = gkl_TC_awake_fix*0.5;
      fac_gkl     = gkl_awake_fix*0.7;
      fac_gh_TC   = gh_TC_awake*3;
      fac_gkca_cx = gk_cx_slow_awake;
      fac_gkm_cx  = gk_cx_slow_awake;
      fac_gkv_cx  = gk_cx_spike_awake;

      fac_pL_cx_map = pL_cx_map_awake;
      fac_pD_cx_map = pD_cx_map_awake;
    }



    //printf("Running time %lf \n",t);

    if(t>print_count){
        // current_time = time(NULL);
        //double time_from_start = difftime(current_time,start_time);
        // printf("milliseconds of sim = %lf real time from timer start in seconds: %lf\n", t,time_from_start);
        printf("milliseconds of sim = %lf \n", t);
        print_count = print_count + 1;
    }

    ii = ii + 1; //number iterations
	
    if (((ii)/(TAUr))*(TAUr) == (ii))
        runMP(1);
    else
        runMP(0);	  

    //applies field effects
    if(LFP.local_field_effect == 1){
      LFP.apply_field(ii,t,ttime,t,cell_sizes); //careful there is another barrier in here
    }

    //prints all 1D cell Voltages or all x cells for one y in 2d voltages occasionally
    if((t > ttime) && ((ii/(50))*(50) == ii)){ //multiples of 50 print
      print_freq(LFP.cx_base_v_SOMA,LFP.cx5b_base_v_SOMA,cell_sizes,t,Num_Types);
    }

//    //prints all cell voltages occationally
//    if((t > t3D) && (((ii)/(500))*(500) == (ii))){ //multiples of 500 print
//      print_occ(cell_sizes);
//    }

    if(t > t3D){
            if ((t > 1500 ) && ( t < 3500 ) && ((ii/(50))*(50) == ii) ){
                    print_occ(cell_sizes);
            }
            else if     (((ii)/(500))*(500) == (ii)){ 
                    print_occ(cell_sizes);
            }
    }


    //prints strength of connections
    if((ii >= fre_print_cs) && (ii % fre_print_cs) == 0 && print_c_sten == 1){
      int i =0;
      //TODO we are able to enumerate E_CX directly
      for(i=0; i<num_cells; i++){
	if(cells[i]->type == E_CX){
	  ((CXsyn*)cells[i])->print_c_stren(f28);
	}
      }
    }

    //implements homeostatic mechanisms
    if((ii >= homeo.fre_window) && (ii % homeo.fre_window) == 0){
      homeo.spike_fre_calc(cell_sizes,cell_numbers); //figures out frequency of spiking
      int i = 0;
      for(i=0; i<num_cells; i++){
        if(cells[i]->type == E_CX || cells[i]->type == E_CXa || cells[i]->type == E_CX6 || cells[i]->type == E_CX3 || cells[i]->type == E_CX4 || cells[i]->type == E_CX5a || cells[i]->type == E_CX5b ){
          //for openmp just sets num spikes to 0 
	  cells[i]->num_spikes=0;
	}
      }

      homeo.boost_activity(cell_sizes,cells,num_cells); //causes homeostatic changes
    }

    t = t + TAU; //increase time
  } //end while loop

  int computationTime = (int)difftime(time(NULL),start_time);
  printf("Computation time: %d s\n",(int)difftime(time(NULL),start_time));
  FILE * computationTimeFile = fopen((output_location+"computationSeconds.txt").c_str(), "w");
  fprintf(computationTimeFile, "%d\n", computationTime);
  fclose(computationTimeFile);
  close_files(LFP.field_file,LFP.num_field_layers);

  return 0;
}

