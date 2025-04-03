///Created by M.Bazhenov, 1997-2009
//All rights reserved
//Most recently edited by Peter Lonjers(just search the net if you need to find me)
//I added all the MPI stuff but don't know much about the bio
#include <omp.h>
#include <assert.h>
#include "CellSyn.h"    //header for all classes describing currents and basic cell types
#include "io.h"
#include "network.h"
#include <stdio.h>

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
// current stimulation parameters
double stim_on;
double stim_start;
double stim_stop;
double stim_stren;

//TODO encapsulare these into single "network" object
// number of types of cells
int Num_Types; 
int Num_Cells_CX = 10242;
//gives the dimensions of the arrays of cells of each type
Pair *cell_sizes;
//total number of cells
int num_cells;
//total number of cells of each type
int *cell_numbers;
int *clust_ids;
int *cut_ids;
int *stim_ids;
//cells info is an array of all cells containing the color of each cell and the numbers of different types of connections it has.
Cell_Info ****cells_info = NULL;
//array of all the pointers to cells belonging to this process
CellSyn **cells =NULL;

int *mapns, *hhns;

int mapi,hhi;

double ***Dist_Prob;
const char *dist_file = "/bazhlab/ysokolov/meg/newConn/conn_data/v7/dist_inMeters_prob_LH_onesOnDiag.txt";

int awake_end = 300000; //-1;
int stage2_end = -1; // 200000;
int stage3_end = -1; //200000; //stage2_end; // -1;
int recovery_count = 0;
int ramp = 0;


int mriPointsFull = 20500;
int num_dist_param = 2;
int MRI_POINTS = 12500;

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

double fac_yrest;


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


// Load Distance between columns
void load_3D_dist_prob(const char *file){
  //cerr << "load_3D_dist_prob" << endl;

  FILE *f=fopen(file,"r"); assert(f);
  printf("Dist File Opened... \n");

  while(1){
    int from,to;
    double dist, conn_prob;
    if(fscanf(f,"%d %d %lf %lf\n",&from,&to,&dist,&conn_prob)!=4){
      break;
    }
    from--; to--; //MATLAB vs C indexing
    assert(from<MRI_POINTS); assert(to<MRI_POINTS);
    Dist_Prob[from][to][0]=dist;
    Dist_Prob[from][to][1]=conn_prob;
  }

  fclose(f);
}


//gets index for a particular cell in a group(cells on the same process)
//ATM can be used as mapping from layer&2D location into global cells array
int get_cell_index(int type, int m, int n){
  // printf("proc_index- %d \n",cells_info[type][m][n]->proc_index);
  return cells_info[type][m][n]->cell_index;
}

double get_connection_length(int from, int to){
  return Dist_Prob[from][to][0];
}


void init_pL_scaler(){
  for(int i=0; i<num_cells; i++){
    cells[i]->base_cell->pL_scaler = 1;
  }
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

  // 250944
 printf("num cells: %d", num_cells);
 //exit(1);
  
 cells = new CellSyn*[num_cells];
 printf("New Cells!");
 for(m = 0; m<num_cells; m++){
  printf("cells %d\n", m);
  cells[m]= NULL;
 }

 printf("cells allocated");

 // Failing here with new additions for synapse length
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

        // SCALE MINIS AGAIN
        // for(int k=0; k < cells[i]->num_syns; k++){
        //     cells[i]->syns[k]->mini_s = cells[i]->syns[k]->mini_s * 0.0;
        //     cells[i]->syns[k]->mini_fre = cells[i]->syns[k]->mini_fre * 0.0;
        // }

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

vector<double> print_receive_allIN(int m , int n, enum Cell_Type type){
  int group_index = get_cell_index(type,m,n);

 vector<double> message = {};

  if(m == 1060){
    message = cells[group_index]->base_cell->get_in_all();
  }

  return message;
}


double print_receive_syn(int m , int n, enum Cell_Type type){

  int group_index = get_cell_index(type,m,n); // ranges from 2k to 12k?

  double sum = 0;
  double avg = 0;
  int counter = 0;
  int k = 0;
  
  for(k=0; k < cells[group_index]->num_syns; ++k){
    if(cells[group_index]->syns[k]->type == E_AMPAMap_D1){   // || cells[i]->syns[k]->type == E_NMDAMap_D1 
     sum = sum + abs(cells[group_index]->syns[k]->strength);
     counter = counter + 1;
    }
  }
  if(sum==0){
    avg = sum;
  }
  else{
    avg = sum / counter; 
  }

  return avg;
}

double print_receive_g(int m , int n, enum Cell_Type type){

  int group_index = get_cell_index(type,m,n); // ranges from 2k to 12k?

  double sum_g = 0;
  double avg_g = 0;
  int counter = 0;
  int k = 0;
  
  // Average over synapses for a given cell
  for(k=0; k < cells[group_index]->num_syns; ++k){
    if(cells[group_index]->syns[k]->type == E_AMPAMap_D1){   // || cells[i]->syns[k]->type == E_NMDAMap_D1 
     sum_g = sum_g + cells[group_index]->syns[k]->g_track;
     counter = counter + 1;
    }
  }

  if(sum_g==0){
    avg_g = sum_g;
  }
  else{
    avg_g = sum_g / counter; 
  }

  return avg_g;
}

double print_receive_g_1cell1syn(int m , int n, enum Cell_Type type){

  int group_index = get_cell_index(type,m,n); // ranges from 2k to 12k?

  double sum_g = 0;
  double avg_g = 0;
  int counter = 0;
  int k = 0;

  for(k=0; k < cells[group_index]->num_syns; ++k){
    if(cells[group_index]->syns[k]->type == E_AMPAMap_D1){ 
      sum_g = sum_g + cells[group_index]->syns[k]->g_track;
      counter = counter + 1;
    }
  }

  if(sum_g==0){
    avg_g = sum_g;
  }
  else{
    avg_g = sum_g / counter; 
  }

  return avg_g;
}


double print_receive_d(int m , int n, enum Cell_Type type){

  int group_index = get_cell_index(type,m,n); // ransges from 2k to 12k?

  double sum_d = 0.0;
  double avg_d = 0.0;
  int counter = 0;
  int k = 0;
  
  // Average over synapses for a given cell
  for(k=0; k < cells[group_index]->num_syns; ++k){
    if(cells[group_index]->syns[k]->type == E_AMPAMap_D1){
      sum_d = sum_d + cells[group_index]->syns[k]->d_track;
      counter = counter +1;
    }
  }

  if(sum_d==0){
    avg_d = sum_d;
  }
  else{
    avg_d = sum_d / counter; 
  }

  return avg_d;
}


double print_receive_minis(int m , int n, enum Cell_Type type){

  int group_index = get_cell_index(type,m,n); // ransges from 2k to 12k?

  double sum_mini = 0;
  double avg_mini = 0;
  int counter = 0;
  int k = 0;
  
  // Average over synapses for a given cell
  for(k=0; k < cells[group_index]->num_syns; ++k){
    if(cells[group_index]->syns[k]->type == E_AMPAMap_D1){   // || cells[i]->syns[k]->type == E_NMDAMap_D1 
      sum_mini = sum_mini + abs(cells[group_index]->syns[k]->mini_s);
      counter = counter +1;
    }
  }

  if(sum_mini==0){
    avg_mini = sum_mini;
  }
  else{
    avg_mini = sum_mini / counter; //cells[group_index]->num_syns; 
  }

  return avg_mini;
}

double print_receive_minifre(int m , int n, enum Cell_Type type){

  int group_index = get_cell_index(type,m,n); // ransges from 2k to 12k?

  double sum_mini = 0;
  double avg_mini = 0;
  int k = 0;
  
  // Average over synapses for a given cell
  for(k=0; k < cells[group_index]->num_syns; ++k){
    if(cells[group_index]->syns[k]->from_type ==  cells[group_index]->syns[k]->to_type && (cells[group_index]->syns[k]->type == E_AMPAMap_D1)){   // || cells[i]->syns[k]->type == E_NMDAMap_D1 
      sum_mini = sum_mini + abs(cells[group_index]->syns[k]->mini_fre);
    }
  }

  if(sum_mini==0){
    avg_mini = sum_mini;
  }
  else{
    avg_mini = sum_mini / k; //cells[group_index]->num_syns; 
  }

  return avg_mini;
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
  if (argc < 6) {
    puts("Bad command: Should be");
    puts("input_file_name output_directory connections_file_name hsp_ids undercut_ids stim_ids"); 
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
  // above prints "cells are set"
    
  //finds total number of cells
  num_cells = 0;
  for(int iter = 0; iter<Num_Types; iter++){
    num_cells = num_cells+cell_numbers[iter];
  }
  printf("number cells: %d\n",num_cells);

//open up the local HSP ID's
  FILE *hsp_ids;
  if (!(hsp_ids=fopen(argv[4],"r"))) {
    printf("%s file for HSP region ids does not exist or something\n",argv[4]);
    exit(1); 
  }

  get_clusters(hsp_ids,Num_Cells_CX,clust_ids);

  //open up the local undercut ID's
  FILE *undercut_ids;
  if (!(undercut_ids=fopen(argv[5],"r"))) {
    printf("%s file for undercut ids does not exist or something\n",argv[5]);
    exit(1); 
  }

  get_undercut(undercut_ids,Num_Cells_CX,cut_ids);

   //open up the stimulation ID's
  FILE *stimulation_ids;
  if (!(stimulation_ids=fopen(argv[6],"r"))) {
    printf("%s file for stimulation ids does not exist or something\n",argv[6]);
    exit(1); 
  }

  get_stimulation(stimulation_ids,Num_Cells_CX,stim_ids);

 } //load_input_data

//+++++++++++++++++++ MAIN PROGRAM +++++++++++++++++++++++++++++++++++++++++++
int main(int argc,char **argv){

  LocalFieldPotential LFP;

  load_input_data(argc,argv);
  load_input_params(argc,argv,tmax,t3D,ttime,num_mp_threads,print_c_sten,fre_print_cs,LFP.local_field_effect,LFP.lfp_scale,LFP.num_field_layers,homeo.boost,homeo.amp_boost,homeo.con_boost,homeo.fre_boost,homeo.target_f,homeo.fre_window,homeo.start,homeo.num_regions, homeo.num_clusters, homeo.undercut, homeo.undercut_start, stim_on, stim_start, stim_stop, stim_stren);
  // seed random number generator with current second of execution
  // srand(time(NULL)); 
  srand(2002); 

  printf("input params loaded... \n");

  homeo.allocate();

  // Initialize Dist Prob
  Dist_Prob = new double**[mriPointsFull];
  for(int i = 0; i < mriPointsFull; i++){
    Dist_Prob[i] = new double*[mriPointsFull];  
    for(int j = 0; j < mriPointsFull; j++){
      Dist_Prob[i][j] = new double[num_dist_param];  
      for(int k = 0; k < num_dist_param; k++){
        Dist_Prob[i][j][k] = -1;
      }
    }
  }


  printf("Loading Dist Matrix... \n");
  load_3D_dist_prob(dist_file);

  printf("Call Initialize cells ...\n");

  initialize_cells(cells);

  init_pL_scaler();

  printf("Scaling synapses ...\n");

  scale_synapses(cells, num_cells);

  LFP.init();
  open_files(output_location,LFP.field_file,LFP.num_field_layers);  
  LFP.allocate_state_save(cell_sizes);

   FILE * SpikeFre = fopen((output_location+"spikefre").c_str(), "w");

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

  double gkl_awake_fix     = 1; //0.8333; // 0.19; //
  double gkl_s2            = gkl_awake_fix*s2_scale; //
  double gkl_s3            = 1; //gkl_awake_fix*s3_scale;

  double gkl_TC_awake_fix  = 0.8; //correct1_func(awake_ach_fac ,0.0);
  double gkl_TC_s2         = gkl_TC_awake_fix*s2_scale; //
  double gkl_TC_s3         = 0.96;

  double gkl_RE_awake_fix  =  0.9; //correct1_func(awake_ach_fac ,0.0);
  double gkl_RE_s2         = gkl_RE_awake_fix*((2-s2_scale/2)-0.5); //
  double gkl_RE_s3         = 0.81;

  double awake_AMPAd2_fix  = 0.35; //0.2; // 0.19; //wake_ach_fac;
  double s2_AMPAd2         =awake_AMPAd2_fix*s2_scale*0.95; //
  double s3_AMPAd2         = 0.365 + (0.365*0.5);

  double gh_TC_awake       =-1.0; //
  double gh_TC_s2          = 0.0; // -2.0; 
  double gh_TC_s3          =-1.0; 

  double awake_GABAd2_fix  = 0.5; // 0.22; //0.6; //awake_gaba_fac;
  double s2_GABAd2         =awake_GABAd2_fix*s2_scale; //
  double s3_GABAd2         =1; //

  double awake_GABA_TC_fix = 0.6; //awake_gaba_fac;
  double s2_GABA_TC        =awake_GABA_TC_fix*s2_scale; //
  double s3_GABA_TC        = 0.72; //

  double awake_AMPA_TC     = 0.0; //0.1;
  double s2_AMPA_TC        = 1.0; //0.5; 
  double s3_AMPA_TC        = 0.1; 

  double gk_cx_slow_awake  = 1.0; 
  double gk_cx_slow_s2     = 1.0; 
  double gk_cx_slow_s3     = 1.0; 

  double gk_cx_spike_awake = 1.0; 
  double gk_cx_spike_s2    = 1.0; 
  double gk_cx_spike_s3    = 1.0; 

  double pL_cx_map_awake   = 0.15;
  double pL_cx_map_s2      = 0.15;
  double pL_cx_map_s3      = 0.5;

  double pD_cx_map_awake   = 3.0;
  double pD_cx_map_s2      = 5.5;
  double pD_cx_map_s3      = 4.0;

  double yrest_awake       = 0.991;
  double yrest_s2          = 1;
  double yrest_s3          = 1;



  print_connectivity(cell_sizes);

  // DO ONCE to get orig HSP parameters
  homeo.orig_strengths(cells, num_cells);

  // Scale within-column connections
  homeo.reduce_column(cells,num_cells);


  while( t < tmax){ 
    //printf("total_rec:%d root\n",total_rec);
    total_rec = 0;

  //printf("pL scaler: %f", cells[100]->base_cell->pL_scaler);

  if ((t<=awake_end)){

      // Fix all values for awake state
      fac_AMPA_D2 = awake_AMPAd2_fix;
      fac_AMPA_TC = awake_AMPA_TC;
      fac_GABA_D2 = awake_GABAd2_fix;
      fac_GABA_TC = awake_GABA_TC_fix;
      fac_gkl_RE  = gkl_RE_awake_fix;
      fac_gkl_TC  = gkl_TC_awake_fix;
      fac_gkl     = gkl_awake_fix;
      fac_gh_TC   = gh_TC_awake;
      fac_gkca_cx = gk_cx_slow_awake;
      fac_gkm_cx  = gk_cx_slow_awake;
      fac_gkv_cx  = gk_cx_spike_awake;

      fac_pL_cx_map = pL_cx_map_awake;
      fac_pD_cx_map = pD_cx_map_awake;

      fac_yrest = yrest_awake;

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

      fac_yrest = yrest_s2;

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

      fac_yrest = yrest_s3;
    }
  else { //awake after sleep
      fac_AMPA_D2 = awake_AMPAd2_fix;
      fac_AMPA_TC = awake_AMPA_TC;
      fac_GABA_D2 = awake_GABAd2_fix;
      fac_GABA_TC = awake_GABA_TC_fix;
      fac_gkl_RE  = gkl_RE_awake_fix;
      fac_gkl_TC  = gkl_TC_awake_fix;
      fac_gkl     = gkl_awake_fix;
      fac_gh_TC   = gh_TC_awake;
      fac_gkca_cx = gk_cx_slow_awake;
      fac_gkm_cx  = gk_cx_slow_awake;
      fac_gkv_cx  = gk_cx_spike_awake;

      fac_pL_cx_map = pL_cx_map_awake;
      fac_pD_cx_map = pD_cx_map_awake;

      fac_yrest = yrest_awake;
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

    // implement undercut if its time, if i say so
    if((t >= homeo.undercut_start) && (homeo.undercut == 1)){
     //printf("do the cut");
     homeo.perform_undercut(cut_ids,cell_sizes,cells,num_cells,1);
     homeo.undercut = -1; // set to -1 so it only cuts once 
     printf("Cut \n");
    }

    // BM hand ramp connections
    // Note: using CLUST ID's will do all of c1
    // using CUT ID's will do traditionally cut cells
    // if((t >= 10000) && (ramp==1)){
    //   homeo.hand_ramp(clust_ids,cell_sizes,cells,num_cells);
    //   ramp = -1;
    // }

    // Double stim stren every 5 seconds
    // if((t>stim_start+1000) && ((ii % 250000) == 0)){
    //   stim_stren = stim_stren * 2;
    // }

    //implements homeostatic mechanisms
    if((ii >= homeo.fre_window) && (ii % homeo.fre_window) == 0){
      //homeo.spike_fre_calc(cell_sizes,cell_numbers); //figures out frequency of spiking
      homeo.spike_fre_calc_local(clust_ids,cell_sizes,cell_numbers, SpikeFre); 
      int i = 0;
      for(i=0; i<num_cells; i++){
        if(cells[i]->type == E_CX || cells[i]->type == E_CXa || cells[i]->type == E_CX6 || cells[i]->type == E_CX3 || cells[i]->type == E_CX4 || cells[i]->type == E_CX5a || cells[i]->type == E_CX5b ){
          //for openmp just sets num spikes to 0 
	  cells[i]->num_spikes=0;
	}
      }

      if(t >= homeo.start){
        //homeo.boost_activity(cell_sizes,cells,num_cells); //causes homeostatic changes
        homeo.boost_local(clust_ids,cell_sizes,cells,num_cells);
      }
    }

    // Cut Cell Recovery; once for every 5 spike free checks (4 seconds)
    // Only run this 10 times since going in increments of 10...
    // if((t>=(homeo.undercut_start+1000)) && ((ii % (2*homeo.fre_window)) == 0) && (recovery_count <=10)){
    //     homeo.recovery(cut_ids,cell_sizes,cells,num_cells,1);
    //     recovery_count = recovery_count+1;
    //     printf("Recovery Count = %d \n", recovery_count);
    // }

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

