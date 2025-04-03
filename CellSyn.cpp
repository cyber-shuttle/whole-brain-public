#include "CellSyn.h"


CellSyn::CellSyn(int i,int j,int b_size,int t_size,Cell_Info *cell_info,int type){    

  rand_buffer = new drand48_data;
  srand48_r(2000,rand_buffer);
  //printf("srand 48 \n");
  //srand48_r(time(NULL),rand_buffer);

  this->old_v=0;

  this->type = type;
  this->cell_info = cell_info;
  this->num_syns = this->cell_info->total_connects;

  num_spikes = 0;
    
  flip = 0;
  flip_for_maps = 0;

  m = i;
  n = j;
  lhs_size = t_size;
  y = new double[lhs_size];
  f = new double[lhs_size];
  s = new double[lhs_size];
  yk = new double[lhs_size];
  for(int iter = 0; iter < lhs_size; iter++){
	y[iter] = 0;
	f[iter] = 0;
	s[iter] = 0;
	yk[iter] = 0;
  }

  current_step = 0;
  cur_time = 0;
  xk = 0;
  this->tau = TAU;
  num_b_eq = b_size;
}

//TODO this function should use enum Cell_Type instead of plain numbers 
//this function creates all the cells(CellSyns) and does some basic initialization
//all arguments are inputs.
CellSyn* CellSyn::initialize(int type,int i, int j, Cell_Info ****cells_info,string output_location ){
  
  CellSyn *me;

  switch(type){    
  case E_RE:
    me = new REsyn(i,j,N_RE,N_RE,cells_info[E_RE][i][j],type);
    return me;

  case E_REa:
    me = new REsyn(i,j,N_RE,N_RE,cells_info[E_REa][i][j],type);
    return me;

  case E_TC:
    me = new TCcore(i,j,N_TC,N_TC+N_GB*cells_info[E_TC][i][j]->num_connects_norm[E_RE],cells_info[E_TC][i][j],type);
    return me;

  case E_TCa: 
    me = new TCmatrix(i,j,N_TCa,N_TCa+N_GB*cells_info[E_TCa][i][j]->num_connects_norm[E_REa],cells_info[E_TCa][i][j],type);
    return me;


  case E_rRE:
    me = new REsyn(i,j,N_RE,N_RE,cells_info[E_rRE][i][j],type);
    return me;

  case E_rREa:
    me = new REsyn(i,j,N_RE,N_RE,cells_info[E_rREa][i][j],type);
    return me;

  case E_rTC:
    me = new TCcore(i,j,N_TC,N_TC+N_GB*cells_info[E_rTC][i][j]->num_connects_norm[E_rRE],cells_info[E_rTC][i][j],type);
    return me;

  case E_rTCa: 
    me = new TCmatrix(i,j,N_TCa,N_TCa+N_GB*cells_info[E_rTCa][i][j]->num_connects_norm[E_rREa],cells_info[E_rTCa][i][j],type);
    return me;


  case E_CX:
	if (I_HH==1){
	  me = new CXsyn(i,j,N_CX,N_CX+N_GB*cells_info[E_CX][i][j]->num_connects_norm[E_IN],cells_info[E_CX][i][j],type);}
	else {
	  // printf("\n Creating maps");
	  me = new CXsyn_Map(i,j,N_CX,N_CX+N_GB*cells_info[E_CX][i][j]->num_connects_norm[E_IN],cells_info[E_CX][i][j],type);}
    return me;

  case E_CXa:
	if (I_HH==1){
	  me = new CXasyn(i,j,N_CXa,N_CXa+N_GB*cells_info[E_CXa][i][j]->num_connects_norm[E_INa],cells_info[E_CXa][i][j],type);}
	else{
	  me = new CXasyn_Map(i,j,N_CXa,N_CXa+N_GB*cells_info[E_CXa][i][j]->num_connects_norm[E_INa],cells_info[E_CXa][i][j],type);}
    return me;

  case E_IN:
	if (I_HH==1){
	  me = new INsynCore(i,j,N_IN,N_IN,cells_info[E_IN][i][j],type);}
	else{
	  me = new INsynCore_Map(i,j,N_IN,N_IN,cells_info[E_IN][i][j],type);}
    return me;

  case E_INa:
	if (I_HH==1){
	  me = new INsynMatrix(i,j,N_INa,N_INa,cells_info[E_INa][i][j],type);}
	else{
	  me = new INsynMatrix_Map(i,j,N_INa,N_INa,cells_info[E_INa][i][j],type); }
    return me;

//  case E_CX6:
//	if (I_HH==1){
//	  me = new CXsyn6(i,j,N_CX6,N_CX6+N_GB*cells_info[E_CX6][i][j]->num_connects_norm[E_IN6],cells_info[E_CX6][i][j],type);}
//	else{
//	  me = new CXsyn6_Map(i,j,N_CX6,N_CX6+N_GB*cells_info[E_CX6][i][j]->num_connects_norm[E_IN6],cells_info[E_CX6][i][j],type);}
//    return me;
//
//  case E_IN6:
//	if (I_HH==1){
//	  me = new INsyn6(i,j,N_IN6,N_IN6,cells_info[E_IN6][i][j],type);}
//	else{
//	  me = new INsyn6_Map(i,j,N_IN6,N_IN6,cells_info[E_IN6][i][j],type);}
//    return me;



  case E_rCX:
	if (I_HH==1){
	  me = new CXsyn(i,j,N_CX,N_CX+N_GB*cells_info[E_rCX][i][j]->num_connects_norm[E_rIN],cells_info[E_rCX][i][j],type);}
	else {
	  // printf("\n Creating maps");
	  me = new CXsyn_Map(i,j,N_CX,N_CX+N_GB*cells_info[E_rCX][i][j]->num_connects_norm[E_rIN],cells_info[E_rCX][i][j],type);}
    return me;

  case E_rCXa:
	if (I_HH==1){
	  me = new CXasyn(i,j,N_CXa,N_CXa+N_GB*cells_info[E_rCXa][i][j]->num_connects_norm[E_rINa],cells_info[E_rCXa][i][j],type);}
	else{
	  me = new CXasyn_Map(i,j,N_CXa,N_CXa+N_GB*cells_info[E_rCXa][i][j]->num_connects_norm[E_rINa],cells_info[E_rCXa][i][j],type);}
    return me;

  case E_rIN:
	if (I_HH==1){
	  me = new INsynCore(i,j,N_IN,N_IN,cells_info[E_rIN][i][j],type);}
	else{
	  me = new INsynCore_Map(i,j,N_IN,N_IN,cells_info[E_rIN][i][j],type);}
    return me;

  case E_rINa:
	if (I_HH==1){
	  me = new INsynMatrix(i,j,N_INa,N_INa,cells_info[E_rINa][i][j],type);}
	else{
	  me = new INsynMatrix_Map(i,j,N_INa,N_INa,cells_info[E_rINa][i][j],type); }
    return me;

//   case E_rCX6:
// 	if (I_HH==1){
// 	  me = new CXsyn6(i,j,N_CX6,N_CX6+N_GB*cells_info[E_rCX6][i][j]->num_connects_norm[E_rIN6],cells_info[E_rCX6][i][j],type);}
// 	else{
// 	  me = new CXsyn6_Map(i,j,N_CX6,N_CX6+N_GB*cells_info[E_rCX6][i][j]->num_connects_norm[E_rIN6],cells_info[E_rCX6][i][j],type);}
//     return me;
// 
//   case E_rIN6:
// 	if (I_HH==1){
// 	  me = new INsyn6(i,j,N_IN6,N_IN6,cells_info[E_rIN6][i][j],type);}
// 	else{
// 	  me = new INsyn6_Map(i,j,N_IN6,N_IN6,cells_info[E_rIN6][i][j],type);}
//     return me;






  case E_CX2:
	if (I_HH==1){
	  me = new CX2syn(i,j,N_CX2,N_CX2+N_GB*cells_info[E_CX2][i][j]->num_connects_norm[E_IN2],cells_info[E_CX2][i][j],type);}
	else{
	  me = new CX2syn_Map(i,j,N_CX2,N_CX2+N_GB*cells_info[E_CX2][i][j]->num_connects_norm[E_IN2],cells_info[E_CX2][i][j],type);}
    return me;

  case E_IN2:
	if (I_HH==1){
	  me = new IN2syn(i,j,N_IN2,N_IN2,cells_info[E_IN2][i][j],type);}
	else{
	  me = new IN2syn_Map(i,j,N_IN2,N_IN2,cells_info[E_IN2][i][j],type);}
    return me;



  case E_CX3:
	if (I_HH==1){
	  me = new CX3syn(i,j,N_CX3,N_CX3+N_GB*cells_info[E_CX3][i][j]->num_connects_norm[E_IN3],cells_info[E_CX3][i][j],type);}
	else{
	  me = new CX3syn_Map(i,j,N_CX3,N_CX3+N_GB*cells_info[E_CX3][i][j]->num_connects_norm[E_IN3],cells_info[E_CX3][i][j],type);}
    return me;

  case E_IN3:
	if (I_HH==1){
	  me = new IN3syn(i,j,N_IN3,N_IN3,cells_info[E_IN3][i][j],type);}
	else{
	  me = new IN3syn_Map(i,j,N_IN3,N_IN3,cells_info[E_IN3][i][j],type);}
    return me;



  case E_CX4:
	if (I_HH==1){
	  me = new CX4syn(i,j,N_CX4,N_CX4+N_GB*cells_info[E_CX4][i][j]->num_connects_norm[E_IN4],cells_info[E_CX4][i][j],type);}
	else{
	  me = new CX4syn_Map(i,j,N_CX4,N_CX4+N_GB*cells_info[E_CX4][i][j]->num_connects_norm[E_IN4],cells_info[E_CX4][i][j],type);}
    return me;

  case E_IN4:
	if (I_HH==1){
	  me = new IN4syn(i,j,N_IN4,N_IN4,cells_info[E_IN4][i][j],type);}
	else{
	  me = new IN4syn_Map(i,j,N_IN4,N_IN4,cells_info[E_IN4][i][j],type);}
    return me;



  case E_CX5a:
	if (I_HH==1){
	  me = new CX5asyn(i,j,N_CX5a,N_CX5a+N_GB*cells_info[E_CX5a][i][j]->num_connects_norm[E_IN5a],cells_info[E_CX5a][i][j],type);}
	else{
	  me = new CX5asyn_Map(i,j,N_CX5a,N_CX5a+N_GB*cells_info[E_CX5a][i][j]->num_connects_norm[E_IN5a],cells_info[E_CX5a][i][j],type);}
    return me;

  case E_IN5a:
	if (I_HH==1){
	  me = new IN5asyn(i,j,N_IN5a,N_IN5a,cells_info[E_IN5a][i][j],type);}
	else{
	  me = new IN5asyn_Map(i,j,N_IN5a,N_IN5a,cells_info[E_IN5a][i][j],type);}
    return me;



  case E_CX5b:
	if (I_HH==1){
	  me = new CX5bsyn(i,j,N_CX5b,N_CX5b+N_GB*cells_info[E_CX5b][i][j]->num_connects_norm[E_IN5b],cells_info[E_CX5b][i][j],type);}
	else{
	  me = new CX5bsyn_Map(i,j,N_CX5b,N_CX5b+N_GB*cells_info[E_CX5b][i][j]->num_connects_norm[E_IN5b],cells_info[E_CX5b][i][j],type);}
    return me;

  case E_IN5b:
	if (I_HH==1){
	  me = new IN5bsyn(i,j,N_IN5b,N_IN5b,cells_info[E_IN5b][i][j],type);}
	else{
	  me = new IN5bsyn_Map(i,j,N_IN5b,N_IN5b,cells_info[E_IN5b][i][j],type);}
    return me;



  case E_CX6:
	if (I_HH==1){
	  me = new CX6syn(i,j,N_CX6,N_CX6+N_GB*cells_info[E_CX6][i][j]->num_connects_norm[E_IN6],cells_info[E_CX6][i][j],type);}
	else{
	  me = new CX6syn_Map(i,j,N_CX6,N_CX6+N_GB*cells_info[E_CX6][i][j]->num_connects_norm[E_IN6],cells_info[E_CX6][i][j],type);}
    return me;

  case E_IN6:
	if (I_HH==1){
	  me = new IN6syn(i,j,N_IN6,N_IN6,cells_info[E_IN6][i][j],type);}
	else{
	  me = new IN6syn_Map(i,j,N_IN6,N_IN6,cells_info[E_IN6][i][j],type);}
    return me;





  case E_rCX2:
	if (I_HH==1){
	  me = new CX2syn(i,j,N_CX2,N_CX2+N_GB*cells_info[E_rCX2][i][j]->num_connects_norm[E_rIN2],cells_info[E_rCX2][i][j],type);}
	else{
	  me = new CX2syn_Map(i,j,N_CX2,N_CX2+N_GB*cells_info[E_rCX2][i][j]->num_connects_norm[E_rIN2],cells_info[E_rCX2][i][j],type);}
    return me;

  case E_rIN2:
	if (I_HH==1){
	  me = new IN2syn(i,j,N_IN2,N_IN2,cells_info[E_rIN2][i][j],type);}
	else{
	  me = new IN2syn_Map(i,j,N_IN2,N_IN2,cells_info[E_rIN2][i][j],type);}
    return me;



  case E_rCX3:
	if (I_HH==1){
	  me = new CX3syn(i,j,N_CX3,N_CX3+N_GB*cells_info[E_rCX3][i][j]->num_connects_norm[E_rIN3],cells_info[E_rCX3][i][j],type);}
	else{
	  me = new CX3syn_Map(i,j,N_CX3,N_CX3+N_GB*cells_info[E_rCX3][i][j]->num_connects_norm[E_rIN3],cells_info[E_rCX3][i][j],type);}
    return me;

  case E_rIN3:
	if (I_HH==1){
	  me = new IN3syn(i,j,N_IN3,N_IN3,cells_info[E_rIN3][i][j],type);}
	else{
	  me = new IN3syn_Map(i,j,N_IN3,N_IN3,cells_info[E_rIN3][i][j],type);}
    return me;



  case E_rCX4:
	if (I_HH==1){
	  me = new CX4syn(i,j,N_CX4,N_CX4+N_GB*cells_info[E_rCX4][i][j]->num_connects_norm[E_rIN4],cells_info[E_rCX4][i][j],type);}
	else{
	  me = new CX4syn_Map(i,j,N_CX4,N_CX4+N_GB*cells_info[E_rCX4][i][j]->num_connects_norm[E_rIN4],cells_info[E_rCX4][i][j],type);}
    return me;

  case E_rIN4:
	if (I_HH==1){
	  me = new IN4syn(i,j,N_IN4,N_IN4,cells_info[E_rIN4][i][j],type);}
	else{
	  me = new IN4syn_Map(i,j,N_IN4,N_IN4,cells_info[E_rIN4][i][j],type);}
    return me;



  case E_rCX5a:
	if (I_HH==1){
	  me = new CX5asyn(i,j,N_CX5a,N_CX5a+N_GB*cells_info[E_rCX5a][i][j]->num_connects_norm[E_rIN5a],cells_info[E_rCX5a][i][j],type);}
	else{
	  me = new CX5asyn_Map(i,j,N_CX5a,N_CX5a+N_GB*cells_info[E_rCX5a][i][j]->num_connects_norm[E_rIN5a],cells_info[E_rCX5a][i][j],type);}
    return me;

  case E_rIN5a:
	if (I_HH==1){
	  me = new IN5asyn(i,j,N_IN5a,N_IN5a,cells_info[E_rIN5a][i][j],type);}
	else{
	  me = new IN5asyn_Map(i,j,N_IN5a,N_IN5a,cells_info[E_rIN5a][i][j],type);}
    return me;



  case E_rCX5b:
	if (I_HH==1){
	  me = new CX5bsyn(i,j,N_CX5b,N_CX5b+N_GB*cells_info[E_rCX5b][i][j]->num_connects_norm[E_rIN5b],cells_info[E_rCX5b][i][j],type);}
	else{
	  me = new CX5bsyn_Map(i,j,N_CX5b,N_CX5b+N_GB*cells_info[E_rCX5b][i][j]->num_connects_norm[E_rIN5b],cells_info[E_rCX5b][i][j],type);}
    return me;

  case E_rIN5b:
	if (I_HH==1){
	  me = new IN5bsyn(i,j,N_IN5b,N_IN5b,cells_info[E_rIN5b][i][j],type);}
	else{
	  me = new IN5bsyn_Map(i,j,N_IN5b,N_IN5b,cells_info[E_rIN5b][i][j],type);}
    return me;



  case E_rCX6:
	if (I_HH==1){
	  me = new CX6syn(i,j,N_CX6,N_CX6+N_GB*cells_info[E_rCX6][i][j]->num_connects_norm[E_rIN6],cells_info[E_rCX6][i][j],type);}
	else{
	  me = new CX6syn_Map(i,j,N_CX6,N_CX6+N_GB*cells_info[E_rCX6][i][j]->num_connects_norm[E_rIN6],cells_info[E_rCX6][i][j],type);}
    return me;

  case E_rIN6:
	if (I_HH==1){
	  me = new IN6syn(i,j,N_IN6,N_IN6,cells_info[E_rIN6][i][j],type);}
	else{
	  me = new IN6syn_Map(i,j,N_IN6,N_IN6,cells_info[E_rIN6][i][j],type);}
    return me;



  default:
    printf("tried to create a non existant cell type in initialize\n");
    exit(1);
  }
}


// Called from memory function - syns[i] = initialize_synapses(syn_info[i]);
Syn* CellSyn::initialize_synapses(Syn_Info *syn_info){

  Syn *syni=NULL;
	  
  if(syn_info->type == E_AMPA){
	syni = new AMPA(syn_info);
  }else if(syn_info->type == E_GiantAMPA){
  	syni = new Syn(syn_info);
	syni->giant = get_GiantSyn(syn_info->from_type,syn_info->type);
  }else if(syn_info->type == E_GiantAMPAMap){
  	syni = new Syn(syn_info);
	syni->giant = get_GiantSyn(syn_info->from_type,syn_info->type);
  }else if(syn_info->type == E_GABA_A){
	syni = new GABA_A(syn_info);
  }else if(syn_info->type == E_AMPA_D3){ 
	syni = new AMPA_D3(syn_info);
  }else if(syn_info->type == E_AMPA_D2){
	syni = new AMPA_D2(syn_info); 	  	
  }else if(syn_info->type == E_AMPA_D1){
	syni = new AMPA_D1(syn_info); 
  }else if(syn_info->type == E_NMDA_D1){
	syni = new NMDA_D1(syn_info); 
  }else if(syn_info->type == E_GABA_B){
	syni = new GB(syn_info); 
  }else if(syn_info->type == E_GABA_A_D2){
	syni = new GABA_A_D2(syn_info); 
  }else if(syn_info->type == E_GABAMap_A_D2){
	syni = new GABAAmapD1(syn_info); 
  }else if(syn_info->type == E_AMPAMapD){
        syni = new AMPAmapD(syn_info);
  }else if(syn_info->type == E_AMPAMap_D1){
	syni = new AMPAmapD1(syn_info);
  }else if(syn_info->type == E_NMDAMap_D1){
	syni = new NMDAmapD1(syn_info); 
  }else if(syn_info->type == E_GAP){
	syni = new GAP(syn_info);
  }else{
	printf("Synapse initialization failed: type %d from neuron: %d to neuron %d \n", syn_info->type, syn_info->from_type, syn_info->to_type);
	exit(1);
  }

  return syni;
}


GiantSyn* CellSyn::define_GiantSyn(enum Cell_Type cell, enum Syn_Type syn){
  GiantSyn *S=NULL;
  //TODO add case statement for more synapses
  if (syn == E_GiantAMPAMap)
	S = new GiantSynAMPAMap();
  else if (syn == E_GiantAMPA)
	S = new GiantSynAMPA();

  if (!S) {printf("Error Giant synapse initialization!\n");exit(1);}
  return S;
}


void CellSyn::signal_spike(){
    
  double message = base_cell->get_v_soma();

  if(this->old_v < 0 && message >= 0){
	this->old_v = message;
	flip = 1;
	flip_for_maps=1;
	num_spikes = num_spikes +1;
	print_spike_time(cur_time,this->type,this->m,this->n);
  }else{
	flip = 0;
	this->old_v = message;
	return;
  }
 
}

void CellSyn::reset_flip_maps(){
    
  if(flip_for_maps == 1){
	flip_for_maps=0;
  }
 
}

void CellSyn::zero_temps(){
  // not sure if it is need 
  int i = 0;
  for(i =0; i<lhs_size; i++){
	s[i]=0;
	yk[i]=0;
  }
}

void CellSyn::reset_y(){
  int i = 0;
  for(i = 0; i < lhs_size; ++i){ 
	y[i] = y[i] + (tau/6.0)*(s[i] + f[i]);
  }
}

void CellSyn::step(){

  if (this->ismap){
	this->calc(cur_time,y,f,2);
	cur_time=cur_time+TAU_Map;
  } 
  else{
	int i = 0;
	if(current_step == 0){
	  this->calc(cur_time,y,f,1);
	  for(i = 0; i < lhs_size; ++i){ 
		s[i] = f[i]; 
		yk[i] = y[i] + (tau/2.0)*f[i]; 
	  }
	  xk = cur_time + tau/2.0;
	}
  
	if(current_step ==1){
	  this->calc(xk,yk,f,0);
	  for(i = 0; i < lhs_size; ++i) { 
		s[i] = s[i] + 2.0*f[i]; 
		yk[i] = y[i] + (tau/2.0)*f[i]; 
	  }
	}

	if(current_step == 2){
	  this->calc(xk,yk,f,0);
	  for(i = 0; i < lhs_size; ++i) { 
		s[i] = s[i] + 2.0*f[i]; 
		yk[i] = y[i] + tau*f[i]; 
	  }
	  xk = cur_time + tau;
	}

	if(current_step == 3){
	  this->calc(xk,yk,f,0);
	  cur_time = cur_time + tau;
	}
    
	current_step = current_step +1;
	//TODO kill this condition
	if(current_step > 3){
	  current_step = 0;
	}
  }

}


extern CellSyn **cells;
extern Cell_Info ****cells_info;
//this function runs thought each synapse for a cell and calculates its current before calculating the cell itself
//step == 0 <=> RK 2,3,4 step;  step == 1 <=> RK 1; step == 2 <=> MAP
void CellSyn::calc(double x, double *y_ini, double *f_ini, int step){
  
  double current = 0.0;
  int input_soma = 0; 	  // input_soma - pre-synaptic spike
  double input_dend = 0;  
  int i = 0;
  int k = 0;

  e_current=0; i_current=0;

  double vd=base_cell->get_v_dend();

  // reset giant synapse when step
  if (step>0)
    reset_GiantSyns();

  // Run through all synapses
  for(i =0; i < num_syns; i++){

	int syn_type=syns[i]->type;

    double tmp_current=0;

	int from=syns[i]->from_cell;
	GiantSyn *gs = syns[i]->giant; // Check if this synapse part of giant synapse
	  
	if(syn_type == E_GAP){
	  input_dend = cells[from]->base_cell->get_v_dend();
      tmp_current = ((GAP*)syns[i])->calc_gap(x,vd, input_dend); // /max_connect_gap[type][ftype];
      current = current - tmp_current;
	  continue; 
	}

	if (step>0) {
	  if (this->ismap)
		input_soma = cells[from]->flip_for_maps;
	  else {
		input_soma = cells[from]->flip;	// check if input neuron spiked
               }

	  // if it is giant synapse add spikes
	  if (gs && input_soma)
	    gs->add_spike();
    } else {
        input_soma=0; // Ensures the spike is used only when in step>0, since the reset occurs at step=1 or 2
    }

	if(gs){}
	else{
	  if(syn_type == E_GABA_B){
		//current = current - ((GB*)syns[i])->calc_gb(x,y+num_b_eq+2*k, f+num_b_eq+2*k, vd, input_soma)/max_connect[type][ftype];
		// scaling now implimented in initialization of network - Prevents division every time calc is called
        tmp_current = ((GB*)syns[i])->calc_gb(x,y+num_b_eq+2*k, f+num_b_eq+2*k, vd, input_soma);
        current = current - tmp_current;
		k = k + 1;
	  }else{ 
		//general case - for all other synapses
		// current = current - (syns[i])->calc(x, vd, input_soma)/max_connect[type][ftype];
        // tmp_current = (syns[i])->calc(x, vd, input_soma);
        // current = current - tmp_current;

        (syns[i])->appliedCalc = false; //mark synapse and not having applied calc function
                                      // due to delay we do not know if the current synapse 
                                      // has a spike at this time or not

        if(input_soma == 1)
        {
            storedSpikes.push(StoredSpike(syns[i], x + syns[i]->delay));
        }
      }
    }
// if ((x>1000 && x<1010) && (syn_type == E_AMPAMap_D1)){
//    printf("Synaptic strenth for neuron %d : %lf", from, syns[i]->strength);
// }

    // if (step>0) {
    //   if ((syn_type == E_AMPA) || (syn_type == E_NMDA_D1) || (syn_type == E_AMPA_D1) || (syn_type == E_AMPA_D2) || (syn_type == E_AMPA_D3) ||  (syn_type == E_AMPAMap_D1) ||  (syn_type == E_NMDAMap_D1) ||  (syn_type == E_GiantAMPA) ||  (syn_type == E_GiantAMPAMap)) {
    //     e_current=e_current-tmp_current;
	//   }
    //   else{
    //     i_current=i_current-tmp_current;
    //   }
    // }

    // if (type==E_CX || type==E_CXa || type==E_CX6) {
    //   y_var= base_cell->get_y_var();
	// }

  }

  // apply synapses with delay that should apply a spike at this timestep
  while(storedSpikes.size() > 0 && storedSpikes.top().timeApplySpike <= x)
  {
    double tmp_current = storedSpikes.top().synapse->calc(x, vd, 1);
    int syn_type = storedSpikes.top().synapse->type;
    storedSpikes.top().synapse->appliedCalc = true;
    storedSpikes.pop();
    current = current - tmp_current;

    // current stuff
    if (step > 0)
    {
      if ((syn_type == E_AMPA) || (syn_type == E_NMDA_D1) || (syn_type == E_AMPA_D1) || (syn_type == E_AMPA_D2) || (syn_type == E_AMPA_D3) || (syn_type == E_AMPAMap_D1) || (syn_type == E_NMDAMap_D1) || (syn_type == E_GiantAMPA) || (syn_type == E_GiantAMPAMap))
      {
        e_current = e_current - tmp_current;
      }
      else
      {
        i_current = i_current - tmp_current;
      }
    }
    // e_current = current;
    // i_current = current;
    if (type == E_CX || type == E_CXa || type == E_CX6 || type == E_CX3 || type == E_CX4 || type == E_CX5a || type == E_CX5b )
    {
      y_var = base_cell->get_y_var();
    }
  }
  
  // apply synapses with delay that should apply a no spike at this timestep
  for(i =0; i < num_syns; i++)
  {
    if( (syns[i])->appliedCalc == false)
    {
      int syn_type = syns[i]->type;
      (syns[i])->appliedCalc = true;
      double tmp_current = (syns[i])->calc(x, vd, 0); // current if there is no spike
      current = current - tmp_current;

      // current stuff
      if (step > 0)
      {
        if ((syn_type == E_AMPA) || (syn_type == E_NMDA_D1) || (syn_type == E_AMPA_D1) || (syn_type == E_AMPA_D2) || (syn_type == E_AMPA_D3) || (syn_type == E_AMPAMap_D1) || (syn_type == E_NMDAMap_D1) || (syn_type == E_GiantAMPA) || (syn_type == E_GiantAMPAMap))
        {
          e_current = e_current - tmp_current;
        }
        else
        {
          i_current = i_current - tmp_current;
        }
      }
      // e_current = current;
      // i_current = current;
      if (type == E_CX || type == E_CXa || type == E_CX6 || type == E_CX3 || type == E_CX4 || type == E_CX5a || type == E_CX5b )
      {
        y_var = base_cell->get_y_var();
      }
    }

  }

  current =  current - calc_GiantSyns();
  base_cell->calc(x, current, y_ini, f_ini);

}

GiantSyn* CellSyn::find_GiantSyn(enum Cell_Type cell, enum Syn_Type syn){
  list<GiantSyn*>::iterator it=giant_syns.begin();
  list<GiantSyn*>::iterator end=giant_syns.end();
  for (;it!=end;it++)
	if ( (*it)->from_type == cell && (*it)->syn_type == syn )
	  return (*it);
  return NULL;
}

GiantSyn* CellSyn::get_GiantSyn(enum Cell_Type cell, enum Syn_Type syn){
  //Here can be placed #ifndef GiantSyn return NULL; 
  GiantSyn *S=NULL;
  S = find_GiantSyn(cell,syn);
  if (!S) {
	S = define_GiantSyn(cell,syn);
	giant_syns.push_back(S);
  }
  return S;
}

void CellSyn::reset_GiantSyns(){
  list<GiantSyn*>::iterator it=giant_syns.begin();
  list<GiantSyn*>::iterator end=giant_syns.end();
  for (;it!=end;it++)
	(*it)->reset_spikes();
}

double CellSyn::calc_GiantSyns(){
  double sum = 0;
  list<GiantSyn*>::iterator it=giant_syns.begin();
  list<GiantSyn*>::iterator end=giant_syns.end();
  for (;it!=end;it++)
	sum += (*it)->calc();

  return sum;
}



