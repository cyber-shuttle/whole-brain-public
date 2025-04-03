#include "io.h"
#include <sstream>
#include "params.h"
#include <string.h>

double print_receive(int m , int n, enum Cell_Type type); //still in main
double print_receive_var(int m , int n, enum Cell_Type type, int index); //still in main

void print_receive_gsynapse(int m , int n, enum Cell_Type type, enum Syn_Type stype, FILE *fp);
void print_receive_gsynapse_index(int m , int n, enum Cell_Type from_type, enum Cell_Type type, enum Syn_Type stype, FILE *fp);
void print_receive_gsynapse_multLayer(int m , int n, enum Cell_Type from_type, enum Cell_Type type, enum Syn_Type stype, FILE *fp);

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
		int&  homeo_num_regions

		){

	add_double_param(  tmax  );
	add_double_param(  t3D   );
	add_double_param(  ttime );
	add_int_param(  num_mp_threads );
	add_int_param(  print_c_sten   );
	add_int_param(  fre_print_cs   );
	add_int_param(  LFP_local_field_effect );
	add_double_param(  LFP_lfp_scale );
	add_int_param(  LFP_num_field_layers );
	add_int_param(  homeo_boost  );
	add_double_param(  homeo_amp_boost  );
	add_double_param(  homeo_con_boost  );
	add_double_param(  homeo_fre_boost  );
	add_double_param(  homeo_target_f   );
	add_int_param(  homeo_fre_window    );
	add_int_param(  homeo_num_regions   );

	assert(load_parameters(argv[1]));
	assert(cmdline_parameters(argc,argv));
	print_parameters();

} //load_input_params

//giant list of all the output files some are probably not used
FILE *flocal, *f2, *f3, *f4, *f6, *f7, *f8, *f9, *f10, *f11, *f12, *f13, *f14, *f15, *f16, *f17, *f18, *f19,*f20,*f21,*f22,*f23,*f24,*f25,*f26,*f27,*f28, *f29, *f30, *f31, *f32; 
FILE *cx3, *cx4, *cx5a, *cx5b, *in3, *in4, *in5a, *in5b; // *cx2, *in2;
FILE *cx3Currents, *cx4Currents, *cx5aCurrents, *cx5bCurrents; 
//FILE *in3Currents, *in4Currents, *in5aCurrents, *in5bCurrents; 

FILE *fconn, *fSpikeTime;
FILE *fconn_cx2_cx2, *fconn_cx3_cx2, *fconn_cx4_cx2, *fconn_cx5a_cx2, *fconn_cx5b_cx2, *fconn_cx6_cx2;
FILE *fconn_cx2_cx3, *fconn_cx3_cx3, *fconn_cx4_cx3, *fconn_cx5a_cx3, *fconn_cx5b_cx3, *fconn_cx6_cx3;


FILE *fconduct_cx2_cx2, *fconduct_cx3_cx2, *fconduct_cx4_cx2, *fconduct_cx5a_cx2, *fconduct_cx5b_cx2, *fconduct_cx6_cx2;
FILE *fconduct_cx2_cx3, *fconduct_cx3_cx3, *fconduct_cx4_cx3, *fconduct_cx5a_cx3, *fconduct_cx5b_cx3, *fconduct_cx6_cx3;



void print_connectivity(Pair *cell_sizes){

  int i =0;
  int j =0;
  
  for(i = 0; i < cell_sizes[E_CX].x; ++i){
    for(j = 0; j < cell_sizes[E_CX].y; ++j){
      print_receive_gsynapse_index(i, j, E_IN, E_CX, E_GABAMap_A_D2, fconn);
    }    
  }
  fprintf(fconn,"\n");
  
//  j=0; 
//  for(i = 0; i < cell_sizes[E_CX].x; ++i){
//      print_receive_gsynapse_index(i, j, E_CX, E_CX, E_AMPAMap_D1, fconn_cx2_cx2);    
//  }
//  fprintf(fconn_cx2_cx2,"\n");
//
//  j=0; 
//  for(i = 0; i < cell_sizes[E_CX].x; ++i){
//      print_receive_gsynapse_index(i, j, E_CX3, E_CX, E_AMPAMap_D1, fconn_cx3_cx2);    
//  }
//  fprintf(fconn_cx3_cx2,"\n");
//
//  j=0; 
//  for(i = 0; i < cell_sizes[E_CX].x; ++i){
//      print_receive_gsynapse_index(i, j, E_CX4, E_CX, E_AMPAMap_D1, fconn_cx4_cx2);    
//  }
//  fprintf(fconn_cx4_cx2,"\n");
//
//  j=0; 
//  for(i = 0; i < cell_sizes[E_CX].x; ++i){
//      print_receive_gsynapse_index(i, j, E_CX5a, E_CX, E_AMPAMap_D1, fconn_cx5a_cx2);    
//  }
//  fprintf(fconn_cx5a_cx2,"\n");
//
//  j=0; 
//  for(i = 0; i < cell_sizes[E_CX].x; ++i){
//      print_receive_gsynapse_index(i, j, E_CX5b, E_CX, E_AMPAMap_D1, fconn_cx5b_cx2);    
//  }
//  fprintf(fconn_cx5b_cx2,"\n");
//
//  j=0; 
//  for(i = 0; i < cell_sizes[E_CX].x; ++i){
//      print_receive_gsynapse_index(i, j, E_CX6, E_CX, E_AMPAMap_D1, fconn_cx6_cx2);    
//  }
//  fprintf(fconn_cx6_cx2,"\n");
//
//
//// //////////////////////
//  j=0; 
//  for(i = 0; i < cell_sizes[E_CX3].x; ++i){
//      print_receive_gsynapse_index(i, j, E_CX, E_CX3, E_AMPAMap_D1, fconn_cx2_cx3);    
//  }
//  fprintf(fconn_cx2_cx3,"\n");
//
//  j=0; 
//  for(i = 0; i < cell_sizes[E_CX3].x; ++i){
//      print_receive_gsynapse_index(i, j, E_CX3, E_CX3, E_AMPAMap_D1, fconn_cx3_cx3);    
//  }
//  fprintf(fconn_cx3_cx3,"\n");
//
//  j=0; 
//  for(i = 0; i < cell_sizes[E_CX3].x; ++i){
//      print_receive_gsynapse_index(i, j, E_CX4, E_CX3, E_AMPAMap_D1, fconn_cx4_cx3);    
//  }
//  fprintf(fconn_cx4_cx3,"\n");
//
//  j=0; 
//  for(i = 0; i < cell_sizes[E_CX3].x; ++i){
//      print_receive_gsynapse_index(i, j, E_CX5a, E_CX3, E_AMPAMap_D1, fconn_cx5a_cx3);    
//  }
//  fprintf(fconn_cx5a_cx3,"\n");
//
//  j=0; 
//  for(i = 0; i < cell_sizes[E_CX3].x; ++i){
//      print_receive_gsynapse_index(i, j, E_CX5b, E_CX3, E_AMPAMap_D1, fconn_cx5b_cx3);    
//  }
//  fprintf(fconn_cx5b_cx3,"\n");
//
//  j=0; 
//  for(i = 0; i < cell_sizes[E_CX3].x; ++i){
//      print_receive_gsynapse_index(i, j, E_CX6, E_CX3, E_AMPAMap_D1, fconn_cx6_cx3);    
//  }
//  fprintf(fconn_cx6_cx3,"\n");

} 

void print_occ(Pair *cell_sizes){

	int i =0;
	int j =0;

	for(i = 0; i < cell_sizes[E_RE].x; ++i){
		for(j = 0; j < cell_sizes[E_RE].y; ++j){
			fprintf(f6,"%lf ", print_receive(i,j,E_RE));
		}
		fprintf(f6,"\n");
	}

	for(i = 0; i < cell_sizes[E_TC].x; ++i){
		for(j = 0; j < cell_sizes[E_TC].y; ++j){
			fprintf(f8,"%lf ", print_receive(i,j,E_TC));
		}
		fprintf(f8,"\n");
	}

	for(i = 0; i < cell_sizes[E_TCa].x; ++i){
		for(j = 0; j < cell_sizes[E_TC].y; ++j){
			fprintf(f16,"%lf ", print_receive(i,j,E_TCa));
		}
		fprintf(f16,"\n");
	}

	for(i = 0; i < cell_sizes[E_CX].x; ++i){
		for(j = 0; j < cell_sizes[E_CX].y; ++j){
			fprintf(f10,"%lf ", print_receive(i,j,E_CX));
		}
		fprintf(f10,"\n");
	}

//	for(i = 0; i < cell_sizes[E_CXa].x; ++i){
//		for(j = 0; j < cell_sizes[E_CXa].y; ++j){
//			fprintf(f14,"%lf ", print_receive(i,j,E_CXa));
//		}
//		fprintf(f14,"\n");
//	}

	for(i = 0; i < cell_sizes[E_IN].x; ++i){
		for(j = 0; j < cell_sizes[E_IN].y; ++j){
			fprintf(f12,"%lf ", print_receive(i,j,E_IN));
		}
		fprintf(f12,"\n");
	}

	// TB IN layers
//	for(i = 0; i < cell_sizes[E_INa].x; ++i){
//		for(j = 0; j < cell_sizes[E_INa].y; ++j){
//			fprintf(f18,"%lf ", print_receive(i,j,E_INa));
//		}
//		fprintf(f18,"\n");
//	}

	for(i = 0; i < cell_sizes[E_CX6].x; ++i){
		for(j = 0; j < cell_sizes[E_CX6].y; ++j){
			fprintf(f24,"%lf ", print_receive(i,j,E_CX6));
		}
		fprintf(f24,"\n");
	}

	for(i = 0; i < cell_sizes[E_IN6].x; ++i){
		for(j = 0; j < cell_sizes[E_IN6].y; ++j){
			fprintf(f26,"%lf ", print_receive(i,j,E_IN6));
		}
		fprintf(f26,"\n");
	}


	//for(i = 0; i < cell_sizes[E_CX].x; ++i){
       //      for(j = 0; j < cell_sizes[E_CX].y; ++j){
       for(i = 100; i < 101; ++i){
               for(j = 0; j < 1; ++j){; // cell_sizes[E_CX].y; ++j){
                       print_receive_gsynapse(i, j, E_CX, E_GABAMap_A_D2, f20);
               }
       }
       fprintf(f20,"\n");







}

void print_spike_time(double time, int type, int m, int n){
  fprintf(fSpikeTime, "%lf %d %d %d \n",time,type,m,n);
}             

void print_freq(double **cx_base_v_SOMA, double **cx5b_base_v_SOMA, Pair *cell_sizes, double const t, int Num_Types){

	int i = 0;
	int j = 0;

	//for(i=0; i <cell_sizes[3].x; i++){
	//for(j=0; j<cell_sizes[3].y; j++){
	//cx_base_v_SOMA[i][j] = print_receive(i,j,E_CX);
	//}
	//}

	//for(i=0; i <cell_sizes[4].x; i++){
	//for(j=0; j<cell_sizes[4].y; j++){
	//cxa_base_v_SOMA[i][j] = print_receive(i,j,E_CXa);
	//}
	//}

	// av=0;
	//for(j = 0; j < cell_sizes[3].y; ++j){
	//for(i = 0; i < cell_sizes[3].x; ++i){
	//av=av + cx_base_v_SOMA[i][j];
	//}
	//}

	//ava=0;
	//for(j = 0; j < cell_sizes[4].y; ++j){
	//for(i = 0; i < cell_sizes[4].x; ++i){
	//  ava=ava + cxa_base_v_SOMA[i][j];
	//}
	//}

	//synAMPA=0;
	//synNMDA=0;
	//synGABA=0;
	//synAMPAtc=0;

	//for(j = 0; j < cell_sizes[3].y; ++j)
	//for(i = 0; i < cell_sizes[3].x; ++i){
	//order of print_receives is important
	//double a_tc_cx_I =print_receive(i,j,E_CX);
	//double ga_in_cx_I =print_receive(i,j,E_CX);
	//double nmda_cx_cx_I =print_receive(i,j,E_CX);
	//double nmda_cxa_cx_I =print_receive(i,j,E_CX);
	//double a_cxa_cx_I =print_receive(i,j,E_CX);
	//double a_cx_cx_I = print_receive(i,j,E_CX);
	//synAMPA=synAMPA+a_cx_cx_I+a_cxa_cx_I;
	//synNMDA=synNMDA+nmda_cx_cx_I+nmda_cxa_cx_I;
	//synGABA=synGABA+ga_in_cx_I;
	//synAMPAtc=synAMPAtc+a_tc_cx_I;
	//}

	fprintf(f7,"%lf ", t);
	for(i = 0; i < cell_sizes[E_RE].x; ++i){
		fprintf(f7,"%lf ",print_receive(i,cell_sizes[E_RE].y/2,E_RE));
	}
	if (Num_Types>E_rRE) for(i = 0; i < cell_sizes[E_rRE].x; ++i){
		fprintf(f7,"%lf ",print_receive(i,cell_sizes[E_rRE].y/2,E_rRE));
	}
	fprintf(f7,"\n");


	fprintf(f29,"%lf ", t);
	for(i = 0; i < cell_sizes[E_REa].x; ++i){
		fprintf(f29,"%lf ",print_receive(i,cell_sizes[E_REa].y/2,E_REa));
	}
	if (Num_Types>E_rREa) for(i = 0; i < cell_sizes[E_rREa].x; ++i){
		fprintf(f29,"%lf ",print_receive(i,cell_sizes[E_rREa].y/2,E_rREa));
	}
	fprintf(f29,"\n");

	fprintf(f9,"%lf ", t);
	for(i = 0; i < cell_sizes[E_TC].x; ++i){
		fprintf(f9,"%lf ", print_receive(i,cell_sizes[E_TC].y/2,E_TC));
	}
	if (Num_Types>E_rTC) for(i = 0; i < cell_sizes[E_rTC].x; ++i){
		fprintf(f9,"%lf ", print_receive(i,cell_sizes[E_rTC].y/2,E_rTC));
	}
	fprintf(f9,"\n");

	fprintf(f17,"%lf ", t);
	for(i = 0; i < cell_sizes[E_TCa].x; ++i){
		fprintf(f17,"%lf ", print_receive(i,cell_sizes[E_TCa].y/2,E_TCa));
	}
	if (Num_Types>E_rTCa) for(i = 0; i < cell_sizes[E_rTCa].x; ++i){
		fprintf(f17,"%lf ", print_receive(i,cell_sizes[E_rTCa].y/2,E_rTCa));
	}
	fprintf(f17,"\n");



	// **** Print for CX cells **** 

	//print time
	fprintf(f11,"%lf ", t);

	//print individual cell data
	for(i = 0; i < cell_sizes[E_CX].x; ++i){
		for(j = 0; j < cell_sizes[E_CX].y; ++j){

			fprintf(f11,"%lf ", print_receive(i,j,E_CX));

			// y variable
			// fprintf(f11,"%lf ", print_receive_var(i,j,E_CX,0));

			// Print currents
			fprintf(f30,"%lf ", print_receive_var(i,j,E_CX,1000));     
			fprintf(f30,"%lf ", print_receive_var(i,j,E_CX,1001));     
		}
	}

	// rCX
	if (Num_Types>E_rCX) for(i = 0; i < cell_sizes[E_rCX].x; ++i){

		// fprintf(f11,"%lf ", print_receive(i,cell_sizes[E_rCX].y/2,E_rCX));

		// y variable
		//fprintf(f11,"%lf ", print_receive_var(i,cell_sizes[E_rCX].y/2,E_rCX,0));

		fprintf(f30,"%lf ", print_receive_var(i,cell_sizes[E_rCX].y/2,E_rCX,1000));     
		fprintf(f30,"%lf ", print_receive_var(i,cell_sizes[E_rCX].y/2,E_rCX,1001));     

	}

	fprintf(f11,"\n");
	fprintf(f30,"\n");

	//fflush(f11);



//	for(i = 0; i < cell_sizes[E_CX2].x; ++i){
//		for(j = 0; j < cell_sizes[E_CX2].y; ++j){
//
//			fprintf(cx2,"%lf ", print_receive(i,j,E_CX));
//
//			// y variable
//			// fprintf(cx2,"%lf ", print_receive_var(i,j,E_CX,0));
//
//			// Print currents
//			// fprintf(f30,"%lf ", print_receive_var(i,j,E_CX,1000));     
//			// fprintf(f30,"%lf ", print_receive_var(i,j,E_CX,1001));     
//		}
//	}
//	fprintf(cx2,"\n");

	fprintf(cx3,"%lf ", t);	

	for(i = 0; i < cell_sizes[E_CX3].x; ++i){
		for(j = 0; j < cell_sizes[E_CX3].y; ++j){

			fprintf(cx3,"%lf ", print_receive(i,j,E_CX3));

			// y variable
			// fprintf(cx3,"%lf ", print_receive_var(i,j,E_CX,0));

			// Print currents
			fprintf(cx3Currents,"%lf ", print_receive_var(i,j,E_CX3,1000));     
			fprintf(cx3Currents,"%lf ", print_receive_var(i,j,E_CX3,1001));     
		}
	}
	fprintf(cx3,"\n");
	fprintf(cx3Currents,"\n");

	/////// TODO: add printing of currents for the right hemisphere

	fprintf(cx4,"%lf ", t);

	for(i = 0; i < cell_sizes[E_CX4].x; ++i){
		for(j = 0; j < cell_sizes[E_CX4].y; ++j){

			fprintf(cx4,"%lf ", print_receive(i,j,E_CX4));

			// y variable
			// fprintf(cx4,"%lf ", print_receive_var(i,j,E_CX,0));

			// Print currents
			fprintf(cx4Currents,"%lf ", print_receive_var(i,j,E_CX4,1000));     
			fprintf(cx4Currents,"%lf ", print_receive_var(i,j,E_CX4,1001));     
		}
	}
	fprintf(cx4,"\n");
	fprintf(cx4Currents,"\n");

	fprintf(cx5a,"%lf ", t);

	for(i = 0; i < cell_sizes[E_CX5a].x; ++i){
		for(j = 0; j < cell_sizes[E_CX5a].y; ++j){

			fprintf(cx5a,"%lf ", print_receive(i,j,E_CX5a));

			// y variable
			// fprintf(cx5a,"%lf ", print_receive_var(i,j,E_CX,0));

			// Print currents
			fprintf(cx5aCurrents,"%lf ", print_receive_var(i,j,E_CX5a,1000));     
			fprintf(cx5aCurrents,"%lf ", print_receive_var(i,j,E_CX5a,1001));     
		}
	}
	fprintf(cx5a,"\n");
	fprintf(cx5aCurrents,"\n");

	fprintf(cx5b,"%lf ", t);

	for(i = 0; i < cell_sizes[E_CX5b].x; ++i){
		for(j = 0; j < cell_sizes[E_CX5b].y; ++j){

			fprintf(cx5b,"%lf ", print_receive(i,j,E_CX5b));

			// y variable
			// fprintf(cx5b,"%lf ", print_receive_var(i,j,E_CX,0));

			// Print currents
			fprintf(cx5bCurrents,"%lf ", print_receive_var(i,j,E_CX5b,1000));     
			fprintf(cx5bCurrents,"%lf ", print_receive_var(i,j,E_CX5b,1001));     
		}
	}
	fprintf(cx5b,"\n");
	fprintf(cx5bCurrents,"\n");

//	for(i = 0; i < cell_sizes[E_IN2].x; ++i){
//		for(j = 0; j < cell_sizes[E_IN2].y; ++j){
//
//			fprintf(in2,"%lf ", print_receive(i,j,E_CX));
//
//			// y variable
//			// fprintf(in2,"%lf ", print_receive_var(i,j,E_CX,0));
//
//			// Print currents
//			// fprintf(f30,"%lf ", print_receive_var(i,j,E_CX,1000));     
//			// fprintf(f30,"%lf ", print_receive_var(i,j,E_CX,1001));     
//		}
//	}
//	fprintf(in2,"\n");

	fprintf(in3,"%lf ", t);

	for(i = 0; i < cell_sizes[E_IN3].x; ++i){
		for(j = 0; j < cell_sizes[E_IN3].y; ++j){

			fprintf(in3,"%lf ", print_receive(i,j,E_IN3));

			// y variable
			// fprintf(in3,"%lf ", print_receive_var(i,j,E_CX,0));

			// Print currents
//			fprintf(in3Currents,"%lf ", print_receive_var(i,j,E_IN3,1000));     
//			fprintf(in3Currents,"%lf ", print_receive_var(i,j,E_IN3,1001));     
		}
	}
	fprintf(in3,"\n");

	fprintf(in4,"%lf ", t);

	for(i = 0; i < cell_sizes[E_IN4].x; ++i){
		for(j = 0; j < cell_sizes[E_IN4].y; ++j){

			fprintf(in4,"%lf ", print_receive(i,j,E_IN4));

			// y variable
			// fprintf(in4,"%lf ", print_receive_var(i,j,E_CX,0));

			// Print currents
//			fprintf(in4Currents,"%lf ", print_receive_var(i,j,E_IN4,1000));     
//			fprintf(in4Currents,"%lf ", print_receive_var(i,j,E_IN4,1001));     
		}
	}
	fprintf(in4,"\n");

	fprintf(in5a,"%lf ", t);

	for(i = 0; i < cell_sizes[E_IN5a].x; ++i){
		for(j = 0; j < cell_sizes[E_IN5a].y; ++j){

			fprintf(in5a,"%lf ", print_receive(i,j,E_IN5a));

			// y variable
			// fprintf(in5a,"%lf ", print_receive_var(i,j,E_CX,0));

			// Print currents
//			fprintf(in5aCurrents,"%lf ", print_receive_var(i,j,E_IN5a,1000));     
//			fprintf(in5aCurrents,"%lf ", print_receive_var(i,j,E_IN5a,1001));     
		}
	}
	fprintf(in5a,"\n");

	fprintf(in5b,"%lf ", t);

	for(i = 0; i < cell_sizes[E_IN5b].x; ++i){
		for(j = 0; j < cell_sizes[E_IN5b].y; ++j){

			fprintf(in5b,"%lf ", print_receive(i,j,E_IN5b));

			// y variable
			// fprintf(in5b,"%lf ", print_receive_var(i,j,E_CX,0));

			// Print currents
//			fprintf(in5bCurrents,"%lf ", print_receive_var(i,j,E_IN5b,1000));     
//			fprintf(in5bCurrents,"%lf ", print_receive_var(i,j,E_IN5b,1001));     
		}
	}
	fprintf(in5b,"\n");

	// **** Print for CXa cells **** 

	
// 	//print time
//	fprintf(f15,"%lf ", t);
//
//	//print individual cell data
//	for(i = 0; i < cell_sizes[E_CXa].x; ++i){
//    for(j = 0; j < cell_sizes[E_CXa].y; ++j){
//
//      //fprintf(f15,"%lf ", print_receive(i,cell_sizes[E_CXa].y/2,E_CXa));
//      fprintf(f15,"%lf ", print_receive(i,j,E_CXa));
//
//      // y variable
//      //fprintf(f15,"%lf ", print_receive_var(i,cell_sizes[E_CXa].y/2,E_CXa,0));
//      //fprintf(f15,"%lf ", print_receive_var(i,j,E_CXa,0));
//
//      // Print currents
//      fprintf(f31,"%lf ", print_receive_var(i,cell_sizes[E_CXa].y/2,E_CXa,1000));     
//      fprintf(f31,"%lf ", print_receive_var(i,cell_sizes[E_CXa].y/2,E_CXa,1001));     
//    }
//	}
//
//	if (Num_Types>E_rCXa) for(i = 0; i < cell_sizes[E_rCXa].x; ++i){
//
//		//fprintf(f15,"%lf ", print_receive(i,cell_sizes[E_rCXa].y/2,E_rCXa));
//
//		// y variable
//		//fprintf(f15,"%lf ", print_receive_var(i,cell_sizes[E_rCXa].y/2,E_rCXa,0));
//
//		// Print currents
//		fprintf(f31,"%lf ", print_receive_var(i,cell_sizes[E_rCXa].y/2,E_rCXa,1000));     
//		fprintf(f31,"%lf ", print_receive_var(i,cell_sizes[E_rCXa].y/2,E_rCXa,1001));     
//
//	}
//
//	fprintf(f15,"\n");
//	fprintf(f31,"\n");

	// **** Print for CX6 cells **** 

  fprintf(f25,"%lf ", t); // TB IN layers
  for(i = 0; i < cell_sizes[E_CX6].x; ++i){
    for(j = 0; j < cell_sizes[E_CX6].y; ++j){

      fprintf(f25,"%lf ", print_receive(i,j,E_CX6));  // TB IN layers

      // y variable
      //fprintf(f25,"%lf ", print_receive_var(i,j,E_CX6,0));

      // Print currents
      fprintf(f32,"%lf ", print_receive_var(i,j,E_CX6,1000));     
      fprintf(f32,"%lf ", print_receive_var(i,j,E_CX6,1001));     
    }
  }
	/*fprintf(f25,"%lf ", t); // TB IN layers
	//if (Num_Types>E_rCX6) for(i = 0; i < cell_sizes[E_rCX6].x; ++i){
	if (Num_Types>E_rCX6) for(i = 0; i < cell_sizes[E_rCX6].x; ++i){

		fprintf(f25,"%lf ", print_receive(i,cell_sizes[E_rCX6].y/2,E_rCX6));  // TB IN layers

		// y variable
		fprintf(f25,"%lf ", print_receive_var(i,cell_sizes[E_rCX6].y/2,E_rCX6,0));

		// Print currents
		fprintf(f32,"%lf ", print_receive_var(i,cell_sizes[E_rCX6].y/2,E_rCX6,1000));     
		fprintf(f32,"%lf ", print_receive_var(i,cell_sizes[E_rCX6].y/2,E_rCX6,1001));     

	}*/
	fprintf(f25,"\n"); // TB IN layers
	fprintf(f32,"\n");



	fprintf(f13,"%lf ", t);
	for(i = 0; i < cell_sizes[E_IN].x; ++i){
		fprintf(f13,"%lf ", print_receive(i,cell_sizes[E_IN].y/2,E_IN));
	}
	fprintf(f13,"\n");

//	fprintf(f19,"%lf ", t); // TB IN layers
//	for(i = 0; i < cell_sizes[E_INa].x; ++i){
//		fprintf(f19,"%lf ", print_receive(i,cell_sizes[E_INa].y/2,E_INa));  // TB IN layers
//	}
//	fprintf(f19,"\n"); // TB IN layers

	/*fprintf(f25,"%lf ", t); // TB IN layers
	for(i = 0; i < cell_sizes[E_CX6].x; ++i){
		fprintf(f25,"%lf ", print_receive(i,cell_sizes[E_CX6].y/2,E_CX6));  // TB IN layers
	}
	fprintf(f25,"%lf ", t); // TB IN layers
	for(i = 0; i < cell_sizes[E_rCX6].x; ++i){
		fprintf(f25,"%lf ", print_receive(i,cell_sizes[E_rCX6].y/2,E_rCX6));  // TB IN layers
	}
	fprintf(f25,"\n"); // TB IN layers
*/

	fprintf(f27,"%lf ", t); // TB IN layers
	for(i = 0; i < cell_sizes[E_IN6].x; ++i){
		fprintf(f27,"%lf ", print_receive(i,cell_sizes[E_IN6].y/2,E_IN6));  // TB IN layers
	}
	fprintf(f27,"\n"); // TB IN layers






//       j=0;
//       for(i = 0; i < cell_sizes[E_CX].x; ++i){
//            print_receive_gsynapse_multLayer(i, j, E_CX, E_CX, E_AMPAMap_D1, fconduct_cx2_cx2);
//       }
//       fprintf(fconduct_cx2_cx2,"\n");
//
//       j=0;
//       for(i = 0; i < cell_sizes[E_CX].x; ++i){
//            print_receive_gsynapse_multLayer(i, j, E_CX3, E_CX, E_AMPAMap_D1, fconduct_cx3_cx2);
//       }
//       fprintf(fconduct_cx3_cx2,"\n");
//
//       j=0;
//       for(i = 0; i < cell_sizes[E_CX].x; ++i){
//            print_receive_gsynapse_multLayer(i, j, E_CX4, E_CX, E_AMPAMap_D1, fconduct_cx4_cx2);
//       }
//       fprintf(fconduct_cx4_cx2,"\n");
//
//       j=0;
//       for(i = 0; i < cell_sizes[E_CX].x; ++i){
//            print_receive_gsynapse_multLayer(i, j, E_CX5a, E_CX, E_AMPAMap_D1, fconduct_cx5a_cx2);
//       }
//       fprintf(fconduct_cx5a_cx2,"\n");
//
//       j=0;
//       for(i = 0; i < cell_sizes[E_CX].x; ++i){
//            print_receive_gsynapse_multLayer(i, j, E_CX5b, E_CX, E_AMPAMap_D1, fconduct_cx5b_cx2);
//       }
//       fprintf(fconduct_cx5b_cx2,"\n");
//
//       j=0;
//       for(i = 0; i < cell_sizes[E_CX].x; ++i){
//            print_receive_gsynapse_multLayer(i, j, E_CX6, E_CX, E_AMPAMap_D1, fconduct_cx6_cx2);
//       }
//       fprintf(fconduct_cx6_cx2,"\n");
//
//
//       j=0;
//       for(i = 0; i < cell_sizes[E_CX3].x; ++i){
//            print_receive_gsynapse_multLayer(i, j, E_CX, E_CX3, E_AMPAMap_D1, fconduct_cx2_cx3);
//       }
//       fprintf(fconduct_cx2_cx3,"\n");
//
//       j=0;
//       for(i = 0; i < cell_sizes[E_CX3].x; ++i){
//            print_receive_gsynapse_multLayer(i, j, E_CX3, E_CX3, E_AMPAMap_D1, fconduct_cx3_cx3);
//       }
//       fprintf(fconduct_cx3_cx3,"\n");
//
//       j=0;
//       for(i = 0; i < cell_sizes[E_CX3].x; ++i){
//            print_receive_gsynapse_multLayer(i, j, E_CX4, E_CX3, E_AMPAMap_D1, fconduct_cx4_cx3);
//       }
//       fprintf(fconduct_cx4_cx3,"\n");
//
//       j=0;
//       for(i = 0; i < cell_sizes[E_CX3].x; ++i){
//            print_receive_gsynapse_multLayer(i, j, E_CX5a, E_CX3, E_AMPAMap_D1, fconduct_cx5a_cx3);
//       }
//       fprintf(fconduct_cx5a_cx3,"\n");
//
//       j=0;
//       for(i = 0; i < cell_sizes[E_CX3].x; ++i){
//            print_receive_gsynapse_multLayer(i, j, E_CX5b, E_CX3, E_AMPAMap_D1, fconduct_cx5b_cx3);
//       }
//       fprintf(fconduct_cx5b_cx3,"\n");
//
//       j=0;
//       for(i = 0; i < cell_sizes[E_CX3].x; ++i){
//            print_receive_gsynapse_multLayer(i, j, E_CX6, E_CX3, E_AMPAMap_D1, fconduct_cx6_cx3);
//       }
//       fprintf(fconduct_cx6_cx3,"\n");

}

void open_files(string output_location,FILE **field_file, int num_field_layers){
	printf("open files\n");
	int i = 0;
	for(i=0; i<num_field_layers; i++){
		stringstream ss;
		ss << i;
		field_file[i] = fopen((output_location+"field_file_" + ss.str()  ).c_str(), "w");
	}

	f2 = fopen((output_location+"dat").c_str(), "w");
	f6 = fopen((output_location+"graf_re").c_str(), "w");
	if (!(f7=fopen((output_location+"time_re").c_str(), "w"))){
		printf("probably out put folder doesn't exist\n");
		exit(1); 
	}
	f29 = fopen((output_location+"time_rea").c_str(), "w");

	//f7 = fopen((output_location+"time_re").c_str(), "w");
	f8 = fopen((output_location+"graf_tc").c_str(), "w");
	f9 = fopen((output_location+"time_tc").c_str(), "w");
	f10 = fopen((output_location+"graf_cx").c_str(), "w");
	f11 = fopen((output_location+"time_cx").c_str(), "w");
	f12 = fopen((output_location+"graf_in").c_str(), "w");
	f13 = fopen((output_location+"time_in").c_str(), "w");
	//TANYA MODIF layers 
	f14 = fopen((output_location+"graf_cxa").c_str(), "w");
	f15 = fopen((output_location+"time_cxa").c_str(), "w");
	f16 = fopen((output_location+"graf_tca").c_str(), "w");
	f17 = fopen((output_location+"time_tca").c_str(), "w");
	f18 = fopen((output_location+"graf_ina").c_str(), "w"); // TB IN layers
	f19 = fopen((output_location+"time_ina").c_str(), "w"); // TB IN layers
	//END TANYA MODIF layers 
	f20 = fopen((output_location+"time_G_AMPA0_CX_CX").c_str(), "w");
	f21 = fopen((output_location+"time_G_AMPA0_CXa_CX").c_str(), "w");
	f22 = fopen((output_location+"time_G_AMPA0_CX_CXa").c_str(), "w"); // TB IN layers
	f23 = fopen((output_location+"time_G_AMPA0_CXa_CXa").c_str(), "w"); // TB IN layers

	f24 = fopen((output_location+"graf_cx6").c_str(), "w"); // TB IN layers
	f25 = fopen((output_location+"time_cx6").c_str(), "w"); // TB IN layers
	f26 = fopen((output_location+"graf_in6").c_str(), "w"); // TB IN layers
	f27 = fopen((output_location+"time_in6").c_str(), "w"); // TB IN layers

	f28 = fopen((output_location+"cx_cx_g_ampa0").c_str(), "w");

	f30 = fopen((output_location+"e_i_currents_cx").c_str(), "w");
//	f31 = fopen((output_location+"e_i_currents_cxa").c_str(), "w");
	f32 = fopen((output_location+"e_i_currents_cx6").c_str(), "w");

 // cx2 = fopen((output_location+"time_cx2").c_str(), "w");
  cx3 = fopen((output_location+"time_cx3").c_str(), "w");
  cx4 = fopen((output_location+"time_cx4").c_str(), "w");
  cx5a = fopen((output_location+"time_cx5a").c_str(), "w");
  cx5b = fopen((output_location+"time_cx5b").c_str(), "w");

 // in2 = fopen((output_location+"time_in2").c_str(), "w");
  in3 = fopen((output_location+"time_in3").c_str(), "w");
  in4 = fopen((output_location+"time_in4").c_str(), "w");
  in5a = fopen((output_location+"time_in5a").c_str(), "w");
  in5b = fopen((output_location+"time_in5b").c_str(), "w");

  cx3Currents = fopen((output_location+"e_i_currents_cx3").c_str(), "w");
  cx4Currents = fopen((output_location+"e_i_currents_cx4").c_str(), "w");
  cx5aCurrents = fopen((output_location+"e_i_currents_cx5a").c_str(), "w");
  cx5bCurrents = fopen((output_location+"e_i_currents_cx5b").c_str(), "w");

//  in3Currents = fopen((output_location+"e_i_currents_in3").c_str(), "w");
//  in4Currents = fopen((output_location+"e_i_currents_in4").c_str(), "w");
//  in5aCurrents = fopen((output_location+"e_i_currents_in5a").c_str(), "w");
//  in5bCurrents = fopen((output_location+"e_i_currents_in5b").c_str(), "w");

  fSpikeTime = fopen((output_location+"spike_time").c_str(),"w");
  fconn = fopen((output_location+"conn_index_cx_cx_G").c_str(), "w");

//  fconn_cx2_cx2 = fopen((output_location+"conn_index_cx2_cx2_G").c_str(), "w");
//  fconn_cx3_cx2 = fopen((output_location+"conn_index_cx3_cx2_G").c_str(), "w");
//  fconn_cx4_cx2 = fopen((output_location+"conn_index_cx4_cx2_G").c_str(), "w");
//  fconn_cx5a_cx2 = fopen((output_location+"conn_index_cx5a_cx2_G").c_str(), "w");
//  fconn_cx5b_cx2 = fopen((output_location+"conn_index_cx5b_cx2_G").c_str(), "w");
//  fconn_cx6_cx2 = fopen((output_location+"conn_index_cx6_cx2_G").c_str(), "w");
//
//  fconn_cx2_cx3 = fopen((output_location+"conn_index_cx2_cx3_G").c_str(), "w");
//  fconn_cx3_cx3 = fopen((output_location+"conn_index_cx3_cx3_G").c_str(), "w");
//  fconn_cx4_cx3 = fopen((output_location+"conn_index_cx4_cx3_G").c_str(), "w");
//  fconn_cx5a_cx3 = fopen((output_location+"conn_index_cx5a_cx3_G").c_str(), "w");
//  fconn_cx5b_cx3 = fopen((output_location+"conn_index_cx5b_cx3_G").c_str(), "w");
//  fconn_cx6_cx3 = fopen((output_location+"conn_index_cx6_cx3_G").c_str(), "w");
//
//
//  fconduct_cx2_cx2 =  fopen((output_location+"conduct_cx2_cx2_G").c_str(), "w");
//  fconduct_cx3_cx2 =  fopen((output_location+"conduct_cx3_cx2_G").c_str(), "w");
//  fconduct_cx4_cx2 =  fopen((output_location+"conduct_cx4_cx2_G").c_str(), "w");
//  fconduct_cx5a_cx2 = fopen((output_location+"conduct_cx5a_cx2_G").c_str(), "w");
//  fconduct_cx5b_cx2 = fopen((output_location+"conduct_cx5b_cx2_G").c_str(), "w");
//  fconduct_cx6_cx2 =  fopen((output_location+"conduct_cx6_cx2_G").c_str(), "w");
//  fconduct_cx2_cx3 =  fopen((output_location+"conduct_cx2_cx3_G").c_str(), "w");
//  fconduct_cx3_cx3 =  fopen((output_location+"conduct_cx3_cx3_G").c_str(), "w");
//  fconduct_cx4_cx3 =  fopen((output_location+"conduct_cx4_cx3_G").c_str(), "w");
//  fconduct_cx5a_cx3 = fopen((output_location+"conduct_cx5a_cx3_G").c_str(), "w");
//  fconduct_cx5b_cx3 = fopen((output_location+"conduct_cx5b_cx3_G").c_str(), "w");
//  fconduct_cx6_cx3 =  fopen((output_location+"conduct_cx6_cx3_G").c_str(), "w");



	printf("files open for write\n");
}


//TODO convert this to something sensible
void close_files(FILE **field_file, int num_field_layers){

	int i = 0;
	for(i=0; i<num_field_layers; i++){
		fclose(field_file[i]);
	}
	
        if(fSpikeTime!=NULL){
           fclose(fSpikeTime);
        } 

	if(f2!=NULL){
		fclose(f2);
	}
	if(f6!=NULL){
		fclose(f6);
	}
	if(f7!=NULL){
		fclose(f7);
	}
	if(f8!=NULL){
		fclose(f8);
	}
	if(f9!=NULL){
		fclose(f9);
	}
	if(f10!=NULL){
		fclose(f10);
	}
	if(f11!=NULL){
		fclose(f11);
	}
	if(f12!=NULL){
		fclose(f12);
	}
	if(f13!=NULL){
		fclose(f13);
	}
	if(f14!=NULL){
		fclose(f14);
	}
	if(f15!=NULL){
		fclose(f15);
	}
	if(f16!=NULL){
		fclose(f16);
	}
	if(f17!=NULL){
		fclose(f17);
	}
	if(f18!=NULL){
		fclose(f18);
	}  // TB IN layers
	if(f19!=NULL){
		fclose(f19);
	}
	if(f20!=NULL){
		fclose(f20);
	}
	if(f21!=NULL){
		fclose(f21);
	}
	if(f22!=NULL){
		fclose(f22);
	}
	if(f23!=NULL){
		fclose(f23);
	}
	if(f24!=NULL){
		fclose(f24);
	}
	if(f25!=NULL){
		fclose(f25);
	}
	if(f26!=NULL){
		fclose(f26);
	}
	if(f27!=NULL){
		fclose(f27);
	}
	if(f28!=NULL){
		fclose(f28);
	}
	if(f29!=NULL){
		fclose(f29);
	}
	if(f30!=NULL){
		fclose(f30);
	}
//	if(f31!=NULL){
//		fclose(f31);
//	}
	if(f32!=NULL){
		fclose(f32);
	}
	if(fconn!=NULL){
                fclose(fconn);
	}
//
//	if(fconn_cx2_cx2 !=NULL){
//                fclose(fconn_cx2_cx2);
//        }
//	if(fconn_cx3_cx2 !=NULL){
//                fclose(fconn_cx3_cx2);
//        }
//	if(fconn_cx4_cx2 !=NULL){
//                fclose(fconn_cx4_cx2);
//        }
//	if(fconn_cx5a_cx2 !=NULL){
//                fclose(fconn_cx5a_cx2);
//        }
//	if(fconn_cx5b_cx2 !=NULL){
//                fclose(fconn_cx5b_cx2);
//        }
//	if(fconn_cx6_cx2 !=NULL){
//                fclose(fconn_cx6_cx2);
//        }
//
//
//	if(fconn_cx2_cx3 !=NULL){
//                fclose(fconn_cx2_cx3);
//        }
//	if(fconn_cx3_cx3 !=NULL){
//                fclose(fconn_cx3_cx3);
//        }
//	if(fconn_cx4_cx3 !=NULL){
//                fclose(fconn_cx4_cx3);
//        }
//	if(fconn_cx5a_cx3 !=NULL){
//                fclose(fconn_cx5a_cx3);
//        }
//	if(fconn_cx5b_cx3 !=NULL){
//                fclose(fconn_cx5b_cx3);
//        }
//	if(fconn_cx6_cx3 !=NULL){
//                fclose(fconn_cx6_cx3);
//        }
//
//	if(fconduct_cx2_cx2 !=NULL){
//                fclose(fconduct_cx2_cx2);
//        }
//	if(fconduct_cx3_cx2 !=NULL){
//                fclose(fconduct_cx3_cx2);
//        }
//	if(fconduct_cx4_cx2 !=NULL){
//                fclose(fconduct_cx4_cx2);
//        }
//	if(fconduct_cx5a_cx2 !=NULL){
//                fclose(fconduct_cx5a_cx2);
//        }
//	if(fconduct_cx5b_cx2 !=NULL){
//                fclose(fconduct_cx5b_cx2);
//        }
//	if(fconduct_cx6_cx2 !=NULL){
//                fclose(fconduct_cx6_cx2);
//        }
//
//
//	if(fconduct_cx2_cx3 !=NULL){
//                fclose(fconduct_cx2_cx3);
//        }
//	if(fconduct_cx3_cx3 !=NULL){
//                fclose(fconduct_cx3_cx3);
//        }
//	if(fconduct_cx4_cx3 !=NULL){
//                fclose(fconduct_cx4_cx3);
//        }
//	if(fconduct_cx5a_cx3 !=NULL){
//                fclose(fconduct_cx5a_cx3);
//        }
//	if(fconduct_cx5b_cx3 !=NULL){
//                fclose(fconduct_cx5b_cx3);
//        }
//	if(fconduct_cx6_cx3 !=NULL){
//                fclose(fconduct_cx6_cx3);
//        }



//	if(cx2!=NULL){
//		fclose(cx2);
//	}
	if(cx3!=NULL){
		fclose(cx3);
	}
	if(cx4!=NULL){
		fclose(cx4);
	}
	if(cx5a!=NULL){
		fclose(cx5a);
	}
	if(cx5b!=NULL){
		fclose(cx5b);
	}
//	if(in2!=NULL){
//		fclose(in2);
//	}
	if(in3!=NULL){
		fclose(in3);
	}
	if(in4!=NULL){
		fclose(in4);
	}
	if(in5a!=NULL){
		fclose(in5a);
	}
	if(in5b!=NULL){
		fclose(in5b);
	}
	if(cx3Currents!=NULL){
		fclose(cx3Currents);
	}
	if(cx4Currents!=NULL){
		fclose(cx4Currents);
	}
	if(cx5aCurrents!=NULL){
		fclose(cx5aCurrents);
	}
	if(cx5bCurrents!=NULL){
		fclose(cx5bCurrents);
	}
//	if(in3Currents!=NULL){
//		fclose(in3Currents);
//	}
//	if(in4Currents!=NULL){
//		fclose(in4Currents);
//	}
//	if(in5aCurrents!=NULL){
//		fclose(in5aCurrents);
//	}
//	if(in5bCurrents!=NULL){
//		fclose(in5bCurrents);
//	}
}

//returns MB
void print_used_mem(){
	FILE* proc = fopen("/proc/self/status", "r");
	string res;
	char line[128];

	while (fgets(line, 128, proc)){
		if (!strncmp(line, "VmRSS:", 6)){
			printf("Resident RAM: %s",line);
			break;
		}
	}
	fclose(proc);
} 
