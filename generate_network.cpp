#include <string>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <random>
#include <map>

using namespace std;
ofstream * ConnSummaryFile;

typedef unsigned char BYTE;
//default_random_engine generator;



#define MAXTYPES 33
class CellTypes {
 public:
  struct size { int x, y; string name;};
  int n; 		//number of all types we have
  size cell[MAXTYPES]; 	//info for each type
  bool right_hemisphere(int cellt){//is the cell in left or right hemisphere (left is default for networks without hemispheric difference)
    return (cell[cellt].name.find("r")==0);
  }//right_hemisphere
  int get_type(const string &t){  //map name into index
    for (int i=0; i<n; i++) if (cell[i].name==t) return i;
    assert(0); //no such neuron type exists
  }//get_type
  void dump(){ for (int i=0; i<n; i++) {dump(i);cerr<<"\n";} }
  void dump(int i){
     assert(i<n);
     (*ConnSummaryFile)<<cell[i].name<<" "<<cell[i].x<<" "<<cell[i].y;
  }
  void output(){ //fixed output format used by simulation
    cout<<n<<"\n";
    for (int i=0; i<n; i++)
     cout<<cell[i].x<<" "<<cell[i].y<<"\n";
  }
  CellTypes():n(0),MRI_ratio(1){};
  double MRI_ratio; //mainly for debugging purposes when we want simulate smaller network from hardcoded MRI input data, we shrink number of neurons by MRU_ratio factor
  void MRI_reduce(){ for (int i=0; i<n; i++) cell[i].x*=MRI_ratio; }
};//CellTypes

CellTypes CT;

//connection types (not particular connections between neurons)
class Connections{
 public:
  //connection types between layers
  struct conn { 
    int from, to; string synapse;
      double radius_min, radius_max, cs; //radius; column size
      double probab, probab_oc; // probability within  and outside of column
    string distribution; 
    double strength,
    mini_strength,   //minimal strength of miniature synaptic potentials
    mini_freq;       //frequency of spontaneous minis (miniature synaptic potentials)
    int range;       //long/short range; this value has different meaning for generated (1==short range) and MRI network
   };//conn
  int n; //number of connections
  conn C[MAXTYPES*MAXTYPES]; //connections
  
  int find(int from,int to,string synapse,int range_type){
    for (int i=0; i<n; i++)
      if (C[i].from == from && C[i].to == to && C[i].synapse == synapse && range_type == C[i].range)
        return i;
   cerr<<"Missing connection "<<from<<"->"<<to<<" "<<synapse<<", range "<<range_type<<"\n";
   //exit(0); //connection type not entered in config file
   return 100;
  }//find 
  void dump(){ for (int i=0; i<n; i++) dump(i); }//dump
  void dump(int i){
     assert(i<n);
     (*ConnSummaryFile)<<C[i].from<<" "<<C[i].to<<" "<<C[i].synapse<<" "<<
     C[i].radius_min<<" "<<C[i].radius_max<<" "<<C[i].probab<<"\n";
  }
  int stats[MAXTYPES*MAXTYPES]; //connection statistics for debug purposes
  void dump_stats(){
    for (int c=0; c<n; c++) {
      (*ConnSummaryFile)<<c+1<<" of "<<n<<" , edges: "<<stats[c]<<" \ttype: (";
      CT.dump(C[c].from);(*ConnSummaryFile)<<") -> (";CT.dump(C[c].to);(*ConnSummaryFile)<<") ";
      dump(c);
    }
  }//dump stats
  Connections():n(0){ for (int i=0;i<MAXTYPES*MAXTYPES;i++) stats[i]=0;};
};//Connections

Connections CN;

//C++11 has std::round
double round(double x){
  return (x<0.0) ? ceil(x-0.5f) : floor(x+0.5f);
}

#define sqr(x) ((x)*(x))
//euclidean distance; we could use manhattan for 2D networks actually
//check carefuly comments in generate_connections routine in case
//you want to change this metric.
inline double dist(double x1, double y1, double x2, double y2){
  return sqrt(sqr(x1-x2)+sqr(y1-y2));
}

inline int is_column(int x1, int y1, int x2, int y2, int CS){
    if ((x1/CS == x2/CS) && (y1/CS==y2/CS)) return 1;
    else return 0;
}


//TODO use mersenne twister instead of this weak guy
inline double rand01(){return ((double) rand() / (RAND_MAX));}

string get_input_line(FILE *f,string &line){
  char buf[4096];
  string tok1;

  while(1) {
    if (!fgets(buf,sizeof(buf),f)) {return line="";}; //EOF
    line=buf;
    if (line.empty()) continue;
    stringstream parser(buf);
    if (!(parser>>tok1)) continue;    //whitespace line
    if (tok1[0]=='#') continue;       //just comment
    break;
  }
  return tok1;
}

bool both_hemispheres_detected=false; //in case we load ucsd data for both hemispheres. signalized in config via existence of "r" prefix in cell type name

int parse_config(char const *file){
  string line,tok;
  double x,y,ratio_x,ratio_y;
  FILE *f=fopen(file,"r");
  assert(f);

  // Types
  tok=get_input_line(f,line);
  assert(tok=="[Types]");
  //ratios
  tok=get_input_line(f,line);
  stringstream parser(line);
  assert(parser >> ratio_x); assert(parser >> ratio_y);
  //types enumeration
  tok=get_input_line(f,line);
  while (tok!="[Connections]"){
    assert(CT.n<MAXTYPES);
    parser.str(line);
    assert(parser >> CT.cell[CT.n].name >> x >> y);

    if (!both_hemispheres_detected && CT.cell[CT.n].name.find("r")==0) { 
      cerr<<"Configuration for both-hemisphere connections detected (rXXX cell type found).\n";
      both_hemispheres_detected=true;
      //only x will be used, because MRI data are 1 dimensionally indexed, and it must be set _after_ loading MRI data
      if (ratio_x!=1) {
        cerr<<"Ratio will be used in different way because both hemispheres detected\n";
	CT.MRI_ratio=ratio_x;
	ratio_x=1;
      }
    }

    CT.cell[CT.n].x=x; CT.cell[CT.n].y=y;
    CT.n++;

    tok=get_input_line(f,line);
  }//types

  //ratios after loop because MRI scenario
  for (int i=0; i<CT.n;i++) {
    CT.cell[CT.n].x=round(CT.cell[CT.n].x*ratio_x); CT.cell[CT.n].y=round(CT.cell[CT.n].y*ratio_y);
  }

  //Short range connections
  string from, to, syn, dist; double rad,cs,p,poc,str,mf,ms;
  assert(tok=="[Connections]");
  tok=get_input_line(f,line);
  while (tok!="[Longrange]"){
    assert(CN.n<MAXTYPES*MAXTYPES);
    parser.str(line);
    assert(parser >> from >> to >> syn >> rad >> cs >> p >> poc >> dist >> str >> ms >> mf);

    Connections::conn &c=CN.C[CN.n];	 //current connection
    c.from=CT.get_type(from); c.to=CT.get_type(to); c.synapse=syn;
    c.radius_min=0; c.radius_max=rad; c.cs=cs; c.probab=p; c.probab_oc=poc;
    c.distribution=dist; c.strength=str; c.mini_strength=ms; c.mini_freq=mf;
    c.range=1;
    CN.n++; 
    
    tok=get_input_line(f,line);
  }//connections

  //Long range connections
  double rad_min,rad_max;
  assert(tok=="[Longrange]");
  tok=get_input_line(f,line);
  while (!feof(f)){
    assert(CN.n<MAXTYPES*MAXTYPES);
    parser.str(line);
    assert(parser >> from >> to >> syn >> rad_min >> rad_max >> p >> dist >> str >> ms >> mf);

    Connections::conn &c=CN.C[CN.n];	 //current connection
    c.from=CT.get_type(from); c.to=CT.get_type(to); c.synapse=syn;
    c.radius_min=rad_min; c.radius_max=rad_max; c.probab=p;
    c.range=0;
    CN.n++; 

    //TODO we should change rad_min for small (1D) networks to be bigger or use specific function
    
    tok=get_input_line(f,line);
  }//connections

  fclose(f);
}

//representing single synaptic connection when generating network
class edge {
 public:
  int from_t, from_x, from_y, to_t, to_x, to_y; //out==from, in==to
  double delay, weightFactor;
  Connections::conn *syntype; 
  edge(int ft,int fx,int fy,int tt,int tx,int ty,Connections::conn *s, double d):
    from_t(ft),from_x(fx),from_y(fy),to_t(tt),to_x(tx),to_y(ty),syntype(s), delay(d),weightFactor(1){};
  edge(int ft,int fx,int fy,int tt,int tx,int ty,Connections::conn *s, double d, double weightF):
    from_t(ft),from_x(fx),from_y(fy),to_t(tt),to_x(tx),to_y(ty),syntype(s), delay(d),weightFactor(weightF){};
  static int ct,cx,cy; //cache for last to_ printed -> we want to print "In: " line sometimes
  int output(){
     int new_neuron=0; //did we enter new neurons while printing connections; for total neurons count.
     if (to_t!=ct || to_x!=cx || to_y!=cy) {
       cout<<"In: "<<to_t<<" "<<to_x<<" "<<to_y<<"\n";
       ct=to_t; cx=to_x; cy=to_y;
       new_neuron=1;
     }
     if (weightFactor>0){
     	cout<<from_t<<" "<<from_x<<" "<<from_y<<" "<<syntype->synapse<<" "<<syntype->strength*weightFactor
        	<<" "<<syntype->mini_strength<<" "<<syntype->mini_freq<<" "<<syntype->range<<" "<<delay<<"\n";
     }
     return new_neuron;
  }
  void dump(){
     cerr<<from_t<<" "<<from_x<<" "<<from_y<< " -> " <<to_t<<" "<<to_x<<" "<<to_y<<" "<<syntype->synapse<<"\n";
  }
  inline static bool cmp(const edge &a, const edge &b){ //for sorting
     if (a.to_t!=b.to_t) return (a.to_t<b.to_t);
     if (a.to_y!=b.to_y) return (a.to_y<b.to_y);
     if (a.to_x!=b.to_x) return (a.to_x<b.to_x);
     if (a.from_t!=b.from_t) return (a.from_t<b.from_t);
     if (a.from_y!=b.from_y) return (a.from_y<b.from_y);
     if (a.from_x!=b.from_x) return (a.from_x<b.from_x);
     return false;
  }//cmp
};//edge
int edge::ct=-1; int edge::cx=-1; int edge::cy=-1;
class Edges{
 public:
  vector<edge> edges;
  //this is rough version, performance can be substantially improved in case we still need it by
  //1) testing/enumerating neurons only in radius - DONE
  //2) Using 3D fixed array for sorting results instead of 1D vector
  void generate_connections(bool useEuclideanDelays){ //generates geometry topology
    cerr<<"No input file given, generating network on our own\n";

    for (int c=0; c<CN.n; c++){ //layer connections
      int FX = CT.cell[CN.C[c].from].x; int FY = CT.cell[CN.C[c].from].y; //input layer dimension
      int TX = CT.cell[CN.C[c].to].x; int TY = CT.cell[CN.C[c].to].y; //output layer dimension
      int ec=0;	//edges counter for statistics
      for (int fx=0; fx<FX; fx++)
        for (int fy=0; fy<FY; fy++) {
          //get scaled position of target neuron in TO layer
          int sc_x = (double)TX/FX * fx; int sc_y = (double)TY/FY * fy;
          //restrict the local neighbourhood where we look for possible connections
	  //by given radius. this speedup trick can fail in case someone switch to
	  //non-eucledian metric.
          //[max(sc_x-radius,0)..min(sc_x+radius,TX)]
	  //[max(sc_y-radius,0)..min(sc_y+radius,TY)]
          int radius = CN.C[c].radius_max + 1;

          for (int tx=max(sc_x-radius,0); tx<min(sc_x+radius,TX); tx++) //particular connections between neurons
            for (int ty=max(sc_y-radius,0); ty<min(sc_y+radius,TY); ty++) { 

             if (sc_x==tx && sc_y==ty && CT.cell[CN.C[c].from].name==CT.cell[CN.C[c].to].name) continue; //never connect neuron to itself
             //distance is measured between 'From' neuron projected into 'TO' layer
             double d=dist(sc_x,sc_y, tx,ty);

             if ( d>=CN.C[c].radius_min && d<=CN.C[c].radius_max ){
                 double delay = 0;
                 if(useEuclideanDelays){
                   delay = dist((double)fx, (double)fy, (double)tx, (double) ty) * 0.5;// 2 millisecond delay for every index the neurons are apart
                   //delay = d; //this distance doesn't correspond to he indexes used to create the synapse (fx, fy, tx, ty) but it may be better to use in other circumstances
                   //           // fx != sc_x and fy != sc_y in all cases 
                 }
                  // check if neurons are in same or different column
                 if (CN.C[c].range==1) {
                     if (is_column(sc_x,sc_y,tx,ty,CN.C[c].cs)) {
                       if (rand01() <= CN.C[c].probab) 
                         edges.push_back(edge(CN.C[c].from,fx,fy,CN.C[c].to,tx,ty,&CN.C[c],delay)),ec++;
                    } else {  //out of column
                        if (rand01() <= CN.C[c].probab_oc) 
                          edges.push_back(edge(CN.C[c].from,fx,fy,CN.C[c].to,tx,ty,&CN.C[c],delay)),ec++;
                      } 
                 } else {  //longrange
                     if (rand01() <= CN.C[c].probab) 
                       edges.push_back(edge(CN.C[c].from,fx,fy,CN.C[c].to,tx,ty,&CN.C[c],delay)),ec++;
                   }
             }//radius
            }//ty
          }//fy

      //dbg
      (*ConnSummaryFile)<<c+1<<" of "<<CN.n<<" , edges: "<<ec<<" \ttype: (";
      CT.dump(CN.C[c].from);(*ConnSummaryFile)<<") -> (";CT.dump(CN.C[c].to);(*ConnSummaryFile)<<") ";
      CN.dump(c);
    }//c
  }//generate_connections
  void generate_3D_connections(bool);
  void generate_MultiLayer_connections(bool);
  void write_connections(bool  right_hemisphere=false){
   cerr<<"Writing output file\n";
   //header
   if (!right_hemisphere) CT.output();
   int te=0,tn=0;	//total edges, total neurons
   vector<edge>::iterator it=edges.begin();
   vector<edge>::iterator end=edges.end();
   for (;it!=end;it++){
    te++; 
    tn += (*it).output(); 
    //(*it).dump();
   }
   cerr<<"Total counts (hemisphere "<<right_hemisphere<<"): "<<tn<<" neurons, "<<te<<" edges.\n";
  }//write_connections

  //assign connection attributes for each edge in network (this can be dynamic) 
  //this routine is separated from topology generation since we load some topology-only networks 
  void apply_synaptic_types(){

   cerr<<"Sorting connections\n";
   sort(edges.begin(),edges.end(),edge::cmp);

   cerr<<"Adjusting synaptic types\n";
   vector<edge>::iterator it=edges.begin();
   vector<edge>::iterator end=edges.end();
   for (;it!=end;it++){
    edge &c=(*it);
    if (c.syntype->distribution=="fixed") 
      continue; //syntype->already contains values we will use in output

    //TODO decide what to do with gaussian or normal distribution
//    if (c.syntype->distribution=="uniform") {}// ..
//    if (c.syntype->distribution=="gauss") {}// ..
   }    
 }//apply

 //dump all output connection for a given neuron (we usually write down the opposite)
 void dump_from_neuron_edges(int type, int x, int y, int to_type){
   vector<edge>::iterator it=edges.begin();
   vector<edge>::iterator end=edges.end();
   for (;it!=end;it++){
    edge &c=(*it);
    if (c.from_t==type && c.from_x == x && c.from_y == y && (to_type==-1 || c.to_t == to_type))
      cerr<<type<<" "<<x<<","<<y<<" -> "<<c.to_t<<" "<<c.to_x<<","<<c.to_y<<"\n";
   }    
 }//dump
 void dump_from_neuron_edges(int type, int x, int y){ dump_from_neuron_edges(type, x, y, -1); }
 void load_MRI_network(const char *file);
 void generate_3D_intrahemispheric_connections();
 int connect_local_neighbours(int c /*rule*/,int fx, int tx,bool neighbourhood_hemisphere);
 bool check(int type,int neuron); //hack for shrinking MRI data based on ratio, we allow adding edges which neurons pass this check  
};//Edges

Edges CE;

//Tables for mapping specific synapses to MAP types for the hybrid model
//Tables are derived from the loading code, but the whole thing can be explained
//much simpler at the end: CX & IN are Map based neurons and _any_ synapse on them
//must be converted into MAP type; lazy to simplify this now.
//The mechanism needs to be checked in case short X long range connection differ
//in synapses used.
//All to all connections, e.g. CX* -> CX*
const string rule_all[][4] =
                       { {"CX","CX",  "AMPA_D2",  "AMPAMap_D1"},
                         {"CX","CX",  "NMDA_D1",  "NMDAMap_D1"},
                         {"TC","CX",  "AMPA_D1",  "AMPAMap_D1"},
                         {"TC","IN",  "AMPA_D1",  "AMPAMap_D1"} };

//Paired connections, eg. CX/a/6 -> IN/a/6 (not CX6->INa)
const string rule_paired[][4]=
                       { {"CX","IN",  "AMPA_D2",  "AMPAMap_D1"},
                         {"CX","IN",  "NMDA_D1",  "NMDAMap_D1"},
                         {"IN","CX","GABA_A_D2","GABAMap_A_D1"},
                         {"CX","CX",  "NMDA_D1",  "NMDAMap_D1"}}; //these should be actually part of the previous

const int rules_all=sizeof(rule_all)/sizeof(rule_all[0]);
const int rules_paired=sizeof(rule_paired)/sizeof(rule_paired[0]);

//does the rule and string match in the beginning in both cases
int match(string const from, string const rule_from, string const to, string const rule_to){
  return !from.compare(0,rule_from.length(),rule_from) && !to.compare(0,rule_to.length(),rule_to);
}

//does the rule and string match in the beginning in both cases, does the suffix matches (e.g. 6 in  CX6)
int match_paired(string const from, string const rule_from, string const to, string const rule_to){ 
   assert(from.length() && to.length());
   if (from==rule_from && to==rule_to) //without suffix
     return true;
   return !from.compare(0,rule_from.length(),rule_from) &&
          !to.compare(0,rule_to.length(),rule_to) &&
	  from[from.length()-1] == to[to.length()-1];          
}

//load network from the file (fixed structure) we
//get from ucsd folks approximating MRI topology
void Edges::load_MRI_network(const char *file){
  cerr<<"Loading network\n";

  //cleanup x,y header information, we will fill it with what we read from here.
  for (int ii=0; ii<CT.n; ii++) CT.cell[ii].x=CT.cell[ii].y=0;

  FILE *f=fopen(file,"r"); assert(f);
  //now we go through the connections file and fill in the info we need
  while(1){
    int type,x_loc,y_loc;
    if(fscanf(f,"Incoming connections for cell: type: %d x: %d y: %d\n",&type,&x_loc,&y_loc) != 3)
      break;

    CT.cell[type].x = max(x_loc,CT.cell[type].x); CT.cell[type].y = max(y_loc,CT.cell[type].y);

    //second pass going through all the connections to the cell
    while(1){
      char *syn_type = new char[256];
      double strength = 0.0;
      double mini_s = 0.0;
      double mini_fre = 0.0;
      int mri_range = 0;
      int in_type,in_x_loc,in_y_loc;

      if(fscanf(f,"type: %d x: %d y: %d Syntype: %s range: %d \n",&in_type,&in_x_loc,&in_y_loc,syn_type,&mri_range) != 5)
	break;
      
      //transform synaptic type from given file to hybrid model with MAP synapse type currently used.
      string syn_type_lit(syn_type);
      string from(CT.cell[in_type].name);
      string to(CT.cell[type].name);

      //Short/long range connection conversion:
      //3 used for inter-hemispheric "short" range connections, 0 for intrahemispheric long range connections, 1 forshort range
      //3 is currently ignored, we don't have connections defined
      if (mri_range == 3) continue;
      int short_range = mri_range; //except for 3 the same meaning as in normal generation

      //transform synapse names according to the tables defined above
      for (int i=0; i<rules_all; i++){
        if (rule_all[i][2]!=syn_type_lit) continue;
	if (match(from,rule_all[i][0],to,rule_all[i][1]))
	  syn_type_lit=rule_all[i][3];
      }
      for (int i=0; i<rules_paired; i++){
        if (rule_paired[i][2]!=syn_type_lit) continue;
	if (match_paired(from,rule_paired[i][0],to,rule_paired[i][1]))
	  syn_type_lit=rule_paired[i][3];
      }

      int conn_type = CN.find(in_type,type,syn_type_lit,short_range);
      edges.push_back(edge(in_type,in_x_loc,in_y_loc,type,x_loc,y_loc,&CN.C[conn_type],0));
      CN.stats[conn_type]++; //for stats printing
    }//while 2nd pass
   }//while file

  CN.dump_stats();
  fclose(f);
}//load_MRI_network

#define MRI_POINTS 12500
float Distance3D[2][MRI_POINTS][MRI_POINTS]; //distance lookup table given for MRI data; [0][][] for left hemispher
float maxDistance = 0;// max distance between two neurons
float max_TC_Distance = 0;
float maxDelay = 40;// max delay between farthest neurons. all delays will be scaled according to this
float maxDelay_TCc = 40; // core TC 
float maxDelay_TCm = 40; // matrix TC
float avgDist = 0;
float numCons = 0;
bool useDelays = false;

double getScaledDelay(double distance){
  double dist = useDelays ? distance : 0; //if you should use delays or not use delays
  return (dist / maxDistance) * maxDelay;
}

double getScaled_TC_Delay(double distance, double maxDelay){
  double dist = useDelays ? distance : 0; //if you should use delays or not use delays
  return (dist / max_TC_Distance) * maxDelay;
}

void load_3D_distances(const char *file, bool right_hemisphere){
  cerr<<"Loading 3D network distances, hemisphere: "<<right_hemisphere<<"\n";

  for (int i=0; i<MRI_POINTS; i++)
    for (int ii=0; ii<MRI_POINTS; ii++)
      Distance3D[right_hemisphere][i][ii]=-1;

  FILE *f=fopen(file,"r"); assert(f);

  //go through the distance matrix
  int i=0;
  while(1){
    int from,to;
    float dist;
    if(fscanf(f,"%d %d %f\n",&from,&to,&dist)!=3)
      break;
    from--; to--; //MATLAB vs C indexing
    assert(from<MRI_POINTS); assert(to<MRI_POINTS);
    Distance3D[right_hemisphere][from][to]=dist;
    //avgDist += dist;
    //numCons = numCons + 1;

    if(dist > maxDistance)// get max distance for scaling later
      maxDistance = dist;
  }//while

  fclose(f);
}

#define MAX_SUBNET_POINTS 12500	//The actual number of smaller net is to be read out from network.cfg (TC cells)
float Subnet[2][MAX_SUBNET_POINTS]; //distance lookup table given for MRI data, [0][] for left hemisphere

//load subset of CX neurons to be flagged as indexes for TC,RE & IN
void load_3D_subnet(const char *file, bool right_hemisphere){
  cerr<<"Loading 3D sub network mapping, hemisphere: "<<right_hemisphere<<"\n";

  for (int i=0; i<MAX_SUBNET_POINTS; i++)
    Subnet[right_hemisphere][i]=-1;

  FILE *f=fopen(file,"r"); assert(f);

  //go through the distance matrix
  int i=0;
  while(1){
    int id,map;
    if(fscanf(f,"%d %d\n",&id,&map)!=2)
      break;
    id--; map--; //MATLAB vs C indexing
//    cerr<<id<<" "<<map<<" " << CT.cell[CT.get_type("TC")].x<<"\n";
    if (right_hemisphere) assert(id< CT.cell[CT.get_type("rTC")].x);
    if (!right_hemisphere) assert(id< CT.cell[CT.get_type("TC")].x);
    assert(id<MRI_POINTS); assert(map<MAX_SUBNET_POINTS); 

    Subnet[right_hemisphere][id]=map;
  }//while

  fclose(f);
}

int map_LH_RH [MAX_SUBNET_POINTS]; //mapping between LH-RH neurons
int map_LH_RH_small [MAX_SUBNET_POINTS]; //mapping between LH-RH neurons, just the subset
void load_3D_LH_RH_correspondence(const char *file, int Map[]){
  cerr<<"Loading 3D LH-RH interhemispheric mapping "<<file<<"\n";

  for (int i=0; i<MRI_POINTS; i++)
    Map[i]=-1;

  FILE *f=fopen(file,"r"); assert(f);

  int i=0;
  while(1){
    int id,map;
    if(fscanf(f,"%d %d\n",&id,&map)!=2)
      break;
    id--; map--; //MATLAB vs C indexing
    assert(id<MRI_POINTS);// assert(map<MAX_SUBNET_POINTS); 

    Map[id]=map;

    // printf("Mapping id: %d to %d \n",id,map);
  }//while

  fclose(f);
}

inline double dist3D(int id_from, int id_to, bool right_hemisphere){
  // printf("Distance3D - %lf from %d to %d \n",Distance3D[right_hemisphere][id_from][id_to],id_from,id_to);

 assert(Distance3D[right_hemisphere][id_from][id_to]!=-1);
 return Distance3D[right_hemisphere][id_from][id_to];
}

#define NUM_TC_POINTS 12500
float TC_Dist[NUM_TC_POINTS];

void load_thalCort_dist(const char *file){
  cerr << "load_thalCort_dist" << endl;

  for (int ii=0; ii<NUM_TC_POINTS; ii++)
    TC_Dist[ii]=-1;

  FILE *f=fopen(file,"r"); assert(f);

  int i=0;
  while(1){
    float dist;
    if(fscanf(f,"%f\n",&dist)!=1)
      break;    

    if (dist>=0){
    	TC_Dist[i] = dist;
    }

    if(dist > max_TC_Distance){// get max distance for scaling later
      max_TC_Distance = dist;
    }

    i++;
  }//while

  fclose(f);

}


#define MRI_POINTS_FULL 20500 // 20500 //  12500
#define NUM_PARAM 2

//float Dist_Prob[MRI_POINTS_FULL][MRI_POINTS_FULL][NUM_PARAM]; //distance lookup table given for MRI data; [0][][] for left hemispher
float ***Dist_Prob;

void load_3D_dist_prob(const char *file){
cerr << "load_3D_dist_prob" << endl;

// for (int i=0; i<MRI_POINTS_FULL; i++){
//    for (int ii=0; ii<MRI_POINTS_FULL; ii++){
//      for (int kk=0; kk<NUM_PARAM; kk++){
//      // cerr << i << " " << ii << " " << kk << endl;
//	      Dist_Prob[i][ii][kk]=-1;
//      }
//    }
//  }

  FILE *f=fopen(file,"r"); assert(f);

  int i=0;
  while(1){
    int from,to;
    float dist, conn_prob;
    if(fscanf(f,"%d %d %f %f\n",&from,&to,&dist,&conn_prob)!=4)
      break;
    from--; to--; //MATLAB vs C indexing
    assert(from<MRI_POINTS); assert(to<MRI_POINTS);
    Dist_Prob[from][to][0]=dist;
    Dist_Prob[from][to][1]=conn_prob;

    if(dist > maxDistance)// get max distance for scaling later
      maxDistance = dist;
  }//while

  fclose(f);
}

float **Dist_intraThalamus;
void load_intraThalamus_dist(const char *file){
	cerr << "load_intraThalamus_dist" << endl;
	FILE *f = fopen(file, "r"); assert(f);
	while (1){
		int from, to;
		float dist;
		if(fscanf(f, "%d %d %f\n", &from,&to,&dist)!=3)
			break;
		from--; to--;
		assert(from<MRI_POINTS); assert(to<MRI_POINTS);
		Dist_intraThalamus[from][to]=dist;
	}
	fclose(f);
}


inline double dist_only(int id_from, int id_to){
  // printf("Distance3D - %lf from %d to %d \n",Distance3D[right_hemisphere][id_from][id_to],id_from,id_to);

// assert(Dist_Prob[id_from][id_to][0]!=-1);
 return Dist_Prob[id_from][id_to][0];
}

inline double prob_only(int id_from, int id_to){
  // printf("Distance3D - %lf from %d to %d \n",Distance3D[right_hemisphere][id_from][id_to],id_from,id_to);

// assert(Dist_Prob[id_from][id_to][1]!=-1);
 return Dist_Prob[id_from][id_to][1];
}

inline double intraThalamicDistOnSphere(int id_from, int id_to){
  return Dist_intraThalamus[id_from][id_to];
}
//////////////// read files //////////////////////

#define NUM_LAYERS 6
//int weight_factor[MRI_POINTS_FULL][MRI_POINTS_FULL][NUM_LAYERS][NUM_LAYERS]; //distance lookup table given for MRI data; [0][][] for left hemispher
//int ****weight_factor;
float ****weight_factor;
    
void load_weight_factors(const char *file){

cerr << "load_weight_factors" << endl;
  
//  for (int i=0; i<MRI_POINTS_FULL; i++){
//    for (int ii=0; ii<MRI_POINTS_FULL; ii++){
//      for (int kk=0; kk<NUM_LAYERS; kk++){
//        for (int ll=0; ll<NUM_LAYERS; ll++){
//              weight_factor[i][ii][kk][ll]=-1;
//	}
//      }
//    }
//  }
    
  FILE *f=fopen(file,"r"); assert(f);

  int i=0;
  while(1){
    int fromNeuron,toNeuron, fromLayer,toLayer; //,factor;    
    float factor;
    if(fscanf(f,"%d %d %d %d %f\n",&fromNeuron,&toNeuron,&fromLayer,&toLayer,&factor)!=5)
      break; 
    fromNeuron--;toNeuron--;fromLayer--;toLayer--; //MATLAB vs C indexing
    assert(fromNeuron<MRI_POINTS_FULL); assert(toNeuron<MRI_POINTS_FULL);
    weight_factor[fromNeuron][toNeuron][fromLayer][toLayer]=factor;
      
  }//while
        
  fclose(f);
}


// float weight_factor_INPY[MRI_POINTS_FULL][MRI_POINTS_FULL];
float **weight_factor_INPY;

void load_weightFactors_INPY(const char *file){
cerr << "load_weightFactors_INPY" << endl;
//  for (int i=0; i<MRI_POINTS_FULL; i++){
//    for (int ii=0; ii<MRI_POINTS_FULL; ii++){
//              weight_factor_INPY[i][ii] =-1;
//    }
//  }

  FILE *f=fopen(file,"r"); assert(f);

  int i=0;
  while(1){
    int fromNeuron,toNeuron; //,factor;
    float factor;
    if(fscanf(f,"%d %d %f \n",&fromNeuron,&toNeuron,&factor)!=3)
      break;
    fromNeuron--;toNeuron--; //MATLAB vs C indexing
    assert(fromNeuron<MRI_POINTS_FULL); assert(toNeuron<MRI_POINTS_FULL);
    weight_factor_INPY[fromNeuron][toNeuron]=factor;

  }//while

  fclose(f);
}


/////////////// end read fiels //////////////////

//////////////// Read binary files /////////////////////

map<std::string, int> mapLayerNameId_IN ;
map<std::string, int> mapLayerNameId_PY ;

map<std::string, int> fillDict_PY(map<string, int> m){
        m.insert(pair<string,int>("CX",0));
        m.insert(pair<string,int>("CX3",1));
        m.insert(pair<string,int>("CX4",2));
        m.insert(pair<string,int>("CX5a",3));
        m.insert(pair<string,int>("CX5b",4));
        m.insert(pair<string,int>("CX6",5));

        m.insert(pair<string,int>("rCX",0));
        m.insert(pair<string,int>("rCX3",1));
        m.insert(pair<string,int>("rCX4",2));
        m.insert(pair<string,int>("rCX5a",3));
        m.insert(pair<string,int>("rCX5b",4));
        m.insert(pair<string,int>("rCX6",5));

        return m;
}

map<std::string, int> fillDict_IN(map<string, int> m){

        m.insert(pair<string,int>("IN",0));
        m.insert(pair<string,int>("IN3",1));
        m.insert(pair<string,int>("IN4",2));
        m.insert(pair<string,int>("IN5a",3));
        m.insert(pair<string,int>("IN5b",4));
        m.insert(pair<string,int>("IN6",5));

        m.insert(pair<string,int>("rIN",0));
        m.insert(pair<string,int>("rIN3",1));
        m.insert(pair<string,int>("rIN4",2));
        m.insert(pair<string,int>("rIN5a",3));
        m.insert(pair<string,int>("rIN5b",4));
        m.insert(pair<string,int>("rIN6",5));

        return m;
}

void intReshape( vector<BYTE> fileData, int dim1, int dim2, int dim3, int dim4 ){

//	for (int i=0; i<MRI_POINTS_FULL; i++){
//	    for (int ii=0; ii<MRI_POINTS_FULL; ii++){
//		 for (int kk=0; kk<NUM_LAYERS; kk++){
//		     for (int ll=0; ll<NUM_LAYERS; ll++){
//             		  weight_factor[i][ii][kk][ll]=-1;
//       		     }
//                 }
//            }
//        }

        for (int ll=0; ll< dim4; ll++){
                for (int kk=0; kk< dim3; kk++){			                        
                        for (int jj=0; jj< dim2; jj++){

//				cerr << "---------------" << endl;
//                        	cerr <<ll << " " << kk << " " << jj<<  endl;

                                for (int ii=0; ii< dim1; ii++){
                                       // weight_factor[ii][jj][kk][ll] = fileData[ ll*(dim3*dim2*dim1) + kk*(dim2*dim1) + jj*dim1 + ii ];
                                  // try
                                  // {
                                  //      weight_factor[ii][jj][kk][ll] = fileData[ ll*(dim3*dim2*dim1) + kk*(dim2*dim1) + jj*dim1 + ii ];

                                  // }
                                  // catch (std::exception& e)
                                  // {
                                  //     cerr << "seg error hit" << endl;
                                  //     cerr << ii << " " << jj << " " << kk << " " << ll << endl;

                                  // }
                                      
				    //  cerr << "---------------" << endl;
                                    //  cerr << ii << " " << jj << " " << kk << " " << ll << endl;
                                    //  cerr << 1 << endl;
                                    //  int value = fileData[ ll*(dim3*dim2*dim1) + kk*(dim2*dim1) + jj*dim1 + ii ];
                                    //  cerr << 2 << endl;
                                    //   weight_factor[ii][jj][kk][ll] = value;
                                    //  cerr << 3 << endl;
				    
				    weight_factor[ii][jj][kk][ll] = fileData[ ll*(dim3*dim2*dim1) + kk*(dim2*dim1) + jj*dim1 + ii ];
                               }

                       }
               }
        }
}

void load_weightFactors_binary(const char *filename){

cerr << "load_weightFactors_binary" << endl;
    streampos fileSize;
    ifstream file(filename, ios::binary);

    // get its size:
    file.seekg(0, ios::end);
    fileSize = file.tellg();
    file.seekg(0, ios::beg);

    // read the data:
    vector<BYTE> fileData(fileSize);
    file.read((char*) &fileData[0], fileSize);

    intReshape( fileData, 20484, 20484, 6, 6);
}

////////////// End - read binary files /////////////////


void Edges::generate_3D_connections(bool right_hemisphere){ //generates connection using external distance matrix
    cerr<<"Generating connection file from pre-set 3D MRI data, hemisphere: "<<right_hemisphere<<"\n";

    int small = CT.cell[CT.get_type("TC")].x; int large = CT.cell[CT.get_type("CX")].x;

    if (both_hemispheres_detected && right_hemisphere)
      small = CT.cell[CT.get_type("rTC")].x, large = CT.cell[CT.get_type("rCX")].x;
      
    cerr<<"Full set(CX): "<<large<<", Subnet(TC,RE,IN...): "<<small<<"\n";

    for (int c=0; c<CN.n; c++){ //layer connections
      if (right_hemisphere == !CT.right_hemisphere(CN.C[c].from)) continue;  //must be intrahemispheric
      if (right_hemisphere == !CT.right_hemisphere(CN.C[c].to)) continue;
      int FX = CT.cell[CN.C[c].from].x; int FY = CT.cell[CN.C[c].from].y; //input layer dimension
      int TX = CT.cell[CN.C[c].to].x; int TY = CT.cell[CN.C[c].to].y; //output layer dimension
      assert(FY==1); assert(TY==1); //MRI data are stored as 1D indexes

      int ec=0;	//edges counter for statistics

      for (int fxi=0, fx= (FX==large) ? 0:Subnet[right_hemisphere][0]; (FX==large) ? fx<FX: Subnet[right_hemisphere][fxi]!=-1; (FX==large)? fx++: fx=Subnet[right_hemisphere][++fxi]) {
          //scaling is actually hidden in Subnet array, identity in large case, map in small case
	  int fy = 0;
          int sc_x = fx; int sc_y = fy;

          for (int txi=0, tx= (TX==large) ? 0:Subnet[right_hemisphere][0]; (TX==large) ? tx<TX: Subnet[right_hemisphere][txi]!=-1; (TX==large)? tx++: tx=Subnet[right_hemisphere][++txi]) { //particular connections between neurons
	    int ty=0;
	    if (sc_x==tx && sc_y==ty && CT.cell[CN.C[c].from].name==CT.cell[CN.C[c].to].name) continue; //never connect neuron to itself
            //distance is measured between 'From' neuron projected into 'TO' layer
             // printf("sc_x: %d, tx: %d, right_hem:%d ", sc_x,tx,right_hemisphere);
             double d=dist3D(sc_x, tx, right_hemisphere);
             double delayTime = getScaledDelay(d);
             if ( d>=CN.C[c].radius_min && d<=CN.C[c].radius_max ){
                // column is no more checked
               if (CN.C[c].range==1) {
                 if (rand01() <= CN.C[c].probab_oc) 
                   check(CN.C[c].from,(FX==large)? fx:fxi) && check(CN.C[c].to,(TX==large)? tx:txi) && (edges.push_back(edge(CN.C[c].from,(FX==large)? fx:fxi,fy,CN.C[c].to,(TX==large)? tx:txi,ty,&CN.C[c],delayTime)),ec++);
               } else {  //longrange
                   if (rand01() <= CN.C[c].probab) 
                       check(CN.C[c].from, (FX==large)? fx:fxi) && check(CN.C[c].to,(TX==large)? tx:txi) && (edges.push_back(edge(CN.C[c].from,(FX==large)? fx:fxi,fy,CN.C[c].to,(TX==large)? tx:txi,ty,&CN.C[c],delayTime)),ec++);
                 } //longrange
             }//radius
            }//tx
          }//fx

      //dbg
      (*ConnSummaryFile)<<c+1<<" of "<<CN.n<<" , edges: "<<ec<<" \ttype: (";
      CT.dump(CN.C[c].from);(*ConnSummaryFile)<<") -> (";CT.dump(CN.C[c].to);(*ConnSummaryFile)<<") ";
      CN.dump(c);
    }//c
  }//generate_3Dconnections



void Edges::generate_MultiLayer_connections(bool right_hemisphere){
cerr<<"Generating connection file from pre-set 3D MRI data, new MultiLayer connections, hemisphere: "<<right_hemisphere<<"\n";

// // uncomment to use a non-deterministic seed:
//    std::random_device rd;
//    std::mt19937 gen(rd());
std::mt19937 gen(1937);

    int small = CT.cell[CT.get_type("TC")].x; int large = CT.cell[CT.get_type("CX")].x;
    int indOffset =0;
    if (both_hemispheres_detected && right_hemisphere){
      small = CT.cell[CT.get_type("rTC")].x, large = CT.cell[CT.get_type("rCX")].x;
      indOffset = CT.cell[CT.get_type("CX")].x;
    }


    cerr<<"Full set(CX): "<<large<<", Subnet(TC,RE,IN...): "<<small<<"\n";

    for (int c=0; c<CN.n; c++){ //layer connections
      if (right_hemisphere == !CT.right_hemisphere(CN.C[c].from)) continue;  //must be intrahemispheric
      if (right_hemisphere == !CT.right_hemisphere(CN.C[c].to)) continue;
      int FX = CT.cell[CN.C[c].from].x; int FY = CT.cell[CN.C[c].from].y; //input layer dimension
      int TX = CT.cell[CN.C[c].to].x; int TY = CT.cell[CN.C[c].to].y; //output layer dimension
      assert(FY==1); assert(TY==1); //MRI data are stored as 1D indexes


      string fromCellType = CT.cell[CN.C[c].from].name;
      string toCellType = CT.cell[CN.C[c].to].name;

      int fromLayer = -1;
      int toLayer = -1;
      int from_IN = 0;
      int to_IN=0;
      if( mapLayerNameId_PY.find(CT.cell[CN.C[c].from].name) != mapLayerNameId_PY.end() ) {
          fromLayer = mapLayerNameId_PY.find(CT.cell[CN.C[c].from].name)->second;
      }else if( mapLayerNameId_IN.find(CT.cell[CN.C[c].from].name) != mapLayerNameId_IN.end() ) {
          from_IN=1;
      }
            
      if( mapLayerNameId_PY.find(CT.cell[CN.C[c].to].name) != mapLayerNameId_PY.end() ) {
          toLayer=mapLayerNameId_PY.find(CT.cell[CN.C[c].to].name)->second;
      }else if( mapLayerNameId_IN.find(CT.cell[CN.C[c].to].name) != mapLayerNameId_IN.end() ) {
          to_IN=1;
      }

      assert(fromLayer<6); assert(toLayer<6);

      int ec=0; //edges counter for statistics

      for (int fxi=0, fx= (FX==large) ? 0:Subnet[right_hemisphere][0]; (FX==large) ? fx<FX: Subnet[right_hemisphere][fxi]!=-1; (FX==large)? fx++: fx=Subnet[right_hemisphere][++fxi]) {
          //scaling is actually hidden in Subnet array, identity in large case, map in small case
          int fy = 0;
          int sc_x = fx; int sc_y = fy;

          for (int txi=0, tx= (TX==large) ? 0:Subnet[right_hemisphere][0]; (TX==large) ? tx<TX: Subnet[right_hemisphere][txi]!=-1; (TX==large)? tx++: tx=Subnet[right_hemisphere][++txi]) { //particular connections between neurons
            int ty=0;
            if (sc_x==tx && sc_y==ty && CT.cell[CN.C[c].from].name==CT.cell[CN.C[c].to].name) continue; //never connect neuron to itself
            //distance is measured between 'From' neuron projected into 'TO' layer
             // printf("sc_x: %d, tx: %d, right_hem:%d ", sc_x,tx,right_hemisphere);
             double d=dist_only(sc_x + indOffset, tx + indOffset); 
	     // if(d<0) continue;
             double delayTime=0.0; // = getScaledDelay(d);

		// use spherical Thalamus
		size_t foundFromTC = fromCellType.find("TC");
                size_t foundFromRE = fromCellType.find("RE");
                size_t foundToTC = toCellType.find("TC");
                size_t foundToRE = toCellType.find("RE");

		if( ((foundFromTC != string::npos)||(foundFromRE != string::npos)) && ((foundToTC != string::npos)||(foundToRE != string::npos))   ){
			d = intraThalamicDistOnSphere(sc_x, tx); // no indOffsetfor the right hemi within thalamus	
		}

	        if(d<0) continue;

		// use dealyes only on CX->CX and CX->IN connections
		size_t foundFromCX = fromCellType.find("CX");
		size_t foundToCX = toCellType.find("CX");
		size_t foundToIN = toCellType.find("IN");
		if ( ! (( foundFromCX != string::npos) && (( foundToCX!= string::npos ) || ( foundToIN != string::npos )))){ 
			delayTime=0.0;
		}else{
 			delayTime = getScaledDelay(d);
		}


	     	double currEdgeProb = prob_only(sc_x + indOffset, tx + indOffset);
		double clampProb = (CN.C[c].probab_oc*currEdgeProb <1) ? CN.C[c].probab_oc*currEdgeProb : 1;
                
		if (fromLayer>-1 && toLayer >-1){ // PY->PY connections            	 
//			if(txi==0){cerr<<"PY-PY connection --- 1.1 \n";}
/////			std::default_random_engine generator;
			std::bernoulli_distribution distribution(clampProb);
/////		   if (distribution(generator)){
		   if (distribution(gen)){
//			if(txi==0){cerr<<"generated PY-PY connection --- 1.2 \n";}
		       double currEdgeWeightFactor = weight_factor[sc_x + indOffset][tx + indOffset][fromLayer][toLayer];
		      //  if (d<= 0.0001 ){
			    //  currEdgeWeightFactor *= 5; // local weights scaled up to test their influence
                      //  }
          // if (d>0.005 ){
			    // currEdgeWeightFactor = currEdgeWeightFactor / 2; // all weights scaled down
          //  }

          if (d>0.01 ){
			      currEdgeWeightFactor = currEdgeWeightFactor / 5; // long weights scaled down
          }

          if ( d > CN.C[c].radius_max ){
		    	     currEdgeWeightFactor = 0;
		       }
		       if (currEdgeWeightFactor >0){
          	          check(CN.C[c].from,(FX==large)? fx:fxi) && check(CN.C[c].to,(TX==large)? tx:txi) && (edges.push_back(edge(CN.C[c].from,(FX==large)? fx:fxi,fy,CN.C[c].to,(TX==large)? tx:txi,ty,&CN.C[c],delayTime,currEdgeWeightFactor)),ec++);
		       }
		   }
		}else if(fromLayer>-1 && to_IN>0){ // PY->IN connections		
//			std::default_random_engine generator;
			std::bernoulli_distribution distribution(clampProb);
//		   if (distribution(generator)){
		   if (distribution(gen)){
		       double currEdgeWeightFactor = weight_factor[sc_x + indOffset][tx + indOffset][fromLayer][fromLayer]; // Since PY->IN are only within the same layer, use 'fromLayer' twice
		       if (currEdgeWeightFactor >0){
        	          check(CN.C[c].from,(FX==large)? fx:fxi) && check(CN.C[c].to,(TX==large)? tx:txi) && (edges.push_back(edge(CN.C[c].from,(FX==large)? fx:fxi,fy,CN.C[c].to,(TX==large)? tx:txi,ty,&CN.C[c],delayTime,currEdgeWeightFactor)),ec++);
		       }
		   }				
		}else if(from_IN>0 && toLayer>-1){ // IN->PY connections
		   if ( d>=CN.C[c].radius_min && d<=CN.C[c].radius_max ){
                       if (CN.C[c].range==1) {
                         if (rand01() <= CN.C[c].probab_oc){
			    double currEdgeWeightFactor = weight_factor_INPY[sc_x + indOffset][tx + indOffset];
			    if (currEdgeWeightFactor >0){
                               check(CN.C[c].from,(FX==large)? fx:fxi) && check(CN.C[c].to,(TX==large)? tx:txi) && (edges.push_back(edge(CN.C[c].from,(FX==large)? fx:fxi,fy,CN.C[c].to,(TX==large)? tx:txi,ty,&CN.C[c],delayTime,currEdgeWeightFactor)),ec++);
			    }
			}
                       } 
                     }//radius		
		}else{ // all other types of connections, all are based on the disc model
	             if ( d>=CN.C[c].radius_min && d<=CN.C[c].radius_max ){
			     

      			    // use TC_Dist only on thalamocortical connections (CX,TC,IN,RE)
  		    	    size_t foundFromCX = fromCellType.find("CX");
  		    	    size_t foundFromIN = fromCellType.find("IN");
  		    	    size_t foundToCX = toCellType.find("CX");
  		    	    size_t foundToIN = toCellType.find("IN");

  		    	    size_t foundFromTC = fromCellType.find("TC");
  		    	    size_t foundFromRE = fromCellType.find("RE");
  		    	    size_t foundToTC = toCellType.find("TC");  		    	    
			    size_t foundToRE = toCellType.find("RE");

			    // different maxDealy on Core and Matrix connections
			    // be careful with condition since there is CX5a which is not matrix (TCa)
  		    	    size_t foundFromMatrix = fromCellType.find("a");
  		    	    size_t foundToMatrix = toCellType.find("a");  		    	    

			    if ( ((foundFromCX != string::npos) || (foundFromIN != string::npos) ) && ( (foundToTC != string::npos) || (foundToRE != string::npos)  )  ){
			    	if ( (foundToMatrix != string::npos) ){				    
					// TODO for both hemis use: sc_x + indOffset, tx + indOffset
 				  delayTime = getScaled_TC_Delay( TC_Dist[sc_x], maxDelay_TCm);
				}else{				
 				  delayTime = getScaled_TC_Delay( TC_Dist[sc_x], maxDelay_TCc);
				}
			    }

			    if ( ((foundToCX != string::npos) || (foundToIN != string::npos) ) && ( (foundFromTC != string::npos) || (foundFromRE != string::npos)  )  ){
			    	if ( (foundFromMatrix != string::npos) ){				    
					// TODO for both hemis use: sc_x + indOffset, tx + indOffset
 				  delayTime = getScaled_TC_Delay( TC_Dist[tx], maxDelay_TCm);
				}else{				
 				  delayTime = getScaled_TC_Delay( TC_Dist[tx], maxDelay_TCc);
				}
			    }

			    if (delayTime<0){
				    delayTime=0;
			    }

			     double currEdgeWeightFactor=1;
			     if ((CT.cell[CN.C[c].from].name == "TC")&& (CT.cell[CN.C[c].to].name=="CX4")){
				     currEdgeWeightFactor=3;			     
			     } 
	                // column is no more checked
	               if (CN.C[c].range==1) {
	                 if (rand01() <= CN.C[c].probab_oc)
	                   check(CN.C[c].from,(FX==large)? fx:fxi) && check(CN.C[c].to,(TX==large)? tx:txi) && (edges.push_back(edge(CN.C[c].from,(FX==large)? fx:fxi,fy,CN.C[c].to,(TX==large)? tx:txi,ty,&CN.C[c],delayTime,currEdgeWeightFactor)),ec++);
	               } else {  //longrange
	                 if (rand01() <= CN.C[c].probab)
	                   check(CN.C[c].from,(FX==large)? fx:fxi) && check(CN.C[c].to,(TX==large)? tx:txi) && (edges.push_back(edge(CN.C[c].from,(FX==large)? fx:fxi,fy,CN.C[c].to,(TX==large)? tx:txi,ty,&CN.C[c],delayTime,currEdgeWeightFactor)),ec++);
	               } //longrange
	             }//radius
		}
             	               	             	  	     
            }//tx
          }//fx

      //dbg
      (*ConnSummaryFile)<<c+1<<" of "<<CN.n<<" , edges: "<<ec<<" \ttype: (";
      CT.dump(CN.C[c].from);(*ConnSummaryFile)<<") -> (";CT.dump(CN.C[c].to);(*ConnSummaryFile)<<") ";
      CN.dump(c);
    }//c
} //generate_MultiLayer_connections


//for given pair fx,tx of homotopic neurons create neighbourhood around tx and connect it to fx (works for both hemispheres)
int Edges::connect_local_neighbours(int c /*rule*/,int fx, int tx,bool neighbourhood_hemisphere){
    int ec=0;//edges counter
    for (int nbtx=0; map_LH_RH[nbtx]!=-1; nbtx++){  //should be the same as largeLH
      // printf("fx: %d, tx: %d, nbtx: %d, right_hem:%d ", fx, tx,nbtx, neighbourhood_hemisphere);
      double d=dist3D(tx, nbtx, neighbourhood_hemisphere);
      if ( d>=CN.C[c].radius_min && d<=CN.C[c].radius_max ){
        assert(CN.C[c].range==1);//no long range done here
        if (rand01() <= CN.C[c].probab_oc) {
	  int from = neighbourhood_hemisphere ? CN.C[c].from : CN.C[c].to;
	  int to   = neighbourhood_hemisphere ? CN.C[c].to   : CN.C[c].from;
    double dist = dist3D(from, to, neighbourhood_hemisphere);
          double delayTime = getScaledDelay(d);
          check(from,fx) && check(to,nbtx) && (edges.push_back(edge(from,fx,0,to,nbtx,0,&CN.C[c], delayTime)),ec++); //ty==fy==0, see assert in caller
        }
      }//radius
    }//for nbtx
    return ec;
}//connect_local_neighbours

bool Edges::check(int type, int neuronx){
  return CT.cell[type].x > neuronx;
}





void Edges::generate_3D_intrahemispheric_connections(){ //generates connection using external distance matrix
    cerr<<"Generating intrahemispheric connections from pre-set 3D MRI data \n";

    int smallLH = CT.cell[CT.get_type("TC")].x; int largeLH = CT.cell[CT.get_type("CX")].x;
    int smallRH = CT.cell[CT.get_type("rTC")].x;int largeRH = CT.cell[CT.get_type("rCX")].x;

    cerr<<"Full set(CX): "<<largeLH<<" (L), "<<largeRH<<" (R), Subnet(TC,RE,IN...): "<<smallLH<<" (L), "<<smallRH<<" (R)\n";

    for (int c=0; c<CN.n; c++){ //layer connections
      if (CT.right_hemisphere(CN.C[c].from)==CT.right_hemisphere(CN.C[c].to)) continue;	//only connections between different hemispheres
      int FX = CT.cell[CN.C[c].from].x; int FY = CT.cell[CN.C[c].from].y; //input layer dimension
      int TX = CT.cell[CN.C[c].to].x; int TY = CT.cell[CN.C[c].to].y; //output layer dimension
      assert(FY==1); assert(TY==1); //MRI data are stored as 1D indexes
      assert(FX==largeLH); assert(TX==largeRH); // at this moment we allow only CX-rCX rules

      int ec=0;	//edges counter for statistics

      //HACK* before we get proper fine mesh file from Eran, at this moment we switch from fine mesh to coarse mesh //TODO remove hack
      FX=smallLH; TX=smallRH;

      for (int fx=0; fx<FX; fx++) {
        int tx=map_LH_RH[fx]; 		//fx-tx pair of homotopic neurons
	
	//HACK* use coarse mesh and remap it to fine mesh //TODO remove hack
        tx=map_LH_RH_small[fx];
        int fx2=Subnet[0][fx]; int tx2=Subnet[1][tx];

        // printf("1 c: %d,fx2: %d,tx2: %d tx: %d \n", c,fx2,tx2, tx);
        ec+=connect_local_neighbours(c,fx2,tx2,1/*right*/);
        // printf("0 c: %d,fx2: %d,tx2: %d\n", c,fx2,tx2);
        ec+=connect_local_neighbours(c,tx2,fx2,0/*left*/);
      }//fx

      //dbg
      (*ConnSummaryFile)<<c+1<<" of "<<CN.n<<" , edges: "<<ec<<" \ttype: (";
      CT.dump(CN.C[c].from);(*ConnSummaryFile)<<") -> (";CT.dump(CN.C[c].to);(*ConnSummaryFile)<<") ";
      CN.dump(c);
    }//c
}//generate_3D_intrahemi_connections

int main(int argc, char *argv[]){
  cerr << "start program" << endl;
  // int * temp = new int[20500];
  // int * tempp = new int[20500];
  // these need to be declared as local valieables so that the compiler doesnt know the sizes
  
  int mriPointsFull = 20500;
  int thalPopulations = 2*642; 
  int numParam = 2;
  int numLayers = 6;

  // These variables take too much memory so we need to initialize them at runtime manually
  // for some reason compiler didnt like initialization in single call
  // Dist_Prob = new float[MRI_POINTS_FULL][MRI_POINTS_FULL][NUM_PARAM];
     


  cerr << 1 << endl;
  Dist_Prob = new float**[mriPointsFull];
  for(int i = 0; i < mriPointsFull; i++){
    Dist_Prob[i] = new float*[mriPointsFull];  
    for(int j = 0; j < mriPointsFull; j++){
      Dist_Prob[i][j] = new float[numParam];  
      for(int k = 0; k < numParam; k++){
        Dist_Prob[i][j][k] = -1;
      }
    }
  }
  cerr << 2 << endl;

  // weight_factor = new int[MRI_POINTS_FULL][MRI_POINTS_FULL][NUM_LAYERS][NUM_LAYERS];
  weight_factor = new float***[mriPointsFull];
  for(int i = 0; i < mriPointsFull; i++){
    weight_factor[i] = new float**[mriPointsFull];  
    for(int j = 0; j < mriPointsFull; j++){
      weight_factor[i][j] = new float*[numLayers];  
      for(int k = 0; k < numLayers; k++){
        weight_factor[i][j][k] = new float[numLayers];  
        for(int l = 0; l < numLayers; l++){
          weight_factor[i][j][k][l] = -1;
        }
      }
    }
  }
//  weight_factor = new int***[mriPointsFull];
//  for(int i = 0; i < mriPointsFull; i++){
//    weight_factor[i] = new int**[mriPointsFull];  
//    for(int j = 0; j < mriPointsFull; j++){
//      weight_factor[i][j] = new int*[numLayers];  
//      for(int k = 0; k < numLayers; k++){
//        weight_factor[i][j][k] = new int[numLayers];  
//        for(int l = 0; l < numLayers; l++){
//          weight_factor[i][j][k][l] = -1;
//        }
//      }
//    }
//  }
  cerr << 3 << endl;
  weight_factor_INPY = new float*[mriPointsFull];
  for(int i = 0; i < mriPointsFull; i++){
    weight_factor_INPY[i] = new float[mriPointsFull];  
    for(int j = 0; j < mriPointsFull; j++){
      weight_factor_INPY[i][j] = -1;
    }
  }
  cerr << 4 << endl;

   Dist_intraThalamus = new float*[thalPopulations];
   for(int i =0; i< thalPopulations; i++){
      Dist_intraThalamus[i] = new float[thalPopulations];
      for(int j =0; j< thalPopulations; j++){
          Dist_intraThalamus[i][j] = -1;
      }
   }

  cerr << 5 << endl;


//     for some reason the code doesnt like when we use the #define variables 
//     printf("1");
//  Dist_Prob = new float**[MRI_POINTS_FULL];
//  for(int i = 0; i < MRI_POINTS_FULL; i++){
//    Dist_Prob[i] = new float*[MRI_POINTS_FULL];  
//    for(int j = 0; j < MRI_POINTS_FULL; j++){
//      Dist_Prob[i][j] = new float[NUM_PARAM];  
//      for(int k = 0; k < NUM_PARAM; k++){
//        Dist_Prob[i][j][k] = 0;
//      }
//    }
//  }
//     printf("2");
//  // weight_factor = new int[MRI_POINTS_FULL][MRI_POINTS_FULL][NUM_LAYERS][NUM_LAYERS];
//  weight_factor = new int***[MRI_POINTS_FULL];
//  for(int i = 0; i < MRI_POINTS_FULL; i++){
//    weight_factor[i] = new int**[MRI_POINTS_FULL];  
//    for(int j = 0; j < MRI_POINTS_FULL; j++){
//      weight_factor[i][j] = new int*[NUM_LAYERS];  
//      for(int k = 0; k < NUM_LAYERS; k++){
//        weight_factor[i][j][k] = new int[NUM_LAYERS];  
//        for(int l = 0; l < NUM_LAYERS; l++){
//          weight_factor[i][j][k][l] = 0;
//        }
//      }
//    }
//  }
//     printf("3");

  //connection summary file
  ofstream connSummaryFile;
  ConnSummaryFile = &connSummaryFile;
  connSummaryFile.open ("connSummaryFileexample.txt");
  //any commandline parameter taken as network connections file
  if (argc==2){
    parse_config("mriNetwork.cfg");
    CE.load_MRI_network(argv[1]); 		//connection file we get verbatim from ucsd group
  }
  else if (argc==3) {				//partly processed data from ucsd folks - only subnet and distances given
    parse_config("mriNetwork.cfg");
    load_3D_subnet(argv[2],0);
    load_3D_distances(argv[1],0);
    CE.generate_3D_connections(0);
   } 
   else if (argc==7) {				//partly processed data from ucsd folks - both hemispheres and their connections given
    useDelays = true;
    parse_config("mriNetwork.cfg");
    assert(both_hemispheres_detected);
    load_3D_subnet(argv[2],0); //left
    load_3D_subnet(argv[4],1); //right
    load_3D_LH_RH_correspondence(argv[5],map_LH_RH);
    load_3D_LH_RH_correspondence(argv[6],map_LH_RH_small);
    load_3D_distances(argv[1],0);
    load_3D_distances(argv[3],1);
    CT.MRI_reduce();  //reduction of number neurons in case ratio used in config file
    //we do this independently for distinct hemipheres, because of memory 
    //(some crazy radius scenarios get up to 200 GB of RAM)
    CE.generate_3D_connections(0);
  //  CE.apply_synaptic_types();  
  //  CE.write_connections(0);
  //  CE.edges.clear();
    
    CE.generate_3D_connections(1);
  //  CE.apply_synaptic_types();
  //  CE.write_connections(1);
  //  CE.edges.clear();

    //intrahemispheric connections
    CE.generate_3D_intrahemispheric_connections();
    CE.apply_synaptic_types();
    CE.write_connections(0);
    CE.edges.clear();
    return 0;
   } 
  //else if(argc==8){
  else if(argc>=8){
  cerr << "Entered 8 parameters function call!! ;)" << endl;

    useDelays = true;
    parse_config("mriNetworkNew.cfg");
    assert(both_hemispheres_detected);
    load_3D_subnet(argv[2],0); //left
    load_3D_subnet(argv[4],1); //right
    load_3D_LH_RH_correspondence(argv[5],map_LH_RH);
    load_3D_LH_RH_correspondence(argv[6],map_LH_RH_small);

    mapLayerNameId_IN= fillDict_IN(mapLayerNameId_IN);
    mapLayerNameId_PY= fillDict_PY(mapLayerNameId_PY);

   load_3D_dist_prob(argv[1]);
    load_weightFactors_INPY(argv[7]);
//    load_weightFactors_binary(argv[3]);
    load_weight_factors(argv[3]);

    load_thalCort_dist(argv[8]);
    load_intraThalamus_dist(argv[9]);

    CT.MRI_reduce();  //reduction of number neurons in case ratio used in config file
    //we do this independently for distinct hemipheres, because of memory
    //(some crazy radius scenarios get up to 200 GB of RAM)
    CE.generate_MultiLayer_connections(0);
  //  CE.apply_synaptic_types();
  //  CE.write_connections(0);
  //  CE.edges.clear();

    CE.generate_MultiLayer_connections(1);
  //  CE.apply_synaptic_types();
  //  CE.write_connections(1);
  //  CE.edges.clear();

    //intrahemispheric connections
//    CE.generate_3D_intrahemispheric_connections();
    CE.apply_synaptic_types();
    CE.write_connections(0);
    CE.edges.clear();
    return 0;

  }else{
     parse_config("mriNetwork.cfg");
     CE.generate_connections(false); 		//our own random connectivity
   }

  CE.apply_synaptic_types();	//needed in case we need sample parameters for strength/mini_X parameters
  //CE.dump_from_neuron_edges(4,20,20);
  //CE.dump_from_neuron_edges(0,0,0,2);
  CE.write_connections();
  //CT.dump();
  //CN.dump();
  return 0;
}
