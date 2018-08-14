/*
 *  Dr. Rohitash Chandra, 2009  (c.rohitash@gmail.com)

   aicrg.softwarefoundationfiji.org


This program build on the G3-PCX C code from here: http://www.iitk.ac.in/kangal/codes.shtml

The paper on the G3-PCX: " Kalyanmoy Deb, Ashish Anand, and Dhiraj Joshi, A Computationally Efficient Evolutionary Algorithm for Real-Parameter Optimization, Evolutionary Computation 2002 10:4, 371-395"

 In this code, CC is used for solving general function optimisation problems (Ellip, Rosenbrock, Schwefel's).


 */
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <ctime>


#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<ctype.h>

time_t TicTime;
time_t TocTime;
using namespace std;

typedef vector<double> Layer;
typedef vector<double> Nodes;
typedef vector<double> Frame;
typedef vector<int> Sizes;
typedef vector<vector<double> > Weight;
typedef vector<vector<double> > Data;




const int maxgen = 200  ; //max number of function eval. (termination criteria)




 #define rosen          // choose the function:


#define EPSILON 1e-40

#define MINIMIZE 1      //set 1 to minimize and -1 to maximize
#define LIMIT 1e-20     //accuracy of best solution fitness desired
#define KIDS 2          //pool size of kids to be formed (use 2,3 or 4)
#define M 1             //M+2 is the number of parents participating in xover (use 1)
#define family 2        //number of parents to be replaced by good individuals(use 1 or 2)
#define sigma_zeta 0.01
#define sigma_eta 0.01   //variances used in PCX (best if fixed at these values)

#define  PopSize 20

#define NPSize KIDS + 2   //new pop size

#define RandParent M+2     //number of parents participating in PCX


double d_not[PopSize];

double seed,basic_seed;






//*************************************************


//-------------------------------------------------------

class Individual{
    friend class GeneticAlgorithm;
    protected:
        Layer Chrome;
        double Fitness;
        Layer BitChrome;
    public:
        Individual(){
        }
        void print();
};

class GeneticAlgorithm{
    protected:
        Individual Population[PopSize];
        int TempIndex[PopSize];
        Individual  NewPop[NPSize];
        int mom[PopSize];
        int list[NPSize];
        int MaxGen;
        int NumVariable;
        double BestFit;
        int BestIndex;
        int NumEval;
        int kids;
        int TotalEval;
        double Error;
        double Cycles;
        bool Sucess;

    public:
        GeneticAlgorithm(int stringSize){
    	   NumVariable = stringSize;
    	   NumEval=0;
        }

        int GetEval(){
            return TotalEval;
        }

        double GetCycle(){
            return Cycles;
        }

        double GetError(){
            return Error;
        }

        bool GetSucess(){
            return Sucess;
        }
        double Fitness(){
            return BestFit;
        }
        double RandomWeights();
        double RandomAddition();
        void PrintPopulation();
        int GenerateNewPCX(int pass);
        double Objective(Layer x);
        void InitilisePopulation();
        void Evaluate();
        double modu(double index[]);
        double innerprod(double Ind1[],double Ind2[]);
        double RandomParents();
        double MainAlgorithm(double RUN, ofstream &out1, ofstream &out2, ofstream &out3);
        double rand_normal(double mean, double stddev) ;
        void my_family();   //here a random family (1 or 2) of parents is created who would be replaced by good individuals
        void find_parents() ;
        void rep_parents() ;  //here the best (1 or 2) individuals replace the family of parents
        void sort();
};

   //-------------------------------


double GeneticAlgorithm::RandomWeights()
{
    int chance;
    double randomWeight;
    double NegativeWeight;
    // srand ( time(NULL) );
    chance =rand()%2;
    if(chance ==0){
        randomWeight =rand()%  100000;
        return randomWeight*0.01;
    }
    if(chance ==1){
        NegativeWeight =rand()% 100000;
        return NegativeWeight*0.01;
    }
}


double GeneticAlgorithm::rand_normal(double mean, double stddev) { //Box Muller Random Numbers
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached) {
        // choose a point x,y in the unit circle uniformly at random
        double x, y, r;
        do {
            //  scale two random integers to doubles between -1 and 1
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;
            r = x*x + y*y;
        } while (r == 0.0 || r > 1.0);
        // Apply Box-Muller transform on x, y
        double d = sqrt(-2.0*log(r)/r);
        double n1 = x*d;
        n2 = y*d;
        // scale and translate to get desired mean and standard deviation
        double result = n1*stddev + mean;
        n2_cached = 1;
        return result;
    }
    else {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

double GeneticAlgorithm::RandomAddition()
{
      int chance;
      double randomWeight;
      double NegativeWeight;
      chance =rand()%2;

      if(chance ==0){
          randomWeight =rand()% 100;
          return randomWeight*0.009;
      }

      if(chance ==1){
          NegativeWeight =rand()% 100;
          return NegativeWeight*-0.009;
     }
}



void GeneticAlgorithm::InitilisePopulation(){
    double x, y;
    for(int row = 0; row < PopSize  ; row++) {
        for(int col = 0; col < NumVariable ; col++){
            //   x=randomperc();    // x is a uniform random number in (0,1)
            //	y=(-10.0)+(5.0*x); // the formula used is y=a+(b-a)*x if y should be a random number in (a,b)
            Population[row].Chrome.push_back(RandomWeights());
        }
    }
    for(int row = 0; row < NPSize  ; row++) {
        for(int col = 0; col < NumVariable ; col++)
            NewPop[row].Chrome.push_back(0);
        }
    }

void GeneticAlgorithm::Evaluate(){
    // solutions are evaluated and best id is computed
    Population[0].Fitness= Objective( Population[0].Chrome);
    BestFit = Population[0].Fitness;
    BestIndex = 0;
    for(int row = 0; row < PopSize  ; row++){
        Population[row].Fitness= Objective( Population[row].Chrome);
        if ((MINIMIZE * BestFit) > (MINIMIZE * Population[row].Fitness)){
            BestFit = Population[row].Fitness;
            BestIndex = row;
        }
        NumEval++;
    }
}


void GeneticAlgorithm::PrintPopulation(){
    for(int row = 0; row < PopSize  ; row++) {
        for(int col = 0; col < NumVariable ; col++){
            cout<< Population[row].Chrome[col]<<" ";
        }
        cout<<endl;
    }

    for(int row = 0; row < PopSize  ; row++)
        cout<< Population[row].Fitness<<endl;
    cout<<" ---"<<endl;

    cout<<BestFit<<"  "<<BestIndex<<endl;
    for(int row = 0; row < NPSize  ; row++){
        for(int col = 0; col < NumVariable ; col++)
            cout<< NewPop[row].Chrome[col]<<" ";
        cout<<endl;
    }
}


double GeneticAlgorithm:: Objective(Layer x){
    int i,j,k;
    double fit, sumSCH;
    fit=0.0;

    #ifdef ellip
        // Ellipsoidal function
        for(j=0;j<NumVariable;j++)
            fit+=((j+1)*(x[j]*x[j]));
    #endif

    #ifdef schwefel
        // Schwefel's function
        for(j=0; j<NumVariable; j++){
            for(i=0,sumSCH=0.0; i<j; i++)
                sumSCH += x[i];
            fit += sumSCH * sumSCH;
        }
    #endif

    #ifdef rosen
        //   Rosenbrock's function
        for(j=0; j<NumVariable-1; j++)
            fit += 100.0*(x[j]*x[j] - x[j+1])*(x[j]*x[j] - x[j+1]) + (x[j]-1.0)*(x[j]-1.0);
    #endif
  return(fit);
}
//------------------------------------------------------------------------
//here a random family (1 or 2) of parents is created who would be replaced by good individuals
void GeneticAlgorithm:: my_family(){
    int i,j,index;
    int swp;
    double u;

    for(i=0;i<PopSize;i++)
    mom[i]=i;

    for(i=0;i<family;i++){
        //   u=randomperc();
        //   0index=(u*(PopSize-i))+i;
        index = (rand()%PopSize) +i;

        cout<<"is index         --------------------------------   +++++++++++++++++++++++++++++++++++++++++++++++++++++++ "<<index<<endl;
        if(index>(PopSize-1))
            index=PopSize-1;
        swp=mom[index];
        mom[index]=mom[i];
        mom[i]=swp;
    }
}

//here the parents to be replaced are added to the temporary sub-population to assess their goodness against the new solutions formed which will be the basis of whether they should be kept or not

void GeneticAlgorithm::find_parents(){
    int i,j,k;
    double u,v;

    my_family();
    //cout<<kids<<endl;
    for(j=0;j<family;j++){
        for(i=0;i<NumVariable;i++)
            NewPop[kids+j].Chrome[i] = Population[mom[j]].Chrome[i];
        NewPop[kids+j].Fitness = Objective(NewPop[kids+j].Chrome);
        NumEval++;
    }
}


//here the best (1 or 2) individuals replace the family of parents
void GeneticAlgorithm::rep_parents() {
    int i,j;
    for(j=0;j<family;j++){
        for(i=0;i<NumVariable;i++)
            Population[mom[j]].Chrome[i]=NewPop[list[j]].Chrome[i];
        Population[mom[j]].Fitness = Objective(Population[mom[j]].Chrome);
        cout<< mom[j]<<" "<< Population[mom[j]].Fitness<<  "                                     @@@@" <<endl;
        NumEval++;
    }
}


void GeneticAlgorithm::sort(){
    int i,j, temp;
    double dbest;
    for (i=0;i<(kids+family);i++)
        list[i] = i;
    if(MINIMIZE){
        for (i=0; i<(kids+family-1); i++){
            dbest = NewPop[list[i]].Fitness;
            for (j=i+1; j<(kids+family); j++){
                if(NewPop[list[j]].Fitness < dbest){
                    dbest = NewPop[list[j]].Fitness;
                    temp = list[j];
                    list[j] = list[i];
                    list[i] = temp;
                }
            }
        }
    }
    else{
        for (i=0; i<(kids+family-1); i++){
            dbest = NewPop[list[i]].Fitness;
            for (j=i+1; j<(kids+family); j++){
                if(NewPop[list[j]].Fitness > dbest){
                    dbest = NewPop[list[j]].Fitness;
                    temp = list[j];
                    list[j] = list[i];
                    list[i] = temp;
                }
            }
        }
    }
}


//---------------------------------------------------------------------
double GeneticAlgorithm::modu(double index[]){
    int i;
    double sum,modul;
    sum=0.0;

    for(i=0;i<NumVariable ;i++)
        sum+=(index[i]*index[i]);

    modul=sqrt(sum);
    return modul;
}

// calculates the inner product of two vectors
double GeneticAlgorithm::innerprod(double Ind1[],double Ind2[]){
    int i;
    double sum = 0.0;
    for(i=0;i<NumVariable ;i++)
        sum+=(Ind1[i]*Ind2[i]);
    return sum;
}

int GeneticAlgorithm::GenerateNewPCX(int pass){
    int i,j,num,k;
    double Centroid[NumVariable];
    double tempvar,tempsum,D_not,dist;
    double tempar1[NumVariable];
    double tempar2[NumVariable];
    double D[RandParent];
    double d[NumVariable];
    double diff[RandParent][NumVariable];
    double temp1,temp2,temp3;
    int temp;

    for(i=0;i<NumVariable;i++)
        Centroid[i]=0.0;

    // centroid is calculated here
    for(i=0;i<NumVariable;i++){
        for(j=0;j<RandParent;j++){
            cout<< TempIndex[j]<<" ";
            Centroid[i]+=Population[TempIndex[j]].Chrome[i];
        }
        Centroid[i]  /= RandParent;
        //for(i=0;i<NumVariable;i++)
        //    cout<<Centroid[i] << " ";
        // cout<< j << "     --   centroid "<<endl;
    }
    // cout<< "  ****                                                     centroid done "<<endl;
    // calculate the distace (d) from centroid to the index parent arr1[0]
    // also distance (diff) between index and other parents are computed
    for(j=1;j<RandParent;j++){
        for(i=0;i<NumVariable;i++){
            if(j == 1)
                d[i]=Centroid[i]-Population[TempIndex[0]].Chrome[i];
            diff[j][i]=Population[TempIndex[j]].Chrome[i]-Population[TempIndex[0]].Chrome[i];
        }
        /*if (modu(diff[j]) < EPSILON){
            cout<< "RUN Points are very close to each other. Quitting this run   " <<endl;
            return (0);
        }*/
    }
    for(i=0;i<NumVariable;i++)
        cout<<Centroid[i] << " ";
    cout<<   "     --   centroid "<<endl;

    for(i=0;i<NumVariable;i++)
        cout<<d[i] << " ";
    cout<<"                              -----------       d "<<endl;
    for(j=1;j<RandParent;j++){
        for(i=0;i<NumVariable;i++)
        cout<<diff[j][i]<< " ";
        cout<<endl;
    }
    cout<<"                              -----------       diff ---------------- "<<endl;
    dist=modu(d); // modu calculates the magnitude of the vector
    if (dist < EPSILON){
      cout<< "RUN Points are very close to each other. Quitting this run    " <<endl;
      return (0);
    }
    // orthogonal directions are computed (see the paper)
    for(i=1;i<RandParent;i++){
        temp1=innerprod(diff[i],d);
        temp2=temp1/(modu(diff[i])*dist);
        temp3=1.0-pow(temp2,2.0);
        D[i]=modu(diff[i])*sqrt(temp3);
    }

    //for(i=1;i<RandParent;i++)
    // cout<<D[i]<< " ";
    // cout<<"                              -----------       D  ----------------  ++++++++++++++ "<<endl;


    D_not=0;
    for(i=1;i<RandParent;i++)
        D_not+=D[i];
    D_not/=(RandParent-1); //this is the average of the perpendicular distances from all other parents (minus the index parent) to the index vector
    // Next few steps compute the child, by starting with a random vector
    //  double rndx[] = {-0.03183216, -0.09051064, -0.11249914,  0.03915154, -0.10039502};
    for(j=0;j<NumVariable;j++){
        tempar1[j]=rand_normal(0, D_not*sigma_eta);
        tempar2[j]=tempar1[j];
    }
    for(j=0;j<NumVariable;j++){
        tempar2[j] = tempar1[j]-((innerprod(tempar1,d)*d[j])/  pow(dist,2.0));
    }
    for(j=0;j<NumVariable;j++)
        tempar1[j]=tempar2[j];


    //for(j=0;j<NumVariable;j++)
        // cout<<tempar1[j]<<" ";

    //  cout<<"                              - ---------------- tempar1 ++++++++++++++ "<<endl;
    for(k=0;k<NumVariable;k++){
        NewPop[pass].Chrome[k]=Population[TempIndex[0]].Chrome[k]+tempar1[k];
        cout<<  NewPop[pass].Chrome[k]<<" ";
    }
    cout<<pass<<"  -  ------------------------->>>> --------------------pcx"<<endl;
    tempvar=rand_normal(0, sigma_zeta);
    for(k=0;k<NumVariable;k++){
        NewPop[pass].Chrome[k] += (tempvar*d[k]);
        cout<<  NewPop[pass].Chrome[k]<<" + ";
        //NewPop[pass].Chrome[k] += (tempvar*d[k]);
    }
    cout<<pass<< "  -  --------------------------------------------------------  "<<endl;
    NewPop[pass].Fitness = Objective( NewPop[pass].Chrome);
    NumEval++;
    return (1);
}




//------------------------------------------------------------------------
double GeneticAlgorithm::  RandomParents(){
    int i,j,index;
    int swp;
    double u;
    int delta;
    for(i=0;i<PopSize;i++)
        TempIndex[i]=i;
    swp=TempIndex[0];
    TempIndex[0]=TempIndex[BestIndex];  // best is always included as a parent and is the index parent
    // this can be changed for solving a generic problem
    TempIndex[BestIndex]=swp;
    // shuffle the other parents
    for(i=1;i<RandParent;i++){
        //u=randomperc();
        index=(rand()%PopSize)+i;
        cout<<index<< " is index random parent"<<endl;
        if(index>(PopSize-1))
        index=PopSize-1;
        swp=TempIndex[index];
        TempIndex[index]=TempIndex[i];
        TempIndex[i]=swp;
    }
}


void TIC( void ){
    TicTime = time(NULL);
}

void TOC( void ){
    TocTime = time(NULL);
}

double StopwatchTimeInSeconds(){
    return difftime(TocTime,TicTime);
}

double GeneticAlgorithm:: MainAlgorithm(double RUN, ofstream &out1, ofstream &out2, ofstream &out3 ){
    double x = 0;
    ifstream inFile;
    inFile.open("pop.txt");
    if (!inFile){
        cout << "Unable to open file";
        exit(1); // terminate with error
    }

    clock_t start = clock();
    double tempfit =0;
    int count =0;
    int tag;
    kids = KIDS;

    InitilisePopulation();

    for(int j=0;j<PopSize;j++){
        for(int i=0;i<NumVariable;i++)
            inFile >> Population[j].Chrome[i];
    }


    for(int j=0;j<PopSize;j++){
        for(int i=0;i<NumVariable;i++)
            cout<<Population[j].Chrome[i]<<" * ";
        cout<<endl;
    }
    Evaluate();
    tempfit=Population[BestIndex].Fitness;
    while( NumEval < maxgen){
        RandomParents();           //random array of parents to do PCX is formed
        //int ind[] = {9, 17, 6, 3, 4, 5, 2, 7, 8, 0 , 10 , 11 , 12, 13, 14 , 15 , 16, 1 , 18,  19};
        for(int j=0;j<PopSize;j++){
            // TempIndex[j] = ind[j];
            cout<<TempIndex[j]<< " ";
            out3<<TempIndex[j]<< " ";
        }
        cout<<endl;
        out3<<endl;
        for(int i=0;i<kids;i++){
            cout<<i<< "                               ---------------- *****  ----------------  ----------------------------------- # # #"<<endl;
        	tag = GenerateNewPCX(i);
       	    if (tag == 0) break;
        }
        //  if (tag == 0) break;
        find_parents();  // form a pool from which a solution is to be
        //   replaced by the created child
        for(int j=0;j<PopSize;j++){
            // TempIndex[j] = ind[j];
            cout<<mom[j]<< " ~ ";
            out2<<mom[j]<< " ";
        }
        cout<<endl;
        out2<<endl;
        for(int j=0;j<NPSize;j++){
            for(int i=0;i<NumVariable;i++)
                cout<<NewPop[j].Chrome[i] << " ";
            cout<<j<< " find parents "<< NewPop[j].Fitness<<endl;
        }
        sort();          // sort the kids+parents by fitness
        for(int j=0;j<NPSize;j++){
            for(int i=0;i<NumVariable;i++)
                cout<<NewPop[j].Chrome[i] << " ";
                cout<<j<< " sort parents "<< NewPop[j].Fitness<<endl;
        }
        rep_parents();   // a chosen parent is replaced by the child
       	//finding the best in the population
        BestIndex=0;
        tempfit=Population[0].Fitness;
        for(int j=0;j<NPSize;j++){
            for(int i=0;i<NumVariable;i++)
                cout<<NewPop[j].Chrome[i] << " ";
            cout<<j<< "  rep_parents "<< NewPop[j].Fitness<<endl;
        }
        //cout<<tempfit<<endl;
        for(int i=1;i<PopSize;i++){
            if((MINIMIZE * Population[i].Fitness) < (MINIMIZE * tempfit)){
                tempfit=Population[i].Fitness;
                BestIndex=i;
            }
        }
        for(int j=0;j<PopSize;j++){
            for(int i=0;i<NumVariable;i++)
                cout<<Population[j].Chrome[i]<<" * ";
            cout<<endl;
        }
        for(int j=0;j<PopSize;j++){
            cout<<Population[j].Fitness << "       ";
            // out2<<Population[j].Fitness<< " ";
        }
        cout<<endl;
        // out2<<endl;
        cout<< NumEval <<"   "  << Population[BestIndex].Fitness<< "  "<< BestIndex<<"   --------------- *** --------------***---------------------*** ---------- "<<endl;
	}
    out1 <<"   ---------------  "<<endl;
    clock_t finish = clock();
    Cycles = ((double)(finish - start))/CLOCKS_PER_SEC;
    cout<<Cycles<<" ----"<<endl;
    Error = tempfit;
    TotalEval= NumEval;

    for(int i=0;i<NumVariable;i++)
        cout<<Population[BestIndex].Chrome[i]<<"  ***    ";
    cout<<"---->  "<< tempfit<<"     "<<count*kids<<"      "<<Cycles<<endl;
    //out2<<"Fitness of this best solution:"<<tempfit<<endl;

    cout<<"Best solution obtained after X function evaluations:"<<count*kids<<" "<<NumEval<<endl;

    for(int i=0;i<NumVariable;i++)
        cout<<Population[BestIndex].Chrome[i]<<" ";
    cout<<endl;

    cout<<"Fitness of this best solution:"<<tempfit<<" "<<StopwatchTimeInSeconds()<<endl;
    inFile.close();
}


int main(void)
{

    int VSize =2; //number of variables (dimension) for problem

    ofstream out1;
    out1.open("out1.txt");
    ofstream out2;
         out2.open("out2.txt");
     	ofstream out3;
     	     out3.open("out3.txt");

     	    Sizes EvalAverage;
     	       	 Layer ErrorAverage;
     	       	 Layer CycleAverage;

     	       	 int MeanEval=0;
     	       	 double MeanError=0;
     	       double MeanCycle=0;

     	       	 int EvalSum=0;
     	       	 double ErrorSum=0;
     	       	double CycleSum=0;

     	       	int maxrun=1;
    for(int RUN=1;RUN<=maxrun;RUN++)
         {
        GeneticAlgorithm GenAlg(VSize);
        GenAlg.MainAlgorithm(RUN,out1,out2, out3);


         EvalAverage.push_back(GenAlg.GetEval());
         MeanEval+=GenAlg.GetEval();

         ErrorAverage.push_back(GenAlg.GetError());
        MeanError+= GenAlg.GetError();
            cout<<GenAlg.GetCycle()<<" ----------------"<<endl;
        CycleAverage.push_back(GenAlg.GetCycle());
             MeanCycle+= GenAlg.GetCycle();


         }


     MeanEval=MeanEval/EvalAverage.size();
     MeanError=MeanError/ErrorAverage.size();
     MeanCycle=MeanCycle/CycleAverage.size();


    //out3<< MeanEval<<"   "<<MeanError<<"    "<<MeanCycle<<  endl;
    EvalAverage.empty();
    ErrorAverage.empty();
    CycleAverage.empty();

    out1.close();
    out2.close();
    out3.close();


 return 0;

};
