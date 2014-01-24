#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/*Define structure to get next event*/
struct nev{
	int i,j,k;
	double rtot;
};

/*Output of each simulation*/
struct outsim{ //PP for each sim we need to record the time to failure, time of exit of sanctuary (?), fcc, fcompf, ssdc, and viral load. 
	double time_failure,time_exit;
	int fcc,fcompf,ssdc;  
	float vltot; //average viral load per day
};

/*Prototype functions*/
void migration_rates();
void r0_calculation();
void mutation_rates();
struct outsim one_run();
struct nev next_event(); //PP to determine what is the next event 
int choose(double rmax,float rates[4][4][3],int jik[]); // to "choose" event when rates are already calculated, using random uniform number.  

/*Define global variables*/
FILE *tfailure, *texit, *firstfilled, *pathsm, *pathsdc, *comp2, *comp3,*mmts, *me, *t3, *mrp, *tt, *path, *pathm, *vload, *fa, *parameters;
float m, mig[4][5]={0}; //PP why 4x5 dimension? We only use the 4x4 part of the matrix, but may use the rest for other things. 
int tcf; // threshold value for compartment infected
float R, uf1, ub1, uf2, ub2, c1, c2, s, dy, dx; // PP R is R0 for WT in sanctuary, uf are forward mut rates, ub backward mut rates. c is strength of the drug.
float R0[4][4]; // R0 of strain i in compartment j
float comp[5][4]={0};  // number of strains in each compartment. i:0-3 strains, j: compartment, i=4 total number of uninfected cells when there is no infection, so comp[4][0] = num uninfected cells in sanctuary. 
double mutrate[4][4]={0}; // probability of mutation from strain of type i to strain of type j
double t; //PP I guess this is time. 
double mmt; // counts how many migration events of a mutant from the sanctuary to the SDC have occurred 
unsigned long int seed; //the seed
int outputmutantinfo=0; //if 1 output info on mutants for debugging
int rununtilfailureOrngens=0;//if 0 run until failure if >0 run this num of generations

/*Random number generation*/
gsl_rng * r; /*global random number generator*/

int main () {	
	/*Random number generation*/
	const gsl_rng_type * T; // define pointer for the type of random number generaton
	T = gsl_rng_default; // use the default random number generator algorithm (Mersenne twister)
	gsl_rng_env_setup(); 
	r = gsl_rng_alloc (T); // allocate memory for generator of type T

	/*Variable definition*/
	int rep,nr; // nr = Number of replicates, rep is current replicate	
	struct outsim results; //PP results will hold the results for one simulation run

	/*Read parameters*/
	parameters=fopen("parameters.txt", "rt"); // "rt" = Open a text file for reading. (The file must exist.)

	fscanf(parameters,"R:%f c1:%f c2:%f s:%f dy:%f dx:%f m:%f uf1:%f ub1:%f uf2:%f ub2:%f x0:%f x1:%f x2:%f x3:%f nr:%d tcf:%d seed:%lu ruf:%d omi:%d", &R, &c1, &c2, &s, &dy, &dx, &m, &uf1, &ub1, &uf2, &ub2, &comp[4][0], &comp[4][1], &comp[4][2], &comp[4][3], &nr, &tcf, &seed, &rununtilfailureOrngens, &outputmutantinfo);
	
	//	gsl_rng_set(r, time(NULL)); // Use clock to set the seed of the random number generator
	gsl_rng_set(r, seed); // Use clock to set the seed of the random number generator
		
	/*Calculate migration rates*/
	migration_rates();
	
	/*Calculate R0 values */
	r0_calculation();
	
	/*Calculate mutation rates */
	mutation_rates();
	
	/*Opening files to print results*/
	tfailure=fopen("tfailure.txt","a+"); // time to treatment failure // "a+" means append to the file 
	texit=fopen("texit.txt","a+"); // time to exit the sanctuary
	firstfilled=fopen("firstfilled.txt", "a+"); // first compartment that is colonized (values: 0,1,2,3) //PP why "a" and not "a+" ?
	pathsm=fopen("path_sucmut.txt", "a+"); // origin of first cell that colonizes the double-drug compartment (values: -1,0,1,2,10,11,12) 
	//(-1: no treatment failure | 0,1,2: cell present in the DDC that mutates from wt,m1 or m2 to m12 | 10,11,12: m12 cell that migrates from the sanctuary, SDC1 or SDC2 to the DDC) // PP this is called fcc in the function that runs a sim. 
	pathsdc=fopen("path_sdc.txt", "a+"); // state of the SDC when the DDC is colonized (values:-1,0,1,2,3) 
	//(-1: no TF | 0: only the sanctuary | 1: sanctuary + SDC1 | 2: sanctuary + SDC2 | 3: both SDC )
	vload=fopen("vload.txt", "a+"); // average vload per day
	mmts=fopen("mmt.txt", "a+"); // rate of migrantion of SD resistant mutants from sacntuary to SDC 
	
	
	/*Running all the replicates*/
	for(rep=0;rep<nr;rep++){ //PP nr is the number of replicates, rep is the number of the current simulation
		//printf("rep %d\n",rep);
		mmt=0.0;		
		/*Run one simulation*/
		results=one_run(); //PP results will hold the results for one simulation run
		
		/************** Print results of each simulation ****************/
		
		/*total time to treatment failure*/	
		fprintf(tfailure,"%f ", results.time_failure);
		
		/*time to exit the sanctuary*/	
		fprintf(texit,"%f ", results.time_exit);
		
		/*first compartment filled*/
		fprintf(firstfilled,"%d  ",results.fcompf);
		
		/*first cell that colonizes*/	
		fprintf(pathsm, "%d ", results.fcc);
		
		/*state of the SDC when the DDC is colonized*/	
		fprintf(pathsdc,"%d  ",results.ssdc);
		
		/*average viral load per day*/
		fprintf(vload, "%f ",results.vltot); 
		
		/*rate of migrantion of SD resistant mutants from sacntuary to SDC */
		fprintf(mmts, "%f ",mmt/results.time_failure);
	}
	
	/*Close files with results*/
	fclose(tfailure);
	fclose(texit);
	fclose(firstfilled);
	fclose(pathsm);
	fclose(pathsdc);
	fclose(vload);
	fclose(mmts);
	
	gsl_rng_free(r); // free memory associated to the random number generator
	return 0;
}

/****************************FUNCTION DEFINITIONS************************************/


/*FUNCTION TO CALCULATE MIGRATION RATES*/
void migration_rates(){
	float sumc; //PP HERE sumc is the sum of carrying capacities in all compartments except the one for which we are calculating the mig rate. 
	int i,j;
	
	/*Migration independent of carrying capacity*/
	/*
	 for (i=0; i<4; i++){ 
		for (j=0; j<4; j++)	
			if(i!=j) 
				mig[i][j]=m;
	 };
	 */
	
	/*Migration dependent of carrying capacity*/
	//PP: this fills a 4x4 matrix, but mig is a 4x5 matrix. why is that? 
	sumc=0;
	for (i=0; i<4; i++){ 
		sumc=0;
		for (j=0; j<4; j++)
			if(i!=j) //PP to make the sims fit with the analytical results, we can comment out this line and the one 3 lines further
				sumc+=comp[4][j];//comp[4][j] is the total carrying capacity of compartment j. 
		for (j=0; j<4; j++){
			if(i!=j) //PP to make the sims fit with the analytical results, we can comment out this line and the one 3 lines back	
				mig[i][j]=m*(comp[4][j]/sumc);
			else
				mig[i][j]=0;
		}
	}
	/*
	printf("migrat\n");
	for (i=0; i<4; i++) {
		for (j=0; j<4; j++) 
			printf("%f  ",mig[i][j]);
		printf("\n");
	}*/

}

/*FUNCTION TO CALCULATE R0 VALUES*/
void r0_calculation(){
	int i,j; 
	for (i=0; i<4; i++) {//i loops over the strains (i=0 WT, i=1 Res1, i=2 Res2, i=3 Res12)
		for (j=0; j<4; j++) {//j loops over the compartments (j=0 Sanctuary, j=1 Drug1, j=2 Drug2, j=3 Drug12)
			//PP set the R0 values as if there are no drugs
			if (i==0)
				R0[i][j]=R;
			else{
				if(i==1 || i==2)
					R0[i][j]=R*(1-s);
				else
					R0[i][j]=R*(1-s)*(1-s);
			}
			//PP add the effect of drugs in comp 1, 2, and 3 for strain 0, 1, and 2. 
			if((j==1 && (i==0 || i==2)) || (j==3 && (i==0 || i==2))){
				R0[i][j]*=(1-c1);
				if (j==3 && i==0) 
					R0[i][j]*=(1-c2);
			}
			else {
				if((j==2 && (i==0 || i==1)) || (j==3 &&  i==1)){
					R0[i][j]*=(1-c2);
				}
			}
		}
	};
	// !!! PP set R0 of strain 1 in comp 1 to 0 so that I can test mig rate
	if (rununtilfailureOrngens>0) R0[1][1]=0.0; //to test mig rate   
	if (rununtilfailureOrngens==0) R0[1][1]=1.081; //to test mig rate   
	
	/*
	printf("r0\n");
	for (i=0; i<4; i++) {
		for (j=0; j<4; j++) 
			printf("%f  ",R0[i][j]);
		printf("\n");
	} */
	
}


/*FUNCTION TO MAKE MATRIX OF MUTATION RATES*/
void mutation_rates(){
	int i,j;
	mutrate[0][1]=uf1;
	mutrate[0][2]=uf2;
	mutrate[0][3]=uf1*uf2;
	mutrate[0][0]=1-(mutrate[0][1]+mutrate[0][2]+mutrate[0][3]); //PP so mutation is parent-independent
	mutrate[1][0]=ub1;
	mutrate[1][2]=ub1*uf2;
	mutrate[1][3]=uf2;
	mutrate[1][1]=1-(mutrate[1][0]+mutrate[1][2]+mutrate[1][3]);
	mutrate[2][0]=ub2;
	mutrate[2][1]=ub2*uf1;
	mutrate[2][3]=uf1;
	mutrate[2][2]=1-(mutrate[2][0]+mutrate[2][1]+mutrate[2][3]);
	mutrate[3][0]=ub1*ub2;
	mutrate[3][1]=ub2;
	mutrate[3][2]=ub1;
	mutrate[3][3]=1-(mutrate[3][0]+mutrate[3][1]+mutrate[3][2]);	
	/*
	printf("mutrate\n");
	for (i=0; i<4; i++) {
		for (j=0; j<4; j++) 
			printf("%f  ",mutrate[i][j]);
		printf("\n");
	}*/
}

/*FUNCTION TO RUN ONE SIMULATION*/
struct outsim one_run(){ 
	float sumc,vltot; //vltot average viral load per day // sumc is used to add mut and mig rates. 
	int i,j,k,l,Ntotal,ii,jj,first,ddc; //PP Ntotal is the total pop size
	double tplus,random_rep; //tplus is the time till the next event. random_rep is a random uniform number
	struct nev nev_rec; //PP next event structure
	struct outsim results; //PP results of this simulation
	
	/*Initial conditions for every simulation*/
	for (i=0; i<4; i++)
		for (j=0; j<4; j++) 
			comp[i][j]=0; //PP all compartments are empty except the sanctuary (next line)
	comp[0][0]= comp[4][0]*dx*(R-1)/R; //Is this the starting pop size in the sanctuary? PP this is x0*dx*R0-1/R0*dy, so that the sim starts with the equilibrium pop size in the sanctuary. Maybe shoould be replaced by comp[0][0] = x0*dx*(R-1)/(R*dy);
//	printf("%lf ", comp[0][0]);
//	random_rep=gsl_rng_uniform (r); 
//if (random_rep<0.5) 
//		comp[1][0]=1;
//	else
//		comp[1][0]=2;
	comp[1][0] = gsl_ran_poisson (r, comp[0][0]*uf1/s);
//	printf("%lf ", comp[1][0]);
//	printf("\n");

	if (outputmutantinfo==1){
	printf("%d\t%d\t%d\t",i,j,k);
	printf("%f\t",t);
	printf("%f\t",comp[0][0]); // output num wt
	printf("%f\n",comp[1][0]);
	}
	
	Ntotal=comp[0][0]+comp[1][0];	//PP Ntotal is the total pop size	

	vltot=0.0; //average viral load per day
	t=0.0; //time
	first=0; 
	results.time_exit=0.0;
	results.fcompf=0;
	results.ssdc=-1;
	results.fcc=-1;
	
	/*START SIMULATION*/
	
//	while (t<1000){ // !! PP run for fixed amount of time, in stead of until invasion. 
	while ((comp[1][1]<tcf && rununtilfailureOrngens==0)|| t<rununtilfailureOrngens &&rununtilfailureOrngens>0){ //as long as the fourth compartment has fewer individuals than tcf (the threshhold). 
			/*
		printf("comp\n");
		for (ii=0; ii<4; ii++) {
			for (jj=0; jj<4; jj++) 
				printf("%f  ",comp[ii][jj]);
			printf("\n");
		}*/
		
		if (comp[3][3]==0) // used to look for the first succesful colonizer of the DDC
			ddc=0;
		
		/*CHOOSE NEXT EVENT*/
		nev_rec=next_event();
		j=nev_rec.j; //PP I reverted i and j 
		i=nev_rec.i; 
		k=nev_rec.k;
		if (i==1 & j==0 & outputmutantinfo==1){ //if something happens to a mutant
			printf("%d\t%d\t%d\t",i,j,k);
			printf("%f\t",t);
			printf("%f\t",comp[0][0]); // output num wt
			printf("%f\n",comp[1][0]);
		}
//			printf(" cc %f %d %d %d\n",t, i,j,k);
		
		/*EXECUTE EVENT*/
		/*Growth*/
		if(k==0){
			random_rep=gsl_rng_uniform (r);
			sumc=0.0;
			for (l=0; l<4; l++) { //loops over all strains that it can mutate to. 
				sumc+=mutrate[i][l]; 
				if (random_rep<sumc) { // PP the mutation rates add up to 1, so random_rep doesn't need to be multiplied with total mut rate. 
					if (j==0 & l==1 & i!=1 & outputmutantinfo==1){ //if new mutant born
						printf("%d\t%d\t%d\t",i,j,k);
						printf("%f\t",t);
						printf("%f\t",comp[0][0]); // output num wt
						printf("%f\n",comp[1][0]);
					}					
					
					comp[l][j]+=1; //add the mutant
					Ntotal+=1; //PP update total pop size
										
					if (comp[3][3]==1 && ddc==0) { // used to look for the first succesful colonizer of the DDC 
						//PP should we take into account that colonization can happen by more than 1 ind at a time, by using >= 1? 
						results.fcc=i;	
						// origin of first cell that colonizes the double-drug compartment (values: -1,0,1,2,10,11,12) 
						//(-1: no treatment failure | 0,1,2: cell present in the DDC that mutates from wt,m1 or m2 to m12 | 10,11,12: m12 cell that migrates from the sanctuary, SDC1 or SDC2 to the DDC) // written to path_sucmut.txt 
						ddc=1; 
					};
				
					break;
				}
			}
		}			
		
		else {
			/*Death*/
			if (k==1){ 
				Ntotal-=1; //PP update total pop size
				comp[i][j]-=1;//PP remove the dead ind
			}
			else{ 
				/* Migration */
				random_rep=gsl_rng_uniform (r);
				sumc=0.0;
				for (l=0; l<4; l++) {//loops over all possible compartments the migrant can go to. 
					sumc+=mig[j][l]/m; // normalize the migration rates by dividing by total mig rate. 
					if (random_rep<sumc) {
						comp[i][j]-=1;
						comp[i][l]+=1;
						
						if (comp[3][3]==1 && ddc==0) { // used to look for the first succesful colonizer of the DDC
							results.fcc=10+j;	
							// origin of first cell that colonizes the double-drug compartment (values: -1,0,1,2,10,11,12) 
							//(-1: no treatment failure | 0,1,2: cell present in the DDC that mutates from wt,m1 or m2 to m12 | 10,11,12: m12 cell that migrates from the sanctuary, SDC1 or SDC2 to the DDC) // written to path_sucmut.txt 
							ddc=1;
							}	
						if (i==1 && j==0) {// used to count migration events of a mutant from the sanctuary j==0 is in sanctuary, i ==1 is mutant strain
								mmt+=1;	
							if (outputmutantinfo==0) {printf("%lf  ", t); } 
							}
						break;
					}
				}
			}
		}
		
		/*UPDATE TIME*/
		//random_rep=gsl_rng_uniform (r); //PP I don't understand why we need a new random_rep? 
		tplus=gsl_ran_exponential(r, 1/(nev_rec.rtot)); //PP determine a random time when the event happens 
		t+=tplus; //PP update t with the random time. 
		
		/*Average viral load*/
		vltot+=(float) Ntotal*tplus; //average viral load per day (will later be divided by time)
		
		/*Exit of the sanctuary*/
		if (first==0 && comp[1][1]>tcf) {
			first=1;
			results.fcompf=1;
			results.time_exit=t;
		}
		else {
			if (first==0 && comp[2][2]>tcf){ 
				first=1;
				results.fcompf=2;
				results.time_exit=t;
			}
			else {
				if (first==0 && comp[3][3]>tcf){ 
					first=1;
					results.fcompf=3;
					results.time_exit=t;
				}
			}
			
		};
		
	
	};
	
	/*END SIMULATION*/
	
	
	/****Save and return results*****/
	
	/*time to treatment failure*/
	results.time_failure=t;
	
	/*average viral load*/
	results.vltot= (float) vltot/t; //average viral load per day
	
	/*state of the SDC when the DDC is colonized*/	
	if (comp[3][3]>=tcf-1) {
		if (comp[1][1]>tcf) {
			if (comp[2][2]>tcf)
				results.ssdc=3;
			else
				results.ssdc=1;
		}
		else {
			if (comp[2][2]>tcf) 
				results.ssdc=2;
			else 
				results.ssdc=0;
		}	
	};	
	
	return results;
}


/*FUNCTION TO DETERMINE RATE OF ALL POSSIBLE EVENTS*/
struct nev next_event(){ // returns type of event that will happen and total rate of events. 
	float rates[4][4][3]={0}; // i: type of strain, j: compartment, k: event (0:growth, 1:death, 2:migration)
	float sumr[4]={0},x0; //PP sum of rates for each compartment  
	int i,j,k;
	int jik[3]; 
	struct nev nev_result;
	double rtot; 
	rtot=0.0;
	
	
	/*Create matrix with the rate of all the possible events*/
	for (j=0; j<4; j++){ //PP loop over the compartments
		sumr[j]=0; 
		/*Sum term for calculating rate of growth in compartment j*/
		x0=comp[4][j]; //total num inds in comp j
		for (i=0; i<4; i++) // loop over all strains 
			sumr[j]+=R0[i][j]*dy*comp[i][j]; //add the R0s* dy* numinds
		/*Fill in matrix with rates*/
		for (i=0; i<4; i++){ //loop again over all strains
			if (comp[i][j]!=0) //if the strain is present 
				for(k=0; k<3; k++){ // loop over possible life events
					if (k==0)  // growth 
						rates[i][j][k]=comp[i][j]*((dy*x0*dx*R0[i][j])/(x0*dx + sumr[j])); //PP growth rate is numinds * dy*x0*dx*R0/(x0*dx+totalrate)
					else {
						if (k==1)  // death
							rates[i][j][k]=dy*comp[i][j]; //PP death rate is just numinds*dy
						else
							rates[i][j][k]=m*comp[i][j]; //migration
					}
				rtot+=rates[i][j][k]; //PP keep track of total rate
				}	
		}

	};

	/*Choose event that will happen*/
	choose(rtot,rates,jik); //PP call function choose to determine which event will happen next 

	nev_result.j=jik[0]; //PP should be j = jik[0], no? I reverted j and i // which compartment is the event? 
	nev_result.i=jik[1]; //PP which strain? 
	nev_result.k=jik[2]; //PP what kind of event? growth, death, migration. 
	nev_result.rtot=rtot; //total rate
	
	return nev_result;
}


/*FUNCTION TO CHOOSE NEXT EVENT*/
int choose(double rtot,float rates[4][4][3],int jik[]){
	double random_rep; // PP will be random uniform number
	double cumr; //PP cumulative rate, add rates and stop when cumr > random_rep*rtot 
	int i,j,k;
	
	/*
	int ii,jj; 
	printf("rates for cells of type 0\n");
	for (ii=0; ii<4; ii++) {
		for (jj=0; jj<3; jj++) 
			printf("%f  ",rates[ii][0][jj]);
		printf("\n");
	}
	*/

	random_rep=gsl_rng_uniform (r);
	cumr=0.0;
	for (j=0; j<4; j++)
		for (i=0; i<4; i++)
			for (k=0; k<3; k++){
				cumr+=rates[i][j][k];
				//printf("cumr: %f, %d %d %d\n",cumr,j,i,k);
				if (random_rep*rtot<cumr){ 
					jik[0]=j;
					jik[1]=i;
					jik[2]=k;
					//printf("%d %d %d \n",i,j,k);
					return 0; //PP once the function sees the return statement it goes out of the function
				}
			} // PP does this function not need an exit statement? Or is that implicit in c? It exits when it hits the return statement  				
	return 0;
}


