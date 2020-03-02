/************************************************************************************************************
 * Individual-based model for coevolution of coordination and proportion of helpers                         *
 * For manuscript "Dividing labour in social organisms: coordinated or random specialisation?" by Guy       *
 * Alexander Cooper, Jorge Peña, Ming Liu, Stuart Andrew West                                               *
 * This file runs simulations with repeated time series                                                     *
 * Average genotypes and phenotypes across repeats are logged every 10 generations in the output files      *
 ************************************************************************************************************
 * Key parameters                                                                                           *
 * s1: Number of evolving trait. 0: no coordination evolving, 1: coordination evolves, 2: coordination and  *
 *  target propotion of helping evolve                                                                      *
 * s2: Specifies the function of coordination cost. 0: function I (as in Fig. S6), 1: function IIa,         *
 *  2: function IIb, 3: function IIc.                                                                       *
 * n: group size                                                                                            *
 * e: essentiality of helping                                                                               *
 * T, mut_rate: Number of generations and mutation rate                                                     *
 * N_total, N_rep: Population size and number of repeats for each parameter settings                        *
 ************************************************************************************************************/

// Include standard tool files
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include "../dSFMT-src-2.2.3/dSFMT.c"

// Global functions
    double Mean_array(double p[], int length);                                  // Average of an array
    double sum(double p[], int length);                                         // Sum of an array
    int RandFromProb(double p[], int length, double RandNum);                   // Random seqence generator wrighted by the values of elements in the array
    double normal_dist_BM (double mean, double sd, double u1, double u2);       // Random number generator for normal distribution
    double **d_array2d(long size_1, long size_2);                               // Create a 2d array (stack memory)
    void free_d_array2d(double **x);                                            // Clear the 2d array
    double ****d_array4d(long size_1, long size_2, long size_3, long size_4);   // Create a 4d array (stack memory)
    void free_d_array4d(double ****x);                                          // Clear the 4d array

int main(){
    // General parameters
        int T= 30000;                   // Number of generations
        int N_rep= 10;                  // Number of repetitions
        int n= 18;                      // Number of individual in each social group
        int N_total= 10000;             // Total number of individuals
        int G= ceil((double)N_total/n); // Number of groups
        double e= 0.95;                 // Essentiality of helping, related to fitness function
        double p_star= (2*e-1)/(2*e);   // Optimal probability of being a helper
        int num_alloc=  5*n;            // Number of allocating types
        double mut_rate= 0.001;         // Mutation rate of s
        double MutStep= 0.1;            // Size of mutation for S
        double d= 0.0002;               // Scaling coefficient of coordination cost (for s2==0, 2, 3)
        double alpha= 0.1;              // Scaling coefficient of coordination cost (for s2==1)
        double beta= 0.01;              // Nonlinearity coefficient for s2==1
        double psi= 0.1;                // Nonlinearity coefficient for s2==2
        double phi= 0.001;              // Nonlinearity coefficient for s2==3

    // Switches
        int s1= 2;              
            // 0: no coordination (s==0), 1: coordination level (s) evolves
            // 2: coordination level and target proportion of helpers (s and p) coevolve
            // 4: full coodination (s==1)
        int s2= 1;              
            // 0: Linear function, cost= d*y, y=s*n*(n-1)/2 
            // 1: Exponential function, cost= alpha*(1-exp(-beta*y))
            // 2: Another exponential function, cost= psi*(d*y)^psi
            // 3: Type II functional response, cost= d*y/(1+phi*y) 

    // Log variables
        double p_i;             // Proportion of interacting helpers
        int n_i;                // Number of interacting neughbours
        int k_i;                // Number of interacting helpers

    // Output settings
        FILE *out, *out2, *out3, *out4;
        out= fopen("out_s_10.txt","w");
        fprintf(out, "T\ts_1\ts_2\ts_3\ts_4\ts_5\ts_6\ts_7\ts_8\ts_9\ts_10\n");
        out2= fopen("out_p_10.txt","w");
        fprintf(out2, "T\tp_1\tp_2\tp_3\tp_4\tp_5\tp_6\tp_7\tp_8\tp_9\tp_10\n");
        out3= fopen("out_porp_10.txt","w");
        fprintf(out3, "T\tp_1\tp_2\tp_3\tp_4\tp_5\tp_6\tp_7\tp_8\tp_9\tp_10\n");
        out4= fopen("out_cord_10.txt","w");
        fprintf(out4, "T\tp_1\tp_2\tp_3\tp_4\tp_5\tp_6\tp_7\tp_8\tp_9\tp_10\n");

    // Random number genertor
        int seed;
        dsfmt_t dsfmt;
        seed= time(NULL);
        if(seed==0)seed= 1;
        dsfmt_init_gen_rand(&dsfmt,seed);
        int i,j,k,l;
        int t, rep, idx, k_num, g_samp, int_tmp;
        double pp, pp2, pp3, u1, u2, social_temp, rho, s_temp, p_temp, cost, mut_size;
        double idx_list[n];
        for(i=0; i<n; i++) idx_list[i]= i;
    
    // Temporal space
        double ****Popn= d_array4d(N_rep, G, n, 3);     // The current population
        double ****PopnNext= d_array4d(N_rep, G, n, 3); // The population of next generation, sampled from Popn according to group fitness
        double **ProbGrup= d_array2d(N_rep,G);          // Relative fitness of each group
        double **PorpHelp= d_array2d(N_rep, G);         // Proportion of helpers present in the population (phenotypic)
        double **PorpCord= d_array2d(N_rep, G);         // Level of coordination in the population (phenotypic)
        double **link= d_array2d(n,n);                  // The presence of interaction
        double tmp[G];                  // Temporary array to store a property of the group members within the focal social group
        double idvl_tmp[n];             // Temporary array for individuals in a group
        double globalS[N_rep];          // Population average of S, the genotype of coordination
        double globalP[N_rep];          // Population average of P, the genotype of becoming a helper
        double globalPorp[N_rep];       // Population average of Porp, the phenotype of helping
        double globalCord[N_rep];       // Population average of Cord, the phenotype of coordination
    // Indexes of the individual
        int id_s=   0;                  // Genotype, probability of interacting with group members
        int id_p=   1;                  // Genotype, probability of becoming a helper
        int id_type=2;                  // Phenotype, actual division of labor, 1= helper, 0= pure reproductive

    // Initialization
    t=0;
    for(rep=0; rep<N_rep; rep++){
        globalS[rep]= 0.0;
        for(i=0; i<G; ++i){
            for(j=0; j<n; ++j) {
                Popn[rep][i][j][id_s]=      0.0;
                if(s1==4) Popn[rep][i][j][id_s]= 1.0;
                Popn[rep][i][j][id_p]=      0.0;
                if(s1==1) Popn[rep][i][j][id_p]= p_star;
                Popn[rep][i][j][id_type]=   0.0;
    }}}
    // Start of simulation
    for(t=1; t<=T; ++t){
        for(rep=0; rep<N_rep; rep++){
            // Group level events
            for(i=0; i< G; i++){
                // Reset the links
                for(j=0; j<n; j++){
                    for(k=0; k<n; k++){
                        if(j==k) link[j][k]= NAN;
                        if(j>k){
                            pp= dsfmt_genrand_open_close(&dsfmt);
                            if(pp< Popn[rep][i][0][id_s]) link[j][k]= 1.0;  // link present
                            else link[j][k]= 0.0;                           // link absent
                            link[k][j]= link[j][k];
                }}}
                // Determine division (being a helper or reproductive)
                // Evolving both S and P
                if(s1==2 || s1==4){
                    for(j=0; j<num_alloc; j++){
                        n_i= k_i= 0;
                        // pick an individual randomly
                        pp= dsfmt_genrand_open_close(&dsfmt);
                        idx= floor(pp*n);
                        // calculate the probability of being a helper from all connected neughbours
                        for(k=0; k<n; k++){
                            if(idx!=k && link[idx][k]== 1.0){
                                n_i+= 1;
                                if(Popn[rep][i][k][id_type]== 1.0) k_i+= 1;
                            }
                        }
                        // Conditional allocation
                        if(n_i> 0){
                            p_i= (double)k_i/n_i;
                        // If actual proportion is higher than its genotype -> be a reproductive
                            if(Popn[rep][i][idx][id_p]< p_i) Popn[rep][i][idx][id_type]= 0.0;
                            // If actual proportion is lower than its genotype -> be a helper
                            if(Popn[rep][i][idx][id_p]> p_i) Popn[rep][i][idx][id_type]= 1.0;
                            if(Popn[rep][i][idx][id_p]== p_i){
                                if(dsfmt_genrand_open_close(&dsfmt)< 0.5) \
                                    Popn[rep][i][idx][id_type]= 1.0;
                                else Popn[rep][i][idx][id_type]= 0.0;
                            }
                        }
                        else{
                            if(dsfmt_genrand_open_close(&dsfmt)< Popn[rep][i][idx][id_p]) \
                                Popn[rep][i][idx][id_type]= 1.0;
                            else Popn[rep][i][idx][id_type]= 0.0;
                        }
                    }
                }
                // Evolving S based on optimal P (p*)
                if(s1==1){
                    for(j=0; j<num_alloc; j++){
                        n_i= k_i= 0;
                        // pick an individual randomly
                        pp= dsfmt_genrand_open_close(&dsfmt);
                        idx= floor(pp*n);
                        // calculate the probability of being a helper from all connected neughbours
                        for(k=0; k<n; k++){
                            if(idx!=k && link[idx][k]== 1.0){
                                n_i+= 1;
                                if(Popn[rep][i][k][id_type]== 1.0) k_i+= 1;
                            }
                        }
                        // Conditional allocation
                        if(n_i>0){
                            p_i= (double)k_i/n_i;
                            // If actual proportion is higher than the optimal -> be a reproductive
                            if(p_star< p_i) Popn[rep][i][idx][id_type]= 0.0;
                            // If actual proportion is lower than the optimal -> be a helper
                            if(p_star> p_i) Popn[rep][i][idx][id_type]= 1.0;
                            if(p_star== p_i){
                                if(dsfmt_genrand_open_close(&dsfmt)< 0.5) Popn[rep][i][idx][id_type]= 0.0;
                                else Popn[rep][i][idx][id_type]= 1.0;
                            }
                        }
                        else{
                            if(dsfmt_genrand_open_close(&dsfmt)< p_star) Popn[rep][i][idx][id_type]= 1.0;
                            else Popn[rep][i][idx][id_type]= 0.0;
                        }
                }}
                // No coordination at all (S fixed at 0)
                if(s1==0){
                    for(j=0; j<n; j++){
                        pp= dsfmt_genrand_open_close(&dsfmt);
                        if(pp< Popn[rep][i][j][id_p]) Popn[rep][i][j][id_type]= 1.0;
                        else Popn[rep][i][j][id_type]= 0.0;
                    }
                }
                // Fitness calculation
                    k_num= 0;
                    for(j=0; j<n; j++){
                        if(Popn[rep][i][j][id_type]== 1.0) k_num+= 1;
                    }
                    // Define the cost of coordination
                    if(s2==0) cost= d*Popn[rep][i][0][id_s]*n*(n-1)/2;
                    if(s2==1) cost= alpha*(1- exp(Popn[rep][i][0][id_s]*(-1)*beta*n*(n-1)/2));
                    if(s2==2) cost= psi*pow(d*Popn[rep][i][0][id_s]*n*(n-1)/2, psi);
                    if(s2==3) cost= d*(Popn[rep][i][0][id_s]*n*(n-1)/2)/ \
                        (1+ phi*Popn[rep][i][0][id_s]*n*(n-1)/2);
                    // Boundary condition for s2!= 1
                    if(cost>1) cost=1;
                    // Calculation
                    ProbGrup[rep][i]= (n- k_num)*(1-cost)*((1-e)+e*k_num/n);
                    PorpHelp[rep][i]= (double)k_num/n;
                // Actual level of coordination
                    for(j=0; j<n; j++){
                        n_i= 0;
                        for(k=0; k<n; k++){
                            if(k!=j && link[j][k]== 1.0) n_i+= 1;
                        }
                        idvl_tmp[j]= (double)n_i/(n-1);
                    }
                    PorpCord[rep][i]= Mean_array(idvl_tmp, n);
            }
            // Mutation
            for(i=0; i<G; i++){
                // Sample an individual to form a group
                g_samp= RandFromProb(ProbGrup[rep], G, dsfmt_genrand_open_close(&dsfmt));
                // Mutation of s
                if(s1==0 || s1==4) mut_size= 0.0;
                else{
                    if(dsfmt_genrand_open_close(&dsfmt)< mut_rate){
                        mut_size= MutStep*normal_dist_BM(0,1,dsfmt_genrand_open_close(&dsfmt),dsfmt_genrand_open_close(&dsfmt));
                    }
                    else mut_size= 0.0;
                }
                s_temp= Popn[rep][g_samp][0][id_s]+ mut_size;
                // Mutation of p
                if(dsfmt_genrand_open_close(&dsfmt)< mut_rate){
                    mut_size= MutStep*normal_dist_BM(0,1,dsfmt_genrand_open_close(&dsfmt),dsfmt_genrand_open_close(&dsfmt));
                }
                else mut_size= 0.0;
                p_temp= Popn[rep][g_samp][0][id_p]+ mut_size;
                // Dealing with the extreme conditions
                    if(s_temp> 1) s_temp= 1.0;
                    if(s_temp< 0) s_temp= 0.0;
                    if(p_temp> 1) p_temp= 1.0;
                    if(p_temp< 0) p_temp= 0.0;
                // Clonal group, filling all group member for the next generation with the same values
                for(j=0; j<n; j++){
                    PopnNext[rep][i][j][id_s]= s_temp;
                    if(s1==4) PopnNext[rep][i][j][id_s]= 0.0;
                    PopnNext[rep][i][j][id_p]= p_temp;
                    PopnNext[rep][i][j][id_type]= 0.0;
                }
            // End of a group
            }
        // End of a population
        }
        for(rep=0; rep<N_rep; rep++){
            for(i=0; i<G; i++){
                for(j=0; j<n; j++){
                    for(k=0; k<3; k++) Popn[rep][i][j][k]= PopnNext[rep][i][j][k];
            }}
            for(i=0; i<G; i++) tmp[i]= Popn[rep][i][0][id_s];
            globalS[rep]= Mean_array(tmp, G);
            for(i=0; i<G; i++) tmp[i]= Popn[rep][i][0][id_p];
            globalP[rep]= Mean_array(tmp, G);
            for(i=0; i<G; i++) tmp[i]= PorpHelp[rep][i];
            globalPorp[rep]= Mean_array(tmp, G);
            for(i=0; i<G; i++) tmp[i]= PorpCord[rep][i];
            globalCord[rep]= Mean_array(tmp, G);
        }
        // Print
        if(t% 10==0){
            fprintf(out,"%d\t",t);
            for(rep=0; rep<N_rep; rep++){
                fprintf(out, "%lf\t", globalS[rep]);
                if(rep==(N_rep-1)) fprintf(out, "\n");
            }
            fprintf(out2,"%d\t",t);
            for(rep=0; rep<N_rep; rep++){
                fprintf(out2, "%lf\t", globalP[rep]);
                if(rep==(N_rep-1)) fprintf(out2, "\n");
            }
            fprintf(out3,"%d\t",t);
            for(rep=0; rep<N_rep; rep++){
                fprintf(out3, "%lf\t", globalPorp[rep]);
                if(rep==(N_rep-1)) fprintf(out3, "\n");
            }
            fprintf(out4,"%d\t",t);
            for(rep=0; rep<N_rep; rep++){
                fprintf(out4, "%lf\t", globalCord[rep]);
                if(rep==(N_rep-1)) fprintf(out4, "\n");
            }
        }
    }
    free_d_array2d(link);
    free_d_array2d(ProbGrup);
    free_d_array2d(PorpHelp);
    free_d_array2d(PorpCord);
    free_d_array4d(Popn);
    free_d_array4d(PopnNext);
    fclose(out);
    fclose(out2);
    fclose(out3);
    fclose(out4);
    return 0;
}

////////////////////// Functions are defined in below //////////////////////
double Mean_array(double p[], int length){
    int i;
    double temp, sum;
    temp= sum= 0.0;
    for (i=0; i<length; ++i) sum+= p[i];
    if(length>0) temp= sum/length;
    else temp= NAN;
    return temp;
}

double sum(double p[], int length){
    int i;
    double sum=0.0;
    for(i=0; i< length; ++i){
        if(p[i]> 0.0 && sum> (DBL_MAX- p[i])) printf("Double overflow detected\n");
        else sum+= p[i];
    }
    return sum;
}

int RandFromProb(double p[], int length, double RandNum){
    int i,temp,check;
    double sum, pp;
    sum= pp= 0.0;
    for(i=0; i< length; ++i) sum+= p[i];
    if(sum> 0.0){
        for(i=0; i< length; ++i){
            if(RandNum> pp/sum) check=1;
            pp+= p[i];
            if(RandNum< pp/sum&& check==1){
                temp= i;
                i= length;
            }
            check=0;
        }
        return temp;
    }
    else{
        return floor(length*RandNum);
    }
}

double normal_dist_BM (double mean, double sd, double u1, double u2){
    // Using Box-Muller method to generate pseudo-normal distributed numbers in [0,1]
    double z1;
	z1= sqrt(-2* log(u1))* cos(2* M_PI* u2);
	return z1*sd+ mean;
}

double **d_array2d(long size_1, long size_2){
    double **x;
    long i;
    x= (double **) malloc((size_t)(size_1*sizeof(double)));
    for(i=0; i< size_1; i++) x[i]= (double *) malloc((size_t)(size_2*sizeof(double)));
    return x;
}

void free_d_array2d(double **x){
    free(x[0]);
    free(x);
}

double ****d_array4d(long size_1, long size_2, long size_3, long size_4){
    double ****x;
    long i,j,k;
    x= (double ****) malloc((size_t)(size_1*sizeof(double ***)));
    for(i=0; i< size_1; i++) {
        x[i]= (double ***) malloc((size_t)(size_2*sizeof(double **)));
        for(j=0; j< size_2; j++){
            x[i][j]= (double **) malloc((size_t)(size_3*sizeof(double *)));
            for(k=0; k< size_3; k++) x[i][j][k]= (double *) malloc((size_t)(size_4*sizeof(double)));
        }
    }
    return x;
}

void free_d_array4d(double ****x){
    free(x[0][0][0]);
    free(x[0][0]);
    free(x[0]);
    free(x);
}
