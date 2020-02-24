/************************************************************************************************************
 * Individual-based model for coevolution of coordination and proportion of helpers                         *
 * For manuscript "Dividing labour in social organisms: coordinated or random specialisation?" by Guy       *
 * Alexander Cooper, Jorge Peña, Ming Liu, Stuart Andrew West                                               *
 * This file runs simulations with fixed group size, fixed number of repetitions, and variable              *
 * essentialities of helping.                                                                               *
 * Average genotypes and phenotypes across repeats are logged in the output file (summary.txt)              *
 ************************************************************************************************************
 * Key parameters                                                                                           *
 * s1: Number of evolving trait. 0: no coordination evolving, 1: coordination evolves, 2: coordination and  *
 *  target propotion of helping evolve                                                                      *
 * s2: Specifies the function of coordination cost. 0: function I (as in Fig. S6), 1: function IIa,         *
 *  2: function IIb, 3: function IIc.                                                                       *
 * n: group size                                                                                            *
 * T, mut_rate: Number of generations and mutation rate                                                     *
 * N_total, N_rep: Population size and number of repeats for each parameter settings                        *
 ************************************************************************************************************/

// Include standard tool files
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "../dSFMT-src-2.2.3/dSFMT.c"

// Global functions
    double Mean_array(double p[], int length);                                  // Average of an array
    double SD_array(double p[], int length);                                    // SD of an array
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
        int N_rep= 1;                  // Number of repetitions
        // int G= 1000;                 // Number of groups
        int n= 20;                      // Number of individual in each social group
        int N_total= 10000;             // Total number of individuals
        int G= ceil((double)N_total/n); // Number of groups
        double e= 1.0;                  // Trait essentiality, related to fitness return function
        double p_star= (2*e-1)/(2*e);   // Optimal probability of being a helper
        int num_alloc=  5*n;            // Number of allocating divisions
        double mut_rate= 0.001;         // Mutation rate
        double MutStep= 0.1;            // Size of mutation
        double d= 0.0002;               // Scaling coefficient of coordination cost (for s2==0, 2, 3)
        double alpha= 0.1;              // Scaling coefficient of coordination cost (for s2==1)
        double beta= 0.01;              // Nonlinearity coefficient for s2==1
        double psi= 0.1;                // Nonlinearity coefficient for s2==2
        double phi= 0.001;              // Nonlinearity coefficient for s2==3

    // Switches
        int s1= 2;              
            // 0: no coordination (s==0), 1: coordination level (s) evolves
            // 2: coordination level and target proportion of helpers (s and p) coevolve
        int s2= 1;              
            // 0: Linear function, cost= d*y, y=s*n*(n-1)/2 
            // 1: Exponential function, cost= alpha*(1-exp(-beta*y))
            // 2: Another exponential function, cost= psi*(d*y)^psi
            // 3: Type II functional response, cost= d*y/(1+phi*y) 

    // Generate the list of essentialities
        int i_e;
        int length_param= 20;
        double esseniality[length_param];
        for(i_e=0; i_e< length_param; i_e++) esseniality[i_e]= 0.5+ i_e*(0.5)/(double)(length_param-1);

    // Output temporal space and settings
        int length_log= 0.1*T;
        // Temporal space for the last length_log generaions (the averages)
        double log_avg_s[length_log], log_asd_s[length_log], log_avg_p[length_log], log_asd_p[length_log], log_avg_porp[length_log], log_asd_porp[length_log], log_avg_cord[length_log], log_asd_cord[length_log];
        // Temporal space for the last length_log generaions (the standard deviaitons)
        double log_sd_avg_s[length_log], log_sd_asd_s[length_log], log_sd_avg_p[length_log], log_sd_asd_p[length_log], log_sd_avg_porp[length_log], log_sd_asd_porp[length_log], log_sd_avg_cord[length_log], log_sd_asd_cord[length_log];
        // Creating the output file
        FILE *out;
        out= fopen("summary.txt", "w");
        fprintf(out, "essentiality\tgroup_size\tavg_s\tasd_s\tavg_p\tasd_p\tavg_porp\tasd_porp\tavg_cord\tasd_cord\tsd_avg_s\tsd_asd_s\tsd_avg_p\tsd_asd_p\tsd_avg_porp\tsd_asd_porp\tsd_avg_cord\tsd_asd_cord\n");
        // average s over repeats, average sd of s over repeats, average p over repeats, average sd of p over repeats, sd of average s, sd of sd of s, sd of average p, sd of sd of p.

    // Random number genertor initialization
        int seed;
        dsfmt_t dsfmt;
        seed= time(NULL);
        if(seed==0)seed= 1;
        dsfmt_init_gen_rand(&dsfmt,seed);

    // Variables
        int i,j,k;
        int t, rep, idx, k_num, g_samp;
        double pp, pp2, pp3, u1, u2, social_temp, rho, s_temp, p_temp, cost, avg_S, asd_S, avg_P, asd_P, avg_porp, asd_porp, avg_cord, asd_cord, sd_avg_S, sd_avg_P, sd_asd_S, sd_asd_P, sd_avg_porp, sd_asd_porp, sd_avg_cord, sd_asd_cord, mut_size;
        double p_i;                     // proportion of interacting helpers
        int n_i;                        // number of interacting neughbours
        int k_i;                        // number of interacting helpers
    
    // Temporal space
        double ****Popn= d_array4d(N_rep, G, n, 3);     // The current population
        double ****PopnNext= d_array4d(N_rep, G, n, 3); // The population of next generation, sampled from Popn according to group fitness
        double **ProbGrup= d_array2d(N_rep,G);          // Relative fitness of each group
        double **PorpHelp= d_array2d(N_rep, G);         // Proportion of helpers present in the population (phenotypic)
        double **PorpCord= d_array2d(N_rep, G);         // Level of coordination in the population (phenotypic)
        double link[n][n];              // Interaction matrix of a group
        double tmp[G];                  // Temporary array for groups in a population
        double idvl_tmp[n];             // Temporary array for individuals in a group
        double globalS[N_rep];          // Population average of S, the genotype of coordination
        double globalS_sd[N_rep];       // Individual-level standard deivation of S of the population
        double globalP[N_rep];          // Population average of P, the genotype of becoming a helper
        double globalP_sd[N_rep];       // Individual-level standard deivation of P of the population
        double globalCord[N_rep];       // Population average of Cord, the phenotype of coordination
        double globalCord_sd[N_rep];    // Individual-level standard deivation of Cord of the population
        double globalPorp[N_rep];       // Population average of Porp, the phenotype of helping
        double globalPorp_sd[N_rep];    // Individual-level standard deivation of Porp of the population

    // Indexes of genotypes and phenotype
        int id_s=   0;                  // Genotype, probability of interacting with group members
        int id_p=   1;                  // Genotype, probability of becoming a helper
        int id_type=2;                  // Phenotype, actual division of labor, 1= helper, 0= pure reproductive

    // Initialization
    for(i_e=0; i_e< length_param; i_e++){
        e= esseniality[i_e];
        p_star= (2*e-1)/(2*e);          // Optimal proportion of helpers, used when s1==1
        t=0;
        // Initializing population
        for(rep=0; rep<N_rep; rep++){
            globalS[rep]= 0.0;
            for(i=0; i<G; ++i){
                for(j=0; j<n; ++j) {
                    Popn[rep][i][j][id_s]=      0.0;
                    Popn[rep][i][j][id_p]=      0.0;
                    Popn[rep][i][j][id_type]=   0.0;
        }}}
        // Start of simulation 
        for(t=1; t<T; ++t){
            for(rep=0; rep<N_rep; rep++){
                // Group level events: linkage, allocation, fitness, and the actual propotion of helpers
                for(i=0; i< G; i++){
                    // Reset the links
                    for(j=0; j<n; j++){
                        for(k=0; k<n; k++){
                            if(j==k) link[j][k]= NAN;
                            if(j>k){
                                pp= dsfmt_genrand_open_close(&dsfmt);
                                if(pp< Popn[rep][i][0][id_s]) link[j][k]= 1.0;   // link present
                                else link[j][k]= 0.0;                       // link absent
                                link[k][j]= link[j][k];
                    }}}
                    // Determine the divisions (being a reproductive or a helper)
                    // Evolving both S and P
                    if(s1==2){
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
                    }}
                    // Evolving S based on optimal P (p*)
                    if(s1==1){
                        for(j=0; j<num_alloc; j++){
                            n_i= k_i= 0;
                            pp= dsfmt_genrand_open_close(&dsfmt);
                            idx= floor(pp*n);
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
                                if(dsfmt_genrand_open_close(&dsfmt)< p_star)    \
                                    Popn[rep][i][idx][id_type]= 1.0;
                                else Popn[rep][i][idx][id_type]= 0.0;
                            }
                    }}
                    // No coordination at all (S fixed at 0)
                    if(s1==0){
                        for(j=0; j<n; j++){
                            if(dsfmt_genrand_open_close(&dsfmt)< Popn[rep][i][j][id_p]) \
                                Popn[rep][i][j][id_type]= 1.0;
                    }}
                    // Fitness calculation
                        // Get the number of helpers in the focal group
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
                    // Actual level of coordination
                        if(t> (T-length_log-1)){
                            PorpHelp[rep][i]= (double)k_num/n;
                            for(j=0; j<n; j++){
                                n_i= 0;
                                for(k=0; k<n; k++){
                                    if(k!=j && link[j][k]== 1.0) n_i+= 1;
                                }
                                idvl_tmp[j]= (double)n_i/(n-1);
                            }
                            PorpCord[rep][i]= Mean_array(idvl_tmp, n);
                        }
                }
                // Mutation
                for(i=0; i<G; i++){
                    // Sample an individual to form a group
                    g_samp= RandFromProb(ProbGrup[rep], G, dsfmt_genrand_open_close(&dsfmt));
                    // Mutation of s
                    if(s1==0) mut_size= 0.0;
                    else{
                        if(dsfmt_genrand_open_close(&dsfmt)< mut_rate){
                            mut_size= MutStep*normal_dist_BM(0,1, \
                                dsfmt_genrand_open_close(&dsfmt),dsfmt_genrand_open_close(&dsfmt));
                        }
                        else mut_size= 0.0;
                    }
                    s_temp= Popn[rep][g_samp][0][id_s]+ mut_size;
                    // Mutation of p
                    if(dsfmt_genrand_open_close(&dsfmt)< mut_rate){
                        mut_size= MutStep*normal_dist_BM(0,1, \
                            dsfmt_genrand_open_close(&dsfmt),dsfmt_genrand_open_close(&dsfmt));
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
                        PopnNext[rep][i][j][id_p]= p_temp;
                        PopnNext[rep][i][j][id_type]= 0.0;
                    }
                // End of a group
                }
            // End of a population
            }
            // Calculate population values
            for(rep=0; rep<N_rep; rep++){
                for(i=0; i<G; i++){
                    for(j=0; j<n; j++){
                        for(k=0; k<3; k++) Popn[rep][i][j][k]= PopnNext[rep][i][j][k];
                }}
                if(t> (T-length_log-1)){
                    for(i=0; i<G; i++) tmp[i]= Popn[rep][i][0][id_s];
                    globalS[rep]= Mean_array(tmp, G);
                    globalS_sd[rep]= SD_array(tmp, G);
                    for(i=0; i<G; i++) tmp[i]= Popn[rep][i][0][id_p];
                    globalP[rep]= Mean_array(tmp, G);
                    globalP_sd[rep]= SD_array(tmp, G);
                    // Also calculate the actual divisions
                    for(i=0; i<G; i++) tmp[i]= PorpHelp[rep][i];
                    globalPorp[rep]= Mean_array(tmp, G);
                    globalPorp_sd[rep]= SD_array(tmp, G);
                    for(i=0; i<G; i++) tmp[i]= PorpCord[rep][i];
                    globalCord[rep]= Mean_array(tmp, G);
                    globalCord_sd[rep]= SD_array(tmp, G);
                }
            }
            if(t> (T-length_log-1)){
            // Average over repetitions of current generation
                log_avg_s[t-T+length_log]= Mean_array(globalS, rep);            // Average of population S
                log_asd_s[t-T+length_log]= Mean_array(globalS_sd, rep);         // Average of population S_sd
                log_avg_p[t-T+length_log]= Mean_array(globalP, rep);            // Average of population P
                log_asd_p[t-T+length_log]= Mean_array(globalP_sd, rep);         // Average of population P_sd
                log_avg_porp[t-T+length_log]= Mean_array(globalPorp, rep);      // Average of population Porp
                log_asd_porp[t-T+length_log]= Mean_array(globalPorp_sd, rep);   // Average of population Porp_sd
                log_avg_cord[t-T+length_log]= Mean_array(globalCord, rep);      // Average of population Cord
                log_asd_cord[t-T+length_log]= Mean_array(globalCord_sd, rep);   // Average of population Cord_sd
            // SD over repetitions of current genertaion
                log_sd_avg_s[t-T+length_log]= SD_array(globalS, rep);           // SD of population S
                log_sd_asd_s[t-T+length_log]= SD_array(globalS_sd, rep);        // SD of population S_sd
                log_sd_avg_p[t-T+length_log]= SD_array(globalP, rep);           // SD of population P
                log_sd_asd_p[t-T+length_log]= SD_array(globalP_sd, rep);        // SD of population P_sd
                log_sd_avg_porp[t-T+length_log]= SD_array(globalPorp, rep);     // SD of population Porp
                log_sd_asd_porp[t-T+length_log]= SD_array(globalPorp_sd, rep);  // SD of population Porp_sd
                log_sd_avg_cord[t-T+length_log]= SD_array(globalCord, rep);     // SD of population Cord
                log_sd_asd_cord[t-T+length_log]= SD_array(globalCord_sd, rep);  // SD of population Cord_sd
            }
            // Print final results
            if(t== (T-1)){
            // Average (of average over repetitions) over generations
                avg_S= Mean_array(log_avg_s, length_log);       // Average level of coordination 
                asd_S= Mean_array(log_asd_s, length_log);       // Average SD of coordination (SD over individuals)
                avg_P= Mean_array(log_avg_p, length_log);       // Average target proportion of helpers
                asd_P= Mean_array(log_asd_p, length_log);       // Average SD of target proportion
                avg_porp= Mean_array(log_avg_porp, length_log); // Average actual proportion of helpers
                asd_porp= Mean_array(log_asd_porp, length_log); // Average SD of actual proportion
                avg_cord= Mean_array(log_avg_cord, length_log); // Average actual level of coordination
                asd_cord= Mean_array(log_asd_cord, length_log); // Average SD of actual level of coordination
            // Average (of SD over repetitions) over genertaions
                sd_avg_S= Mean_array(log_sd_avg_s, length_log);         // Average SD of average level of coordination (SD over repetitions)
                sd_asd_S= Mean_array(log_sd_asd_s, length_log);         // Average SD of SD of coordination (SD over repetitions)
                sd_avg_P= Mean_array(log_sd_avg_p, length_log);         // Average SD of average target proportion 
                sd_asd_P= Mean_array(log_sd_asd_p, length_log);         // Average SD of SD of target proportion
                sd_avg_porp= Mean_array(log_sd_avg_porp, length_log);   // Average SD of actual propotion of helpers
                sd_asd_porp= Mean_array(log_sd_asd_porp, length_log);   // Average SD of SD of actual helper proportion
                sd_avg_cord= Mean_array(log_sd_avg_cord, length_log);   // Average SD of actual level of coordination
                sd_asd_cord= Mean_array(log_sd_asd_cord, length_log);   // Average SD of SD of actual coordination level
            // Print
                fprintf(out, "%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", e, n, avg_S, asd_S, avg_P, asd_P, avg_porp, asd_porp, avg_cord, asd_cord, sd_avg_S, sd_asd_S, sd_avg_P, sd_asd_P, sd_avg_porp, sd_asd_porp, sd_avg_cord, sd_asd_cord);
            }
    }}
    free_d_array4d(Popn);
    free_d_array4d(PopnNext);
    free_d_array2d(ProbGrup);
    free_d_array2d(PorpHelp);
    free_d_array2d(PorpCord);
    fclose(out);
    return 0;
}

////////////////////// Functions are defined in below //////////////////////
double Mean_array(double p[], int length){
    int i;
    double temp;
    temp= 0.0;
    if(length> 0){
        for (i=0; i<length; ++i) temp+= p[i]/length;
        return temp;
    }
    else return NAN;
}

double SD_array(double p[], int length){
    if(length> 1){
        int i;
        double avg, temp;
        avg= temp= 0.0;
        for (i=0; i<length; ++i) avg+= p[i]/length;
        for (i=0; i<length; ++i) temp+= pow(p[i]-avg, 2);
        return sqrt(temp/(length-1));
    }
    else return NAN;
}

double sum(double p[], int length){
    int i;
    double sum=0.0;
    for(i=0; i< length; ++i) sum+= p[i];
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
