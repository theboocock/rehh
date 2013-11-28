#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gperftools/profiler.h>
#include <R.h>
#include "ehh_utils.h"

	

void r_scan_hh(int *Rdata,
			   int *number_SNPs,
			   int *number_chromosomes,
			   double *map,
			   double *IHH,
			   double *IES,
			   int *min_number_haplotypes,
			   double *min_EHH,
			   double *min_EHHS,
				 int    *big_gap,
			   int    *small_gap,
				 int *small_gap_penalty
)

{
	int i,j;
	int focal_SNP;
	int allele;
	int n;
	int **data;
	int **number_haplotypes;
	double ehh,ehhs;
	double *tmp;
	double **EHH;
	double *EHHS;
	/* Change map positions based on the distance */
    ProfilerStart("/home/smilefreak/profiles/new_profile.prof");
	double *gap_list;
	int map_distance = 0;
	int kilobase = 1000;
	*small_gap_penalty = *small_gap_penalty * kilobase;
	*small_gap = *small_gap * kilobase;	
	*big_gap = *big_gap * kilobase;
	gap_list = (double *) malloc((*number_SNPs - 1)*sizeof(double));
	for(i = 0;i < *number_SNPs - 1;i++){
		int gap = map[i+1] - map[i];
        //Rprintf("gap %d \n",gap);
		if( *big_gap != 0 ){
			if(gap >= *big_gap ){
				gap_list[i] = 0;
			}else{
                gap_list[i] =1.0; 
            }
		
		}else if(*small_gap != 0){
			if(gap >= *small_gap){
				gap_list[i] = *small_gap_penalty/(double)gap;
			}else{
			    gap_list[i] = 1.0;
            }
		}else{
			    gap_list[i] = 1.0;
		}
	}
	data = (int **) malloc(*number_chromosomes * sizeof(int *));
	for (i = 0; i < *number_chromosomes; i++) {
		data[i] = (int *) malloc(*number_SNPs * sizeof(int));
	}
	for (i = 0; i < *number_chromosomes; i++) {
		for (j = 0; j < *number_SNPs; j++) {	
			if (Rdata[(j * *number_chromosomes) + i] == 0) data[i][j] = 9;
			if (Rdata[(j * *number_chromosomes) + i] == 1) data[i][j] = 0;
			if (Rdata[(j * *number_chromosomes) + i] == 2) data[i][j] = 1;
		}
	}
	EHHS = (double *) malloc(*number_SNPs * sizeof(double));
	number_haplotypes = (int **) malloc(*number_SNPs * sizeof(int *));
	EHH = (double **) malloc(*number_SNPs * sizeof(double *));
	for (j = 0; j < *number_SNPs; j++) {	
		number_haplotypes[j] = (int *) malloc(2 * sizeof(int));
		EHH[j] = (double *) malloc(2 * sizeof(double));
	}
	tmp = (double *) malloc(*number_SNPs * sizeof(double));
	for (focal_SNP = 0; focal_SNP < *number_SNPs; focal_SNP++) {
		for (j = 0; j < *number_SNPs; j++) {	
			for (allele = 0; allele < 2; allele++) {
				EHH[j][allele] = 0.0;
				number_haplotypes[j][allele] = 0;
			}
		}
		EHH[focal_SNP][0] = EHH[focal_SNP][1] = 1.0;
//		number_haplotypes[focal_SNP][0] = number_haplotypes[focal_SNP][1] = 0;
		for (i = 0; i < *number_chromosomes; i++) {
			if (data[i][focal_SNP] == 0) number_haplotypes[focal_SNP][0] += 1;
			if (data[i][focal_SNP] == 1) number_haplotypes[focal_SNP][1] += 1;
		}
		for (allele = 0; allele < 2; allele++) {
			if (number_haplotypes[focal_SNP][allele] < *min_number_haplotypes) {
				EHH[focal_SNP][allele] = 0.0;
			} else {
				for (j = (focal_SNP - 1); j >= 0; j--) {
					compute_EHH(data,*number_chromosomes,*number_SNPs,focal_SNP,j,LEFT,allele,&ehh,&n);
					if (n > *min_number_haplotypes) {
						EHH[j][allele] = ehh;
						number_haplotypes[j][allele] = n;
							if(ehh < *min_EHH){ 
									break;
									}
					} else {
						break;
					}

				}
				for (j = (focal_SNP + 1); j < *number_SNPs; j++) {
					compute_EHH(data,*number_chromosomes,*number_SNPs,focal_SNP,j,RIGHT,allele,&ehh,&n);
					if (n > *min_number_haplotypes) {
						EHH[j][allele] = ehh;
						number_haplotypes[j][allele] = n;
							if(ehh < *min_EHH){ 
									break;
									}
					} else {
						break;
					}

				}
			}
		}
		for (allele = 0; allele < 2; allele++) {
			for (j = 0; j < *number_SNPs; j++) {
				tmp[j] = EHH[j][allele];
			}			
			IHH[(allele * *number_SNPs) + focal_SNP] = integrate_gaps(map,tmp,*number_SNPs,*min_EHH,gap_list);
		}	
		for (j = 0; j < *number_SNPs; j++) {
			EHHS[j] = 0.0;
			number_haplotypes[j][0] = 0;
		}
		EHHS[focal_SNP] = 1.0;
		number_haplotypes[focal_SNP][0] = 0;
		for (i = 0; i < *number_chromosomes; i++) {
			if (data[i][focal_SNP] == 0) number_haplotypes[focal_SNP][0] += 1;
			if (data[i][focal_SNP] == 1) number_haplotypes[focal_SNP][0] += 1;
		}
		if (number_haplotypes[focal_SNP][0] < *min_number_haplotypes) {
			EHHS[focal_SNP] = 0.0;
		} else {
			for (j = (focal_SNP - 1); j >= 0; j--) {
				compute_EHHS(data,*number_chromosomes,*number_SNPs,focal_SNP,j,LEFT,&ehhs,&n);
					if (n > *min_number_haplotypes) {
						EHHS[j] = ehhs;
						number_haplotypes[j][0] = n;
							if(ehhs < *min_EHHS){ 
									break;
									}
					} else {
						break;
					}


			}
			for (j = (focal_SNP + 1); j < *number_SNPs; j++) {
				compute_EHHS(data,*number_chromosomes,*number_SNPs,focal_SNP,j,RIGHT,&ehhs,&n);
					if (n > *min_number_haplotypes) { 
						EHHS[j] = ehhs;
						number_haplotypes[j][0] = n;
							if(ehhs < *min_EHHS){ 
									break;
									}
					} else {
						break;
					}

			}
		}
		IES[focal_SNP] = integrate_gaps(map,EHHS,*number_SNPs,*min_EHHS,gap_list);
    //Rprintf("SNP = %d \tout of\t %d SNPs\n",focal_SNP+1,*number_SNPs);
	}
	for (i = 0; i < *number_chromosomes; i++) {
		free(data[i]);
	}
	free(data);	
	free(gap_list);
	free(EHHS);
	for (j = 0; j < *number_SNPs; j++) {
		free(number_haplotypes[j]);
		free(EHH[j]);
	}
	free(number_haplotypes);
	free(EHH);
	free(tmp);
    ProfilerStop();
}
