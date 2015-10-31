#include <Rcpp.h>
using namespace Rcpp;

// Declare functions
std::vector<float> calcFijPairwiseCpp(NumericMatrix alfreq1,
                                      NumericMatrix alfreq2,
                                      NumericMatrix ref_al_freq,
                                      int n_loci,
                                      int n_alleles_per_loci,
                                      int n_gene_copies);

NumericMatrix calcAlleleFreqCpp(NumericMatrix gen_data,
                                int n_loci,
                                int n_alleles_per_loci);


std::vector<float> calcAlleleFreqCppInd(int alleles_1,
                              int alleles_2,
                              int n_alleles_per_loci);

// ***************************************************
// ***************************************************
// *** CALCULATE PAIRWISE FIJ BETWEEN A PAIR OF INDIVIDUALS
// ***************************************************
// ***************************************************

// Returns a vector with Fij estimate at each locus

// [[Rcpp::export]]
std::vector<float> calcFijPairwiseCpp(NumericMatrix alfreq1,
                                      NumericMatrix alfreq2,
                                      NumericMatrix ref_al_freq,
                                      int n_loci,
                                      int n_alleles_per_loci,
                                      int n_gene_copies){
  
  float denom[n_loci];
  float numer[n_loci];
  std::vector<float> fij(n_loci);

  
  for(int locus = 0;  locus < n_loci; ++locus){ // Loop through loci
    
    denom[locus] = 0; // Initialize numerator and denominator
    numer[locus] = 0;
    
      for(int allele = 0; allele < n_alleles_per_loci; ++allele){ // Loop through alleles
     
        // Calculate Numerator and Denominator of Loiselle et al. 1995
        numer[locus]  +=  (alfreq1(locus, allele) - ref_al_freq(locus, allele)) *
                          (alfreq2(locus, allele) - ref_al_freq(locus, allele)) +
          (ref_al_freq(locus, allele)*(1 - ref_al_freq(locus, allele))) / (n_gene_copies- 1);
        
        denom[locus] += ref_al_freq(locus, allele) * (1 - ref_al_freq(locus, allele));
        
      } // End allele loop
    
    fij[locus] = numer[locus] / denom[locus]; // Calculate Fij per locus.. should equal -999 for missing data
    
  } // End loci loop
  
  return(fij);
}



// [[Rcpp::export]]
std::vector<float> calcFijPopCpp(NumericMatrix gen_data,
                                  NumericMatrix ref_al_freq,
                                  NumericVector offspr_gen,
                                  int n_loci,
                                  int n_alleles_per_loci,
                                  int n_gene_copies){
  
  int n_ad = gen_data.nrow(); // Save number of adults
  
  // Initalize a 3d array that will save allele frequency for each individual
  float ind_al_freq[n_loci][n_alleles_per_loci][n_ad+1]; // Add one - first index will be offspring
  int row = 0; // Needed for loop through 3d array
  std::vector<float> fij(n_ad); // Length equal to number of adults - will be returned by function
  std::vector<float> fij_by_locus(n_loci);
  NumericMatrix  alfreq1(n_loci, n_alleles_per_loci);
  NumericMatrix  alfreq2(n_loci, n_alleles_per_loci);
  
  std::vector<float> freq_table(n_alleles_per_loci);

  // Loop through individuals and calculate allele frequency
  // First index of in_al_freq will be offspring allele frequency
  for(int indi = 0; indi < (n_ad + 1); ++indi){
    
    row = 0; // Row of individual allele frequency 3d array - corresponds to locus
    
    for(int locus = 0; locus < (n_loci * 2); locus += 2){ // Loop through loci
      
      // If on first individual, save it as offspring allele frequency
      if(indi == 0){
        freq_table = calcAlleleFreqCppInd(offspr_gen[locus],
                                                             offspr_gen[locus + 1],
                                                             n_alleles_per_loci);
      } else { 
       freq_table = calcAlleleFreqCppInd(gen_data(indi-1, locus),
                                                    gen_data(indi-1, locus + 1),
                                                    n_alleles_per_loci);
      }
     
      // Code to print out individual frequency tables
      /*
      std::cout << "Indi: " << indi <<" | Locus - "<< locus << "|:";
      for (std::vector<float>::iterator i = freq_table.begin(); i != freq_table.end(); ++i)
          std::cout << *i << ' ';
        Rcpp::Rcout << "\n";
       */
      
      // Have to then loop through allele frequency table to fill in individual part 
      // of the array
      int allele = 0;
      for (std::vector<float>::iterator i = freq_table.begin(); i != freq_table.end(); ++i){
        ind_al_freq[row][allele][indi] = *i;
        allele += 1;
      } // End allele loop
      
      row += 1; // Jump to next row in array, which is the next locus
      
    } // End loci loop
  } // End individuals loop
  
  
  
  // Calculate pairwise Fij between all pairs of individuals
  
for(int ind = 1; ind < (n_ad + 1); ++ind){ 
    
    // Loop through and save allele frequency for each individual quickly
    for(int locus = 0; locus < n_loci; ++locus){
      for(int allele = 0; allele < n_alleles_per_loci; ++allele){
        
        if(ind == 1){  
          alfreq1(locus, allele) = ind_al_freq[locus][allele][0]; // Save offspring al freq on first loop
        //Rcout << ind_al_freq[locus][allele][id1];
        } 
        alfreq2(locus, allele) = ind_al_freq[locus][allele][ind];
        //Rcout << ind_al_freq[locus][allele][id2];
      }
     }
    
    // Returns a vector with fij estimate by locus  
    fij_by_locus =   calcFijPairwiseCpp(alfreq1,
                                        alfreq2,
                                        ref_al_freq,
                                        n_loci, 
                                        n_alleles_per_loci,
                                        n_gene_copies);
  
  
  // Take average of fij by locus and save as mean estimate for each adult
  for(int locus = 0; locus < n_loci; ++locus){
    fij[ind-1] += fij_by_locus[locus];
  }
  
  fij[ind-1] = fij[ind-1] / n_loci;
  
  Rcout << "At end " << ind << "\n";
  
 } // End Adult loop
  
return(fij);

}



// ********************************************

// Calcs allele frequency for each locus and returns a n_loci * n_alleles matrix
// [[Rcpp::export]]
NumericMatrix calcAlleleFreqCpp(NumericMatrix gen_data,
                                int n_loci,
                                int n_alleles_per_loci){

  // Initialize matrix that will hold allele frequencies
  NumericMatrix freq_table(n_loci, n_alleles_per_loci);
  
  
  for(int row = 0; row < gen_data.nrow(); ++row){ // Loop through individuals
    
    int col_index = 0;
    for(int locus = 0; locus < n_loci; ++locus){ // Loop +2 over loci for column index
      
       freq_table(locus, (gen_data(row, col_index) - 1))     += 1; // First column of locus
       freq_table(locus, (gen_data(row, col_index + 1) - 1)) += 1; // Second column of locus
       col_index += 2;
       
    } // End locus loop
  } // End individual loop
    
  // Divide by total number of alleles to get frequency
  // Loop over rows and columns of freq table
  
  int n_tot_alleles = gen_data.nrow() * 2; // Number of total alleles is n_ind * 2
  
  for(int i = 0; i < freq_table.nrow(); i++){
    for(int j = 0; j < freq_table.ncol(); j++){
      freq_table(i, j) = freq_table(i, j) / n_tot_alleles;
    }
  }

  return(freq_table);
  
}




// Calculate allele frequency of individuals
// Needs integer inputs
std::vector<float> calcAlleleFreqCppInd(int alleles_1,
                              int alleles_2,
                              int n_alleles_per_loci){
  
  
  std::vector<float> freq_table(n_alleles_per_loci);
  float AlleleTotal = 2; // FOr looking at individual alleles

  
  // Loop through number of alleles and add them to frequency table
  for(int allele = 0; allele < 1; ++allele){

    freq_table[alleles_1 - 1] += 1; // Minus 1 is because of zero indexing in Cpp
    freq_table[alleles_2 - 1] += 1;
  }
  
  // Divide by total number of alleles to get frequency
  for(int allele = 0; allele < n_alleles_per_loci; ++allele){
    freq_table[allele] = freq_table[allele] / AlleleTotal;
  }

  return(freq_table);
  
}



