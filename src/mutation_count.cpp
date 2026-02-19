#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_set>

// Check for OpenMP support
#ifdef _OPENMP
#include <omp.h>
#define PARALLEL_FOR _Pragma("omp parallel for")
#else
#define PARALLEL_FOR
#pragma message("OpenMP not supported. Compilation will proceed without parallel execution.")
#endif


/**
 * Count or locate mutations between germline and input sequences.
 *
 * This function compares each pair of germline and input sequences, counting or locating positions where the bases differ (mutations),
 * optionally skipping the first X positions. Non-mismatch characters (N, ., -) are ignored. Can run in parallel if OpenMP is available.
 *
 * @param germs Vector of germline sequences (character).
 * @param inputs Vector of input sequences (character), same length as germs.
 * @param X Integer. Number of leading positions to skip (default 0).
 * @param parallel Logical. Whether to use parallel processing (default FALSE).
 * @param return_count Logical. If TRUE, return the number of mutations per sequence; if FALSE, return the positions of mutations (1-based).
 *
 * @return Integer vector (if return_count=TRUE) or list of integer vectors (if return_count=FALSE), one per sequence.
 *
 * @details
 * - Sequences are padded with Ns to equal length if needed.
 * - Only positions after X are considered.
 * - Ignores mismatches involving N, ., or -.
 * - If OpenMP is available and parallel=TRUE, runs in parallel.
 *
 * @examples
 * mutation_count(c("ACGT", "ACGT"), c("ACGA", "ACCT"), X=0, parallel=FALSE, return_count=TRUE)
 * mutation_count(c("ACGT", "ACGT"), c("ACGA", "ACCT"), X=0, parallel=FALSE, return_count=FALSE)
 */
// [[Rcpp::export]]
Rcpp::RObject mutation_count(std::vector<std::string> germs, 
                             std::vector<std::string> inputs, 
                             int X = 0, 
                             bool parallel = false, 
                             bool return_count = false) {
  // Set default non-mismatch characters
  std::unordered_set<char> non_mismatch_chars = {'N', '.', '-'};
  
 if (germs.size() != inputs.size()) {
   Rcpp::stop("The size of germs and inputs must be the same.");
 }
 
 size_t num_sequences = germs.size();
 auto pad_with_ns = [](std::string& seq, size_t target_length) {
   if (seq.size() < target_length) {
     seq.append(target_length - seq.size(), 'N');
   }
 };
 
 if (!parallel) {
   if (return_count) {
     std::vector<int> mutation_counts(num_sequences);
     for (size_t i = 0; i < num_sequences; ++i) {
       std::string germ = germs[i];
       std::string input = inputs[i];
       size_t max_length = std::max(germ.size(), input.size());
       pad_with_ns(germ, max_length);
       pad_with_ns(input, max_length);
       int count = 0;
       for (size_t j = 0; j < max_length; ++j) {
         if (j >= static_cast<size_t>(X) && germ[j] != input[j] &&
             non_mismatch_chars.find(input[j]) == non_mismatch_chars.end() &&
             non_mismatch_chars.find(germ[j]) == non_mismatch_chars.end()) {
           count++;
         }
       }
       mutation_counts[i] = count;
     }
     return Rcpp::wrap(mutation_counts);
   } else {
     Rcpp::List snp_list(num_sequences);
     for (size_t i = 0; i < num_sequences; ++i) {
       std::string germ = germs[i];
       std::string input = inputs[i];
       size_t max_length = std::max(germ.size(), input.size());
       pad_with_ns(germ, max_length);
       pad_with_ns(input, max_length);
       std::vector<int> snp_indices;
       for (size_t j = 0; j < max_length; ++j) {
         if (j >= static_cast<size_t>(X) && germ[j] != input[j] &&
             non_mismatch_chars.find(input[j]) == non_mismatch_chars.end() &&
             non_mismatch_chars.find(germ[j]) == non_mismatch_chars.end()) {
           snp_indices.push_back(j + 1);
         }
       }
       snp_list[i] = Rcpp::wrap(snp_indices);
     }
     return snp_list;
   }
 }
 
 if (parallel) {
   if (return_count) {
     std::vector<int> mutation_counts(num_sequences);
     PARALLEL_FOR
     for (size_t i = 0; i < num_sequences; ++i) {
       std::string germ = germs[i];
       std::string input = inputs[i];
       size_t max_length = std::max(germ.size(), input.size());
       pad_with_ns(germ, max_length);
       pad_with_ns(input, max_length);
       int count = 0;
       for (size_t j = 0; j < max_length; ++j) {
         if (j >= static_cast<size_t>(X) && germ[j] != input[j] &&
             non_mismatch_chars.find(input[j]) == non_mismatch_chars.end() &&
             non_mismatch_chars.find(germ[j]) == non_mismatch_chars.end()) {
           count++;
         }
       }
       mutation_counts[i] = count;
     }
     return Rcpp::wrap(mutation_counts);
   } else {
     Rcpp::List snp_list(num_sequences);
     PARALLEL_FOR
     for (size_t i = 0; i < num_sequences; ++i) {
       std::string germ = germs[i];
       std::string input = inputs[i];
       size_t max_length = std::max(germ.size(), input.size());
       pad_with_ns(germ, max_length);
       pad_with_ns(input, max_length);
       std::vector<int> snp_indices;
       for (size_t j = 0; j < max_length; ++j) {
         if (j >= static_cast<size_t>(X) && germ[j] != input[j] &&
             non_mismatch_chars.find(input[j]) == non_mismatch_chars.end() &&
             non_mismatch_chars.find(germ[j]) == non_mismatch_chars.end()) {
           snp_indices.push_back(j + 1);
         }
       }
       snp_list[i] = Rcpp::wrap(snp_indices);
     }
     return snp_list;
   }
 }
 
 return Rcpp::wrap(R_NilValue);
}
