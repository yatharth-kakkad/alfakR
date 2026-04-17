// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <random>
#include <numeric>   // std::accumulate, std::iota
#include <algorithm> // std::sample, std::min/max, std::sort
#include <iterator>  // std::begin, std::end, std::back_inserter
#include <cmath>     // std::round
#include <stdexcept> // For exception catching
// Removed: const int N_CHROMOSOME_TYPES = 22;

// --- Hash function for using std::vector<int> as map key ---
struct VectorHasher {
  std::size_t operator()(const std::vector<int>& v) const {
    std::size_t seed = v.size();
    for(int i : v) {
      seed ^= std::hash<int>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

typedef std::unordered_map<std::vector<int>, long long, VectorHasher> PopulationMap;
typedef std::unordered_map<std::vector<int>, double, VectorHasher> FitnessMap;

// --- Helper Function Definitions ---
// Parse "1.2.3..." string to vector<int>
// expected_n_chr: if > 0, validates length. If <=0, length is determined by string.
std::vector<int> parse_karyotype_string_abm(const std::string& k_str, int expected_n_chr = -1) {
  std::vector<int> cn;
  std::stringstream ss(k_str);
  std::string segment;
  
  if (k_str.empty()) {
    if (expected_n_chr == 0) return cn; // Valid empty karyotype if explicitly expected
    Rcpp::warning("Empty karyotype string provided when non-empty expected.");
    return {}; // Return empty to signify error
  }
  
  while(std::getline(ss, segment, '.')) {
    try {
      int val = std::stoi(segment);
      if (val <= 0) { // Assuming karyotype counts must be positive
        Rcpp::warning("Non-positive count (%d) in karyotype string: %s. Skipping.", val, k_str.c_str());
        return {}; 
      }
      cn.push_back(val);
    } catch (const std::invalid_argument& ia) {
      Rcpp::warning("Non-integer segment '%s' in karyotype string: %s. Skipping.", segment.c_str(), k_str.c_str());
      return {};
    } catch (const std::out_of_range& oor) {
      Rcpp::warning("Integer segment '%s' out of range in karyotype string: %s. Skipping.", segment.c_str(), k_str.c_str());
      return {};
    }
  }
  
  if (expected_n_chr > 0 && static_cast<int>(cn.size()) != expected_n_chr) {
    Rcpp::warning("Karyotype string '%s' has %d types, expected %d. Skipping.", k_str.c_str(), cn.size(), expected_n_chr);
    return {};
  }
  if (cn.empty() && !k_str.empty()){ // e.g. string was "..."
    Rcpp::warning("Karyotype string '%s' parsed to zero elements. Skipping.", k_str.c_str());
    return {};
  }
  return cn;
}

// Convert vector<int> to "1.2.3..." string
std::string karyotype_to_string_abm(const std::vector<int>& cn) { // Renamed to avoid conflict if linked
  if (cn.empty()) return ""; 
  std::stringstream ss;
  for(size_t i = 0; i < cn.size(); ++i) {
    ss << cn[i] << (i == cn.size() - 1 ? "" : ".");
  }
  return ss.str();
}

// ---- GRF fitness -----------------------------------------------------------
double grf_fitness(const std::vector<int>& cn,
                   const std::vector<std::vector<double>>& centroids,
                   double lambda)
{
  if(lambda <= 0.0) return 0.0;                     // safety guard
  double acc = 0.0;
  for(const auto& r : centroids){
    double d2 = 0.0;
    for(size_t k = 0; k < cn.size(); ++k){
      double diff = static_cast<double>(cn[k]) - r[k];
      d2 += diff * diff;
    }
    acc += std::sin(std::sqrt(d2) / lambda);
  }
  
  double fitness_scalar = 1/(3.14159*sqrt(centroids.size()));
  
  return acc*fitness_scalar;                                       // can be <0, leave as is
}


// Get fitness from LUT, return default if not found
double get_fitness_abm( // Renamed
    const std::vector<int>& cn,
    const FitnessMap& fitness_lut,
    double default_fitness_value = 0.0 // Provide a default for novel types
) {
  auto it = fitness_lut.find(cn);
  if (it != fitness_lut.end()) {
    return it->second; 
  } else {
    // Rcpp::warning("Karyotype %s not in LUT, using default fitness %f", 
    //               karyotype_to_string_abm(cn).c_str(), default_fitness_value); // Can be noisy
    return default_fitness_value; 
  }
}

std::pair<std::vector<int>, std::vector<int>> generate_misseg_daughters(
    int n_errors, 
    const std::vector<int>& parent_cn,
    int n_chrom_total, // Sum of elements in parent_cn
    std::mt19937& rng_engine 
) {
  if (parent_cn.empty() || n_errors < 0 || n_chrom_total <= 0 || n_errors > n_chrom_total) { // n_errors can be 0 for faithful
    if (n_errors > 0) { // Only warn if errors were expected but inputs are bad for misseg
      // Rcpp::warning("generate_misseg_daughters: invalid inputs for missegregation.");
    }
    if (n_errors == 0) { // Faithful division
      return {parent_cn, parent_cn};
    }
    return {{},{}}; // Invalid for missegregation attempt
  }
  if (n_errors == 0) { // Faithful division, return two copies of parent
    return {parent_cn, parent_cn};
  }
  
  // Sample which k distinct chromosomes (not pairs) missegregate
  std::vector<int> individual_chrom_indices;
  individual_chrom_indices.reserve(n_chrom_total);
  for(size_t type_idx = 0; type_idx < parent_cn.size(); ++type_idx) {
    for(int k_copy = 0; k_copy < parent_cn[type_idx]; ++k_copy) {
      individual_chrom_indices.push_back(type_idx); // Store the type index
    }
  }
  // Now individual_chrom_indices contains one entry for each individual chromosome, valued by its type index.
  // e.g., for {2,1}, it's {0,0,1}
  
  if(static_cast<int>(individual_chrom_indices.size()) != n_chrom_total) {
    Rcpp::stop("Internal error in generate_misseg_daughters: n_chrom_total doesn't match expanded indices.");
  }
  
  std::vector<int> sampled_positions_for_misseg; // These are positions in the list of all chromosomes
  sampled_positions_for_misseg.resize(n_chrom_total);
  std::iota(sampled_positions_for_misseg.begin(), sampled_positions_for_misseg.end(), 0);
  
  std::vector<int> missegregating_positions; // Store the *positions* of chromosomes that will missegregate
  missegregating_positions.reserve(n_errors);
  std::sample(sampled_positions_for_misseg.begin(), sampled_positions_for_misseg.end(),
              std::back_inserter(missegregating_positions),
              n_errors, rng_engine);
  // No sort needed for positions if we just iterate through them.
  
  std::vector<int> daughter1_cn = parent_cn;
  std::vector<int> daughter2_cn = parent_cn;
  std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
  
  for (int pos : missegregating_positions) {
    int chrom_type_idx = individual_chrom_indices[pos]; // Get the chromosome type that missegregates
    if (uniform_dist(rng_engine) < 0.5) {
      daughter1_cn[chrom_type_idx]++;
      daughter2_cn[chrom_type_idx]--;
    } else {
      daughter1_cn[chrom_type_idx]--;
      daughter2_cn[chrom_type_idx]++;
    }
  }
  
  bool d1_valid = true;
  bool d2_valid = true;
  for(size_t i = 0; i < parent_cn.size(); ++i) { // Loop up to parent_cn.size()
    if(daughter1_cn[i] <= 0) d1_valid = false; // Counts cannot be negative
    if(daughter2_cn[i] <= 0) d2_valid = false; // Changed from <=0 to <0, assuming 0 is viable if some other chrom present
  }
  
  std::pair<std::vector<int>, std::vector<int>> result;
  if (d1_valid) result.first = daughter1_cn;
  if (d2_valid) result.second = daughter2_cn;
  return result;
}


// [[Rcpp::export]]
Rcpp::List run_karyotype_abm(
    Rcpp::List  initial_population_r,
    Rcpp::List  fitness_map_r,                 
    double      p_missegregation,
    double      dt,
    int         n_steps,
    long long   max_population_size,
    double      culling_survival_fraction = 0.1,
    int         record_interval             = 1,
    int         seed                        = -1,
    Rcpp::NumericMatrix grf_centroids       = Rcpp::NumericMatrix(0,0),
    double      grf_lambda                  = NA_REAL
)
{
  if (!R_finite(p_missegregation) || p_missegregation < 0.0 || p_missegregation > 1.0) {
    Rcpp::stop("`p_missegregation` must be finite and in [0, 1].");
  }
  if (!R_finite(dt) || dt <= 0.0) {
    Rcpp::stop("`dt` must be a positive finite value.");
  }
  if (n_steps < 0) {
    Rcpp::stop("`n_steps` must be non-negative.");
  }
  if (!R_finite(culling_survival_fraction) || culling_survival_fraction < 0.0 || culling_survival_fraction > 1.0) {
    Rcpp::stop("`culling_survival_fraction` must be finite and in [0, 1].");
  }
  if (record_interval == 0) {
    Rcpp::stop("`record_interval` must not be zero.");
  }
  if (grf_centroids.nrow() > 0 && (!R_finite(grf_lambda) || grf_lambda <= 0.0)) {
    Rcpp::stop("`grf_lambda` must be a positive finite value when GRF centroids are supplied.");
  }

  PopulationMap population;
  FitnessMap fitness_map; 
  std::mt19937 rng_engine;
  if (seed == -1) {
    std::random_device rd;
    rng_engine.seed(rd());
  } else {
    rng_engine.seed(static_cast<unsigned int>(seed));
  }
  int n_chr_types_sim = -1; // To be determined from first valid karyotype
  
  // ---------------------------------------------------------------------------
  // recording logic
  const bool record_on_interval = (record_interval >= 1);   // specified timesteps
  const bool record_on_cull     = (record_interval <  0);   // only at culling
  const int  interval_every     = record_on_interval ? record_interval : 1;
  
  
  // ---- Decide fitness mode ---------------------------------------------------
  const bool use_grf = (grf_centroids.nrow() > 0) && !Rcpp::NumericVector::is_na(grf_lambda);
  
  std::vector<std::vector<double>> centroids_std;
  if(use_grf){
    centroids_std.reserve(grf_centroids.nrow());
    for(int i = 0; i < grf_centroids.nrow(); ++i){
      Rcpp::NumericVector r_row = grf_centroids.row(i);           
      centroids_std.emplace_back(r_row.begin(), r_row.end());     
    }
  }
  
  // ---------------------------------------------------------------------------
  // Dimension‑sanity: derive or verify n_chr_types_sim in GRF mode
  if(use_grf){
    const int K = static_cast<int>(grf_centroids.ncol());
    if(n_chr_types_sim == -1){
      n_chr_types_sim = K;                         // no dimension yet → take GRF’s K
    }else if(n_chr_types_sim != K){
      Rcpp::stop("Dimension mismatch: initial data use %d chromosome types, "
                   "but GRF centroids have %d columns.", n_chr_types_sim, K);
    }
  }
  
  
  Rcpp::CharacterVector initial_k_names = initial_population_r.names();
  for(int i = 0; i < initial_k_names.size(); ++i) {
    std::string k_str = Rcpp::as<std::string>(initial_k_names[i]);
    double count_r = 0;
    try { count_r = Rcpp::as<double>(initial_population_r[k_str]); }
    catch (...) { /* Handle or skip */ continue; }
    
    std::vector<int> cn;
    if (n_chr_types_sim == -1) {
      cn = parse_karyotype_string_abm(k_str); 
      if (!cn.empty()) {
        n_chr_types_sim = cn.size();
      } else { Rcpp::warning("Failed to parse first karyotype '%s' to determine chromosome count.", k_str.c_str()); continue; }
    } else {
      cn = parse_karyotype_string_abm(k_str, n_chr_types_sim); 
      if (cn.empty()) { continue; }
    }
    
    if(!cn.empty()) {
      if(count_r >= 0) {
        population[cn] = static_cast<long long>(std::round(count_r));
      } else {
        population[cn] = 0LL;
      }
    } 
  }
  
  if (n_chr_types_sim == -1 && initial_k_names.size() > 0) {
    Rcpp::stop("Could not determine a consistent number of chromosome types from initial_population_r.");
  }
  if (initial_k_names.size() == 0) { // If initial_population_r was empty
    // This state could be an error, or an empty simulation.
    // For now, if fitness_map_r is also empty, it might proceed to return an empty result.
    // If fitness_map_r is not empty, n_chr_types_sim needs to be derived from it.
    if (fitness_map_r.size() > 0) {
      Rcpp::CharacterVector fm_names = fitness_map_r.names();
      std::string fm_first_k_str = Rcpp::as<std::string>(fm_names[0]);
      std::vector<int> fm_cn_first = parse_karyotype_string_abm(fm_first_k_str);
      if(!fm_cn_first.empty()) n_chr_types_sim = fm_cn_first.size();
      else Rcpp::stop("Initial population is empty and cannot determine chromosome types from fitness map's first key.");
    } else {
      Rcpp::warning("Initial population and fitness map are both empty. Returning empty results.");
      return Rcpp::List(); // Return empty list
    }
  }
  
  // record initial population:
  Rcpp::CharacterVector fitness_k_names = fitness_map_r.names();
  for(int i = 0; i < fitness_k_names.size(); ++i) {
    std::string k_str = Rcpp::as<std::string>(fitness_k_names[i]);
    double fitness_r = 0;
    try { fitness_r = Rcpp::as<double>(fitness_map_r[k_str]); }
    catch (...) { /* Handle or skip */ continue; }
    
    std::vector<int> cn = parse_karyotype_string_abm(k_str, n_chr_types_sim); // Validate against determined length
    if(!cn.empty()) {
      fitness_map[cn] = fitness_r;
    } 
  }
  

  Rcpp::List results_over_time;
  if ((record_on_interval || record_on_cull) && n_chr_types_sim > 0) {
    Rcpp::NumericVector counts_r_init;
    Rcpp::CharacterVector names_r_init;
    long long initial_total = 0;
    PopulationMap cleaned_initial_population; // To store valid initial types
    
    for(const auto& pair : population) {
      if (use_grf || fitness_map.count(pair.first)) { 
        if (pair.second > 0) { 
          counts_r_init.push_back(static_cast<double>(pair.second));
          names_r_init.push_back(karyotype_to_string_abm(pair.first));
          initial_total += pair.second;
          cleaned_initial_population[pair.first] = pair.second; // Keep this one
        }
      } else if(!use_grf) {     // only warn in LUT mode
        Rcpp::warning("Initial population contains karyotype '%s' not in fitness map. It will be ignored.",
                      karyotype_to_string_abm(pair.first).c_str());
      }
    }
    population = std::move(cleaned_initial_population); // Update population to only valid types
    
    if (names_r_init.size() > 0) {
      counts_r_init.names() = names_r_init;
    }
    results_over_time.push_back(counts_r_init, "0"); 
    if(initial_total == 0 && population.empty()) {
      Rcpp::warning("Initial population is empty after filtering against fitness map or all counts were zero.");
    }
  } else if (record_interval >= 1 && n_chr_types_sim <= 0) {
    Rcpp::warning("Cannot record initial state because number of chromosome types is undetermined.");
  }
  
  std::poisson_distribution<long long> poisson_dist;
  std::binomial_distribution<long long> binomial_dist_ll;
  bool warned_large_dt = false;
  
  for (int step = 1; step <= n_steps; ++step) {
    if (population.empty()) {
      Rcpp::Rcout << "Population extinct at step " << step << std::endl;
      break;
    }
    if (n_chr_types_sim <=0) { // Should have been caught earlier
      Rcpp::stop("Internal error: n_chr_types_sim not properly set before simulation loop.");
    }
    
    PopulationMap net_changes_this_step;
    std::vector<std::vector<int>> current_karyotypes_vec; // Renamed
    current_karyotypes_vec.reserve(population.size());
    for(auto const& [cn_map_key, count_map_val] : population) { // Renamed iteration vars
      if (count_map_val > 0) current_karyotypes_vec.push_back(cn_map_key);
    }
    
    for (const auto& parent_cn : current_karyotypes_vec) {
      auto parent_it = population.find(parent_cn);
      if(parent_it == population.end() || parent_it->second <=0) continue;
      long long parent_count = parent_it->second;
      
      double fitness = use_grf
      ? grf_fitness(parent_cn, centroids_std, grf_lambda)
        : get_fitness_abm(parent_cn, fitness_map);
      
      
      if (fitness <= 0) {
        continue; // Non-positive fitness means no division in the current ABM.
      }
      if (!warned_large_dt && fitness * dt > 0.1) {
        Rcpp::warning("ABM time step may be too large: fitness * dt > 0.1; growth may be underestimated.");
        warned_large_dt = true;
      }
      
      double expected_divisions = static_cast<double>(parent_count) * fitness * dt;
      if (expected_divisions < 0) expected_divisions = 0; // Should not happen if fitness > 0
      
      poisson_dist.param(typename std::poisson_distribution<long long>::param_type(expected_divisions));
      long long n_divs = poisson_dist(rng_engine);
      n_divs = std::min(n_divs, parent_count); // Cannot divide more cells than exist
      
      if (n_divs <= 0) continue;
      
      int n_chrom_total_parent = 0; // Renamed from n_chrom_pairs
      for(int count_val : parent_cn) n_chrom_total_parent += count_val;
      if(n_chrom_total_parent == 0 && n_divs > 0) { // Dividing an empty karyotype?
        net_changes_this_step[parent_cn] -= n_divs; // These divisions effectively lead to loss
        continue;
      }
      
      // Each completed division consumes one parent exactly once. Daughters are the
      // only source of gains, whether the division is faithful or missegregating.
      net_changes_this_step[parent_cn] -= n_divs; // All dividing parents are removed first
      
      for (long long div_idx = 0; div_idx < n_divs; ++div_idx) {
        int k_errors_this_division = 0;
        if (n_chrom_total_parent > 0 && p_missegregation > 0) { // Only attempt binomial if possible errors
          std::binomial_distribution<long long> division_error_count_dist(n_chrom_total_parent, p_missegregation);
          k_errors_this_division = static_cast<int>(division_error_count_dist(rng_engine));
        }
        
        std::pair<std::vector<int>, std::vector<int>> daughters =
          generate_misseg_daughters(k_errors_this_division, parent_cn, n_chrom_total_parent, rng_engine);
        
        if (!daughters.first.empty() && (use_grf || fitness_map.count(daughters.first))) {
          net_changes_this_step[daughters.first]++;
        }
        if (!daughters.second.empty() && (use_grf || fitness_map.count(daughters.second))) {
          net_changes_this_step[daughters.second]++;
        }
      }
    } 
    
    for(const auto& pair_change : net_changes_this_step) { // Renamed pair
      if (use_grf || fitness_map.count(pair_change.first)) {
        population[pair_change.first] += pair_change.second;
      }
    }
    
    long long current_total_pop = 0; // Renamed from next_total_pop
    for (auto it = population.begin(); it != population.end(); ) {
      if (it->second <= 0 || (!use_grf && !fitness_map.count(it->first))) {
        it = population.erase(it);
      } else {
        current_total_pop += it->second;
        ++it;
      }
    }
    
    if (max_population_size > 0 && current_total_pop > max_population_size) {
      
      if (record_on_cull) {
        Rcpp::NumericVector counts_cull;
        Rcpp::CharacterVector names_cull;
        for (const auto& pair_rec : population) {
          counts_cull.push_back(static_cast<double>(pair_rec.second));
          names_cull.push_back(karyotype_to_string_abm(pair_rec.first));
        }
        if (counts_cull.length() > 0) counts_cull.names() = names_cull;
        results_over_time.push_back(counts_cull, std::to_string(step));
      }
      
      double sampling_fraction = culling_survival_fraction; 
      // Rcpp::Rcout << "Step " << step << ": Population " << current_total_pop // Optional verbose logging
      //             << " exceeded cap " << max_population_size
      //             << ". Culling to approx. " << static_cast<long long>(round(current_total_pop * sampling_fraction)) << " cells." << std::endl;
      
      PopulationMap sampled_population;
      if (sampling_fraction > 0.0 && sampling_fraction <= 1.0) {
        for (auto const& [cn_sample, count_sample] : population) { // Renamed iteration vars
          binomial_dist_ll.param(typename std::binomial_distribution<long long>::param_type(count_sample, sampling_fraction));
          long long sampled_count = binomial_dist_ll(rng_engine);
          if (sampled_count > 0) {
            sampled_population[cn_sample] = sampled_count;
          }
        }
        population = std::move(sampled_population); 
      } else { 
        // Rcpp::Rcout << "Step " << step << ": Sampling fraction is zero or invalid. Population culled entirely." << std::endl; // Optional
        population.clear(); 
      }
      // Recalculate current_total_pop after culling for accurate reporting if needed immediately
      current_total_pop = 0;
      for(const auto& pair_recalc : population) current_total_pop += pair_recalc.second;
    } 
    
    // const int report_freq = std::max(1, n_steps / 10); // Reporting logic can be kept if desired
    // if (step % report_freq == 0 || step == n_steps) { ... }
    
    if (record_on_interval && (step % interval_every == 0 || step == n_steps)) {
      Rcpp::NumericVector counts_r_step; // Renamed
      Rcpp::CharacterVector names_r_step; // Renamed
      if (n_chr_types_sim > 0) { // Only proceed if n_chr_types_sim is valid
        for(const auto& pair_record : population) { // Renamed
          counts_r_step.push_back(static_cast<double>(pair_record.second));
          names_r_step.push_back(karyotype_to_string_abm(pair_record.first));
        }
        if(counts_r_step.length() > 0) counts_r_step.names() = names_r_step;
      }
      results_over_time.push_back(counts_r_step, std::to_string(step));
    }
    Rcpp::checkUserInterrupt();
  } 
  return results_over_time;
}
