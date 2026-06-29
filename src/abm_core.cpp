// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <numeric>   // std::accumulate, std::iota
#include <algorithm> // std::sample, std::min/max, std::sort
#include <iterator>  // std::begin, std::end, std::back_inserter
#include <cmath>     // std::round
#include <limits>
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
typedef std::unordered_map<std::vector<int>, std::vector<std::vector<int>>, VectorHasher> AdjacencyMap;
typedef std::unordered_set<std::vector<int>, VectorHasher> KaryotypeSet;

namespace {

constexpr double ABM_MAX_EXACT_INTEGER = 9007199254740991.0;

bool is_integer_valued_double(double x) {
  return std::floor(x) == x;
}

long long validate_initial_population_count(double count_r, const std::string& k_str) {
  if (!R_finite(count_r)) {
    Rcpp::stop("`initial_population_r[['%s']]` must be finite.", k_str.c_str());
  }
  if (count_r < 0.0) {
    Rcpp::stop("`initial_population_r[['%s']]` must be non-negative.", k_str.c_str());
  }
  if (!is_integer_valued_double(count_r)) {
    Rcpp::stop("`initial_population_r[['%s']]` must be integer-valued.", k_str.c_str());
  }
  if (count_r > ABM_MAX_EXACT_INTEGER) {
    Rcpp::stop("`initial_population_r[['%s']]` exceeds the supported exact-integer range.", k_str.c_str());
  }
  if (count_r > static_cast<double>(std::numeric_limits<long long>::max())) {
    Rcpp::stop("`initial_population_r[['%s']]` exceeds the supported long long range.", k_str.c_str());
  }
  return static_cast<long long>(count_r);
}

} // namespace

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
      if (val < 0) {
        Rcpp::warning("Negative count (%d) in karyotype string: %s. Skipping.", val, k_str.c_str());
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
    if(daughter1_cn[i] < 0) d1_valid = false;
    if(daughter2_cn[i] < 0) d2_valid = false;
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
    try {
      count_r = Rcpp::as<double>(initial_population_r[k_str]);
    } catch (...) {
      Rcpp::stop("`initial_population_r[['%s']]` must be coercible to a finite numeric count.", k_str.c_str());
    }
    
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
      population[cn] = validate_initial_population_count(count_r, k_str);
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

namespace {

double clamp_probability(double x) {
  if (!R_finite(x)) {
    Rcpp::stop("Probability calculation produced a non-finite value.");
  }
  return std::max(0.0, std::min(1.0, x));
}

PopulationMap parse_population_list_abm(Rcpp::List population_r, int& n_chr_types) {
  PopulationMap out;
  Rcpp::CharacterVector names = population_r.names();
  if (names.size() == 0) {
    Rcpp::stop("`initial_population_r` must be a named non-empty list.");
  }
  for (int i = 0; i < names.size(); ++i) {
    std::string k_str = Rcpp::as<std::string>(names[i]);
    std::vector<int> cn = n_chr_types < 0
      ? parse_karyotype_string_abm(k_str)
      : parse_karyotype_string_abm(k_str, n_chr_types);
    if (cn.empty()) {
      Rcpp::stop("Invalid karyotype in `initial_population_r`: '%s'.", k_str.c_str());
    }
    if (n_chr_types < 0) {
      n_chr_types = static_cast<int>(cn.size());
    }
    double count_r = Rcpp::as<double>(population_r[k_str]);
    long long count = validate_initial_population_count(count_r, k_str);
    if (count > 0) {
      out[cn] += count;
    }
  }
  if (out.empty()) {
    Rcpp::stop("`initial_population_r` must contain at least one positive count.");
  }
  return out;
}

FitnessMap parse_fitness_list_abm(Rcpp::List fitness_r, int expected_n_chr, const char* arg_name) {
  FitnessMap out;
  Rcpp::CharacterVector names = fitness_r.names();
  if (names.size() == 0) {
    Rcpp::stop("`%s` must be a named non-empty list.", arg_name);
  }
  for (int i = 0; i < names.size(); ++i) {
    std::string k_str = Rcpp::as<std::string>(names[i]);
    std::vector<int> cn = parse_karyotype_string_abm(k_str, expected_n_chr);
    if (cn.empty()) {
      Rcpp::stop("Invalid karyotype in `%s`: '%s'.", arg_name, k_str.c_str());
    }
    double fitness = Rcpp::as<double>(fitness_r[k_str]);
    if (!R_finite(fitness)) {
      Rcpp::stop("`%s[['%s']]` must be finite.", arg_name, k_str.c_str());
    }
    out[cn] = fitness;
  }
  return out;
}

AdjacencyMap parse_adjacency_list_abm(Rcpp::List adjacency_r, int expected_n_chr) {
  AdjacencyMap out;
  Rcpp::CharacterVector names = adjacency_r.names();
  if (names.size() == 0) {
    Rcpp::stop("`adjacency_r` must be a named non-empty list.");
  }
  for (int i = 0; i < names.size(); ++i) {
    std::string source_str = Rcpp::as<std::string>(names[i]);
    std::vector<int> source = parse_karyotype_string_abm(source_str, expected_n_chr);
    if (source.empty()) {
      Rcpp::stop("Invalid source karyotype in `adjacency_r`: '%s'.", source_str.c_str());
    }
    Rcpp::CharacterVector neighbor_names = Rcpp::as<Rcpp::CharacterVector>(adjacency_r[source_str]);
    std::vector<std::vector<int>> neighbors;
    neighbors.reserve(neighbor_names.size());
    for (int j = 0; j < neighbor_names.size(); ++j) {
      std::string neighbor_str = Rcpp::as<std::string>(neighbor_names[j]);
      std::vector<int> neighbor = parse_karyotype_string_abm(neighbor_str, expected_n_chr);
      if (neighbor.empty()) {
        Rcpp::stop("Invalid neighbor karyotype '%s' in `adjacency_r[['%s']]`.",
                   neighbor_str.c_str(), source_str.c_str());
      }
      neighbors.push_back(neighbor);
    }
    out[source] = std::move(neighbors);
  }
  return out;
}

KaryotypeSet parse_karyotype_set_abm(Rcpp::CharacterVector karyotypes_r, int expected_n_chr, const char* arg_name) {
  KaryotypeSet out;
  for (int i = 0; i < karyotypes_r.size(); ++i) {
    std::string k_str = Rcpp::as<std::string>(karyotypes_r[i]);
    std::vector<int> cn = parse_karyotype_string_abm(k_str, expected_n_chr);
    if (cn.empty()) {
      Rcpp::stop("Invalid karyotype in `%s`: '%s'.", arg_name, k_str.c_str());
    }
    out.insert(cn);
  }
  if (out.empty()) {
    Rcpp::stop("`%s` must contain at least one karyotype.", arg_name);
  }
  return out;
}

Rcpp::NumericVector population_to_named_vector_abm(const PopulationMap& population) {
  Rcpp::NumericVector counts;
  Rcpp::CharacterVector names;
  for (const auto& pair : population) {
    if (pair.second > 0) {
      counts.push_back(static_cast<double>(pair.second));
      names.push_back(karyotype_to_string_abm(pair.first));
    }
  }
  if (counts.length() > 0) {
    counts.names() = names;
  }
  return counts;
}

long long population_total_abm(const PopulationMap& population) {
  long long total = 0;
  for (const auto& pair : population) {
    if (pair.second > 0) {
      total += pair.second;
    }
  }
  return total;
}

long long population_diversity_abm(const PopulationMap& population) {
  long long diversity = 0;
  for (const auto& pair : population) {
    if (pair.second > 0) {
      ++diversity;
    }
  }
  return diversity;
}

long long transition_total_abm(const PopulationMap& population, const KaryotypeSet& transition_set) {
  long long total = 0;
  for (const auto& pair : population) {
    if (pair.second > 0 && transition_set.count(pair.first) > 0) {
      total += pair.second;
    }
  }
  return total;
}

double transition_fraction_abm(const PopulationMap& population, const KaryotypeSet& transition_set) {
  long long total = population_total_abm(population);
  if (total <= 0) {
    return 0.0;
  }
  return static_cast<double>(transition_total_abm(population, transition_set)) /
    static_cast<double>(total);
}

void add_metric_row_abm(std::vector<std::string>& conditions,
                        std::vector<int>& steps,
                        std::vector<double>& times,
                        std::vector<double>& totals,
                        std::vector<double>& diversities,
                        std::vector<double>& transition_totals,
                        std::vector<double>& transition_fractions,
                        const std::string& condition,
                        int step,
                        double dt,
                        const PopulationMap& population,
                        const KaryotypeSet& transition_set) {
  conditions.push_back(condition);
  steps.push_back(step);
  times.push_back(static_cast<double>(step) * dt);
  totals.push_back(static_cast<double>(population_total_abm(population)));
  diversities.push_back(static_cast<double>(population_diversity_abm(population)));
  transition_totals.push_back(static_cast<double>(transition_total_abm(population, transition_set)));
  transition_fractions.push_back(transition_fraction_abm(population, transition_set));
}

PopulationMap transition_treatment_step_abm(
    const PopulationMap& population,
    const FitnessMap& untreated_fitness,
    const FitnessMap& treated_fitness,
    const AdjacencyMap& adjacency,
    const KaryotypeSet& transition_set,
    int step,
    int tau1_step,
    int tau2_step,
    bool use_second_treatment,
    double p_missegregation,
    double base_death_rate,
    double base_birth_rate,
    double fitness_birth_scale,
    double second_treatment_strength,
    std::mt19937& rng_engine) {
  PopulationMap after_death;
  PopulationMap after_move;
  PopulationMap next_population;
  std::binomial_distribution<long long> binomial_dist;

  for (const auto& pair : population) {
    const std::vector<int>& k = pair.first;
    long long count = pair.second;
    if (count <= 0) {
      continue;
    }
    bool second_active = use_second_treatment && tau2_step >= 0 &&
      step >= tau2_step && transition_set.count(k) > 0;
    double death_probability = base_death_rate +
      (second_active ? second_treatment_strength : 0.0);
    death_probability = clamp_probability(death_probability);
    binomial_dist.param(typename std::binomial_distribution<long long>::param_type(
      count, 1.0 - death_probability));
    long long survivors = binomial_dist(rng_engine);
    if (survivors > 0) {
      after_death[k] += survivors;
    }
  }

  for (const auto& pair : after_death) {
    const std::vector<int>& k = pair.first;
    long long count = pair.second;
    if (count <= 0) {
      continue;
    }
    binomial_dist.param(typename std::binomial_distribution<long long>::param_type(
      count, p_missegregation));
    long long movers = binomial_dist(rng_engine);
    long long stayers = count - movers;
    if (stayers > 0) {
      after_move[k] += stayers;
    }

    auto adj_it = adjacency.find(k);
    if (movers <= 0 || adj_it == adjacency.end() || adj_it->second.empty()) {
      if (movers > 0) {
        after_move[k] += movers;
      }
      continue;
    }

    const std::vector<std::vector<int>>& neighbors = adj_it->second;
    std::uniform_int_distribution<std::size_t> neighbor_dist(0, neighbors.size() - 1);
    for (long long i = 0; i < movers; ++i) {
      after_move[neighbors[neighbor_dist(rng_engine)]] += 1;
    }
  }

  const FitnessMap& active_fitness = step < tau1_step ? untreated_fitness : treated_fitness;
  for (const auto& pair : after_move) {
    const std::vector<int>& k = pair.first;
    long long count = pair.second;
    if (count <= 0) {
      continue;
    }
    auto fit_it = active_fitness.find(k);
    double fitness = fit_it == active_fitness.end() ? 0.0 : fit_it->second;
    double birth_probability = clamp_probability(base_birth_rate + fitness_birth_scale * fitness);
    binomial_dist.param(typename std::binomial_distribution<long long>::param_type(
      count, birth_probability));
    long long births = binomial_dist(rng_engine);
    next_population[k] += count + births;
  }

  return next_population;
}

// When `track_tau2` is set, the transition fraction of this same run's
// trajectory is monitored at every step (not just recorded steps) and the
// step of its maximum is written to `detected_tau2_step`/
// `detected_best_transition_fraction`. This lets condition 1's own realized
// trajectory determine when the second treatment should trigger, rather than
// inferring it from an unrelated random draw.
Rcpp::List simulate_transition_condition_abm(
    const PopulationMap& initial_population,
    const FitnessMap& untreated_fitness,
    const FitnessMap& treated_fitness,
    const AdjacencyMap& adjacency,
    const KaryotypeSet& transition_set,
    const std::string& condition,
    int tau1_step,
    int tau2_step,
    bool use_second_treatment,
    double p_missegregation,
    double base_death_rate,
    double base_birth_rate,
    double fitness_birth_scale,
    double second_treatment_strength,
    double dt,
    int n_steps,
    int record_interval,
    std::mt19937& rng_engine,
    std::vector<std::string>& metric_conditions,
    std::vector<int>& metric_steps,
    std::vector<double>& metric_times,
    std::vector<double>& metric_totals,
    std::vector<double>& metric_diversities,
    std::vector<double>& metric_transition_totals,
    std::vector<double>& metric_transition_fractions,
    bool track_tau2 = false,
    int* detected_tau2_step = nullptr,
    double* detected_best_transition_fraction = nullptr) {
  PopulationMap population = initial_population;
  Rcpp::List records;
  bool record_on_interval = record_interval >= 1;
  int interval_every = record_on_interval ? record_interval : 1;
  int best_tau2_step = tau1_step;
  double best_transition_fraction = -1.0;

  records.push_back(population_to_named_vector_abm(population), "0");
  add_metric_row_abm(metric_conditions, metric_steps, metric_times, metric_totals,
                     metric_diversities, metric_transition_totals,
                     metric_transition_fractions, condition, 0, dt,
                     population, transition_set);

  for (int step = 1; step <= n_steps; ++step) {
    if (!population.empty()) {
      population = transition_treatment_step_abm(
        population, untreated_fitness, treated_fitness, adjacency, transition_set,
        step, tau1_step, tau2_step, use_second_treatment, p_missegregation,
        base_death_rate, base_birth_rate, fitness_birth_scale,
        second_treatment_strength, rng_engine);
    }
    if (track_tau2 && step >= tau1_step) {
      double fraction = transition_fraction_abm(population, transition_set);
      if (fraction > best_transition_fraction) {
        best_transition_fraction = fraction;
        best_tau2_step = step;
      }
    }
    if (record_on_interval && (step % interval_every == 0 || step == n_steps)) {
      records.push_back(population_to_named_vector_abm(population), std::to_string(step));
      add_metric_row_abm(metric_conditions, metric_steps, metric_times, metric_totals,
                         metric_diversities, metric_transition_totals,
                         metric_transition_fractions, condition, step, dt,
                         population, transition_set);
    }
    Rcpp::checkUserInterrupt();
  }
  if (track_tau2) {
    if (detected_tau2_step != nullptr) *detected_tau2_step = best_tau2_step;
    if (detected_best_transition_fraction != nullptr) *detected_best_transition_fraction = best_transition_fraction;
  }
  return records;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List run_transition_treatment_abm(
    Rcpp::List initial_population_r,
    Rcpp::List untreated_fitness_map_r,
    Rcpp::List treated_fitness_map_r,
    Rcpp::List adjacency_r,
    Rcpp::CharacterVector transition_karyotypes_r,
    double p_missegregation,
    double base_death_rate,
    double base_birth_rate,
    double fitness_birth_scale,
    double second_treatment_strength,
    int tau1_step,
    double dt,
    int n_steps,
    int record_interval = 1,
    int seed = -1) {
  if (!R_finite(p_missegregation) || p_missegregation < 0.0 || p_missegregation > 1.0) {
    Rcpp::stop("`p_missegregation` must be finite and in [0, 1].");
  }
  if (!R_finite(base_death_rate) || base_death_rate < 0.0 || base_death_rate > 1.0) {
    Rcpp::stop("`base_death_rate` must be finite and in [0, 1].");
  }
  if (!R_finite(base_birth_rate) || base_birth_rate < 0.0 || base_birth_rate > 1.0) {
    Rcpp::stop("`base_birth_rate` must be finite and in [0, 1].");
  }
  if (!R_finite(fitness_birth_scale)) {
    Rcpp::stop("`fitness_birth_scale` must be finite.");
  }
  if (!R_finite(second_treatment_strength) || second_treatment_strength < 0.0 || second_treatment_strength > 1.0) {
    Rcpp::stop("`second_treatment_strength` must be finite and in [0, 1].");
  }
  if (tau1_step < 0) {
    Rcpp::stop("`tau1_step` must be non-negative.");
  }
  if (!R_finite(dt) || dt <= 0.0) {
    Rcpp::stop("`dt` must be a positive finite value.");
  }
  if (n_steps < 0) {
    Rcpp::stop("`n_steps` must be non-negative.");
  }
  if (record_interval <= 0) {
    Rcpp::stop("`record_interval` must be positive.");
  }

  std::mt19937 condition1_rng;
  std::mt19937 condition2_rng;
  if (seed == -1) {
    std::random_device rd;
    unsigned int base_seed = rd();
    condition1_rng.seed(base_seed);
    condition2_rng.seed(base_seed + 1U);
  } else {
    condition1_rng.seed(static_cast<unsigned int>(seed));
    condition2_rng.seed(static_cast<unsigned int>(seed) + 1U);
  }

  int n_chr_types = -1;
  PopulationMap initial_population = parse_population_list_abm(initial_population_r, n_chr_types);
  FitnessMap untreated_fitness = parse_fitness_list_abm(untreated_fitness_map_r, n_chr_types, "untreated_fitness_map_r");
  FitnessMap treated_fitness = parse_fitness_list_abm(treated_fitness_map_r, n_chr_types, "treated_fitness_map_r");
  AdjacencyMap adjacency = parse_adjacency_list_abm(adjacency_r, n_chr_types);
  KaryotypeSet transition_set = parse_karyotype_set_abm(transition_karyotypes_r, n_chr_types, "transition_karyotypes_r");

  std::vector<std::string> metric_conditions;
  std::vector<int> metric_steps;
  std::vector<double> metric_times;
  std::vector<double> metric_totals;
  std::vector<double> metric_diversities;
  std::vector<double> metric_transition_totals;
  std::vector<double> metric_transition_fractions;

  int tau2_step = tau1_step;
  double best_transition_fraction = -1.0;
  Rcpp::List condition1 = simulate_transition_condition_abm(
    initial_population, untreated_fitness, treated_fitness, adjacency, transition_set,
    "condition1_first_treatment_only", tau1_step, -1, false, p_missegregation,
    base_death_rate, base_birth_rate, fitness_birth_scale, 0.0, dt, n_steps,
    record_interval, condition1_rng, metric_conditions, metric_steps, metric_times,
    metric_totals, metric_diversities, metric_transition_totals,
    metric_transition_fractions, /*track_tau2=*/true, &tau2_step, &best_transition_fraction);

  Rcpp::List condition2 = simulate_transition_condition_abm(
    initial_population, untreated_fitness, treated_fitness, adjacency, transition_set,
    "condition2_second_treatment", tau1_step, tau2_step, true, p_missegregation,
    base_death_rate, base_birth_rate, fitness_birth_scale, second_treatment_strength,
    dt, n_steps, record_interval, condition2_rng, metric_conditions, metric_steps,
    metric_times, metric_totals, metric_diversities, metric_transition_totals,
    metric_transition_fractions);

  Rcpp::DataFrame metrics = Rcpp::DataFrame::create(
    Rcpp::Named("condition") = metric_conditions,
    Rcpp::Named("step") = metric_steps,
    Rcpp::Named("time") = metric_times,
    Rcpp::Named("total_population") = metric_totals,
    Rcpp::Named("diversity") = metric_diversities,
    Rcpp::Named("transition_population") = metric_transition_totals,
    Rcpp::Named("transition_fraction") = metric_transition_fractions,
    Rcpp::Named("stringsAsFactors") = false);

  double c1_final_pop = 0.0;
  double c2_final_pop = 0.0;
  double c1_final_diversity = 0.0;
  double c2_final_diversity = 0.0;
  for (std::size_t i = 0; i < metric_conditions.size(); ++i) {
    if (metric_steps[i] == n_steps && metric_conditions[i] == "condition1_first_treatment_only") {
      c1_final_pop = metric_totals[i];
      c1_final_diversity = metric_diversities[i];
    }
    if (metric_steps[i] == n_steps && metric_conditions[i] == "condition2_second_treatment") {
      c2_final_pop = metric_totals[i];
      c2_final_diversity = metric_diversities[i];
    }
  }

  double population_reduction_pct = c1_final_pop > 0.0
    ? (1.0 - c2_final_pop / c1_final_pop) * 100.0
    : NA_REAL;
  double diversity_reduction_pct = c1_final_diversity > 0.0
    ? (1.0 - c2_final_diversity / c1_final_diversity) * 100.0
    : NA_REAL;

  Rcpp::DataFrame endpoint_summary = Rcpp::DataFrame::create(
    Rcpp::Named("tau2_step") = tau2_step,
    Rcpp::Named("tau2_time") = static_cast<double>(tau2_step) * dt,
    Rcpp::Named("tau2_transition_fraction") = best_transition_fraction,
    Rcpp::Named("condition1_final_population") = c1_final_pop,
    Rcpp::Named("condition2_final_population") = c2_final_pop,
    Rcpp::Named("population_reduction_pct") = population_reduction_pct,
    Rcpp::Named("condition1_final_diversity") = c1_final_diversity,
    Rcpp::Named("condition2_final_diversity") = c2_final_diversity,
    Rcpp::Named("diversity_reduction_pct") = diversity_reduction_pct);

  return Rcpp::List::create(
    Rcpp::Named("tau2_step") = tau2_step,
    Rcpp::Named("tau2_time") = static_cast<double>(tau2_step) * dt,
    Rcpp::Named("tau2_transition_fraction") = best_transition_fraction,
    Rcpp::Named("condition1") = condition1,
    Rcpp::Named("condition2") = condition2,
    Rcpp::Named("metrics") = metrics,
    Rcpp::Named("endpoint_summary") = endpoint_summary);
}
