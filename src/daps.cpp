#define ARMA_64BIT_WORD 1
#define ARMA_PRINT_CONFIG
#include <RcppArmadillo.h>
#include <unordered_set>
#include <algorithm>
#include <map>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

struct VectorHash {
  size_t operator()(const std::vector<int>& v) const {
    size_t seed = v.size();
    for (auto& i : v) {
      seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

class DAPSObject {
private:
  struct DAPSParam {
    int p;
    int L;
    arma::mat proposal;
    double proposal_thresh;
    double min_abs_corr;
    bool twas_weight;
  } par;

  struct MLRData {
    arma::mat XtX;
    arma::vec Xty;
    double yty;
    int n;
    arma::vec prior;
    arma::vec grid;
  } dat;

  struct ModelResult {
    vector<vector<int>> models;
    Rcpp::CharacterVector model_configs;
    int m;
    arma::vec log_prior;
    arma::vec log_bf;
    arma::vec log_posterior;
    arma::vec posterior_prob;
    vector<map<int, double>> reg_weights;
  } mod;

  struct SummaryResult {
    double log_nc;
    arma::vec pip;
    arma::vec twas_weights;
    arma::mat effect_pip;
  } rst;

  struct SignalCluster {
    vector<set<int>> clusters;
    vector<vector<int>> signal_clusters;
    arma::sp_mat LD;
    arma::vec cpip;
  } sc;


public:
  void parse_params(const arma::vec &prior,
                    const arma::mat &proposal,
                    double proposal_thresh,
                    const arma::vec &grid,
                    bool twas_weight,
                    double min_abs_corr,
                    int n) {
    dat.prior = prior;
    par.proposal = proposal;
    par.proposal_thresh = proposal_thresh;
    dat.grid = grid;
    par.twas_weight = twas_weight;
    par.min_abs_corr = min_abs_corr;
    dat.n = n;
    par.p = proposal.n_rows;
    par.L = proposal.n_cols;
    sc.clusters.resize(par.L);
  }


  void load_data(const arma::mat &X,
                 const arma::vec &y) {
    dat.XtX = X.t() * X;
    dat.Xty = X.t() * y;
    dat.yty = arma::dot(y, y);
  }
  void load_data_ss(const arma::mat &XtX,
                    const arma::vec &Xty,
                    const double yty) {
    dat.XtX = XtX;
    dat.Xty = Xty;
    dat.yty = yty;
  }


  void sampling() {
    unordered_set<vector<int>, VectorHash> unique_models; // set of unique models
    int max_models = static_cast<int>(1.0 / par.proposal_thresh) + par.p;
    unique_models.reserve(max_models);
    mod.models.reserve(max_models);

    double log_thresh = log10(par.proposal_thresh);
    arma::mat log_proposal = arma::log10(par.proposal);
    arma::umat sorted_indices(par.p, par.L); // sorted indices for each effect

    for (int l = 0; l < par.L; l++) {
      sorted_indices.col(l) = arma::sort_index(log_proposal.col(l), "descend");
    }

    vector<int> v(par.L, 0); // initialize model indicator vector
    int k = 0;
    while (k < par.L) {
      k = 0;
      vector<int> gamma(par.L);
      double log_q = 0.0;
      for (int l = 0; l < par.L; l++) {
        gamma[l] = sorted_indices(v[l], l);
        log_q += log_proposal(gamma[l], l);
        if (log_q < log_thresh) { // early stopping if below threshold
          break;
        }
      }

      if (log_q >= log_thresh) { // model passes proposal threshold
        vector<int> gamma_set;
        bool valid = true;
        for (int l = 0; l < par.L; l++) {
          if (gamma[l] != par.p - 1) { //skip null indicator
            if (find(gamma_set.begin(), gamma_set.end(), gamma[l]) != gamma_set.end()) {
              valid = false; // duplicate SNP in model
              break;
            } else {
              gamma_set.push_back(gamma[l]);
            }
          }
        }
        if (valid) { // no duplicate SNPs in model and unique model
          sort(gamma_set.begin(), gamma_set.end());
          if (unique_models.insert(gamma_set).second) {
            mod.models.push_back(gamma);
            for (int l = 0; l < par.L; l++) {
              if (gamma[l] != par.p - 1) {sc.clusters[l].insert(gamma[l]);}
            }
          }
        }
        while (k < par.L) { // increment
          v[k]++;
          if (v[k] < par.p) {
            break;
          } else {
            v[k] = 0;
            k++;
          }
        }
      } else {
        bool exit = false; //backtrack to find next valid model
        while (k < par.L && !exit) {
          if (v[k] > 0) {
            v[k] = 0;
            k++;
            while (k < par.L) {
              v[k]++;
              if (v[k] < par.p) {
                exit = true;
                break;
              } else {
                v[k] = 0;
                k++;
              }
            }
          } else {
            k++;
          } 
        }
      }
    }
    for (int i = 0; i < par.p; i++) { // add single SNP models and null model
      vector<int> gamma_set;
      if (i != par.p - 1) {
        gamma_set.push_back(i);
      }
      if (unique_models.insert(gamma_set).second) {
        vector<int> gamma(1, i);
        mod.models.push_back(gamma);
      }
    }
  }


  void posterior() {
    mod.m = mod.models.size();
    mod.log_prior.set_size(mod.m);
    mod.log_bf.set_size(mod.m);
    if (par.twas_weight) {mod.reg_weights.resize(mod.m); rst.twas_weights.set_size(par.p - 1);}

    arma::vec log_prior_diff = arma::log10(dat.prior) - arma::log10(1.0 - dat.prior);
    double log_prior_null = arma::sum(arma::log10(1.0 - dat.prior));
    arma::vec weights = arma::ones(dat.grid.n_elem) / dat.grid.n_elem;

    for (int i = 0; i < mod.m; i++) {
      // compute prior
      vector<int> snps;
      double lp = log_prior_null;
      for (int j : mod.models[i]) {
        if (j != par.p - 1) {
          snps.push_back(j);
          lp += log_prior_diff[j];
        }
      }
      mod.log_prior[i] = lp;

      // extract snp data
      arma::uvec idx = arma::conv_to<arma::uvec>::from(snps);
      arma::mat GtG = dat.XtX.submat(idx, idx);
      arma::vec Gty = dat.Xty.elem(idx);

      // compute BF with decomposition
      arma::vec eigval;
      arma::mat eigvec;
      arma::eig_sym(eigval, eigvec, GtG);
      arma::vec proj_y = eigvec.t() * Gty;
      arma::vec proj_y2 = arma::pow(proj_y, 2);
      // compute BF over grid
      arma::vec log_bf_vec(dat.grid.n_elem);
      for (int k = 0; k < dat.grid.n_elem; k++) {
        double phi2 = dat.grid(k);
        double log_det = arma::sum(arma::log10(1.0 + phi2 * eigval));
        double quad_form = phi2 * arma::sum(proj_y2 / (1.0 + phi2 * eigval));
        log_bf_vec(k) = -0.5 * log_det - 0.5 * dat.n * log10(1.0 - quad_form / dat.yty);
      }
      // average over grid
      double max_log_bf = arma::max(log_bf_vec);
      double weighted_sum = arma::sum(weights % arma::exp(log(10.0) * (log_bf_vec - max_log_bf)));
      mod.log_bf[i] = max_log_bf + log10(weighted_sum);
      // regression weights
      if (par.twas_weight) {
        arma::vec beta_hat = eigvec * (proj_y / (eigval + 1e-12));
        for (int j = 0; j < snps.size(); j++) {
          mod.reg_weights[i][snps[j]] = beta_hat(j);
        }
      }
    }
    mod.log_posterior = mod.log_prior + mod.log_bf;
    // compute posterior model probabilities
    double max_log_posterior = arma::max(mod.log_posterior);
    arma::vec exp_val = arma::exp(log(10.0) * (mod.log_posterior - max_log_posterior));
    double sum_exp_val = arma::sum(exp_val);
    mod.posterior_prob = exp_val / sum_exp_val;
    rst.log_nc = max_log_posterior + log10(sum_exp_val);
  }

  void summarize() {
    rst.pip.zeros(par.p - 1);
    rst.effect_pip.zeros(par.p - 1, par.L);
    mod.model_configs = Rcpp::CharacterVector(mod.m);

    // summarize the effect PIPs in L effects
    for (int i = 0; i < mod.m; i++) {
      vector<int>& model = mod.models[i];
      double post_prob = mod.posterior_prob(i);
      // save model configuration
      string config = "";

      bool first = true;
      for (int snp : model) {
        if (snp != par.p - 1) { // skip null indicator
          if (!first) {config += " + ";}
          config += std::to_string(snp + 1); first = false;
        }
      }
      mod.model_configs[i] = first ? "null" : config;
      // compute PIPs
      if (model.size() > 1) {
        for (int l = 0; l < par.L; l++) {
          int snp = model[l];
          if (snp != par.p - 1) { // skip null indicator
            rst.effect_pip(snp, l) += post_prob;
            rst.pip(snp) += post_prob;
            if (par.twas_weight) {rst.twas_weights(snp) += post_prob * mod.reg_weights[i][snp];}
          }
        }
      } else {
        int snp = model[0];
        if (snp != par.p - 1) { // skip null indicator
          for (int l = 0; l < par.L; l++) {
            if (sc.clusters[l].find(snp) != sc.clusters[l].end()) {
              rst.effect_pip(snp, l) += post_prob;
            }
          }
          rst.pip(snp) += post_prob;
          if (par.twas_weight) {rst.twas_weights(snp) += post_prob * mod.reg_weights[i][snp];}
        }
      }
    }
    rst.pip = arma::clamp(rst.pip, 0.0, 1.0);
    rst.effect_pip = arma::clamp(rst.effect_pip, 0.0, 1.0);
  }

  void signal_cluster() {
    sc.signal_clusters.clear();
    sc.cpip.zeros(par.L);
    sc.LD.set_size(par.p - 1, par.p - 1);

    // for each effect, construct signal cluster based on effect PIPs and LD
    for (int l = 0; l < par.L; l++) {
      // order SNPs in cluster by effect PIP
      vector<pair<int, double>> snp_pip_pairs;
      for (int snp : sc.clusters[l]) {
        snp_pip_pairs.push_back(make_pair(snp, rst.effect_pip(snp, l)));
      }
      sort(snp_pip_pairs.begin(), snp_pip_pairs.end(),
           [](const pair<int, double> &a, const pair<int, double> &b) {
             return a.second > b.second;
           });
      if(snp_pip_pairs.empty()) {
        sc.signal_clusters.push_back(vector<int>());
        continue;
      }
      // build signal cluster
      int lead_snp = snp_pip_pairs[0].first;
      vector<int> current_cluster = {lead_snp};
      double cluster_cpip = snp_pip_pairs[0].second;
      // remaining SNPs to consider
      for (size_t j = 1; j < snp_pip_pairs.size(); j++) {
        int next_snp = snp_pip_pairs[j].first;
        bool add_to_cluster = true;
        for (int cluster_snp : current_cluster) {
          double abs_corr = fabs(dat.XtX(next_snp, cluster_snp)) / sqrt(dat.XtX(next_snp, next_snp) * dat.XtX(cluster_snp, cluster_snp));
          if (abs_corr < par.min_abs_corr) {
            add_to_cluster = false;
            break;
          } else {sc.LD(cluster_snp, next_snp) = abs_corr; sc.LD(next_snp, cluster_snp) = abs_corr;}
        }
        if (add_to_cluster) {
          current_cluster.push_back(next_snp);
          cluster_cpip += snp_pip_pairs[j].second;
        }
      }
      sc.signal_clusters.push_back(current_cluster);
      sc.cpip(l) = min(cluster_cpip, 1.0);
    }
  }

  Rcpp::List get_output() const {
    return Rcpp::List::create(
      Named("model_index")   = mod.models,
      Named("model_configs") = mod.model_configs,
      Named("log_prior")     = mod.log_prior,
      Named("log_bf")        = mod.log_bf,
      Named("log_posterior") = mod.log_posterior,
      Named("posterior_prob")= mod.posterior_prob,
      Named("pip")           = rst.pip,
      Named("reg_weights")   = mod.reg_weights,
      Named("twas_weights")  = rst.twas_weights,
      Named("effect_pip")    = rst.effect_pip,
      Named("signal_clusters") = sc.signal_clusters,
      Named("cpip")          = sc.cpip,
      Named("LD")            = sc.LD,
      Named("min_abs_corr")  = par.min_abs_corr,
      Named("log_nc")        = rst.log_nc
    );
  }
};



//' DAPS main function in Rcpp
//'
//' @param X Genotype matrix (n x p)
//' @param y Phenotype vector (n x 1)
//' @param n Sample size (scalar)
//' @param prior Prior inclusion probabilities (p x 1)
//' @param proposal Proposal matrix (p x L) from SuSiE
//' @param proposal_thresh Proposal threshold (scalar)
//' @param grid Grid of scaled prior variance values (k x 1)
//' @param twas_weight Logical indicating whether to use TWAS weights
//' @param min_abs_corr Minimum absolute correlation for TWAS weights (scalar)
//' @return A list containing models and input data
//'
//' @export
// [[Rcpp::export]]
Rcpp::List daps_main(const arma::mat &X,
                     const arma::vec &y,
                     int n,
                     const arma::vec &prior,
                     const arma::mat &proposal,
                     double proposal_thresh,
                     const arma::vec &grid,
                     bool twas_weight,
                     double min_abs_corr) {

  DAPSObject daps;

  daps.parse_params(prior, proposal, proposal_thresh, grid, twas_weight, min_abs_corr, n);

  daps.load_data(X, y);

  daps.sampling();

  daps.posterior();

  daps.summarize();

  daps.signal_cluster();

  return daps.get_output();
}



//' DAPS_SS main function in Rcpp
//'
//' @param XtX Cross-product of genotype matrix (p x p)
//' @param Xty Cross-product of genotype and phenotype vector (p x 1)
//' @param yty Cross-product of phenotype vector (scalar)
//' @param n Sample size (scalar)
//' @param prior Prior inclusion probabilities (p x 1)
//' @param proposal Proposal matrix (p x L) from SuSiE
//' @param proposal_thresh Proposal threshold (scalar)
//' @param grid Grid of scaled prior variance values (k x 1)
//' @param twas_weight Logical indicating whether to use TWAS weights
//' @param min_abs_corr Minimum absolute correlation for TWAS weights (scalar)
//' @return A list containing models and input data
//'
//' @export
// [[Rcpp::export]]
Rcpp::List daps_ss_main(const arma::mat &XtX,
                        const arma::vec &Xty,
                        const double yty,
                        int n,
                        const arma::vec &prior,
                        const arma::mat &proposal,
                        double proposal_thresh,
                        const arma::vec &grid,
                        bool twas_weight,
                        double min_abs_corr) {

  DAPSObject daps;

  daps.parse_params(prior, proposal, proposal_thresh, grid, twas_weight, min_abs_corr, n);

  daps.load_data_ss(XtX, Xty, yty);
  
  daps.sampling();

  daps.posterior();

  daps.summarize();

  daps.signal_cluster();

  return daps.get_output();
}