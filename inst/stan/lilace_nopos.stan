functions {
  real my_dirichlet_multinomial_lpmf(array[,] int y_arr, matrix alpha) {
    int N = dims(y_arr)[1];
    int K = dims(y_arr)[2];
    matrix[N, K] y = to_matrix(y_arr);
    vector[K] v = rep_vector(1.0, K);
    vector[N] alpha_0 = alpha * v;
    vector[N] n = y * v;
    vector[N] y_rowsum_1  = lgamma(y+1) * v;
    vector[N] y_rowsum_alpha  = lgamma(y+alpha) * v;
    vector[N] alpha_rowsum = lgamma(alpha) * v;
    return sum(lgamma(alpha_0) - lgamma(n + alpha_0) + lgamma(n + 1) -
        y_rowsum_1 + y_rowsum_alpha - alpha_rowsum);
  }
}
data {
  int<lower=0> V; // # of baseline variants
  int<lower=0> S; // # of synonymous variants
  int<lower=0> N; // # of observation
  int<lower=0> N_syn; // # of syn observation
  int<lower=0> P; // # of positions
  int<lower=0> R; // # of reps
  array[N] int nMAPv;
  array[N] int nMAPr;
  array[V] int vMAPp;
  array[V] int vMAPs;
  vector[N] n_counts; // variant counts vector
  int<lower=0> K; // # of bins
  array[N, K] int<lower=0> y; // FACS bin count

  array[N_syn] int sMAPr;
  vector[N_syn] n_syn_counts;
  array[N_syn, K] int<lower=0> y_syn; // FACS bin count
}
parameters {
  array[R] simplex[K] q; // simplex for every rep
  vector<lower=0> [V] sigma;
  vector[V] theta;
  vector[V] mu;
  vector<lower=0> [R] a;
  vector<lower=0> [R] b;
}
model { 
  // priors
  theta ~ std_normal();
  sigma ~ inv_gamma(1, 1);
  mu ~ normal(theta, sigma);
  for (r in 1:R) {
    q[r] ~ uniform(0,1);
  }

  // overdispersion
  vector[N] phi;
  phi = a[nMAPr] .* n_counts + b[nMAPr];
  matrix [N, K] phi_mat = rep_matrix(phi, K);

  vector[N_syn] phi_syn;
  phi_syn = a[sMAPr] .* n_syn_counts + b[sMAPr];
  matrix [N_syn, K] phi_syn_mat = rep_matrix(phi_syn, K);
  // baseline 
  array[R] row_vector[K] q_copy;
  for (r in 1:R) {
    q_copy[r] = to_row_vector(q[r]);
  }
  y_syn ~ my_dirichlet_multinomial(phi_syn_mat .* to_matrix(q_copy[sMAPr]));
  array[R] row_vector[K - 1] t0; // latent cutoff
  for (r in 1:R) {
    t0[r] = to_row_vector(inv_Phi(head(cumulative_sum(q[r]), K - 1)));
  }
  
  // latent normal transform
  matrix[N, K-1] t2cdf;
  matrix[N, K] p0;
  t2cdf = Phi(to_matrix(t0[nMAPr]) - rep_matrix(mu[nMAPv], K-1));
  p0[,1] = t2cdf[,1];
  for (k in 2:(K - 1)) {
    p0[,k] = t2cdf[,k] - t2cdf[,k - 1];
  }
  p0[,K] = 1 - t2cdf[,K - 1];

  // sample
  y ~ my_dirichlet_multinomial(phi_mat .* p0);
}
