data {
  int<lower=1> N; int<lower=2> K; int<lower=1> P; int<lower=1> Q; int<lower=1> Pw;
  vector<lower=0>[N] time; array[N] int<lower=0,upper=1> status;
  matrix[N,P] X; matrix[N,Q] D; matrix[N,Pw] W;
  real<lower=0> intercept_scale; real<lower=0> shape_scale; real<lower=0> design_scale;
  real<lower=0> gate_scale; real<lower=0> tau_0; real<lower=0> slab_scale; real<lower=0> slab_df;
}
transformed data { vector[N] log_t = log(time); }
parameters {
  vector[K] alpha_z; vector[K] log_k_z; matrix[Pw,K] gamma_z;
  array[K] vector[P] z; array[K] vector<lower=0>[P] lambda; array[K] real<lower=0> tau;
  real<lower=0> caux; array[K] vector[Q] bdes_z;
}
transformed parameters {
  vector[K] alpha; vector<lower=0>[K] k_shape;
  array[K] vector[P] beta; array[K] vector[Q] beta_des;
  matrix[Pw,K] gamma_full;
  real<lower=0> c2 = square(slab_scale) * caux;
  for (k in 1:K) {
    alpha[k]   = alpha_z[k] * intercept_scale;
    k_shape[k] = exp(log_k_z[k] * shape_scale);
    beta_des[k] = bdes_z[k] * design_scale;
    {
      vector[P] lam2 = square(lambda[k]);
      beta[k] = z[k] .* sqrt(c2 * lam2 ./ (c2 + square(tau[k]) * lam2)) * tau[k];
    }
  }
  for (v in 1:Pw) {
    row_vector[K] g = gamma_z[v] * gate_scale;
    gamma_full[v] = g - mean(g);
  }
}
model {
  alpha_z ~ std_normal(); log_k_z ~ std_normal(); to_vector(gamma_z) ~ std_normal();
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  for (k in 1:K) {
    z[k] ~ std_normal(); lambda[k] ~ cauchy(0, 1); tau[k] ~ cauchy(0, tau_0);
    bdes_z[k] ~ std_normal();
  }
  {
    matrix[N,K] LP; matrix[N,K] gate = W * gamma_full;
    for (k in 1:K) LP[, k] = alpha[k] + X * beta[k] + D * beta_des[k];
    for (i in 1:N) {
      vector[K] logpi = log_softmax(gate[i]');
      vector[K] lpc;
      for (k in 1:K) {
        real e = k_shape[k] * (log_t[i] - LP[i, k]);
        lpc[k] = logpi[k] - exp(e) + (status[i] == 1 ? log(k_shape[k]) + e - log_t[i] : 0);
      }
      target += log_sum_exp(lpc);
    }
  }
}
generated quantities {
  vector[N] log_lik;
  {
    matrix[N,K] LP; matrix[N,K] gate = W * gamma_full;
    for (k in 1:K) LP[, k] = alpha[k] + X * beta[k] + D * beta_des[k];
    for (i in 1:N) {
      vector[K] logpi = log_softmax(gate[i]');
      vector[K] lpc;
      for (k in 1:K) {
        real e = k_shape[k] * (log_t[i] - LP[i, k]);
        lpc[k] = logpi[k] - exp(e) + (status[i] == 1 ? log(k_shape[k]) + e - log_t[i] : 0);
      }
      log_lik[i] = log_sum_exp(lpc);
    }
  }
}
