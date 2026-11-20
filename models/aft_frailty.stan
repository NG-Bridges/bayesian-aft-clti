data {
  int<lower=1> N; int<lower=1> P; int<lower=1> Q; int<lower=1> G;
  vector<lower=0>[N] time; array[N] int<lower=0,upper=1> status;
  matrix[N,P] X; matrix[N,Q] D;
  vector[G] gh_x; vector[G] gh_logw;
  real<lower=0> intercept_scale; real<lower=0> shape_scale; real<lower=0> design_scale;
  real<lower=0> tau_0; real<lower=0> slab_scale; real<lower=0> slab_df;
  real<lower=0> frailty_scale;
}
transformed data { vector[N] log_t = log(time); }
parameters {
  real alpha_z; real log_k_z; vector[Q] bdes_z;
  vector[P] z; vector<lower=0>[P] lambda; real<lower=0> tau; real<lower=0> caux;
  real<lower=0> sigma_u;
}
transformed parameters {
  real alpha = alpha_z * intercept_scale;
  real<lower=0> k_shape = exp(log_k_z * shape_scale);
  vector[Q] beta_des = bdes_z * design_scale;
  real<lower=0> c2 = square(slab_scale) * caux;
  vector[P] lam2 = square(lambda);
  vector[P] beta = z .* sqrt(c2 * lam2 ./ (c2 + square(tau) * lam2)) * tau;
}
model {
  alpha_z ~ std_normal(); log_k_z ~ std_normal(); bdes_z ~ std_normal();
  z ~ std_normal(); lambda ~ cauchy(0, 1); tau ~ cauchy(0, tau_0);
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  sigma_u ~ normal(0, frailty_scale);
  {
    vector[N] r = log_t - (alpha + X * beta + D * beta_des);
    vector[G] cg = sqrt2() * sigma_u * gh_x;
    for (i in 1:N) {
      vector[G] lq;
      for (g in 1:G) {
        real e = k_shape * (r[i] - cg[g]);
        lq[g] = gh_logw[g] - exp(e) + (status[i] == 1 ? log(k_shape) + e - log_t[i] : 0);
      }
      target += log_sum_exp(lq);
    }
  }
}
generated quantities {
  vector[N] log_lik;
  {
    vector[N] r = log_t - (alpha + X * beta + D * beta_des);
    vector[G] cg = sqrt2() * sigma_u * gh_x;
    for (i in 1:N) {
      vector[G] lq;
      for (g in 1:G) {
        real e = k_shape * (r[i] - cg[g]);
        lq[g] = gh_logw[g] - exp(e) + (status[i] == 1 ? log(k_shape) + e - log_t[i] : 0);
      }
      log_lik[i] = log_sum_exp(lq);
    }
  }
}
