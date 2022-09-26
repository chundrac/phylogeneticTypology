functions {
    vector gen_relaxed_dists(real[] rho, matrix ancestor_lens, int[] child, int B, int N) {
    vector[N] relaxed_dists = rep_vector(0,N);
    for (b in 1:B) {
      relaxed_dists[child[b]] = dot_product(ancestor_lens[b,],to_vector(rho));
    }
    return(relaxed_dists);
  }
  real geo_llik(int[,] mrca22, vector relaxed_dists, vector tip_vals, real root_val, int T, real sigma) {
    matrix[T,T] cov22;
    vector[T] mu2;for (i in 1:(T-1)) {
      cov22[i,i] = relaxed_dists[mrca22[i,i]];
      for (j in (i+1):T) {
        cov22[i,j] = relaxed_dists[mrca22[i,j]];
        cov22[j,i] = cov22[i,j];
      }
    }
    cov22[T,T] = relaxed_dists[mrca22[T,T]];
    mu2 = rep_vector(root_val,T);
    return(multi_normal_cholesky_lpdf(tip_vals|mu2,cholesky_decompose(sigma*cov22)));
  }
  matrix evprob(real z, real alpha, real beta) {
    matrix[2,2] P;
    P[1,1] = (beta/(alpha+beta)) + (alpha/(alpha+beta)*exp(-(alpha+beta)*z));
    P[1,2] = (alpha/(alpha+beta)) - (alpha/(alpha+beta)*exp(-(alpha+beta)*z));
    P[2,1] = (beta/(alpha+beta)) - (beta/(alpha+beta)*exp(-(alpha+beta)*z));
    P[2,2] = (alpha/(alpha+beta)) + (beta/(alpha+beta)*exp(-(alpha+beta)*z));
    return P;
  }
  //compute likelihood via Felsenstein's Pruning Algorithm
  vector pruning_vec(int N, int B, int[] child, int[] parent, real[] brlen, matrix tiplik, vector alpha, vector beta, int T, int[] roots) {
    vector[2] pi;                         //stationary probability
    matrix[N,2] lambda = log(tiplik);     //likelihoods at tips+nodes
    vector[T] lliks;                       //log likelihoods for each family for feature d
    for (b in 1:B) {
      matrix[2,2] P = evprob(brlen[b], alpha[b], beta[b]); //via matrix exponentiation
      for (d in 1:2) {
        lambda[parent[b],d] += log(dot_product(P[d],exp(lambda[child[b]])));
      }
    }
    for (t in 1:T) {
      pi[1] = -log(2) + lambda[parent[roots[t]],1];
      pi[2] = -log(2) + lambda[parent[roots[t]],2];
      lliks[t] = log_sum_exp(pi);
    }
    return(lliks);
  }
}
data {
  int<lower=1> N;                          //number of nodes in large families
  int<lower=1> P;                          //number of present geo values
  int<lower=1> M;                          //number of missing geo values, minus root
  int<lower=1> T;                          //number of large families
  int<lower=1> B;                          //number of branches
  int<lower=1> F;                          //number of families
  int<lower=1> D;                          //number of features
  int<lower=1> n_present[T];               //number of present lon values per large fam
  int<lower=1> n_missing[T];               //number of absent lon values per large fam, minus root
  int<lower=1> child[B];                //child of each branch
  int<lower=1> parent[B];               //parent of each branch
  real<lower=0> brlen[B];               //length of each branch
  int<lower=1,upper=F> fam_ID[B];       //family ID of each branch
  matrix[N,D*2] tiplik;                 //likelihoods for data at tips+internal nodes in of each tree
  vector[P] lon;
  vector[P] lat;
  matrix[B,B] ancestor_lens;
  int mrca22[P,P];
  int roots[T];                         //indices of root for each tree in parent list
}
parameters {
  real<lower=0> sigma;
  real<lower=0> rho[B];
  real lon_root[T];
  real lat_root[T];
  real<lower=0> tau;
  real alpha_s[D];
  real alpha_p[D];
  real beta_s_fam[F,D];
  real beta_p_fam[F,D];
  real beta_s_geo;
}
transformed parameters {
  vector[D*T] log_lik;
  vector[B] s[D];
  vector[B] p[D];
  vector[N] relaxed_dists = gen_relaxed_dists(rho, ancestor_lens, child, B, N);
  vector[T] tip_lon_loglik;
  vector[T] tip_lat_loglik;
  { int j = 1;
    for (t in 1:T) {
      tip_lon_loglik[t] = geo_llik(mrca22[j:(j+n_present[t]-1),j:(j+n_present[t]-1)],relaxed_dists,lon[j:(j+n_present[t]-1)],lon_root[t],n_present[t],sigma);
      tip_lat_loglik[t] = geo_llik(mrca22[j:(j+n_present[t]-1),j:(j+n_present[t]-1)],relaxed_dists,lat[j:(j+n_present[t]-1)],lat_root[t],n_present[t],sigma);
      j += n_present[t];
    }
  }
  { int k = 1;
  for (d in 1:D) {
	  for (b in 1:B) {
      s[d,b] = inv_logit(alpha_s[d] + beta_s_fam[fam_ID[b],d] + beta_s_geo*rho[b])*tau;
      p[d,b] = inv_logit(alpha_p[d] + beta_p_fam[fam_ID[b],d]);
    }
    log_lik[k:(k+T-1)] = pruning_vec(N,B,child,parent,brlen,tiplik[,((2*d)-1):(2*d)],p[d,].*s[d,],(1-p[d,]).*s[d,],T,roots);
    k += T;
  }
  }
}
model {
  rho ~ gamma(1,1);
  sigma ~ gamma(1,1);
  lon_root ~ uniform(-180,180);
  lat_root ~ uniform(-90,90);
  tau ~ uniform(0,10);
  alpha_s ~ normal(0,1);
  alpha_p ~ normal(0,1);
  for (i in 1:F) {
    beta_s_fam[i,] ~ normal(0,1);
    beta_p_fam[i,] ~ normal(0,1);
  }
  beta_s_geo ~ normal(0,1);
  target += sum(tip_lon_loglik);
  target += sum(tip_lat_loglik);
  target += sum(log_lik);
}
