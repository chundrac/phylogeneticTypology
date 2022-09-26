functions {
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
  real<lower=0> tau;
  real alpha_s[D];
  real alpha_p[D];
  real beta_s_fam[F,D];
  real beta_p_fam[F,D];
}
transformed parameters {
  vector[D*T] log_lik;
  vector[B] s[D];
  vector[B] p[D];
  { int k = 1;
  for (d in 1:D) {
	  for (b in 1:B) {
      s[d,b] = inv_logit(alpha_s[d] + beta_s_fam[fam_ID[b],d])*tau;
      p[d,b] = inv_logit(alpha_p[d] + beta_p_fam[fam_ID[b],d]);
    }
    log_lik[k:(k+T-1)] = pruning_vec(N,B,child,parent,brlen,tiplik[,((2*d)-1):(2*d)],p[d,].*s[d,],(1-p[d,]).*s[d,],T,roots);
    k += T;
  }
  }
}
model {
  tau ~ uniform(0,10);
  alpha_s ~ normal(0,1);
  alpha_p ~ normal(0,1);
  for (i in 1:F) {
    beta_s_fam[i,] ~ normal(0,1);
    beta_p_fam[i,] ~ normal(0,1);
  }
  target += sum(log_lik);
}
