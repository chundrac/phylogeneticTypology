require(rstan)
require(bayesplot)
require(loo)

model.loglik <- NULL
model.fam.loglik <- NULL
model.fam.geo.loglik <- NULL

beta.s.geo <- NULL
lon.all <- NULL
lat.all <- NULL

for (iteration in 1:10) {
  load(paste('output/fit_model.stan_',iteration,'.Rdata',sep=''))
  model.loglik <- rbind(model.loglik,loo::extract_log_lik(fit))
  
  load(paste('output/fit_model-fam.stan_',iteration,'.Rdata',sep=''))
  model.fam.loglik <- rbind(model.fam.loglik,loo::extract_log_lik(fit))
  
  load(paste('output/fit_model-fam-geo.stan_',iteration,'.Rdata',sep=''))
  model.fam.geo.loglik <- rbind(model.fam.geo.loglik,loo::extract_log_lik(fit))
  
  beta.s.geo <- rbind(beta.s.geo,as.matrix(fit,'beta_s_geo'))
  
  lon.all <- rbind(lon.all,as.matrix(fit,'lon_root'))
  lat.all <- rbind(lat.all,as.matrix(fit,'lat_root'))
  
}

colnames(lon.all) <- families
colnames(lat.all) <- families

rm(list=ls()[!(ls() %in% c('model.loglik','model.fam.loglik','model.fam.geo.loglik','beta.s.geo','lon.all','lat.all'))])

save.image(file='posterior-params.Rdata')

loo_compare(loo(model.loglik),loo(model.fam.loglik),loo(model.fam.geo.loglik))

wts = loo_model_weights(list(model.loglik,model.fam.loglik,model.fam.geo.loglik))
