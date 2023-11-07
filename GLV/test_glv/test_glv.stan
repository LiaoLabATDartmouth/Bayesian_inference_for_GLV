data {
	int<lower=0> N;
	vector[N] dlogX;
	vector[N] growth_Taxon_1;
	vector[N] interaction_Taxon_1_Taxon_1;
	vector[N] interaction_Taxon_1_Taxon_2;
	vector[N] interaction_Taxon_1_Taxon_3;
	vector[N] perturbation_Taxon_1_Perturbation_Periodic;
	vector[N] perturbation_Taxon_1_Perturbation_Square;
	vector[N] growth_Taxon_2;
	vector[N] interaction_Taxon_2_Taxon_1;
	vector[N] interaction_Taxon_2_Taxon_2;
	vector[N] interaction_Taxon_2_Taxon_3;
	vector[N] perturbation_Taxon_2_Perturbation_Periodic;
	vector[N] perturbation_Taxon_2_Perturbation_Square;
	vector[N] growth_Taxon_3;
	vector[N] interaction_Taxon_3_Taxon_1;
	vector[N] interaction_Taxon_3_Taxon_2;
	vector[N] interaction_Taxon_3_Taxon_3;
	vector[N] perturbation_Taxon_3_Perturbation_Periodic;
	vector[N] perturbation_Taxon_3_Perturbation_Square;
}
parameters {
	real<lower=0,upper=14.04> sigma;
	real alpha_Taxon_1;
	real<upper=0> beta_Taxon_1_Taxon_1;
	real beta_Taxon_1_Taxon_2;
	real beta_Taxon_1_Taxon_3;
	real epsilon_Taxon_1_Perturbation_Periodic;
	real epsilon_Taxon_1_Perturbation_Square;
	real alpha_Taxon_2;
	real beta_Taxon_2_Taxon_1;
	real<upper=0> beta_Taxon_2_Taxon_2;
	real beta_Taxon_2_Taxon_3;
	real epsilon_Taxon_2_Perturbation_Periodic;
	real epsilon_Taxon_2_Perturbation_Square;
	real alpha_Taxon_3;
	real beta_Taxon_3_Taxon_1;
	real beta_Taxon_3_Taxon_2;
	real<upper=0> beta_Taxon_3_Taxon_3;
	real epsilon_Taxon_3_Perturbation_Periodic;
	real epsilon_Taxon_3_Perturbation_Square;
}
model {
	sigma ~ uniform(0,14.04);
	alpha_Taxon_1 ~ normal(0,100.00);
	beta_Taxon_1_Taxon_1 ~ normal(0,100.00);
	beta_Taxon_1_Taxon_2 ~ normal(0,100.00);
	beta_Taxon_1_Taxon_3 ~ normal(0,100.00);
	epsilon_Taxon_1_Perturbation_Periodic ~ normal(0,100.00);
	epsilon_Taxon_1_Perturbation_Square ~ normal(0,100.00);
	alpha_Taxon_2 ~ normal(0,100.00);
	beta_Taxon_2_Taxon_1 ~ normal(0,100.00);
	beta_Taxon_2_Taxon_2 ~ normal(0,100.00);
	beta_Taxon_2_Taxon_3 ~ normal(0,100.00);
	epsilon_Taxon_2_Perturbation_Periodic ~ normal(0,100.00);
	epsilon_Taxon_2_Perturbation_Square ~ normal(0,100.00);
	alpha_Taxon_3 ~ normal(0,100.00);
	beta_Taxon_3_Taxon_1 ~ normal(0,100.00);
	beta_Taxon_3_Taxon_2 ~ normal(0,100.00);
	beta_Taxon_3_Taxon_3 ~ normal(0,100.00);
	epsilon_Taxon_3_Perturbation_Periodic ~ normal(0,100.00);
	epsilon_Taxon_3_Perturbation_Square ~ normal(0,100.00);
	dlogX ~ normal(alpha_Taxon_1*growth_Taxon_1+beta_Taxon_1_Taxon_1*interaction_Taxon_1_Taxon_1+beta_Taxon_1_Taxon_2*interaction_Taxon_1_Taxon_2+beta_Taxon_1_Taxon_3*interaction_Taxon_1_Taxon_3+epsilon_Taxon_1_Perturbation_Periodic*perturbation_Taxon_1_Perturbation_Periodic+epsilon_Taxon_1_Perturbation_Square*perturbation_Taxon_1_Perturbation_Square+alpha_Taxon_2*growth_Taxon_2+beta_Taxon_2_Taxon_1*interaction_Taxon_2_Taxon_1+beta_Taxon_2_Taxon_2*interaction_Taxon_2_Taxon_2+beta_Taxon_2_Taxon_3*interaction_Taxon_2_Taxon_3+epsilon_Taxon_2_Perturbation_Periodic*perturbation_Taxon_2_Perturbation_Periodic+epsilon_Taxon_2_Perturbation_Square*perturbation_Taxon_2_Perturbation_Square+alpha_Taxon_3*growth_Taxon_3+beta_Taxon_3_Taxon_1*interaction_Taxon_3_Taxon_1+beta_Taxon_3_Taxon_2*interaction_Taxon_3_Taxon_2+beta_Taxon_3_Taxon_3*interaction_Taxon_3_Taxon_3+epsilon_Taxon_3_Perturbation_Periodic*perturbation_Taxon_3_Perturbation_Periodic+epsilon_Taxon_3_Perturbation_Square*perturbation_Taxon_3_Perturbation_Square, sigma);
}