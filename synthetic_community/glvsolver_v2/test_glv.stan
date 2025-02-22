functions {
    vector row_wise_to_vector(matrix A) {
        int num_rows = rows(A);
        int num_cols = cols(A);
        vector[num_rows * num_cols] result;

        // Fill the result vector in row-wise order
        for (i in 1:num_rows) {
            for (j in 1:num_cols) {
                result[(i - 1) * num_cols + j] = A[i, j];
            }
        }
        return result;
    }
}
data {
	int<lower=0> N;
	vector[N] dlogX;
	matrix[N, 3] growth;
	matrix[N, 9] interaction;
	matrix[N, 6] perturbation;
}
parameters {
	real<lower=0, upper=100.00> sigma;
	vector[3] alpha;
	matrix[3, 3] beta;
	matrix[3, 2] epsilon;
}
transformed parameters {
    matrix[3, 3] constrained_beta;
    for (i in 1:3) {
        for (j in 1:3) {
            if (i == j) {
                constrained_beta[i, j] = -abs(beta[i, j]);
            } else {
                constrained_beta[i, j] = beta[i, j];
            }
        }
    }
}
model {
	sigma ~ uniform(0,100.00);
	alpha ~ normal(0,100.00);
	beta[1, 1] ~ normal(0,100.00);
	beta[1, 2] ~ normal(0,100.00);
	beta[1, 3] ~ normal(0,100.00);
	beta[2, 1] ~ normal(0,100.00);
	beta[2, 2] ~ normal(0,100.00);
	beta[2, 3] ~ normal(0,100.00);
	beta[3, 1] ~ normal(0,100.00);
	beta[3, 2] ~ normal(0,100.00);
	beta[3, 3] ~ normal(0,100.00);
	epsilon[1, 1] ~ normal(0,100.00);
	epsilon[1, 2] ~ normal(0,100.00);
	epsilon[2, 1] ~ normal(0,100.00);
	epsilon[2, 2] ~ normal(0,100.00);
	epsilon[3, 1] ~ normal(0,100.00);
	epsilon[3, 2] ~ normal(0,100.00);
	vector[N] mean = growth * alpha + interaction * row_wise_to_vector(beta) + perturbation * row_wise_to_vector(epsilon);
	dlogX ~ normal(mean, sigma);
}