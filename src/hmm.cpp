// Implementation of the forwards-backwards algorithm for
// parameter estimation of Gaussian density HMM.

#include <R.h>
#include <Rmath.h>
#include <vector>

using namespace std;

const double really_small_num = -1000000000.;

double
calc_observed_likelihood_iter(vector< vector<double> > &yll, 
							  vector< vector<double> > &tpm,
							  double *pi, bool do_filter, double *filter,
							  vector< vector<double> > &alpha,
							  vector< vector<double> > &beta,
							  bool print_info)
{

	// Initialization

	int T = yll[0].size(), k = yll.size();
	double sum_log_lik = 0, ll = 0;

	for(int i = 0; i < k; i++)
	{
		
		ll += exp(alpha[i][0] = yll[i][0] + pi[i]);
		beta[i][T - 1] = 0;
		if (print_info)
			Rprintf("yll[%d][0] = %f\tpi[%d] = %f\n", i, yll[i][0], i, pi[i]);

	}
	sum_log_lik = log(ll);
	for(int i = 0; i < k; i++)
	{
		
		alpha[i][0] -= sum_log_lik;
		if (print_info)
			Rprintf("\talpha[%d][0] = %f\tbeta[%d][%d] = %f\n",
					i, alpha[i][0], i, T - 1, beta[i][T - 1]);
		if (do_filter)
			filter[i] = exp(alpha[i][0]);

	}
	
	// Calculate the observed log-likelihood using scaling variables and
	// a dynamic programming table.

	for(int m = 1, t = T - 2, offs = k; m < T; m++, t--, offs += k)
	{

		double avf = 0, avb = 0;
		
		for(int i = 0; i < k; i++)
		{

			beta[i][t] = alpha[i][m] = 0;
			for(int j = 0; j < k; j++)
			{
				
				alpha[i][m] +=
					exp(yll[i][m] + alpha[j][m - 1] + tpm[j][i]);
				beta[i][t] +=
					exp(yll[i][t] + beta[j][t + 1] + tpm[i][j]);

			}
			avf += alpha[i][m], avb += beta[i][t];

		}
		avf = log(avf), avb = log(avb), ll = 0;
		for(int i = 0; i < k; i++)
		{
			
			alpha[i][m] = log(alpha[i][m]) - avf;
			beta[i][t] = log(beta[i][t]) - avb;
			if (do_filter)
				filter[offs + i] = exp(alpha[i][m]);
			if (print_info)
				Rprintf("\talpha[%d][%d] = %f\tbeta[%d][%d] = %f",
						i, m, alpha[i][m], i, t, beta[i][t]);
			
		}
		sum_log_lik += avf;
		if (print_info)
			Rprintf("\tavf = %f\t%f\n", avf, sum_log_lik);

	}
	
	// compute the log-likelihood including correction factor
	
// 	ll = 0;
// 	for(int i = 0; i < k; i++)
// 		ll += exp(f_tilde_prev[i]);
	
// 	return sum_log_lik + log(ll);
	return sum_log_lik;

}

double
calc_complete_likelihood(vector< vector<double> > &yll, vector< vector<double> > &tpm,
						 double *pi, int *states)
{

	int T = yll[0].size(), k = yll.size();
	vector< vector<double> > log_lik_table(k, vector<double>(T));
	vector< vector<unsigned int> > ptrs(k, vector<unsigned int>(T));
	
	// Initialization

//	Calculate the complete log-likelihood using a dynamic programming
//	table

	vector<double> v(k);
	
	for(int i = 0; i < k; i++)
		log_lik_table[i][0] = yll[i][0] + pi[i];
	for(int m = 1; m < T; m++)
		for(int i = 0; i < k; i++)
		{

			int best_ind = 0;
			double s = really_small_num;
			
			for(int j = 0; j < k; j++)
			{
				
				v[j] = log_lik_table[j][m - 1] + tpm[j][i];
				if (v[j] > s)
					s = v[j], best_ind = j;

			}
			log_lik_table[i][m] = yll[i][m] + v[ptrs[i][m] = best_ind];

		}

//	Track back the Viterbi path using the pointers stored in the
//	previous step
	
	double s = log_lik_table[0][T - 1];
	
	states[T - 1] = 0;
	for(int i = 1; i < k; i++)
		if (log_lik_table[i][T - 1] > s)
			s = log_lik_table[i][T - 1], states[T - 1] = i;
// 	Rprintf("%d\n", hidden_states[T - 1]);
	for(int m = T - 2; m >= 0; m--)
		states[m] = ptrs[states[m + 1]][m + 1];
	
	return log_lik_table[states[T - 1]][T - 1];

}

extern "C"
{

	void
	calc_observed_likelihood(int *_T, double *_y, int *_k, double *means,
							 double *sigma, double *TPM, double *pi, int *maxiter,
							 double *eps, double *_log_lik, double *filter,
							 int *hidden_states, bool *print_info)
	{

		// Initialization
		
		int k = *_k, T = *_T;
		vector< vector<double> > tpm(k, vector<double>(k)),
			yll(k, vector<double>(T)), alpha(k, vector<double>(T)),
			beta(k, vector<double>(T)), gamma(k, vector<double>(T)),
			ksi(k, vector<double>(k));
		vector<double> old_pi(k), old_means(k);

		for(int i = 0; i < k; i++)
			for(int j = 0, offs = 0; j < k; j++, offs += k)
				tpm[i][j] = log(TPM[offs + i]);

		*_log_lik = really_small_num;
		int iter = 0;
		
		while(iter < *maxiter)
		{

			double log_sigma = log(*sigma);
			
			for(int m = 0; m < T; m++)
				for(int i = 0; i < k; i++)
				{

					double x = (_y[m] - means[i]) / *sigma;
					
					yll[i][m] = -(M_LN_SQRT_2PI + 0.5 * x * x + log_sigma);
// 					Rprintf("\tyll[%d][%d] = %f", i, m, yll[i][m]);
					
				}
// 			Rprintf("\n");
			
			// Calculate alpha and beta
			
			double log_lik =
				calc_observed_likelihood_iter(yll, tpm, pi, false, filter,
											  alpha, beta, *print_info);
// 			Rprintf("\t*_log_lik = %f\tlog_lik = %f\n", *_log_lik, log_lik);

			if (log_lik < *_log_lik + *eps)
				break; // time to compute the filtered probs and exit
			else
				*_log_lik = log_lik;
			
			// Calculate gamma and ksi

			vector<double> sum_gamma(k, 0);
			vector< vector<double> > sum_ksi(k, vector<double>(k, 0));
			
			for(int m = 0; m < T; m++)
			{

				double avfb = 0, avksi = 0;
				
				for(int i = 0; i < k; i++)
				{
					
					avfb += exp(gamma[i][m] = alpha[i][m] + beta[i][m]);
					if (m < T - 1)
						for(int j = 0; j < k; j++)
							avksi +=
								exp(ksi[i][j] = alpha[i][m] + tpm[i][j] +
									yll[j][m + 1] + beta[j][m + 1]);

				}
				avfb = log(avfb), avksi = log(avksi);
				for(int i = 0; i < k; i++)
				{
					
					sum_gamma[i] += (gamma[i][m] = exp(gamma[i][m] - avfb));
					if (m < T - 1)
						for(int j = 0; j < k; j++)
							sum_ksi[i][j] += exp(ksi[i][j] - avksi);
// 					Rprintf("\tgamma[%d][%d] = %f", i, m, gamma[i][m]);
					
				}
// 				Rprintf("\n");
				
			}
			
			// Re-estimate the parameters
			// 1. pi
			
			for(int i = 0; i < k; i++)
			{
				
				old_pi[i] = pi[i];
				pi[i] = gamma[i][0];

			}
			
			// 2. tpm (a in Rabiner)
			
// 			Rprintf("tpm:\n");
			for(int i = 0; i < k; i++)
			{
				
				for(int j = 0; j < k; j++)
				{
					
					tpm[i][j] =
						log(sum_ksi[i][j] / (sum_gamma[i] - gamma[i][T - 1]));
// 					Rprintf("\t%f", tpm[i][j]);
					
				}
// 				Rprintf("\n");
				old_means[i] = means[i];
				means[i] = 0; // reset the means

			}
			
			// 3. emission prob. parameters

			double sum_sigma = 0;
			
			for(int m = 0; m < T; m++)
				for(int i = 0; i < k; i++)
				{

					double x = _y[m] - old_means[i];
					
					means[i] += gamma[i][m] * _y[m];
					sum_sigma += gamma[i][m] * x * x;

				}
			for(int i = 0; i < k; i++)
			{
				
				means[i] /= sum_gamma[i];
// 				Rprintf("\tmeans[%d] = %f", i, means[i]);

			}
			*sigma = sqrt(sum_sigma / T);
// 			Rprintf("\tsigma = %f\tpi[0] = %f\tpi[1] = %f\n", sigma, pi[0], pi[1]);

			iter++;
			
		}
		*maxiter = iter + 1;
		for(int i = 0; i < k; i++)
			for(int j = 0, offs = 0; j < k; j++, offs += k)
				TPM[offs + i] = tpm[i][j];

		// Compute the filtered probs and exit

// 		*_log_lik =
// 			calc_observed_likelihood_iter(yll, tpm, pi, true, filter, alpha, beta, *print_info);
		for(int m = 0, offs = 0; m < T; m++, offs += k)
			for(int i = 0; i < k; i++)
				filter[offs + i] = gamma[i][m];
		
		// Decode the hidden states using Viterbi alg.
		
		calc_complete_likelihood(yll, tpm, pi, hidden_states);
		
	}

}
