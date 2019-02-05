#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>

#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>

#include "correlation_function.h"
#include "Stopwatch.h"
#include "Arsenal.h"

using namespace std;

void Correlation_function::Fit_correlation_function()
{
	err << "--> Getting HBT radii by Gaussian fit method" << endl;
	for (int iKT = 0; iKT < n_KT_bins; ++iKT)
	for (int iKphi = 0; iKphi < n_Kphi_bins; ++iKphi)
	for (int iKL = 0; iKL < n_KL_bins; ++iKL)
	{
		out << "  --> Fitting in pair-momentum bin: "
			<< "KT=[" << KT_pts[iKT] << ", " << KT_pts[iKT+1] << "] GeV, "
			<< "Kphi=[" << Kphi_pts[iKphi] << ", " << Kphi_pts[iKphi+1] << "], "
			<< "KL=[" << KL_pts[iKL] << ", " << KL_pts[iKL+1] << "] GeV"
			<< endl;
		find_minimum_chisq_correlationfunction_full( iKT, iKphi, iKL );
	}

	err << "--> Finished getting HBT radii by Gaussian fit method" << endl;

	return;
}


void Correlation_function::find_minimum_chisq_correlationfunction_full( int iKT, int iKphi, int iKL )
{

	const size_t data_length = n_qo_bins*n_qs_bins*n_ql_bins;  // # of points

    double lambda, R_o, R_s, R_l, R_os, R_ol, R_sl;
    int dim = 7;
    int s_gsl;

    double * V = new double [dim];
    double * qweight = new double [dim];
    double ** T = new double * [dim];
    for(int i = 0; i < dim; i++)
    {
        V[i] = 0.0;
        T[i] = new double [dim];
        for(int j = 0; j < dim; j++)
            T[i][j] = 0.0;
    }

	gsl_matrix * T_gsl = gsl_matrix_alloc (dim, dim);
	gsl_matrix * T_inverse_gsl = gsl_matrix_alloc (dim, dim);
	gsl_permutation * perm = gsl_permutation_alloc (dim);

	for (int i = 0; i < n_qo_bins; i++)
	for (int j = 0; j < n_qs_bins; j++)
	for (int k = 0; k < n_ql_bins; k++)
    {
		int idx = indexer(iKT, iKphi, iKL, i, j, k);

        double q_out_local = 0.5*(qo_pts[i]+qo_pts[i+1]);
        double q_side_local = 0.5*(qs_pts[j]+qs_pts[j+1]);
        double q_long_local = 0.5*(ql_pts[k]+ql_pts[k+1]);

		double correl_local = correlation_function[idx]-1.0;
		double correl_err_local = correlation_function_error[idx];

		if(correl_local < 1.0e-15) continue;

		bool ignore_central_point = true;
		if ( 	ignore_central_point
				and i==(n_qo_bins-1)/2
				and j==(n_qs_bins-1)/2
				and k==(n_ql_bins-1)/2)
			continue;
//			correl_err_local = 1.0e10;	//ignore central point
        double sigma_k_prime = correl_err_local/correl_local;

        double inv_sigma_k_prime_sq = 1./(sigma_k_prime*sigma_k_prime);
        double log_correl_over_sigma_sq = log(correl_local)*inv_sigma_k_prime_sq;

		qweight[0] = - 1.0;
		qweight[1] = q_out_local*q_out_local;
		qweight[2] = q_side_local*q_side_local;
		qweight[3] = q_long_local*q_long_local;
		qweight[4] = q_out_local*q_side_local;
		qweight[5] = q_out_local*q_long_local;
		qweight[6] = q_side_local*q_long_local;

        for(int ij = 0; ij < dim; ij++)
        {
            V[ij] += qweight[ij]*log_correl_over_sigma_sq;
            T[0][ij] += qweight[ij]*inv_sigma_k_prime_sq;
        }

        for(int ij = 1; ij < dim; ij++)
            T[ij][0] = T[0][ij];
            
        for(int ij = 1; ij < dim; ij++)
        {
            for(int lm = 1; lm < dim; lm++)
                T[ij][lm] += -qweight[ij]*qweight[lm]*inv_sigma_k_prime_sq;
        }
    }

    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
            gsl_matrix_set(T_gsl, i, j, T[i][j]);

    // Make LU decomposition of matrix T_gsl
    gsl_linalg_LU_decomp (T_gsl, perm, &s_gsl);
    // Invert the matrix m
    gsl_linalg_LU_invert (T_gsl, perm, T_inverse_gsl);

    double ** T_inverse = new double * [dim];
    for(int i = 0; i < dim; i++)
    {
        T_inverse[i] = new double [dim];
        for(int j = 0; j < dim; j++)
            T_inverse[i][j] = gsl_matrix_get(T_inverse_gsl, i, j);
    }
    double * results = new double [dim];
    for(int i = 0; i < dim; i++)
    {
        results[i] = 0.0;
        for(int j = 0; j < dim; j++)
            results[i] += T_inverse[i][j]*V[j];
    }

	// compute results
	lambda_Correl[indexerK(iKT, iKphi, iKL)] = exp(results[0]);
	R2_out[indexerK(iKT, iKphi, iKL)] = results[1]*hbarC*hbarC;
	R2_side[indexerK(iKT, iKphi, iKL)] = results[2]*hbarC*hbarC;
	R2_long[indexerK(iKT, iKphi, iKL)] = results[3]*hbarC*hbarC;
	R2_outside[indexerK(iKT, iKphi, iKL)] = results[4]*hbarC*hbarC;
	R2_outlong[indexerK(iKT, iKphi, iKL)] = results[5]*hbarC*hbarC;
	R2_sidelong[indexerK(iKT, iKphi, iKL)] = results[6]*hbarC*hbarC;

	// compute chi^2
    double chi_sq = 0.0;
	for (int i = 0; i < n_qo_bins; i++)
	for (int j = 0; j < n_qs_bins; j++)
	for (int k = 0; k < n_ql_bins; k++)
    {
		int idx = indexer(iKT, iKphi, iKL, i, j, k);

        double q_out_local = 0.5*(qo_pts[i]+qo_pts[i+1]);
        double q_side_local = 0.5*(qs_pts[j]+qs_pts[j+1]);
        double q_long_local = 0.5*(ql_pts[k]+ql_pts[k+1]);

		double correl_local = correlation_function[idx]-1.0;
		double correl_err_local = correlation_function_error[idx];

		if(correl_local < 1.0e-15) continue;

		bool ignore_central_point = true;
		if ( 	ignore_central_point
				and i==(n_qo_bins-1)/2
				and j==(n_qs_bins-1)/2
				and k==(n_ql_bins-1)/2)
			continue;
//			correl_err_local = 1.0e10;	//ignore central point
        double sigma_k_prime = correl_err_local/correl_local;

        chi_sq += pow( ( log(correl_local) - results[0] 
						+ results[1]*q_out_local*q_out_local 
						+ results[2]*q_side_local*q_side_local
						+ results[3]*q_long_local*q_long_local
						+ results[4]*q_out_local*q_side_local
						+ results[5]*q_out_local*q_long_local
						+ results[6]*q_side_local*q_long_local )
						, 2 )
                  /sigma_k_prime/sigma_k_prime;
    }

    double chi_sq_per_dof = chi_sq/(data_length - dim);
    cout << "chi_sq = " << chi_sq << endl;
    cout << "Number d.o.f. = " << data_length - dim << endl;
    cout << "chi_sq/d.o.f. = " << chi_sq_per_dof << endl;
	cout << "Goodness-of-fit parameter Q = "
			<< gsl_sf_gamma_inc_Q (0.5*(data_length - dim), 0.5*chi_sq) << endl;

	//============================================
	// compute curvature and covariance matrices
    double ** curvature_mat = new double * [dim];
    double ** covariance_mat = new double * [dim];
    for(int i = 0; i < dim; i++)
    {
        curvature_mat[i] = new double [dim];
        covariance_mat[i] = new double [dim];
        for(int j = 0; j < dim; j++)
		{
            curvature_mat[i][j] = 0.0;
            covariance_mat[i][j] = 0.0;
		}
    }

	//=============================
	for (int i = 0; i < n_qo_bins; i++)
	for (int j = 0; j < n_qs_bins; j++)
	for (int k = 0; k < n_ql_bins; k++)
    {
		int idx = indexer(iKT, iKphi, iKL, i, j, k);

        double q_out_local = 0.5*(qo_pts[i]+qo_pts[i+1]);
        double q_side_local = 0.5*(qs_pts[j]+qs_pts[j+1]);
        double q_long_local = 0.5*(ql_pts[k]+ql_pts[k+1]);

		double correl_local = correlation_function[idx]-1.0;
		double correl_err_local = correlation_function_error[idx];

		if(correl_local < 1.0e-15) continue;

		bool ignore_central_point = true;
		if ( 	ignore_central_point
				and i==(n_qo_bins-1)/2
				and j==(n_qs_bins-1)/2
				and k==(n_ql_bins-1)/2)
			continue;
		//	correl_err_local = 1.0e10;	//ignore central point

        double sigma_k_prime = correl_err_local/correl_local;
        double inv_sigma_k_prime_sq = 1./(sigma_k_prime*sigma_k_prime);

		qweight[0] = - 1.0;
		qweight[1] = q_out_local*q_out_local;
		qweight[2] = q_side_local*q_side_local;
		qweight[3] = q_long_local*q_long_local;
		qweight[4] = q_out_local*q_side_local;
		qweight[5] = q_out_local*q_long_local;
		qweight[6] = q_side_local*q_long_local;

		for(int ij = 0; ij < dim; ij++)
		for(int lm = 0; lm < dim; lm++)
            curvature_mat[ij][lm] += qweight[ij]*qweight[lm]
										* inv_sigma_k_prime_sq;
    }

	//=============================
	// covariance matrix is inverse of curvaure matrix
	gsl_matrix_invert( curvature_mat, covariance_mat, dim );

	//=============================
	// determine errors on fit parameters
	// (i.e., diagonal elements of covariance matrix)
	lambda_Correl_err[indexerK(iKT, iKphi, iKL)] = lambda_Correl[indexerK(iKT, iKphi, iKL)]*(exp(sqrt(covariance_mat[0][0]))-1);
	R2_out_err[indexerK(iKT, iKphi, iKL)] = sqrt(covariance_mat[1][1])*hbarC*hbarC;
	R2_side_err[indexerK(iKT, iKphi, iKL)] = sqrt(covariance_mat[2][2])*hbarC*hbarC;
	R2_long_err[indexerK(iKT, iKphi, iKL)] = sqrt(covariance_mat[3][3])*hbarC*hbarC;
	R2_outside_err[indexerK(iKT, iKphi, iKL)] = sqrt(covariance_mat[4][4])*hbarC*hbarC;
	R2_outlong_err[indexerK(iKT, iKphi, iKL)] = sqrt(covariance_mat[5][5])*hbarC*hbarC;
	R2_sidelong_err[indexerK(iKT, iKphi, iKL)] = sqrt(covariance_mat[6][6])*hbarC*hbarC;

	//=============================
    // clean up
    gsl_matrix_free (T_gsl);
    gsl_matrix_free (T_inverse_gsl);
    gsl_permutation_free (perm);

    delete [] qweight;
    delete [] V;
    for(int i = 0; i < dim; i++)
    {
        delete [] T[i];
        delete [] T_inverse[i];
		delete [] curvature_mat[i];
		delete [] covariance_mat[i];
    }
    delete [] T;
    delete [] T_inverse;
	delete [] curvature_mat;
	delete [] covariance_mat;
    delete [] results;

	return;
}


//End of file
