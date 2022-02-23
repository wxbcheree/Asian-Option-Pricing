
//  main.cpp
//  project
//
//  Created by WXB on 2021/5/2.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>

#include <sstream>
#include <stdexcept>
#include "sobol.hpp"
using namespace std;


double T = 1.0;
double K = 100;
double S_0 = 100;
double r = 0.1;
double q = 0.0;
double sigma = 0.2;
double m = 50.0;
int M = 50;
int L = 10;

double sigma_t = sqrt(((2*m+1)*sigma*sigma)/(3*m));
double q_t = q + 0.5*sigma*sigma - 0.5*sigma_t*sigma_t;
double delta = T/m;
double T_t = 0.5*(T+delta);
double se,se_control,se_qmc, se_qmc_control;
int row;
int col;
double **sobol;



double N(const double& z) {
    if (z > 6.0) { return 1.0; }; // this guards against overflow
    if (z < -6.0) { return 0.0; };
    double b1 = 0.31938153;
    double b2 = -0.356563782;
    double b3 = 1.781477937;
    double b4 = -1.821255978;
    double b5 = 1.330274429;
    double p = 0.2316419;
    double c2 = 0.3989423;
    double a=fabs(z);
    double t = 1.0/(1.0+a*p);
    double b = c2*exp((-z)*(z/2.0));
    double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
    n = 1.0-b*n;
    if ( z < 0.0 ) n = 1.0 - n;
    return n;
}

double max(float a, float b) {
    return (b < a) ? a : b;
}

double getgaussrand()
{
    static double V1, V2, S;
    static int phase = 0;
    double X;

    if (phase == 0) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;

            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while (S >= 1 || S == 0);

        X = V1 * sqrt(-2 * log(S) / S);
    }
    else
        X = V2 * sqrt(-2 * log(S) / S);

    phase = 1 - phase;

    return X;
}


double calculate_expected_x(double S,
                                    double K,
                                    double r,
                                    double q,
                                    double sigma,
                                    double time) {

    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+((r-q+0.5*pow(sigma, 2))*time))/(sigma*time_sqrt);
    double d2 = (log(S/K)+((r-q-0.5*pow(sigma, 2))*time))/(sigma*time_sqrt);
    return exp(r*T_t)*(S*exp(-q*time)*N(d1) - K*exp(-r*time)*N(d2));
}


double generate_St(double S_0, double r, double q, double sigma, double T, double rdm){

    return S_0*exp((r-q-0.5*pow(sigma,2))*T + sigma*pow(T,0.5)*rdm);
}

double generate_arithmetic(){
    double sum = 0;
    double temp_rdm;
    temp_rdm = getgaussrand();

    double S_t = generate_St(S_0, r, q, sigma, delta, temp_rdm);
    for (int i = 0; i < m; i++) {
        temp_rdm = getgaussrand();
        sum += S_t;

        S_t = generate_St(S_t, r, q, sigma, delta, temp_rdm);
    }
    return sum/m;

}


double option_price_call_arithmetic_direct(long sample_size){
    double S_t_bar = generate_arithmetic();
    double X_bar = max(S_t_bar - K, 0);
    double Y_bar = pow(X_bar, 2);
    double call_price;

    for (int i = 2; i<= sample_size; i++) {
        S_t_bar = generate_arithmetic();
        if (S_t_bar > K) {
        call_price = max(S_t_bar - K, 0);
        X_bar = (1.0-1.0/i)*X_bar + (1.0/i)*call_price;
        Y_bar= (1.0-1.0/i)*Y_bar + (1.0/i)*pow(call_price,2);
        }
        else{
        X_bar = (1.0-1.0/i)*X_bar ;
        Y_bar= (1.0-1.0/i)*Y_bar ;
        }
    }
    X_bar = exp(-r*T)*X_bar;
    Y_bar = exp(-2*r*T)*Y_bar;

    se = sqrt((1.0/(sample_size-1))*(Y_bar-pow(X_bar,2)));

    return X_bar;
}


double *generate_geometric_arithmetic(){
    double sum = 0;
    double product = 1;
    double temp_rdm = getgaussrand();
    double S_t = generate_St(S_0, r, q, sigma, delta, temp_rdm);
    for (int i = 0; i < m; i++) {
        temp_rdm = getgaussrand();
        
        product *= S_t;
        //cout << product<<endl;
        sum += S_t;
        S_t = generate_St(S_t, r, q, sigma, delta, temp_rdm);
    }
    product = pow(product, 1/m);
    sum = sum/m;
    static double S_t_average[2];
    S_t_average[0] = product;
    S_t_average[1] = sum;
    //cout << S_t_average[0]<<endl;
    //cout << S_t_average[1]<<endl;

    return S_t_average;
}

double calculate_bj(){
    // Y: Arithmetic
    // X: geometric
    int no_of_sample = 1000;
    double X[1000];
    double Y[1000];
    double var_x = 0;
    double var_y = 0;
    double cov = 0;

    double *S_t_average;
    for (int i = 0; i<no_of_sample; i++) {
        S_t_average = generate_geometric_arithmetic();


        X[i] = max(S_t_average[0]-K,0);
        Y[i] = max(S_t_average[1]-K,0);
        //cout<<X[i]<<endl;

    }



    double sum_x = 0;
    double sum_y = 0;


    for(int i = 0; i < no_of_sample; i++){

        sum_x += X[i];
        sum_y += Y[i];
    }
    double mean_x = sum_x / no_of_sample;
    double mean_y = sum_y / no_of_sample;


    for(int i = 0; i < no_of_sample; i++){

        cov += (X[i] - mean_x) * (Y[i] - mean_y);
        var_x += (X[i] - mean_x) * (X[i] - mean_x);
        var_y += (Y[i] - mean_y) * (Y[i] - mean_y);
    }
    return cov/var_x;
    //return cov/(sqrt(var_x*var_y));

}


double option_price_call_arithmetic_control(long sample_size){
    double b_j = calculate_bj();
    double expected_x = calculate_expected_x(S_0, K, r, q_t, sigma_t, T_t);
    double *S_t_average;
    S_t_average = generate_geometric_arithmetic();
    double Y_i = max(S_t_average[1]-K,0);
    double X_i = max(S_t_average[0]-K,0);
    double Y_i_b = Y_i + b_j*(expected_x - X_i);
    double Y_i_b_2 = pow(Y_i_b, 2);
    double temp;
    for (int i = 2; i <= sample_size; i++) {
        S_t_average = generate_geometric_arithmetic();
        double Y_i = max(S_t_average[1]-K,0);
        double X_i = max(S_t_average[0]-K,0);
        temp = Y_i + b_j*(expected_x - X_i);
        Y_i_b = (1.0-1.0/i)*Y_i_b + (1.0/i)*temp;
        Y_i_b_2= (1.0-1.0/i)*Y_i_b_2 + (1.0/i)*pow(temp,2);

    }
    Y_i_b = exp(-r*T)*Y_i_b;
    Y_i_b_2 = exp(-2*r*T)*Y_i_b_2;
    se_control = sqrt((1.0/(sample_size-1))*(Y_i_b_2-pow(Y_i_b,2)));
    return Y_i_b;
}


double RationalApproximation(double t)
{
    // https://www.johndcook.com/blog/normal_cdf_inverse/
    
    
    double c[] = {2.515517, 0.802853, 0.010328};
    double d[] = {1.432788, 0.189269, 0.001308};
    return t - ((c[2]*t + c[1])*t + c[0]) /
               (((d[2]*t + d[1])*t + d[0])*t + 1.0);
}

double NormalCDFInverse(double p)
{
    if (p < 0.5)
    {
        return -RationalApproximation( sqrt(-2.0*log(p)) );
    }
    else
    {
        return RationalApproximation( sqrt(-2.0*log(1-p)) );
    }
}


double quasi_Monte_Carlo_simulation_with_sobol(int sample_size) {
    double x_bar_average = 0;
    double y_bar_average = 0;
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator;// (seed);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    for (int i = 1; i <= L; i++) {
        double* sobol_seq = i8_sobol_generate(50, sample_size, 10000);
        double uniform[50];
        for (int p = 0; p < M; p++) {
            uniform[p] = distribution(generator);
        }
        double sum = 0;
        double S_t = S_0;
        double mod = 0;
        for (int o = 0; o < M; o++) {
            mod = sobol_seq[o] + uniform[o] - floor(sobol_seq[o] + uniform[o]);
            S_t = S_t * exp((r - q - 0.5 * pow(sigma, 2)) * delta + sigma * sqrt(delta) * NormalCDFInverse(mod));
            sum += S_t;
        }
        double average_St = sum / m;
        double value = max(average_St - K, 0);
        double x_bar = value;
        for (int n = 2; n <= sample_size; n++) {
            sum = 0;
            S_t = S_0;
            for (int o = 0; o < M; o++) {
                mod = sobol_seq[(n-1)*50 + o] + uniform[o] - floor(sobol_seq[(n - 1) * 50 + o]
                    + uniform[o]);
                S_t = S_t * exp((r - q - 0.5 * pow(sigma, 2)) * delta + sigma * sqrt(delta) * NormalCDFInverse(mod));
                sum += S_t;
            }
            average_St = sum / m;
            value = max(average_St - K, 0);
            x_bar = (1.0 - 1.0 / n) * x_bar + (1.0 / n) * value;
        }
        x_bar = exp(-r * T) * x_bar;
        x_bar_average = (1 - 1.0 / i) * x_bar_average + (1.0 / i) * x_bar;
        y_bar_average = (1 - 1.0 / i) * y_bar_average + (1.0 / i) * x_bar*x_bar;
    }
    se_qmc = sqrt((1.0 / 9.0) * (y_bar_average- pow(x_bar_average, 2)));
    return x_bar_average;
}



double calculate_bj_rqmc_control(){
    // Y: Arithmetic
    // X: geometric
    int no_of_sample = 1000;
    double* sobol_bj = i8_sobol_generate(50, no_of_sample, 0);

    double X[1000];
    double Y[1000];
    double var_x = 0;
    double var_y = 0;
    double cov = 0;
    double S_t;
    double sum, product;
    double ave_sum, ave_product;


    for (int i = 0; i < no_of_sample; i++){
        sum = 0;
        product = 1;
        for (int j = 0; j < M; j++){
            S_t = S_0;
            S_t = S_t * exp((r - q - 0.5 * pow(sigma, 2)) * delta + sigma * sqrt(delta) * NormalCDFInverse(sobol_bj[i*M+j]));
            sum += S_t;
            product *= S_t;
        }
        ave_sum = sum / 50;
        ave_product = pow(product, 1 / m);
        Y[i] = max(ave_sum - K, 0);
        X[i] = max(ave_product - K, 0);

    }


    double sum_x = 0;
    double sum_y = 0;


    for(int i = 0; i < no_of_sample; i++){

        sum_x += X[i];
        sum_y += Y[i];
    }
    double mean_x = sum_x / no_of_sample;
    double mean_y = sum_y / no_of_sample;


    for(int i = 0; i < no_of_sample; i++){

        cov += (X[i] - mean_x) * (Y[i] - mean_y);
        var_x += (X[i] - mean_x) * (X[i] - mean_x);
        var_y += (Y[i] - mean_y) * (Y[i] - mean_y);
    }
    return cov/var_x;
    //return cov/(sqrt(var_x*var_y));

}






double RMCQ_with_control(int sample_size){
    double average_yib = 0;
    double average_yib2 = 0;
    double b_j = calculate_bj_rqmc_control();
    std::default_random_engine generator;// (seed);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double expected_x = calculate_expected_x(S_0, K, r, q_t, sigma_t, T_t);

    for (int i = 1; i <= L; i++) {

        double* sobol_seq = i8_sobol_generate(m, sample_size, 100);

        double uniform[50];
        for (int p = 0; p < M; p++) {
            uniform[p] = distribution(generator);
        }

        double sum = 0;
        double product = 1;
        double S_t = S_0;
        double mod = 0;
        for (int o = 0; o < M; o++) {
            mod = sobol_seq[o] + uniform[o] - floor(sobol_seq[o] + uniform[o]);
            S_t = S_t * exp((r - q - 0.5 * pow(sigma, 2)) * delta + sigma * sqrt(delta) * NormalCDFInverse(mod));
            sum += S_t;
            product *= S_t;
        }
        double average_St_sum = sum / m;
        double average_St_product = pow(product, 1.0/m);

        double value_sum = max(average_St_sum - K, 0);
        double value_product = max(average_St_product - K, 0);
        double Y_i_b = value_sum + b_j*(expected_x - value_product);
        //double Y_i_b_2 = pow(Y_i_b, 2);


        for (int n = 2; n <= sample_size; n++) {
            sum = 0;
            product = 1;
            S_t = S_0;
            double temp;
            for (int o = 0; o < M; o++) {
                mod = sobol_seq[M*(n-1)+o] + uniform[o] - floor(sobol_seq[M * (n - 1) + o] + uniform[o]);
                S_t = S_t * exp((r - q - 0.5 * pow(sigma, 2)) * delta + sigma * sqrt(delta) * NormalCDFInverse(mod));
                sum += S_t;
                product *= S_t;
            }
            average_St_sum = sum/m;
            average_St_product = pow(product, 1.0/m);

            value_sum = max(average_St_sum - K, 0);
            value_product = max(average_St_product - K, 0);
            temp = value_sum + b_j*(expected_x - value_product);
            Y_i_b = (1.0-1.0/n)*Y_i_b + (1.0/n)*temp;


        }
        Y_i_b = exp(-r*T)*Y_i_b;
        average_yib = (1 - 1.0 / i) * average_yib + (1.0 / i) * Y_i_b;
        average_yib2 = (1 - 1.0 / i) * average_yib2 + (1.0 / i) * Y_i_b*Y_i_b;


    }
    se_qmc_control = sqrt((1.0 / 9.0) * (average_yib2 - pow(average_yib, 2)));

    return average_yib;
}



int main(int argc, const char * argv[]){
        double start, end, time;
        double start_control, end_control, time_control;
        double start_qmc, end_qmc, time_qmc;
        double start_qmc_control, end_qmc_control, time_qmc_control;



        start = clock();
        std::cout << "Estimated price(Monte_Carlo direct approach): " <<option_price_call_arithmetic_direct(100000)<<std::endl;
        end = end = clock();
        time = double(end - start) / CLOCKS_PER_SEC;
        std::cout << "The standard error " << se<< std::endl;
        std::cout << "The computational time in seconds = " << time << "s" << std::endl;
        std::cout << "The efficiency = " << time*se*se  << std::endl;

        start_control = clock();
        std::cout << "Estimated price(with control variate): " <<option_price_call_arithmetic_control(10000000)<<std::endl;
        end_control = clock();
        time_control = double(end_control - start_control) / CLOCKS_PER_SEC;
        std::cout << "The standard error " << se_control << std::endl;
        std::cout << "The computational time in seconds = " << time_control << "s" << std::endl;
        std::cout << "The efficiency = " << time_control*se_control*se_control  << std::endl;

        start_qmc = clock();
        std::cout << "Estimated price(with qmc): " <<quasi_Monte_Carlo_simulation_with_sobol(1000000)<<std::endl;
        end_qmc = clock();
        time_qmc = double(end_qmc - start_qmc) / CLOCKS_PER_SEC;
        std::cout << "The standard error " << se_qmc << std::endl;
        std::cout << "The computational time in seconds = " << time_qmc << "s" << std::endl;
        std::cout << "The efficiency = " << time_qmc * se_qmc * se_qmc << std::endl;

        start_qmc_control = clock();
        cout << RMCQ_with_control(1000000) << endl;
        end_qmc_control = clock();
        time_qmc_control = double(end_qmc_control - start_qmc_control) / CLOCKS_PER_SEC;
        std::cout << "The standard error " << se_qmc_control << std::endl;
        std::cout << "The computational time in seconds = " << time_qmc_control << "s" << std::endl;
        std::cout << "The efficiency = " << time_qmc_control * se_qmc_control * se_qmc_control << std::endl;




    }

