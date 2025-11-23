#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "json.hpp"

#define M_PI 3.14159265358979323846

using json = nlohmann::json;

void calculate_well(
    double h, double fi, double k, double mu, double c, 
    double p_0, double r_w, const double* pressure_profile,
    double* time_days, double* flow_rate_m3_day, double* kf_result,
    int N
) {
    double dt = 86400.0;

    double kf = k / (fi * mu * c);
    *kf_result = kf;
    
    double B = -(r_w * r_w) / (4.0 * kf * dt);
    
    auto expi_approx = [](double x) -> double {
        if (x >= 0) return 0.0;
        
        double result = 0.57721566 + std::log(-x);
        double term = x;
        result += term;
        
        for (int n = 2; n <= 10; n++) {
            term *= x * (n - 1) / (n * n);
            result += term;
        }
        return result;
    };
    
    double Ei0 = -expi_approx(B);
    double A = 4.0 * kf * M_PI * fi * c * h / Ei0;
    
	std::vector<double> time_seconds(N);
    for (int i = 0; i < N; i++) {
        time_seconds[i] = (i + 1) * dt;
    }
    
    std::vector<double> E(N);
    for (int n = 0; n < N; n++) {
        double current_k = n + 1;
        double Ei_n = -expi_approx(B / current_k);
        double Ei_n1 = -expi_approx(B / (current_k + 1));
        E[n] = (Ei_n - Ei_n1) / Ei0;
    }
    
    std::vector<double> flow_rate_m3_sec(N);
    flow_rate_m3_sec[0] = A * (p_0 - pressure_profile[0]);
    
    for (int n = 1; n < N; n++) {
        double sum_term = 0.0;
        for (int i = 0; i < n; i++) {
            sum_term += E[i] * flow_rate_m3_sec[n - 1 - i];
        }
        flow_rate_m3_sec[n] = A * (p_0 - pressure_profile[n]) + sum_term;
    }
    
    for (int i = 0; i < N; i++) {
        time_days[i] = time_seconds[i] / 86400.0;  
        flow_rate_m3_day[i] = flow_rate_m3_sec[i] * 86400.0; 
    }
}

int main() {
    try {
       std::ifstream input("input_params.json");
        if (!input) {
            std::cerr << "n/a" << std::endl;
            return 1;
        }
        
        json params;
        input >> params;
        input.close();

        double h = params["height"];
        double fi = params["porosity"];
        double k = params["permeability"];
        double mu = params["viscosity"];
        double c = params["compressibility"];
        double p_0 = params["initial_pressure"];
        double r_w = params["well_radius"];
        double p_zab = params["bottomhole_pressure"];
        int N = params["points_count"];
        
        std::cout << "ПРОЧИТАНЫ ПАРАМЕТРЫ ИЗ JSON:" << std::endl;
        std::cout << "k=" << k << " м²" << std::endl;
        
        std::vector<double> pressure_profile(N);
        for(int i = 0; i < 50; i++) pressure_profile[i] = p_0;
        for(int i = 50; i < N; i++) pressure_profile[i] = p_zab;
        
        std::vector<double> time_days(N);
        std::vector<double> flow_rate_m3_day(N);
        double kf_result;
        
        calculate_well(h, fi, k, mu, c, p_0, r_w, pressure_profile.data(), 
                       time_days.data(), flow_rate_m3_day.data(), &kf_result, N);
        
		json results;
        results["piezometric_coefficient"] = kf_result;
        results["time_days"] = time_days;
        results["flow_rate_m3_day"] = flow_rate_m3_day;
        
        std::ofstream out("result.json");
        out << results.dump(2);
        out.close();
        
        std::cout << "result.json is ready" << std::endl;
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}