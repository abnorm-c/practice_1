#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#define M_PI 3.14159265358979323846

void calculate_well(
    double h, double fi, double k, double mu, double c, 
    double p_0, double r_w, const double* pressure_profile,
    double* time_result, double* flow_result, double* kf_result,
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
    
    for (int i = 0; i < N; i++) {
        time_result[i] = (i + 1) * dt;
    }
    
    std::vector<double> E(N);
    for (int n = 0; n < N; n++) {
        double current_k = n + 1;
        double Ei_n = -expi_approx(B / current_k);
        double Ei_n1 = -expi_approx(B / (current_k + 1));
        E[n] = (Ei_n - Ei_n1) / Ei0;
    }
    

    flow_result[0] = A * (p_0 - pressure_profile[0]);
    
    for (int n = 1; n < N; n++) {
        double sum_term = 0.0;
        for (int i = 0; i < n; i++) {
            sum_term += E[i] * flow_result[n - 1 - i];
        }
        flow_result[n] = A * (p_0 - pressure_profile[n]) + sum_term;
    }
}

int main() {
    try {
        double h = 20.0, fi = 0.15, k = 1e-15, mu = 1e-3, c = 0.5e-9;
        double p_0 = 250e5, r_w = 0.12;
        int N = 300;
        
        std::vector<double> pressure_profile(N);
        for(int i = 0; i < 50; i++) pressure_profile[i] = p_0;
        for(int i = 50; i < N; i++) pressure_profile[i] = 150e5;
        
        std::vector<double> time_result(N);
        std::vector<double> flow_result(N);
        double kf_result;
        
        calculate_well(h, fi, k, mu, c, p_0, r_w, pressure_profile.data(), 
                       time_result.data(), flow_result.data(), &kf_result, N);
        std::ofstream out("result.txt");
        out << kf_result << std::endl;
        
        for(int i = 0; i < N; i++) {
            out << time_result[i] << " " << flow_result[i] << std::endl;
        }
        out.close();
        
        std::cout << "Done" << std::endl;
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Неизвестная ошибка" << std::endl;
        return 1;
    }
}