#include <pybind11/pybind11.h>
#include <cmath>

double calculate_kf_cpp(double k, double fi, double mu, double c) {
    if (fi == 0 || mu == 0 || c == 0) {
        return 0;
    }
    double kf = k / (fi* mu * c); 
    return kf;
}

double calculate_b_cpp(double r_w, double kf, double dt) {
    if (kf == 0 || dt == 0) {
        return 0;
    } 
    double B = - (r_w * r_w) / (4 * kf * dt);
    return B;
}

double calculate_kf_cpp(double k, double fi, double mu, double c);
double calculate_b_cpp(double r_w, double kf, double dt);

PYBIND11_MODULE(cpp_core, m) {
    m.def("calculate_kf", &calculate_kf_cpp,
        "Расчет коэффициента пьезопроводности",
        pybind11::arg("k"), pybind11::arg("fi"), pybind11::arg("mu"), pybind11::arg("c"));

     m.def("calculate_B", &calculate_b_cpp,
          "Расчет параметра B",
          pybind11::arg("r_w"), pybind11::arg("kf"), 
          pybind11::arg("dt"));

}