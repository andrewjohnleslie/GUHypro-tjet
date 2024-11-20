#include "totalT.h"
#include "/home/robert/repos/fast-cpp-csv-parser/csv.h"

int main(){

//    io::CSVReader<10> in("/home/robert/Downloads/exit_totalT(3).csv");
    io::CSVReader<14> in("/home/robert/repos/Cantera_HyPro/HyPro/test/scramjet_combustors/TotalT/BK_EDM_averages_len.csv");
    double x;
    double y;
    double z;
    double rho;
    double t;
    double p;
    double u;
    double Y_O2;
    double Y_N2;
    double Y_H2;
    double Y_H2O;
    double Y_OH;
    double Y_O;
    double Y_H;
    in.read_header(io::ignore_extra_column, "x(m)","y(m)","z(m)", "rho(kg/m3)", "T(K)","p(Pa)","u(m/s)","Y_O2(-)","Y_N2(-)","Y_H2(-)","Y_H2O(-)","Y_OH(-)","Y_O(-)","Y_H(-)");
    std::cout << "x(m), y(m), z(m), rho(kg/m3), T(K), p(Pa), u(m/s), T0(K), p0(K), Y_O2(-), Y_N2(-), Y_H2(-), Y_H2O(-), Y_OH(-), Y_O(-), Y_H(-)" << std::endl;
    while(in.read_row(x,y,z,rho,t,p,u,Y_O2,Y_N2,Y_H2,Y_H2O,Y_OH,Y_O,Y_H)){
//        std::cout << t << std::endl;

        Cantera::compositionMap inletComposition;
        inletComposition.emplace("O2", Y_O2);
        inletComposition.emplace("N2", Y_N2);
        inletComposition.emplace("H2", Y_H2);
        inletComposition.emplace("H2O", Y_H2O);
        inletComposition.emplace("OH", Y_OH);
        inletComposition.emplace("O", Y_O);
        inletComposition.emplace("H", Y_H);

        std::cout << x << ", " << y << ", " << z << ", " << rho << ", " << t << ", " << p << ", " << u;
        find_total_t(t,p,inletComposition,u);
        std::cout << ", " << Y_O2 << ", " << Y_N2 << ", " << Y_H2 << ", " << Y_H2O << ", " << Y_OH << ", " << Y_O << ", " << Y_H << std::endl;;
    }

//    Cantera::compositionMap inletComposition;
////    inletComposition.emplace("O2", 0.0631917247);
////    inletComposition.emplace("N2", 0.7499926192);
////    inletComposition.emplace("H2", 0.0033875664);
////    inletComposition.emplace("H2O", 0.1514206995);
////    inletComposition.emplace("OH", 0.0194478701);
////    inletComposition.emplace("O", 0.0109971729);
////    inletComposition.emplace("H", 0.0015137379);
//    inletComposition.emplace("O2", 0.2342262316);
//    inletComposition.emplace("N2", 0.4780626059);
//    inletComposition.emplace("H2", 0.01385757092);
//    inletComposition.emplace("H2O", 0.2722714673);
//    inletComposition.emplace("OH", 0.001187885184);
//    inletComposition.emplace("O", 0.0001155936797);
//    inletComposition.emplace("H", 0.0002786455224);
//
//
////    find_total_t(2189.723914,127727.4617,inletComposition,2221.165291);
//    find_total_t(1328.723321,98126.12144,inletComposition,1665.411877);
}