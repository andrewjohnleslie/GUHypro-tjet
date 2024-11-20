//
// Created by mark on 25/11/19.
//

#ifndef HYPRO_OFFDESIGNDATA_H
#define HYPRO_OFFDESIGNDATA_H
#include "systemModule.h"

namespace hypro {
    struct plotData {

        double m2A2;
        double FT;
        double Tratio;
        double pratio;
        double minod;
        double throttle;
        double mbar_;
		systemModule *sys;

    };

    struct outputData {
      	double thrust;
		vector<double> Mfr;
		Collection::OutType out;
		double dynamic_pressure;
		vector<double> emissions;
		double total_flow;
		systemModule *sys;
		double kH;
		double Aratio;
		double T04;

		double kH_HP;
		double Aratio_HP;
		double kH_LP;
		double Aratio_LP;
		double T045;
		double T02;
		double p9;
		double p05;
		double p04;
		double p03;
		double p02;
		double V9;
		double Ufs;
		double f;
		double cp2;
		double cp4;
		double y5;
		double y9;

		double throttle {0};
		double Mach {0};

		double min;
		double A4;
		double mbar9;
		double m2A2;
    };

    class offDesignData {

    public:

		offDesignData();
		virtual ~offDesignData();

		void disp(plotData outputs[], int numbOfPoints, int cppaste);
		void analysis(plotData &output, outputData &DP, outputData &offDP);
		double findmbar(double p0, double p, double y);

    };
}
#endif //HYPRO_OFFDESIGNDATA_H
