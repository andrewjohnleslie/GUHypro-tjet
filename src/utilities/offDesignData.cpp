#include "offDesignData.h"

namespace hypro {


	offDesignData::offDesignData(){}
	offDesignData::~offDesignData(){}


			double offDesignData::findmbar(double p0, double p, double y){

					double mbar;
					double pratio =  p0/p;
					double pchoke = std::pow(((y+1)/2),(y/(y-1)));

					if (pratio < pchoke) {
						mbar = (y/(y-1)) * sqrt(2*((std::pow((1/pratio),(2/y)))   -    (std::pow((1/pratio),((y+1)/y)))));

					}

					else {
						mbar = (y/(y-1)) * sqrt(2*((std::pow((1/pchoke),(2/y)))   -    (std::pow((1/pchoke),((y+1)/y)))));

					}

					return mbar;
			}

			void offDesignData::analysis(plotData &output, outputData &DP, outputData &offDP){

					output.sys = offDP.sys;
					output.throttle = offDP.throttle;

					if (DP.p03 == offDP.p03){

						double mbar = findmbar(DP.p05, DP.p9, DP.y9);
						output.FT = DP.thrust;
						output.pratio = DP.p03 / DP.p02;
						output.Tratio = DP.T04 / DP.T02;
						output.minod = DP.min;
						output.m2A2 = DP.min * (sqrt(DP.T02 * DP.cp2) / DP.p02);
						output.mbar_ = mbar;


					}

					else {

						double mbar = findmbar(DP.p05, DP.p9, DP.y9);
						double mbarOff = findmbar(offDP.p05, offDP.p9, offDP.y9);
						double A4 = (((1+DP.f)*DP.min) * sqrt(DP.cp4*DP.T04))/(DP.p04*mbar);


						output.minod =  ((mbarOff * A4 * offDP.p04)/
										((1+offDP.f)*(sqrt(offDP.cp4*offDP.T04))));

						fmt::print(std::cout, "\n{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12}\n\n","minod", "p04", "A4", "mbaroff ", "f", "cp4", "T04");
						fmt::print(std::cout, "{:<12.2f} {:<12.2f} {:<12.5f} {:<12.5f} {:<12.5f} {:<12.2f} {:<12.2f} \n",DP.min, DP.p04, A4, mbar, DP.f, DP.cp4, DP.T04);
						fmt::print(std::cout, "{:<12.2f} {:<12.2f} {:<12.5f} {:<12.5f} {:<12.5f} {:<12.2f} {:<12.2f} \n\n",output.minod, offDP.p04, A4, mbarOff, offDP.f, offDP.cp4, offDP.T04);

						output.pratio = offDP.p03/offDP.p02;
						output.Tratio = offDP.T04/offDP.T02;
						output.FT = (output.minod*(1+offDP.f)*offDP.V9) - (output.minod*offDP.Ufs);
						output.m2A2 =(output.minod/offDP.p02) * sqrt(offDP.T02 * offDP.cp2);
						output.mbar_ = mbarOff;
						
						
						fmt::print(std::cout, "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} \n","min od", "A4", "mbaroff ", "Tratio", "FT", "m2A2", "gamma_9", "cp2");
						fmt::print(std::cout, "{:<12.7f} {:<12.7f} {:<12.7f} {:<12.5f} {:<12.7f} {:<12.7f} {:<12.7f} {:<12.7f} \n",output.minod, A4, mbarOff, output.Tratio, output.FT, output.m2A2, offDP.y9, offDP.cp2);


					}



			}


			void offDesignData::disp(plotData outputs[], int numbOfPoints , int cppaste){

				std::cout <<  "\n" << "*************************** " << "Post Processing" << " ***************************" << "\n\n";

				if (cppaste == 1){

					for (std::size_t index {0}; index<(7); index++){

						for (std::size_t coutt {0}; coutt<(numbOfPoints); coutt++){

							if (index == 0){
								fmt::print(std::cout, "{:<12.7f}\n", outputs[coutt].throttle);
							}
							else if (index == 1){
								fmt::print(std::cout, "{:<12.7f}\n", outputs[coutt].minod);
							}
							else if (index == 2){
								fmt::print(std::cout, "{:<12.7f}\n", outputs[coutt].Tratio);
							}
							else if (index == 3){
								fmt::print(std::cout, "{:<12.5f}\n", outputs[coutt].FT);
							}
							else if (index == 4){
								fmt::print(std::cout, "{:<12.7f}\n", outputs[coutt].pratio);
							}
							else if (index == 5){
								fmt::print(std::cout, "{:<12.7f}\n", outputs[coutt].m2A2);
							}
							else if (index == 6){
								fmt::print(std::cout, "{:<12.7f}\n", outputs[coutt].mbar_);
							}
						}
					}

				}

				else {

					fmt::print(std::cout, "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12}\n","throttle", "min od", "pratio ", "Tratio", "FT", "m2A2", "mbar");

					for (std::size_t coutt = 0; coutt<(numbOfPoints); coutt++){
						fmt::print(std::cout, "{:<12.6f} {:<12.7f} {:<12.7f} {:<12.7f} {:<12.5f} {:<12.7f} {:<12.5f}\n",outputs[coutt].throttle, outputs[coutt].minod, outputs[coutt].pratio, outputs[coutt].Tratio, outputs[coutt].FT, outputs[coutt].m2A2,  outputs[coutt].mbar_);

					}
				}

			}


	};

