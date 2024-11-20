//
// Created by robert on 13/11/17.
//

#include "../lapcat24.h"

// int main() {
int run_burrows(int combustionMode, double combeff, double splitFrac,
				std::string chemistryModel = "gri30")
{
	/* Modes:
   *      Combustion:
   *          0: EffComb
   *          1: EquilComb
   *          2: Split EffComb
   *          3: Split EquilComb
   */

	// Initial Conditions

	// Geometry Modes
	bool par_mixer = false;
	bool injector = true;
	bool premixed = false;

	double T;
	double U;
	double p;
	Cantera::compositionMap inletComposition;
	if (par_mixer == true || injector == true)
	{
		T = 1236.383329;
		U = 1689.308478;
		p = 97154.45662;
		inletComposition.emplace("O2", 0.258);
		inletComposition.emplace("N2", 0.486);
		inletComposition.emplace("H2O", 0.256);
	}
	else if (premixed)
	{
		T = 1193.896487;
		U = 1671.915823;
		p = 96972.47395;
		inletComposition.emplace("H2", 0.01281971967);
		inletComposition.emplace("O2", 0.2546925123);
		inletComposition.emplace("N2", 0.4797696161);
		inletComposition.emplace("H2O", 0.2527181517);
	}

	// Initiate vector for Areas
	std::vector<double> A;

	// Calculate the global equivalence ratio, used by injectionPhi
	double phi = ((1.0 * 0.4661) / (0.258 * 35.889)) / (1.0 / 8.0);

	if (combustionMode == 0)
	{
		// EffComb
		// Post Injector
		A = {0.09376,
			 0.09376,
			 0.09376,
			 0.1048};
	}
	else if (combustionMode == 1)
	{
		// EquilComb
		A = {0.09376,
			 0.09376,
			 0.09376,
			 0.1048};
	}
	else if (combustionMode == 2)
	{
		// Isolator, Split EffComb
		// // Injector
		if (injector == true)
		{
			A = {0.089,
				 0.089,
				 0.09376,
				 0.09376,
				 0.09376 - 0.09376 * splitFrac,
				 0.09376 * splitFrac,
				 0.09376 - 0.09376 * splitFrac,
				 0.09376 * splitFrac,
				 0.09376,
				 0.1048};
		}
		else if (par_mixer == true)
		{
			A = {0.089,
				 0.004,
         0.093,
				 0.09376,
				 0.09376,
				 0.09376 - 0.09376 * splitFrac,
				 0.09376 * splitFrac,
				 0.09376 - 0.09376 * splitFrac,
				 0.09376 * splitFrac,
				 0.09376,
				 0.1048};
		}
		else if (premixed == true)
		{
		}
	}
	else if (combustionMode == 3)
	{
		// Isolator Split EquilComb
		if (injector == true)
		{
			A = {0.09376,
				 0.09376,
				 0.09376,
				 0.09376 - 0.09376 * splitFrac,
				 0.09376 * splitFrac,
				 0.09376 - 0.09376 * splitFrac,
				 0.09376 * splitFrac,
				 0.09376,
				 0.1048};
		}
		else if (par_mixer == true)
		{
			A = {0.08976,
				 0.004,
				 0.09376,
				 0.09376,
				 0.09376 - 0.09376 * splitFrac,
				 0.09376 * splitFrac,
				 0.09376 - 0.09376 * splitFrac,
				 0.09376 * splitFrac,
				 0.09376,
				 0.1048};
		}
		else if (premixed == true)
		{
		}
	}
	else
	{
		std::cout << "Error" << std::endl;
	}

	std::cout << combustionMode << ", " << combeff << ", " << splitFrac
			  << ", ";
	burrows(A, combustionMode, p, T, U, inletComposition, phi, combeff, chemistryModel);
	return 1;
}

void print_columns()
{
	std::cout << "Engine Mode, Combustion Efficiency, Separated "
				 "Ratio, Chemical Model, Density, Velocity, Temperature, "
				 "Pressure, Mass Flow, T0, Y_H2, Y_O2, Y_N2, Y_H2O, Y_OH, Y_O, "
				 "Y_H, Y_NO, Y_NO2"
			  << std::endl;
}

void journal_data()
{
	// Journal Paper Settings
	// Split Complete for SR = 0.1, 0.116, 0.185, 0.4, 1
	std::cout << "############## Split Complete for SR = 0.1, 0.116, 0.185, 0.4 "
				 "##############\n"
			  << std::endl;
	print_columns();
	std::vector<double> combEffs = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35,
									0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7,
									0.75, 0.8, 0.85, 0.9, 0.95, 1.0}; //{0.5};

	for (auto &combEff : combEffs)
		run_burrows(2, combEff, 0.1, "edm");
	for (auto &combEff : combEffs)
		run_burrows(2, combEff, 0.116, "edm");
	for (auto &combEff : combEffs)
		run_burrows(2, combEff, 0.116, "frc");
	for (auto &combEff : combEffs)
		run_burrows(2, combEff, 0.4, "edm");
	for (auto &combEff : combEffs)
		run_burrows(0, combEff, 1, "edm");

	// Split Equilibrium for SR = 0.116, 0.185, 0.1, 0.4, 1
	std::cout
		<< "\n\n############## Split Equilibrium for SR = 0.116, 0.185, 0.1, "
		   "0.4, 1 ##############"
		<< std::endl;
	print_columns();
	run_burrows(3, 0, 0.1, "gri30");
	run_burrows(3, 0, 0.116, "gri30");
	run_burrows(3, 0, 0.116, "edm");
	run_burrows(3, 0, 0.185, "gri30");
	run_burrows(3, 0, 0.185, "frc");
	run_burrows(3, 0, 0.4, "gri30");
	run_burrows(1, 0, 1, "gri30");
}

int main()
{
	/*
  Examples of runtime use

  ** Run through multiple engine modes, combustion efficiencies, and split
  ratios
  std::vector<double> engineModes = {2, 4};
  std::vector<double> combEffs = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
  0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};//{0.5};
  std::vector<double> splits = {0.1, 0.116, 0.185, 0.4};

  for (auto &engineMode : engineModes ) {
      for (auto &combEff : combEffs ){
          for (auto & split : splits)
              run_burrows(engineMode, combEff, split);
      }
  }

  ** Run through multiple engine modes and split ratios
  std::vector<double> engineModes = {3};
  std::vector<double> splits = {0.1, 0.116, 0.185, 0.4};
      for (auto &engineMode : engineModes ) {
      for (auto & split : splits){
          run_burrows(engineMode, 0, split);
      }
  }

  ** Run through multiple engine modes
  std::vector<double> engineModes = {1};
  for (auto &engineMode : engineModes ) {
      run_burrows(engineMode, 0, 0);
  }

  ** Run through multiple engine modes and combustion efficiencies
  std::vector<double> engineModes = {0};
  std::vector<double> combEffs = {1};
      for (auto &engineMode : engineModes ) {
      for (auto & combeff : combEffs){
          run_burrows(engineMode, combeff, 0);
      }
  }
  */
	run_burrows(2, 1, 0.116, "edm");
	//   journal_data();
	return 0;
}