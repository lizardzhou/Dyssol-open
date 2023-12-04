/* Copyright (c) 2020, Dyssol Development Team. All rights reserved. This file is part of Dyssol. See LICENSE file for license information. */

#define DLL_EXPORT
#include "VibratedFluidizedBedDryer.h"
#include "Helpers.h"
#include <sstream>

extern "C" DECLDIR CBaseUnit* DYSSOL_CREATE_MODEL_FUN()
{
	return new CVibratedFluidizedBedDryer();
}

//////////////////////////////////////////////////////////////////////////
/// Unit
void CVibratedFluidizedBedDryer::CreateBasicInfo()
{
	/// Basic unit's info ///
	SetUnitName  ("Vibrated Fluidized Bed Dryer (steady-state)");
	SetAuthorName("Buchholz (based on Zhengyu Lu's Master's Thesis)");
	SetUniqueID  ("359EE781CDA64DB8ACFFD92D294E7225");
}

void CVibratedFluidizedBedDryer::CreateStructure()
{
	/// Add ports ///
	AddPort("InGas"      , EUnitPort::INPUT);
	AddPort("InParticle" , EUnitPort::INPUT);
	AddPort("OutGas"     , EUnitPort::OUTPUT);
	AddPort("OutParticle", EUnitPort::OUTPUT);

	// Operational parameters
	AddStringParameter("Operational parameters", "", "");
	AddConstRealParameter("eps_0", 0, "-", "Initial bed porosity.", 0, 1);
	AddConstRealParameter("H_0", 0, "m", "Initial bed height (fixed bed).", 0, 1000);
	AddConstRealParameter("Vibration_f", 0, "1/s", "Frequency of vibration.", 0, 1000);
	AddConstRealParameter("Vibration_A", 0, "m", "Amplitude of vibration.", 0, 1000);
	AddConstRealParameter("T_Wall", 293.15, "K", "Temperature of the inner wall of the dryer.", 273.15, 1000);

	// Geometric parameters
	AddStringParameter("Geometric parameters", "", "");
	AddComboParameter("Shape", E2I(EShape::CYLINDRICAL), E2I({ EShape::CYLINDRICAL, EShape::RECTANGULAR }), { "Cylindrical", "Rectangular" }, "Geometry of the fluidized bed");
	AddConstRealParameter("Diameter", 0.18, "m", "Diameter of the fluidized bed", 0.0);
	AddConstRealParameter("Length"  , 0.18, "m", "Length of the fluidized bed"  , 0.0);
	AddConstRealParameter("Width"   , 0.18, "m", "Width of the fluidized bed"   , 0.0);
	AddParametersToGroup("Shape", "Cylindrical", { "Diameter" });
	AddParametersToGroup("Shape", "Rectangular", { "Length", "Width"});
	AddConstUIntParameter("nor", 100, "-", "Number of orifices of the distributor plate. \n nor = 0 corresponds to porous plate distributors.", 0, 1000000);

	// Material parameters
	AddStringParameter("Material parameters", "", "");
	AddConstRealParameter("k_dc", 3.5, "-", "k for normalized drying curve. \nthe normalized drying curve (nu) of the material is represented by: \n nu = k* eta / (1. + eta*(k - 1.)), \n with eta representing the particle moiture content.", -10000, 10000);
	AddConstRealParameter("Xcr", 0.2, "kg/kg", "Critical water content of the particles (transition between 1st and 2nd drying period).", 0, 100);
	AddConstRealParameter("phi_p", 0.8, "-", "Relative humidity in equilibrium state at particle surface. \nThis assumption is used to calculated the mass transfer between particles and supsension gas", 0, 1);
	AddComboParameter("Geldart", E2I(EGeldart::D), E2I({ EGeldart::A, EGeldart::B, EGeldart::C, EGeldart::D, EGeldart::AUTO }), { "A", "B", "C", "D", "Auto" }, "Geldart class of particles");
	AddConstRealParameter("u_mf", 0, "m/s", "Manual input value for minimum fluidization velocity, especially important for vibration and/or cohesive particles. \nIf 0, the u_mf is calculated, using the correlation by Wen&Yu.", 0, 1000);

	// Compound parameters
	AddStringParameter("Compound names", "", "");
	AddCompoundParameter("Particles (s)", "Name of the particle compound in the material database.");
	AddCompoundParameter("H2O (l)", "Name of the H2O (l) compound in the material database.");
	AddCompoundParameter("H2O (g)", "Name of the H2O (g) compound in the material database.");
	AddCompoundParameter("Air (g)", "Name of the air compound in the material database.");

	// Solver parameters
	AddStringParameter("Solver parameters", "", "");
	AddConstUIntParameter("MaxIterations", 2000, "-", "Maximum number of allowed solver iterations.", 1, 100000);
	AddConstUIntParameter("AndersonAcceleration", 0, "-", "Fixed point solver parameter for adjusting Anderson Acceleration, to improve solver stability.", 0, 10000);
	AddConstUIntParameter("HeightClasses", 15, "-", "Number of height classes.", 1, 10000);
	AddConstUIntParameter("RTClasses", 200, "-", "Number of residence time classes.", 1, 10000);
	AddCheckBoxParameter("Switch_HT", true, "Switch for heat transfer. Heat transfer will NOT be calculated for input value of '0'. Otherwise heat transfer is considered in the model. \nThis option may be used to only calculate the hydrodynamic model.");
	AddCheckBoxParameter("Switch_MT", true, "Switch for mass transfer. Mass transfer will NOT be calculated for input value of '0'. Otherwise mass transfer is considered in the model. \nThis option may be used to only calculate the hydrodynamic model.");

	// Paths for writing out csv data
	AddStringParameter("Paths", "", "");
	AddStringParameter("csv_path", "C:\\", "Directory where model results are exported to (as .csv files).");

	// Model settings
	AddStringParameter("Model setting parameters", "", "");
	AddConstRealParameter("C_Qpp", 0, "-", "Switch for Qpp (heat transfer between particles of different classes). If zero, Qpp will not be considered.", 0, 1);
	AddConstRealParameter("C_Qw", 0, "-", "Switch for Qw (heat transfer from the bed to the wall). If zero, Qw will not be considered.", 0, 1);
	AddConstRealParameter("C_Ntaumax", 0.99, "-", "Fraction of particles, considered inside the residence time classes. The default value 0.99 means that 99 % of the particles are evaluated.", 0, 1);
	AddComboParameter("C_ApAbed", E2I(EArea::SAUTER_DIAMETER), E2I({ EArea::SAUTER_DIAMETER, EArea::MEAN_AND_FRACTIONS, EArea::MEAN }), { "Sauter diameter", "Mean and fractions", "Mean" }, "Setting parameter for the calculation of ApAbed (surface area of all particles in the bed). \n 0: Use Sauter diameter. 1: Use mean particle diameters and mass fractions. \n 2: Only use mean particle diameters");
	AddConstRealParameter("C_Ktis", 1, "-", "Setting parameter for the calculation of residence time distribution (Number of tanks in series model) \n 1: ideally mixed (CSTR) \n >1:  tanks in series model (no integer required)", 1, 100);

	// Adjustment parameters / Assumptions
	AddStringParameter("Adjustment parameters", "", "");
	AddConstRealParameter("H_NTU", 0.05, "m", "Height at which NTU between bubble and suspension phase is assumed to = 1 \n Default is 0.05 m acc. to Groenewold & Tsotsas 1999. \nThis parameter is used to calculate the mass transfer coefficient between bubbles and supsension", 0, 1000);
	AddConstRealParameter("OverwriteBubbleVolFraction", 2, "-", "Input value for bubble volume fraction ,epsilon_B, (0 < epsilon_B < 1). \nIf value >1, bubble volume fraction will be calculated according to Lehmann et al., Powder Technology 2019", 0, 2);

	m_FDNLModel.SetUserData(this);
	m_TDNLModel.SetUserData(this);
}

void CVibratedFluidizedBedDryer::Initialize(double _time)
{
	//// CHECK CALCULATION SETUP ////
	//Check if phases are defined
	if (!IsPhaseDefined(EPhase::VAPOR))
		RaiseError("Vapor phase has not been defined.");
	if (!IsPhaseDefined(EPhase::SOLID))
		RaiseError("Solid phase has not been defined.");
	//Check if distributions are defined
	if (!IsDistributionDefined(EDistrTypes::DISTR_SIZE))
		RaiseError("Size distribution has not been defined.");
	if (!IsDistributionDefined(EDistrTypes::DISTR_MOISTURE))
		RaiseError("Moisture distribution has not been defined.");

	// Stream pointers
	m_inGas       = GetPortStream("InGas");
	m_inParticle  = GetPortStream("InParticle");
	m_outGas      = GetPortStream("OutGas");
	m_outParticle = GetPortStream("OutParticle");
	SPorts inPorts(m_inGas, m_inParticle);
	// Discretization parameters
	size_t heightGridPoints   = GetConstUIntParameterValue("HeightClasses") + 1;
	size_t RTGridPoints       = GetConstUIntParameterValue("RTClasses") + 1;
	size_t sizeClasses        = GetClassesNumber(EDistrTypes::DISTR_SIZE);
	size_t moistureClasses    = GetClassesNumber(EDistrTypes::DISTR_MOISTURE);
	size_t FDHeightGridPoints = GetConstUIntParameterValue("HeightClasses") * 10 + 1;

	if (sizeClasses < 1)
		RaiseError("No classes for size distribution defined.");
	if (moistureClasses < 1)
		RaiseError("No classes for moisture distribution defined.");
	SDiscretizationParameters discrParams(heightGridPoints, RTGridPoints, sizeClasses, moistureClasses, FDHeightGridPoints);

	// Compound parameters
	std::string compKeyParticles = GetCompoundParameterValue("Particles (s)");
	std::string compKeyH2Ol      = GetCompoundParameterValue("H2O (l)");
	std::string compKeyH2Og      = GetCompoundParameterValue("H2O (g)");
	std::string compKeyAir       = GetCompoundParameterValue("Air (g)");
	SCompoundKeys compKeys(compKeyParticles, compKeyH2Ol, compKeyH2Og, compKeyAir);
	m_compKeys = compKeys;
	// Calculation parameters
	const size_t maxIter     = GetConstUIntParameterValue("MaxIterations");
	const size_t andersonAcc = GetConstUIntParameterValue("AndersonAcceleration");

	// Material parameters
	const double k_dc     = GetConstRealParameterValue("k_dc");
	const double X_cr     = GetConstRealParameterValue("Xcr");
	const double phi_eq_p = GetConstRealParameterValue("phi_p");
	const auto geldart = static_cast<EGeldart>(GetComboParameterValue("Geldart"));
	if (geldart != EGeldart::A && geldart != EGeldart::B && geldart != EGeldart::C && geldart != EGeldart::D)
		RaiseWarning("Unknown Geldart class. Determine Geldart class automatically...");
	double T_0 = 273.15;
	double hv_H2O = GetCompoundProperty(compKeyH2Ol, ECompoundConstProperties::HEAT_OF_VAPORIZATION_AT_NORMAL_BOILING_POINT) / GetCompoundProperty(compKeyH2Ol, ECompoundConstProperties::MOLAR_MASS);
	SUnit unit = SUnit(T_0, hv_H2O, phi_eq_p, X_cr, k_dc);

	// Geometric parameters
	SGeometry geometry = ReadGeometryParameter();

	// Bed parameters
	double eps_0    = GetConstRealParameterValue("eps_0");
	double H_fb_0   = GetConstRealParameterValue("H_0");
	double u_mf_man = GetConstRealParameterValue("u_mf");
	double f        = GetConstRealParameterValue("Vibration_f");
	double Af       = GetConstRealParameterValue("Vibration_A");
	double T_Wall   = GetConstRealParameterValue("T_Wall");
	double rho_s    = GetCompoundProperty(compKeyParticles, ECompoundTPProperties::DENSITY, T_0, 1e5);
	SBed bedParams(eps_0, H_fb_0, u_mf_man, rho_s, f, Af, T_Wall, geometry);

	SConstantParameters constantParameters(unit, geometry, bedParams);

	// Paths
	const std::string csvPath = GetStringParameterValue("csv_path");

	// Setting parameters
	const double C_Qpp = GetConstRealParameterValue("C_Qpp");
	if (C_Qpp == 0.0)
		ShowInfo("Unit parameter C_Qpp = 0: No heat transfer between particles is considered.");
	const double C_Qw = GetConstRealParameterValue("C_Qw");
	if (C_Qw == 0.0)
		ShowInfo("Unit parameter C_Qw = 0: No heat transfer to wall considered.");
	double C_Ntaumax = GetConstRealParameterValue("C_Ntaumax");
	if (C_Ntaumax < 0.5)
	{
		RaiseWarning("Unit parameter Ntaumax smaller than 0.5 leads to poor results. Using default value 0.99 instead.");
		C_Ntaumax = 0.99;
	}
	const auto C_ApAbed = static_cast<EArea>(GetComboParameterValue("C_ApAbed"));
	double C_Ktis = GetConstRealParameterValue("C_Ktis");
	bool isHT = GetCheckboxParameterValue("Switch_HT");
	bool isMT = GetCheckboxParameterValue("Switch_MT");

	double C_H_NTU = GetConstRealParameterValue("H_NTU");
	double C_BubbleVolFraction = GetConstRealParameterValue("OverwriteBubbleVolFraction");

	SSettings settings = SSettings(csvPath, C_Qpp, C_Qw, C_Ntaumax, C_ApAbed, C_Ktis, isHT, isMT, C_H_NTU, C_BubbleVolFraction);

	/// Simulation setup infos ///
	ShowInfo("\nFluidized bed model - simulation settings");
	ShowInfo("Heat transfer: " + BoolToOnOff(isHT));
	ShowInfo("Mass transfer: " + BoolToOnOff(isMT));
	ShowInfo("\n");

	/// Add state variables ///
	AddStateVariable("Bed mass [kg]", 0);
	AddStateVariable("u_mf [m/s]", 0);
	AddStateVariable("Geldart type", 0);
	AddStateVariable("Deviation H2O-Balance [%]", 0);
	AddStateVariable("Deviation Enthalpy-Balance [%]", 0);
	AddStateVariable("Tau_max [s]", 0);
	AddStateVariable("T_Particle [K]", 0);
	AddStateVariable("X_Particle [kg/kg]", 0);
	AddStateVariable("nu_mean [-]", 0);
	AddStateVariable("Q_Wall [W]", 0);
	AddStateVariable("alpha_SW_mean [W/(m^2*K)]", 0);
	AddStateVariable("alpha_PW_mean [W/(m^2*K)]", 0);

	/// Add plots ///
	// Fluid dynamics
	AddPlot("Bed height", "Time [s]", "Height [m]");
	AddCurveOnPlot("Bed height", "Fluidized bed reactor");
	AddPlot("d_bubble", "Height [m]", "Diameter [m]", "Time [s]");
	// Gas
	AddPlot("Y_g", "Height [m]", "Y [kg/kg]");
	AddPlot("phi_g", "Height [m]", "phi [%]");
	AddPlot("h_g", "Height [m]", "h [J/kg]");
	AddPlot("T_g", "Height [m]", "T [degC]");
	// Particles
	AddPlot("Number density distribution", "Residence time [s]", "Density [1/s]", "Time [s]");
	std::vector<double> vX_p_in = GetClassesMeans(DISTR_MOISTURE);
	for (size_t i = 0; i < moistureClasses; ++i)
	{
		AddPlot("X_p (X_in- " + std::to_string(vX_p_in[i]) + ")", "Residence time [s]", "X [kg/kg]", "Diameter [m]");
		AddPlot("h_p (X_in- " + std::to_string(vX_p_in[i]) + ")", "Residence time [s]", "h [J/kg]", "Diameter [m]");
		AddPlot("T_p (X_in- " + std::to_string(vX_p_in[i]) + ")", "Residence time [s]", "T [degC]", "Diameter [m]");

		AddPlot("Q_3_X (X_in- " + std::to_string(vX_p_in[i]) + ")", "X [kg/kg]", "Q3 [%]", "Diameter [m]");
	}
	// Heat transfer coefficient
	AddPlot("Heat transfer", "Height [m]", "Heat transfer coefficient [W/(m^2 K)]", "Diameter [m]");
	// Ratio surface of particles to cross-sectional area of the bed
	AddPlot("ApAbed", "Diameter [m]", "ApAbed [-]", "Time [s]");
	// Gas distribution coefficient
	AddPlot("Bubble volume fraction nu", "Height [m]", "nu [-]", "Time [s]");
	//////////////////////////////////////////////////////////////////////////////////////////////////
	//==============================
	// Set NLModel 1
	//==============================
	// Part 1: Calculating local bubble size
	m_FDNLModel.Setup(this, compKeys, constantParameters, inPorts, discrParams, settings);

	STDParameters sTDParams = m_FDNLModel.CalculateTDParameters(_time, geldart);
	// Initialize model parameters
	m_FDNLModel.Initialize(_time, sTDParams);

	m_FDNLModel.AddVariables();
	m_FDNLSolver.SetFixedPointSolverParameters(0, 1);
	m_FDNLSolver.SetStrategy(ENLSolverStrategy::Fixedpoint);

	//Set model to solver

	if (!m_FDNLSolver.SetModel(&m_FDNLModel))
		RaiseError(m_FDNLSolver.GetError());

	//////////////////////////////////////////////////////////////////////////////////////////////////
	//==============================
	// Set NLModel 2
	//==============================
	// Part 2: Calculating parameters of heat and mass transfer
	m_TDNLModel.Setup(this, compKeys, constantParameters, inPorts, discrParams, settings);
	m_TDNLModel.AddVariables(_time);

	m_TDNLSolver.SetStrategy(ENLSolverStrategy::Fixedpoint);
	m_TDNLSolver.SetMaxIter(maxIter);
	m_TDNLSolver.SetFixedPointSolverParameters(andersonAcc, 1);
	// Set the maximum iteration times of the NLSolver
	//m_TDNLSolver.SetSolverMaxIter(nMaxIter);
	if (!m_TDNLSolver.SetModel(&m_TDNLModel))
		RaiseError(m_TDNLSolver.GetError());
	//////////////////////////////////////////////////////////////////////////////////////////////////
}

SGeometry CVibratedFluidizedBedDryer::ReadGeometryParameter()
{
	const auto form = static_cast<EShape>(GetComboParameterValue("Shape"));
	const size_t numOr = GetConstUIntParameterValue("nor");

	switch (form)
	{
	case EShape::CYLINDRICAL:
	{
		const double dd = GetConstRealParameterValue("Diameter");
		return SGeometry{ dd, numOr };
	}
	case EShape::RECTANGULAR:
	{
		const double dl = GetConstRealParameterValue("Length");
		const double db = GetConstRealParameterValue("Width");
		return SGeometry{ dl, db, numOr };
	}
	}

	RaiseError("Invalid input for geometrical parameter " + StringFunctions::Quote("Shape") + ".");
	return {};
}

void CVibratedFluidizedBedDryer::Simulate(double _time)
{
	////////////////////////////////////////////////////////////////////////////////////
	//                  Calculate minimum fluidization point						  //
	////////////////////////////////////////////////////////////////////////////////////
	const auto geldart = static_cast<EGeldart>(GetComboParameterValue("Geldart"));
	STDParameters TDParams = m_FDNLModel.CalculateTDParameters(_time, geldart);

	////////////////////////////////////////////////////////////////////////////////////
	//                  Solving fluid dynamics model                                  //
	////////////////////////////////////////////////////////////////////////////////////
	// Initialize model parameters
	m_FDNLModel.Initialize(_time, TDParams);
	// Run model calculations
	if (!m_FDNLSolver.Calculate(_time))
		RaiseError(m_FDNLSolver.GetError());

	// Results of FD model
	SResultsNLModelFD ResultsFDModel = m_FDNLModel.GetResultsNLModelFD();

	/// Plotting ///
	// Fluid dynamics
	if (ResultsFDModel.valid)
	{
		AddPointOnCurve("Bed height", "Fluidized bed reactor", _time, ResultsFDModel.h_fb);
		AddCurveOnPlot("d_bubble", _time);
		AddPointsOnCurve("d_bubble", _time, ResultsFDModel.h_z, ResultsFDModel.d_v_z);
	}
	else
	{
		m_outGas->CopyFromStream(_time, m_inGas);
		m_outParticle->CopyFromStream(_time, m_outParticle);
		RaiseError("No valid results available for fluid dynamics model.");
		return;
	}
	/// State variables ///
	SetStateVariable("Bed mass [kg]", m_FDNLModel.GetInitialBedMass(), _time);
	SetStateVariable("u_mf [m/s]", m_FDNLModel.GetUmf(), _time);
	SetStateVariable("Geldart type", E2D(m_FDNLModel.GetGeldart().type), _time);

	////////////////////////////////////////////////////////////////////////////////////
	//                    Solving heat and mass transfer                              //
	////////////////////////////////////////////////////////////////////////////////////
	// Initialize model parameters
	m_TDNLModel.Initialize(_time, TDParams, ResultsFDModel);
	// Run model calculations
	if (!m_TDNLSolver.Calculate(_time))
		RaiseError(m_TDNLSolver.GetError());

	// Results of TD model
	SResultsNLModelTD ResultsTDModel = m_TDNLModel.GetTDModelResults();

	if (ResultsTDModel.valid)
	{
		/// TD model attributes ///
		size_t RTGridPoints = m_TDNLModel.GetRTGridPoints();
		size_t heightGridPoints = m_TDNLModel.GetHeightGridPoints();
		size_t sizeClasses = m_TDNLModel.GetSizeClasses();
		size_t moistureClasses = m_TDNLModel.GetMoistureClasses();
		std::vector<double> h_z = m_TDNLModel.GetHeightVector();
		std::vector<double> tau_i = m_TDNLModel.GetRTTauVector();
		std::vector<double> n_i = m_TDNLModel.GetRTNumberDensityVector();
		/// Plotting ///
		// Gas
		AddCurveOnPlot("Y_g"  , "Y_s_t-"   + StringFunctions::Double2String(_time));
		AddCurveOnPlot("Y_g"  , "Y_b_t-"   + StringFunctions::Double2String(_time));
		AddCurveOnPlot("phi_g", "phi_s_t-" + StringFunctions::Double2String(_time));
		AddCurveOnPlot("phi_g", "phi_b_t-" + StringFunctions::Double2String(_time));
		AddCurveOnPlot("h_g"  , "h_s_t-"   + StringFunctions::Double2String(_time));
		AddCurveOnPlot("h_g"  , "h_b_t-"   + StringFunctions::Double2String(_time));
		AddCurveOnPlot("T_g"  , "T_s_t-"   + StringFunctions::Double2String(_time));
		AddCurveOnPlot("T_g"  , "T_b_t-"   + StringFunctions::Double2String(_time));
		for (size_t z = 0; z < heightGridPoints; ++z)
		{
			double Y_s_z = ResultsTDModel.Y_s_z[z];
			double Y_b_z = ResultsTDModel.Y_b_z[z];
			double T_s_z = ResultsTDModel.T_s_z[z];
			double T_b_z = ResultsTDModel.T_b_z[z];
			AddPointOnCurve("Y_g"  , "Y_s_t-"   + StringFunctions::Double2String(_time), h_z[z], Y_s_z);
			AddPointOnCurve("Y_g"  , "Y_b_t-"   + StringFunctions::Double2String(_time), h_z[z], Y_b_z);
			AddPointOnCurve("phi_g", "phi_s_t-" + StringFunctions::Double2String(_time), h_z[z], m_TDNLModel.GetRelativeHumidity(Y_s_z, T_s_z, 101325) * 100);
			AddPointOnCurve("phi_g", "phi_b_t-" + StringFunctions::Double2String(_time), h_z[z], m_TDNLModel.GetRelativeHumidity(Y_s_z, T_s_z, 101325) * 100);
			AddPointOnCurve("h_g"  , "h_s_t-"   + StringFunctions::Double2String(_time), h_z[z], ResultsTDModel.h_s_z[z]);
			AddPointOnCurve("h_g"  , "h_b_t-"   + StringFunctions::Double2String(_time), h_z[z], ResultsTDModel.h_b_z[z]);
			AddPointOnCurve("T_g"  , "T_s_t-"   + StringFunctions::Double2String(_time), h_z[z], T_s_z - 273.15);
			AddPointOnCurve("T_g"  , "T_b_t-"   + StringFunctions::Double2String(_time), h_z[z], T_b_z - 273.15);
		}

		AddCurveOnPlot("Number density distribution", _time);
		AddPointsOnCurve("Number density distribution", _time, tau_i, n_i);

		AddCurveOnPlot("Bubble volume fraction nu", _time);
		AddPointsOnCurve("Bubble volume fraction nu", _time, ResultsFDModel.h_z, ResultsFDModel.nu_z);

		// Particles
		if (_time == 0)
		{
			std::vector<double> X_p_in = GetClassesMeans(DISTR_MOISTURE);
			std::vector<double> d_d = GetClassesMeans(DISTR_SIZE);
			std::vector<double> n_i_dv = m_TDNLModel.GetRTNumberDensityVector();
			std::vector<double> Q3_X(RTGridPoints, 100);
			for (size_t i = 1; i < RTGridPoints; ++i)
			{
				Q3_X[i] = Q3_X[i - 1] - 0.5 * (n_i_dv[i] + n_i_dv[i - 1]) * (tau_i[i] - tau_i[i - 1]) * 100;
			}

			for (size_t d = 0; d < sizeClasses; ++d)
			{
				for (size_t f = 0; f < moistureClasses; ++f)
				{
					AddCurveOnPlot("X_p (X_in- " + std::to_string(X_p_in[f]) + ")", d_d[d]);
					AddCurveOnPlot("h_p (X_in- " + std::to_string(X_p_in[f]) + ")", d_d[d]);
					AddCurveOnPlot("T_p (X_in- " + std::to_string(X_p_in[f]) + ")", d_d[d]);
					AddPointsOnCurve("X_p (X_in- " + std::to_string(X_p_in[f]) + ")", d_d[d], tau_i, ResultsTDModel.X_p_dfi[d][f]);
					AddPointsOnCurve("h_p (X_in- " + std::to_string(X_p_in[f]) + ")", d_d[d], tau_i, ResultsTDModel.h_p_dfi[d][f]);
					std::vector<double> T_p_i(RTGridPoints);
					for (size_t i = 0; i < RTGridPoints; ++i)
					{
						T_p_i[i] = ResultsTDModel.T_p_dfi[d][f][i] - 273.15;
					}
					AddPointsOnCurve("T_p (X_in- " + std::to_string(X_p_in[f]) + ")", d_d[d], tau_i, T_p_i);

					AddCurveOnPlot("Q_3_X (X_in- " + std::to_string(X_p_in[f]) + ")", d_d[d]);
					AddPointsOnCurve("Q_3_X (X_in- " + std::to_string(X_p_in[f]) + ")", d_d[d], ResultsTDModel.X_p_dfi[d][f], Q3_X);
				}
			}
		}

		// Heat transfer coefficient
		if (_time == 0)
		{
			std::vector<double> d_d = GetClassesMeans(DISTR_SIZE);

			for (size_t d = 0; d < sizeClasses; ++d)
			{
				AddCurveOnPlot("Heat transfer", d_d[d]);
				AddPointsOnCurve("Heat transfer", d_d[d], h_z, ResultsTDModel.alpha_ps_dz[d]);
			}
		}

		AddCurveOnPlot("ApAbed", _time);
		AddPointsOnCurve("ApAbed", _time, m_TDNLModel.GetMeanDiameterVector(), m_TDNLModel.GetRatioApAbedVector());

		/// State variables ///
		SetStateVariable("Deviation H2O-Balance [%]", ResultsTDModel.BalanceH2ODeviation * 100, _time);
		SetStateVariable("Deviation Enthalpy-Balance [%]", ResultsTDModel.BalanceEnthalpyDeviation * 100, _time);
		SetStateVariable("Tau_max [s]", tau_i[RTGridPoints - 1], _time);
		SetStateVariable("T_Particle [K]", ResultsTDModel.T_p_out, _time);
		SetStateVariable("X_Particle [kg/kg]", ResultsTDModel.X_p_out, _time);
		// Calculate mean value of gas distribution over the bed height
		double nu_mean = 0;
		for (size_t z = 1; z < ResultsFDModel.zeta_z.size(); ++z)
		{
			nu_mean += 0.5 * (ResultsFDModel.nu_z[z] + ResultsFDModel.nu_z[z - 1]) * (ResultsFDModel.zeta_z[z] - ResultsFDModel.zeta_z[z - 1]);
		}
		SetStateVariable("nu_mean [-]", nu_mean, _time);
		SetStateVariable("Q_Wall [W]", ResultsTDModel.Q_pw, _time);

		double T_s_mean = 0;
		for (size_t z = 1; z < ResultsTDModel.T_s_z.size(); ++z)
		{
			T_s_mean += 0.5 * (ResultsTDModel.T_s_z[z] + ResultsTDModel.T_s_z[z - 1]);
		}
		T_s_mean /= static_cast<double>(ResultsTDModel.T_s_z.size() - 1);
		double v_s = m_FDNLModel.GetU0();
		double alpha_SW_mean = m_TDNLModel.CalculateAlpha_SW(T_s_mean, 0, 101325, ResultsFDModel.h_fb, v_s);
		SetStateVariable("alpha_SW_mean [W/(m^2*K)]", alpha_SW_mean, _time);
		double alpha_PW_mean = m_TDNLModel.CalculateAlpha_PW(ResultsTDModel.T_p_out, ResultsTDModel.X_p_out, T_s_mean, 101325);
		SetStateVariable("alpha_PW_mean [W/(m^2*K)]", alpha_PW_mean, _time);

		/// Set outlet streams ///
		// Gas flow
		m_outGas->CopyFromStream(_time, m_inGas);
		m_outGas->SetMassFlow(_time, 0);
		double mflow_g_in = m_inGas->GetPhaseMassFlow(_time, EPhase::VAPOR);
		double Y_g_in = ResultsTDModel.Y_s_z[0];
		double Y_g_out = ResultsTDModel.Y_s_z[heightGridPoints - 1];
		double mflow_g_out = mflow_g_in / (1 + Y_g_in) * (1 + Y_g_out);
		m_outGas->SetPhaseMassFlow(_time, EPhase::VAPOR, mflow_g_out);
		m_outGas->SetCompoundFraction(_time, m_compKeys.air, EPhase::VAPOR, 1 / (1 + Y_g_out));
		m_outGas->SetCompoundFraction(_time, m_compKeys.H2Og, EPhase::VAPOR, Y_g_out / (1 + Y_g_out));
		m_outGas->SetTemperature(_time, ResultsTDModel.T_g_out);
		// Particle flow
		m_outParticle->CopyFromStream(_time, m_inParticle);
		m_outParticle->SetTemperature(_time, ResultsTDModel.T_p_mean);

		if (_time > 0)
			WriteCSV();
	}
	else
	{
		m_outGas->CopyFromStream(_time, m_inGas);
		m_outParticle->CopyFromStream(_time, m_outParticle);
		RaiseWarning("No valid results available for thermodynamics model.");
	}
}

void CVibratedFluidizedBedDryer::SaveState()
{
	m_FDNLSolver.SaveState();
	m_TDNLSolver.SaveState();
}

void CVibratedFluidizedBedDryer::LoadState()
{
	m_FDNLSolver.LoadState();
	m_TDNLSolver.LoadState();
}

/// Property functions ///
double CVibratedFluidizedBedDryer::GetMolarMass(const std::string& _compKey) const
{
	return GetCompoundProperty(_compKey, ECompoundConstProperties::MOLAR_MASS);
}
double CVibratedFluidizedBedDryer::GetCompoundGasDensity(double _t_eval, double _p_eval, const std::string& _compKey) const
{
	const double MM_air = GetCompoundProperty(_compKey, ECompoundConstProperties::MOLAR_MASS);
	return _p_eval * MM_air / (MOLAR_GAS_CONSTANT * _t_eval);
}
double CVibratedFluidizedBedDryer::GetCompoundViscosity(double _t_eval, double _p_eval, const std::string& _compKey) const
{
	return GetCompoundProperty(_compKey, ECompoundTPProperties::VISCOSITY, _t_eval, _p_eval);
}
double CVibratedFluidizedBedDryer::GetCompoundHeatCapacity(double _t_eval, double _p_eval, const std::string& _compKey) const
{
	return GetCompoundProperty(_compKey, ECompoundTPProperties::HEAT_CAPACITY_CP);
}
double CVibratedFluidizedBedDryer::GetCompoundThermalConductivity(double _t_eval, double _p_eval, const std::string& _compKey) const
{
	return GetCompoundProperty(_compKey, ECompoundTPProperties::THERMAL_CONDUCTIVITY, _t_eval, _p_eval);
}

double CVibratedFluidizedBedDryer::GetMoistGasDensity(double _t, double _p, double _y, const std::string& _keyAir) const
{
	const double MM_air = GetCompoundProperty(_keyAir, ECompoundConstProperties::MOLAR_MASS);
	const double rho_air = _p * MM_air / (MOLAR_GAS_CONSTANT * _t);
	return rho_air * (1 + _y);
}

double CVibratedFluidizedBedDryer::GetMoistGasEnthalpy(double _t, double _p, double _y, double _t_0, double _hv_H2O, const std::string& _keyAir, const std::string& _keyH2Ol) const
{
	const double cp_air = GetCompoundProperty(_keyAir, ECompoundTPProperties::HEAT_CAPACITY_CP);
	const double cp_H2O_l = GetCompoundProperty(_keyH2Ol, ECompoundTPProperties::HEAT_CAPACITY_CP);

	const double h_air_dry = cp_air * (_t - _t_0);
	const double h_H2O_g = cp_H2O_l * (_t - _t_0) + _hv_H2O;

	const double h_air_moist = h_air_dry + _y * h_H2O_g;
	return h_air_moist;
}
double CVibratedFluidizedBedDryer::GetMoistGasTemperature(double _h, double _y, double _t_eval, double _p_eval, double _t_0, double _hv_H2O, const std::string& _keyAir, const std::string& _keyH2Ol) const
{
	const double cp_air = GetCompoundProperty(_keyAir, ECompoundTPProperties::HEAT_CAPACITY_CP);
	const double cp_H2O_l = GetCompoundProperty(_keyH2Ol, ECompoundTPProperties::HEAT_CAPACITY_CP);

	const double T_air_moist = _t_0 + (_h - _y * _hv_H2O) / (cp_air + _y * cp_H2O_l);
	return T_air_moist;
}
double CVibratedFluidizedBedDryer::GetRelativeHumidity(double _y, double _t, double _p, const std::string& _keyAir, const std::string& _keyH2Ol) const
{
	const double M_air = 1.;				// Reference mass of air [kg]
	const double M_wg = M_air * _y;	// Mass of gaseous water [kg]

	// TO DO
	const double MM_air = GetCompoundProperty(_keyAir, ECompoundConstProperties::MOLAR_MASS); // Molar mass of air [kg/mol]
	const double MM_H2Og = GetCompoundProperty(_keyH2Ol, ECompoundConstProperties::MOLAR_MASS); // Molar mass of water [kg/mol]

	const double n_air = M_air / MM_air;	// Amount of substance of air [mol]
	const double n_wg = M_wg / MM_H2Og;		// Amount of substance of gaseous water [mol]

	const double P_wg = n_wg / (n_wg + n_air) * _p; // Partial pressure of gaseous water [Pa]
	const double P_sat = GetCompoundProperty(_keyH2Ol, ECompoundTPProperties::VAPOR_PRESSURE, _t, _p); // Saturated vapor pressure of water [Pa]

	const double phi = P_wg / P_sat; // Relative humidity [-]

	return phi;
}
double CVibratedFluidizedBedDryer::GetGasEquilibriumMoistureContent(double _t_p, double _phi_eq, double _p, const std::string& _keyAir, const std::string& _keyH2Ol, const std::string& _keyH2Og) const
{
	// TO DO hygroscopic/not hygroscopic materials
	double dP_sat = GetCompoundProperty(_keyH2Ol, ECompoundTPProperties::VAPOR_PRESSURE, _t_p, _p); // Saturation partial pressure of gaseous water at particle temperature [Pa]

	// TO DO
	double dMM_air = GetCompoundProperty(_keyAir, ECompoundConstProperties::MOLAR_MASS); // Molar mass of air [kg/mol]
	double dMM_H2Og = GetCompoundProperty(_keyH2Og, ECompoundConstProperties::MOLAR_MASS); // Molar mass of water [kg/mol]

	double dP_eq = dP_sat * _phi_eq;							// Equilibrium partial pressure [Pa]

	// TO DO warning
	// Check non positive pressure difference
	if (_p - dP_eq <= 0)
		return 0;

	double dY_eq = dMM_H2Og / dMM_air * dP_eq / (_p - dP_eq);	// Eqilibrium moisture content [kg/kg]
	return dY_eq;
}
double CVibratedFluidizedBedDryer::GetDiffusionCoefficientOld(double _t, double _p, const std::string& _keyAir, const std::string& _keyH2Og) const
{
	// References:
	// [1] Poling, Bruce E., et al. The properties of gases and liquids. Vol. 5. New York: Mcgraw-hill, 2001.
	// [2] Fuller, E. N., and J. C. Giddings. "A comparison of methods for predicting gaseous diffusion coefficients." Journal of Chromatographic Science 3.7 (1965): 222-227.
	// [3] Fuller, Edward N., Paul D. Schettler, and J. Calvin Giddings. "New method for prediction of binary gas-phase diffusion coefficients." Industrial & Engineering Chemistry 58.5 (1966): 18-27.
	// [4] Fuller, Edward N., Keith Ensley, and J. Calvin Giddings. "Diffusion of halogenated hydrocarbons in helium. The effect of structure on collision cross sections." The Journal of Physical Chemistry 73.11 (1969): 3679-3685.

	double dP_bar = _p / 1e5;	// Pressure [bar]
	// TO DO
	double dMM_air = GetCompoundProperty(_keyAir, ECompoundConstProperties::MOLAR_MASS); // Molar mass of air [kg/mol]
	double dMM_H2Og = GetCompoundProperty(_keyH2Og, ECompoundConstProperties::MOLAR_MASS); // Molar mass of water [kg/mol]

	// TO DO check
	double dM_wg = 2. / 1. / dMM_air + 1. / dMM_H2Og;		// Mean molar mass (?) [kg/mol]
	double dSigma_nu_w = 2. * 2.31 + 6.11;	// Summation of atomic diffusion volumes of water
	double dSigma_nu_g = 19.7;				// Summation of atomic diffusion volumes of air
	double ddelta_wg_cm2s = 0.00143 * pow(_t, 1.75) / dP_bar / pow(dM_wg, 0.5) / pow(pow(dSigma_nu_w, 1. / 3.) + pow(dSigma_nu_g, 1. / 3.), 2.); // Diffusion coefficient D_wg in [cm^2/s]
	double ddelta_wg_m2s = ddelta_wg_cm2s * 1e-4; // Diffusion coefficient delta_wg [m^2/s]

	return ddelta_wg_m2s;
}
double CVibratedFluidizedBedDryer::GetDiffusionCoefficient(double _t, double _p) const
{
	// Model used in dissertation of Burgschweiger (2000)
	// From Appendix A 2.3
	return 2.252 / _p * pow((_t / 273), 1.81);
}

double CVibratedFluidizedBedDryer::GetWaterVaporEnthalpy(double _t, double _p, double _t_0, double _hv_H2O, const std::string& _keyH2Ol) const
{
	double dcp_H2O_l = GetCompoundProperty(_keyH2Ol, ECompoundTPProperties::HEAT_CAPACITY_CP);
	return dcp_H2O_l * (_t - _t_0) + _hv_H2O;
}

// Particle properties
double CVibratedFluidizedBedDryer::GetDryParticleDensity(double _t_eval, double _p_eval, const std::string& _keyParticle) const
{
	return GetCompoundProperty(_keyParticle, ECompoundTPProperties::DENSITY, _t_eval, _p_eval);
}
double CVibratedFluidizedBedDryer::GetMoistParticleDensity(double _x_p, double _t_eval, double _p_eval, const std::string& _keyParticle) const
{
	double drho_s = GetCompoundProperty(_keyParticle, ECompoundTPProperties::DENSITY, _t_eval, _p_eval);
	return drho_s * (1 + _x_p);
}

double CVibratedFluidizedBedDryer::GetMoistParticleHeatCapacity(double _x_p, double _t, double _p, const std::string& _keyParticle, const std::string& _keyH2Ol) const
{
	double dcp_s = GetCompoundProperty(_keyParticle, ECompoundTPProperties::HEAT_CAPACITY_CP);
	double dcp_H2O_l = GetCompoundProperty(_keyH2Ol, ECompoundTPProperties::HEAT_CAPACITY_CP);
	return dcp_s + _x_p * dcp_H2O_l;
}
double CVibratedFluidizedBedDryer::GetMoistParticleEnthalpy(double _t, double _p, double _x, double _t_0, const std::string& _keyParticle, const std::string& _keyH2Ol) const
{
	double dcp_s = GetCompoundProperty(_keyParticle, ECompoundTPProperties::HEAT_CAPACITY_CP);
	double dcp_H2O_l = GetCompoundProperty(_keyH2Ol, ECompoundTPProperties::HEAT_CAPACITY_CP);

	double dh_p_dry = dcp_s * (_t - _t_0);
	double dh_H2O_l = dcp_H2O_l * (_t - _t_0);

	double dh_p_moist = dh_p_dry + _x * dh_H2O_l;
	return dh_p_moist;
}
double CVibratedFluidizedBedDryer::GetMoistParticleTemperature(double _h, double _x, double _t_eval, double _p_eval, double _t_0, const std::string& _keyParticle, const std::string& _keyH2Ol) const
{
	double dcp_s = GetCompoundProperty(_keyParticle, ECompoundTPProperties::HEAT_CAPACITY_CP);
	double dcp_H2O_l = GetCompoundProperty(_keyH2Ol, ECompoundTPProperties::HEAT_CAPACITY_CP);

	double dT_p_moist = _t_0 + _h / (dcp_s + _x * dcp_H2O_l);
	return dT_p_moist;
}
std::vector<double> CVibratedFluidizedBedDryer::GetMeanDiameterVector() const
{
	return GetClassesMeans(EDistrTypes::DISTR_SIZE);
}
std::vector<double> CVibratedFluidizedBedDryer::GetMeanMoistureVector() const
{
	return GetClassesMeans(EDistrTypes::DISTR_MOISTURE);
}

/// Auxiliary functions ///

void CVibratedFluidizedBedDryer::WriteCSV()
{
	const std::string csvPath = m_TDNLModel.GetCSVPath();
	if (csvPath.empty()) return;

	std::filesystem::create_directories(csvPath);
	if (!std::filesystem::exists(csvPath))
		RaiseWarning("Directory for writing out results " + StringFunctions::Quote(csvPath) + " was not found.");
	else
	{
		unsigned count = 1;
		char buffer[4];

		while (true)
		{
			snprintf(buffer, 4, "%03d", count);
			std::string resFolderPath = csvPath + "/Results_Run" + std::string(buffer);
			if (std::filesystem::exists(resFolderPath))
			{
				count++;
			}
			else
			{
				std::filesystem::create_directory(resFolderPath);
				WriteParameters(resFolderPath);
				WriteStateVariables(resFolderPath);
				WritePlots(resFolderPath);
				break;
			}
		}
	}
}

void CVibratedFluidizedBedDryer::WriteParameters(const std::string& _resultsPath)
{
	std::vector<std::vector<std::string>> params;

	for (const auto& param : GetUnitParametersManager().GetParameters())
	{
		std::string name = param->GetName();
		const EUnitParameter type = param->GetType();

		std::string val;
		switch (type)
		{
		case EUnitParameter::CONSTANT:
			val = std::to_string(GetConstRealParameterValue(name));
			break;
		case EUnitParameter::STRING:
			val = GetStringParameterValue(name);
			break;
		case EUnitParameter::COMPOUND:
			val = GetCompoundParameterValue(name);
			break;
		default:
			val = "";
			break;
		}
		params.push_back(std::vector<std::string>{ name, val });
	}

	const std::string filePath = _resultsPath + "/Parameters.csv";
	if (!std::filesystem::exists(filePath))
	{
		std::ofstream fileCSV;
		fileCSV.open(filePath);
		writerCSV2D(params, "Parameters", fileCSV);
		fileCSV.close();
	}
}

void CVibratedFluidizedBedDryer::WriteStateVariables(const std::string& _resultsPath)
{
	const auto& manager = GetStateVariablesManager();
	const size_t numberSV = manager.GetStateVariablesNumber();

	std::vector<std::vector<std::string>> vars;

	for (const auto& sv : GetStateVariablesManager().GetAllStateVariables())
	{
		const std::string name = sv->GetName();
		const std::string val = std::to_string(sv->GetValue());
		vars.push_back(std::vector<std::string>{ name , val });
	}

	const std::string filePath = _resultsPath + "/StateVariables.csv";
	if (!std::filesystem::exists(filePath))
	{
		std::ofstream fileCSV;
		fileCSV.open(filePath);
		writerCSV2D(vars, "StateVariables", fileCSV);
		fileCSV.close();
	}
}

void CVibratedFluidizedBedDryer::WritePlots(const std::string& _resultsPath) const
{
	for (const auto& plot : GetPlotsManager().GetAllPlots())
	{
		std::string nPlotname = plot->GetName();

		std::vector<std::vector<std::string>> vals;

		vals.push_back(strVector(plot->GetLabelX(), plot->GetCurve(static_cast<size_t>(0))->GetXValues()));

		if (plot->GetCurvesNumber() == 0)
			continue;

		for (const auto& curve : plot->GetAllCurves())
		{
			vals.push_back(strVector(curve->GetName() + "_" + plot->GetLabelY(), curve->GetYValues()));
		}

		transpose(vals);

		const std::string filePath = _resultsPath + "/Plot_" + nPlotname + ".csv";
		if (!std::filesystem::exists(filePath))
		{
			std::ofstream fileCSV;
			fileCSV.open(filePath);
			writerCSV2D(vals, nPlotname, fileCSV);
			fileCSV.close();
		}
	}
}


//////////////////////////////////////////////////////////////////////////
/// Base Model

void CNLModelBase::Setup(void* _unit, const SCompoundKeys& _compoundKeys, const SConstantParameters& _constants, const SPorts& _inPort, const SDiscretizationParameters& _discretization, const SSettings& _settings)
{
	m_unit           = _unit;
	m_compoundKeys   = _compoundKeys;
	m_constants      = _constants;
	m_inPorts        = _inPort;
	m_discretization = _discretization;
	m_settings       = _settings;
}


/// Calculation functions ///
STDParameters CNLModelBase::CalculateTDParameters(double _time, EGeldart _geldart)
{
	const SMYTP MYTPInGas{ m_inPorts.gas, _time, EPhase::VAPOR, m_compoundKeys.air, m_compoundKeys.H2Og };
	const SMYTP MYTPInParticle{ m_inPorts.particle, _time, EPhase::SOLID, m_compoundKeys.particles, m_compoundKeys.H2Ol };

	const double mflow_g_dry = MYTPInGas.m;
	const double A_fb = GetBedArea();
	const double d_p_sauter = GetInParticleSauter(_time);
	const double rho_p = GetInParticleDensity(_time);
	const double rho_g = GetInGasDensity(_time);
	const double eta_g = GetInGasViscosity(_time);
	const double nu_g = eta_g / rho_g;
	const double Ar = STANDARD_ACCELERATION_OF_GRAVITY * (d_p_sauter * d_p_sauter * d_p_sauter) * (rho_p - rho_g) / (nu_g * nu_g * rho_g);	// Archimedes number [-]
	double Re_mf = 33.7 * (sqrt(1. + 3.6e-5 * Ar) - 1.);		// Reynolds number, Wen & Yu [-]
	const double u_mf_man = GetManualMinimumVelocity();
	double u_mf = Re_mf * nu_g / d_p_sauter;				// Minimum fluidization velocity at inlet conditions [m/s]
	if (u_mf_man > 0)
	{
		u_mf = u_mf_man;
		Re_mf = u_mf * d_p_sauter / nu_g;
	}

	const double u_0 = mflow_g_dry / rho_g / A_fb;			// Gas velocity at inlet conditions[m/s]

	SGeldart geldart;
	if (_geldart == EGeldart::A || _geldart == EGeldart::B || _geldart == EGeldart::C || _geldart == EGeldart::D)
		geldart = SGeldart(_geldart);
	else
		geldart = SGeldart(d_p_sauter, rho_p, rho_g);


	const double theta = CalculateTheta(geldart);
	const SFluidDynamics FD = SFluidDynamics(Re_mf, u_mf, u_0, theta);

	return STDParameters(MYTPInGas, MYTPInParticle, FD, geldart);
}
SBoundaryCondition CNLModelBase::CalculateBoundaryConditions(double _time)
{
	double dT_in_g = GetInGasMYTP(_time).T;
	double dh_in_g = GetInGasEnthalpy(_time);
	double dY_in_g = GetInGasMoisture(_time);

	double dT_in_p = GetInParticleMYTP(_time).T;
	std::vector<double> vh_in_p_f = GetInParticleEnthalpy(_time);
	std::vector<double> vX_in_p_f = GetMeanMoistureVector();

	return SBoundaryCondition(dT_in_g, dh_in_g, dY_in_g, dT_in_p, vh_in_p_f, vX_in_p_f);
}

SBoundaryCondition CNLModelBase::CalculateInitialConditions(SBoundaryCondition _bc)
{
	double dT_in_g = _bc.T_g;
	double dh_in_g = _bc.h_g;
	double dY_in_g = _bc.Y_g;

	double dT_in_p = _bc.T_p;
	std::vector<double> vh_in_p_f = _bc.h_p;
	std::vector<double> vX_in_p_f = _bc.X_p;

	double dT_mean = 0.5 * (dT_in_g + dT_in_p);
	double dY = 1.2 * dY_in_g;
	double dh_g = GetMoistGasEnthalpy(dT_mean, 101325, dY);
	std::vector<double> dX = std::vector<double>(_bc.X_p.size(), 0);
	std::vector<double> dh_p = std::vector<double>(_bc.h_p.size(), 0);
	for (size_t f = 0; f < _bc.X_p.size(); ++f)
	{
		dX[f] = 0.8 * _bc.X_p[f];
		dh_p[f] = GetMoistParticleEnthalpy(dT_mean, 101325, dX[f]);
	}

	return SBoundaryCondition(dT_mean, dh_g, dY, dT_mean, dh_p, dX);
}

double CNLModelBase::CalculateTheta(const SGeldart& _geldart) const
{
	double dd_fb = GetBedDiameter();
	double dtheta = 0;
	switch (_geldart.type)
	{
	case EGeldart::A:
		if (dd_fb <= 0.05)
			dtheta = 1.18;
		else if ((dd_fb >= 0.05) && (dd_fb <= 1.))
			dtheta = 3.2 * pow(dd_fb, 1. / 3.);
		else if (dd_fb > 1.)
			dtheta = 3.2;
		break;
	case EGeldart::B:
		if (dd_fb <= 0.1)
			dtheta = 0.63;
		else if ((dd_fb >= 0.1) && (dd_fb <= 1.))
			dtheta = 2 * pow(dd_fb, 0.5);
		else if (dd_fb > 1.)
			dtheta = 2.;
		break;
	case EGeldart::C:
		// TO DO
		if (dd_fb <= 0.05)
			dtheta = 1.18;
		else if ((dd_fb >= 0.05) && (dd_fb <= 1.))
			dtheta = 3.2 * pow(dd_fb, 1. / 3.);
		else if (dd_fb > 1.)
			dtheta = 3.2;
		break;
	case EGeldart::D:
		dtheta = 0.87;
		break;
	default:
		break;
	}
	return dtheta;
}
double CNLModelBase::CalculateParticleMeanMoisture(double _time, const std::vector<double>& _mfrac_f)
{
	std::vector<double> vX_p = GetMeanMoistureVector();
	double dX_p_mean = 0;
	for (size_t f = 0; f < vX_p.size(); ++f)
	{
		dX_p_mean += _mfrac_f[f] * vX_p[f];
	}
	return dX_p_mean;
}

/// Stream functions ///
double CNLModelBase::GetInGasMoisture(double _time) const
{
	std::string sKeyAir = m_compoundKeys.air;
	std::string sKeyH2Og = m_compoundKeys.H2Og;
	double dmflow_Air = m_inPorts.gas->GetCompoundMassFlow(_time, sKeyAir, EPhase::VAPOR);
	double dmflow_H2Og = m_inPorts.gas->GetCompoundMassFlow(_time, sKeyH2Og, EPhase::VAPOR);
	return dmflow_H2Og / dmflow_Air;
}
double CNLModelBase::GetInGasEnthalpy(double _time)
{
	double dY = GetInGasMoisture(_time);
	double dT = m_inPorts.gas->GetTemperature(_time);
	double dP = m_inPorts.gas->GetPressure(_time);
	return GetMoistGasEnthalpy(dT, dP, dY);
}
std::vector<std::vector<double>> CNLModelBase::GetInParticleMassfractions(double _time) const
{
	std::vector<double> vmfrac_d = m_inPorts.particle->GetDistribution(_time, EDistrTypes::DISTR_SIZE);
	std::vector<double> vmfrac_f = m_inPorts.particle->GetDistribution(_time, EDistrTypes::DISTR_MOISTURE);

	std::vector<std::vector<double>> vmfrac_df(vmfrac_d.size(), std::vector<double>(vmfrac_f.size()));
	for (size_t d = 0; d < vmfrac_d.size(); ++d)
	{
		for (size_t f = 0; f < vmfrac_f.size(); ++f)
		{
			vmfrac_df[d][f] = vmfrac_d[d] * vmfrac_f[f];
		}
	}
	return vmfrac_df;
}
double CNLModelBase::GetInParticleMeanMoisture(double _time)
{
	std::vector<double> vmfrac_f = m_inPorts.particle->GetDistribution(_time, EDistrTypes::DISTR_MOISTURE);
	return CalculateParticleMeanMoisture(_time, vmfrac_f);
}
std::vector<double> CNLModelBase::GetInParticleEnthalpy(double _time)
{
	size_t nMoistureClasses = GetMoistureClasses();
	std::vector<double> vX_f = GetMeanMoistureVector();

	double dT = m_inPorts.particle->GetTemperature(_time);
	double dP = m_inPorts.particle->GetPressure(_time);
	std::vector<double> vh_f(nMoistureClasses);
	for (size_t f = 0; f < nMoistureClasses; ++f)
	{
		vh_f[f] = GetMoistParticleEnthalpy(dT, dP, vX_f[f]);
	}
	return vh_f;
}
double CNLModelBase::GetInParticleDensity(double _time) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	std::string sParticles = m_compoundKeys.particles;
	double dT = m_inPorts.gas->GetTemperature(_time);
	double dP = m_inPorts.gas->GetTemperature(_time);
	return unit->GetCompoundProperty(sParticles, ECompoundTPProperties::DENSITY, dT, dP);
}
double CNLModelBase::GetInParticleSauter(double _time) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	std::vector<double> vSizeGrid = unit->GetNumericGrid(DISTR_SIZE);
	std::vector<double> vq3_p = m_inPorts.particle->GetPSD(_time, PSD_q3);
	return GetSauterDiameter(vSizeGrid, vq3_p);
}
double CNLModelBase::GetInGasDensity(double _time) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	std::string sKeyAir = m_compoundKeys.air;
	double dT = m_inPorts.gas->GetTemperature(_time);
	double dP = m_inPorts.gas->GetTemperature(_time);
	return unit->GetCompoundProperty(sKeyAir, ECompoundTPProperties::DENSITY, dT, dP);
}
double CNLModelBase::GetInGasViscosity(double _time) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	std::string sKeyAir = m_compoundKeys.air;
	double dT = m_inPorts.gas->GetTemperature(_time);
	double dP = m_inPorts.gas->GetTemperature(_time);
	double deta = unit->GetCompoundProperty(sKeyAir, ECompoundTPProperties::VISCOSITY, dT, dP);
	return deta;
}

/// Material property functions ///
// Constant properties
double CNLModelBase::GetGasMolarMass() const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	std::string sKeyAir = m_compoundKeys.air;
	return unit->GetMolarMass(sKeyAir);
}
double CNLModelBase::GetDiffusionCoefficient(double _t, double _p) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	return unit->GetDiffusionCoefficient(_t, _p);
};

// Gas properties
double CNLModelBase::GetDryGasDensity(double _t, double _p) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	std::string sKeyAir = m_compoundKeys.air;
	return unit->GetCompoundGasDensity(_t, _p, sKeyAir);
}
double CNLModelBase::GetMoistGasDensity(double _t, double _p, double _y) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	std::string sKeyAir = m_compoundKeys.air;
	return unit->GetMoistGasDensity(_t, _p, _y, sKeyAir);
}
double CNLModelBase::GetDryGasViscosity(double _t, double _p) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	std::string sKeyAir = m_compoundKeys.air;
	return unit->GetCompoundViscosity(_t, _p, sKeyAir);
}
double CNLModelBase::GetDryGasHeatCapacity(double _t, double _p) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	std::string sKeyAir = m_compoundKeys.air;
	return unit->GetCompoundHeatCapacity(_t, _p, sKeyAir);
}
double CNLModelBase::GetWaterHeatCapacity(double _t, double _p) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	std::string sKeyH2O = m_compoundKeys.H2Ol;
	return unit->GetCompoundHeatCapacity(_t, _p, sKeyH2O);
}
double CNLModelBase::GetParticleHeatCapacity(double _t, double _p) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	std::string sKeyParticle = m_compoundKeys.particles;
	return unit->GetCompoundHeatCapacity(_t, _p, sKeyParticle);
}
double CNLModelBase::GetDryGasThermalConductivity(double _t, double _p) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	std::string sKeyAir = m_compoundKeys.air;
	return unit->GetCompoundThermalConductivity(_t, _p, sKeyAir);
}
double CNLModelBase::GetGasEquilibriumMoistureContent(double _t_p, double _phi_eq, double _p) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	std::string sKeyAir = m_compoundKeys.air;
	std::string sKeyH2Ol = m_compoundKeys.H2Ol;
	std::string sKeyH2Og = m_compoundKeys.H2Og;
	return unit->GetGasEquilibriumMoistureContent(_t_p, _phi_eq, _p, sKeyAir, sKeyH2Ol, sKeyH2Og);
}
double CNLModelBase::GetMoistGasEnthalpy(double _t, double _p, double _y) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	double dT_0 = GetReferenceTemperature();
	double dhv_H2O = GetEvaporationEnthalpy();
	std::string sKeyAir = m_compoundKeys.air;
	std::string sKeyH2Ol = m_compoundKeys.H2Ol;
	return unit->GetMoistGasEnthalpy(_t, _p, _y, dT_0, dhv_H2O, sKeyAir, sKeyH2Ol);
}
double CNLModelBase::GetMoistGasTemperature(double _h, double _y, double _t_prop, double _p_prop) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	double dT_0 = GetReferenceTemperature();
	double dhv_H2O = GetEvaporationEnthalpy();
	std::string sKeyAir = m_compoundKeys.air;
	std::string sKeyH2Ol = m_compoundKeys.H2Ol;
	return unit->GetMoistGasTemperature(_h, _y, _t_prop, _p_prop, dT_0, dhv_H2O, sKeyAir, sKeyH2Ol);
}
double CNLModelBase::GetRelativeHumidity(double _y, double _t, double _p) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	std::string sKeyAir = m_compoundKeys.air;
	std::string sKeyH2Og = m_compoundKeys.H2Og;
	return unit->GetRelativeHumidity(_y, _t, _p, sKeyAir, sKeyH2Og);
}
double CNLModelBase::GetParticleEquilibriumMoistureContent(double _phi, double _t) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	std::string sKeyParticle = m_compoundKeys.particles;
	return unit->GetParticleEquilibriumMoistureContent(_phi, _t, sKeyParticle);
}
double CNLModelBase::GetWaterVaporEnthalpy(double _t, double _p) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	double dT_0 = GetReferenceTemperature();
	double dhv_H2O = GetEvaporationEnthalpy();
	std::string sKeyAir = m_compoundKeys.air;
	std::string sKeyH2Ol = m_compoundKeys.H2Ol;
	return unit->GetWaterVaporEnthalpy(_t, _p, dT_0, dhv_H2O, sKeyH2Ol);
}

// Particle properties
double CNLModelBase::GetDryParticleDensity(double _t, double _p) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	std::string sKeyParticle = m_compoundKeys.particles;
	return unit->GetDryParticleDensity(_t, _p, sKeyParticle);
}
double CNLModelBase::GetMoistParticleDensity(double _x_p, double _t, double _p) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	std::string sKeyParticle = m_compoundKeys.particles;
	return unit->GetMoistParticleDensity(_x_p, _t, _p, sKeyParticle);
}
double CNLModelBase::GetMoistParticleHeatCapacity(double _x_p, double _t, double _p) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	std::string sKeyParticle = m_compoundKeys.particles;
	std::string sKeyH2Ol = m_compoundKeys.H2Ol;
	return unit->GetMoistParticleHeatCapacity(_x_p, _t, _p, sKeyParticle, sKeyH2Ol);
}
double CNLModelBase::GetMoistParticleEnthalpy(double _t, double _p, double _x) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	double dT_0 = GetReferenceTemperature();
	std::string sKeyH2Ol = m_compoundKeys.H2Ol;
	std::string sKeyParticle = m_compoundKeys.particles;
	return unit->GetMoistParticleEnthalpy(_t, _p, _x, dT_0, sKeyParticle, sKeyH2Ol);
}
double CNLModelBase::GetMoistParticleTemperature(double _h, double _x, double _t_prop, double _p_prop) const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	double dT_0 = GetReferenceTemperature();
	std::string sKeyH2Ol = m_compoundKeys.H2Ol;
	std::string sKeyParticle = m_compoundKeys.particles;
	return unit->GetMoistParticleTemperature(_h, _x, _t_prop, _p_prop, dT_0, sKeyParticle, sKeyH2Ol);
}

/// Unit functions ///
std::vector<double> CNLModelBase::GetParticleSizeGrid() const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	return unit->GetPSDGridDiameters();
}

std::vector<double> CNLModelBase::GetMeanDiameterVector() const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	return unit->GetMeanDiameterVector();
}
std::vector<double> CNLModelBase::GetMeanMoistureVector() const
{
	auto* unit = static_cast<CVibratedFluidizedBedDryer*>(m_unit);
	return unit->GetMeanMoistureVector();
}

//////////////////////////////////////////////////////////////////////////
/// Fluid dynamics - Model

CNLModelFD::~CNLModelFD()
{
	ClearModel();
}

/// Declaration ///
void CNLModelFD::AddVariables()
{
	ClearModel();

	// Get attributes
	size_t nFDHeightGridPoints = GetFDHeightGridPoints();
	double dh_0 = GetInitialBedHeight();

	double dd_fb = GetBedDiameter();
	double dA_fb = GetBedArea();
	double du_0 = GetU0();
	size_t nNum_or = GetNumberOfOrifices();

	double dVflow_or = dA_fb * du_0;
	if (nNum_or > 0)
		dVflow_or /= nNum_or;								// V_or, gas volume flow rate through a distributor, Eq. 2.2

	double dd_v_0 = 1.3 * pow(pow(dVflow_or, 2.0) / STANDARD_ACCELERATION_OF_GRAVITY, 0.2); //dv,0, theoretical diameter of a bubble, Eq. 2.1


	// Add NL variables
	for (unsigned z = 0; z < (nFDHeightGridPoints); z++)
	{
		//m_vnd_v_z.push_back(AddNLVariable(dd_v_0, 0.0));
		m_nd_v_z.push_back(AddNLVariable(0, 0.0));
	}
	m_h_fb = AddNLVariable(dh_0, 0.0);
}

void CNLModelFD::Initialize(double _time, const STDParameters& _tdParameters)
{
	SetTDParameters(_tdParameters);
}

/// Solver functions
void CNLModelFD::CalculateFunctions(double* _vars, double* _func, void* _unit)
{
	////////////////////////////////////////////////////////////////////////////////////
	//                              ATTRIBUTES		                                  //
	////////////////////////////////////////////////////////////////////////////////////
	// Discretization
	size_t nFDHeightGridPoints = GetFDHeightGridPoints();
	// Bed input
	double dh_0 = GetInitialBedHeight();
	double deps_0 = GetInitialBedPorosity();
	double dLambda = GetVibrationIntensity();
	// Geometry
	double dL_char_fb = GetBedCharacteristicLength();
	double dA_fb = GetBedArea();
	size_t nNum_or = GetNumberOfOrifices();
	// Operating parameters
	SGeldart sGeldart = GetGeldart();
	double du_0 = GetU0();
	double du_mf = GetUmf();
	double dA_0 = CalculatedA_0_GeldartC();
	////////////////////////////////////////////////////////////////////////////////////
	//                              VECTOR DECLARATION                                //
	////////////////////////////////////////////////////////////////////////////////////
	std::vector<double> vd_v_z(nFDHeightGridPoints);
	std::vector<double> vd_v_z_update(nFDHeightGridPoints);
	std::vector<double> vh_z(nFDHeightGridPoints);
	std::vector<double> vVflow_b_z(nFDHeightGridPoints);
	std::vector<double> veps_b_z(nFDHeightGridPoints);

	std::vector<double> vu_b_z(nFDHeightGridPoints);
	std::vector<double> vPsi_z(nFDHeightGridPoints);


	////////////////////////////////////////////////////////////////////////////////////
	//                              GET VARIABLE VALUES                               //
	////////////////////////////////////////////////////////////////////////////////////
	for (size_t z = 0; z < nFDHeightGridPoints; ++z)
	{
		vd_v_z[z] = _vars[m_nd_v_z[z]];
	}
	double dH_fb = _vars[m_h_fb];
	////////////////////////////////////////////////////////////////////////////////////
	//                              CALCULATIONS                                      //
	////////////////////////////////////////////////////////////////////////////////////
	double dDelta_z = dH_fb / (nFDHeightGridPoints - 1);
	// Height discretization
	for (size_t z = 0; z < nFDHeightGridPoints; ++z)
	{
		vh_z[z] = z * dDelta_z;
	}

	double dtheta = GetTheta();

	double dd_v_0 = 0;
	double dVflow_or = dA_fb * du_0;

	// Calculate initial bubble diameter
	switch (sGeldart.type)
	{

		// Geldart C
	case EGeldart::C:

		dd_v_0 = 0.21 * (pow((du_0 - du_mf), 0.49) * pow((0 + 4 * sqrt(dA_0)), 0.48)) / (pow(STANDARD_ACCELERATION_OF_GRAVITY, 0.2));

		break;

		// Geldart A,B,D

	default:
		if (nNum_or > 0) {
			dVflow_or /= nNum_or;								// V_or, gas volume flow rate through a distributor, Eq. 2.2
			//dd_v_0 = 0.0001;
			dd_v_0 = 1.3 * pow(pow(dVflow_or, 2.0) / STANDARD_ACCELERATION_OF_GRAVITY, 0.2); //dv,0, theoretical diameter of a bubble, Eq. 2.1
		}
		else {	//(nNum_or <= 0)
			dd_v_0 = 0.001;
		}

		break;
	}
	// Modified Hilligardt & Werther 1986 model

	double deps_b_mean = 0; // Mean value of bubble volume fraction
	for (size_t z = 0; z < nFDHeightGridPoints; ++z)
	{
		double dh_z = vh_z[z];
		double dd_v_z = vd_v_z[z];
		if (dd_v_z == 0)
			dd_v_z = dd_v_0;

		double dPsi = CalculatePsi(dh_z, dL_char_fb, sGeldart);
		double dVflow_b_z = dPsi * (du_0 - du_mf) / (1. + dLambda);
		double du_b_z = dVflow_b_z + 0.71 * dtheta * sqrt(STANDARD_ACCELERATION_OF_GRAVITY * dd_v_z);
		double deps_b_z = dVflow_b_z / du_b_z;

		vPsi_z[z] = dPsi;
		vVflow_b_z[z] = dVflow_b_z;
		vu_b_z[z] = du_b_z;
		veps_b_z[z] = deps_b_z;

		double dd_v_z_update;
		if (z == 0)
		{
			dd_v_z_update = dd_v_0;
		}
		else
		{
			double dDelta_h_z = vh_z[z] - vh_z[z - 1];
			dd_v_z_update = Calculated_v_z_update(vd_v_z_update[z - 1], vd_v_z[z], deps_b_z, du_b_z, dDelta_h_z, z, sGeldart, nNum_or, dLambda);
			// Calculation of mean bubble volume fraction
			deps_b_mean += dDelta_h_z / vh_z[nFDHeightGridPoints - 1] * deps_b_z;
		}
		// Store new value of bubble diameter
		vd_v_z_update[z] = dd_v_z_update;

		if (vd_v_z_update[z] < 0)			//z == nFDHeightGridPoints - 1 ||
			unsigned nDebug = 1;
	}

	/* Set mean bubble volume fraction to 0 for Geldart D, according to Chen, Bachmann et al. 2017 .. nice try, but didnt work out -- solver issues
	switch (sGeldart.nType)
	{
		// Geldart D
	case 4:
		deps_b_mean = 0;
	}
	*/

	// Gas velocity in suspension
	double du_d = 0;
	switch (sGeldart.type)
	{
		// Geldart C
	case EGeldart::C:
		du_d = du_mf * pow((1 + 1.5 * deps_b_mean), 2. / 3); // Assumption: Usage of mean bubble volume fraction
		break;
		// Geldart A,B and D (TO DO: Check for correlation of Geldart D)
	default:
		if (nNum_or > 0)		// industrial distributor
		{
			du_d = (1. / 4. * (du_0 - du_mf) + du_mf); // / (1. + dLambda);
		}
		else					// porous plate distributor
		{
			du_d = (1. / 3. * (du_0 - du_mf) + du_mf);// /(1. + dLambda);
		}

		break;
	}


	double deps_d = deps_0 * pow(du_d / du_mf, 1.0 / (4.65 * (1. + dLambda)));
	//double deps_fb = 0;
	//for (size_t z = 1; z < nFDHeightGridPoints; ++z)
	//{
	//	double dleft = (1 - (1 - veps_b_z[z - 1]) * (1 - deps_d));
	//	double dright = (1 - (1 - veps_b_z[z]) * (1 - deps_d));

	//	deps_fb += (dleft + dright) / 2 * dDelta_z;		// TO DO
	//}
	//deps_fb /= dH_fb;


	/* // To be corrected IMPORTANT !!!!
	// Richardson Zaki for bed expansion (after Burgschwieger 2002)
	// temperature dependency will be critical
	double drho_g = GetDryGasDensity(_dT_ps, _dP);				// Density [kg/m^3]
	double deta_g = GetDryGasViscosity(_dT_ps, _dP);			// Dynamic viscosity [Pa s]
	double dnu_g = deta_g / drho_g;								// Kinematic viscosity [m^2 / s]

	double Re_0 = du_0 * dd_p / dnu_g;	// add GetSauter in FD model!						// Reynolds number
	double Ar = CalculateArchimedes();	// To add!!						// Archimedes number
	double Re_elu = sqrt(4 / 3 * Ar);							// Reynolds elutriation
	double Re_mf_RZ = 42.9*(1 - deps_mf)* (sqrt(1 + pow(deps_mf, 3) *Ar / pow(1 - deps_mf, 2) * 3214) - 1);

	double n_RZ = log(Re_mf_RZ / Re_elu) / log(deps_mf);

	double deps_fb =  pow((Re_0 / Re_mf_RZ), 1 / n_RZ* (1. + dLambda));
	*/

	double deps_fb;
	if (sGeldart.type == EGeldart::B || sGeldart.type == EGeldart::D)
		deps_fb = deps_d;
	else
		deps_fb = (1 - (1 - deps_b_mean) * (1 - deps_d));	// Assumption: Usage of mean bubble volume fraction

	double dH_fb_update = (1 - deps_0) / (1 - deps_fb) * dh_0;

	////////////////////////////////////////////////////////////////////////////////////
	//                              SET FUNCTION VALUES                               //
	////////////////////////////////////////////////////////////////////////////////////
	for (size_t z = 0; z < nFDHeightGridPoints; ++z)
	{
		//_pFunc[m_vnd_v_z[z]] = vd_v_z_update[z] - vd_v_z[z];
		_func[m_nd_v_z[z]] = vd_v_z_update[z];
	}
	//_pFunc[m_nh_fb] = dH_fb_update - dH_fb;
	_func[m_h_fb] = dH_fb_update;
}

void CNLModelFD::ResultsHandler(double _time, double* _vars, void* _unit)
{
	size_t nFDHeightGridPoints = GetFDHeightGridPoints();
	double dd_fb = GetBedDiameter();
	double du_0 = GetU0();
	double du_mf = GetUmf();
	SGeldart sGeldart = GetGeldart();
	double debs_Bmanuel = GetSettingsBubbleVolFractOverwrite();
	////////////////////////////////////////////////////////////////////////////////////
	//                             GET VARIABLE VALUES                                //
	////////////////////////////////////////////////////////////////////////////////////
	// Bed height
	double dH_fb = _vars[m_h_fb];
	// Bubble diameters
	std::vector<double> vd_v_z;
	for (size_t z = 0; z < nFDHeightGridPoints; ++z)
		vd_v_z.push_back(_vars[m_nd_v_z[z]]);
	// Height discretization
	std::vector<double> vh_z;
	double dDelta_z = dH_fb / (nFDHeightGridPoints - 1);
	for (size_t z = 0; z < nFDHeightGridPoints; ++z)
		vh_z.push_back(z * dDelta_z);
	////////////////////////////////////////////////////////////////////////////////////
	//                              CALCULATIONS                                      //
	////////////////////////////////////////////////////////////////////////////////////
	std::vector<double> vnu_z(nFDHeightGridPoints); //nu, ratio of gas in the bubble phase to the total gas flow, height relevant (because of psi), Eq. 3.8
	double dLambda = GetVibrationIntensity();
	for (unsigned z = 0; z < nFDHeightGridPoints; z++)
	{
		double dh_z = vh_z[z];
		double dd_v_z = vd_v_z[z];

		double dtheta = GetTheta();
		double dPsi = CalculatePsi(dh_z, dd_fb, sGeldart);
		double dVflow_b_z = dPsi * (du_0 - du_mf) / (1. + dLambda);
		double du_b_z = dVflow_b_z + 0.71 * dtheta * sqrt(STANDARD_ACCELERATION_OF_GRAVITY * dd_v_z);
		double deps_b_z = dVflow_b_z / du_b_z;
		vnu_z[z] = deps_b_z;

		// Set bubble volume fraction via user input
		if (debs_Bmanuel < 1)
		{
			vnu_z[z] = debs_Bmanuel;
		}

		// Set mean bubble volume fraction to almost 0 for Geldart D, according to Chen, Bachmann et al. 2017
		// the TD model connot handle 0, as it will calculates/returns NaNs in the property function
		/*switch (sGeldart.nType)
		{
			// Geldart D
		case 4:
			vnu_z[z] = 0.0001 ;
		}
		*/

		//double dpsi_z = CalculatePsi(dh_z, dd_fb, sGeldart); //Psi Fit for all reactor types
		//vnu_z[z] = dpsi_z*(du_0 - du_mf) / du_0; //Eq. 3.8
	}
	////////////////////////////////////////////////////////////////////////////////////
	//                              STORE RESULTS                                     //
	////////////////////////////////////////////////////////////////////////////////////
	// Check results
	SResultsNLModelFD sResults;
	if (dH_fb > 0 && dH_fb == dH_fb)
		sResults = SResultsNLModelFD(dH_fb, vh_z, vd_v_z, vnu_z);
	SetResultsNLModelFD(sResults);
}

/// Additional functions ///
double CNLModelFD::Calculated_v_z_update(double _d_v_z_1, double _d_v_z, double _eps_b_z, double _u_b_z, double _delta_h_z, size_t _heighIndex, const SGeldart& _geldart, size_t _num_or, double _lambda)
{
	double dd_v_z = 0;
	switch (_geldart.type)
	{
		// Case type "D": Hilligardt and Werther (1986)
	case EGeldart::D:
	{
		double du_0 = GetU0();		// Superficial velocity [m/s]
		double du_mf = GetUmf();	// Minimum fluidization velocity [m/s]
		double deps_mf = GetInitialBedPorosity();		// Porosity at minium fluidization velocity [-]

		// better way than Alaathar: Using the dependent bubble diameter here as well, not his randomly estimeated bubble diameter from (Alaathar Dissertation pages 30 and 31)
		//double du_b_z_single = 0.71 * sqrt(STANDARD_ACCELERATION_OF_GRAVITY * _dd_v_z);	// Velocity of single bubble [m/s]

		// Velocity in suspension phase [m/s]
		double du_d = 0;
		if (_num_or > 0)		// industrial distributor
		{
			du_d = (1. / 4. * (du_0 - du_mf) + du_mf); // / (1. + _dLambda);
		}
		else					// porous plate distributor
		{
			du_d = (1. / 3. * (du_0 - du_mf) + du_mf); // / (1. + _dLambda);
		}

		double du_b = _d_v_z; // fmax(du_b_z_single, du_0 - 2.71 * du_mf + 1.36 * du_b_z_single);	// Velocity of bubbles [m/s]


		double deps_s = deps_mf * pow((du_d / du_mf), 1. / (4.65 * (1 + _lambda))); // Porosity of suspension phase [-]

		double dalpha = du_b * deps_s / du_d;
		double dzeta = dalpha > 1 ? 0 : 1 - pow(dalpha, 3.0);
		dd_v_z = _d_v_z_1 + _delta_h_z * pow((2. * _eps_b_z / (9. * MATH_PI)), 1. / 3.) / (1 - dzeta * pow((6. * _eps_b_z / MATH_PI), 1. / 3.));


		if (du_b <= 0 || dd_v_z < 0)
			unsigned nDEBUG = 1;

		break;
	}
	// Case type "C": Zou et al. (2011)
	case EGeldart::C:
	{
		double du_0 = GetU0();		// Superficial velocity [m/s]
		double du_mf = GetUmf();	// Minimum fluidization velocity [m/s]
		double dA_fb = GetBedArea();
		size_t nNum_or = GetNumberOfOrifices();
		double dh = _delta_h_z * _heighIndex;
		//double dA_0 = 0;
		//if(nNum_or > 0)
		//	dA_0 = dA_fb / nNum_or; // Assumption: Crosssectional area of orifice in one opening of the distributor
		//if (nNum_or == 0)
		//	dA_0 = 0.0001;	 // Assumption for porous plate distributor
		double dA_0 = CalculatedA_0_GeldartC();

		dd_v_z = 0.21 * (pow((du_0 - du_mf), 0.49) * pow((dh + 4 * sqrt(dA_0)), 0.48)) / (pow(STANDARD_ACCELERATION_OF_GRAVITY, 0.2));

		if (dd_v_z < 0)
			unsigned nDEBUG = 1;

		break;

	}
	// Case type "A" and "B": Hilligardt and Werther (1987)

	default:
	{
		double dA_fb = GetBedArea();
		double du_mf = GetUmf();
		double dlambda = 280. * du_mf / STANDARD_ACCELERATION_OF_GRAVITY;		// Average life time of a bubble [s]
		dd_v_z = _d_v_z_1 + _delta_h_z * (pow((2. * _eps_b_z / 9. / MATH_PI), 1. / 3.) - _d_v_z / 3. / dlambda / _u_b_z);
		if (dd_v_z < 0)
			unsigned nDEBUG = 1;

		break;
	}
	}
	return dd_v_z;
}
double CNLModelFD::CalculatePsi(double _z, double _l_char_fb, const SGeldart& _geldart) const
{
	double dpsi = -1.; //psi, deviation from the apparent bubble volume flow rate
	double zdbed = _z / _l_char_fb; //z/d_bed

	switch (_geldart.type)
	{
		// TO DO
		//Geldart Type A
	case EGeldart::A:
		if (zdbed <= 1.)
			dpsi = 0.8;
		else
			dpsi = 0.8; //No suitable data for z/dbed. Lbed/dbed should be checked before calling this function, or use RaiseError here.
		break;
		//Geldart Type B
	case EGeldart::B:
		if (zdbed <= 1.7)
			dpsi = 0.67;
		else if ((zdbed >= 1.7) && (zdbed <= 4.))
			dpsi = 0.51 * pow(zdbed, 0.5);
		else if (zdbed > 4.)
			dpsi = 1.;
		else
		{
			return -1;
		}
		break;
		// TO DO
		//Geldart Type C
	case EGeldart::C:
		if (zdbed <= 1.)
			dpsi = 1;
		else
			dpsi = 1; // Burgschweiger 2000
		break;
		//Geldart Type D
	case EGeldart::D:
		if (zdbed < 0.55)
			dpsi = 0.26;
		else if ((zdbed >= 0.55) && (zdbed <= 8))
			dpsi = 0.35 * pow(zdbed, 0.5);
		else if (zdbed > 8)
			dpsi = 1;
		else
		{
			return -1;
		}
		break;
	default:
		break;
	}
	return dpsi;
}

double CNLModelFD::CalculatedA_0_GeldartC() const
{
	double dA_fb = GetBedArea();
	size_t nNum_or = GetNumberOfOrifices();

	double dA_0 = 0;
	if (nNum_or > 0)
		dA_0 = dA_fb / nNum_or; // Assumption: Crosssectional area of orifice in one opening of the distributor
	if (nNum_or == 0)
		dA_0 = 0.0001;	// Assumption for porous plate distributor

	return dA_0;
}

/// Declaration ///
void CNLModelFD::ClearModel()
{
	ClearVariables();

	Clear1D(m_nd_v_z);
}

//////////////////////////////////////////////////////////////////////////
/// Thermodynamics - Model

/// Destructor ///
CNLModelTD::~CNLModelTD()
{
	ClearModel();
}

/// Declaration / Initialization ///
void CNLModelTD::AddVariables(double _time)
{
	ClearModel();

	// Get attributes
	size_t nSizeClasses = GetSizeClasses();
	size_t nMoistureClasses = GetMoistureClasses();
	size_t nRTGridPoints = GetRTGridPoints();
	size_t nHeightGridPoints = GetHeightGridPoints();


	// Boundary conditions
	m_HeightSubGridPoints = 0;
	m_RTSubGridPoints = 0;
	SBoundaryCondition sBC = CalculateBoundaryConditions(_time);

	//Add and set X, h, T in the particle phase as NL Variable
	m_iX_p_dfi.resize(nSizeClasses);
	m_ih_p_dfi.resize(nSizeClasses);

	SBoundaryCondition sInit = CalculateInitialConditions(sBC);

	for (size_t d = 0; d < nSizeClasses; d++)
	{
		m_iX_p_dfi[d].resize(nMoistureClasses);
		m_ih_p_dfi[d].resize(nMoistureClasses);
		for (size_t f = 0; f < nMoistureClasses; f++)
		{
			m_iX_p_dfi[d][f].resize(nRTGridPoints);
			m_ih_p_dfi[d][f].resize(nRTGridPoints);
			for (size_t i = 0; i < nRTGridPoints; i++)
			{
				//Calculate and set the initial value of the parameters to be solved by the solver
				m_iX_p_dfi[d][f][i] = AddNLVariable(sInit.X_p[f], 0.0, 100, 100);
				m_ih_p_dfi[d][f][i] = AddNLVariable(sInit.h_p[f], 0.0, 1e-5, 1e-5);
			}
		}
	}

	//Add Y, h and T of the suspension and bubble phase as NL Variable
	m_iY_b_z.resize(nHeightGridPoints);
	m_ih_b_z.resize(nHeightGridPoints);
	m_iY_s_z.resize(nHeightGridPoints);
	m_ih_s_z.resize(nHeightGridPoints);

	for (size_t z = 0; z < (nHeightGridPoints); z++)
	{
		m_ih_s_z[z] = AddNLVariable(sInit.h_g, 0.0, 100, 100);
		m_iY_s_z[z] = AddNLVariable(sInit.Y_g, 0.0, 1e-5, 1e-5);
		m_ih_b_z[z] = AddNLVariable(sInit.h_g, 0.0, 100, 100);
		m_iY_b_z[z] = AddNLVariable(sInit.Y_g, 0.0, 1e-5, 1e-5);
	}
}

void CNLModelTD::Initialize(double _time, const STDParameters& _tdParameters, const SResultsNLModelFD& _resultsNLModelFD)
{
	// TD Parameters
	SetTDParameters(_tdParameters);
	// Boundary Conditions
	SBoundaryCondition sBC = CalculateBoundaryConditions(_time);
	SetBoundaryConditions(sBC);
	// TD Model parameters
	double dh_fb = _resultsNLModelFD.h_fb;
	std::vector<std::vector<double>> vmfrac_df = CalculateMassFractionVector(_time);

	std::vector<double> vApAbed_d = CalculateRatioApAbed(_time, dh_fb, vmfrac_df);
	double ddp_sauter_in = GetInParticleSauter(_time);
	double deps_fb = CaclulateBedPorosity(_time, dh_fb);
	double dC_Ntaumax = GetSettingsNtaumax();
	double dtau_max = CalculateTauMax(dC_Ntaumax);

	size_t nHeightGridPoints = GetHeightGridPoints();
	// Height and dimensionless height vectors
	std::vector<double> vh_z(nHeightGridPoints);
	std::vector<double> vzeta_z(nHeightGridPoints);
	for (unsigned z = 0; z < nHeightGridPoints; ++z)
	{
		vh_z[z] = z / ((double)nHeightGridPoints - 1) * dh_fb;
		vzeta_z[z] = z / ((double)nHeightGridPoints - 1);
	}
	// Assumption: calculate mean value of gas distribution over the bed height
	double dnu_mean = 0;
	for (unsigned z = 1; z < nHeightGridPoints; ++z)
	{
		dnu_mean += 0.5 * (_resultsNLModelFD.nu_z[z] + _resultsNLModelFD.nu_z[z - 1]) * (vzeta_z[z] - vzeta_z[z - 1]);
	}
	std::vector<double> vtau_i = CalculateRTTauVector(dtau_max);
	std::vector<double> vn_i = CalculateRTDensityVector(dtau_max);
	std::vector<double> vd_d = GetMeanDiameterVector();
	std::vector<std::vector<double>> vAps_df = CalculateParticleSurfaceVector(_time, vd_d);
	SetTDModelParameters(STDModelParameters(vApAbed_d, ddp_sauter_in, deps_fb, dtau_max, dh_fb, vh_z, vzeta_z, dnu_mean, vtau_i, vn_i, vmfrac_df, vd_d, vAps_df));
}


/// Solver functions ///
void CNLModelTD::CalculateFunctions(double* _vars, double* _func, void* _unit)
{
	////////////////////////////////////////////////////////////////////////////////////
	//                         CONSTANT MODEL ATTRIBUTES                              //
	////////////////////////////////////////////////////////////////////////////////////
	// Discretization
	size_t nRTGridPoints = GetRTGridPoints();
	size_t nHeightGridPoints = GetHeightGridPoints();
	size_t nSizeClasses = GetSizeClasses();
	size_t nMoistureClasses = GetMoistureClasses();

	////////////////////////////////////////////////////////////////////////////////////
	//                              VECTOR DECLARATION                                //
	////////////////////////////////////////////////////////////////////////////////////
	std::vector<double> vQ_sb_z;  // 1D-vector for storage of heat flow from suspension to bubble phase at height class z
	std::vector<double> vM_sb_z;  // 1D-vector for storage of mass flow from suspension to bubble phase at height class z
	std::vector<double> vH_sb_z;  // 1D-vector for storage of enthalpy flow from suspension to bubble phase at height class z
	std::vector<double> vQ_ps_z;  // 1D-vector for storage of total heat flow from all particles to suspension phase at height class z
	std::vector<double> vM_ps_z;  // 1D-vector for storage of total mass flow from all particles to suspension phase at height class z
	std::vector<double> vH_ps_z;  // 1D-vector for storage of total enthalpy flow from all particles to suspension phase at height class z
	std::vector<std::vector<std::vector<double>>> vvvQ_ps_dfi;  // 3D-vector for storage of total heat flow from particles of size class d, initial moisture class f and residence class i to suspension phase
	std::vector<std::vector<std::vector<double>>> vvvM_ps_dfi;  // 3D-vector for storage of total mass flow from particles of size class d, initial moisture class f and residence class i to suspension phase
	std::vector<std::vector<std::vector<double>>> vvvH_ps_dfi;  // 3D-vector for storage of total enthalpy flow from particles of size class d, initial moisture class f and residence class i to suspension phase
	std::vector<std::vector<std::vector<double>>> vvvQ_pp_dfi;	// 3D-vector for storage of heat flow between particles of size class d, initial moisture class f and residence class i
	std::vector<std::vector<std::vector<double>>> vvvQ_pw_dfi;	// 3D-vector for storage of heat flow between particles and wall of size class d, initial moisture class f and residence class i

	std::vector<double> vT_s_z;	// 1D-vector for storage of temperatures of suspension phase at height class z
	std::vector<double> vT_b_z;	// 1D-vector for storage of temperatures of bubble phase at height class z
	std::vector<std::vector<std::vector<double>>> vvvT_p_dfi;	// 3D-vector for storage of temperatures of particle phase at size class d, initial moisture class f and residence class i

	////////////////////////////////////////////////////////////////////////////////////
	//                              GET VARIABLE VALUES                               //
	////////////////////////////////////////////////////////////////////////////////////
	//	Particle phase
	GetValues3D(_vars, m_X_p_dfi, m_iX_p_dfi);
	GetValues3D(_vars, m_h_p_dfi, m_ih_p_dfi);
	CalculateParticleTemperatures(m_h_p_dfi, m_X_p_dfi, vvvT_p_dfi);
	//	Suspension phase
	GetValues1D(_vars, m_Y_s_z, m_iY_s_z);
	GetValues1D(_vars, m_h_s_z, m_ih_s_z);
	CalculateGasTemperatures(m_h_s_z, m_Y_s_z, vT_s_z);
	//	Bubble phase
	GetValues1D(_vars, m_Y_b_z, m_iY_b_z);
	GetValues1D(_vars, m_h_b_z, m_ih_b_z);
	CalculateGasTemperatures(m_h_b_z, m_Y_b_z, vT_b_z);

	////////////////////////////////////////////////////////////////////////////////////
	//                              UPDATE FUNCTION VALUES                            //
	////////////////////////////////////////////////////////////////////////////////////
	std::vector<double> vn_i = GetRTNumberDensityVector();
	std::vector<double> vtau_i = GetRTTauVector();
	std::vector<double> vzeta_z = GetDimensionlessHeightVector();

	/// Heat and mass transfer ///
	//  Particle -> Suspension (ps)
	ExchangeParticleSuspension(vzeta_z, vtau_i, vn_i, m_Y_s_z, vT_s_z, m_X_p_dfi, vvvT_p_dfi, vM_ps_z, vH_ps_z, vQ_ps_z, vvvM_ps_dfi, vvvH_ps_dfi, vvvQ_ps_dfi);
	//	Suspension -> Bubble (sb)
	ExchangeSuspensionBubble(m_Y_s_z, vT_s_z, m_Y_b_z, vT_b_z, vM_sb_z, vH_sb_z, vQ_sb_z);
	//	Particle -> Particle (pp)
	ExchangeParticleParticle(vtau_i, vn_i, m_X_p_dfi, vvvT_p_dfi, vvvQ_pp_dfi);
	//	Particle -> Particle (pp)
	ExchangeParticleWall(vtau_i, vn_i, m_X_p_dfi, vvvT_p_dfi, vvvQ_pw_dfi);

	/// Balance equations ///
	//	Particle phase
	BalanceParticlePhase(vtau_i, m_X_p_dfi, m_h_p_dfi, vvvT_p_dfi, m_X_p_dfi_update, m_h_p_dfi_update, vvvM_ps_dfi, vvvQ_ps_dfi, vvvH_ps_dfi, vvvQ_pp_dfi, vvvQ_pw_dfi);
	//	Suspension phase
	BalanceSuspensionPhase(vzeta_z, m_Y_s_z, m_h_s_z, vT_s_z, m_Y_s_z_update, m_h_s_z_update, vM_ps_z, vH_ps_z, vQ_ps_z, vM_sb_z, vH_sb_z, vQ_sb_z);
	//	Bubble phase
	BalanceBubblePhase(vzeta_z, m_Y_b_z, m_h_b_z, vT_b_z, m_Y_b_z_update, m_h_b_z_update, vM_sb_z, vH_sb_z, vQ_sb_z);

	/// Set function values ///
	//	Particle phase
	SetValues3D(_func, m_X_p_dfi_update, m_iX_p_dfi);
	SetValues3D(_func, m_h_p_dfi_update, m_ih_p_dfi);
	//	Suspension phase
	SetValues1D(_func, m_Y_s_z_update, m_iY_s_z);
	SetValues1D(_func, m_h_s_z_update, m_ih_s_z);
	//	Bubble phase
	SetValues1D(_func, m_Y_b_z_update, m_iY_b_z);
	SetValues1D(_func, m_h_b_z_update, m_ih_b_z);
}

void CNLModelTD::ResultsHandler(double _time, double* _vars, void* _unit)
{
	/// Discretization ///
	size_t nRTGridPoints = GetRTGridPoints();
	size_t nHeightGridPoints = GetHeightGridPoints();
	size_t nSizeClasses = GetSizeClasses();
	size_t nMoistureClasses = GetMoistureClasses();
	/// Vector declaration ///
	// Particle phase (p)
	std::vector<std::vector<std::vector<double>>> vvvX_p_dfi, vvvh_p_dfi, vvvT_p_dfi;			// 3D-vector for storage of moisture, enthalpy and temperature of particles of size class d, moisture class f and residence class i
	// Suspension phase (s)
	std::vector<double> vY_s_z, vh_s_z, vT_s_z;			// 1D-vector for storage of moisture content, enthalpy and temperature of suspension phase at height class z
	// Bubble phase (b)
	std::vector<double> vY_b_z, vh_b_z, vT_b_z;			// 1D-vector for storage of moisture content, enthalpy and temperature of bubble phase at height class z
	// Heat flow to wall
	std::vector<std::vector<std::vector<double>>> vvvQ_pw_dfi;	// 3D-vector for storage of heat flow between particles and wall of size class d, initial moisture class f and residence class i

	/// Variable vectors
	//	Particle phase (p)
	GetValues3D(_vars, vvvX_p_dfi, m_iX_p_dfi);
	GetValues3D(_vars, vvvh_p_dfi, m_ih_p_dfi);
	CalculateParticleTemperatures(vvvh_p_dfi, vvvX_p_dfi, vvvT_p_dfi);
	//	Suspension phase (s)
	GetValues1D(_vars, vY_s_z, m_iY_s_z);
	GetValues1D(_vars, vh_s_z, m_ih_s_z);
	CalculateGasTemperatures(vh_s_z, vY_s_z, vT_s_z);
	//	Bubble phase (b)
	GetValues1D(_vars, vY_b_z, m_iY_b_z);
	GetValues1D(_vars, vh_b_z, m_ih_b_z);
	CalculateGasTemperatures(vh_b_z, vY_b_z, vT_b_z);

	//	Particle -> Particle (pp)
	std::vector<double> vn_i = GetRTNumberDensityVector();
	std::vector<double> vtau_i = GetRTTauVector();
	ExchangeParticleWall(vtau_i, vn_i, m_X_p_dfi, vvvT_p_dfi, vvvQ_pw_dfi);
	double dQ_pw = 0;
	for (size_t d = 0; d < nSizeClasses; ++d)
	{
		for (size_t f = 0; f < nMoistureClasses; ++f)
		{
			for (size_t i = 1; i < nRTGridPoints; ++i)
			{
				dQ_pw += 0.5 * (vn_i[i] * vvvQ_pw_dfi[d][f][i] + vn_i[i - 1] * vvvQ_pw_dfi[d][f][i - 1]) * (vtau_i[i] - vtau_i[i - 1]);
			}
		}
	}
	/// Water mass balance ///
	// Mass flows
	// Dry air flow
	double dmflow_g = GetInGasMYTP(_time).m;	// Mass flow of dry air [kg/s]
	double dnu_0 = GetTDGasSplitFactor();		// Split factor of gas flow between suspension and bubble phase [-]
	// Dry particle flow
	double dC_mflow_p = 1;		// obsolete residual from earlier version, not used anymore
	double dmflow_p = dC_mflow_p * GetInParticleMYTP(_time).m;		// Mass flow of dry particles [kg/s]
	std::vector<std::vector<double>> vvmfrac_df = GetMassFractionVector();	// Mass fractions of size and moisture classes [-]
	// Inlet water mass flows
	// Gas phase
	double dY_g_in = GetBoundaryConditions().Y_g;
	double dmflow_H2O_g_in = dY_g_in * dmflow_g;
	// Particle phase
	std::vector<double> vX_p_in = GetBoundaryConditions().X_p;
	double dmflow_H2O_p_in = 0;
	for (size_t d = 0; d < nSizeClasses; ++d)
	{
		for (size_t f = 0; f < nMoistureClasses; ++f)
		{
			dmflow_H2O_p_in += vX_p_in[f] * vvmfrac_df[d][f] * dmflow_p;
		}
	}
	// Overall inlet mass flow
	double dmflow_H2O_in = dmflow_H2O_p_in + dmflow_H2O_g_in;
	// Outlet water mass flows
	// Gas phase
	double dmflow_H2O_b_out = vY_b_z[nHeightGridPoints - 1] * dnu_0 * dmflow_g;
	double dmflow_H2O_s_out = vY_s_z[nHeightGridPoints - 1] * (1 - dnu_0) * dmflow_g;
	double dmflow_H2O_g_out = dmflow_H2O_b_out + dmflow_H2O_s_out;
	// Particle phase
	double dX_p_out = CalculateMoistureMeanValue(vvvX_p_dfi, 0);
	double dT_p_out = CalculateTemperatureMeanValue(vvvT_p_dfi, 0);

	double dmflow_H2O_p_out = dX_p_out * dmflow_p;
	// Overall outlet mass flow
	double dmflow_H2O_out = dmflow_H2O_p_out + dmflow_H2O_g_out;
	// Balance equation
	double dmflow_H2O_diff = dmflow_H2O_out - dmflow_H2O_in;
	// Relative deviation
	double dBalanceH2ODeviation = dmflow_H2O_diff / dmflow_H2O_in;

	/// Enthalpy balance ///
	// Inlet enthalpy flows
	// Gas phase
	double dh_g = GetBoundaryConditions().h_g;
	double dhflow_g_in = dmflow_g * dh_g;
	// Particle phase
	double dhflow_p_in = 0;
	std::vector<double> vh_p_in = GetBoundaryConditions().h_p;
	for (size_t d = 0; d < nSizeClasses; ++d)
	{
		for (size_t f = 0; f < nMoistureClasses; ++f)
		{
			dhflow_p_in += vh_p_in[f] * vvmfrac_df[d][f] * dmflow_p;
		}
	}
	// Overall inlet enthalpy flow
	double dhflow_in = dhflow_g_in + dhflow_p_in;
	// Outlet enthalpy flows
	// Gas phase
	double dhflow_b_out = dnu_0 * dmflow_g * vh_b_z[nHeightGridPoints - 1];
	double dhflow_s_out = (1 - dnu_0) * dmflow_g * vh_s_z[nHeightGridPoints - 1];
	double dhflow_g_out = dhflow_b_out + dhflow_s_out;
	// Particle flow
	double dh_p_out = CalculateTemperatureMeanValue(vvvh_p_dfi, 0);
	double dhflow_p_out = dh_p_out * dmflow_p;
	// Overall outlet enthalpy flow
	double dhflow_out = dhflow_b_out + dhflow_s_out + dhflow_p_out + dQ_pw;
	// Balance equation
	double dhflow_diff = dhflow_out - dhflow_in;
	// Relative deviation
	double dBalanceEnthalpyDeviation = dhflow_diff / dhflow_in;

	// Heat transfer coefficient
	std::vector<std::vector<double>> vvalpha_ps_dz(nSizeClasses, std::vector<double>(nHeightGridPoints, 0));
	double dP = GetTDInGasMYTP().P;													// Gas pressure at inlet [Pa]
	std::vector<double> vApAbed_d = GetRatioApAbedVector();							// 1D-vector with mean particle diameters of size classes d [m]
	std::vector<double> vdp_d = GetTDMeanDiameterVector();							// 1D-vector with mean particle diameters of size classes d [m]
	std::vector<std::vector<double>> vvw_df = GetMassFractionVector();				// 2D-vector with mass fractions of size classes d and moisture classes f [kg/kg]

	std::vector<std::vector<std::vector<std::vector<double>>>> vvvvalpha_ps_dfiz;
	Resize4D(vvvvalpha_ps_dfiz, nSizeClasses, nMoistureClasses, nRTGridPoints, nHeightGridPoints);	// 4D-vector for storage of heat transfer coefficients  from particles of size class d and residence class i to suspension gas at height class z

	ParallelFor(nHeightGridPoints, [&](size_t z)
		{
			for (size_t d = 0; d < nSizeClasses; d++)
			{
				for (size_t f = 0; f < nMoistureClasses; f++)
				{
					for (size_t i = 0; i < nRTGridPoints; i++)
					{
						double dT_p_dfi = vvvT_p_dfi[d][f][i];		// Temperature of particles of size class d and moisture class f at residence time grid point i [K]
						double dT_s_z = vT_s_z[z];		// Temperature of suspension gas at height grid point z [K]
						double dp_d = vdp_d[d];						// Particle diameter of size class d [m]
						double dT_ps = 0.5 * (dT_p_dfi + dT_s_z);		// Mean temperature between suspension and particle phase

						double dApAbed_d = vApAbed_d[d];				// Ratio of particle surface to cross-sectional area of dryer or size class d [-]
						double dalpha_ps = CalculateAlpha_PS(dT_ps, dP, dp_d, dApAbed_d);		// Heat transfer coefficient [W/(m^2 K)]

						// Heat transfer coefficient
						vvvvalpha_ps_dfiz[d][f][i][z] = dalpha_ps;
					}
				}
			}
		});
	ParallelFor(nHeightGridPoints, [&](size_t z)
		{
			for (size_t d = 0; d < nSizeClasses; d++)
			{
				for (size_t f = 0; f < nMoistureClasses; f++)
				{
					for (size_t i = 1; i < nRTGridPoints; i++)
					{
						double dtau_i = vtau_i[i];
						double dtau_i_1 = vtau_i[i - 1];
						double dn_i = vn_i[i];
						double dn_i_1 = vn_i[i - 1];
						vvalpha_ps_dz[d][z] += 0.5 * (vvvvalpha_ps_dfiz[d][f][i][z] + vvvvalpha_ps_dfiz[d][f][i - 1][z]) * 0.5 * (dn_i + dn_i_1) * (dtau_i - dtau_i_1);
					}
					vvalpha_ps_dz[d][z] *= vvw_df[d][f];
				}
			}
		});


	/// Storage of results ///
	double dT_p_mean = CalculateTemperatureMeanValue(vvvT_p_dfi, 0);
	double dh_g_out = dhflow_g_out / dmflow_g;
	double dY_g_out = dnu_0 * vY_b_z[nHeightGridPoints - 1] + (1 - dnu_0) * vY_s_z[nHeightGridPoints - 1];
	double dT_g_est = dnu_0 * vT_b_z[nHeightGridPoints - 1] + (1 - dnu_0) * vT_s_z[nHeightGridPoints - 1];
	double dT_g_out = GetMoistGasTemperature(dh_g_out, dY_g_out, dT_g_est, dP);
	SResultsNLModelTD sResults = SResultsNLModelTD(vvvX_p_dfi, vvvh_p_dfi, vvvT_p_dfi, dT_p_mean, vY_s_z, vh_s_z, vT_s_z, vY_b_z, vh_b_z, vT_b_z, dT_g_out, dX_p_out, dT_p_out, vvalpha_ps_dz, dQ_pw, dBalanceH2ODeviation, dBalanceEnthalpyDeviation);
	SetTDModelResults(sResults);
}


/// Exchange flow functions ///
void CNLModelTD::ExchangeParticleSuspension(const std::vector<double>& _zeta_z, const std::vector<double>& _tau_i, const std::vector<double>& _n_i, const std::vector<double>& _y_s_z, const std::vector<double>& _t_s_z, const std::vector<std::vector<std::vector<double>>>& _x_p_dfi, const std::vector<std::vector<std::vector<double>>>& _t_p_dfi, std::vector<double>& _m_ps_z, std::vector<double>& _h_ps_z, std::vector<double>& _q_ps_z, std::vector<std::vector<std::vector<double>>>& _m_ps_dfi, std::vector<std::vector<std::vector<double>>>& _h_ps_dfi, std::vector<std::vector<std::vector<double>>>& _q_ps_dfi)
{
	/// Discretization ///
	size_t nRTGridPoints = GetRTGridPoints();
	size_t nHeightGridPoints = GetHeightGridPoints();
	size_t nSizeClasses = GetSizeClasses();
	size_t nMoistureClasses = GetMoistureClasses();

	/// Declaration of temporary vectors ///
	// Transfer flows for each size, moisture, residence time and height class
	std::vector<std::vector<std::vector<std::vector<double>>>> vvvvQ_ps_dfiz;
	std::vector<std::vector<std::vector<std::vector<double>>>> vvvvM_ps_dfiz;
	std::vector<std::vector<std::vector<std::vector<double>>>> vvvvH_ps_dfiz;

	/// Resize exchange flow vectors ///
	_q_ps_z.resize(nHeightGridPoints);  // 1D-vector for storage of total heat flow from all particles to suspension phase at height class z
	_m_ps_z.resize(nHeightGridPoints);  // 1D-vector for storage of total mass flow from all particles to suspension phase at height class z
	_h_ps_z.resize(nHeightGridPoints);  // 1D-vector for storage of total enthalpy flow from all particles to suspension phase at height class z
	Resize3D(_q_ps_dfi, nSizeClasses, nMoistureClasses, nRTGridPoints);  // 3D-vector for storage of total heat flow from particles of size class d and residence class i to suspension phase
	Resize3D(_m_ps_dfi, nSizeClasses, nMoistureClasses, nRTGridPoints);  // 3D-vector for storage of total mass flow from particles of size class d and residence class i to suspension phase
	Resize3D(_h_ps_dfi, nSizeClasses, nMoistureClasses, nRTGridPoints);  // 3D-vector for storage of total enthalpy flow from particles of size class d and residence class i to suspension phase

	Resize4D(vvvvQ_ps_dfiz, nSizeClasses, nMoistureClasses, nRTGridPoints, nHeightGridPoints);	// 4D-vector for storage of heat flow from particles of size class d and residence class i to suspension gas at height class z
	Resize4D(vvvvM_ps_dfiz, nSizeClasses, nMoistureClasses, nRTGridPoints, nHeightGridPoints);	// 4D-vector for storage of mass flow from particles of size class d and residence class i to suspension gas at height class z
	Resize4D(vvvvH_ps_dfiz, nSizeClasses, nMoistureClasses, nRTGridPoints, nHeightGridPoints);	// 4D-vector for storage of enthalpy flow from particles of size class d and residence class i to suspension gas at height class z

	/// Model attributes ///
	std::vector<double> vdp_d = GetTDMeanDiameterVector();							// 1D-vector with mean particle diameters of size classes d [m]
	std::vector<double> vApAbed_d = GetRatioApAbedVector();							// 1D-vector with mean particle diameters of size classes d [m]

	std::vector<std::vector<double>> vvAps_df = GetParticleSurfaceVector();			// 2D-vector with total particle surfaces of size classes d and moisture classes f [m^2]

	double dP = GetTDInGasMYTP().P;													// Gas pressure at inlet [Pa]
	double dphi_eq = GetPhiEq();													// Relative humidity in equilibrium state of the particle. [-]

	bool bSwitchHT = GetSettingsSwitchHT();				// Variable to switch of heat transfer
	bool bSwitchMT = GetSettingsSwitchMT();				// Variable to switch of mass transfer

	// Transfer flows for each size, residence time and height class
	for (size_t d = 0; d < nSizeClasses; d++)
	{
		for (size_t f = 0; f < nMoistureClasses; f++)
		{
			ParallelFor(nRTGridPoints, [&](size_t i)
			{
				for (size_t z = 0; z < nHeightGridPoints; z++)
				{
					// Get state of particle and suspension phase
					double dT_p_dfi = _t_p_dfi[d][f][i];		// Temperature of particles of size class d and moisture class f at residence time grid point i [K]
					double dX_p_dfi = _x_p_dfi[d][f][i];		// Moisture content of particles of size class d and moisture class f at residence time grid point i [kg/kg]

					double dT_s_z = _t_s_z[z];		// Temperature of suspension gas at height grid point z [K]
					double dY_s_z = _y_s_z[z];		// Moisture content of suspension gas at height grid point z [kg/kg]

					double dn_i = _n_i[i];						// Number density of particles at residence time grid point i [1/s]
					double dA_ps_df = vvAps_df[d][f];			// Surface of particles of size class d and moisture class f [m^2]
					double dp_d = vdp_d[d];						// Particle diameter of size class d [m]
					double dApAbed_d = vApAbed_d[d];				// Ratio of particle surface to cross-sectional area of dryer or size class d [-]

					double dT_ps = 0.5 * (dT_p_dfi + dT_s_z);		// Mean temperature between suspension and particle phase

					double dalpha_ps = CalculateAlpha_PS(dT_ps, dP, dp_d, dApAbed_d);		// Heat transfer coefficient [W/(m^2 K)]
					double dphi_s = GetRelativeHumidity(dY_s_z, dT_s_z, dP);					// Relative humidity [-]
					double dXeq = GetParticleEquilibriumMoistureContent(dphi_s, dT_p_dfi);		// Particle equilibrium moisture content [kg/kg]

					double dDryingCurve = CalculateNormalizedDryingCurve(dX_p_dfi, dXeq);		// Normalized drying curve [-]
					double drho_air_ps = GetDryGasDensity(dT_ps, dP);							// Dry air density [kg/m^3]
					double dYeq = GetGasEquilibriumMoistureContent(dT_p_dfi, dphi_eq, dP);		// Equilibrium moisture content of the gas [kg/kg]

					double dbeta_ps = CalculateBeta_PS(dT_ps, dY_s_z, dP, dp_d, dApAbed_d);		// Mass transfer coefficient [m/s]
					double dh_H2O_p = GetWaterVaporEnthalpy(dT_p_dfi, dP);						// Specific water vapor enthalpy [J/kg]

					// Calculate exchange streams
					double dM_ps_dfiz = dn_i * dDryingCurve * drho_air_ps * dbeta_ps * dA_ps_df * (dYeq - dY_s_z);		// Water vapor mass flow [kg/s]
					dM_ps_dfiz = dM_ps_dfiz > 0 ? dM_ps_dfiz : 0;														// Water vapor mass flow has to be positive: Going from particle to suspension phase
					dM_ps_dfiz = bSwitchMT ? dM_ps_dfiz : 0;
					double dQ_ps_dfiz = dn_i * dA_ps_df * dalpha_ps * (dT_p_dfi - dT_s_z);								// Heat flow due to heat exchange [J/s]
					dQ_ps_dfiz = bSwitchHT ? dQ_ps_dfiz : 0;
					double dH_ps_dfiz = dM_ps_dfiz * dh_H2O_p;															// Enthalpy flow due to mass transfer [J/s]

					// Set calculated exchange streams to storage vectors
					vvvvQ_ps_dfiz[d][f][i][z] = dQ_ps_dfiz;
					vvvvM_ps_dfiz[d][f][i][z] = dM_ps_dfiz;
					vvvvH_ps_dfiz[d][f][i][z] = dH_ps_dfiz;
				}
			});
		}
	}
	// Integration
	// Summation over all height classes for each size and residence time class
	for (size_t d = 0; d < nSizeClasses; d++)
	{
		for (size_t f = 0; f < nMoistureClasses; f++)
		{
			ParallelFor(nRTGridPoints, [&](size_t i)
				{
					if (i == 0)
						return;

					double dSumM_ps_dfi = 0;
					double dSumH_ps_dfi = 0;
					double dSumQ_ps_dfi = 0;
					double dSumM_ps_dfi_1 = 0;
					double dSumH_ps_dfi_1 = 0;
					double dSumQ_ps_dfi_1 = 0;
					for (size_t z = 1; z < nHeightGridPoints; z++)
					{
						// TO DO Possible optimization/Parallelization
						double dzeta_z = _zeta_z[z];
						double dzeta_z_1 = _zeta_z[z - 1];
						dSumM_ps_dfi += (vvvvM_ps_dfiz[d][f][i][z] + vvvvM_ps_dfiz[d][f][i][z - 1]) * 0.5 * (dzeta_z - dzeta_z_1);
						dSumH_ps_dfi += (vvvvH_ps_dfiz[d][f][i][z] + vvvvH_ps_dfiz[d][f][i][z - 1]) * 0.5 * (dzeta_z - dzeta_z_1);
						dSumQ_ps_dfi += (vvvvQ_ps_dfiz[d][f][i][z] + vvvvQ_ps_dfiz[d][f][i][z - 1]) * 0.5 * (dzeta_z - dzeta_z_1);
						dSumM_ps_dfi_1 += (vvvvM_ps_dfiz[d][f][i - 1][z] + vvvvM_ps_dfiz[d][f][i - 1][z - 1]) * 0.5 * (dzeta_z - dzeta_z_1);
						dSumH_ps_dfi_1 += (vvvvH_ps_dfiz[d][f][i - 1][z] + vvvvH_ps_dfiz[d][f][i - 1][z - 1]) * 0.5 * (dzeta_z - dzeta_z_1);
						dSumQ_ps_dfi_1 += (vvvvQ_ps_dfiz[d][f][i - 1][z] + vvvvQ_ps_dfiz[d][f][i - 1][z - 1]) * 0.5 * (dzeta_z - dzeta_z_1);
					}
					double dn_i = _n_i[i];
					double dn_i_1 = _n_i[i - 1];
					_m_ps_dfi[d][f][i] = (dSumM_ps_dfi / dn_i + dSumM_ps_dfi_1 / dn_i_1) * 0.5;
					_h_ps_dfi[d][f][i] = (dSumH_ps_dfi / dn_i + dSumH_ps_dfi_1 / dn_i_1) * 0.5;
					_q_ps_dfi[d][f][i] = (dSumQ_ps_dfi / dn_i + dSumQ_ps_dfi_1 / dn_i_1) * 0.5;
				});
		}
	}
	// Summation over all size and residence time classes for each height class
	ParallelFor(nHeightGridPoints, [&](size_t z)
	{
		if (z == 0)
			return;

		double dSumM_ps_z = 0;
		double dSumH_ps_z = 0;
		double dSumQ_ps_z = 0;
		double dSumM_ps_z_1 = 0;
		double dSumH_ps_z_1 = 0;
		double dSumQ_ps_z_1 = 0;
		for (size_t d = 0; d < nSizeClasses; d++)
		{
			for (size_t f = 0; f < nMoistureClasses; f++)
			{
				for (size_t i = 1; i < nRTGridPoints; i++)
				{
					double dtau_i = _tau_i[i];
					double dtau_i_1 = _tau_i[i - 1];
					dSumM_ps_z += (vvvvM_ps_dfiz[d][f][i][z] + vvvvM_ps_dfiz[d][f][i - 1][z]) * 0.5 * (dtau_i - dtau_i_1);
					dSumH_ps_z += (vvvvH_ps_dfiz[d][f][i][z] + vvvvH_ps_dfiz[d][f][i - 1][z]) * 0.5 * (dtau_i - dtau_i_1);
					dSumQ_ps_z += (vvvvQ_ps_dfiz[d][f][i][z] + vvvvQ_ps_dfiz[d][f][i - 1][z]) * 0.5 * (dtau_i - dtau_i_1);
					dSumM_ps_z_1 += (vvvvM_ps_dfiz[d][f][i][z - 1] + vvvvM_ps_dfiz[d][f][i - 1][z - 1]) * 0.5 * (dtau_i - dtau_i_1);
					dSumH_ps_z_1 += (vvvvH_ps_dfiz[d][f][i][z - 1] + vvvvH_ps_dfiz[d][f][i - 1][z - 1]) * 0.5 * (dtau_i - dtau_i_1);
					dSumQ_ps_z_1 += (vvvvQ_ps_dfiz[d][f][i][z - 1] + vvvvQ_ps_dfiz[d][f][i - 1][z - 1]) * 0.5 * (dtau_i - dtau_i_1);
				}
			}
		}
		_m_ps_z[z] = (dSumM_ps_z + dSumM_ps_z_1) * 0.5;
		_h_ps_z[z] = (dSumH_ps_z + dSumH_ps_z_1) * 0.5;
		_q_ps_z[z] = (dSumQ_ps_z + dSumQ_ps_z_1) * 0.5;
	});
}

void CNLModelTD::ExchangeSuspensionBubble(const std::vector<double>& _y_s_z, const std::vector<double>& _t_s_z, const std::vector<double>& _y_b_z, const std::vector<double>& _t_b_z, std::vector<double>& _m_sb_z, std::vector<double>& _h_sb_z, std::vector<double>& _q_sb_z)
{
	/// Discretization ///
	size_t nRTGridPoints = GetRTGridPoints();
	size_t nHeightGridPoints = GetHeightGridPoints();

	/// Resize exchange flow vectors ///
	_q_sb_z.resize(nHeightGridPoints);  // 1D-vector for storage of heat flow from suspension to bubble phase at height class z
	_m_sb_z.resize(nHeightGridPoints);  // 1D-vector for storage of mass flow from suspension to bubble phase at height class z
	_h_sb_z.resize(nHeightGridPoints);  // 1D-vector for storage of enthalpy flow from suspension to bubble phase at height class z

	/// Model attributes ///
	double dmflow_g = GetTDInGasMYTP().m;	// Gas flow at inlet [kg/s]
	double dP = GetTDInGasMYTP().P;			// Gas pressure at inlet [Pa]
	double dh_fb = GetBedHeight();			// Height of fluidized bed [m]

	// Transfer flows for each height grid point
	std::vector<double> vtempQ_sb_z(nHeightGridPoints, 0);  // 1D-vector for temporary storage of heat flow from suspension to bubble phase at height class z
	std::vector<double> vtempM_sb_z(nHeightGridPoints, 0);  // 1D-vector for temporary storage of mass flow from suspension to bubble phase at height class z
	std::vector<double> vtempH_sb_z(nHeightGridPoints, 0);  // 1D-vector for temporary storage of enthalpy flow from suspension to bubble phase at height class z

	bool bSwitchHT = GetSettingsSwitchHT();				// Variable to switch of heat transfer
	bool bSwitchMT = GetSettingsSwitchMT();				// Variable to switch of mass transfer

	ParallelFor(nHeightGridPoints, [&](size_t z)
	{
		// Get state of particle and suspension phase
		double dT_s_z = _t_s_z[z];
		double dY_s_z = _y_s_z[z];
		double dT_b_z = _t_b_z[z];
		double dY_b_z = _y_b_z[z];

		// Mean temperature
		double dT_sb = 0.5 * (dT_s_z + dT_b_z);
		double dY_sb = 0.5 * (dY_s_z + dY_b_z);

		// TO DO comments
		// TO DO moist air flow ?
		double drho_air = GetMoistGasDensity(dT_sb, dP, dY_sb);
		double dbeta_A_sb = CalculateBetaA_SB(dT_sb, dP, dh_fb, dmflow_g);
		double dalpha_sb = CalculateAlphaA_SB(dT_sb, dP, dbeta_A_sb);

		double dh_H2O_sb = GetWaterVaporEnthalpy(dT_sb, dP);

		double dM_sb_z = drho_air * dbeta_A_sb * (_y_s_z[z] - _y_b_z[z]);
		dM_sb_z = bSwitchMT ? dM_sb_z : 0;
		double dH_sb_z = dM_sb_z * dh_H2O_sb;
		double dQ_sb_z = dalpha_sb * (_t_s_z[z] - _t_b_z[z]);
		dQ_sb_z = bSwitchHT ? dQ_sb_z : 0;

		vtempM_sb_z[z] = dM_sb_z;
		vtempH_sb_z[z] = dH_sb_z;
		vtempQ_sb_z[z] = dQ_sb_z;
	});
	// Trapezoidal rule
	for (size_t z = 1; z < nHeightGridPoints; z++)
	{
		_m_sb_z[z] = (vtempM_sb_z[z] + vtempM_sb_z[z - 1]) * 0.5;
		_h_sb_z[z] = (vtempH_sb_z[z] + vtempH_sb_z[z - 1]) * 0.5;
		_q_sb_z[z] = (vtempQ_sb_z[z] + vtempQ_sb_z[z - 1]) * 0.5;
	}
}

void CNLModelTD::ExchangeParticleParticle(const std::vector<double>& _tau_i, const std::vector<double>& _n_i, const std::vector<std::vector<std::vector<double>>>& _x_p_dfi, const std::vector<std::vector<std::vector<double>>>& _t_p_dfi, std::vector<std::vector<std::vector<double>>>& _q_pp_dfi)
{
	/// Discretization ///
	size_t nRTGridPoints = _n_i.size();
	size_t nSizeClasses = GetSizeClasses();
	size_t nMoistureClasses = GetMoistureClasses();

	/// Resize exchange flow vectors ///
	Resize3D(_q_pp_dfi, nSizeClasses, nMoistureClasses, nRTGridPoints);	// 3D-vector for storage of heat flow between particles of size class d and residence class i

	/// Model attributes ///
	double dP = GetTDInGasMYTP().P;			// Gas pressure at inlet [Pa]
	std::vector<std::vector<double>> vvAps_df = GetParticleSurfaceVector();

	double dC_Qpp = GetSettingsQpp();

	if (dC_Qpp == 0)
		return;

	std::vector<std::vector<std::vector<double>>> vvvtempQ_pp_dfi(nSizeClasses, std::vector<std::vector<double>>(nMoistureClasses, std::vector<double>(nRTGridPoints, 0)));	// 3D-vector for temporary storage of heat flow between particles of size class d and residence class i

	bool bSwitchHT = GetSettingsSwitchHT();				// Variable to switch of heat transfer

	double dT_p_mean = CalculateTemperatureMeanValue(_t_p_dfi, m_RTSubGridPoints);			// Mean particle temperature [K]
	double dX_p_mean = CalculateMoistureMeanValue(_x_p_dfi, m_RTSubGridPoints);			// Mean particle moisture content [kg/kg]
	double dalpha_pp = CalculateAlpha_PP(dT_p_mean, dX_p_mean, dT_p_mean, dP);	// Surface for interparticle heat transfer [m^2]

	for (size_t d = 0; d < nSizeClasses; d++)
	{
		for (size_t f = 0; f < nMoistureClasses; f++)
		{
			ParallelFor(nRTGridPoints, [&](size_t i)
			{
				// Get state of particle phase
				double dT_p_di = _t_p_dfi[d][f][i];
				double dX_p_di = _x_p_dfi[d][f][i];

				double dA_pp = vvAps_df[d][f];

				double dn_i = _n_i[i];
				double dQ_pp_dfi = dalpha_pp * dn_i * dC_Qpp * dA_pp * (dT_p_mean - dT_p_di);
				dQ_pp_dfi = bSwitchHT ? dQ_pp_dfi : 0;

				vvvtempQ_pp_dfi[d][f][i] = dQ_pp_dfi;
			});
		}
	}

	// Trapezoidal rule
	for (size_t d = 0; d < nSizeClasses; d++)
	{
		for (size_t f = 0; f < nMoistureClasses; f++)
		{
			ParallelFor(nRTGridPoints, [&](size_t i)
			{
				if (i == 0)
					return;
				double dn_i = _n_i[i];
				double dn_i_1 = _n_i[i - 1];
				_q_pp_dfi[d][f][i] = 0.5 * (vvvtempQ_pp_dfi[d][f][i - 1] / dn_i_1 + vvvtempQ_pp_dfi[d][f][i] / dn_i);
			});
		}
	}

	return;
}

void CNLModelTD::ExchangeParticleWall(const std::vector<double>& _tau_i, const std::vector<double>& _n_i, const std::vector<std::vector<std::vector<double>>>& _x_p_dfi, const std::vector<std::vector<std::vector<double>>>& _t_p_dfi, std::vector<std::vector<std::vector<double>>>& _q_pp_dfi)
{
	/// Discretization ///
	size_t nRTGridPoints = _n_i.size();
	size_t nSizeClasses = GetSizeClasses();
	size_t nMoistureClasses = GetMoistureClasses();

	/// Resize exchange flow vectors ///
	Resize3D(_q_pp_dfi, nSizeClasses, nMoistureClasses, nRTGridPoints);	// 3D-vector for storage of heat flow between particles of size class d and residence class i

	/// Model attributes ///
	double dP = GetTDInGasMYTP().P;			// Gas pressure at inlet [Pa]
	std::vector<std::vector<double>> vvAps_df = GetParticleSurfaceVector();

	double dC_Qw = GetSettingsQw();

	if (dC_Qw == 0)
		return;

	std::vector<std::vector<std::vector<double>>> vvvtempQ_pw_dfi(nSizeClasses, std::vector<std::vector<double>>(nMoistureClasses, std::vector<double>(nRTGridPoints, 0)));	// 3D-vector for temporary storage of heat flow between particles of size class d and residence class i

	bool bSwitchHT = GetSettingsSwitchHT();				// Variable to switch of heat transfer

	double dT_p_mean = CalculateTemperatureMeanValue(_t_p_dfi, m_RTSubGridPoints);			// Mean particle temperature [K]
	double dX_p_mean = CalculateMoistureMeanValue(_x_p_dfi, m_RTSubGridPoints);			// Mean particle moisture content [kg/kg]
	double dalpha_pw = CalculateAlpha_PW(dT_p_mean, dX_p_mean, dT_p_mean, dP);	// Surface for interparticle heat transfer [m^2]
	double dT_w = GetWallTemperature();
	double dh_fb = GetBedHeight();
	double dWall_circum = GetDryerCircumference();

	for (size_t d = 0; d < nSizeClasses; d++)
	{
		for (size_t f = 0; f < nMoistureClasses; f++)
		{
			ParallelFor(nRTGridPoints, [&](size_t i)
			{
				// Get state of particle phase
				double dT_p_di = _t_p_dfi[d][f][i];
				double dX_p_di = _x_p_dfi[d][f][i];

				// double dA_pw_pp = vvAps_df[d][f];	// Using surface of all paticles in class

				double dA_pw = dh_fb * dWall_circum;		// using surface area of dryer wall, in contact with particles

				double dn_i = _n_i[i];
				double dQ_pw_dfi = dalpha_pw * dn_i * dC_Qw * dA_pw * (dT_p_di - dT_w);
				dQ_pw_dfi = bSwitchHT ? dQ_pw_dfi : 0;

				vvvtempQ_pw_dfi[d][f][i] = dQ_pw_dfi;
			});
		}
	}

	// Trapezoidal rule
	for (size_t d = 0; d < nSizeClasses; d++)
	{
		for (size_t f = 0; f < nMoistureClasses; f++)
		{
			ParallelFor(nRTGridPoints, [&](size_t i)
			{
				if (i == 0)
					return;
				double dn_i = _n_i[i];
				double dn_i_1 = _n_i[i - 1];
				_q_pp_dfi[d][f][i] = 0.5 * (vvvtempQ_pw_dfi[d][f][i - 1] / dn_i_1 + vvvtempQ_pw_dfi[d][f][i] / dn_i);
			});
		}
	}
}

/// Balance functions ///
void CNLModelTD::BalanceParticlePhase(const std::vector<double>& _vtau_i, const std::vector<std::vector<std::vector<double>>>& _x_p_dfi, const std::vector<std::vector<std::vector<double>>>& _h_p_dfi, const std::vector<std::vector<std::vector<double>>>& _t_p_dfi, std::vector<std::vector<std::vector<double>>>& _x_p_dfi_update, std::vector<std::vector<std::vector<double>>>& _h_p_dfi_update, const std::vector<std::vector<std::vector<double>>>& _m_ps_dfi, const std::vector<std::vector<std::vector<double>>>& _q_ps_dfi, const std::vector<std::vector<std::vector<double>>>& _h_ps_dfi, const std::vector<std::vector<std::vector<double>>>& _q_pp_dfi, const std::vector<std::vector<std::vector<double>>>& _q_pw_dfi)
{
	/// Discretization ///
	size_t nRTGridPoints = _vtau_i.size();		// Number of residence time grid points [-]
	size_t nSizeClasses = GetSizeClasses();			// Number of particle size classes [-]
	size_t nMoistureClasses = GetMoistureClasses();	// Number of moisture classes [-]

	/// Resize solution vectors ///
	Resize3D(_x_p_dfi_update, nSizeClasses, nMoistureClasses, nRTGridPoints);	// 3D-vector for storage of updated moisture of particles of size class d, moisture class f and residence class i
	Resize3D(_h_p_dfi_update, nSizeClasses, nMoistureClasses, nRTGridPoints);	// 3D-vector for storage of updated enthalpy of particles of size class d, moisture class f and residence class i

	/// Model attributes ///
	double dP = GetTDInGasMYTP().P;			// Gas pressure at inlet [Pa]
	double dM_p = GetInitialBedMass();				// Initial dry bed mass [kg]
	std::vector<std::vector<double>> vvmfrac_df = GetMassFractionVector();	// Mass fraction of particles of size class d and moisture class f [-]

	/// Boundary conditions ///
	SBoundaryCondition sBC = GetBoundaryConditions();

	/// Balance equations ///
	// Iteration over size classes d and residence time classes i
	for (size_t d = 0; d < nSizeClasses; d++)
	{
		for (size_t f = 0; f < nMoistureClasses; f++)
		{
			// Set values for i = 0 (residence times boundary conditions)
			m_X_p_dfi_update[d][f][0] = sBC.X_p[f];
			m_h_p_dfi_update[d][f][0] = sBC.h_p[f];
			// Iterate over remaining residence time grid points
			for (size_t i = 1; i < nRTGridPoints; i++)
			{
				double dDelta_tau = _vtau_i[i] - _vtau_i[i - 1];	// Residence time class size [s]

				double dm_p_d = dM_p * vvmfrac_df[d][f];	// Mass of particles of size class d and moisture class f [kg]


				double dX_p_dfi = _x_p_dfi[d][f][i];
				double dh_p_dfi = _h_p_dfi[d][f][i];
				double dT_p_dfi = _t_p_dfi[d][f][i];

				// Calculate Xp(d,i)
				double dX_p_dfi_update = _x_p_dfi[d][f][i - 1] + dDelta_tau / dm_p_d * (-_m_ps_dfi[d][f][i]);	// Partice moisture [kg/kg]

				// Calculate hp(d,i)
				double dh_p_dfi_update = _h_p_dfi[d][f][i - 1] + dDelta_tau / dm_p_d * (-_q_ps_dfi[d][f][i] - _h_ps_dfi[d][f][i] + _q_pp_dfi[d][f][i] - _q_pw_dfi[d][f][i]);	// Particle enthalpy [J/kg]

				// Save values to storage vectors
				_x_p_dfi_update[d][f][i] = dX_p_dfi_update;
				_h_p_dfi_update[d][f][i] = dh_p_dfi_update;
			}
		}
	}
}

void CNLModelTD::BalanceSuspensionPhase(const std::vector<double>& _zeta_z, const std::vector<double>& _y_s_z, const std::vector<double>& _h_s_z, const std::vector<double>& _t_s_z, std::vector<double>& _y_s_z_update, std::vector<double>& _h_s_z_update, const std::vector<double>& _m_ps_z, const std::vector<double>& _h_ps_z, const std::vector<double>& _q_ps_z, const std::vector<double>& _m_sb_z, const std::vector<double>& _h_sb_z, const std::vector<double>& _q_sb_z) const
{
	/// Discretization ///
	size_t nHeightGridPoints = _zeta_z.size();	// Number of height grid points [-]

	/// Resize solution vectors ///
	_y_s_z_update.resize(nHeightGridPoints, 0.0);	// 1D-vector for storage of updated water loading of suspension phase at height class z
	_h_s_z_update.resize(nHeightGridPoints, 0.0);	// 1D-vector for storage of updated enthalpy of suspension phase at height class z

	/// Model attributes ///
	double dmflow_g = GetTDInGasMYTP().m;	// Gas flow at inlet [kg/s]
	double dP = GetTDInGasMYTP().P;			// Gas pressure at inlet [Pa]
	double dnu_0 = GetTDGasSplitFactor();	// Split factor vector of suspension and bubble phase
	double dmflow_s = (1 - dnu_0) * dmflow_g;	// Suspension gas flow [kg/s]

	/// Boundary conditions ///
	SBoundaryCondition sBC = GetBoundaryConditions();

	/// Balance equations ///
	// Boundary condition at z = 0
	_y_s_z_update[0] = sBC.Y_g;
	_h_s_z_update[0] = sBC.h_g;
	// Iteration over height elements z
	for (size_t z = 1; z < nHeightGridPoints; z++)
	{
		/* REWORK NEEDED */
		// Step-width of zeta
		double dDelta_zeta = _zeta_z[z] - _zeta_z[z - 1];

		// Mass flow of suspension gas
		// Calculate Ys(z)
		double dY_s_z_update = _y_s_z[z - 1] + dDelta_zeta / dmflow_s * (_m_ps_z[z] - _m_sb_z[z]);	// Moisture content at height class z [kg/kg]

		// Calculate hs(z)
		double dh_s_z_update = _h_s_z[z - 1] + dDelta_zeta / dmflow_s * (_h_ps_z[z] + _q_ps_z[z] - _h_sb_z[z] - _q_sb_z[z]); // Enthalpy at height class z [J/kg]

		// Save values to storage vectors
		_y_s_z_update[z] = dY_s_z_update;
		_h_s_z_update[z] = dh_s_z_update;
	}
}

void CNLModelTD::BalanceBubblePhase(const std::vector<double>& _zeta_z, const std::vector<double>& _y_b_z, const std::vector<double>& _h_b_z, const std::vector<double>& _t_b_z, std::vector<double>& _y_b_z_update, std::vector<double>& _h_b_z_update, const std::vector<double>& _m_sb_z, const std::vector<double>& _h_sb_z, const std::vector<double>& _q_sb_z) const
{
	/// Discretization ///
	size_t nHeightGridPoints = _zeta_z.size();	// Number of height grid points [-]

	/// Resize solution vectors ///
	_y_b_z_update.resize(nHeightGridPoints, 0.0);	// 1D-vector for storage of updated water loading of bubble phase at height class z
	_h_b_z_update.resize(nHeightGridPoints, 0.0);	// 1D-vector for storage of updated enthalpy of bubble phase at height class z

	/// Model attributes ///
	double dmflow_g = GetTDInGasMYTP().m;	// Gas flow at inlet [kg/s]
	double dP = GetTDInGasMYTP().P;			// Gas pressure at inlet [Pa]
	double dnu_0 = GetTDGasSplitFactor();	// Split factor vector of suspension and bubble phase
	double dmflow_b = (dnu_0)*dmflow_g;	// Bubble gas flow [kg/s]

	/// Boundary conditions ///
	SBoundaryCondition sBC = GetBoundaryConditions();

	/// Balance equations ///
	// Boundary condition at z = 0
	_y_b_z_update[0] = sBC.Y_g;
	_h_b_z_update[0] = sBC.h_g;
	// Iteration over height elements z
	for (size_t z = 1; z < nHeightGridPoints; z++)
	{
		// Step-width of zeta
		double dDelta_zeta = _zeta_z.at(z) - _zeta_z.at(z - 1);

		// Mass flow of suspension gas
		// Calculate Ys(z)
		double dY_b_z_update = 0;
		double dh_b_z_update = 0;
		// check bubble volume fraction
		if (dnu_0 >= 0)
		{
			// Calculate Ys(z)
			dY_b_z_update = _y_b_z[z - 1] + dDelta_zeta / dmflow_b * (_m_sb_z[z]);	// Moisture content at height class z [kg/kg]

			// Calculate hs(z)
			dh_b_z_update = _h_b_z[z - 1] + dDelta_zeta / dmflow_b * (_h_sb_z[z] + _q_sb_z[z]);	// Enthalpy at height class z [J/kg]
		}
		else    // if Geldart D, no bubble phase assumed - late fix SL
		{
			// do nothing and leave bubble proberites 0, as pre allocated -- does not work, becuase of NaNs in TD model and other solver issues for very small epsilo_B...
		}

		// Save values to storage vectors
		_y_b_z_update[z] = dY_b_z_update;
		_h_b_z_update[z] = dh_b_z_update;

	}
}

/// Additional functions ///
double CNLModelTD::IntegrateTrapez(double _val1, double _val2, double _step)
{
	return 0.5 * (_val1 + _val2) * _step;
}

void CNLModelTD::CalculateParticleTemperatures(const std::vector<std::vector<std::vector<double>>>& _h_p_dfi, const std::vector<std::vector<std::vector<double>>>& _x_p_dfi, std::vector<std::vector<std::vector<double>>>& _t_p_dfi) const
{
	double dT_Ref = GetReferenceTemperature();

	if (_h_p_dfi.empty())
		return;
	size_t nSizeClasses = _h_p_dfi.size();
	if (_h_p_dfi[0].empty())
		return;
	size_t nMoistureClasses = _h_p_dfi[0].size();
	if (_h_p_dfi[0][0].empty())
		return;
	size_t nRTGridPoints = _h_p_dfi[0][0].size();

	_t_p_dfi.resize(nSizeClasses, std::vector<std::vector<double>>(nMoistureClasses, std::vector<double>(nRTGridPoints, 0)));
	for (size_t d = 0; d < nSizeClasses; ++d)
		for (size_t f = 0; f < nMoistureClasses; ++f)
			for (size_t i = 0; i < nRTGridPoints; ++i)
				_t_p_dfi[d][f][i] = GetMoistParticleTemperature(_h_p_dfi[d][f][i], _x_p_dfi[d][f][i], dT_Ref, 101325);
}

void CNLModelTD::CalculateGasTemperatures(const std::vector<double>& _h_g_z, const std::vector<double>& _y_g_z, std::vector<double>& _t_g_z) const
{
	double dT_Ref = GetReferenceTemperature();
	if (_h_g_z.empty())
		return;
	size_t nHeightGridPoints = _h_g_z.size();

	_t_g_z.resize(nHeightGridPoints, 0);
	for (size_t z = 0; z < nHeightGridPoints; ++z)
	{
		_t_g_z[z] = GetMoistGasTemperature(_h_g_z[z], _y_g_z[z], dT_Ref, 101325);
	}
}

std::vector<double> CNLModelTD::CalculateRatioApAbed(double _time, double _h_fb, const std::vector<std::vector<double>>& _mfrac_df)
{
	double deps_fb = CaclulateBedPorosity(_time, _h_fb);		// Porosity of the fluidized bed [-]
	std::vector<double> vd_p_d = GetMeanDiameterVector();		// Vector with mean particle diameters of size classes d [m]
	std::vector<double> vApAbed_d(vd_p_d.size());				// Ap/Abed, ratio of interface area of particles to cross sectional area of the bed [-]

	// Get setting parameter for the calculation of ApAbed
	const auto C_ApAbed = GetSettingsApAbed();
	// If 0: Calculate ApAbed using the Sauter diameter
	// If 1: Calculate ApAbed for each diameter class using the corresponding mean diameter and the respective mass fraction
	switch (C_ApAbed)
	{
	case EArea::MEAN_AND_FRACTIONS:
	{
		for (size_t d = 0; d < vd_p_d.size(); ++d)
		{
			double dmfrac_d = 0;
			for (size_t f = 0; f < _mfrac_df[d].size(); ++f)
			{
				dmfrac_d += _mfrac_df[d][f];
			}
			vApAbed_d[d] = 6. / vd_p_d[d] * (1. - deps_fb) * dmfrac_d * _h_fb;		// TO DO: Check equation: does it make sense to include dmfrac??
		}
		break;
	}
	case EArea::MEAN:
	{
		for (size_t d = 0; d < vd_p_d.size(); ++d)
		{
			vApAbed_d[d] = 6. / vd_p_d[d] * (1. - deps_fb) * _h_fb;
		}
		break;
	}
	case EArea::SAUTER_DIAMETER:
	{
		double ddp_sauter = GetInParticleSauter(_time);
		for (size_t d = 0; d < vd_p_d.size(); ++d)
		{
			vApAbed_d[d] = 6. / ddp_sauter * (1. - deps_fb) * _h_fb;
		}
		break;
	}
	}

	return vApAbed_d;
}

double CNLModelTD::CaclulateBedPorosity(double _time, double _h_fb) const
{
	double dh_0 = GetInitialBedHeight();
	double deps_0 = GetInitialBedPorosity();
	double deps_fb = 1 - (1 - deps_0) * dh_0 / _h_fb;
	return deps_fb;
}

double CNLModelTD::CalculateMoistureMeanValue(const std::vector<std::vector<std::vector<double>>>& _vals, size_t _rtSubGridPoints) const
{
	std::vector<double> vtau_i = GetRTTauVector();
	std::vector<double> vn_i = GetRTNumberDensityVector();

	size_t nSizeClasses = GetSizeClasses();
	size_t nMoistureClasses = GetMoistureClasses();
	size_t nRTGridPoints = vtau_i.size();

	std::vector<std::vector<double>> vvmfrac_df = GetMassFractionVector();

	double dVal_p_mean = 0;
	for (size_t d = 0; d < nSizeClasses; d++)
	{
		for (size_t f = 0; f < nMoistureClasses; f++)
		{
			double dVal_p_mean_df = 0;
			double dN_df = 0;
			for (size_t i = 1; i < nRTGridPoints; i++)
			{
				dN_df += 0.5 * (vn_i[i - 1] + vn_i[i]) * (vtau_i[i] - vtau_i[i - 1]);
				dVal_p_mean_df += 0.5 * (vn_i[i - 1] * _vals[d][f][i - 1] + vn_i[i] * _vals[d][f][i]) * (vtau_i[i] - vtau_i[i - 1]);
			}
			dVal_p_mean_df /= dN_df;
			dVal_p_mean += vvmfrac_df[d][f] * dVal_p_mean_df;
		}
	}
	return dVal_p_mean;
}

double CNLModelTD::CalculateTemperatureMeanValue(const std::vector<std::vector<std::vector<double>>>& _vals, size_t _rtSubGridPoints) const
{
	std::vector<double> vtau_i = GetRTTauVector();
	std::vector<double> vn_i = GetRTNumberDensityVector();

	size_t nSizeClasses = GetSizeClasses();
	size_t nMoistureClasses = GetMoistureClasses();
	size_t nRTGridPoints = vtau_i.size();

	std::vector<std::vector<double>> vvmfrac_df = GetMassFractionVector();
	std::vector<std::vector<double>> vvAps_df = GetParticleSurfaceVector();

	// Total surface
	double dAps = 0;
	for (size_t d = 0; d < nSizeClasses; d++)
		for (size_t f = 0; f < nMoistureClasses; f++)
			dAps += vvAps_df[d][f];

	double dVal_p_mean = 0;
	for (size_t d = 0; d < nSizeClasses; d++)
	{
		for (size_t f = 0; f < nMoistureClasses; f++)
		{
			double dVal_p_mean_df = 0;
			double dN_df = 0;
			for (size_t i = 1; i < nRTGridPoints; i++)
			{
				dN_df += 0.5 * (vn_i[i - 1] + vn_i[i]) * (vtau_i[i] - vtau_i[i - 1]);
				dVal_p_mean_df += 0.5 * (vn_i[i - 1] * _vals[d][f][i - 1] + vn_i[i] * _vals[d][f][i]) * (vtau_i[i] - vtau_i[i - 1]);
			}
			dVal_p_mean_df /= dN_df;
			dVal_p_mean += vvAps_df[d][f] / dAps * dVal_p_mean_df;
		}
	}
	return dVal_p_mean;
}

double CNLModelTD::CalculateTauMax(double _n_tau_max) const
{
	double dmlfow_p_in = GetTDInParticleMYTP().m;
	double dM_p = GetInitialBedMass();
	double dtau_m = dM_p / dmlfow_p_in;
	// Determine tau_max such that a determined fraction of particles are considered (value: "dN_tau_max"), through integration of the number density function
	double dK = GetSettingsKtis();
	double dtau_max = 0;
	if (dK == 1)
	{
		dtau_max = -log(1 - _n_tau_max) / (dmlfow_p_in / dM_p);
	}
	else
	{
		double dtau = 100;
		while (IntegrateTIS(dtau, dtau_m, dK) < _n_tau_max)
		{
			dtau += 100;
		}
		dtau_max = dtau;
	}

	return dtau_max;
}

std::vector<double> CNLModelTD::CalculateRTTauVector(double _tau_max) const
{
	size_t nRTGridPoints = GetRTGridPoints();
	std::vector<double> vtau_i(nRTGridPoints);

	double dDelta_tau = _tau_max / (nRTGridPoints - 1);
	for (size_t i = 0; i < nRTGridPoints; ++i)
	{
		vtau_i[i] = dDelta_tau * i;
	}
	return vtau_i;
}

std::vector<double> CNLModelTD::CalculateRTDensityVector(double _tau_max) const
{
	double dmlfow_p_in = GetTDInParticleMYTP().m;
	double dM_p = GetInitialBedMass();
	size_t nRTGridPoints = GetRTGridPoints();
	std::vector<double> vn_i(nRTGridPoints, 0);
	std::vector<double> vtau_i = CalculateRTTauVector(_tau_max);
	double dtau_m = dM_p / dmlfow_p_in;

	double dK = GetSettingsKtis();
	// Ideally mix reactor
	if (dK == 1)
	{

		for (size_t i = 0; i < nRTGridPoints; ++i)
		{
			vn_i[i] = 1.0 / dtau_m * exp(-vtau_i[i] / dtau_m);
		}
	}
	// TIS model
	else
	{
		double dGamma = Gamma(dK);
		vn_i = TISModel(vtau_i, dtau_m, dK, dGamma);
		/*for (size_t i = 0; i < nRTGridPoints; ++i)
		{
			vn_i[i] = 1.0 / dtau_m * pow((vtau_i[i] / dtau_m), (dK - 1)) * pow(dK, dK) / Gamma(dK) * exp(-vtau_i[i] * dK / dtau_m);
		}*/
	}

	return vn_i;
}

std::vector<std::vector<double>> CNLModelTD::CalculateParticleSurfaceVector(double _time, const std::vector<double>& _d_d)
{
	std::vector<std::vector<double>> vvmfrac_df = CalculateMassFractionVector(_time);
	size_t nSizeClasses = GetSizeClasses();
	size_t nMoistureClasses = GetMoistureClasses();

	double drho_s = GetDryParticleDensity(273.15, 101325);
	double dM_p = GetInitialBedMass();

	std::vector<std::vector<double>> vvA_ps_df(nSizeClasses, std::vector<double>(nMoistureClasses)); //Total particle surface in each particle size classes [m2]
	for (size_t d = 0; d < nSizeClasses; d++)
		for (size_t f = 0; f < nMoistureClasses; f++)
			vvA_ps_df[d][f] = 6. * dM_p * vvmfrac_df[d][f] / (drho_s * _d_d[d]);

	return vvA_ps_df;
}

std::vector<std::vector<double>> CNLModelTD::CalculateMassFractionVector(double _time) const
{
	return GetInParticleMassfractions(_time);		// Proportion of each size classes, based on mass fraction. Copied from CFluidizedBedDryer
}

double CNLModelTD::CalculateAlpha_PS(double _t_ps, double _p_g, double _d_p, double _ratioApAbed_d) const
{
	// Get unit parameters
	double deps_mf = GetInitialBedPorosity();				// Bed porosity at minimum fluidization, equal to INITIAL bed porosity [-]
	double dRe_mf = GetRemf();								// Reynolds number at minimum fluidization [-]
	double du_0 = GetU0();									// Superficial velocity [m/s]

	// Get material property
	// TO DO: Moist Gas
	double drho_g = GetDryGasDensity(_t_ps, _p_g);				// Density [kg/m^3]
	double deta_g = GetDryGasViscosity(_t_ps, _p_g);			// Dynamic viscosity [Pa s]
	double dnu_g = deta_g / drho_g;								// Kinematic viscosity [m^2 / s]
	double dcp_g = GetDryGasHeatCapacity(_t_ps, _p_g);			// Heat capacity [kg /(kg K)]
	double dlambda_g = GetDryGasThermalConductivity(_t_ps, _p_g);		// Thermal conductivity [W/(m K)]

	// Nusselt number calculation
	double dRe = dRe_mf / deps_mf;								// Renolds number in the fluidized bed [-]
	double dPr = dnu_g * dcp_g * drho_g / dlambda_g;				// Prandtl number [-]
	double dNu_lam = 0.644 * pow(dRe, 0.5) * pow(dPr, 1. / 3.);		// Laminar Nusselt number [-]
	double dNu_tur = 0.037 * pow(dRe, 0.8) * dPr / (1. + 2.443 * pow(dRe, -0.1) * (pow(dPr, 2. / 3.) - 1.)); // Turbulent Nusselt number [-]
	double dNu_p = 2. + pow(pow(dNu_lam, 2.) + pow(dNu_tur, 2.), 0.5);	// Nusselt number for single particle [-]
	double dNu_ps = (1. + 1.5 * (1. - deps_mf)) * dNu_p;			// Nusselt number between particles and suspension gas [-]

	double dd_p = GetInSauter();
	double dRe_0 = du_0 * dd_p / dnu_g;	//	Renolds number at superficial velocity [-]

	//double dRe_0 = du_0 * _dd_p / dnu_g;	//	Renolds number at superficial velocity [-]
	double dNu_ps_apparent = dRe_0 * dPr / _ratioApAbed_d * log(1. + dNu_ps * _ratioApAbed_d / dRe_0 / dPr); // Nu_ps', apparent Nusselt number [-]

	// Heat transfer coefficient
	double dalpha_ps = dNu_ps_apparent * dlambda_g / dd_p;				//alpha_ps, heat transfer coefficient between particles and suspension gas [W / (m^2 K)]
	//double dalpha_ps = dNu_ps_apparent * dlambda_g / _dd_p;				//alpha_ps, heat transfer coefficient between particles and suspension gas [W / (m^2 K)]

	// TO DO error
	if (dalpha_ps != dalpha_ps)
		dalpha_ps = 0.;

	return dalpha_ps;
}

double CNLModelTD::CalculateBeta_PS(double _t_ps, double _y_s, double _p, double _d_p, double _ratioApAbed_d) const
{
	// Get unit parameters
	double deps_mf = GetInitialBedPorosity();				// Bed porosity at minimum fluidization, equal to INITIAL bed porosity [-]
	double dRe_mf = GetRemf();								// Reynolds number at minimum fluidization [-]
	double du_0 = GetU0();									// Superficial velocity [m/s]

	// Get material property
	double drho_g = GetMoistGasDensity(_t_ps, _p, _y_s);				// Density [kg/m^3]
	// TO DO: Moist Gas
	double deta_g = GetDryGasViscosity(_t_ps, _p);			// Dynamic viscosity [Pa s]
	double dnu_g = deta_g / drho_g;								// Kinematic viscosity [m^2 / s]

	double ddelta_g = GetDiffusionCoefficient(_t_ps, _p);		// Diffusion coefficient of water in gas [m2/s]

	// Model of Gnielinski (1980)
	double dRe = dRe_mf / deps_mf;						// Renolds number in the fluidized bed [-]
	double dSc = dnu_g / ddelta_g;						// Schmidt number [-]
	double dSh_lam = 0.644 * pow(dRe, 0.5) * pow(dSc, 1. / 3.); // Laminar Sherwood number [-]
	double dSh_tur = 0.037 * pow(dRe, 0.8) * dSc / (1. + 2.443 * pow(dRe, -0.1) * (pow(dSc, 2. / 3.) - 1.)); // Turbulent Sherwood number [-]
	double dSh_p = 2. + pow(pow(dSh_lam, 2.) + pow(dSh_tur, 2.), 0.5); // Sherwood number for single particle [-]
	double dSh_bed = (1. + 1.5 * (1. - deps_mf)) * dSh_p; //Sherwood number [-]

	// Model of Groenewold and Tsotsas (1997)
	double dd_p = GetInSauter();
	double dRe_0 = du_0 * dd_p / dnu_g; // Renolds number at superficial velocity [-]

	//double dRe_0 = du_0*_dd_p / dnu_g; // Renolds number at superficial velocity [-]
	double dSh_ps_apparent = dRe_0 * dSc / _ratioApAbed_d * log(1. + dSh_bed * _ratioApAbed_d / dRe_0 / dSc); //Sh_ps', apparent Sherwood-Zahl

	// Mass transfer coefficient
	double dbeta_ps = dSh_ps_apparent * ddelta_g / dd_p; // Mass transfer coefficient between particles and suspension gas [m/s]
	//double dbeta_ps = dSh_ps_apparent*ddelta_g / _dd_p; // Mass transfer coefficient between particles and suspension gas [m/s]

	// TO DO error
	if (dbeta_ps != dbeta_ps)
		dbeta_ps = 0.;

	return dbeta_ps;
}

double CNLModelTD::CalculateNormalizedDryingCurve(double _x, double _x_eq) const
{
	// Get Parameter
	double dk = GetDryingCurveK();		// Parameter for normalized drying curve [-]
	double dX_cr = GetXcr();			// Critical moisture content [kg/kg]

	// Normalized drying curve
	double dnu = 1.;
	// If the moisture content X is larger than critical moisture content X_cr, nu is limited to 1 (first drying period)
	if (_x >= dX_cr)
		dnu = 1.;
	else
	{
		// TO DO check
		// Neglect case for which X is smaller than X_eq, particle would take moisture from suspension
		if (_x <= _x_eq)
			return 0;

		double deta = (_x - _x_eq) / (dX_cr - _x_eq); // Normalized moistrue content [-]

		 // Normalized Drying Curve based model of van Meel (1958)
		dnu = dk * deta / (1. + deta * (dk - 1.)); // Normalized Drying Curve
	}
	// TO DO check for nan
	return dnu;
}

double CNLModelTD::CalculateAlphaA_SB(double _t_sb, double _p, double _beta_sb_A_sb) const
{
	// Get material property
	// TO DO: Moist Gas
	double drho_g = GetDryGasDensity(_t_sb, _p);					// Density [kg/m^3]
	double dcp_g = GetDryGasHeatCapacity(_t_sb, _p);				// Heat capacity [kg /(kg K)]
	double dlambda_g = GetDryGasThermalConductivity(_t_sb, _p);		// Thermal conductivity [W/(m K)]
	double ddelta_wg = GetDiffusionCoefficient(_t_sb, _p);		// Diffusion coefficient of water in gas [m2/s]

	double dLe = dlambda_g / dcp_g / drho_g / ddelta_wg;	//	Lewis number [-]
	double dexponent_m = 1. / 3.;							// Exponent m from Martin (1980) [-]

	// TO DO error
	if (dLe < 0)
		return 0.;

	double dalpha_sb_A_sb = drho_g * dcp_g * _beta_sb_A_sb * pow(dLe, 1. - dexponent_m); //alpha_SB * A_SB, kinetic heat transfer coefficient, A.3.74

	return dalpha_sb_A_sb;
}

double CNLModelTD::CalculateBetaA_SB(double _t_sb, double _p, double _h_fb, double _m_g_dot) const
{
	// Get material property
	// TODO: Moist Gas
	double drho_g = GetDryGasDensity(_t_sb, _p);			// Density [kg/m^3]

	double d_H_NTU = GetSettingsH_NTU();

	// TODO: Check model
	double dNTU_sb = _h_fb / d_H_NTU; // default value = 0.05;							// Number of transfer units [-]
	double dbeta_sb_A_sb = dNTU_sb * _m_g_dot / drho_g;	// Kinetic mass transfer coefficient times Area

	return  dbeta_sb_A_sb;
}

// TODO comments
double CNLModelTD::CalculateAlpha_PW(double _t_p, double _x_p, double _t_g, double _p) const
{
	/// Get unit parameters
	double deps_mf = GetInitialBedPorosity();				// Bed porosity at minimum fluidization, equal to INITIAL bed porosity [-]
	double deps_fb = GetBedPorosity();						// Bed porosity at operating conditions [-]
	double ddp_p = GetInSauter();							// Sauter diameter of particles [m]

	// Material data of dry air
	double dMM_g = GetGasMolarMass();								// Molar mass [kg/mol]
	double dcp_g = GetDryGasHeatCapacity(_t_g, _p);				// Heat capacity [J/kgK]
	double dlambda_g = GetDryGasThermalConductivity(_t_g, _p);	// Thermal conductivity [W/(m*K)]

	// Material data of moist particles
	double dcp_p = GetMoistParticleHeatCapacity(_x_p, _t_p, _p);	// Heat capacitiy of moist particles
	double drho_p = GetMoistParticleDensity(_x_p, _t_p, _p);		// Density of moist particles [kg/m^3]

	// Calculate modified free path of the gas molecule
	double dgamma = 1. / (pow(10., 0.6 - (1000. / _t_g + 1.) / 2.8) + 1); // gamma, A.3.85

	double dl = 2. * (2. / dgamma - 1.) * sqrt(2. * MATH_PI * MOLAR_GAS_CONSTANT * _t_g / dMM_g) * dlambda_g / _p / (2. * dcp_g - MOLAR_GAS_CONSTANT / dMM_g); //modified free path of the gas molecule, A.3.84

	//Calculate parameter Z and N
	double dC_k = 2.6; //C_k, A.3.82
	double dNu_pw_max = 4. * ((1. + 2. * dl / ddp_p) * log(1. + ddp_p / 2. / dl) - 1.); //Nu_PW,max, A.3.83
	double dZ = 1. / 6. * drho_p * dcp_p / dlambda_g * sqrt(STANDARD_ACCELERATION_OF_GRAVITY * pow(ddp_p, 3.) * (deps_fb - deps_mf) / 5. / (1 - deps_mf) / (1 - deps_fb)); //Parameter Z, A.3.80
	double dN = dNu_pw_max / dC_k / dZ; //Parameter N, A.3.81

	//Calculate Nusselt Number and heat transfer coefficient
	double dNu_pw = (1. - deps_fb) * dZ * (1. - exp(-dN));
	double dalpha_pw = dNu_pw * dlambda_g / ddp_p;

	return dalpha_pw;
}

double CNLModelTD::CalculateAlpha_SW(double _t_s, double _y_s, double _p, double _h_fb, double _v_s) const
{
	// Material data of dry air
	double dMM_g = GetGasMolarMass();								// Molar mass [kg/mol]
	double dcp_s = GetDryGasHeatCapacity(_t_s, _p);				// Heat capacity [J/kgK]
	double dlambda_s = GetDryGasThermalConductivity(_t_s, _p);	// Thermal conductivity [W/(m*K)]
	double deta_s = GetDryGasViscosity(_t_s, _p);					// Dynamic viscosity [Pa s]
	double drho_s = GetDryGasDensity(_t_s, _p);					// Gas density [kg/m^3]

	double dPr_s = deta_s * dcp_s / dlambda_s;
	double dRe_s = _v_s * _h_fb * drho_s / deta_s;

	//Calculate Nusselt Number and heat transfer coefficient
	double dNu_sw = 0.037 * pow(dRe_s, 0.8) * dPr_s / (1 + 2.443 * pow(dRe_s, -0.1) * (pow(dPr_s, 2. / 3.) - 1));
	double dalpha_sw = dNu_sw * dlambda_s / _h_fb;

	return dalpha_sw;
}

double CNLModelTD::CalculateA_PW() const
{
	double dh_fb = GetBedHeight();
	double dd_fb = GetBedDiameter();
	double dA_pw = MATH_PI * dd_fb * dh_fb;
	return dA_pw;
}

double CNLModelTD::CalculateAlpha_PP(double _t_p, double _x_p, double _t_g, double _p) const
{
	double dalpha_pw = CalculateAlpha_PW(_t_p, _x_p, _t_g, _p);
	// Assumption: alpha_PP equals alpha_PW
	double dalpha_pp = dalpha_pw;

	if (dalpha_pp != dalpha_pp)
		dalpha_pp = 0.;

	return dalpha_pp;
}

double CNLModelTD::CalculateA_PP(double _t_p_mean, double _p) const
{
	double ddp_d = GetInSauter();							// Sauter diameter of particles [m]
	double dM_p = GetInitialBedMass();
	double drho_p = GetDryParticleDensity(_t_p_mean, _p);
	double dA_ps_d = 6 * dM_p / (drho_p * ddp_d);
	return dA_ps_d;
}

/// Clear functions ///
void CNLModelTD::ClearModel()
{
	ClearVariables();

	m_HeightSubGridPoints = 0;
	m_RTSubGridPoints = 0;

	Clear1D(m_iY_b_z);
	Clear1D(m_ih_b_z);
	Clear1D(m_Y_b_z);
	Clear1D(m_h_b_z);
	Clear1D(m_Y_b_z_update);
	Clear1D(m_h_b_z_update);

	Clear1D(m_iY_s_z);
	Clear1D(m_ih_s_z);
	Clear1D(m_Y_s_z);
	Clear1D(m_h_s_z);
	Clear1D(m_Y_s_z_update);
	Clear1D(m_h_s_z_update);

	Clear3D(m_iX_p_dfi);
	Clear3D(m_ih_p_dfi);
	Clear3D(m_X_p_dfi);
	Clear3D(m_h_p_dfi);
	Clear3D(m_X_p_dfi_update);
	Clear3D(m_h_p_dfi_update);
}