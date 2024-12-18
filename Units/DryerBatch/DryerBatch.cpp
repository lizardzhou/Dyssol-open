/* Copyright (c) 2020, Dyssol Development Team. All rights reserved. This file is part of Dyssol. See LICENSE file for license information. */

#define DLL_EXPORT
#include "DryerBatch.h"
#include <sstream>
#include <set>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
//#include <jsoncpp/json/json.h>

extern "C" DECLDIR CBaseUnit* DYSSOL_CREATE_MODEL_FUN()
{
	return new CDryerBatch();
}

//////////////////////////////////////////////////////////////////////////
/// Batch dryer with a water inlet and a fluidization gas inlet, and an outlet for exhaust gas
/// - Heat transfer based on model from Maksym
/// - Mass transfer of water from particle to fluidization gas
/// - Calculation of enthalpy in holdup: using reference temperature 273.15 K and set reference energy as zero. Reason: "absolute enthalpy" isn't meaningful by itself because enthalpy is a RELATIVE MEASUREMENT ON ENERGY, it's a state function depending on changes, not absolute values.
/// - Drying kinetics (1st and 2nd drying stages) - UNDER CONSTRUCTION
/// - Height discretization for gas phase moisture & temperature - UNDER CONSTRUCTION

void CDryerBatch::CreateBasicInfo()
{
	/// Set basic unit info ///
	SetUnitName("Dryer (batch with heat and mass transfer)");
	SetAuthorName("Xiye Zhou, Alexander Hanke");
	SetUniqueID("9771F953-C2F7-417C-95CC-83C503433FB8"); 
	//test for environment
}

void CDryerBatch::CreateStructure()
{
	/// Add ports ///
	AddPort("InletLiquid", EUnitPort::INPUT);
	AddPort("InletFluidizationGas", EUnitPort::INPUT);
	AddPort("InletNozzleAir", EUnitPort::INPUT);
	AddPort("OutletExhaustGas", EUnitPort::OUTPUT);

	/// Add holdups for each phase ///
	AddHoldup("HoldupSolid");
	AddHoldup("HoldupLiquid");
	AddHoldup("HoldupGas");

	/// Add unit parameters ///
	std::vector<size_t>	items;
	std::vector<std::string> itemNames;
	// particle properties
	AddStringParameter("Particle properties", "", "");
	AddConstRealParameter("A_P", 0, "m2", "Total surface of all particles in the granulator.\nIf = 0, PSD is used to calculate surface area", 0);
	AddConstRealParameter("Delta_f", 100, "mum", "Thickness of the water film on particles in micrometer", 1);
	// environment
	AddConstRealParameter("theta_env", 20, "deg C", "Temperature of the environment, assumed to be constant", 0, 100);
	// Process chamber information
	AddStringParameter("Chamber information", "", "");
	AddConstRealParameter("H_tempProbe", 0.07, "m", "Height of temperature probe for chamber temperature over distributor.", 0, 0.35);
	AddConstRealParameter("H_nozzle", 0.25, "m", "Height of two-fluid nozzle over gas distributor.", 0, 0.35);
	AddConstRealParameter("d_bed", 0.18, "m", "Bed diameter", 0.1, 0.2);
	AddConstRealParameter("H_bedFix", 0.08, "m", "Bed height without fluidization", 0.08, 0.2);
	AddConstRealParameter("H_plant", 0.35, "m", "Plant height", 0, 30);
	// Bed properties, use for development of further models
	AddStringParameter("Bed properties", "", "");
	AddConstRealParameter("eps_0", 0.35, "-", "Fixed bed porosity", 0, 1);
	AddConstRealParameter("u_mf", 0, "m/s", "Minimal fluidization velocity\nIf 0, calculate use Wen&Yu correlation", 0, 1);
	items = { 0, 1 };
	itemNames = { "Martin", "Lehmann" };
	AddComboParameter("Bed porosity calculation", 0, items, itemNames, "Methods for calculating the bed porosity during (homogeneous) fluidization. Choose between Martin(VDI-Waermeatlas, chapter M5) and Lehmann dissertation(2021). \nMartin: porosity is a function of Re_mf and Re_elu. \nLehmann: porosity is a function of suspension gas velocity.");
	//AddConstIntParameter ("N_el", 1, "", "Number of hight discretization layers", 1); // # of height discretization layers

	// Heat and mass transfer 
	AddStringParameter("Calculation of heat and mass transfer", "", "");
	AddConstRealParameter("f_a", 0, "-", "Quotient alpha_PF/alpha_GF. Used to calculate heat transfer coefficient particle to water film from heat transfer coeeficient gas to water film.", 0, 1);
	items = { 0, 1, 2 };
	itemNames = { "Martin", "Groenewold & Tsostas", "Self-defined" };
	AddComboParameter("Heat & mass transfer methods", 0, items, itemNames, "Methods for calculating heat and mass transfer coefficient, choose between Martin(VDI-Waermeatlas, chapter M5), Groenewolds & Tsostas (see Rieck dissertation (2020)) and self-defined values.");
	AddConstRealParameter("alpha_GP", 0, "W/(m2*K)", "Heat transfer coefficient between gas and particle.", 0);
	AddConstRealParameter("alpha_GF", 0, "W/(m2*K)", "Heat transfer coefficient between gas and water film.", 0);
	AddConstRealParameter("alpha_PF", 0, "W/(m2*K)", "Heat transfer coefficient between particle and water film.", 0);
	AddParametersToGroup("Heat & mass transfer methods", "Self-defined", { "alpha_GP", "alpha_GF", "alpha_PF" });
	AddConstRealParameter("beta_GP", 0, "m/s", "Mass transfer coefficient for liquid from gas to particle\nIf negative, calculating using methods reported in Martin(VDI-Waermeatlas, chapter M5).");
	AddConstRealParameter("Tolerance temperature", 0.1, "K", "Tolerance for temperature difference for calculating the heat transfer between phases. Temperature difference smaller than the tolerance will be set to zero.", 1e-4, 2);
	AddConstRealParameter("Shrink factor alpha GP", 0.1, "-", "Shrink factor for heat transfer coefficient GP at small temperature difference.");
	AddConstRealParameter("Shrink factor alpha GF", 0.1, "-", "Shrink factor for heat transfer coefficient GF at small temperature difference.");
	AddConstRealParameter("Shrink factor beta", 0.1, "-", "Shrink factor for mass transfer coefficient GF.");

	// Drying kinetics calculation
	AddStringParameter("Drying kinetics", "", "");
	items = {0, 1};
	itemNames = {"REA", "NCDC"};
	AddComboParameter("Methods", 0, items, itemNames, "Method for calculating relative drying rate, choose between Normalized Characteristic Drying Curve (NCDC) or Reaction Engineering Approach (REA, Lehmann dissertation, Fig.4.14)");
	AddConstRealParameter("A", 0.96, "-", "Coefficient for f=1-A*exp(B*(X-X_eq)^C).", 0, 1);
	AddConstRealParameter("B", -19.63, "-", "Coefficient for f=1-A*exp(B*(X-X_eq)^C).", -100, 100);
	AddConstRealParameter("C", 0.73, "-", "Coefficient for f=1-A*exp(B*(X-X_eq)^C).", 0, 1);
	AddConstRealParameter("X_eq", 0.01, "kg/kg", "Equilibrium water content of the particles (transition between 2nd and 3rd drying period).", 0, 0.1);
	AddConstRealParameter("X_cr", 0.06, "kg/kg", "Critical water content of the particles (transition between 1st and 2nd drying period).", 0, 1.0);
	AddConstRealParameter("k_dc", 2, "-", "k for normalized drying curve. \nthe normalized drying curve of the material describes the relative drying rate f: \n f = k*normX/(1+normX*(k-1)), \n normX=(X-X_eq)/(X_cr-X_eq).", 0.1, 10);
	AddParametersToGroup("Methods", "REA", { "X_eq", "A", "B", "C" });
	AddParametersToGroup("Methods", "NCDC", { "X_cr", "X_eq", "k_dc" });
	itemNames = {"Y_eq", "Y_sat"};
	AddComboParameter("Use Y_eq or Y_sat?", 0, items, itemNames, "Choose between Y_eq (equilibrium gas moisture calculated from desorption isotherm) and Y_sat (saturation gas moisture) to calculate the water evaporation rate.");
	AddCheckBoxParameter("Use relative drying rate?", false, "Tick this box to use relative drying rate to calculate the water evaporation rate.");

	/// Read desorption isothermal information
	//AddStringParameter("Path to X_eq data", "C:\\", "Location of equilibrium moisture content with temperature, must be a csv file.");
	//AddCompoundParameter("X_eq compound", "Compound for with the Xeq values are contained in Path Xeq.");
	//AddConstRealParameter("x_l,eq,min", 0, "mass %", "Minimum measured equilibirum moisture fraction of particles.\nIf 0, materials database will be used to calculate euqilibirum moisture content. (Further see Path X_eq)", 0, 100);
	//AddConstRealParameter("theta_eq,min", 0, "degree celsius", "Temperature correponding to w_l,eq,min", 0, 100);	

	//AddStringParameter("Tolerance for solver", "", "");
	//AddConstRealParameter("Relative tolerance", 0.0, "-", "Solver relative tolerance. Set to 0 to use flowsheet-set value", 0);
	//AddConstRealParameter("Absolute tolerance Y", 0.0, "-", "Solver absolute tolerance for gas moisture content Y.\nSet to 0 to use flowsheet-set value", 0);
	//AddConstRealParameter("Absolute tolerance T", 0.0, "-", "Solver absolute tolerance for eqData.temperatures T.\nSet to 0 to use flowsheet-set value", 0);

	/// Tolerance ///
	AddConstRealParameter("Relative tolerance", 1e-3, "-", "Relative tolerance for the solver");
	AddConstRealParameter("Absolute tolerance", 1e-5, "-", "Absolute tolerance for the solver");
	AddConstRealParameter("Maximal time step", 1e-3, "s", "Maximal time step for iteration");

	/// Debug information ///
	AddCheckBoxParameter("Print intermediate results", true, "Tick this box to print intermediate results on simulation window.");

	/// Set this unit as user data of model ///
	m_model.SetUserData(this);
	m_model.m_unit = this; // set this unit as unit for internal calculation
}

void CDryerBatch::Initialize(double _time)
{
	/// Check flowsheet parameters ///
	if (!IsPhaseDefined(EPhase::SOLID))		RaiseError("Solid phase has not been defined.");
	if (!IsPhaseDefined(EPhase::LIQUID))	RaiseError("Liquid phase has not been defined.");
	if (!IsPhaseDefined(EPhase::VAPOR))		RaiseError("Gas phase has not been defined.");
	if (!IsDistributionDefined(DISTR_SIZE))	RaiseError("Size distribution has not been defined.");

	const bool printResult = GetCheckboxParameterValue("Print intermediate results");

	/// Tolerance for solver
	const auto rtol = GetConstRealParameterValue("Relative tolerance");
	const auto atol = GetConstRealParameterValue("Absolute tolerance");
	m_model.SetTolerance(rtol != 0.0 ? rtol : GetRelTolerance(), atol != 0.0 ? atol : GetAbsTolerance());
	/// Set maximal time step for solver
	const double t = GetConstRealParameterValue("Maximal time step");
	//m_solver.SetMaxStep(t);

	/// Settings
	//bool calcNdc = GetCheckboxParameterValue("calcNdc");
	//ignoreVaporInput = GetCheckboxParameterValue("ignoreVaporInput");
	//useREA = GetCheckboxParameterValue("useREA");
	//dryingCurveSetting = GetComboParameterValue("DryingCurve");

	/// Get holdup ///
	m_holdupSolid = GetHoldup("HoldupSolid");
	m_holdupLiquid = GetHoldup("HoldupLiquid");
	m_holdupGas = GetHoldup("HoldupGas");
	const mass mSolidHoldup = m_holdupSolid->GetPhaseMass(_time, EPhase::SOLID); 
	if (mSolidHoldup == 0) // empty run
	{
		particlesGlobal = false;
		RaiseWarning("No particles in system.\nSolids and liquids will be ignored.");
	}
	const std::vector<double> Grid = GetNumericGrid(DISTR_SIZE);
	const std::vector<double> q_3 = m_holdupSolid->GetPSD(_time, PSD_q3, EPSDGridType::DIAMETER);
	const length d32 = CalculateHoldupSauter(_time);

	/// Get pointers to streams ///
	m_inLiquidStream = GetPortStream("InletLiquid");
	m_inGasStream = GetPortStream("InletFluidizationGas");
	m_inNozzleAirStream = GetPortStream("InletNozzleAir");
	m_outExhaustGasStream = GetPortStream("OutletExhaustGas");

	/// Pull compound data ///
	//this->PullCompoundDataFromDatabase(_time);

	/// Setup vapor (gas phase) stream ///
	//m_VaporStream = AddStream("Vapor");
	//std::vector<double> vaporStreamCompoundSetup(size(compoundKeys), 0);
	//m_VaporStream->SetCompoundsFractions(_time, vaporStreamCompoundSetup);
	//m_VaporStream->SetPhaseFraction(_time, EPhase::SOLID, 0);
	//m_VaporStream->SetPhaseFraction(_time, EPhase::LIQUID, 0);
	//m_VaporStream->SetPhaseFraction(_time, EPhase::VAPOR, 1);
	//m_VaporStream->SetCompoundFraction(_time, compoundKeys.at(indicesOfVaporOfPhaseChangingCompound.second), EPhase::VAPOR,1);
	//m_VaporStream->SetMassFlow(_time, 0);

	/// Expansion part and working holdup
	//m_Expander = AddHoldup("Expander");
	//m_Expander->SetPhaseFraction(_time, EPhase::SOLID, 0);
	//m_Expander->SetPhaseFraction(_time, EPhase::LIQUID, 0);
	//m_Expander->SetPhaseFraction(_time, EPhase::GAS, 1);
	//m_Expander->SetCompoundsFractions(_time, m_inGasStream->GetCompoundsFractions(_time));
	//m_Expander->SetMass(_time, 2);
	//workingHoldup = AddHoldup("workingHoldup");
	//workingHoldup->SetPhaseFraction(_time, EPhase::SOLID, 0);
	//workingHoldup->SetPhaseFraction(_time, EPhase::LIQUID, 1);
	//workingHoldup->SetPhaseFraction(_time, EPhase::GAS, 0);
	//workingHoldup->SetCompoundsFractions(_time, m_inSuspensionStream->GetCompoundsFractions(_time));
	//workingHoldup->SetMass(_time, 0);
	
	/////////////////////////////
	/// Read input parameters ///
	/////////////////////////////
	std::ostringstream os;
/// Get particle properties
	const mass mLiquidHoldup = m_holdupLiquid->GetPhaseMass(_time, EPhase::LIQUID);
	const moistureContent initX = mLiquidHoldup / mSolidHoldup; // initial moisture content of holdup particle
	const massFraction x_wInit = ConvertMoistContentToMassFrac(initX); // initial water mass fraction in particle
/// PSD, currently not in use
	/*
	/// Get initial mass in holdup ///
	m_initMass = m_holdup->GetMass(_time);

	/// Calculate initial diameter ratio for calculating G ///
	m_diamRatio.clear();
	m_diamRatio.push_back(0);
	for (size_t i = 1; i < m_classesNum; ++i)
	{
		m_diamRatio.push_back(std::pow((m_sizeGrid[i] + m_sizeGrid[i + 1]) / (m_sizeGrid[i - 1] + m_sizeGrid[i]), 3)); //at least 3 size grids, i.e. 2 classNums
	}
	*/
	const area A_Calculated = CalculateParticleSurfaceArea(_time);//GetSpecificSurface(Grid, q_3) * (mSolidHoldup / rhoParticle);
	//bool updateA = GetCheckboxParameterValue("updateA");
	const area A_Input = GetConstRealParameterValue("A_P");
	double A_P = 0;
	if (A_Input == 0)
	{
		A_P = A_Calculated;
	}
	else
	{
		A_P = A_Input;
	}
	const length Delta_f = GetConstRealParameterValue("Delta_f") * 1e-6; // convert into [m]		

	//if (particlesGlobal) // set particle moisture for non-empty run
	//{
	//	double w_eqMinTemperature = GetConstRealParameterValue("w_l,eq,min temperature") + T_ref;
	//	double X_eqMinBase = CalcuateSolidEquilibriumMoistureContent(_time, w_eqMinTemperature, GetRelativeHumidity(Y_inGas, w_eqMinTemperature));
	//	double w_leqmin = GetConstRealParameterValue("w_l,eq,min") / 100;
	//	//minMoistureContent = w_leqmin / (1. - w_leqmin);//0.0271;
	//	//moistureScaler = minMoistureContent / X_eqMinBase;
	//}
	// Show particle properties in DEBUG mode

/// environment temperature
	T_env = T_ref + GetConstRealParameterValue("theta_env"); // T_ref = 0 degreeC = 273.15 K

/// inlet fluidization gas conditions
	const temperature T_inGas = m_inGasStream->GetTemperature(_time);
	const temperature theta_inGas = m_inGasStream->GetTemperature(_time) - T_ref;
	const moistureContent y_inGas = m_inGasStream->GetCompoundFraction(_time, keyVapor);
	const massFraction Y_inGas = ConvertMassFracToMoistContent(y_inGas);
	const double RH_inGas = CalculateGasRelativeHumidity(Y_inGas, theta_inGas, m_inGasStream->GetPressure(_time));
	const massFlow mFlowInGas = m_inGasStream->GetMassFlow(_time);
	const massFlow mFlowInGasDry = m_inGasStream->GetCompoundMassFlow(_time, keyGas);
	const specificLatentHeat h_inGas = C_PGas * theta_inGas + Y_inGas * (C_PWaterVapor * theta_inGas + Delta_h0);

/// nozzle gas condition
	const massFlow mFlowInNozzleGas = m_inNozzleAirStream->GetMassFlow(_time);
	const moistureContent y_nozzle = m_inNozzleAirStream->GetCompoundFraction(_time, keyVapor);
	const massFraction Y_nozzle = ConvertMassFracToMoistContent(y_nozzle);
	const massFlow mFlowInNozzleGasDry = m_inNozzleAirStream->GetCompoundMassFlow(_time, keyGas);
	const temperature T_nozzleGas = m_inNozzleAirStream->GetTemperature(_time);
	const temperature thetaNozzleGas = m_inNozzleAirStream->GetTemperature(_time) - T_ref; // in degreeC
	const specificLatentHeat h_nozzleGas = C_PGas * thetaNozzleGas + Y_nozzle * (C_PWaterVapor * thetaNozzleGas + Delta_h0);

/// Spray liquid condition
	const massFlow mFlowSprayLiquid = m_inLiquidStream->GetMassFlow(_time);
	const massFraction x_wSusp = m_inLiquidStream->GetPhaseFraction(_time, EPhase::LIQUID);
	const temperature T_sprayLiquid = m_inLiquidStream->GetTemperature(_time);
	const temperature thetaSprayLiquid = m_inLiquidStream->GetTemperature(_time) - T_ref; // in degreeC
	const specificLatentHeat h_susp = thetaSprayLiquid * (C_PParticle * (1 - x_wSusp) + C_PWaterLiquid * x_wSusp); // for granulation: C_PSolid should be for coating material

/// get process chamber properties
	/*heightOfChamberTemperatureProbe = GetConstRealParameterValue("H_tempProbe");
	heightOfNozzle = GetConstRealParameterValue("H_nozzle");
	heightOfChamber = GetConstRealParameterValue("H_chamber");
	diamOfBed = GetConstRealParameterValue("d_bed");
	SetupChamber();*/
	//this->CheckHeightDiscretizationLayers(_time);

/// get bed properties
	const length heightOfFixedBed = GetConstRealParameterValue("H_bedFix");
	const double eps_0 = GetConstRealParameterValue("eps_0");
	const double u_mf = GetConstRealParameterValue("u_mf");
	//suspLayer = DetermineLayersInSectionFilledWithBed(0, heightOfNozzle / chamber.at(0).height) - 1;

	//////////////////////////////////////////////////////////////////////////////////
	// Calculate gas mass from initial conditions and volume // CURRENLY NOT IN USE //
	//volume volumeSolids = m_holdup->GetPhaseMass(_time, EPhase::SOLID) / rhoParticle;
	//volume volumeChamber = 0;
	//for (size_t i = 0; i < chamber.size(); i++)
	//{
	//	volumeChamber += CalculateSectionVolume(i);
	//}		
	//volume volumeGas = volumeChamber - volumeSolids;
	//mass massGas = volumeGas * rhoGas;
	//mass orignalMassGas = m_holdup->GetPhaseMass(_time, EPhase::GAS);
	//if (massGas > orignalMassGas)
	//{
	//	m_holdup->SetPhaseMass(_time, EPhase::GAS, massGas);
	//	os << "Adjusted mass of gas phase in holdup from " << orignalMassGas << " kg to " << massGas << " kg.";
	//	ShowInfo(os.str());
	//	os.str("");
	//}
	//////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////
	// parameters about drying kinetics, currently not in use //
	//eqData.compoundKey = GetCompoundParameterValue("Xeq compound");
	//if (particlesGlobal)
	//{
	//	InitializeMoistureContentDatabase(GetStringParameterValue("Path Xeq")); //"E:\\Dyssol\\Xeq.csv"
	//}
	////phi_eq = GetConstRealParameterValue("phi_eq");
	//k_dc = GetConstRealParameterValue("k_dc");
	//X_cr = GetConstRealParameterValue("X_cr");
	////////////////////////////////////////////////////////////

	//
	//double compoundFractionOfSolidInLiquidInput = m_inSuspensionStream->GetPhase(EPhase::SOLID)->GetCompoundFraction(_time, GetCompoundIndex(compoundKeys[indicesOfVaporOfPhaseChangingCompound.second]));
	//if (compoundFractionOfSolidInLiquidInput > 0 && !ignoreVaporInput)
	//{
	//	w_l = compoundFractionOfVaporOfPhaseChangingCompound / (1 - compoundFractionOfVaporOfPhaseChangingCompound);
	//	ShowInfo("Y_in was calculated using inGasStream as: " + std::to_string(Y_in) + "kg/kg");
	//}
	//else
	//	w_l = 1;*/
	//mass massGas2 = 0;
	//for (size_t i = 0; i < N_total; i++)
	//	massGas2 += GetGasMassOfLayer(_time, i, STANDARD_CONDITION_T, STANDARD_CONDITION_T);

	/// Clear all state variables in model ///
	m_model.ClearVariables();

	///////////////////////////
	/// Add state variables ///
	///////////////////////////	
/// gas in holdup ///
	const mass mGasHoldup = m_holdupGas->GetMass(_time);
	const mass mGasHoldupDry = m_holdupGas->GetCompoundMass(_time, keyGas);
	const temperature T_holdupGas = m_holdupGas->GetTemperature(_time);
	const temperature theta_holdupGas = T_holdupGas - T_ref;
	const moistureContent y_holdupGas = m_holdupGas->GetCompoundFraction(_time, keyVapor);
	const moistureContent Y_holdupGas = ConvertMassFracToMoistContent(y_holdupGas);
	const double RH_holdupGas = CalculateGasRelativeHumidity(Y_holdupGas, theta_holdupGas, m_holdupGas->GetPressure(_time));
	const pressure pressureGasHoldup = m_holdupGas->GetPressure(_time);
	AddStateVariable("Holdup air mass [kg]", mGasHoldupDry);
	AddStateVariable("Holdup water vapor mass [kg]", mGasHoldup * y_holdupGas);
/// gas in outlet ///
	m_model.m_iYOutGas = m_model.AddDAEVariable(true, Y_holdupGas, 0.0, 0.0); // initial outlet gas absolut humidity == holdup gas 
	m_model.m_iTempOutGas = m_model.AddDAEVariable(true, T_holdupGas, 0.0, 0.0); // initial outlet gas temperature == holdup gas 
	AddStateVariable("Outlet gas Y [g/kg]", Y_holdupGas * 1e3);
	AddStateVariable("Outlet gas RH [%]", RH_holdupGas * 100);

/// particle ///
	temperature T_ParticleInit = 0;
	temperature theta_holdupSolid = 0;
	temperature T_LiquidInit = 0;
	temperature theta_holdupLiquid = 0;
	T_ParticleInit = m_holdupSolid->GetTemperature(_time); // in [K]
	theta_holdupSolid = T_ParticleInit - T_ref;
	T_LiquidInit = m_holdupLiquid->GetTemperature(_time); 
	theta_holdupLiquid = T_LiquidInit - T_ref;
	double initPhi = mLiquidHoldup / rhoWater / Delta_f / A_P;
	m_model.m_iTempParticle = m_model.AddDAEVariable(true, T_ParticleInit, 0.0, 0.0); // Particle temperature in [K]
	m_model.m_iPhi = m_model.AddDAEVariable(true, initPhi, 0.0, 0.0); 
	// m_model.miA_P: A_P is constant in case of water spray, A_P as DAE variable will be used for granulation
	AddStateVariable("Particle moisture content [%]", initX * 100);
	AddStateVariable("Particle water mass fraction [%]", x_wInit * 100);
	AddStateVariable("Particle wetness degree [%]", initPhi * 100);
	const bool use_f = GetCheckboxParameterValue("Use relative drying rate?");
	const double fInit = use_f ? CalculateRelativeDryingRate(initX) : 1;
	const size_t Y_maxMethod = GetComboParameterValue("Use Y_eq or Y_sat?");
	const moistureContent Y_satInit = CalculateGasSaturationMoistureContent(theta_holdupGas, pressureGasHoldup);
	const double RH_eqInit = CalculateGasEquilibriumRelativeHumidity(initX);
	const pressure P_satInit = CalculateGasSaturationPressure(theta_holdupGas, pressureGasHoldup);
	const moistureContent Y_eqInit = CalculateGasEquilibriumMoistureContent(pressureGasHoldup, P_satInit, RH_eqInit);
	double Y_maxInit = 0;
	if (Y_maxMethod == 0) // use Y_eq
	{
		Y_maxInit = Y_eqInit;
	}
	else // use Y_sat
	{
		Y_maxInit = Y_satInit;
	}
	double betaInit = CalculateBeta(_time, d32, theta_holdupGas, CalculateDiffusionCoefficient(theta_holdupGas));
	double MFlowInit = betaInit * A_P * initPhi * rhoGas * (Y_maxInit - Y_holdupGas) * fInit;
	AddStateVariable("Water evaporation rate [g/s]", MFlowInit);
	AddStateVariable("Relative drying rate [-]", fInit);
	AddStateVariable("Heat loss gas to wall [J/s]", 0);
		
	os << "\nInitial particle moisture content: " << initX * 1e3 << " g/kg dry solid \nWater mass fraction: " << x_wInit * 100 << " %\n";
	os << "Initial particle wetness degree: " << initPhi * 100 << " %\n";
	ShowInfo(os.str());
	os.str("");
		
/// liquid (film) properties ///
	m_model.m_iTempFilm = m_model.AddDAEVariable(true, T_LiquidInit, 0.0, 0.0); 

///  Add plots ///
	// temperatures
	AddPlot("TEMPERATURES", "Time [s]", "Temperature [degC]");
	AddCurveOnPlot("TEMPERATURES", "Holdup solid [degC]");
	AddCurveOnPlot("TEMPERATURES", "Holdup liquid film [degC]");
	AddCurveOnPlot("TEMPERATURES", "Holdup gas [degC]");
	AddCurveOnPlot("TEMPERATURES", "Outlet gas [degC]");
	// water vapor mass or mass flow 
	AddPlot("MASS water vapor", "Time [s]", "Mass flow [kg/s] or mass [kg]");
	AddCurveOnPlot("MASS water vapor", "Inlet vapor mass flow [kg/s]");
	AddCurveOnPlot("MASS water vapor", "Outlet vapor mass flow [kg/s]");
	AddCurveOnPlot("MASS water vapor", "Vapor mass in holdup [kg]");
	// water liquid mass or mass flow 
	AddPlot("MASS water liquid", "Time [s]", "Mass flow [kg/s] or mass [kg]");
	AddCurveOnPlot("MASS water liquid", "Inlet liquid mass flow [kg/s]");
	AddCurveOnPlot("MASS water liquid", "Outlet liquid mass flow [kg/s]");
	AddCurveOnPlot("MASS water liquid", "Liquid mass in holdup [kg]");
	// energy streams
	AddPlot("HEAT TRANSFER", "Time [s]", "Enthalpy flow [J/s] or coefficients");
	AddCurveOnPlot("HEAT TRANSFER", "Heat transfer GP [J/s]");
	AddCurveOnPlot("HEAT TRANSFER", "alpha_GP [W/(m2*K)]");
	AddCurveOnPlot("HEAT TRANSFER", "Heat transfer GF [J/s]");
	AddCurveOnPlot("HEAT TRANSFER", "alpha_GF [W/(m2*K)]");
	AddCurveOnPlot("HEAT TRANSFER", "Heat transfer PF [J/s]");
	AddCurveOnPlot("HEAT TRANSFER", "alpha_PF [W/(m2*K)]");
	AddCurveOnPlot("HEAT TRANSFER", "Water vapor flow enthalpy [J/s]");
	AddCurveOnPlot("HEAT TRANSFER", "beta [m/s]");
	AddCurveOnPlot("HEAT TRANSFER", "Heat loss gas to wall [J/s]");
	// total energy flow / holdup
	AddPlot("ENERGY total", "Time [s]", "Enthalpy flow [J/s] or Enthalpy [J]");
	AddCurveOnPlot("ENERGY total", "Inlet enthalpy flow [J/s]");
	AddCurveOnPlot("ENERGY total", "Outlet enthalpy flow [J/s]");
	AddCurveOnPlot("ENERGY total", "Enthalpy in holdup [J]");
	// water vapor energy flow / holdup 
	AddPlot("ENERGY water vapor", "Time [s]", "Enthalpy flow [J/s] or Enthalpy [J]");
	AddCurveOnPlot("ENERGY water vapor", "Inlet enthalpy flow [J/s]");
	AddCurveOnPlot("ENERGY water vapor", "Outlet enthalpy flow [J/s]");
	AddCurveOnPlot("ENERGY water vapor", "Enthalpy in holdup [J]");
	// water liquid energy flow / holdup 
	AddPlot("ENERGY water liquid", "Time [s]", "Enthalpy flow [J/s] or Enthalpy [J]");
	AddCurveOnPlot("ENERGY water liquid", "Inlet enthalpy flow [J/s]");
	AddCurveOnPlot("ENERGY water liquid", "Outlet enthalpy flow [J/s]");
	AddCurveOnPlot("ENERGY water liquid", "Enthalpy in holdup [J]");
	// dry air energy flow / holdup
	AddPlot("ENERGY dry air", "Time [s]", "Enthalpy flow [J/s] or Enthalpy [J]");
	AddCurveOnPlot("ENERGY dry air", "Inlet enthalpy flow [J/s]");
	AddCurveOnPlot("ENERGY dry air", "Outlet enthalpy flow [J/s]");
	AddCurveOnPlot("ENERGY dry air", "Enthalpy in holdup [J]");
	// solid Cellets energy flow / holdup
	AddPlot("ENERGY particle", "Time [s]", "Enthalpy flow [J/s] or Enthalpy [J]");
	AddCurveOnPlot("ENERGY particle", "Inlet enthalpy flow [J/s]");
	AddCurveOnPlot("ENERGY particle", "Outlet enthalpy flow [J/s]");
	AddCurveOnPlot("ENERGY particle", "Enthalpy in holdup [J]");

	/// height discretization, CURRENTLY NOT IN USE
	//std::vector<moistureContent> Y_inInit(N_total, Y_in); // vector in case of height discretization
	//std::vector<temperature> temperatureGasInit(N_total, m_holdup->GetTemperature(_time) - T_ref); // vector in case of height discretization
	//for (int i = 0; i < Y_init.size(); i++)
	//{
	//	AddStateVariable("Gas temperature_i [degreeC] " + std::to_string(i), m_holdup->GetTemperature(_time) - T_ref); // [°C]
	//	AddStateVariable("Y_i [g/kg] " + std::to_string(i), Y_in) * 1e3; // [g/kg]
	//}
	//AddStateVariable("Avg. Y [g/kg]", Y_in * 1e3); //  [g/kg]
	//ShowInfo("RH to Y: " + std::to_string(RH_in * GetGasSaturationMoistureContent(T_ref + 20.5)*1000) + " g/kg");
	//ShowInfo("Massflow nozzle: " + std::to_string((5.91 /*Volume flow gas m^3/h*/ * rhoGas / 60 * 1000 + 20 /*Mass flow water g/min*/) / 60 / 1000) + " kg/s");
	//ShowInfo("Phase fraction water: " + std::to_string(20. /*Mass flow water g/min*/ / 60 / 1000 / ((5.91 /*Volume flow gas m^3/h*/ * rhoGas / 60 * 1000 + 20 /*Mass flow water g/min*/) / 60 / 1000)));
	//temperature T_out = T_ref + 48.2;
	//ShowInfo("RH - in vs out: " + std::to_string(RH_in*100) + " at " + std::to_string(T_inf-T_ref) + " deg C turn into " + std::to_string(100 * GetRelativeHumidity(RH_in * GetGasSaturationMoistureContent(T_inf), T_out)) + " at " + std::to_string(T_out-T_ref) + " deg C.");
	//ShowInfo(std::to_string(_y / GetGasSaturationMoistureContent(_t, _p)));
		
	// Check mass and energy balance in DEBUG mode
	/*if (debugToggle)
	{
		AddStateVariable("Mass balance", 0);
		AddStateVariable("Energy balance", 0);
		for (int i = 0; i < 10; i++)
		{
			AddStateVariable("debug" + std::to_string(i), 0);
		}			
	} */

	/// checkForSmallBiot
	//if (particlesGlobal)
	//	if (!CheckForSmallBiot(_time))
	//		RaiseWarning("Biot number: "+std::to_string(CalcBiotNumber(_time)) + ">" + std::to_string(SmallBiotNumber) + "\nAssumption of uniform temperature distribution throughout material is invalid.");
	
	/*if (DiffCoefWarning) {
		RaiseWarning("Critical eqData.temperatures and critical specific molar volumes used, are hard coded for air and water.");
		DiffCoefWarning = false;
	}*/

	//m_model.InitCounterVariables();
	//ShowInfo(std::to_string(m_model.CalculateBetaPA(_time, m_holdup->GetTemperature(_time), m_holdup->GetTemperature(_time)) / m_model.CalculateBetaPA(_time, m_holdup->GetTemperature(_time), m_holdup->GetTemperature(_time) + 4)*100));
	//debugHasBeenShown = false;
	//DiffCoefWarning = true;

	//m_model.unit = static_cast<CDryerBatch*>(this);
	//if (particlesGlobal)
	//{
	//	this->EnergyLiquidPhaseOld = this->m_holdup->GetPhaseMass(_time, EPhase::LIQUID) * this->C_PWaterLiquid * T_LiquidInit;
	//	this->EnergySolidPhaseOld = this->m_holdup->GetPhaseMass(_time, EPhase::SOLID) * this->C_PParticle * T_ParticleInit;
	//}
	//this->EnergyGasPhaseOld = this->m_holdup->GetPhaseMass(_time, EPhase::GAS) * ((this->C_PGas + this->C_PWaterVapor * Y_inGas) * m_holdup->GetTemperature(_time) + this->Delta_h0 * Y_inGas);

	//GetStringParameterValue("Path")

	/// Set tolerances to model ///
	//const auto rtol = GetConstRealParameterValue("Relative tolerance");
	//// separate absolute tolerance for temperature and gas moisture content
	//const auto atolT = GetConstRealParameterValue("Absolute tolerance T");
	//const auto atolY = GetConstRealParameterValue("Absolute tolerance Y");
	//std::vector<double> absolutTolerances;

	/// height discretization
	//for (int i = 0; i < Y_inInit.size(); i++)
	//	absolutTolerances.push_back(atolY != 0.0 ? atolY : GetAbsTolerance()); // 0.0001
	//for (int i = 0; i < temperatureGasInit.size() + (particlesGlobal ? 2 : 0); i++)
	//	absolutTolerances.push_back(atolT != 0.0 ? atolT : GetAbsTolerance()); // 0.01
	//size_t variables = m_model.GetVariablesNumber();
	//m_model.SetTolerance(rtol != 0.0 ? rtol : GetRelTolerance(), rtol != 0.0 ? rtol : GetRelTolerance());
	//m_model.SetTolerance(rtol != 0.0 ? rtol : GetRelTolerance(), absolutTolerances); // 0.01
	m_model.SetTolerance(GetRelTolerance(), GetAbsTolerance());
	//m_solver.SetMaxStep(1);

	/// Set model to a solver ///
	if (!m_solver.SetModel(&m_model))
		RaiseError(m_solver.GetError());

	if (printResult)
	{	
		/// Get number of classes for PSD ///
		size_t m_classesNum = GetClassesNumber(DISTR_SIZE); //n
		/// Get grid of PSD ///
		const std::vector<double> avgClassDiam = GetClassesMeans(DISTR_SIZE); //d_m,i
		const std::vector<double> classSize = GetClassesSizes(DISTR_SIZE); //Delta d
		// output of read data in simulation window
		os << "PSD info\n" << "\tNumber of classes: " << m_classesNum << "\n\tSize Grid\n";
		for (int i = 0; i < Grid.size(); i++)
			os << "\t\tSize Grid: " << i << " = " << Grid[i] << " m\n";
		os << "\tAvg diam\n";
		for (int i = 0; i < avgClassDiam.size(); i++)
			os << "\t\tAvg diam: " << i << " = " << avgClassDiam[i] << " m\n";
		os << "\tClass size\n";
		for (int i = 0; i < classSize.size(); i++)
			os << "\t\tClass size: " << i << " = " << classSize[i] << " m\n";
		os << "\tSauter diameter:" << d32 << " m\n";
		os << "\tToal surface area: " << A_Calculated << " m^2\n" << "\tTotal particle mass: " << mSolidHoldup << " kg\n";
		if (particlesGlobal)
		{
			ShowInfo(os.str());
		}
		os.str("");
		/// Phase mass ///
		os << "Phase mass:\n\t Solid: " << m_holdupSolid->GetPhaseMass(_time, EPhase::SOLID) << " kg";
		os << "\n\t Liquid: " << m_holdupLiquid->GetPhaseMass(_time, EPhase::LIQUID) << " kg";
		os << "\n\t Gas: " << m_holdupGas->GetPhaseMass(_time, EPhase::VAPOR) << " kg\n";
		ShowInfo(os.str());
		os.str("");
	}

	os << "Initialization completed.";
	ShowInfo(os.str());
}

void CDryerBatch::SaveState()
{
	m_solver.SaveState();
	ShowInfo("State saved.");
}

void CDryerBatch::LoadState()
{
	m_solver.LoadState();
	ShowInfo("State loaded.");
}

void CDryerBatch::Simulate(double _timeBeg, double _timeEnd)
{
	ShowInfo("Start simulation.");
	if (!m_solver.Calculate(_timeBeg, _timeEnd))
		RaiseError(m_solver.GetError());
}


///////////////////////////////////
///  Main calculation functions ///
///////////////////////////////////
void CUnitDAEModel::CalculateResiduals(double _time, double* _vars, double* _ders, double* _res, void* _unit)
{
/// Define unit, streams and holdup - data must be read-only for calculating residual, set as const ///
	const auto* unit = static_cast<CDryerBatch*>(_unit);
	const CStream* inGasStream = unit->GetPortStream("InletFluidizationGas");
	const CStream* inLiquidStream = unit->GetPortStream("InletLiquid");
	const CStream* inNozzleAirStream = unit->GetPortStream("InletNozzleAir");
	const CHoldup* holdupSolid = unit->GetHoldup("HoldupSolid");
	const CHoldup* holdupLiquid = unit->GetHoldup("HoldupLiquid");
	const CHoldup* holdupGas = unit->GetHoldup("HoldupGas");
	
/// Read input parameters ///
	const bool printResult = unit->GetCheckboxParameterValue("Print intermediate results");
	/// Phase properties as constant
	const heatCapacity C_PGas = unit->C_PGas;
	const thermalConductivity lambdaGas = unit->lambdaGas;
	const heatCapacity C_PWaterLiquid = unit->C_PWaterLiquid;
	const heatCapacity C_PWaterVapor = unit->C_PWaterVapor;
	const specificLatentHeat Delta_h0 = unit->Delta_h0;
	const heatCapacity C_PParticle = unit->C_PParticle;
	const density rhoGas = unit->rhoGas;
	const density rhoVapor = unit->rhoVapor;
	const density rhoLiquid = unit->rhoWater;
	const density rhoParticle = unit->rhoParticle;
	/// Particle in holdup
	const mass mHoldupSolid = holdupSolid->GetPhaseMass(_time, EPhase::SOLID); // Solid mass holdup [kg] - constant for fluidization with water
	const mass mHoldupLiquid = holdupLiquid->GetPhaseMass(_time, EPhase::LIQUID);	
	const moistureContent X_wP = mHoldupLiquid / mHoldupSolid;
	const massFraction x_wP = unit->ConvertMoistContentToMassFrac(X_wP);
	const length Delta_f = unit->GetConstRealParameterValue("Delta_f") * 1e-6; // convert in [m]
	const area A_P = unit->CalculateParticleSurfaceArea(_time);
	//const double derA_p = 0; // Time change of total particle surface area, constant for fluidization with water
	const length d32 = unit->CalculateHoldupSauter(_time);
	/// Hydrodynamics
	const double eps_0 = unit->GetConstRealParameterValue("eps_0");
	const double u_Gas = unit->CalculateGasVel(_time, d32); 
	/// Liquid in holdup
	const temperature T_holdupLiquid = holdupLiquid->GetTemperature(_time); // in [K]
	const temperature theta_holdupLiquid = T_holdupLiquid - unit->T_ref; // in degreeC
	/// Inlet fluidization gas
	const massFraction y_inGas = inGasStream->GetCompoundFraction(_time, unit->keyVapor);
	const moistureContent Y_inGas = unit->ConvertMassFracToMoistContent(y_inGas);
	const massFlow mFlowInGas = inGasStream->GetMassFlow(_time); // Gas mass flow [kg/s]
	const massFlow mFlowInGasDry = inGasStream->GetCompoundMassFlow(_time, unit->keyGas);
	const temperature T_inGas = inGasStream->GetTemperature(_time); // temperature in [K]
	const temperature theta_inGas = T_inGas - unit->T_ref; // temperature in [degreeC]
	const specificLatentHeat h_inGas = C_PGas * theta_inGas + Y_inGas * (C_PWaterVapor * theta_inGas + Delta_h0);
	/// Inlet nozzle gas
	const massFlow mFlowInNozzleGas = inNozzleAirStream->GetMassFlow(_time);
	const massFraction y_nozzle = inNozzleAirStream->GetCompoundFraction(_time, unit->keyVapor);
	const moistureContent Y_nozzle = unit->ConvertMassFracToMoistContent(y_nozzle);
	const massFlow mFlowInNozzleGasDry = inNozzleAirStream->GetCompoundMassFlow(_time, unit->keyGas);
	const temperature T_nozzleGas = inNozzleAirStream->GetTemperature(_time); // in [K]
	const temperature thetaNozzleGas = T_nozzleGas - unit->T_ref;
	const specificLatentHeat h_nozzleGas = C_PGas * thetaNozzleGas + Y_nozzle * (C_PWaterVapor * thetaNozzleGas + Delta_h0);
	/// Spray liquid
	const massFlow mFlowSprayLiquid = inLiquidStream->GetMassFlow(_time);
	const massFraction x_wSusp = inLiquidStream->GetPhaseFraction(_time, EPhase::LIQUID);
	const temperature T_sprayLiquid = inLiquidStream->GetTemperature(_time);
	const temperature thetaSprayLiquid = inLiquidStream->GetTemperature(_time) - unit->T_ref;
	const specificLatentHeat h_susp = thetaSprayLiquid * (C_PParticle /*coating material*/ * (1 - x_wSusp) + C_PWaterLiquid * x_wSusp);
	
/// DAE system: 5 equations (all ture) ///
	// _vars: determined by the solver & should not be changed, set as const!
	// gas in holdup
	const mass mHoldupGas = holdupGas->GetMass(_time);
	const mass mHoldupGasDry = holdupGas->GetCompoundMass(0, unit->keyGas); // set intial dry air mass as dry air mass during whole process
	const temperature varT_gasHoldup = holdupGas->GetTemperature(_time);
	const temperature varTheta_gasHoldup = varT_gasHoldup - unit->T_ref;
	const moistureContent varYHoldupGas = _vars[m_iYOutGas]; // current simplification without descretization: Y_holdup == Y_out
	const pressure pressureGasHoldup = holdupGas->GetPressure(_time); // Pressure holdup [Pa]	
	// gas outlet
	const moistureContent varYOutGas = _vars[m_iYOutGas];
	const temperature varTempOutGas = _vars[m_iTempOutGas];
	const temperature varThetaOutGas = varTempOutGas - unit->T_ref;
	const double varHFlowOutGasFormula = (mFlowInGasDry + mFlowInNozzleGasDry) * (C_PGas * varThetaOutGas + varYOutGas * (C_PWaterVapor * varThetaOutGas + Delta_h0));
	const double varP_sat_Formula = unit->CalculateGasSaturationPressure(varTheta_gasHoldup, pressureGasHoldup);
	// particle 
	const temperature varTempParticle = _vars[m_iTempParticle];
	const temperature varThetaParticle = varTempParticle - unit->T_ref;
	const double varPhi = _vars[m_iPhi];
	const moistureContent varX = varPhi * A_P * Delta_f * unit->rhoWater / mHoldupSolid;
	// liquid film
	const temperature varTempFlim = _vars[m_iTempFilm];
	const temperature varThetaFilm = varTempFlim - unit->T_ref;
	const moistureContent varY_sat_Formula = unit->CalculateGasSaturationMoistureContent(varThetaFilm, pressureGasHoldup);
	const double varRH_eq_Formula = unit->CalculateGasEquilibriumRelativeHumidity(varX);
	const moistureContent varY_eq_Formula = unit->CalculateGasEquilibriumMoistureContent(pressureGasHoldup, varP_sat_Formula, varRH_eq_Formula);
	// Heat transfer
	const double varAlpha_GP_Formula = unit->CalculateAlpha_GP(_time, varTheta_gasHoldup, d32);
	const double varAlpha_GF_Formula = unit->CalculateAlpha_GF(_time, varTheta_gasHoldup, d32);
	const double varAlpha_PF_Formula = unit->CalculateAlpha_PF(/*varTempFlim, pressureGasHoldup, d32*/ varAlpha_GP_Formula);
	// Mass transfer
	const double varD_a_Formula = unit->CalculateDiffusionCoefficient(varTheta_gasHoldup);
	const double shrinkFactorBeta = unit->GetConstRealParameterValue("Shrink factor beta");
	const double varBeta_FG_Formula = _time < 10? shrinkFactorBeta * unit->CalculateBeta(_time, d32, varThetaOutGas, varD_a_Formula) : unit->CalculateBeta(_time, d32, varThetaOutGas, varD_a_Formula);
	// water vapor
	const bool use_f = unit->GetCheckboxParameterValue("Use relative drying rate?");
	const double f = use_f ? unit->CalculateRelativeDryingRate(varX) : 1;
	const size_t Y_maxMethod = unit->GetComboParameterValue("Use Y_eq or Y_sat?");
	double Y_max = 0;
	if (Y_maxMethod == 0) // use Y_eq
	{
		Y_max = varY_eq_Formula;
	}
	else // use Y_sat
	{
		Y_max = varY_sat_Formula;
	}
	double varMFlowVaporFormula = varBeta_FG_Formula * A_P * varPhi * rhoGas * (Y_max - varYHoldupGas) * f; 
	const double varHFlowVaporFormula = varMFlowVaporFormula * (C_PWaterVapor * varTheta_gasHoldup + Delta_h0);
	// heat flow
	const double tolTemp = unit->GetConstRealParameterValue("Tolerance temperature");
	const double shrinkFactorGP = unit->GetConstRealParameterValue("Shrink factor alpha GP");
	const double shrinkFactorGF = unit->GetConstRealParameterValue("Shrink factor alpha GF");
	double varAlpha_GP = varAlpha_GP_Formula;
	double varAlpha_GF = varAlpha_GF_Formula;
	double varAlpha_PF = varAlpha_PF_Formula;
	double varAlpha_GP_modify = shrinkFactorGP * varAlpha_GP_Formula;
	double varAlpha_GF_modify = shrinkFactorGF * varAlpha_GF_Formula;
	double varAlpha_PF_modify = shrinkFactorGF * varAlpha_PF_Formula;
	double diffTempGP = varTheta_gasHoldup - varThetaParticle;
	double diffTempGF = varTheta_gasHoldup - varThetaFilm;
	double diffTempPF = varThetaParticle - varThetaFilm;
	double varQFlow_GP_Formula = 0;
	double varQFlow_GF_Formula = 0;
	double varQFlow_PF_Formula = 0;
	/// TODO: make the tolerance range differenciable ///
	/*
	modify diffTemp: so that its change is not sharp (jumping between tol and 0)
	using the tanh function
	idea: diffTemp = tol * fun + zero * (1 - fun) -> see Vasyl's functions	
	*/

	if (varPhi < 1)
	{
		if (abs(diffTempGP) < tolTemp)
		{
			varAlpha_GP = varAlpha_GP_modify;			
		}		
		varQFlow_GP_Formula = varAlpha_GP * A_P * (1 - varPhi) * diffTempGP;
	}
	if (abs(diffTempGF) < tolTemp)
	{
		varAlpha_GF = varAlpha_GF_modify;		
	}
	varQFlow_GF_Formula = varAlpha_GF * A_P * varPhi * diffTempGF;
	if (abs(diffTempPF) < tolTemp)
	{
		varAlpha_PF = varAlpha_PF_modify;	
	}
	varQFlow_PF_Formula = varAlpha_PF * A_P * varPhi * diffTempPF;
	temperature T_surface = unit->IterateSurfaceTemp(_time, varT_gasHoldup, d32);
	const double QFlow_GW_Chamber = unit->CalculateHeatLossWall(_time, unit->wallThickness, unit->GetConstRealParameterValue("H_plant"), unit->GetConstRealParameterValue("d_bed"), varTheta_gasHoldup, T_surface, unit->lambdaWall, d32);

	/// _ders: determined by the solver & should not be changed, set as const!
	// gas phase (outlet gas)
	const double derTempOutGas = _ders[m_iTempOutGas];
	const double derYOutGas = _ders[m_iYOutGas];
	const double derTempOutGasFormula = (mFlowInGasDry * h_inGas + mFlowInNozzleGasDry * h_nozzleGas - varHFlowOutGasFormula - varQFlow_GP_Formula - varQFlow_GF_Formula - QFlow_GW_Chamber + varHFlowVaporFormula) / (mHoldupGasDry * (C_PGas + C_PWaterVapor * varYOutGas)) - derYOutGas * (C_PWaterVapor * varThetaOutGas + Delta_h0) / (C_PGas + C_PWaterVapor * varYOutGas);
	const double derYOutGasFormula = (mFlowInGasDry * Y_inGas + mFlowInNozzleGasDry * Y_nozzle - (mFlowInGasDry + mFlowInNozzleGasDry) * varYOutGas) / mHoldupGasDry + varMFlowVaporFormula / mHoldupGasDry; // currently Y_in - Y_out
	// particle (solid) phase
	const double derTempParticle = _ders[m_iTempParticle];
	const double derTempParticleFormula = (varQFlow_GP_Formula - varQFlow_PF_Formula /* - varQ_PW_Formula in the future*/) / (mHoldupSolid * C_PParticle); 
	// liquid phase (water film)
	const double derTempFlim = _ders[m_iTempFilm];
	const double derPhi = _ders[m_iPhi];
	const double derTempFilmFormula = (varQFlow_PF_Formula + varQFlow_GF_Formula + mFlowSprayLiquid * h_susp - varHFlowVaporFormula) / (C_PWaterLiquid * Delta_f * rhoLiquid * A_P * varPhi) - varThetaFilm * derPhi / varPhi;
	const double derPhiFormula = (mFlowSprayLiquid - varMFlowVaporFormula /** (1 + varMFlowVaporFormula * (1 - x_wSusp) / (mFlowSprayLiquid * x_wSusp))*/) / (rhoLiquid * Delta_f * A_P); // need to be adjusted in case of granulation

	/// define _res
	// outlet gas
	_res[m_iYOutGas] = derYOutGas - derYOutGasFormula;
	_res[m_iTempOutGas] = derTempOutGas - derTempOutGasFormula;
	// particle (solid) phase
	_res[m_iTempParticle] = derTempParticle - derTempParticleFormula;
	_res[m_iPhi] = derPhi - derPhiFormula;
	// liquid phase (water film)
	_res[m_iTempFilm] = derTempFlim - derTempFilmFormula;

/// Codes for checking _vars, _ders and _res ///
	//if (printResult)
	//{
	//	std::vector vars(_vars, _vars + GetVariablesNumber());
	//	std::vector ders(_ders, _ders + GetVariablesNumber());
	//	std::vector  res(_res, _res + GetVariablesNumber());
	//}

	// size_t sectionsWithParticles = unit->particlesGlobal ? std::ceil(unit2->DetermineSectionsFilledWithBed(_time, varTempParticle2)) : 0;
	// size_t sectionsWithParticles = unit->particlesGlobal ? std::ceil(unit2->CalculateBedHeightOrDetermineSectionsFilledWithBed(_time, varTempParticle2, false)) : 0;
	//size_t layersWithParticles = 0;
	//for (size_t section = 0; section < sectionsWithParticles; section++)
	//	layersWithParticles += unit->chamber.at(section).layers;
	// Height averaged temperature of the gas phase [K]
	// Dosta 2010 eq. 14a
	//const temperature varAvTempGas = CalculateAverage(_vars, m_iTempOutGas, unit->N_particle);
	// Height averaged moisture content of the gas phase [kg/kg]
	// Dosta 2010 eq. 14b
	//const moistureContent Y_av = CalculateAverage(_vars, m_iYOutGas, unit->N_particle);
	//double rhoGas = unit->rhoGas; // Average density of gas phase for average gas temperature and system pressure [kg/m^3]
	//rhoGas = unit->GetAvgTPCompoundProperty(_time, EPhase::GAS, ECompoundTPProperties::DENSITY, varAvTempGas, pressureGasHoldup);
	//double rhoLiquid = unit->rhoWater; // Average density of gas phase for average gas temperature and system pressure [kg/m^3]
	//rhoLiquid = unit->GetAvgTPCompoundProperty(_time, EPhase::LIQUID, ECompoundTPProperties::DENSITY, varTempFilm2, pressureGasHoldup);
	//double rhoSolid = unit->rhoParticle; // Average density of solid phase for average particle temperature and system pressure [kg/m^3]
	//rhoSolid = unit->GetAvgTPCompoundProperty(_time, EPhase::SOLID, ECompoundTPProperties::DENSITY, varTempParticle2, pressureGasHoldup);
	//double cPLiquidGaseous = unit->C_PWaterVapor;
	//double cPSolid = unit->C_PParticle;
	//cPSolid = unit->GetAvgTPCompoundProperty(_time, EPhase::SOLID, ECompoundTPProperties::HEAT_CAPACITY_CP, varTempParticle2, pressureGasHoldup);
	//const double D = unit->CalculateDiffusionCoefficient(_time, varAvTempGas, varTempFilm2, pressureGasHoldup);
	//const double beta = !unit->particlesGlobal ? 0 : unit->CalculateBeta(_time, varAvTempGas, varTempParticle2, D); // Mass transfer coefficient for film to gas [m/s]
	//if (beta == -1)
	//	unit2->RaiseError("Mass transfer coefficient calculation error. Provided temperature = nan or < 0 K");
	//double Y_sat = unit2->GetGasSaturationMoistureContent(varAvTempGas, pressureGasHoldup); // Saturation moisture content of the gas (average) [kg/kg]	

	/// for height discretization, CURRENTLY NOT IN USE ///
	//ParallelFor(unit->N_total, [&](size_t i)
	//	{
	//		double HFlowInGasFormula = 0; // Enthaly stream entering height layer Dosta
	//		if (i == 0)
	//			HFlowInGasFormula = mFlowInGas * (unit->C_PGas * T_in + Y_in * (unit->C_PWaterVapor * T_in + unit->Delta_h0));
	//		// Dosta 2010 eq. 18 - addjusted for first layer with gas input
	//		else
	//			HFlowInGasFormula = (mFlowInGas + ((i > unit->suspLayer) ? mFlowSuspensionGas : 0)) * (unit->C_PGas * _vars[m_iTempOutGas + (i - 1)] + _vars[m_iYOutGas + (i - 1)] * (unit->C_PWaterVapor * _vars[m_iTempOutGas + (i - 1)] + unit->Delta_h0));
	//		// Dosta 2010 eq. 18
	//		if (i == unit->suspLayer)
	//			HFlowInGasFormula += mFlowSuspensionGas * (unit->C_PGas * mSuspensionTemperature + Y_in * (unit->C_PWaterVapor * mSuspensionTemperature + unit->Delta_h0));
	//		HFlowInGas2[i] = HFlowInGasFormula;
	//		const double HFlowOutGasFormula = (mFlowInGas + ((i >= unit->suspLayer) ? mFlowSuspensionGas : 0)) * (unit->C_PGas * _vars[m_iTempOutGas + i] + _vars[m_iYOutGas + i] * (unit->C_PWaterVapor * _vars[m_iTempOutGas + i] + unit->Delta_h0)); // Enthaly stream exiting height layer
	//		// Dosta 2010 eq. 19
	//		HFlowOutGas2[i] = HFlowOutGasFormula;
	//		if (i < unit->N_particle)
	//			sumEnergyStreamsGas[i] = HFlowInGas2[i] - HFlowOutGas2[i] - Q_AP[i] - Q_AF[i] - Q_GW[i] + (HFlowVapor2/* + H_gasSupp */) / unit->N_particle;
	//		else
	//			sumEnergyStreamsGas[i] = HFlowInGas2[i] - HFlowOutGas2[i] - Q_GW[i];
	//	});
	//std::vector<double> layerGasMasses = unit2->GetGasMassOfLayers(unit->heighestFlowTimepoint, varAvTempGas, varTempParticle2);
	/// for height discretization, CURRENTLY NOT IN USE ///
	//ParallelFor(unit->N_total, [&](size_t i)
	//	{
	//		std::vector vars(_vars, _vars + GetVariablesNumber());
	//		std::vector ders(_ders, _ders + GetVariablesNumber());
	//		std::vector  res(_res, _res + GetVariablesNumber());
	//		double varYOutFormula2 = 0; // Change of gas moisture content
	//		if (i == 0)
	//			varYOutFormula2 = -(mFlowGas + ((i >= unit->suspLayer) ? mFlowSuspensionGas : 0)) * unit->N_total / mHoldupGas * (_vars[m_iYOutGas + i] - Y_in) + MFlowVapor / mHoldupGas;
	//		// Dosta 2010 eq. 12 addjusted for first layer with gas input
	//		else if (i < unit->N_particle)
	//			varYOutFormula2 = -(mFlowGas + ((i >= unit->suspLayer) ? mFlowSuspensionGas : 0)) * unit->N_total / mHoldupGas * (_vars[m_iYOutGas + i] - _vars[m_iYOutGas + (i - 1)]) + MFlowVapor / mHoldupGas;
	//		// Dosta 2010 eq. 12
	//		else
	//			varYOutFormula2 = -(mFlowGas + ((i >= unit->suspLayer) ? mFlowSuspensionGas : 0)) * unit->N_total / mHoldupGas * (_vars[m_iYOutGas + i] - _vars[m_iYOutGas + (i - 1)]);
	//		if (i == unit->suspLayer)
	//			varYOutFormula2 += mFlowSuspensionGas * Y_in * unit->N_total / mHoldupGas;
	//		_res[m_iYOutGas + i] = _ders[m_iYOutGas + i] - varYOutFormula2;
	//		// Change of gas temperature
	//		const double derTempGasFormula2 = unit->N_total / (mHoldupGas * (unit->C_PGas + unit->C_PWaterVapor * _vars[m_iYOutGas + i])) * sumEnergyStreamsGas[i]
	//			- _ders[m_iYOutGas + i] * ((unit->C_PWaterVapor * _vars[m_iTempOutGas + i] + unit->Delta_h0) / (unit->C_PGas + unit->C_PWaterVapor * _vars[m_iYOutGas + i]));
	//		// Dosta 2010 eq. 10
	//		_res[m_iTempOutGas + i] = _ders[m_iTempOutGas + i] - derTempGasFormula2;
	//	}
	//);
	/// Calculation for granulation with particle size change, CURRENLY NOT IN USE ///
	//const double H_granulation = 0; // granulation enthalpy
	//const double deltaHParticleGrowth = 0; // sum of all variables concerning particle growth
	//const double sumEnergyStreamsParticle = Q_AP3 - Q_PW - Q_PF2 - deltaHParticleGrowth + (H_susp2 - HFlowVapor2) * (1 - varPhi2); // Sum of energy streams to particle
	//const double derTempParticleFormula2 = sumEnergyStreamsParticle / (mHoldupSolid * (unit->C_PParticle + (liquidFilm ? XCutt : X) * unit->C_PWaterLiquid)) + (liquidFilm ? 0 : -varTempParticle2 * unit->C_PWaterLiquid / (unit->C_PParticle + unit->C_PWaterLiquid * X) * ((mFlowSuspension - MFlowVapor) / mHoldupSolid));// Change of particel temperature
	//// Dosta 2010 eq. 11 with addjustments: no sum necessary as all values are not height dependent, division not applicable as no summation occures
	//if (unit->particlesGlobal)
	//	_res[m_iTempParticle] = _ders[m_iTempParticle] - derTempParticleFormula2;
	//doubleTransfer = derTempParticleFormula2;
	//const double sumEnergyStreamsFilm = Q_PF2 + Q_AF2 + (H_susp2 - HFlowVapor2) * varPhi2; // Sum of energy streams to film
	// Dosta 2010 eq. 9 - no sum necessary as all values are not height dependent, division not applicable as no summation occures
	// _DEBUG
	//if (/*_time > 7929|| _time > 6000 || _time > 5948 || _time > 300 ||*/ _time > 2797.52)
	//	bool breakint = true;
	// Codes for checking _vars, _ders and _res
	//std::vector vars2(_vars, _vars + GetVariablesNumber());
	//std::vector ders2(_ders, _ders + GetVariablesNumber());
	//std::vector  res2(_res, _res + GetVariablesNumber());
}

void CUnitDAEModel::ResultsHandler(double _time, double* _vars, double* _ders, void* _unit)
{
/// Define unit, streams and holdup ///
	auto* unit = static_cast<CDryerBatch*>(_unit); 
	const CStream* inGasStream = unit->GetPortStream("InletFluidizationGas");
	const CStream* inLiquidStream = unit->GetPortStream("InletLiquid");
	const CStream* inNozzleAirStream = unit->GetPortStream("InletNozzleAir");
	CStream* outGasStream = unit->GetPortStream("OutletExhaustGas");
	CHoldup* holdupSolid = unit->GetHoldup("HoldupSolid");
	CHoldup* holdupLiquid = unit->GetHoldup("HoldupLiquid");
	CHoldup* holdupGas = unit->GetHoldup("HoldupGas");

/// Read parameters ///
	/// Phase properties as constant
	const heatCapacity C_PGas = unit->C_PGas;
	const thermalConductivity lambdaGas = unit->lambdaGas;
	const heatCapacity C_PWaterLiquid = unit->C_PWaterLiquid;
	const heatCapacity C_PWaterVapor = unit->C_PWaterVapor;
	const specificLatentHeat Delta_h0 = unit->Delta_h0;
	const heatCapacity C_PParticle = unit->C_PParticle;
	const density rhoGas = unit->rhoGas;
	const density rhoVapor = unit->rhoVapor;
	const density rhoLiquid = unit->rhoWater;
	const density rhoParticle = unit->rhoParticle;
	/// Particle in holdup
	const mass mHoldupSolid = holdupSolid->GetPhaseMass(_time, EPhase::SOLID); // Solid mass holdup [kg] - constant for fluidization with water
	const mass mHoldupLiquid = holdupLiquid->GetPhaseMass(_time, EPhase::LIQUID);
	const length Delta_f = unit->GetConstRealParameterValue("Delta_f") * 1e-6; // convert in [m]
	const area A_P = unit->CalculateParticleSurfaceArea(_time);
	const length d32 = unit->CalculateHoldupSauter(_time);
	/// Hydrodynamics
	const double eps_0 = unit->GetConstRealParameterValue("eps_0");
	const double u_Gas = unit->CalculateGasVel(_time, d32);
	/// Liquid in holdup
	const temperature T_holdupLiquid = holdupLiquid->GetTemperature(_time); // in [K]
	const temperature theta_holdupLiquid = T_holdupLiquid - unit->T_ref; // in degreeC
	/// Inlet fluidization gas
	const massFraction y_inGas = inGasStream->GetCompoundFraction(_time, unit->keyVapor);
	const moistureContent Y_inGas = unit->ConvertMassFracToMoistContent(y_inGas);
	const massFlow mFlowInGas = inGasStream->GetMassFlow(_time); // Gas mass flow [kg/s]
	const massFlow mFlowInGasDry = inGasStream->GetCompoundMassFlow(_time, unit->keyGas);
	const temperature T_inGas = inGasStream->GetTemperature(_time); // temperature in [K]
	const temperature theta_inGas = T_inGas - unit->T_ref; // temperature in [degreeC]
	const specificLatentHeat h_inGas = C_PGas * theta_inGas + Y_inGas * (C_PWaterVapor * theta_inGas + Delta_h0);
	const specificLatentHeat h_inGasDry = C_PGas * theta_inGas;
	const specificLatentHeat h_inGasVapor = Y_inGas * (C_PWaterVapor * theta_inGas + Delta_h0);
	/// Inlet nozzle gas
	const massFlow mFlowInNozzleGas = inNozzleAirStream->GetMassFlow(_time);
	const massFraction y_nozzle = inNozzleAirStream->GetCompoundFraction(_time, unit->keyVapor);
	const moistureContent Y_nozzle = unit->ConvertMassFracToMoistContent(y_nozzle);
	const massFlow mFlowInNozzleGasDry = inNozzleAirStream->GetCompoundMassFlow(_time, unit->keyGas);
	const temperature T_nozzleGas = inNozzleAirStream->GetTemperature(_time); // in [K]
	const temperature thetaNozzleGas = T_nozzleGas - unit->T_ref;
	const specificLatentHeat h_nozzleGas = C_PGas * thetaNozzleGas + Y_nozzle * (C_PWaterVapor * thetaNozzleGas + Delta_h0);
	const specificLatentHeat h_nozzleGasDry = C_PGas * thetaNozzleGas;
	const specificLatentHeat h_nozzleGasVapor = Y_nozzle * (C_PWaterVapor * thetaNozzleGas + Delta_h0);
	/// Spray liquid
	const massFlow mFlowSprayLiquid = inLiquidStream->GetPhaseMassFlow(_time, EPhase::LIQUID);
	const massFraction x_wSusp = inLiquidStream->GetPhaseFraction(_time, EPhase::LIQUID);
	const temperature T_sprayLiquid = inLiquidStream->GetTemperature(_time);
	const temperature thetaSprayLiquid = inLiquidStream->GetTemperature(_time) - unit->T_ref;
	const specificLatentHeat h_susp = thetaSprayLiquid * (C_PParticle /*coating material*/ * (1 - x_wSusp) + C_PWaterLiquid * x_wSusp);
	/// gas in holdup
	const mass mHoldupGas = holdupGas->GetMass(_time);
	const mass mHoldupGasDry = holdupGas->GetCompoundMass(0, unit->keyGas); // set intial dry air mass as dry air mass during whole process
	const temperature varT_gasHoldup = holdupGas->GetTemperature(_time);
	const temperature varTheta_gasHoldup = varT_gasHoldup - unit->T_ref;
	const pressure pressureGasHoldup = holdupGas->GetPressure(_time); // Pressure holdup [Pa]
	/// gas outlet
	moistureContent varYOutGas = _vars[m_iYOutGas];
	temperature varTempOutGas = _vars[m_iTempOutGas];
	temperature varThetaOutGas = varTempOutGas - unit->T_ref;
	const double varHFlowOutGasFormula = (mFlowInGasDry + mFlowInNozzleGasDry) * (C_PGas * varThetaOutGas + varYOutGas * (C_PWaterVapor * varThetaOutGas + Delta_h0));
	const double varP_sat_Formula = unit->CalculateGasSaturationPressure(varTheta_gasHoldup, pressureGasHoldup);
	/// particle 
	const temperature varTempParticle = _vars[m_iTempParticle];
	const temperature varThetaParticle = varTempParticle - unit->T_ref;
	const double varPhi = _vars[m_iPhi];
	const moistureContent varX = varPhi * A_P * Delta_f * unit->rhoWater / mHoldupSolid; 
	/// liquid film
	const temperature varTempFlim = _vars[m_iTempFilm];
	const temperature varThetaFilm = varTempFlim - unit->T_ref;
	const moistureContent varY_sat_Formula = unit->CalculateGasSaturationMoistureContent(varThetaFilm, pressureGasHoldup);
	const double varRH_eq_Formula = unit->CalculateGasEquilibriumRelativeHumidity(varX);
	const moistureContent varY_eq_Formula = unit->CalculateGasEquilibriumMoistureContent(pressureGasHoldup, varP_sat_Formula, varRH_eq_Formula);
	/// Heat transfer
	const double varAlpha_GP_Formula = unit->CalculateAlpha_GP(_time, varTheta_gasHoldup, d32);
	const double varAlpha_GF_Formula = unit->CalculateAlpha_GF(_time, varTheta_gasHoldup, d32);
	const double varAlpha_PF_Formula = unit->CalculateAlpha_PF(/*varTempFlim, pressureGasHoldup, d32*/ varAlpha_GP_Formula);
	/// Mass transfer
	const double varD_a_Formula = unit->CalculateDiffusionCoefficient(varTheta_gasHoldup);
	const double shrinkFactorBeta = unit->GetConstRealParameterValue("Shrink factor beta");
	const double varBeta_FG_Formula = _time < 10 ? shrinkFactorBeta * unit->CalculateBeta(_time, d32, varThetaOutGas, varD_a_Formula) : unit->CalculateBeta(_time, d32, varThetaOutGas, varD_a_Formula);
	/// water vapor
	const bool usef = unit->GetCheckboxParameterValue("Use relative drying rate?");
	const double f = usef ? unit->CalculateRelativeDryingRate(varX) : 1;
	const size_t Y_maxMethod = unit->GetComboParameterValue("Use Y_eq or Y_sat?");
	double Y_max = 0;
	if (Y_maxMethod == 0) // use Y_eq
	{
		Y_max = varY_eq_Formula;
	}
	else // use Y_sat
	{
		Y_max = varY_sat_Formula;
	}
	double varMFlowVaporFormula = varBeta_FG_Formula * A_P * varPhi * rhoGas * (Y_max - varYOutGas) * f;
	const double varHFlowVaporFormula = varMFlowVaporFormula * (C_PWaterVapor * varTheta_gasHoldup + Delta_h0);
	// heat flow
	const double tolTemp = unit->GetConstRealParameterValue("Tolerance temperature");
	const double shrinkFactorGP = unit->GetConstRealParameterValue("Shrink factor alpha GP");
	const double shrinkFactorGF = unit->GetConstRealParameterValue("Shrink factor alpha GF");
	double varAlpha_GP = varAlpha_GP_Formula;
	double varAlpha_GF = varAlpha_GF_Formula;
	double varAlpha_PF = varAlpha_PF_Formula;
	double varAlpha_GP_modify = shrinkFactorGP * varAlpha_GP_Formula;
	double varAlpha_GF_modify = shrinkFactorGF * varAlpha_GF_Formula;
	double varAlpha_PF_modify = varAlpha_PF_Formula;
	double diffTempGP = varTheta_gasHoldup - varThetaParticle;
	double diffTempGF = varTheta_gasHoldup - varThetaFilm;
	double diffTempPF = varThetaParticle - varThetaFilm;
	double varQFlow_GP_Formula = 0;
	double varQFlow_GF_Formula = 0;
	double varQFlow_PF_Formula = 0;
	/// TODO: make the tolerance range differenciable ///
	/*
	modify diffTemp: so that its change is not sharp (jumping between tol and 0)
	using the tanh function
	idea: diffTemp = tol * fun + zero * (1 - fun) -> see Vasyl's functions
	*/


	if (varPhi < 1)
	{
		if (abs(diffTempGP) < tolTemp)
		{
			varAlpha_GP = varAlpha_GP_modify;
		}
		varQFlow_GP_Formula = varAlpha_GP * A_P * (1 - varPhi) * diffTempGP;
	}
	if (abs(diffTempGF) < tolTemp)
	{
		varAlpha_GF = varAlpha_GF_modify;
	}
	varQFlow_GF_Formula = varAlpha_GF * A_P * varPhi * diffTempGF;
	if (abs(diffTempPF) < tolTemp)
	{
		varAlpha_PF = varAlpha_PF_modify;
	}
	varQFlow_PF_Formula = varAlpha_PF * A_P * varPhi * diffTempPF;
	temperature T_surface = unit->IterateSurfaceTemp(_time, varT_gasHoldup, d32);
	const double QFlow_GW_Chamber = unit->CalculateHeatLossWall(_time, unit->wallThickness, unit->GetConstRealParameterValue("H_plant"), unit->GetConstRealParameterValue("d_bed"), varTheta_gasHoldup, T_surface, unit->lambdaWall, d32);

/// Set holdup properties ///
	// time points
	holdupSolid->AddTimePoint(_time);
	holdupLiquid->AddTimePoint(_time);
	holdupGas->AddTimePoint(_time);
	// set compound masses
	holdupSolid->SetCompoundMass(_time, unit->keySolid, EPhase::SOLID, mHoldupSolid);
	holdupLiquid->SetCompoundMass(_time, unit->keyLiquid, EPhase::LIQUID, mHoldupSolid * varX);
	holdupGas->SetCompoundMass(_time, unit->keyGas, EPhase::VAPOR, mHoldupGasDry);
	holdupGas->SetCompoundMass(_time, unit->keyVapor, EPhase::VAPOR, mHoldupGasDry * varYOutGas); // set holdup gas moisture as outlet
	// set temperature
	holdupSolid->SetTemperature(_time, varTempParticle);
	holdupLiquid->SetTemperature(_time, varTempFlim);
	holdupGas->SetTemperature(_time, varTempOutGas);

/// Set outlet properties ///
	outGasStream->CopyFromStream(_time, inGasStream); // initialization of outlet stream is important!!!
	outGasStream->SetPhaseFraction(_time, EPhase::SOLID, 0);
	outGasStream->SetPhaseFraction(_time, EPhase::LIQUID, 0);
	outGasStream->SetPhaseFraction(_time, EPhase::VAPOR, 1);
	outGasStream->SetTemperature(_time, varTempOutGas);
	const massFlow mFlowOutgas = (mFlowInGasDry + mFlowInNozzleGasDry) * (1 + varYOutGas);
	outGasStream->SetCompoundMass(_time, unit->keyVapor, EPhase::VAPOR, mFlowOutgas* unit->ConvertMoistContentToMassFrac(varYOutGas));
	outGasStream->SetCompoundMass(_time, unit->keyGas, EPhase::VAPOR, mFlowOutgas* (1 - unit->ConvertMoistContentToMassFrac(varYOutGas)));

/// Set state variables ///
	unit->SetStateVariable("Holdup air mass [kg]", mHoldupGasDry, _time);
	unit->SetStateVariable("Holdup water vapor mass [kg]", mHoldupGas - mHoldupGasDry, _time);
	unit->SetStateVariable("Outlet gas Y [g/kg]", varYOutGas * 1e3, _time);
	unit->SetStateVariable("Outlet gas RH [%]", 100 * (unit->CalculateGasRelativeHumidity(varYOutGas, varThetaOutGas, pressureGasHoldup)), _time);
	unit->SetStateVariable("Particle moisture content [%]", varX * 100, _time);
	unit->SetStateVariable("Particle water mass fraction [%]", unit->ConvertMoistContentToMassFrac(varX) * 100, _time);
	unit->SetStateVariable("Particle wetness degree [%]", varPhi * 100, _time);
	unit->SetStateVariable("Water evaporation rate [g/s]", varMFlowVaporFormula * 1e3, _time);
	unit->SetStateVariable("Relative drying rate [-]", f, _time);

/// Plotting ///
	// temperatures
	unit->AddPointOnCurve("TEMPERATURES", "Holdup solid [degC]", _time, varThetaParticle);
	unit->AddPointOnCurve("TEMPERATURES", "Holdup liquid film [degC]", _time, varThetaFilm);
	unit->AddPointOnCurve("TEMPERATURES", "Holdup gas [degC]", _time, varTheta_gasHoldup);
	unit->AddPointOnCurve("TEMPERATURES", "Outlet gas [degC]", _time, varThetaOutGas);
	// water vapor mass or mass flow
	unit->AddPointOnCurve("MASS water vapor", "Inlet vapor mass flow [kg/s]", _time, inGasStream->GetCompoundMassFlow(_time, unit->keyVapor));
	unit->AddPointOnCurve("MASS water vapor", "Outlet vapor mass flow [kg/s]", _time, outGasStream->GetCompoundMassFlow(_time, unit->keyVapor));
	unit->AddPointOnCurve("MASS water vapor", "Vapor mass in holdup [kg]", _time, holdupGas->GetCompoundMass(_time, unit->keyVapor));
	// water liquid mass or mass flow 
	unit->AddPointOnCurve("MASS water liquid", "Inlet liquid mass flow [kg/s]", _time, inLiquidStream->GetCompoundMassFlow(_time, unit->keyLiquid));
	unit->AddPointOnCurve("MASS water liquid", "Outlet liquid mass flow [kg/s]", _time, outGasStream->GetCompoundMassFlow(_time, unit->keyLiquid));
	unit->AddPointOnCurve("MASS water liquid", "Liquid mass in holdup [kg]", _time, holdupLiquid->GetCompoundMass(_time, unit->keyLiquid));
	// energy streams
	unit->AddPointOnCurve("HEAT TRANSFER", "Heat transfer GP [J/s]", _time, varQFlow_GP_Formula);
	unit->AddPointOnCurve("HEAT TRANSFER", "alpha_GP [W/(m2*K)]", _time, varAlpha_GP);
	unit->AddPointOnCurve("HEAT TRANSFER", "Heat transfer GF [J/s]", _time, varQFlow_GF_Formula);
	unit->AddPointOnCurve("HEAT TRANSFER", "alpha_GF [W/(m2*K)]", _time, varAlpha_GF);
	unit->AddPointOnCurve("HEAT TRANSFER", "Heat transfer PF [J/s]", _time, varQFlow_PF_Formula);
	unit->AddPointOnCurve("HEAT TRANSFER", "alpha_PF [W/(m2*K)]", _time, varAlpha_PF);
	unit->AddPointOnCurve("HEAT TRANSFER", "Water vapor flow enthalpy [J/s]", _time, varHFlowVaporFormula);
	unit->AddPointOnCurve("HEAT TRANSFER", "beta [m/s]", _time, varBeta_FG_Formula);
	unit->AddPointOnCurve("HEAT TRANSFER", "Heat loss gas to wall [J/s]", _time, QFlow_GW_Chamber);
	// total energy flow / holdup
	unit->AddPointOnCurve("ENERGY total", "Inlet enthalpy flow [J/s]", _time, mFlowInNozzleGasDry * h_nozzleGas + mFlowInGasDry * h_inGas + mFlowSprayLiquid * h_susp);
	unit->AddPointOnCurve("ENERGY total", "Outlet enthalpy flow [J/s]", _time, varHFlowOutGasFormula);
	unit->AddPointOnCurve("ENERGY total", "Enthalpy in holdup [J]", _time, mHoldupSolid* C_PParticle* varThetaParticle +
																		   mHoldupLiquid * C_PWaterLiquid * varThetaFilm +
																	       mHoldupGasDry * (C_PGas * varTheta_gasHoldup + varYOutGas * (C_PWaterVapor * varTheta_gasHoldup + Delta_h0)));
	// water vapor energy flow / holdup
	unit->AddPointOnCurve("ENERGY water vapor", "Inlet enthalpy flow [J/s]", _time, inGasStream->GetCompoundMassFlow(_time, unit->keyVapor) * (C_PWaterVapor* theta_inGas + Delta_h0) + inNozzleAirStream->GetCompoundMassFlow(_time, unit->keyVapor) * (C_PWaterVapor * thetaNozzleGas + Delta_h0));
	unit->AddPointOnCurve("ENERGY water vapor", "Outlet enthalpy flow [J/s]", _time, outGasStream->GetCompoundMassFlow(_time, unit->keyVapor) * (C_PWaterVapor * varThetaOutGas + Delta_h0));
	unit->AddPointOnCurve("ENERGY water vapor", "Enthalpy in holdup [J]", _time, holdupGas->GetCompoundMass(_time, unit->keyVapor) * (C_PWaterVapor * varTheta_gasHoldup + Delta_h0));
	// water liquid energy flow / holdup
	unit->AddPointOnCurve("ENERGY water liquid", "Inlet enthalpy flow [J/s]", _time, inLiquidStream->GetCompoundMassFlow(_time, unit->keyLiquid) * C_PWaterLiquid * thetaSprayLiquid);
	unit->AddPointOnCurve("ENERGY water liquid", "Outlet enthalpy flow [J/s]", _time, outGasStream->GetCompoundMassFlow(_time, unit->keyLiquid)* C_PWaterLiquid * varThetaOutGas);
	unit->AddPointOnCurve("ENERGY water liquid", "Enthalpy in holdup [J]", _time, holdupLiquid->GetCompoundMass(_time, unit->keyLiquid) * C_PWaterLiquid * varThetaFilm);
	// dry air energy flow / holdup
	unit->AddPointOnCurve("ENERGY dry air", "Inlet enthalpy flow [J/s]", _time, inGasStream->GetCompoundMassFlow(_time, unit->keyGas) * C_PGas * theta_inGas + inNozzleAirStream->GetCompoundMassFlow(_time, unit->keyGas) * C_PGas * thetaNozzleGas);
	unit->AddPointOnCurve("ENERGY dry air", "Outlet enthalpy flow [J/s]", _time, outGasStream->GetCompoundMassFlow(_time, unit->keyGas) * C_PGas * varThetaOutGas);
	unit->AddPointOnCurve("ENERGY dry air", "Enthalpy in holdup [J]", _time, holdupGas->GetCompoundMass(_time, unit->keyGas) * C_PGas * varTheta_gasHoldup);
	// solid Cellets energy flow / holdup
	unit->AddPointOnCurve("ENERGY particle", "Inlet enthalpy flow [J/s]", _time, 0);
	unit->AddPointOnCurve("ENERGY particle", "Outlet enthalpy flow [J/s]", _time, 0);
	unit->AddPointOnCurve("ENERGY particle", "Enthalpy in holdup [J]", _time, holdupSolid->GetCompoundMass(_time, unit->keySolid) * C_PParticle * varThetaParticle);

/// Print results for each simulation time point ///
	const bool printResult = unit->GetCheckboxParameterValue("Print intermediate results");
	if (printResult)
	{
		/// Print simulation results for each time step ///
		unit->ShowInfo("Time = " + std::to_string(_time) + " s:");
		unit->ShowInfo("\tGas temperature = " + std::to_string(varThetaOutGas) + " degreeC");
		unit->ShowInfo("\tParticle temperature = " + std::to_string(varThetaParticle) + " degreeC");
		unit->ShowInfo("\tWater film temperature = " + std::to_string(varThetaFilm) + " degreeC");
		unit->ShowInfo("\tFilm:");
		unit->ShowInfo("\t\tQFlow_GF in = " + std::to_string(varQFlow_GF_Formula) + " J/s");
		unit->ShowInfo("\t\tQFlow_PF in = " + std::to_string(varQFlow_PF_Formula) + " J/s");
		unit->ShowInfo("\t\tHFlowSusp in = " + std::to_string(mFlowSprayLiquid * h_susp) + " J/s");
		unit->ShowInfo("\t\tHFlowVapor out = " + std::to_string(varMFlowVaporFormula * (unit->C_PWaterVapor * varThetaOutGas + unit->Delta_h0)) + " J/s");
		unit->ShowInfo("\tParticle:");
		unit->ShowInfo("\t\tQFlow_GP in = " + std::to_string(varQFlow_GP_Formula) + " J/s");
		unit->ShowInfo("\t\tQFlow_PF out = " + std::to_string(varQFlow_PF_Formula) + " J/s");
		unit->ShowInfo("\tGas:");
		unit->ShowInfo("\t\tHFlowVapor in = " + std::to_string(varMFlowVaporFormula * (unit->C_PWaterVapor * varThetaOutGas + unit->Delta_h0)) + " J/s");
		unit->ShowInfo("\t\tHFlowIn in = " + std::to_string(mFlowInGasDry * h_inGas) + " J/s");
		unit->ShowInfo("\t\tHFlowNozzle in = " + std::to_string(mFlowInNozzleGasDry * h_nozzleGas) + " J/s");
		unit->ShowInfo("\t\tHFlowOut out = " + std::to_string((mFlowInGasDry + mFlowInNozzleGasDry) * (unit->C_PGas * varThetaOutGas + varYOutGas * (unit->C_PWaterVapor * varThetaOutGas + unit->Delta_h0))) + " J/s");
		unit->ShowInfo("\t\tQFlow_GP out = " + std::to_string(varQFlow_GP_Formula) + " J/s");
		unit->ShowInfo("\t\tQFlow_GF out = " + std::to_string(varQFlow_GF_Formula) + " J/s");
		unit->ShowInfo("\t\tQFlow_GW out = " + std::to_string(QFlow_GW_Chamber) + " J/s");
		//unit->ShowInfo("\tsum in&out = " + std::to_string(varMFlowVaporFormula * (unit->C_PWaterVapor * varThetaOutGas + unit->Delta_h0) + mFlowInGasDry * h_inGas + mFlowInNozzleGasDry * h_nozzleGas - (mFlowInGasDry + mFlowInNozzleGasDry) * (unit->C_PGas * varThetaOutGas + _vars[m_iYOutGas] * (unit->C_PWaterVapor * varThetaOutGas + unit->Delta_h0)) - varQFlow_GP_Formula - varQFlow_GF_Formula) + " J/s");
		unit->ShowInfo("\tHeat & mass transfer:");
		unit->ShowInfo("\t\talpha_GF = " + std::to_string(varAlpha_GF) + " W/(m2*K)");
		unit->ShowInfo("\t\talpha_GP = " + std::to_string(varAlpha_GP) + " W/(m2*K)");
		unit->ShowInfo("\t\talpha_PF = " + std::to_string(varAlpha_PF) + " W/(m2*K)");
		unit->ShowInfo("\t\tbeta = " + std::to_string(varBeta_FG_Formula) + " m/s");
	}

	//massFlow mFlowInGas = unit->m_inGasStream->GetMassFlow(_time);
	//const massFlow mFlowInLiquid = unit->m_inLiquidStream->GetPhaseMassFlow(_time, EPhase::LIQUID);
	//const mass mHoldup = unit->m_holdup->GetMass(0);
	//const mass mHoldupSolid = unit->m_holdup->GetPhaseMass(_time, EPhase::SOLID);
	////const mass mHoldupGas = unit->m_holdup->GetPhaseMass(_time, EPhase::GAS);
	//const pressure mPressure = unit->m_holdup->GetPressure(_time);
	//mass mHoldupLiquid = unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID);
	//const area A_P = unit->A_P;
	//const length Delta_f = unit->Delta_f;
	//const moistureContent Y_in = unit->Y_inGas;
	//const massFlow mFlowVapor = VaporFlowStorage;
	//const temperature mTempParticle = _vars[m_iTempParticle];
	//const moistureContent mY_gOut = _vars[m_iYOutGas + (unit->N_total - 1)];
	//const double mTempGasOut = _vars[m_iTempOutGas + (unit->N_total - 1)];
	//const temperature mTempFilm = _vars[m_iTempFilm];
	
	//const mass mWaterNozzle = (_time - prevTime) * mFlowInLiquid;
	//const mass mWaterAirIn = (_time - prevTime) * mFlowInGas * Y_in;
	//const mass mVapor = (_time - prevTime) * mFlowVapor;

	//unit->m_holdup->AddStream(prevTime, _time, unit->m_inLiquidStream);
	////mHoldupLiquid = unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID);
	//unit->m_holdup->RemoveTimePointsAfter(_time);

	//unit->m_holdup->SetPhaseMass(_time, EPhase::LIQUID, unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID) - mVapor);
	//mHoldupLiquid = unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID);
	//unit->m_holdup->SetTemperature(_time, mTempParticle);

	//
	//const massFlow mFLowGasOut = unit->m_outExhaustGasStream->GetMassFlow(_time);

	//// Height averaged temperature of the gas phase [K]
	//// Dosta 2010 eq. 14a
	////const double avVarTempGas = CalculateAverage(_vars, m_iTempOutGas, unit->N_total);
	////unit->SetStateVariable("Avg. temperature gas", avVarTempGas - unit->T_ref, _time);
	//// Output of average gas temperature in celsius
	////
	//// Height averaged moisture content of the gas phase [kg/kg]
	//// Dosta 2010 eq. 14b
	////const double avY = CalculateAverage(_vars, m_iYOutGas, unit->N_total);
	////unit->SetStateVariable("Avg. Y", avY, _time);

	////const double avRH = unit->GetRelativeHumidity(avY, avVarTempGas, mPressure);
	////const double Xeq = !unit->particlesGlobal ? 0 : unit->CalcuateSolidEquilibriumMoistureContent(_time, mTempParticle, avRH);
	////const double XCutt = unit->REAinv(Xeq, 1. - unit->phiCuttOff);
	////double mPhi = unit->particlesGlobal ? (mHoldupLiquid - XCutt * mHoldupSolid) / (unit->rhoWater * Delta_f * unit->A_P) : 0;
	////
	////if (mPhi < 0)
	////	mPhi = 0;
	////if (mPhi > 1)
	////	unit->RaiseWarning("Attention: Overwetting in the dryer at " + std::to_string(_time) + " s");
	//// Give warning if degree of wetness (phi) exceeds 1
	//// Warning texts in GUI

	//*if (unit->particlesGlobal) {
	//	unit->SetStateVariable("X", X, _time);
	//	unit->SetStateVariable("w%", X / (1.0 + X) * 100, _time);
	//	unit->SetStateVariable("phi", mPhi, _time);
	//	unit->SetStateVariable("solid mass in holdup", mHoldupSolid, _time);
	//	unit->SetStateVariable("liquid mass in holdup", unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID), _time);
	//	unit->SetStateVariable("evaporation rate [kg/s]", mFlowVapor, _time);
	//	unit->SetStateVariable("Temperature film", mTempFilm - unit->T_ref, _time);
	//	unit->SetStateVariable("Temperature particle", mTempParticle - unit->T_ref, _time);
	//	unit->SetStateVariable("Moisture precentage", unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID) * 100 / (unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID) + unit->m_holdup->GetPhaseMass(_time, EPhase::SOLID)), _time);
	//	unit->SetStateVariable("RH particles", unit->GetRelativeHumidity(avY, avVarTempGas, mPressure) * 100, _time);
	//}
	//unit->SetStateVariable("T_out", mTempGasOut - unit->T_ref, _time);
	//unit->SetStateVariable("Y_out", mY_gOut, _time);
	//unit->SetStateVariable("RH_out [%]", unit->GetRelativeHumidity(mY_gOut, mTempGasOut, mPressure) * 100, _time);
	//unit->SetStateVariable("gas mass in holdup", mHoldupGas, _time);
	//size_t probeLayer = unit->DetermineLayersInSectionFilledWithBed(0, unit->heightOfChamberTemperatureProbe / unit->chamber.at(0).height) - 1;
	//unit->SetStateVariable("Temperature chamber", _vars[m_iTempOutGas + probeLayer] - unit->T_ref, _time);
	//ParallelFor(unit->N_total, [&](size_t i)
	//	{
	//		unit->SetStateVariable("Temperature gas height layer " + std::to_string(i), _vars[m_iTempOutGas + i] - unit->T_ref, _time);
	//		unit->SetStateVariable("Y height layer " + std::to_string(i), _vars[m_iYOutGas + i], _time);
	//	});*/

	//if (debugToggle) {
	//	unit->SetStateVariable("Mass balance", (_time - prevTime == 0 ? 0 : (-unit->m_holdup->GetMass(_time) + unit->m_holdup->GetMass(prevTime)) / (_time - prevTime)) + ((unit->m_inGasStream->GetMassFlow(_time) + unit->m_inNozzleAirStream->GetMassFlow(_time)) * (1 + Y_in) - unit->m_outExhaustGasStream->GetMassFlow(_time) + unit->m_inLiquidStream->GetMassFlow(_time)), _time);
	//	//unit->SetStateVariable("Energy balance", unit->m_inGasStream->GetMassFlow(_time) - unit->m_outExhaustGasStream->GetMassFlow(_time),_time);
	//	//unit->SetStateVariable("Mass balance new", -unit->m_holdup->GetMass(_time) + unit->m_holdup->GetMass(prevTime) + (_time - prevTime) * (unit->m_inGasStream->GetMassFlow(_time)*(1+Y_in) - unit->m_outExhaustGasStream->GetMassFlow(_time) + unit->m_inSuspensionStream->GetMassFlow(_time)), _time);

	//	double EnergyLiquidPhase = unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID) * unit->C_PWaterLiquid * mTempFilm;
	//	double EnergySolidPhase = unit->m_holdup->GetPhaseMass(_time, EPhase::SOLID) * unit->C_PParticle * mTempParticle;
	//	double EnergyGasPhase = unit->m_holdup->GetPhaseMass(_time, EPhase::GAS) * ((unit->C_PGas + unit->C_PWaterVapor * (unit->YavgOld - avY)) * avVarTempGas/* + unit->Delta_h0 * (unit->YavgOld - Y_av)*/);

	//	double EnergyInputGas = (_time - prevTime) * unit->m_inGasStream->GetMassFlow(_time) * ((unit->C_PGas + unit->Y_inGas * unit->C_PWaterVapor) * unit->m_inGasStream->GetTemperature(_time)/* + unit->Y_in * unit->Delta_h0*/);
	//	double EnergyInputLiquid = (_time - prevTime) * (unit->m_inLiquidStream->GetPhaseMassFlow(_time, EPhase::LIQUID) * unit->C_PWaterLiquid + unit->m_inLiquidStream->GetPhaseMassFlow(_time, EPhase::SOLID) * unit->C_PParticle) * unit->m_inLiquidStream->GetTemperature(_time);
	//	double EnergyOutputGas = (_time - prevTime) * (unit->m_outExhaustGasStream->GetMassFlow(_time) * ((unit->C_PGas + mY_gOut * unit->C_PWaterVapor) * mTempGasOut/* + mY_out * unit->Delta_h0*/));

	//	double HeatLoss = heatLossTransfer;
	//	/*
	//	const double EnergyChangeLiquidPhase = (-unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID) * mTempFilm + unit->m_holdup->GetPhaseMass(prevTime, EPhase::LIQUID) * unit->TempLiquidOld)* unit->C_PWaterLiquid;
	//	const double EnergyChangeSolidPhase = unit->m_holdup->GetPhaseMass(_time, EPhase::SOLID) * unit->C_PParticle * (-mTempParticle + unit->TempSolidOld);
	//	const double EnergyChangeGasPhase = unit->m_holdup->GetPhaseMass(_time, EPhase::GAS) * (unit->C_PGas + unit->C_PWaterVapor * (unit->YavgOld-Y_av)) * (-varAvTempGas + unit->TempGasOld);
	//	const double EnergyInputGas = (_time - prevTime) * unit->m_inGasStream->GetMassFlow(_time) * ((unit->C_PGas + unit->Y_in * unit->C_PWaterVapor) * unit->m_inGasStream->GetTemperature(_time));
	//	const double EnergyInputLiquid = (_time - prevTime) * unit->m_inSuspensionStream->GetMassFlow(_time) * (unit->C_PWaterLiquid * unit->w_l + unit->C_PParticle * (1 - unit->w_l))* unit->m_inSuspensionStream->GetTemperature(_time);
	//	const double EnergyOutputGas = (_time - prevTime) * (-unit->m_outExhaustGasStream->GetMassFlow(_time) * (
	//		(unit->C_PGas + mY_out * unit->C_PWaterVapor) * T_out
	//		+ (mY_out - unit->Y_in) * unit->Delta_h0));
	//	const double EnergyBalance = EnergyChangeLiquidPhase + EnergyChangeSolidPhase + EnergyChangeGasPhase + EnergyInputGas + EnergyInputLiquid + EnergyOutputGas;
	//	unit->TempLiquidOld = mTempFilm;
	//	unit->TempSolidOld = mTempParticle;
	//	unit->TempGasOld = varAvTempGas;
	//	unit->YavgOld = Y_av;
	//	*/
	//	double DeltaEnergyLiquid = EnergyLiquidPhase - unit->EnergyLiquidPhaseOld;
	//	double DeltaEnergySolid = EnergySolidPhase - unit->EnergySolidPhaseOld;
	//	double DeltaEnergyGas = EnergyGasPhase - unit->EnergyGasPhaseOld;
	//	double DeltaGasStreams = EnergyInputGas - EnergyOutputGas;
	//	double DeltaEnergy = DeltaEnergyLiquid + DeltaEnergySolid + DeltaEnergyGas;
	//	double EnergyBalance = -DeltaEnergy + DeltaGasStreams + EnergyInputLiquid - unit->HeatLossOld;
	//	/*unit->EnergyLiquidPhaseOld = EnergyLiquidPhase;
	//	unit->EnergyGasPhaseOld = EnergyGasPhase;
	//	unit->EnergySolidPhaseOld = EnergySolidPhase;
	//	unit->HeatLossOld = HeatLoss;
	//	unit->SetStateVariable("Energy balance new", EnergyBalance, _time);*/
	//	EnergyGasPhase = unit->m_holdup->GetPhaseMass(_time, EPhase::GAS) * ((unit->C_PGas + unit->C_PWaterVapor * avY) * avVarTempGas + unit->Delta_h0 * avY);
	//	EnergyInputGas = unit->m_inGasStream->GetMassFlow(_time) * ((unit->C_PGas + unit->Y_inGas * unit->C_PWaterVapor) * unit->m_inGasStream->GetTemperature(_time) + unit->Y_inGas * unit->Delta_h0);
	//	EnergyInputLiquid = (unit->m_inLiquidStream->GetPhaseMassFlow(_time, EPhase::LIQUID) * unit->C_PWaterLiquid + unit->m_inLiquidStream->GetPhaseMassFlow(_time, EPhase::SOLID) * unit->C_PParticle) * unit->m_inLiquidStream->GetTemperature(_time);
	//	EnergyOutputGas = (unit->m_outExhaustGasStream->GetMassFlow(_time) * ((unit->C_PGas + mY_gOut * unit->C_PWaterVapor) * mTempGasOut + mY_gOut * unit->Delta_h0));
	//	EnergyBalance = (_time - prevTime == 0 ? 0 : -(EnergyGasPhase + EnergyLiquidPhase + EnergySolidPhase) + (unit->EnergyGasPhaseOld + unit->EnergyLiquidPhaseOld + unit->EnergySolidPhaseOld) / (_time - prevTime)) + EnergyInputGas + EnergyInputLiquid - EnergyOutputGas - unit->HeatLossOld;
	//	unit->SetStateVariable("Energy balance", EnergyBalance, _time);
	//	unit->EnergyLiquidPhaseOld = EnergyLiquidPhase;
	//	unit->EnergyGasPhaseOld = EnergyGasPhase;
	//	unit->EnergySolidPhaseOld = EnergySolidPhase;

	//	//unit->SetStateVariable("debug0", EnergyBalance!=0?mPhi/(-EnergyBalance):0, _time);// ~mPhi ~X
	//	/*
	//	unit->SetStateVariable("debug1",unit->GetGasSaturationMoistureContent(varAvTempGas, mPressure)/unit->Y_sat,_time);

	//	//unit->SetStateVariable("debug15", (mPressure - CachePressureRange <= inputCache[1][1] && inputCache[1][1] <= mPressure + CachePressureRange)?1:0, _time);
	//	//unit->SetStateVariable("debug17",resultCache[0],_time);
	//	unit->SetStateVariable("debug18", unit->GetRelativeHumidity(Y_av + Y_in, varAvTempGas, mPressure),_time);
	//	double RH = unit->GetRelativeHumidity(Y_av, varAvTempGas, mPressure);
	//	unit->SetStateVariable("debug19", RH, _time);
	//	double Xeq = unit->CalcuateSolidEquilibriumMoistureContent(_time, mTempParticle, RH);
	//	double deltaX = unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID) / mHoldupSolid - Xeq;
	//	unit->SetStateVariable("debug20", deltaX, _time);
	//	double REA = unit->REA(deltaX);
	//	unit->SetStateVariable("debug21", REA, _time);
	//	double curve = 1 - REA;
	//	unit->SetStateVariable("debug22", curve, _time);
	//	unit->SetStateVariable("debug24", Xeq-XeqTransfer, _time);
	//	unit->SetStateVariable("debug25", X-curveTransfer, _time);
	//	if (_time >= 134)
	//	{
	//		std::vector vars(_vars, _vars + GetVariablesNumber());
	//		std::vector ders(_ders, _ders + GetVariablesNumber());
	//		//std::vector  res(_res, _res + GetVariablesNumber());
	//		double deltaTime = _time - unit->m_holdup->GetPreviousTimePoint(_time);
	//		std::vector<double> DebugVector;
	//		DebugVector.reserve(vars.size() + ders.size()+3); // preallocate memory
	//		DebugVector.insert(DebugVector.end(), vars.begin(), vars.end());
	//		DebugVector.insert(DebugVector.end(), ders.begin(), ders.end());
	//		DebugVector.push_back(deltaX);
	//		DebugVector.push_back(MFlowVaporLiquid);
	//		DebugVector.push_back(mFlowVapor);
	//		//DebugCache.push_back(DebugVector);
	//		if (_time >= 136)
	//			double test = 1;
	//	}
	//	*/
	//	/*
	//	std::vector ders(_ders, _ders + GetVariablesNumber());
	//	std::pair<double,size_t> output = std::make_pair(0,-1);
	//	for (size_t i = 0; i < ders.size();i++ )
	//		if (abs(ders[i])>output.first)
	//			output = std::make_pair(ders[i], i);
	//	std::ofstream myfile;
	//	std::string fileLocation = "E:\\Dyssol\\output.csv";
	//	myfile.open(fileLocation, std::fstream::app);
	//	myfile << _time << "," << output.first << "," << output.second << "\n";
	//	myfile.close();*/
	//	unit->SetStateVariable("debug0", progressCounterTotal, _time);
	//	unit->SetStateVariable("debug1", liquidSideLimitedGlobal, _time);
	//	unit->SetStateVariable("debug2", TransferBool, _time);
	//	unit->SetStateVariable("debug3", doubleTransfer, _time);
	//	double X_eq = !unit->particlesGlobal ? 0 : unit->CalcuateSolidEquilibriumMoistureContent(_time, mTempParticle, unit->GetRelativeHumidity(avY, avVarTempGas));
	//	unit->SetStateVariable("debug4", X_eq, _time);
	//	unit->SetStateVariable("debug5", X - X_eq, _time);
	//	unit->SetStateVariable("debug6", HeatLoss, _time);
	//	unit->SetStateVariable("debug7", (mTempFilm - unit->T_ref) - (mTempFilm - unit->T_ref) * TransferBool, _time);
	//	if (_time >= 99)
	//		bool breakPoint = true;
	//	//unit->SetStateVariable("debug3", RHTransfer, _time);
	//	//unit->SetStateVariable("debug3", unit->GetRelativeHumidity(Y_av, varAvTempGas, mPressure), _time);
	//	std::vector vars(_vars, _vars + GetVariablesNumber());
	//	std::vector ders(_ders, _ders + GetVariablesNumber());
	//} // DEBUG
	//	//Temp particle climbs while temps of gas and film fall
	//	/*
	//	std::pair<double, std::vector<double>> result = CalculateTotalHeatLoss(_time, _unit,X,_vars, true);
	//	unit->ShowInfo("Heat loss - particles");
	//	unit->ShowInfo(std::to_string(result.first));
	//	unit->ShowInfo(std::to_string(result.second[0]));
	//	result = CalculateTotalHeatLoss(_time, _unit, X, _vars, false);
	//	unit->ShowInfo("Heat loss - gas");
	//	unit->ShowInfo(std::to_string(result.first));
	//	unit->ShowInfo(std::to_string(result.second[0]));
	//	unit->ShowInfo("");
	//	unit->ShowInfo(std::to_string(XeqTransfer));
	//	unit->ShowInfo("");*/
	//	//unit->ShowInfo(std::to_string(CalculateAlpha_PW(mTempParticle, X, varAvTempGas, mPressure, _unit, _time)* unit->AreaGasParticleWallHeatLoss* (mTempParticle - unit->TemperatureWall)));
	//	//unit->ShowInfo(std::to_string(unit->CalckA(CalculateAlpha_PW(mTempParticle, X, varAvTempGas, mPressure, _unit, _time), 0, 0.35, { 0.225 ,0.230 }, { 15 }) * (mTempParticle - unit->T_inf)));
	//	//unit->ShowInfo(std::to_string(unit->CalcAlphaOutside(_time, h, d.front(), /*unit->TemperatureWall*/0.4 * (std::max(unit->m_inGasStream->GetTemperature(_time), unit->m_holdup->GetTemperature(_time)) + unit->T_inf), unit->T_inf)));
}


///////////////////////////////////////////////////////////////////
/// Read, calculate and check material properties from databank ///
///////////////////////////////////////////////////////////////////

//void CDryerBatch::PullCompoundDataFromDatabase(double _time)
//{
//	// Get general data for compounds
//	this->compoundKeys = this->GetAllCompounds();
//	std::vector<std::pair< EPhase, int>> compoundsKeyIndexPhasePartnerIndex(compoundKeys.size(), std::make_pair(EPhase::UNDEFINED, -2));
//	std::vector<density> CompoundDensities(compoundKeys.size()); // Compound densities
//	std::vector<heatCapacity> CompoundHeatCapacities(compoundKeys.size()); // Compound heat capacities
//	std::vector<thermalConductivity> CompoundThermalConductivities(compoundKeys.size()); // Compound thermal conductivities
//	std::vector<dynamicViscosity> CompoundDynViscosities(compoundKeys.size()); // Compound dynamic viscosities
//	std::vector<molarMass> CompoundMolarMasses(compoundKeys.size()); // Compound molar masses
//	//std::vector<double> CompoundCriticaleqData.temperatures(compoundKeys.size()); // Compound critical eqData.temperatures
//	for (int i = 0; i < compoundKeys.size(); i++)
//	{
//		CompoundDensities[i] = GetCompoundProperty(compoundKeys[i], ECompoundTPProperties::DENSITY);
//		CompoundHeatCapacities[i] = GetCompoundProperty(compoundKeys[i], ECompoundTPProperties::HEAT_CAPACITY_CP);
//		CompoundThermalConductivities[i] = GetCompoundProperty(compoundKeys[i], ECompoundTPProperties::THERMAL_CONDUCTIVITY);
//		CompoundDynViscosities[i] = GetCompoundProperty(compoundKeys[i], ECompoundTPProperties::VISCOSITY);
//		CompoundMolarMasses[i] = GetCompoundProperty(compoundKeys[i], ECompoundConstProperties::MOLAR_MASS);
//		//CompoundCriticaleqData.temperatures[i] = GetCompoundProperty(compoundKeys[i], ECompoundConstProperties::CRITICAL_TEMPERATURE);
//	}
//
//	/// Properties of solids ///
//	if (m_holdup->GetPhaseFraction(_time, EPhase::SOLID) == 0)
//	{
//		RaiseWarning("Dryer holdup contains no solid.");
//	}
//	std::vector<massFraction> SolidCompoundsDistribution = m_holdup->GetPhase(EPhase::SOLID)->GetCompoundsDistribution(_time);
//	density tempRhoParticle = 0;
//	heatCapacity tempHeatCapacityParticle = 0;
//	thermalConductivity tempThermalConductivityParticle = 0;
//	for (int i = 0; i < compoundKeys.size(); i++)
//	{
//		if (SolidCompoundsDistribution[i] > 0)
//		{
//			compoundsKeyIndexPhasePartnerIndex[i] = std::make_pair(EPhase::SOLID, -1);
//			tempRhoParticle += SolidCompoundsDistribution[i] * CompoundDensities[i];
//			tempHeatCapacityParticle += SolidCompoundsDistribution[i] * CompoundHeatCapacities[i];
//			tempThermalConductivityParticle += SolidCompoundsDistribution[i] * CompoundThermalConductivities[i];
//		}
//	}		
//	rhoParticle = tempRhoParticle;
//	C_PParticle = tempHeatCapacityParticle;
//	lambdaParticle = tempThermalConductivityParticle;
//	//beta_GP = this->GetConstRealParameterValue("beta_GP");
//
//	/// Properties of liquids ///
//	if (m_inLiquidStream->GetPhaseFraction(_time, EPhase::LIQUID) == 0)
//	{
//		RaiseWarning("Port: InletSuspension contains no liquid.");
//	}
//	std::vector<massFraction> LiqudCompoundsDistribution = m_inLiquidStream->GetPhase(EPhase::LIQUID)->GetCompoundsDistribution(_time);
//	density tempRhoLiquid = 0;
//	heatCapacity tempHeatCapacityLiquid = 0;
//	thermalConductivity tempLambdaLiquid = 0;
//	for (int i = 0; i < compoundKeys.size(); i++)
//	{
//		tempRhoLiquid = tempRhoLiquid + LiqudCompoundsDistribution[i] * CompoundDensities[i];
//		tempHeatCapacityLiquid = tempHeatCapacityLiquid + LiqudCompoundsDistribution[i] * CompoundHeatCapacities[i];
//		tempLambdaLiquid = tempLambdaLiquid + LiqudCompoundsDistribution[i] * CompoundThermalConductivities[i];
//	}
//	rhoWater = tempRhoLiquid;
//	C_PWaterLiquid = tempHeatCapacityLiquid;
//	lambdaWater = tempLambdaLiquid;
//
//	/// Propterties of vapors (gas phase) ///
//	if (m_inGasStream->GetPhaseFraction(_time, EPhase::GAS) != 1)
//	{
//		RaiseError("Port: InletFluidizationGas contains non gaseous phases.");
//	}
//	std::vector<double> VaporCompoundsDistribution = m_inGasStream->GetPhase(EPhase::GAS)->GetCompoundsDistribution(_time);
//	density tempRhoGas = 0;
//	heatCapacity tempHeatCapacityGas = 0;
//	thermalConductivity templambdaGas = 0;
//	dynamicViscosity tempEtaGas = 0;
//	molarMass tempMolarMassGas = 0;
//	//temperature tempCriticaleqData.temperatures = 0;
//	for (int i = 0; i < compoundKeys.size(); i++)
//		if (VaporCompoundsDistribution[i] > 0 && compoundsKeyIndexPhasePartnerIndex[i].first == EPhase::UNDEFINED)
//		{
//			compoundsKeyIndexPhasePartnerIndex[i] = std::make_pair(EPhase::GAS, -1);
//			tempRhoGas += VaporCompoundsDistribution[i] * CompoundDensities[i];
//			tempHeatCapacityGas += VaporCompoundsDistribution[i] * CompoundHeatCapacities[i];
//			templambdaGas += VaporCompoundsDistribution[i] * CompoundThermalConductivities[i];
//			tempEtaGas += VaporCompoundsDistribution[i] * CompoundDynViscosities[i];
//			tempMolarMassGas += VaporCompoundsDistribution[i] * CompoundMolarMasses[i];
//			//tempCriticaleqData.temperatures += VaporCompoundsDistribution[i] * CompoundCriticaleqData.temperatures[i];
//		}
//	rhoGas = tempRhoGas;
//	C_PGas = tempHeatCapacityGas;
//	lambdaGas = templambdaGas;
//	etaGas = tempEtaGas;
//	molarMassGas = tempMolarMassGas;
//	//T_critGas = tempCriticaleqData.temperatures;
//	for (int i = 0; i < compoundsKeyIndexPhasePartnerIndex.size(); i++)
//	{
//		if (compoundsKeyIndexPhasePartnerIndex[i].first == EPhase::UNDEFINED)
//		{
//			RaiseError("Compound " + GetCompoundName(compoundKeys[i]) + "could not be fitted into solid, liquid or vapor/gas phase");
//		}
//	}
//	CompoundsKeyIndexPhasePartnerIndex = compoundsKeyIndexPhasePartnerIndex;
//
//	/// Vector to store indices of liquid (first) and vapor (second) form of compound ///
//	std::vector<std::pair<int, int>> liquidVaporIndices;
//	for (int i = 0; i < compoundKeys.size(); i++) // Liquid form
//	{
//		if (LiqudCompoundsDistribution[i] > 0 && compoundsKeyIndexPhasePartnerIndex[i].first == EPhase::UNDEFINED)
//		{
//			double heatOfVaporizationLiquid = GetCompoundProperty(compoundKeys[i], ECompoundConstProperties::HEAT_OF_VAPORIZATION_AT_NORMAL_BOILING_POINT);
//			double heatOfVaporizationVapor;
//			for (int j = 0; j < compoundKeys.size(); j++) // Vapor form
//			{
//				if (i != j)
//				{
//					heatOfVaporizationVapor = GetCompoundProperty(compoundKeys[j], ECompoundConstProperties::HEAT_OF_VAPORIZATION_AT_NORMAL_BOILING_POINT);
//					if (heatOfVaporizationLiquid == heatOfVaporizationVapor)
//					{
//						liquidVaporIndices.push_back(std::make_pair(i, j));
//					}
//				}
//			}
//		}
//	}
//	if (liquidVaporIndices.size() == 0)
//		RaiseError("No compounds for phase transition found. Please add vapor or liquid form of compound or check material database entries for heat of vaporiation.");
//	else
//	{
//		for (int i = 0; i < liquidVaporIndices.size(); i++)
//			ShowInfo("Compounds found for phase change: " + GetCompoundName(compoundKeys[liquidVaporIndices[i].first]) + "," + GetCompoundName(compoundKeys[liquidVaporIndices[i].second]));
//		compoundsKeyIndexPhasePartnerIndex[liquidVaporIndices[0].first] = std::make_pair(EPhase::LIQUID, liquidVaporIndices[0].second);
//		compoundsKeyIndexPhasePartnerIndex[liquidVaporIndices[0].second] = std::make_pair(EPhase::VAPOR, liquidVaporIndices[0].first);
//		indicesOfVaporOfPhaseChangingCompound = liquidVaporIndices[0];
//
//		for (int i = 0; i < compoundKeys.size(); i++) // Liquid form
//			if (LiqudCompoundsDistribution[i] > 0 && compoundsKeyIndexPhasePartnerIndex[i].first == EPhase::UNDEFINED)
//				compoundsKeyIndexPhasePartnerIndex[i] = std::make_pair(EPhase::LIQUID, -1);
//	}
//	if (liquidVaporIndices.size() > 1)
//	{
//		RaiseWarning("Only the first pair will be considered for phase change calculations.");
//	}
//	C_PWaterVapor = GetCompoundProperty(compoundKeys[liquidVaporIndices[0].second], ECompoundTPProperties::HEAT_CAPACITY_CP); // Water vapor heat capacity
//	Delta_h0 = GetCompoundProperty(compoundKeys[liquidVaporIndices[0].first], ECompoundConstProperties::HEAT_OF_VAPORIZATION_AT_NORMAL_BOILING_POINT) / GetCompoundProperty(compoundKeys[liquidVaporIndices[0].first], MOLAR_MASS); // Latent heat (evaporation heat) at 0 degree
//	molarMassPhaseChangingLiquid = CompoundMolarMasses[liquidVaporIndices[0].first];
//
//	/// print material properties on simulation window in DEBUG mode ///
//	if (debugToggle){
//		std::ostringstream  os;
//		os << "Debug info\n"
//			<< "\trhoGas: " << rhoGas << "\n"
//			<< "\tetaGas: " << etaGas << "\n"
//			<< "\tC_PGas: " << C_PGas << "\n"
//			<< "\tlambdaGas: " << lambdaGas << "\n"
//			<< "\tmolarMassGas: " << molarMassGas << "\n" << "\n"
//
//			<< "\trhoWater: " << rhoWater << "\n"
//			<< "\tC_PWaterLiquid: " << C_PWaterLiquid << "\n"
//			<< "\tC_PWaterVapor: " << C_PWaterVapor << "\n"
//			<< "\tDelta_h0: " << Delta_h0 << "\n"
//			<< "\tlambdaWater: " << lambdaWater << "\n"
//			<< "\tmolarMassPhaseChangingLiquid: " << molarMassPhaseChangingLiquid << "\n" << "\n"
//
//			<< "\trhoParticle: " << rhoParticle << "\n"
//			<< "\tbeta_GP: ";
//		/*if (beta_GP < 0)
//			os << "automatic calculation";
//		else
//			os << beta_GP;*/
//		os  << "\n"
//			<< "\tC_PParticle: " << C_PParticle << "\n";
//		ShowInfo(os.str());
//		os.str("");
//	} // DEBUG
//}


/////////////////////////////////////////////////////////////////
/// Functions to calculate/return moisture-related properties ///
///						for gas and solid					  ///
/////////////////////////////////////////////////////////////////
massFraction CDryerBatch::ConvertMoistContentToMassFrac(moistureContent Y) const
{
	return Y / (1. + Y);
}

moistureContent CDryerBatch::ConvertMassFracToMoistContent(massFraction y) const
{
	return y / (1. - y);
}

pressure CDryerBatch::CalculateGasSaturationPressure(temperature theta_Gas, pressure pressureGas) const // P_sat
{
	const double A = 8.07131; // in[mmHg]
	const double B = 1730.61;
	const double C = 233.426; // in[degC]
	const double ratioMM = molarMassPhaseChangingLiquid / molarMassGas;
	return pow(10, (A - B / (C + theta_Gas))) * (pressureGas / 760); //Antoine equation (Source: Springer-Verlag Berlin Heidelberg 2014, E. Drioli, L. Giorno (eds.), Encyclopedia of Membranes, DOI 10.1007/978-3-642-40872-4_26-1). convert[mmHg] to[Pa]
}

moistureContent CDryerBatch::CalculateGasSaturationMoistureContent(temperature theta_Gas, pressure pressureGas) const  // Y_sat
{
	const double ratioMM = molarMassPhaseChangingLiquid / molarMassGas;
	pressure P_sat = CalculateGasSaturationPressure(theta_Gas, pressureGas); //Antoine equation (Source: Springer-Verlag Berlin Heidelberg 2014, E. Drioli, L. Giorno (eds.), Encyclopedia of Membranes, DOI 10.1007/978-3-642-40872-4_26-1). convert[mmHg] to[Pa]
	return ratioMM * P_sat / (pressureGas - P_sat); 
}

double CDryerBatch::CalculateGasRelativeHumidity(moistureContent Y, temperature thetaGas, pressure pressure) //const
{
	const moistureContent Y_sat = CalculateGasSaturationMoistureContent(thetaGas, pressure);
	double RH = Y / Y_sat; // relative humidity
	if (RH > 1)
		RH = 1;
	if (RH < 0)
		RH = 0;
	return RH;
}

double CDryerBatch::CalculateGasEquilibriumRelativeHumidity(/*double _time, temperature temperature,*/ moistureContent X) const // RH_eq
{
	double RH_eq = -17.256 * pow(X, 2.0) + 10.207 * X; // from measurement at 25 degreeC, still not consider temperature->TODO
	if (RH_eq >= 1)
	{
		return 1;
	}
	else
	{
		return RH_eq;
	}
}

moistureContent CDryerBatch::CalculateGasEquilibriumMoistureContent(pressure pressureGas, pressure P_sat, double RH_eq) const // Y_eq
{
	const double ratioMM = molarMassPhaseChangingLiquid / molarMassGas;
	return ratioMM * RH_eq * P_sat / (pressureGas - RH_eq * P_sat); // RH_eq * P_sat = P_eq

	//if (RH == 0)
	//{
	//	return 0;
	//}
	//pressure P_sat = GetCompoundProperty(compoundKeys[indicesOfVaporOfPhaseChangingCompound.first], ECompoundTPProperties::VAPOR_PRESSURE, temperatureParticle, pressureGas); // Saturation partial pressure of phase changing compound at particle temperature [Pa] // Equilibrium partial pressure [Pa]
	//const pressure P_eq = P_sat * RH;
	//const moistureContent Y_eq = ratioMM * P_eq / (pressureGas - P_eq); // Moisture content in kg per kg dry gas
	//return Y_eq;
}


///////////////////////////////////
/// Bed and particle properties ///
///////////////////////////////////
length CDryerBatch::CalculateHoldupSauter(double _time) const
{
	std::vector<double> sizeGrid = GetNumericGrid(DISTR_SIZE);
	std::vector<double> q3_holdup = m_holdupSolid->GetPSD(_time, PSD_q3);
	double d32 = GetSauterDiameter(sizeGrid, q3_holdup);
	return d32;
}

area CDryerBatch::CalculateParticleSurfaceArea(double _time) const
{
	std::vector<double> Grid = GetNumericGrid(DISTR_SIZE);
	const double M_tot = m_holdupSolid->GetPhaseMass(_time, EPhase::SOLID);
	std::vector<double> q_3 = m_holdupSolid->GetPSD(_time, PSD_q3, EPSDGridType::DIAMETER);
	const double A = GetSpecificSurface(Grid, q_3) / rhoParticle * M_tot;
	return A;
}

length CDryerBatch::CalculateFluidizedBedHeight(length H_fix, double eps) const // Soeren diss. eq. (4.18)
{
	const double eps_0 = GetConstRealParameterValue("eps_0");
	return H_fix * (1 - eps_0) / (1 - eps);

	//const double massSolid = m_holdup->GetPhaseMass(_time, EPhase::SOLID);
	//const double densitySolid = GetAvgTPCompoundProperty(_time, EPhase::SOLID, ECompoundTPProperties::DENSITY, particleTemperature, m_holdup->GetPressure(_time));
	//const double volumeSolid = massSolid / densitySolid;
	//const double eps = CalculateBedPorosity(_time);
	//const double volumeBed = volumeSolid / (1. - eps);
	//double volumeOfFilledSections = 0;
	//size_t numberOfFilledSection = 0;
	//double heightOfFilledSection = 0;
	//while (volumeBed > (CalculateSectionVolume(numberOfFilledSection) + volumeOfFilledSections))
	//{
	//	volumeOfFilledSections += CalculateSectionVolume(numberOfFilledSection);
	//	heightOfFilledSection += chamber.at(numberOfFilledSection).height;
	//	numberOfFilledSection++;
	//	if (numberOfFilledSection >= chamber.size())
	//		RaiseError("Bed expands out of chamber.");
	//}
	//double r = chamber.at(numberOfFilledSection).dimensionsInternal.at(0).first / 2;
	//double R = chamber.at(numberOfFilledSection).dimensionsInternal.at(0).second / 2;
	//double h = chamber.at(numberOfFilledSection).height;
	//const double tanWallAngle = (R - r) / h;
	//const double heightBed = heightOfFilledSection - std::cbrt(-MATH_PI * pow(r, 3) * pow(tanWallAngle, 3) - 3 * pow(tanWallAngle, 4) * volumeBed) / (std::cbrt(MATH_PI) * pow(tanWallAngle, 2)) - r / tanWallAngle;
	//return heightBed;
}

///double CDryerBatch::CalculateBedHeightOrDetermineSectionsFilledWithBed(double _time, double particleTemperature, bool outputHeight)
//{
//	const double massSolid = m_holdup->GetPhaseMass(_time, EPhase::SOLID);
//	const double densitySolid = GetAvgTPCompoundProperty(_time, EPhase::SOLID, ECompoundTPProperties::DENSITY, particleTemperature, m_holdup->GetPressure(_time));
//	const double volumeSolid = massSolid / densitySolid;
//	const double eps = CalculateBedPorosity(_time);
//	const double volumeBed = volumeSolid / (1. - eps);
//	double volumeOfFilledSections = 0;
//	size_t numberOfFilledSection = 0;
//	double heightOfFilledSection = 0;
//	while (volumeBed > (CalculateSectionVolume(numberOfFilledSection) + volumeOfFilledSections))
//	{
//		volumeOfFilledSections += CalculateSectionVolume(numberOfFilledSection);
//		heightOfFilledSection += chamber.at(numberOfFilledSection).height;
//		numberOfFilledSection++;
//		if (numberOfFilledSection >= chamber.size())
//			RaiseError("Bed expands out of chamber.");
//	}
//
//	double r = chamber.at(numberOfFilledSection).dimensionsInternal.at(0).first / 2;
//	double R = chamber.at(numberOfFilledSection).dimensionsInternal.at(0).second / 2;
//	double h = chamber.at(numberOfFilledSection).height;
//	const double tanWallAngle = (R - r) / h;
//	const double bedHeight = -std::cbrt(-MATH_PI * pow(r, 3) * pow(tanWallAngle, 3) - 3 * pow(tanWallAngle, 4) * volumeBed) / (std::cbrt(MATH_PI) * pow(tanWallAngle, 2)) - r / tanWallAngle;
//	const double totalBedHeight = heightOfFilledSection + bedHeight;
//	if (outputHeight)
//		return totalBedHeight;
//	else
//		return numberOfFilledSection + bedHeight / h;
//}

///void CDryerBatch::SetupChamber()
//{
//	if (chamber.size() == 0)
//	{
//		chamberSection section;
//		section.name = "main";
//		section.shape = EShape::CYLINDRICAL;
//		section.dimensionsInternal = { std::make_pair(diamOfBed, 0.250) };
//		section.height = { heightOfChamber };
//		section.wallThicknesses = { 0.005 };
//		section.thermalConductivities = { 15 };
//		section.layers = 1;
//		chamber.push_back(section);
//
//		section.name = "adapter";
//		section.dimensionsInternal = { std::make_pair(0.250, 0.500) };
//		section.height = { 0.345 };
//		chamber.push_back(section);
//
//		section.name = "expander";
//		section.dimensionsInternal = { std::make_pair(0.500, 0.500) };
//		section.height = { 0.950 }; // 0.500 from Glatt, 0.95 measured
//		chamber.push_back(section);
//	}
//}

/////////////////////////////
/// Dimensionless numbers ///
/////////////////////////////
dimensionlessNumber CDryerBatch::CalculateReynolds(double _time, length d32) const
{
	double u_Gas = CalculateGasVel(_time, d32);
	return (d32 * u_Gas * rhoGas) / etaGas;
}

dimensionlessNumber CDryerBatch::CalculatePrandtl(temperature avgGasTheta) const // source fo Maksym???
{
	return 0.702 - 2.63e-4 * avgGasTheta - 1.05e-6 * pow(avgGasTheta, 2) - 1.52e-9 * pow(avgGasTheta, 3);
	//return C_PGas * etaGas / lambdaGas;
}

dimensionlessNumber CDryerBatch::CalculateSchmidt(double D_a) const
{
	return etaGas / (rhoGas * D_a);
}

dimensionlessNumber CDryerBatch::CalculateArchimedes(length d32) const
{
	return STANDARD_ACCELERATION_OF_GRAVITY * pow(d32, 3) * (rhoParticle - rhoGas) * rhoGas / pow(etaGas, 2);
}

dimensionlessNumber CDryerBatch::CalculateNusseltFree(double _time, temperature T_surface) const
{
	const temperature theta_env = GetConstRealParameterValue("theta_env");
	const temperature T_env = theta_env + T_ref;
	const dimensionlessNumber Pr_env = CalculatePrandtl(theta_env);
	dimensionlessNumber Gr = CalculateGrashof(T_env, T_surface);
	const dimensionlessNumber f_Pr = pow(1 + pow(0.492 / Pr_env, 9 / 16), -16 / 9); // source: VDI-Warmeatlas, Kap. F2
	const dimensionlessNumber Nu_free = pow(0.825 + 0.387 * pow(Pr_env * Gr * f_Pr, 1 / 6), 2); // source: VDI-Warmeatlas, Kap. F2
	return Nu_free;
}

/// CalculateReynolds for discretized height, CURRENTLY NOT IN USE
//dimensionlessNumber CDryerBatch::CalculateReynolds(double _time, size_t section) const
//{
//	double d = (chamber.at(section).dimensionsInternal.at(0).first + chamber.at(section).dimensionsInternal.at(0).second) / 2;
//	return CalculateGasVel(_time,section)* d* rhoGas / etaGas;
//}

/// CheckForSmallBiot, CURRENTLY NOT IN USE
//bool CDryerBatch::CheckForSmallBiot(double _time) const // Needed for uniform temperature distribution (Diss Soeren, chapter 2.2.2, page 27)
//{
//	const double Bi = CalcBiotNumber(_time);
//	if (Bi < SmallBiotNumber)
//		return true;
//	else
//		return false;
//}
//bool CDryerBatch::LoadDimensions(std::string path)
//{
//	return false;
//}

/////////////////////
/// Mass transfer ///
/////////////////////
massTransferCoefficient CDryerBatch::CalculateBeta(double _time, length d32, double avgGasTheta, double D_a) const
{
	const double betaVal = GetConstRealParameterValue("beta_GP");
	if (betaVal < 0)
	{
		const dimensionlessNumber Ar = CalculateArchimedes(d32);
		const dimensionlessNumber Pr = CalculatePrandtl(avgGasTheta);
		const dimensionlessNumber Sc = CalculateSchmidt(D_a);
		const dimensionlessNumber Re_s = 18 * pow(sqrt(1. + sqrt(Ar) / 9) - 1, 2.0); // Martin (VDI-Waermeatlas, chapter M5): homogeneous fluidization: sink velocity = gas velocity
		const size_t methodIdx = GetComboParameterValue("Heat & mass transfer methods");
		switch (methodIdx)
		{
			case 0: // Martin (VDI-Waermeatlas, chapter M5)
			{
				const dimensionlessNumber Nu_lam = CalculateNusseltSherwoodLam(Re_s, Pr);
				const dimensionlessNumber Nu_turb = CalculateNusseltSherwoodTurb(Re_s, Pr);
				const dimensionlessNumber Nu = CalculateNusseltSherwood(Nu_lam, Nu_turb);
				const dimensionlessNumber Sh = Nu * pow(Sc / Pr, 1. / 3.); // Lewis number = Sc / Pr
				return Sh * D_a / d32;
			}			
			case 1: // Groenewold & Tsostas (see Rieck dissertation (2020), page 150-151)
			{
				const dimensionlessNumber Re_mf = CalculateReynoldsMF(_time, d32);
				const dimensionlessNumber Re = CalculateReynolds(_time, d32); // based on superficial gas velocity
				const dimensionlessNumber eps_mf = CalculateBedPorosityMF(wadellFactor); // from Wen & Yu. eps_mf != eps_0
				const double eps = CalculateBedPorosity(_time, d32);
				const dimensionlessNumber Re_eps = Re_mf / eps_mf;
				const dimensionlessNumber Nu_lam = CalculateNusseltSherwoodLam(Re_s, Pr); // Soeren diss.
				const dimensionlessNumber Nu_turb = CalculateNusseltSherwoodTurb(Re_eps, Pr);
				const dimensionlessNumber Nu = CalculateNusseltSherwood(Nu_lam, Nu_turb);
				const dimensionlessNumber Nu_app = CalculateNusseltSherwoodApp(Nu, eps_mf);
				const dimensionlessNumber Sh_app = Nu_app * pow(Sc / Pr, 1. / 3.); // Lewis number = Sc / Pr
				const area A_P = CalculateParticleSurfaceArea(_time);
				const length d_bed = GetConstRealParameterValue("d_bed");
				const length H_fix = GetConstRealParameterValue("H_bedFix");
				const length H_fb = CalculateFluidizedBedHeight(H_fix, eps);
				const dimensionlessNumber AvH = CalculateAtoF(H_fb, d32, eps);
				dimensionlessNumber Sh_modify = CalculateNusseltSherwoodModify(Re, Sc, Sh_app, AvH);
				return Sh_modify * D_a / d32;
			}
		}
	}
	else
	{
		return betaVal;
	}
}

/////////////////////
/// Heat transfer ///
/////////////////////
double CDryerBatch::CalculateAlpha_GP(double _time, temperature avgGasTheta, length d32) const
{
	const dimensionlessNumber Ar = CalculateArchimedes(d32);
	const dimensionlessNumber Pr = CalculatePrandtl(avgGasTheta);  
	const dimensionlessNumber Re_s = 18 * pow(sqrt(1. + sqrt(Ar) / 9) - 1, 2.0); // Martin (VDI-Waermeatlas, chapter M5): homogeneous fluidization: sink velocity = gas velocity
	const size_t methodIdx = GetComboParameterValue("Heat & mass transfer methods");
	switch (methodIdx)
	{
		case 0: // Martin (VDI-Waermeatlas, chapter M5)
		{
			const dimensionlessNumber Nu_lam = CalculateNusseltSherwoodLam(Re_s, Pr);
			const dimensionlessNumber Nu_turb = CalculateNusseltSherwoodTurb(Re_s, Pr);
			const dimensionlessNumber Nu = CalculateNusseltSherwood(Nu_lam, Nu_turb);
			return Nu * lambdaGas / d32;
		}
		case 1: // Groenewold & Tsostas (see Rieck dissertation (2020), page 150-151)
		{
			const dimensionlessNumber Re_mf = CalculateReynoldsMF(_time, d32);
			const dimensionlessNumber Re = CalculateReynolds(_time, d32); // based on superficial gas velocity
			const dimensionlessNumber eps_mf = CalculateBedPorosityMF(wadellFactor); // from Wen & Yu. eps_mf != eps_0
			const double eps = CalculateBedPorosity(_time, d32);
			const dimensionlessNumber Re_eps = Re_mf / eps_mf;
			const dimensionlessNumber Nu_lam = CalculateNusseltSherwoodLam(Re_s, Pr); // Soeren diss.
			const dimensionlessNumber Nu_turb = CalculateNusseltSherwoodTurb(Re_eps, Pr);
			const dimensionlessNumber Nu = CalculateNusseltSherwood(Nu_lam, Nu_turb);
			const dimensionlessNumber Nu_app = CalculateNusseltSherwoodApp(Nu, eps_mf);
			const area A_P = CalculateParticleSurfaceArea(_time);
			const length d_bed = GetConstRealParameterValue("d_bed");
			const length H_fix = GetConstRealParameterValue("H_bedFix");
			const length H_fb = CalculateFluidizedBedHeight(H_fix, eps);
			const dimensionlessNumber AvH = CalculateAtoF(H_fb, d32, eps);
			dimensionlessNumber Nu_modify = CalculateNusseltSherwoodModify(Re, Pr, Nu_app, AvH);
			return Nu_modify * lambdaGas / d32;
		}		
		case 2: // self-defined
		{
			return GetConstRealParameterValue("alpha_GP");
		}
	}
}

double CDryerBatch::CalculateAlpha_GF(double _time, temperature avgGasTheta, length d32) const
{
	const dimensionlessNumber Ar = CalculateArchimedes(d32);
	const dimensionlessNumber Pr = CalculatePrandtl(avgGasTheta);
	const dimensionlessNumber Re_s = 18 * pow(sqrt(1. + sqrt(Ar) / 9) - 1, 2.0); // Martin (VDI-Waermeatlas, chapter M5): homogeneous fluidization: sink velocity = gas velocity
	const size_t methodIdx = GetComboParameterValue("Heat & mass transfer methods");
	switch (methodIdx)
	{
		case 0: // Martin (VDI-Waermeatlas, chapter M5)
		{
			const dimensionlessNumber Nu_lam = CalculateNusseltSherwoodLam(Re_s, Pr);
			const dimensionlessNumber Nu_turb = CalculateNusseltSherwoodTurb(Re_s, Pr);
			const dimensionlessNumber Nu = CalculateNusseltSherwood(Nu_lam, Nu_turb);
			return Nu * lambdaGas / d32;
		}
		case 1: // Groenewold & Tsostas (see Rieck dissertation (2020), page 150-151)
		{
			const dimensionlessNumber Re_mf = CalculateReynoldsMF(_time, d32);
			const dimensionlessNumber Re = CalculateReynolds(_time, d32); // based on superficial gas velocity
			const dimensionlessNumber eps_mf = CalculateBedPorosityMF(wadellFactor); // from Wen & Yu. eps_mf != eps_0
			const double eps = CalculateBedPorosity(_time, d32);
			const dimensionlessNumber Re_eps = Re_mf / eps_mf;
			const dimensionlessNumber Nu_lam = CalculateNusseltSherwoodLam(Re_s, Pr); // Soeren diss.
			const dimensionlessNumber Nu_turb = CalculateNusseltSherwoodTurb(Re_eps, Pr);
			const dimensionlessNumber Nu = CalculateNusseltSherwood(Nu_lam, Nu_turb);
			const dimensionlessNumber Nu_app = CalculateNusseltSherwoodApp(Nu, eps_mf);
			const area A_P = CalculateParticleSurfaceArea(_time);
			const length d_bed = GetConstRealParameterValue("d_bed");
			const length H_fix = GetConstRealParameterValue("H_bedFix");
			const length H_fb = CalculateFluidizedBedHeight(H_fix, eps);
			const dimensionlessNumber AvH = CalculateAtoF(H_fb, d32, eps);
			dimensionlessNumber Nu_modify = CalculateNusseltSherwoodModify(Re, Pr, Nu_app, AvH);
			return Nu_modify * lambdaGas / d32;
		}
		case 2: // self-defined
		{
			return GetConstRealParameterValue("alpha_GF");
		}
	}
}

double CDryerBatch::CalculateAlpha_PF(/*temperature tempWater, pressure pressureHoldup, length d32*/ double alpha_GP) const
{
	const double f_a = GetConstRealParameterValue("f_a");
	const size_t methodIdx = GetComboParameterValue("Heat & mass transfer methods");
	switch (methodIdx)
	{
	case 0: // Martin
		return alpha_GP * f_a;
	case 1: // Groenewold & Tsostas
		return alpha_GP * f_a;
	case 2: // self-defined
		return GetConstRealParameterValue("alpha_PF");
	}
		
	// CALCULATE FROM DISS RIECK BASED ON Nu = 2
	//const double lambdaParticle = 0.04;
	//const double Nu = 2.0;
	//const double alpha_PF = Nu * lambdaParticle / 300e-6/*d32*/;
	//return alpha_PF;
}

temperature CDryerBatch::IterateSurfaceTemp(double _time, temperature T_gasHoldup, length d32) const
{
	const double tol = 0.01;
	const temperature theta_env = GetConstRealParameterValue("theta_env");
	const temperature T_env = theta_env + T_ref;
	temperature T_surfaceInit = (T_gasHoldup + T_env) / 2; // initial guess of surface temperature
	temperature T_surfaceOld = T_surfaceInit;
	temperature T_surfaceNew = GetNewTempSurface(_time, T_surfaceOld, (T_gasHoldup - T_ref), d32);
	int current_try = 0;
	int max_try = 100;
	while (abs(T_surfaceNew - T_surfaceOld) > tol && current_try < max_try) 
	{
		T_surfaceOld = T_surfaceNew;
		T_surfaceNew = GetNewTempSurface(_time,T_surfaceOld, T_gasHoldup - T_ref, d32);
		current_try++;
		if (current_try == max_try)
		{
			std::cout << "test4\n";
		}
	}
	return T_surfaceNew;
}

temperature CDryerBatch::GetNewTempSurface(double _time, temperature TempSurfaceOld, temperature thetaInside, length d32) const
{
	temperature thetaSurfaceOld = TempSurfaceOld - T_ref;
	double alpha_out = CalculateNusseltFree(_time, TempSurfaceOld) * lambdaGas / lengthChamber;
	double QFlow_Wall = CalculateHeatLossWall(_time, wallThickness, lengthChamber, radiusChamber, thetaInside, thetaSurfaceOld, lambdaWall, d32);
	area A_out = 2 * MATH_PI * lengthChamber * (radiusChamber + wallThickness);
	return (T_env + QFlow_Wall / (alpha_out * A_out));
}

double CDryerBatch::CalculateHeatLossWall(double _time, length wallThickness, length height, length radiusInner, temperature thetaInside, temperature thetaSurface, double lambdaWall, length d32) const
{
	// heat loss through wall = (thetaGasHoldup - thetaEnv) / (1/(kA)) | 1/(kA): total heat resistance
	// 1/(kA) = 1/(alpha_Holdup * A_in) + wallThickness/(lambdaWall * A_m) + 1/(alphaEnv * A_out) = resistIn + resistWall + resistOut	
	area A_in = 2 * MATH_PI * height * radiusInner;
	double alpha_in = CalculateAlpha_GF(_time, thetaInside, d32); // alpha_Holdup == alphaGP == alphaGF
	double resistIn = 1 / (alpha_in * A_in);
	area A_m = 2 * MATH_PI * height * wallThickness / log((radiusInner + wallThickness) / radiusInner);
	double resistWall = wallThickness / (lambdaWall * A_m);
	area A_out = 2 * MATH_PI * height * (radiusInner + wallThickness);
	dimensionlessNumber alpha_out = CalculateNusseltFree(_time, thetaSurface + T_ref) * lambdaGas / lengthChamber;
	double resistOut = 1 / (alpha_out * A_out);
	double resistTotal = resistIn + resistWall + resistOut;
	const temperature theta_env = GetConstRealParameterValue("theta_env");
	return (thetaInside - theta_env) / resistTotal;
}

////////////////////
/// Fluiddynamic ///
////////////////////
dimensionlessNumber CDryerBatch::CalculateReynoldsMF(double _time, length d32) const
{
	dimensionlessNumber Ar = CalculateArchimedes(d32);
	return 33.7 * (sqrt(1. + 3.6e-5 * Ar) - 1); // correlation Wen&Yu
}

double CDryerBatch::CalculateMinFluidizeVel(double _time, length d32) const
{
	double u_mf = this->GetConstRealParameterValue("u_mf");
	if (u_mf == 0)
	{
		dimensionlessNumber Re_mf = CalculateReynoldsMF(_time, d32);
		u_mf = Re_mf * etaGas / (d32 * rhoGas);
	}
	return u_mf;
}

double CDryerBatch::CalculateGasVel(double _time, length d32) const // superficial gas velocity
{
	const length d_bed = this->GetConstRealParameterValue("d_bed");
	const massFlow mFlow_gasIn = m_inGasStream->GetMassFlow(_time);
	double VFlow_gasIn = mFlow_gasIn / rhoGas;
	const area area_bed = MATH_PI * pow(d_bed, 2.) / 4.;
	return VFlow_gasIn / area_bed;
}

double CDryerBatch::CalculateBedPorosity(double _time, length d32, bool homogeneousFluidization) const 
{
	const double eps_0 = GetConstRealParameterValue("eps_0");
	const double eps_mf = CalculateBedPorosityMF(wadellFactor);
	const size_t methodIdx = GetComboParameterValue("Bed porosity calculation");
	dimensionlessNumber Ar = 0;
	dimensionlessNumber Re_mf = 0;
	dimensionlessNumber Re_elu = 0;
	dimensionlessNumber Re = 0;
	dimensionlessNumber n = 0;
	switch (methodIdx) 
	{
		case 0: // Martin (VDI-Waermeatlas, Kap. M5)
		{
			Ar = CalculateArchimedes(d32);
			Re_mf = 42.9 * (1. - eps_0) * (sqrt(1. + pow(eps_0, 3) * Ar / (3214 * pow(1. - eps_0, 2))) - 1);
			Re_elu = homogeneousFluidization ? 18 * pow(sqrt(1. + sqrt(Ar) / 9) - 1, 2) : sqrt(4 * Ar / 3);
			Re = CalculateReynolds(_time, d32);
			n = log(Re_mf / Re_elu) / log(eps_0);
			return pow(Re / Re_elu, 1. / n);
		}
		case 1: // Soeren diss
		{
			double u_gas = CalculateGasVel(_time, d32);
			double u_mf = CalculateMinFluidizeVel(_time, d32);
			double u_suspGas = (u_gas - u_mf) / 3 + u_mf;
			return eps_mf * pow(u_suspGas / u_mf, 1. / 4.65);
		}
	}
}

//double CDryerBatch::DetermineSectionsFilledWithBed(double _time, double particleTemperature)
//{
//	if (particlesGlobal == false)
//		return 0;
//	double bedHeight = CalculateBedHeight(_time, particleTemperature);
//	int sectionFilled = 0;
//	double chamberHeight = 0;
//	while (chamberHeight < bedHeight)
//	{
//		chamberHeight += chamber.at(sectionFilled).height;
//		sectionFilled++;
//		if (sectionFilled > chamber.size())
//			return -1; // Bed larger than chamber
//	}
//	double filledHeight = chamberHeight - chamber.at(sectionFilled).height;
//	return sectionFilled -1 + (bedHeight - filledHeight) / chamber.at(sectionFilled).height;
//}


//----------------------------------:)---- take a break ----:)-----------------------------------------//

/////////////////////////////////////////////////////////////////
/// function related to drying kinetics, CURRENTLY NOT IN USE ///
/////////////////////////////////////////////////////////////////

double CDryerBatch::CalculateRelativeDryingRate(moistureContent X) const
{
	const double k_dc = GetConstRealParameterValue("k_dc");
	const double X_cr = GetConstRealParameterValue("X_cr");
	const double X_eq = GetConstRealParameterValue("X_eq");
	const double REA_A = GetConstRealParameterValue("A");
	const double REA_B = GetConstRealParameterValue("B");
	const double REA_C = GetConstRealParameterValue("C");
	const size_t methodIdx = GetComboParameterValue("Methods");
	switch (methodIdx)
	{
	case 0: // REA
	{
		return 1 - REA_A * exp(REA_B * pow((X - X_eq), REA_C));
	}
	case 1: // NCDC
	{
		const double normX = (X - X_eq) / (X_cr - X_eq);
		if (X <= X_eq)
		{
			return 0;
		}
		else if (X >= X_cr)
		{
			return 1;
		}
		else
		{
			return k_dc * normX / (1. + normX * (k_dc - 1.)); // for spray drying
		}
	}
	}
}

//double CDryerBatch::GetParticleEquilibriumMoistureContent(double temperature, double RH) const
//{
//	if (eqData.RHs.empty() || eqData.temperatures.empty() || eqData.equilibriumMoistureContents.empty() || temperature != temperature || RH != RH)
//		return -1; // Error
//	std::set<double>::iterator itRHmax = eqData.RHs.end();
//	std::advance(itRHmax, -1);
//	// RH bounds [0,1] check
//	if (RH > 1)
//		RH = 1;
//	if (RH < 0)
//		RH = 0;
//
//	// Temperature bounds [0,] check
//	if (temperature < 0)
//		temperature = 0;
//
//	std::set<double>::iterator itT, itRH, itlowT, itupT, itlowRH, itupRH;
//
//	itlowT = eqData.temperatures.lower_bound(temperature);
//	if (itlowT == eqData.temperatures.end() || *itlowT != temperature)
//		std::advance(itlowT, -1);
//	itupT = eqData.temperatures.upper_bound(temperature);
//	if (itlowT == eqData.temperatures.end())
//	{
//		itlowT = itupT;
//		std::advance(itupT, 1);
//	}
//	if (itupT == eqData.temperatures.end())
//	{
//		itupT = itlowT;
//		std::advance(itlowT, -1);
//	}
//
//	itlowRH = eqData.RHs.lower_bound(RH);
//	if (itlowRH == eqData.RHs.end() || *itlowRH != RH)
//		std::advance(itlowRH, -1);
//	itupRH = eqData.RHs.upper_bound(RH);
//	if (itlowRH == eqData.RHs.end())
//	{
//		itlowRH = itupRH;
//		std::advance(itupRH, 1);
//	}
//	if (itupRH == eqData.RHs.end())
//	{
//		itupRH = itlowRH;
//		std::advance(itlowRH, -1);
//	}
//
//	double Xeq;
//	if (itlowT == eqData.temperatures.end() || itupT == eqData.temperatures.end() || itlowRH == eqData.RHs.end() || itupRH == eqData.RHs.end())
//	{
//
//		if (itlowT == eqData.temperatures.end() || itupT == eqData.temperatures.end())
//		{
//			itT = itlowT != eqData.temperatures.end() ? itlowT : itupT;
//			if (itlowRH == eqData.RHs.end() || itupRH == eqData.RHs.end())
//			{
//				itRH = itlowRH != eqData.RHs.end() ? itlowRH : itupRH;
//				Xeq = eqData.equilibriumMoistureContents.at(std::make_pair(*itT, *itRH));
//				return Xeq;
//			}
//			else
//			{
//				double XeqTRHl = eqData.equilibriumMoistureContents.at(std::make_pair(*itT, *itlowRH));
//				double XeqTRHu = eqData.equilibriumMoistureContents.at(std::make_pair(*itT, *itupRH));
//				Xeq = lerp(XeqTRHl, XeqTRHu, (RH - *itlowRH) / (*itupRH - *itlowRH));
//				return Xeq;
//			}
//
//		}
//		else
//		{
//			itRH = itlowRH != eqData.RHs.end() ? itlowRH : itupRH;
//			double XeqTlRH = eqData.equilibriumMoistureContents.at(std::make_pair(*itlowT, *itRH));
//			double XeqTuRH = eqData.equilibriumMoistureContents.at(std::make_pair(*itupT, *itRH));
//			Xeq = lerp(XeqTlRH, XeqTuRH, (temperature - *itlowT) / (*itupT - *itlowT));
//			return Xeq;
//		}
//		return -1;
//	}
//	else
//	{
//		double XeqTlRHl, XeqTuRHl, XeqTlRHu, XeqTuRHu, XeqTl, XeqTu;
//
//		XeqTlRHl = eqData.equilibriumMoistureContents.at(std::make_pair(*itlowT, *itlowRH));
//		XeqTuRHl = eqData.equilibriumMoistureContents.at(std::make_pair(*itupT, *itlowRH));
//
//		/*if (eqData.equilibriumMoistureContents.contains(std::make_pair(*itlowT, *itupRH)))
//			XeqTlRHu = eqData.equilibriumMoistureContents.at(std::make_pair(*itlowT, *itupRH));
//		else
//			XeqTlRHu = eqData.equilibriumMoistureContents.at(std::make_pair(*itlowT, *eqData.maxRH.at(*itlowT)));*/
//
//		try
//		{
//			XeqTlRHu = eqData.equilibriumMoistureContents.at(std::make_pair(*itlowT, *itupRH));
//		}
//		catch (std::out_of_range const& e)
//		{
//			XeqTlRHu = eqData.equilibriumMoistureContents.at(std::make_pair(*itlowT, *eqData.maxRH.at(*itlowT)));
//		}
//
//		/*if (eqData.equilibriumMoistureContents.contains(std::make_pair(*itupT, *itupRH)))
//			XeqTuRHu = eqData.equilibriumMoistureContents.at(std::make_pair(*itupT, *itupRH));
//		else
//			XeqTuRHu = eqData.equilibriumMoistureContents.at(std::make_pair(*itupT, *eqData.maxRH.at(*itupT)));*/
//
//		try
//		{
//			XeqTuRHu = eqData.equilibriumMoistureContents.at(std::make_pair(*itupT, *itupRH));
//		}
//		catch (std::out_of_range const& e)
//		{
//			XeqTuRHu = eqData.equilibriumMoistureContents.at(std::make_pair(*itupT, *eqData.maxRH.at(*itupT)));
//		}
//
//		XeqTl = lerp(XeqTlRHl, XeqTlRHu, (RH - *itlowRH) / (*itupRH - *itlowRH));
//		XeqTu = lerp(XeqTuRHl, XeqTuRHu, (RH - *itlowRH) / (*itupRH - *itlowRH));
//		Xeq = lerp(XeqTl, XeqTu, (temperature - *itlowT) / (*itupT - *itlowT));
//		return Xeq;
//	}
//	return -1;
//}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Function related to particle X from materials database (sorption isotherm), CURRENTLY UNDER CONSTRUCTION ///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//double CDryerBatch::GetEquilibriumRelativeHumidity(double temperature, double X) const
//{
//	if (eqData.RHs.empty() || eqData.temperatures.empty() || eqData.equilibriumMoistureContents.empty() || temperature != temperature || X != X)
//		return -1; // Error
//	//X = X / moistureScaler;
//
//	if (X < 0)
//		return 0;
//
//	std::set<double>::iterator itlowT, itupT, itT, itRH;
//	double XlowT, XupT, Xpolated, XT;
//
//	itRH = eqData.RHs.end();
//	std::advance(itRH, -1);
//
//	itlowT = eqData.temperatures.lower_bound(temperature);
//	if (itlowT == eqData.temperatures.end() || *itlowT != temperature)
//		std::advance(itlowT, -1);
//	itupT = eqData.temperatures.upper_bound(temperature);
//	if (itlowT == eqData.temperatures.end())
//	{
//		itlowT = itupT;
//		std::advance(itupT, 1);
//	}
//	if (itupT == eqData.temperatures.end())
//	{
//		itupT = itlowT;
//		std::advance(itlowT, -1);
//	}
//
//	if (itlowT == eqData.temperatures.end() || itupT == eqData.temperatures.end())
//	{
//		itT = itlowT != eqData.temperatures.end() ? itlowT : itupT;
//		//if (*itRH > *eqData.maxRH.at(*itT))
//		//	itRH = eqData.maxRH.at(*itT);
//		XT = eqData.equilibriumMoistureContents.at(std::make_pair(*itT, *itRH));
//		if (X > XT)
//			return 1;
//		while (X < XT)
//		{
//			std::advance(itRH, -1);
//			if (itRH == eqData.RHs.end())
//				return -1;
//			XT = eqData.equilibriumMoistureContents.at(std::make_pair(*itT, *itRH));
//			if (X > XT)
//				return *itRH;
//		}
//	}
//	else
//	{
//		XlowT = eqData.equilibriumMoistureContents.at(std::make_pair(*itlowT, *itRH));
//		XupT = eqData.equilibriumMoistureContents.at(std::make_pair(*itupT, *itRH));
//		Xpolated = unit->lerp(XlowT, XupT, (temperature - *itlowT) / (*itupT - *itlowT));
//		if (X > Xpolated)
//			return 1;
//		while (X < XlowT)
//		{
//			std::advance(itRH, -1);
//			if (itRH == eqData.RHs.end())
//				return -1;
//			XlowT = eqData.equilibriumMoistureContents.at(std::make_pair(*itlowT, *itRH));
//			XupT = eqData.equilibriumMoistureContents.at(std::make_pair(*itupT, *itRH));
//			Xpolated = lerp(XlowT, XupT, (temperature - *itlowT) / (*itupT - *itlowT));
//			if (X > Xpolated)
//				return *itRH;
//		}
//	}
//	return -1;
//}

//bool CDryerBatch::InitializeMoistureContentDatabase(std::string path)
//{
//	//CCorrelation::GetParameters();
//	std::ifstream  file(path);
//	if (!file.is_open())
//		RaiseError("Moisture content data file has incorrect path.");
//	std::string line;
//	std::vector<std::vector<std::string>> parsedCsv;
//	while (std::getline(file, line))
//	{
//		std::stringstream lineStream(line);
//		std::string cell;
//		std::vector<std::string> parsedRow;
//		while (std::getline(lineStream, cell, ','))
//			parsedRow.push_back(cell);
//
//		parsedCsv.push_back(parsedRow);
//	}
//
//	for (auto& element : parsedCsv.at(0))
//		if (element != "")
//			eqData.RHs.insert(std::stod(element));
//
//	for (auto& element : parsedCsv)
//		if (element.at(0) != "")
//			eqData.temperatures.insert(std::stod(element.at(0)));
//
//	for (int i = 1; i < parsedCsv.size(); i++)
//	{
//		for (int j = 1; j < parsedCsv.at(i).size(); j++)
//			eqData.equilibriumMoistureContents.insert(std::make_pair(std::make_pair(std::stod(parsedCsv.at(i).at(0)), std::stod(parsedCsv.at(0).at(j))), std::stod(parsedCsv.at(i).at(j))));
//		std::set<double>::iterator itRHmax = eqData.RHs.end();
//		int fromMaxRH = eqData.RHs.size() - parsedCsv.at(i).size() + 2;
//		std::advance(itRHmax, -fromMaxRH);
//		eqData.maxRH.insert(std::make_pair(std::stod(parsedCsv.at(i).at(0)), itRHmax));
//	}
//	return true;
//}

//moistureContent CDryerBatch::CalcuateSolidEquilibriumMoistureContent(double _time, temperature temperature, double RH)
//{
//	std::vector<moistureContent> EquilibriumMoistures(compoundKeys.size());
//	moistureContent particleEquilibriumMoistureContent = 0;
//	std::vector<massFraction> SolidCompoundsDistribution = m_holdup->GetPhase(EPhase::SOLID)->GetCompoundsDistribution(_time);
//	for (int i = 0; i < compoundKeys.size(); i++)
//	{
//		if (compoundKeys.at(i) == eqData.compoundKey)
//			EquilibriumMoistures[i] = GetParticleEquilibriumMoistureContent(temperature, RH);
//		else
//			EquilibriumMoistures[i] = GetCompoundProperty(compoundKeys[i], ECompoundTPProperties::EQUILIBRIUM_MOISTURE_CONTENT, temperature, RH);
//	}		
//	// Equilibrium moisture content depends on temperatur and humidity
//	for (int i = 0; i < compoundKeys.size(); i++)
//	{
//		if (EquilibriumMoistures.at(i) < 0)
//		{
//			std::stringstream os;
//			os << "Moisture content of " << GetCompoundName(compoundKeys.at(i)) << "at " << _time << "s returned a negativ moisture content for " << temperature << "K at " << RH * 100 << "%.";
//			RaiseError(os.str());
//			os.str() = "";
//		}
//	}
//		
//	for (int i = 0; i < compoundKeys.size(); i++)
//	{
//		particleEquilibriumMoistureContent += SolidCompoundsDistribution[i] * EquilibriumMoistures[i];
//	}
//	particleEquilibriumMoistureContent *= moistureScaler;
//	return std::max(particleEquilibriumMoistureContent, minMoistureContent);
//}

//////////////////////////////////////////////////////////////
/// Heat loss to environment, CURRENTLY UNDER CONSTRUCTION ///
//////////////////////////////////////////////////////////////

//double CDryerBatch::CalculateAlpha_PW(double _t_p, double _t_g, double _p, double _time) const
//{
//	/// Get unit parameters
//	double eps_mf = eps_0;				// Bed porosity at minimum fluidization, equal to INITIAL bed porosity [-]
//	double eps_fb = CalculateBedPorosity(_time);						// Bed porosity at operating conditions [-]
//	double dp_p = CalculateHoldupSauter(_time);							// Sauter diameter of particles [m]
//
//	// Calculate modified free path of the gas molecule
//	double C_A = 2.8;
//	double gamma = 1. / (pow(10., 0.6 - (1000. / _t_g + 1.) / C_A) + 1); // gamma, A.3.85
//
//	double l = 2. * (2. / gamma - 1.) * sqrt(2. * MATH_PI * MOLAR_GAS_CONSTANT * _t_g / molarMassGas) * lambdaGas / _p / (2. * C_PGas - MOLAR_GAS_CONSTANT / molarMassGas); //modified free path of the gas molecule, A.3.84
//
//	//Calculate parameter Z and N
//	double C_k = 2.6; //C_k, A.3.82
//	double Nu_pw_max = 4. * ((1. + 2. * l / dp_p) * log(1. + dp_p / 2. / l) - 1.); //Nu_PW,max, A.3.83
//	double Z = 1. / 6. * rhoParticle * C_PParticle / lambdaGas * sqrt(STANDARD_ACCELERATION_OF_GRAVITY * pow(dp_p, 3.) * (eps_fb - eps_mf) / 5. / (1 - eps_mf) / (1 - eps_fb)); //Parameter Z, A.3.80
//	double N = Nu_pw_max / C_k / Z; //Parameter N, A.3.81
//
//	//Calculate Nusselt Number and heat transfer coefficient
//	double Nu_pw = (1. - eps_fb) * Z * (1. - exp(-N));
//	double alpha_pw = Nu_pw * lambdaGas / dp_p;
//
//	return alpha_pw;
//}

//double CDryerBatch::CalculateAlpha_GW(double _time, temperature avgGasTemperature) const
//{
//	const double Pr = CalculatePrandtl(_time, avgGasTemperature);
//	const double Ar = CalculateArchimedes(_time, d32);
//	const double NuG = 0.009 * pow(Pr, 1. / 3) * sqrt(Ar);
//	double d = CalculateHoldupSauter(_time);
//	const double alphaGW = NuG * lambdaGas / d;
//	return alphaGW;
//}

//double CDryerBatch::CalckAc(double alphaIn, double alphaOut, double L, std::vector<double> d/* Inner to outer diameter*/, std::vector<double> lambda) const
//{
//	double wall = 0;
//	if (lambda.size() == d.size() - 1)
//		for (int i = 0; i < lambda.size(); i++)
//			wall += log(d[i + 1] / d[i]) / lambda[i];
//	double aIn = 0;
//	if ((alphaIn * d.front()) != 0)
//		aIn = 1. / (alphaIn * d.front());
//	double aOut = 0;
//	if ((alphaOut * d.back()) != 0)
//		aIn = 1. / (alphaOut * d.back());
//	double sum = aIn + 1. / 2 * wall + aOut;
//	const double kA = 1. / (1. / MATH_PI / L * sum);
//	return kA;
//	/*
//		lambda of 1.4404 steel = 15 W/m/K at 20°C
//							T°C rho	al	cp	lam
//		S31603				-100	15.0		 [73, 217]
//		Wst-Nr. 1.4404		0		16.0 466
//		X2CrNiMoi17-12-2	20 7956	16.2 470 12.7
//		ASTM/AISI 316L		100		16.7 486 13.8
//							200		17.1 501 15.5
//							400		18.1 518 18.6
//							600		18.8 539 21.7
//							800		19.3 557 24.8
//
//		https://doi.org/10.1007/978-3-662-57572-7_3
//	*/
//}

//double CDryerBatch::CalckAp(double alphaIn, double alphaOut, std::vector<double> A, std::vector<double> delta, std::vector<double> lambda) const
//{
//	double wall = 0;
//	if (lambda.size() == A.size() - 2 && lambda.size() == delta.size())
//		for (int i = 0; i < lambda.size(); i++)
//			wall += delta[i] / (lambda[i] * A[i + 1]);
//	double aIn = 0;
//	if ((alphaIn * A.front()) != 0)
//		aIn = 1. / (alphaIn * A.front());
//	double aOut = 0;
//	if ((alphaOut * A.back()) != 0)
//		aIn = 1. / (alphaOut * A.back());
//	double sum = aIn + wall + aOut;
//	const double kA = 1. / sum;
//	return kA;
//}

//double CDryerBatch::CalcAlphaOutside(double _time, const double h, const double D, const double Ts, EShape shape) const
//{
//	const double Tstar = 0.5 * (Ts + T_inf); // Temperature for properties
//	const double C_PGas = GetAvgTPCompoundProperty(_time, EPhase::GAS, ECompoundTPProperties::HEAT_CAPACITY_CP, Tstar);
//	const double etaGas = GetAvgTPCompoundProperty(_time, EPhase::GAS, ECompoundTPProperties::VISCOSITY, Tstar);
//	const double lambdaGas = GetAvgTPCompoundProperty(_time, EPhase::GAS, ECompoundTPProperties::THERMAL_CONDUCTIVITY, Tstar);
//	const double rhoGas = GetAvgTPCompoundProperty(_time, EPhase::GAS, ECompoundTPProperties::DENSITY, Tstar);
//
//	const double nyGas = etaGas / rhoGas;
//	const double Pr = C_PGas * etaGas / lambdaGas;
//	const double beta = 1 / (Tstar);
//	const double deltaT = abs(Ts - T_inf);
//	const double Gr = beta * STANDARD_ACCELERATION_OF_GRAVITY * deltaT * pow(h, 3) / pow(nyGas, 2);
//	const double Ra = Gr * Pr;
//
//	// ToDo - Add Ra under angle
//	// ToDo - Implement surface line instead of height
//
//	const double f1 = pow(1. + pow(0.492 / Pr, 9. / 16), -16. / 9);
//	const double NuP = pow(0.825 + 0.387 * pow(Ra * f1, 1. / 6), 2);
//	const double Nu = shape == EShape::CYLINDRICAL ? NuP + 0.435 * h / D : NuP;
//	const double alpha = Nu * lambdaGas / h;
//	// Expected range 2.5-25 https://www.sciencedirect.com/topics/engineering/convection-heat-transfer-coefficient
//	return alpha;
//}

//std::pair<double, std::vector<double>> CUnitDAEModel::CalculateChamberHeatLoss(double _time, void* _unit, double* _vars)
//{
//	auto* unit = static_cast<CDryerBatch*>(_unit);
//	const double sectionsFilledWithBed = unit->DetermineSectionsFilledWithBed(_time, _vars[m_iTempParticle]);
//	//	const double sectionsFilledWithBed = unit->CalculateBedHeightOrDetermineSectionsFilledWithBed(_time, _vars[m_iTempParticle],false);
//	double heightUsage, fullSections;
//	heightUsage = modf(sectionsFilledWithBed, &fullSections);
//
//	double Q_PW = 0;
//	std::vector<double> Q_GW;
//
//	// Fully filled sections
//	for (int section = 0; section < fullSections; section++)
//	{
//		std::vector<double> Q_GWsection(unit->chamber.at(section).layers, 0);
//		Q_PW += CalculateSectionHeatLoss(_time, _unit, _vars, section, &Q_GWsection);
//		Q_GW.insert(Q_GW.end(), Q_GWsection.begin(), Q_GWsection.end());
//	}
//
//	// Partilly filled section
//
//	std::vector<double> Q_GWsection(unit->chamber.at(fullSections).layers, 0);
//	Q_PW += CalculateSectionHeatLoss(_time, _unit, _vars, fullSections, &Q_GWsection, heightUsage);
//	Q_GW.insert(Q_GW.end(), Q_GWsection.begin(), Q_GWsection.end());
//
//	// Fully empty sections
//	for (size_t section = std::max(std::ceil(sectionsFilledWithBed), 1.); section < unit->chamber.size(); section++)
//	{
//		std::vector<double> Q_GWsection(unit->chamber.at(section).layers, 0);
//		Q_PW += CalculateSectionHeatLoss(_time, _unit, _vars, section, &Q_GWsection, 0);
//		Q_GW.insert(Q_GW.end(), Q_GWsection.begin(), Q_GWsection.end());
//	}
//
//	// Top plate
//	Q_GW.back() += CalculateTopPlateHeatLoss(_time, _unit, _vars);
//	return std::make_pair(Q_PW, Q_GW);
//}

/// calculate alpha_GW for discretized height, CURRENLTY NOT IN USE
//double CDryerBatch::CalculateAlpha_GW(double _time, size_t section) const
//{
//	const double Re = CalculateReynolds(_time, section);
//	const double Pr = CalculatePrandtl(_time);
//	double di = (chamber.at(section).dimensionsInternal.at(0).first + chamber.at(section).dimensionsInternal.at(0).second) / 2;
//	double Nu_m2 = 1.616 * std::cbrt(Re * Pr * di / chamber.at(section).height);
//	double Nu2300 = std::cbrt(pow(3.66, 3) + pow(0.7, 3) + pow(Nu_m2 - 0.7, 3));
//	double zeta = pow(1.8 * log10(Re) - 1.5, -2) / 8;
//	double Nu1e4 = zeta * Re * Pr * (1 + pow(di / chamber.at(section).height, 2. / 3)) / (1. + 12.7 * sqrt(zeta) * (pow(Pr, 2. / 3) - 1));
//	double Nu = 0;
//	if (Re < 2300)
//		Nu = Nu2300;
//	else if (Re > 1e4)
//		Nu = Nu1e4;
//	else
//	{
//		double gamma = (Re - 2300) / (1e4 - 2300);
//		Nu = (1 - gamma) * Nu2300 + gamma * Nu1e4;
//	}
//	const double alphaGW = Nu * lambdaGas / di;
//	return alphaGW;
//}


///////////////////////////////////////////////////////////////////////////////////
/// Calculation about sections, for height discretization, CURRENTLY NOT IN USE ///
///////////////////////////////////////////////////////////////////////////////////
/// CheckHeightDiscretizationLayers
//void CDryerBatch::CheckHeightDiscretizationLayers(double _time)
//{
//	N_particle = this->GetConstIntParameterValue("N_el");
//	std::vector<double> sizeGrid = GetNumericGrid(DISTR_SIZE);
//	double maxDistributionSize = sizeGrid.back();
//	if (chamber.at(0).height / N_particle < maxDistributionSize)
//	{
//		this->RaiseWarning("Height of discretization layer smaller than biggest particles.");
//		int64_t tempN_el = N_particle -  1;
//		while (chamber.at(0).height / tempN_el < maxDistributionSize)
//		{
//			tempN_el--;
//			if (tempN_el == 0)
//				break;
//		}
//		N_particle = tempN_el;
//		if (N_particle < 1)
//			this->RaiseError("Size of larges particles " + std::to_string(maxDistributionSize) + " m is larger than total bed height " + std::to_string(chamber.at(0).height) + " m.");
//		this->RaiseWarning("Adjusted number of discretization layers to: " + std::to_string(N_particle));
//	}
//	for (auto& element : chamber)
//	{
//		if (element.layers < 1)
//		{
//			element.layers = 1;
//		}
//			
//	}
//		
//	chamber.at(0).layers = N_particle;
//	double sectionsFilledWithBed = 0;
//	if (particlesGlobal)
//	{
//		std::vector<double> inGasStreamTimePoints = m_inGasStream->GetAllTimePoints();
//		double maxInGasFlow = 0;
//		heighestFlowTimepoint = 0;
//		for (auto& element : inGasStreamTimePoints)
//			if (m_inGasStream->GetMassFlow(element) > maxInGasFlow)
//			{
//				maxInGasFlow = m_inGasStream->GetMassFlow(element);
//				heighestFlowTimepoint = element;
//			}
//		sectionsFilledWithBed = DetermineSectionsFilledWithBed(heighestFlowTimepoint, m_holdup->GetTemperature(heighestFlowTimepoint));
//		// sectionsFilledWithBed = CalculateBedHeightOrDetermineSectionsFilledWithBed(heighestFlowTimepoint, std::max(m_holdup->GetTemperature(0), m_inGasStream->GetTemperature(heighestFlowTimepoint)), false);
//		for (int i = 1; i < std::ceil(sectionsFilledWithBed); i++)
//		{
//			chamber.at(i).layers = std::ceil(chamber.at(0).layers * chamber.at(i).height / chamber.at(0).height);
//			N_particle += chamber.at(i).layers;
//		}
//	}
//	N_total = 0;
//	for (auto& section : chamber)
//		N_total += section.layers;
//}

//double CUnitDAEModel::CalculateAverage(double* _vars, size_t variableKey, int64_t end, int64_t start) const
//{
//	double tempAvg = 0;
//	for (int64_t i = start; i < end; i++)
//		tempAvg += _vars[variableKey + i];
//	if (end - start != 0)
//		tempAvg = tempAvg / (end - start);
//	return tempAvg;
//}

//double CDryerBatch::GetAvgConstCompoundProperty(double _time, EPhase phase, ECompoundConstProperties  property) const
//{
//	std::vector<double> CompoundProperties(compoundKeys.size());
//	for (int i = 0; i < compoundKeys.size(); i++)
//		CompoundProperties[i] = GetCompoundProperty(compoundKeys[i], property);
//
//	std::vector<double> PhaseCompoundsDistribution = m_holdup->GetPhase(phase)->GetCompoundsDistribution(_time);
//
//	double avgProperty = 0;
//	for (int i = 0; i < compoundKeys.size(); i++)
//		avgProperty += PhaseCompoundsDistribution[i] * CompoundProperties[i];
//	return avgProperty;
//}

//double CDryerBatch::GetAvgTPCompoundProperty(double _time, EPhase phase, ECompoundTPProperties  property, double temperature, double pressure) const
//{
//	std::vector<double> CompoundProperties(compoundKeys.size());
//	for (int i = 0; i < compoundKeys.size(); i++)
//		CompoundProperties[i] = GetCompoundProperty(compoundKeys[i], property, temperature, pressure);
//
//	std::vector<double> PhaseCompoundsDistribution = m_holdup->GetPhase(phase)->GetCompoundsDistribution(_time);
//
//	double avgProperty = 0;
//	for (int i = 0; i < compoundKeys.size(); i++)
//		avgProperty += PhaseCompoundsDistribution[i] * CompoundProperties[i];
//	return avgProperty;
//}

//double CDryerBatch::CalculateSectionVolume(size_t section)
//{
//	double sectionVolume = 0;
//	if (chamber.at(section).shape == EShape::CYLINDRICAL)
//		sectionVolume = chamber.at(section).height * MATH_PI / 3 * (
//			pow(chamber.at(section).dimensionsInternal.at(0).first / 2, 2)
//			+ (chamber.at(section).dimensionsInternal.at(0).first * chamber.at(section).dimensionsInternal.at(0).second) / 4
//			+ pow(chamber.at(section).dimensionsInternal.at(0).second / 2, 2)
//			);
//	else if (chamber.at(section).shape == EShape::RECTANGULAR)
//	{
//		std::stringstream os;
//		os << section << " has an rectengular shape.\nWhich volume has not been implemented.";
//		RaiseError(os.str());
//	}
//	else
//	{
//		std::stringstream os;
//		os << section << " has an undefind shape.";
//		RaiseError(os.str());
//	}
//
//	return sectionVolume;
//}

/// testing
//void CDryerBatch::Testing()
//{
//	/*std::vector<double> inGasStreamTimePoints = m_inGasStream->GetAllTimePoints();
//	std::vector<double> layerGasMass(N_total,0);
//	std::vector<std::vector<double>> masses(inGasStreamTimePoints.size(), layerGasMass);
//	for (int i = 0; i < inGasStreamTimePoints.size();i++)
//		masses.at(i) = GetGasMassOfLayers(inGasStreamTimePoints.at(i), STANDARD_CONDITION_T, STANDARD_CONDITION_T);
//	std::vector<double> diff1(N_total, 0);
//	std::vector<double> diff2(N_total, 0);
//	for (int i = 0; i < N_total; i++)
//	{
//		diff1.at(i) = masses.at(2).at(i) - masses.at(0).at(i);
//		diff2.at(i) = masses.at(4).at(i) - masses.at(2).at(i);
//	}*/
//}

//double CDryerBatch::CalculateOverallHeatTransferCoefficientCylinder(size_t section, double alphaInternal, double alphaExternal, double heightUsage)
//{
//	if (chamber.at(section).thermalConductivities.size() != chamber.at(section).wallThicknesses.size() || heightUsage > 1 || heightUsage <= -1)
//	{
//		std::ostringstream  os;
//		os << "CalculateOverallHeatTransferCoefficientCylinder has encountered an error in section: ";
//		os << section;
//		os << ".\nCheck if numer of wall thicknesses is equal to thermal conductivities.";
//		os << "\nOr heightUsage was out of bounds.";
//		RaiseError(os.str());
//	}
//	const double tanWallAngle = (chamber.at(section).dimensionsInternal.at(0).second - chamber.at(section).dimensionsInternal.at(0).first) / (2 * chamber.at(section).height);
//
//	double h, R, r, R_Plus_r;
//
//	if (heightUsage > 0)
//	{
//		h = heightUsage * chamber.at(section).height;
//		r = chamber.at(section).dimensionsInternal.at(0).first / 2;
//		R = h * tanWallAngle + r;
//		R_Plus_r = (R + r);
//	}
//	else
//	{
//		h = (1 - heightUsage) * chamber.at(section).height;
//		R = chamber.at(section).dimensionsInternal.at(0).second / 2;
//		r = R - h * tanWallAngle;
//		R_Plus_r = (R + r);
//	}
//
//	double wall = 0;
//	for (int i = 0; i < chamber.at(section).thermalConductivities.size(); i++)
//	{
//		double sumWallThicknesses = 0;
//		for (int j = 0; j < chamber.at(section).wallThicknesses.size() - 1; j++)
//			sumWallThicknesses += chamber.at(section).wallThicknesses.at(j);
//		wall += log((R_Plus_r + 2 * (chamber.at(section).wallThicknesses.at(i) + sumWallThicknesses)) / (R_Plus_r + sumWallThicknesses)) / chamber.at(section).thermalConductivities.at(i);
//	}
//	double sumWallThicknesses = 0;
//	for (int i = 0; i < chamber.at(section).wallThicknesses.size(); i++)
//		sumWallThicknesses += chamber.at(section).wallThicknesses.at(i);
//	double aIn = 0;
//	if (alphaInternal * R_Plus_r != 0)
//		aIn = 1. / (alphaInternal * R_Plus_r);
//	double aOut = 0;
//	if (alphaExternal * (R_Plus_r + sumWallThicknesses) != 0)
//		aIn = 1. / (alphaExternal * (R_Plus_r + sumWallThicknesses));
//	double sum = aIn + 1. / 2 * wall + aOut;
//
//	const double m = sqrt(pow(R - r, 2) + pow(h, 2));
//	const double kA = 1. / (sum / (MATH_PI * m));
//	return kA;
//}

//double CDryerBatch::CalculateOverallHeatTransferCoefficientTopPlate(double alphaInternal, double alphaExternal)
//{
//	if (chamber.back().thermalConductivities.size() != chamber.back().wallThicknesses.size())
//		bool error = true;
//	double A = 0;
//	if (chamber.back().shape == EShape::CYLINDRICAL)
//		A = pow(chamber.back().dimensionsInternal.at(0).second, 2) * MATH_PI / 4;
//	else
//		A = chamber.back().dimensionsInternal.at(0).second * chamber.back().dimensionsInternal.at(1).second;
//
//	double wall = 0;
//	for (int i = 0; i < chamber.back().thermalConductivities.size(); i++)
//		wall += chamber.back().wallThicknesses.at(i) / (chamber.back().thermalConductivities.at(i) * A);
//	double aIn = 0;
//	if ((alphaInternal * A) != 0)
//		aIn = 1. / (alphaInternal * A);
//	double aOut = 0;
//	if ((alphaExternal * A) != 0)
//		aIn = 1. / (alphaExternal * A);
//	double sum = aIn + wall + aOut;
//	if (sum == 0)
//		bool error = true;
//	const double kA = 1. / sum;
//	return kA;
//}

//std::vector<double> CDryerBatch::GetSectionGasMass(double _time, double gasTemperature, double particleTemperature)
//{
//	std::vector<double> layerMasses(N_total, 0);
//	double chamberVolume = 0;
//	for (size_t section = 0; section < chamber.size(); section++)
//		chamberVolume += CalculateSectionVolume(section);
//
//	double rhoParticle = GetAvgTPCompoundProperty(_time, EPhase::SOLID, ECompoundTPProperties::DENSITY, particleTemperature, m_holdup->GetPressure(_time));
//	double rhoGas = GetAvgTPCompoundProperty(_time, EPhase::GAS, ECompoundTPProperties::DENSITY, gasTemperature, m_holdup->GetPressure(_time));
//
//	const double particlesVolume = m_holdup->GetPhaseMass(_time, EPhase::SOLID) / rhoParticle;
//	const double chamberGasMass = (chamberVolume - particlesVolume) * rhoGas;
//	double sectionsFilledWithBed = DetermineSectionsFilledWithBed(_time, particleTemperature);
//	//	double sectionsFilledWithBed = CalculateBedHeightOrDetermineSectionsFilledWithBed(_time, particleTemperature,false);
//	double heightUsage, fullSections;
//	heightUsage = modf(sectionsFilledWithBed, &fullSections);
//	double sectionVolume = 0;
//	double sectionGasMass = 0;
//	std::vector<double>::reverse_iterator ritLayerMasses = layerMasses.rbegin();
//	std::vector<chamberSection>::reverse_iterator ritChamberSections = chamber.rbegin();
//	double assosiatedGasMass = 0;
//	ritChamberSections++;
//	double test = chamber.rend() - ritChamberSections;
//
//	// Empty section above bed
//	for (size_t section = chamber.size() - 1; section > fullSections; section--)
//	{
//		sectionVolume = CalculateSectionVolume(section);
//		sectionGasMass = sectionVolume * rhoGas;
//		std::vector<double> sectionLayerGasMasses(chamber.at(section).layers, 0);
//		const double r = chamber.at(section).dimensionsInternal.at(0).first / 2;
//		const double R = chamber.at(section).dimensionsInternal.at(0).second / 2;
//		const double h = chamber.at(section).height;
//		const double tanWallAngle = (R - r) / h;
//		for (size_t layer = 0; layer < chamber.at(section).layers; layer++)
//		{
//			const double r1 = r + h * (layer) / chamber.at(section).layers * tanWallAngle;
//			const double R1 = r + h * (layer + 1) / chamber.at(section).layers * tanWallAngle;
//			const double layerVolume = CalculateLayerVolume(section, R1, r1);
//			sectionLayerGasMasses.at(layer) = sectionGasMass * layerVolume / sectionVolume;
//		}
//		std::vector<double>::reverse_iterator ritSectionLayerMasses = sectionLayerGasMasses.rbegin();
//		for (; ritSectionLayerMasses != sectionLayerGasMasses.rend(); ritSectionLayerMasses++)
//		{
//			*ritLayerMasses = *ritSectionLayerMasses;
//			ritLayerMasses++;
//			assosiatedGasMass += *ritSectionLayerMasses;
//		}
//
//	}
//	// Should be working.
//
//
//	size_t layersWithParticles = 0;
//
//
//	// Last section filled with bed
//	{
//		size_t section = fullSections;
//		size_t layersUsed = std::ceil(heightUsage * chamber.at(section).layers);
//		layersWithParticles += layersUsed;
//		for (size_t i = 0; i < fullSections; i++)
//			layersWithParticles += chamber.at(i).layers;
//
//		std::vector<double> sectionLayerGasMasses(chamber.at(section).layers, 0);
//		const double r = chamber.at(section).dimensionsInternal.at(0).first / 2;
//		const double R = chamber.at(section).dimensionsInternal.at(0).second / 2;
//		const double h = chamber.at(section).height;
//		const double tanWallAngle = (R - r) / h;
//
//		for (size_t layer = layersUsed; layer < chamber.at(section).layers; layer++)
//		{
//			const double r1 = r + h * (layer) / chamber.at(section).layers * tanWallAngle;
//			const double R1 = r + h * (layer + 1) / chamber.at(section).layers * tanWallAngle;
//			const double layerVolume = CalculateLayerVolume(section, R1, r1);
//			sectionLayerGasMasses.at(layer) = sectionGasMass * layerVolume / sectionVolume;
//		}
//
//		for (size_t layer = 0; layer < layersUsed; layer++)
//		{
//			const double r1 = r + h * (layer) / chamber.at(section).layers * tanWallAngle;
//			const double R1 = r + h * (layer + 1) / chamber.at(section).layers * tanWallAngle;
//			const double layerVolume = CalculateLayerVolume(section, R1, r1) - (particlesVolume / layersWithParticles);
//			sectionLayerGasMasses.at(layer) = sectionGasMass * layerVolume / sectionVolume;
//		}
//
//		std::vector<double>::reverse_iterator ritSectionLayerMasses = sectionLayerGasMasses.rbegin();
//		for (; ritSectionLayerMasses != sectionLayerGasMasses.rend(); ritSectionLayerMasses++)
//		{
//			*ritLayerMasses = *ritSectionLayerMasses;
//			ritLayerMasses++;
//			assosiatedGasMass += *ritSectionLayerMasses;
//		}
//	}
//
//	// Sections filled completely with particles
//	if (fullSections > 0)
//		for (double section = fullSections - 1; section > -1; section--)
//		{
//			sectionVolume = CalculateSectionVolume(section);
//			sectionGasMass = sectionVolume * rhoGas;
//			std::vector<double> sectionLayerGasMasses(chamber.at(section).layers, 0);
//			const double r = chamber.at(section).dimensionsInternal.at(0).first / 2;
//			const double R = chamber.at(section).dimensionsInternal.at(0).second / 2;
//			const double h = chamber.at(section).height;
//			const double tanWallAngle = (R - r) / h;
//			for (size_t layer = 0; layer < chamber.at(section).layers; layer++)
//			{
//				const double r1 = r + h * (layer) / chamber.at(section).layers * tanWallAngle;
//				const double R1 = r + h * (layer + 1) / chamber.at(section).layers * tanWallAngle;
//				const double layerVolume = CalculateLayerVolume(section, R1, r1) - (particlesVolume / layersWithParticles);
//				sectionLayerGasMasses.at(layer) = sectionGasMass * layerVolume / sectionVolume;
//			}
//			std::vector<double>::reverse_iterator ritSectionLayerMasses = sectionLayerGasMasses.rbegin();
//			for (; ritSectionLayerMasses != sectionLayerGasMasses.rend(); ritSectionLayerMasses++)
//			{
//				*ritLayerMasses = *ritSectionLayerMasses;
//				ritLayerMasses++;
//				assosiatedGasMass += *ritSectionLayerMasses;
//			}
//
//		}
//
//	if (chamberGasMass - assosiatedGasMass > 1e-16)
//		RaiseError("Not all gas mass was accounted for in combined layer gas mass.");
//
//	return layerMasses;
//}

//double CUnitDAEModel::CalculateSectionHeatLoss(double _time, void* _unit, double* _vars, size_t section, std::vector<double>* Q_GW, double heightUsage)
//{
//	auto* unit = static_cast<CDryerBatch*>(_unit);
//
//	bool noTemperatureDelta = unit->particlesGlobal ? ( (_vars[m_iTempParticle] - unit->T_inf) != 0 ? false : true) : true;
//	if (noTemperatureDelta)
//		for (int i = 0; i < Q_GW->size(); i++)
//			if (_vars[m_iTempOutGas + (i + section * Q_GW->size())] - unit->T_inf != 0)
//			{
//				noTemperatureDelta = false;
//				break;
//			}
//	if (noTemperatureDelta)
//		return 0;
//
//	const double alphaGPipe = unit->CalculateAlpha_GW(_time, section);
//	//const double alphaGFB = unit->CalculateAlpha_GW(_time, );
//	size_t layersOfPrevSections = 0;
//	for (size_t i = 0; i < section; i++)
//	{
//		layersOfPrevSections += unit->chamber.at(i).layers;
//	}
//	const double varAvTempGas = CalculateAverage(_vars, m_iTempOutGas, unit->chamber.at(section).layers+layersOfPrevSections, layersOfPrevSections);
//	const double alphaPFB = !unit->particlesGlobal ? 0 : unit->CalculateAlpha_PW(_vars[m_iTempParticle], varAvTempGas, unit->m_holdup->GetPressure(_time), _time);
//
//	double sumWallThicknesses = 0;
//	for (int i = 0; i < unit->chamber.at(section).wallThicknesses.size(); i++)
//		sumWallThicknesses += unit->chamber.at(section).wallThicknesses.at(i);
//	double outerDiameter = sumWallThicknesses + (unit->chamber.at(section).dimensionsInternal.at(0).first + unit->chamber.at(section).dimensionsInternal.at(0).second) / 2;
//	double Q_PW = 0;
//
//	double alphaOut = 0;
//	double kA_P = 0;
//	double kA_G = 0;
//
//	size_t usedLayers = unit->DetermineLayersInSectionFilledWithBed(section, heightUsage);
//	double usedLayerPercentage = 1. / unit->chamber.at(section).layers * usedLayers;
//	size_t layersOfPreveusSections = 0;
//	for (size_t i = 0; i < section; i++)
//		layersOfPreveusSections += unit->chamber.at(i).layers;
//
//	double temperatureWall0 = 0;
//	double temperatureWall = 0;
//
//	if (usedLayerPercentage > 0)
//	{
//		temperatureWall0 = 0;
//		temperatureWall = 0.5 * (std::max(unit->m_inGasStream->GetTemperature(_time), unit->m_holdup->GetTemperature(_time)) - unit->T_inf) + unit->T_inf;
//		while (abs(temperatureWall - temperatureWall0) > 1/*temperatureWall * GetRTol() * 10 + GetATol(unit->N_total)*/)
//		{
//			alphaOut = unit->CalcAlphaOutside(_time, unit->chamber.at(section).height, outerDiameter, temperatureWall);
//
//			kA_P = unit->CalculateOverallHeatTransferCoefficientCylinder(section, alphaPFB, alphaOut, usedLayerPercentage);
//			Q_PW = unit->particlesGlobal ? kA_P * (_vars[m_iTempParticle] - unit->T_inf) : 0;
//
//			//kA_G = unit->CalculateOverallHeatTransferCoefficientCylinder(section, alphaGFB, alphaOut, usedLayerPercentage);
//			for (int i = 0; i < usedLayers; i++)
//				Q_GW->at(i) = kA_G / usedLayers * (_vars[m_iTempOutGas + (i + layersOfPreveusSections)] - unit->T_inf);
//
//			double Q = Q_PW;
//			for (int i = 0; i < usedLayers; i++)
//				Q += Q_GW->at(i);
//
//			temperatureWall0 = temperatureWall;
//			double alphaAout = unit->chamber.at(section).height * outerDiameter * MATH_PI * alphaOut;
//			temperatureWall = Q / alphaAout + unit->T_inf;
//		}
//	}
//	temperatureWall0 = 0;
//	temperatureWall = 0.5 * (std::max(unit->m_inGasStream->GetTemperature(_time), unit->m_holdup->GetTemperature(_time)) - unit->T_inf) + unit->T_inf;
//	while (abs(temperatureWall - temperatureWall0) > /*temperatureWall * GetRTol() * 10 + GetATol(unit->N_total)*/ 1)
//	{
//		alphaOut = unit->CalcAlphaOutside(_time, unit->chamber.at(section).height, outerDiameter, temperatureWall);
//
//		kA_G = unit->CalculateOverallHeatTransferCoefficientCylinder(section, alphaGPipe, alphaOut, -usedLayerPercentage);
//		for (size_t i = usedLayers; i < unit->chamber.at(section).layers; i++)
//			Q_GW->at(i) = kA_G / (unit->chamber.at(section).layers -usedLayers) * (_vars[m_iTempOutGas + (i + layersOfPreveusSections)] - unit->T_inf);
//
//		double Q = 0;
//		for (size_t i = usedLayers; i < unit->chamber.at(section).layers; i++)
//			Q += Q_GW->at(i);
//
//		temperatureWall0 = temperatureWall;
//		double alphaAout = unit->chamber.at(section).height * outerDiameter * MATH_PI * alphaOut;
//		temperatureWall = Q / alphaAout + unit->T_inf;
//	}
//
//	return Q_PW;
//}

//double CUnitDAEModel::CalculateTopPlateHeatLoss(double _time, void* _unit, double* _vars)
//{
//	auto* unit = static_cast<CDryerBatch*>(_unit);
//
//	if (_vars[m_iTempOutGas + (unit->N_total - 1)] - unit->T_inf == 0)
//		return 0;
//
//	const double alphaG = unit->CalculateAlpha_GW(_time, unit->chamber.size() - 1);
//	double temperatureWall0 = 0;
//	double temperatureWall = 0.5 * (std::max(unit->m_inGasStream->GetTemperature(_time), unit->m_holdup->GetTemperature(_time)) - unit->T_inf) + unit->T_inf;
//	double Q_GWTop = 0;
//	double alphaOut = 0;
//	double kA_G = 0;
//
//	while (abs(temperatureWall - temperatureWall0) > 1 /*temperatureWall * GetRTol()*10 + GetATol(unit->N_total)*/)
//	{
//		alphaOut = unit->CalcAlphaOutside(_time, unit->chamber.back().dimensionsInternal.at(0).second, 0, temperatureWall, EShape::RECTANGULAR);
//
//		kA_G = unit->CalculateOverallHeatTransferCoefficientTopPlate(alphaG, alphaOut);
//		Q_GWTop = kA_G  * (_vars[m_iTempOutGas + (unit->N_total-1)] - unit->T_inf);
//
//		temperatureWall0 = temperatureWall;
//		double alphaAout = pow(unit->chamber.back().dimensionsInternal.at(0).second,2)/4 * MATH_PI * alphaOut;
//		temperatureWall = Q_GWTop / alphaAout + unit->T_inf;
//	}
//	return Q_GWTop;
//}