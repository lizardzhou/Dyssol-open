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
/// - With heat transfer based on model from Maksym + Soeren (1st and 2nd drying stages)
/// - Mass transfer of water from particle to fluidization gas
/// - Mass transfer is treated as ideal plug flow reactor and heat exchanger with phase change -> consider the efficiency during mass transfer

void CDryerBatch::CreateBasicInfo()
{
	/// Set basic unit info ///
	SetUnitName("Dryer (batch with heat and mass transfer)");
	SetAuthorName("Alexander Hanke, Xiye Zhou");
	SetUniqueID("E458E673-D296-4CC0-9084-69381811A0A1");
}

void CDryerBatch::CreateStructure()
{
	/// Add ports ///
	AddPort("InletLiquid", EUnitPort::INPUT);
	AddPort("InletFluidizationGas", EUnitPort::INPUT);
	AddPort("InletNozzleAir", EUnitPort::INPUT); 
	AddPort("OutletExhaustGas", EUnitPort::OUTPUT);

	/// Add unit parameters ///
	AddStringParameter("Particle properties", "",	"");
	AddConstRealParameter("A_P", 4, "m2", "Total surface of all particles in the granulator.\nIf <= 0, PSD is used to calculate surface area", 0);
	AddCheckBoxParameter("updateA", true, "Tick this box to enable calculating total particle surface from PSD and mass.");
	AddConstRealParameter("Delta_f"	, 40	, "mum"			, "Thickness of the water film on particles in micrometer"	, 1);
	//AddStringParameter("Path", "C:\\", "Location of dimension file.");
	AddConstRealParameter("w_l (t=0)", 0, "mass %", "Initial moisture fraction of particles.\nIf 0, equalibrium moisture content of ambient conditions will be used.", 0, 100);
	AddConstRealParameter("w_l,eq,min", 0, "mass %", "Minimum measured equilibirum moisture fraction of particles.\nIf 0, materials database will be used to calculate euqilibirum moisture content. (Further see Path X_eq)", 0, 100);
	AddConstRealParameter("theta_eq,min", 0, "degree celsius", "Temperature correponding to w_l,eq,min", 0, 100);
	AddStringParameter("Path to X_eq data", "C:\\", "Location of equilibrium moisture content with temperature, must be a csv file.");
	AddCompoundParameter("X_eq compound", "Compound for with the Xeq values are contained in Path Xeq.");
	AddConstRealParameter("k_dc", 3.5, "-", "k for normalized drying curve. \nthe normalized drying curve (nu) of the material is represented by: \n nu = k* eta / (1. + eta*(k - 1.)), \n with eta representing the particle moiture content.", 0, 5); // Credit to SetUnitName  ("Vibrated Fluidized Bed Dryer (steady-state)");	SetAuthorName("Buchholz (based on Zhengyu Lu's Master's Thesis)");
	AddConstRealParameter("X_cr", 0.025, "kg/kg", "Critical water content of the particles (transition between 1st and 2nd drying period).", 0);
	AddConstRealParameter("theta_env", 20, "deg C", "Temperature of the environment, assumed to be constant", 0, 100);

	// Height of temperature probe and nozzle
	AddStringParameter("Inlet fluidization gas properties", "", "");
	AddCheckBoxParameter("Y_eq", 1, "If unticked, omits reduced saturation vapor pressure.");
	AddConstRealParameter("Y_in"	, 5		, "g/kg dry gas", "Moisture content of the fluidization gas"				, 0);
	AddConstRealParameter("RH_in"	, 48.99	, "%"			, "Relative humidity of the fluidization gas"				, 0, 100);
	AddConstRealParameter("Y_sat"	, 0		, "g/kg dry gas", "Saturation moisture content of the fluidization air\nIf 0, material database will be used to cacluate Y_sat.", 0);
	//AddCheckBoxParameter("calc Y_sat", true, "Tick this box to enable calculating saturated absolute moisture in the air.");
	AddConstRealParameter("beta_GP", -1, "m/s", "Mass transfer coefficient for liquid from gas to particle\nIf 0, calculated using correlation from Gnielinski");
	//AddCheckBoxParameter("calc Beta"	, true, "Tick this box to enable calculating mass transfer coeffcient.");
	
	// Selection of correlation for calculating diffusion coefficient
	std::vector<size_t>	items{ 0, 1, 2 };
	std::vector<std::string> itemNames{ "Dosta", "Tsotsas", "Poos" };
	AddComboParameter("Diff_coeff", 0, items, itemNames, "Correlation for diffusion coefficient from water vapor to air"); // by default use correlation from Dosta (2010)
	// Selection of calculating method for calculating inlet air absolute moisture
	itemNames = { "RH_in","Y_in","Input stream" };
	AddComboParameter("Y_in_Value", 0, items, itemNames, "How Y_in should be calculated.");
	// Selection of calculating method for drying kinetics
	itemNames = { "REA", "Normalized drying curve", "No curve" };
	AddComboParameter("DryingCurve", 0, items, itemNames, "How drying curve should be determined.");
	// Reaction Engineering Approach for drying curve
	AddStringParameter("REA parameters", "", "Values used for funtion fit for Reaction Engineering Approach to drying phases\n Normalized activation energy = REA1 * exp(REA2 * (X - Xeq)^REA3)");
	AddConstRealParameter("REA1", 0.96, "", "");
	AddConstRealParameter("REA2", -19.63, "", "");
	AddConstRealParameter("REA3", 0.73, "", "");
	
	// Process chamber information
	AddStringParameter	 ("Chamber information"	, ""		, "");
	AddConstRealParameter("H_tempProbe", 7, "cm", "Height of temperature probe for chamber temperature over distributor.", 0);
	AddConstRealParameter("H_nozzle", 25, "cm", "Height of two-phase nozzle over distributor.", 0);
	AddConstRealParameter("d_bed"			, 0.2			, "m"	, "Bed diameter", 0.001);
	AddConstRealParameter("H_bed"			, 0.106			, "m"	, "Bed height without fluidization", 0.001);
	AddConstRealParameter("H_chamber"		, 0.25			, "m"	, "Process chamber height", 0.001);
	// Bed properties, use for development of further complexer models
	AddStringParameter("Bed properties", ""	, "Currently not in use, designed for further development");
	AddConstRealParameter("eps_0"			, 0.4			, "-"	, "Bed porosity without fluidization"										, 0, 1);
	AddConstRealParameter("u_mf"			, 0				, "m/s"	, "Minimal fluidization velocity\nIf 0, calculate use Wen&Yu correlation", 0, 1);
	AddConstIntParameter ("N_el"			, 1				, ""	, "Number of hight discretization layers"									, 1);
	//AddConstRealParameter("phi_eq",0.4,"","Targeted relativ humidity for exhaust gas",0.1,1);
	
	// Tolerance for solver
	AddStringParameter("Tolerance for solver", "", "");
	AddConstRealParameter("Relative tolerance", 0.0, "-", "Solver relative tolerance. Set to 0 to use flowsheet-set value", 0);
	AddConstRealParameter("Absolute tolerance Y", 0.0, "-", "Solver absolute tolerance for gas moisture content Y.\nSet to 0 to use flowsheet-set value", 0);
	AddConstRealParameter("Absolute tolerance T", 0.0, "-", "Solver absolute tolerance for eqData.temperatures T.\nSet to 0 to use flowsheet-set value", 0);

	// Debug information
	AddStringParameter("Debug information","","");
	AddCheckBoxParameter("Toggle debug mode", 0, "Tick this box to generate debug console outputs.");

	/// Add holdups ///
	AddHoldup("Holdup");

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

	debugToggle = m_model.debugToggle = GetCheckboxParameterValue("Toggle debug mode");

	// Settings
	/*bool calcBeta = GetCheckboxParameterValue("calcBeta");
	bool calcY_sat = GetCheckboxParameterValue("calcY_sat");
	bool calcNdc = GetCheckboxParameterValue("calcNdc");*/
	////ignoreVaporInput = GetCheckboxParameterValue("ignoreVaporInput");
	//useREA = GetCheckboxParameterValue("useREA");
	//dryingCurveSetting = GetComboParameterValue("DryingCurve");


	/// Get holdup ///
	m_holdup = GetHoldup("Holdup");
	const double M_tot = m_holdup->GetPhaseMass(_time, EPhase::SOLID);
	if (M_tot == 0) // empty run
	{
		particlesGlobal = false;
		RaiseWarning("No particles in system.\nSolids and liquids will be ignored.");
	}
	//HoldupOld = m_holdup->GetPhaseMass(_time, EPhase::LIQUID);
	TempLiquidOld = m_holdup->GetTemperature(_time);
	TempSolidOld = m_holdup->GetTemperature(_time);
	TempGasOld = m_holdup->GetTemperature(_time);

	/// Get pointers to streams ///
	m_inLiquidStream = GetPortStream("InletLiquid");
	m_inGasStream = GetPortStream("InletFluidizationGas");
	m_inNozzleAirStream = GetPortStream("InletNozzleAir");
	m_outExhaustGasStream = GetPortStream("OutletExhaustGas");

	/// Pull compound data ///
	this->PullCompoundDataFromDatabase(_time);

	// Setup vapor (gas phase) stream //
	m_VaporStream = AddStream("Vapor");
	std::vector<double> vaporStreamCompoundSetup(size(compoundKeys), 0);
	m_VaporStream->SetCompoundsFractions(_time, vaporStreamCompoundSetup);
	m_VaporStream->SetPhaseFraction(_time, EPhase::SOLID, 0);
	m_VaporStream->SetPhaseFraction(_time, EPhase::LIQUID, 0);
	m_VaporStream->SetPhaseFraction(_time, EPhase::VAPOR, 1);
	m_VaporStream->SetCompoundFraction(_time, compoundKeys.at(indicesOfVaporOfPhaseChangingCompound.second), EPhase::VAPOR,1);
	m_VaporStream->SetMassFlow(_time, 0);

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
	
	double H_bed = GetConstRealParameterValue("H_bed");
	double H_chamber = GetConstRealParameterValue("H_chamber");
	double d_bed = GetConstRealParameterValue("d_bed");
	T_inf = T_ref + GetConstRealParameterValue("T_env");
	eps_0 = GetConstRealParameterValue("eps_0");
	u_mf = GetConstRealParameterValue("u_mf");
	double H_susp = GetConstRealParameterValue("SuspNozzleHeight") / 100;
	SetupChamber();
	eqData.compoundKey = GetCompoundParameterValue("Xeq compound");
	if (particlesGlobal)
	{
		InitializeMoistureContentDatabase(GetStringParameterValue("Path Xeq")); //"E:\\Dyssol\\Xeq.csv"
	}

	//phi_eq = GetConstRealParameterValue("phi_eq");
	k_dc = GetConstRealParameterValue("k_dc");
	X_cr = GetConstRealParameterValue("X_cr");
	DiffCoeff = GetComboParameterValue("Diff_coeff");
	double RH_in = GetConstRealParameterValue("RH_in") / 100; 
	size_t Y_in_Value = GetComboParameterValue("Y_in_Value");
	Y_sat = GetConstRealParameterValue("Y_sat") * 1e-3; // convert into [kg/kg]
	m_model.Y_eq = GetCheckboxParameterValue("Y_eq");
	Delta_f = GetConstRealParameterValue("Delta_f") * 1e-6; // convert into [m]
	heightOfChamberTemperatureProbe = GetConstRealParameterValue("T-Probe height") * 1e-2; // convert into  [m]
	std::ostringstream os;

	/// Read input parameters ///
	// Inlet air absolute moisture
	switch (Y_in_Value) 
	{
		case 0: // calculate from RH and saturation vapor pressure at the given temperature
		{
			Y_in = RH_in * GetGasSaturationMoistureContent(T_inf);
			os << "Y_in is calculated using RH_in as: " << Y_in << " kg/kg";
			ShowInfo(os.str());
			os.str("");
			break;
		}
			
		case 1: // same as GUI input
		{
			Y_in = GetConstRealParameterValue("Y_in") * 1e-3; // convert in [kg/kg]
			os << "Y_in is the input value: " << Y_in << " kg/kg";
			ShowInfo(os.str());
			os.str("");
			break;
		}
		case 2: // using mass fraction of water defined in GUI input of holdup
		{
			double compoundFractionOfVaporOfPhaseChangingCompound = m_inGasStream->GetPhase(EPhase::GAS)->GetCompoundFraction(_time, GetCompoundIndex(compoundKeys[indicesOfVaporOfPhaseChangingCompound.second]));
			Y_in = compoundFractionOfVaporOfPhaseChangingCompound / (1 - compoundFractionOfVaporOfPhaseChangingCompound);
			os << "Y_in is calculated using defined water vapor fraction in the air as: "<< Y_in << " kg/kg";
			ShowInfo(os.str());
			os.str("");
			break;
		}
	}

	/*
	double compoundFractionOfSolidInLiquidInput = m_inSuspensionStream->GetPhase(EPhase::SOLID)->GetCompoundFraction(_time, GetCompoundIndex(compoundKeys[indicesOfVaporOfPhaseChangingCompound.second]));
	if (compoundFractionOfSolidInLiquidInput > 0 && !ignoreVaporInput)
	{
		w_l = compoundFractionOfVaporOfPhaseChangingCompound / (1 - compoundFractionOfVaporOfPhaseChangingCompound);
		ShowInfo("Y_in was calculated using inGasStream as: " + std::to_string(Y_in) + "kg/kg");
	}
	else
		w_l = 1;*/
	
	// Particle properties
	std::vector<double> Grid = GetNumericGrid(DISTR_SIZE);
	std::vector<double> q_3 = m_holdup->GetPSD(_time, PSD_q3, EPSDGridType::DIAMETER);
	const double A_Calculated = GetSpecificSurface(Grid, q_3) / rhoParticle * M_tot;
	bool updateA = GetCheckboxParameterValue("updateA");
	double A_Input = GetConstRealParameterValue("A_P");
	if (updateA || A_Input <= 0)
	{
		A_P = A_Calculated;
	}
	else
	{
		A_P = A_Input;
	}
	/*if (A_Input > 0)
		A_P = A_Input;
	else
		A_P = A_Calculated;*/

	// Show particle properties in DEBUG mode
	if (debugToggle) 
	{
		/// Get number of classes for PSD ///
		size_t m_classesNum = GetClassesNumber(DISTR_SIZE); //n
		/// Get grid of PSD ///
		std::vector<double> m_sizeGrid = GetNumericGrid(DISTR_SIZE); //d_min
		std::vector<double> m_avgDiam = GetClassesMeans(DISTR_SIZE); //d_m,i
		std::vector<double> m_classSize = GetClassesSizes(DISTR_SIZE); //Delta d
		const double d32 = CalculateHoldupSauter(_time);

		// output of read data in simulation window
		os << "\nPSD info\n" << "\tNumber of classes: " << m_classesNum <<"\n\tSize Grid\n";
		for (int i = 0; i < m_sizeGrid.size(); i++)
			os << "\t\tSize Grid: " << i << " = " << m_sizeGrid[i] << " m\n";
		os << "\tAvg diam\n";
		for (int i = 0; i < m_avgDiam.size(); i++)
			os << "\t\tAvg diam: "<<i<<" = "<<m_avgDiam[i]<<" m\n";
		os <<"\tClass size\n";
		for (int i = 0; i < m_classSize.size(); i++)
			os << "\t\tClass size: "<<i<<" = "<<m_classSize[i]<<" m\n";
		os << "\tSauter diameter:" << d32 << "m\n";
		os << "\tToal surface area: " << A_Calculated << " m^2\n" << "\tTotal particle mass: " << M_tot << " kg\n";
		if (particlesGlobal)
			ShowInfo(os.str());
		os.str("");
	} 

	this->CheckHeightDiscretizationLayers(_time);

	// Calculate gas mass from initial conditions and volume // ToDo: Add more info
	volume volumeSolids = m_holdup->GetPhaseMass(_time, EPhase::SOLID) / rhoParticle;
	volume volumeChamber = 0;
	for (size_t i = 0; i < chamber.size(); i++)
		volumeChamber += CalculateSectionVolume(i);
	volume volumeGas = volumeChamber - volumeSolids;
	mass massGas = volumeGas * rhoGas;

	//mass massGas2 = 0;
	//for (size_t i = 0; i < N_total; i++)
	//	massGas2 += GetGasMassOfLayer(_time, i, STANDARD_CONDITION_T, STANDARD_CONDITION_T);

	mass orignalMassGas = m_holdup->GetPhaseMass(_time, EPhase::GAS);

	if (massGas > orignalMassGas)
	{
		m_holdup->SetPhaseMass(_time, EPhase::GAS, massGas);
		os << "Adjusted mass of gas phase in holdup from " << orignalMassGas << " kg to " << massGas << " kg.";
		ShowInfo(os.str());
		os.str("");
	}

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

	suspLayer = DetermineLayersInSectionFilledWithBed(0, H_susp / chamber.at(0).height) - 1;

	/// Clear all state variables in model ///
	m_model.ClearVariables();
	if (particlesGlobal)
	{
		AddStateVariable("X", 0);
		AddStateVariable("w%", 0);
		AddStateVariable("phi", 0);
	}
	AddStateVariable("T_out", m_holdup->GetTemperature(_time));
	AddStateVariable("Y_out", Y_in);
	AddStateVariable("RH_out [%]", RH_in*100);
	if (particlesGlobal)
		AddStateVariable("solid mass in holdup", m_holdup->GetPhaseMass(_time, EPhase::SOLID));
	if (particlesGlobal)
		AddStateVariable("liquid mass in holdup", m_holdup->GetPhaseMass(_time, EPhase::LIQUID));
	AddStateVariable("gas mass in holdup", m_holdup->GetPhaseMass(_time, EPhase::VAPOR));
	if (particlesGlobal)
		AddStateVariable("evaporation rate [kg/s]", 0);
	//AddStateVariable("temperature exhaust gas", m_holdup->GetTemperature(_time) - 273.15); // convert in Cesius
	//AddStateVariable("temperature water film", m_inSuspensionStream->GetTemperature(_time) - 273.15); // convert in Cesius
	//AddStateVariable("temperature particle", m_holdup->GetTemperature(_time) - 273.15); // convert in Cesius
	

	//ShowInfo("RH to Y: " + std::to_string(RH_in * GetGasSaturationMoistureContent(T_ref + 20.5)*1000) + " g/kg");
	//ShowInfo("Massflow nozzle: " + std::to_string((5.91 /*Volume flow gas m^3/h*/ * rhoGas / 60 * 1000 + 20 /*Mass flow water g/min*/) / 60 / 1000) + " kg/s");
	//ShowInfo("Phase fraction water: " + std::to_string(20. /*Mass flow water g/min*/ / 60 / 1000 / ((5.91 /*Volume flow gas m^3/h*/ * rhoGas / 60 * 1000 + 20 /*Mass flow water g/min*/) / 60 / 1000)));

	//temperature T_out = T_ref + 48.2;
	//ShowInfo("RH - in vs out: " + std::to_string(RH_in*100) + " at " + std::to_string(T_inf-T_ref) + " deg C turn into " + std::to_string(100 * GetRelativeHumidity(RH_in * GetGasSaturationMoistureContent(T_inf), T_out)) + " at " + std::to_string(T_out-T_ref) + " deg C.");
	////ShowInfo(std::to_string(_y / GetGasSaturationMoistureContent(_t, _p)));

	if (particlesGlobal)
	{
		double w_eqMinTemperature = GetConstRealParameterValue("w_l,eq,min temperature") + T_ref;
		double X_eqMinBase = CalcuateSolidEquilibriumMoistureContent(_time, w_eqMinTemperature, GetRelativeHumidity(Y_in, w_eqMinTemperature));
		double w_leqmin = GetConstRealParameterValue("w_l,eq,min") / 100;
		//minMoistureContent = w_leqmin / (1. - w_leqmin);//0.0271;
		//moistureScaler = minMoistureContent / X_eqMinBase;
	}

	// State variables for model of Dosta
	temperature tempParticleInit = 0;
	double phiInit = 0;
	if (particlesGlobal)
	{
		tempParticleInit = m_holdup->GetTemperature(_time);
		if (tempParticleInit <= 0)
			RaiseError("Provided holdup temperature is <= 0 K.");
		const double w_lInit = GetConstRealParameterValue("w_l(t=0)") / 100;
		if (w_lInit < 0)
			InitX = CalcuateSolidEquilibriumMoistureContent(_time, tempParticleInit, GetRelativeHumidity(Y_in, tempParticleInit));
		else
			InitX = w_lInit / (1. - w_lInit);
		os << "Initial particle moisture content: " << InitX * 1000 << " g/kg\nLiquid content: " << InitX / (1. + InitX) * 100 << " %";
		ShowInfo(os.str());
		os.str("");
		m_holdup->SetPhaseMass(_time, EPhase::LIQUID, InitX * M_tot);
		phiInit = InitX * M_tot / rhoWater / Delta_f / A_P;
	}
	//m_model.m_iPhi2 = m_model.AddDAEVariable(true, phiInit, 0, 1.0); // Degree of wetness of the particles
	//InitX = 0;

	//m_model.m_iA_p = m_model.AddDAEVariable(true, A_P, 0,1.0); // Total surface area of particles [m^2]
	//m_model.m_iMFlowVapor2 = m_model.AddDAEVariable(false, 0, 0, 0.0); // Vapor flow [kg/s]
	m_model.ClearVariables();

	std::vector<moistureContent> Y_inInit(N_total,Y_in); // Vector to initilize Y variable for model for all height layers
	//m_model.m_iY_g = m_model.AddDAEVariables(true, Y_inInit,0,1.0); // Moisture content of gas layers in [kg/kg]

	std::vector<temperature> temperatureGasInit(N_total, m_holdup->GetTemperature(_time));
	//m_model.m_iTempGas = m_model.AddDAEVariables(true, temperatureGasInit, 0,1.0); // Temperature of gas layers in [K]

	temperature InitTemperatureLiquid;
	if (particlesGlobal)
	{
		if (m_holdup->GetTemperature(_time) > m_inLiquidStream->GetTemperature(_time) && phiInit == 0)
			InitTemperatureLiquid = m_inLiquidStream->GetTemperature(_time);
		else
			InitTemperatureLiquid = m_holdup->GetTemperature(_time);
		m_model.m_iTempFilm = m_model.AddDAEVariable(true, InitTemperatureLiquid, 0, 1.0); // Temperature of liquid film in [K]
		m_model.m_iTempParticle = m_model.AddDAEVariable(true, tempParticleInit, 0, 1.0); // Temperature of particles in [K]
	}

	const pressure dP_sat = GetCompoundProperty(compoundKeys[indicesOfVaporOfPhaseChangingCompound.first], ECompoundTPProperties::VAPOR_PRESSURE, m_holdup->GetTemperature(_time), m_holdup->GetPressure(_time)); // Saturation partial pressure of phase changing compound at particle temperature [Pa]					// Equilibrium partial pressure [Pa]
	const double ratioMM = molarMassPhaseChangingLiquid / molarMassGas;
	//const moistureContent Y_sat = ratioMM * dP_sat / (m_holdup->GetPressure(_time) - (1 - ratioMM) * dP_sat);
	/*const double phiIn = Y_in / Y_sat;
	if (phi_eq < phiIn) {
		phi_eq = phiIn;
		RaiseWarning("Phi for exaust gas was raised to input levels.");
	}*/

	// State variables for Dosta model varaiables
	AddStateVariable("Avg. temperature gas", m_holdup->GetTemperature(_time)); // [°C]
	if (particlesGlobal)
	{
		AddStateVariable("Temperature film", 0); // [°C]
		AddStateVariable("Temperature particle", m_holdup->GetTemperature(_time)); // [°C]
		AddStateVariable("Moisture precentage", m_holdup->GetPhaseMass(_time, EPhase::LIQUID) / (m_holdup->GetPhaseMass(_time, EPhase::LIQUID) + m_holdup->GetPhaseMass(_time, EPhase::SOLID)));
		AddStateVariable("RH particles", GetRelativeHumidity(Y_in, m_holdup->GetTemperature(_time), m_holdup->GetPressure(_time)));
	}
	AddStateVariable("Avg. Y", Y_in); // [kg/kg]
	AddStateVariable("Temperature chamber", m_holdup->GetTemperature(_time));
	for (int i = 0; i < Y_inInit.size(); i++)
	{
		AddStateVariable("Temperature gas height layer " + std::to_string(i), m_holdup->GetTemperature(_time)); // [°C]
		AddStateVariable("Y height layer " + std::to_string(i), Y_in); // [kg/kg]
	}

	// Check mass and energy balance in DEBUG mode
	if (debugToggle)
	{
		AddStateVariable("Mass balance", 0);
		AddStateVariable("Energy balance", 0);
		for (int i = 0; i < 10; i++)
			AddStateVariable("debug" + std::to_string(i), 0);
	} 

	if (particlesGlobal)
		if (!CheckForSmallBiot(_time))
			RaiseWarning("Biot number: "+std::to_string(CalcBiotNumber(_time)) + ">" + std::to_string(SmallBiotNumber) + "\nAssumption of uniform temperature distribution throughout material is invalid.");
	/*if (DiffCoefWarning) {
		RaiseWarning("Critical eqData.temperatures and critical specific molar volumes used, are hard coded for air and water.");
		DiffCoefWarning = false;
	}*/

	m_model.InitCounterVaraibles();
	//ShowInfo(std::to_string(m_model.CalculateBetaPA(_time, m_holdup->GetTemperature(_time), m_holdup->GetTemperature(_time)) / m_model.CalculateBetaPA(_time, m_holdup->GetTemperature(_time), m_holdup->GetTemperature(_time) + 4)*100));
	debugHasBeenShown = false;
	//DiffCoefWarning = true;

	//m_model.unit = static_cast<CDryerBatch*>(this);
	if (particlesGlobal)
	{
		this->EnergyLiquidPhaseOld = this->m_holdup->GetPhaseMass(_time, EPhase::LIQUID) * this->C_PWaterLiquid * InitTemperatureLiquid;
		this->EnergySolidPhaseOld = this->m_holdup->GetPhaseMass(_time, EPhase::SOLID) * this->C_PParticle * tempParticleInit;
	}
	this->EnergyGasPhaseOld = this->m_holdup->GetPhaseMass(_time, EPhase::GAS) * ((this->C_PGas + this->C_PWaterVapor * Y_in) * m_holdup->GetTemperature(_time) + this->Delta_h0 * Y_in);

	//GetStringParameterValue("Path")

	/// Set tolerances to model ///
	const auto rtol = GetConstRealParameterValue("Relative tolerance");
	const auto atolT = GetConstRealParameterValue("Absolute tolerance T");
	const auto atolY = GetConstRealParameterValue("Absolute tolerance Y");

	std::vector<double> absolutTolerances;
	for (int i = 0; i < Y_inInit.size(); i++)
		absolutTolerances.push_back(atolY != 0.0 ? atolY : GetAbsTolerance()); // 0.0001
	for (int i = 0; i < temperatureGasInit.size() + (particlesGlobal ? 2 : 0); i++)
		absolutTolerances.push_back(atolT != 0.0 ? atolT : GetAbsTolerance()); // 0.01
	size_t variables = m_model.GetVariablesNumber();
	m_model.SetTolerance(rtol != 0.0 ? rtol : GetRelTolerance(), rtol != 0.0 ? rtol : GetRelTolerance());
	m_model.SetTolerance(rtol != 0.0 ? rtol : GetRelTolerance(), absolutTolerances); // 0.01
	m_solver.SetMaxStep(1);

	/// Set model to a solver ///
	if (!m_solver.SetModel(&m_model))
		RaiseError(m_solver.GetError());

	if (debugToggle)
	{
		os << "Phase masses:\nsolid: " << m_holdup->GetPhaseMass(_time, EPhase::SOLID);
		os << "\nliquid: " << m_holdup->GetPhaseMass(_time, EPhase::LIQUID);
		os << "\ngas/vapor: " << m_holdup->GetPhaseMass(_time, EPhase::GAS) << "\n";
		ShowInfo(os.str());
		os.str("");
	}


	os << "\nInitialization completed.\n";
	ShowInfo(os.str());

	//Testing();
}

void CDryerBatch::SaveState()
{
	m_solver.SaveState();
}

void CDryerBatch::LoadState()
{
	m_solver.LoadState();
}

void CDryerBatch::Simulate(double _timeBeg, double _timeEnd)
{
	if (!m_solver.Calculate(_timeBeg, _timeEnd))
		RaiseError(m_solver.GetError());
}

double CDryerBatch::CalculateSectionVolume(size_t section)
{
	double sectionVolume = 0;
	if (chamber.at(section).shape == EShape::CYLINDRICAL)
		sectionVolume = chamber.at(section).height * MATH_PI / 3 * (
			pow(chamber.at(section).dimensionsInternal.at(0).first / 2, 2)
			+ (chamber.at(section).dimensionsInternal.at(0).first * chamber.at(section).dimensionsInternal.at(0).second) / 4
			+ pow(chamber.at(section).dimensionsInternal.at(0).second / 2, 2)
			);
	else if (chamber.at(section).shape == EShape::RECTANGULAR)
	{
		std::stringstream os;
		os << section << " has an rectengular shape.\nWhich volume has not been implemented.";
		RaiseError(os.str());
	}
	else
	{
		std::stringstream os;
		os << section << " has an undefind shape.";
		RaiseError(os.str());
	}
		
	return sectionVolume;
}

// Functions to setup and check the inputs of for the calculation
void CDryerBatch::PullCompoundDataFromDatabase(double _time)
{
	/// Get general data for compounds
	this->compoundKeys = this->GetAllCompounds();
	std::vector<std::pair< EPhase, int>> compoundsKeyIndexPhasePartnerIndex(compoundKeys.size(), std::make_pair(EPhase::UNDEFINED, -2));
	std::vector<density> CompoundDensities(compoundKeys.size()); // Compound densities
	std::vector<heatCapacity> CompoundHeatCapacities(compoundKeys.size()); // Compound heat capacities
	std::vector<thermalConductivity> CompoundThermalConductivities(compoundKeys.size()); // Compound thermal conductivities
	std::vector<dynamicViscosity> CompoundDynViscosities(compoundKeys.size()); // Compound dynamic viscosities
	std::vector<molarMass> CompoundMolarMasses(compoundKeys.size()); // Compound molar masses
	//std::vector<double> CompoundCriticaleqData.temperatures(compoundKeys.size()); // Compound critical eqData.temperatures
	for (int i = 0; i < compoundKeys.size(); i++)
	{
		CompoundDensities[i] = GetCompoundProperty(compoundKeys[i], ECompoundTPProperties::DENSITY);
		CompoundHeatCapacities[i] = GetCompoundProperty(compoundKeys[i], ECompoundTPProperties::HEAT_CAPACITY_CP);
		CompoundThermalConductivities[i] = GetCompoundProperty(compoundKeys[i], ECompoundTPProperties::THERMAL_CONDUCTIVITY);
		CompoundDynViscosities[i] = GetCompoundProperty(compoundKeys[i], ECompoundTPProperties::VISCOSITY);
		CompoundMolarMasses[i] = GetCompoundProperty(compoundKeys[i], ECompoundConstProperties::MOLAR_MASS);
		//CompoundCriticaleqData.temperatures[i] = GetCompoundProperty(compoundKeys[i], ECompoundConstProperties::CRITICAL_TEMPERATURE);
	}

	// Properties of solids
	if (m_holdup->GetPhaseFraction(_time, EPhase::SOLID) == 0)
		RaiseWarning("Dryer holdup contains no solid.");

	std::vector<massFraction> SolidCompoundsDistribution = m_holdup->GetPhase(EPhase::SOLID)->GetCompoundsDistribution(_time);

	density tempRhoParticle = 0;
	heatCapacity tempHeatCapacityParticle = 0;
	thermalConductivity tempThermalConductivityParticle = 0;
	for (int i = 0; i < compoundKeys.size(); i++)
		if (SolidCompoundsDistribution[i] > 0)
		{
			compoundsKeyIndexPhasePartnerIndex[i] = std::make_pair(EPhase::SOLID, -1);
			tempRhoParticle += SolidCompoundsDistribution[i] * CompoundDensities[i];
			tempHeatCapacityParticle += SolidCompoundsDistribution[i] * CompoundHeatCapacities[i];
			tempThermalConductivityParticle += SolidCompoundsDistribution[i] * CompoundThermalConductivities[i];
		}
	rhoParticle = tempRhoParticle;
	C_PParticle = tempHeatCapacityParticle;
	lambdaParticle = tempThermalConductivityParticle;
	// Particle
	beta_GP = this->GetConstRealParameterValue("beta_GP");

	/// Properties of liquids
	if (m_inLiquidStream->GetPhaseFraction(_time, EPhase::LIQUID) == 0)
		RaiseWarning("Port: InletSuspension contains no liquid.");

	std::vector<massFraction> LiqudCompoundsDistribution = m_inLiquidStream->GetPhase(EPhase::LIQUID)->GetCompoundsDistribution(_time);

	density tempRhoLiquid = 0;
	heatCapacity tempHeatCapacityLiquid = 0;
	thermalConductivity tempLambdaLiquid = 0;
	for (int i = 0; i < compoundKeys.size(); i++)
	{
		tempRhoLiquid = tempRhoLiquid + LiqudCompoundsDistribution[i] * CompoundDensities[i];
		tempHeatCapacityLiquid = tempHeatCapacityLiquid + LiqudCompoundsDistribution[i] * CompoundHeatCapacities[i];
		tempLambdaLiquid = tempLambdaLiquid + LiqudCompoundsDistribution[i] * CompoundThermalConductivities[i];
	}
	rhoWater = tempRhoLiquid;
	C_PWaterLiquid = tempHeatCapacityLiquid;
	lambdaWater = tempLambdaLiquid;

	// Vector to store indices of liquid (first) and vapor (second) form of compound
	std::vector<std::pair<int, int>> liquidVaporIndices;

	for (int i = 0; i < compoundKeys.size(); i++) // Liquid form
	{
		if (LiqudCompoundsDistribution[i] > 0 && compoundsKeyIndexPhasePartnerIndex[i].first == EPhase::UNDEFINED)
		{
			double heatOfVaporizationLiquid = GetCompoundProperty(compoundKeys[i], ECompoundConstProperties::HEAT_OF_VAPORIZATION_AT_NORMAL_BOILING_POINT);
			double heatOfVaporizationVapor;
			for (int j = 0; j < compoundKeys.size(); j++) // Vapor form
			{
				if (i != j)
				{
					heatOfVaporizationVapor = GetCompoundProperty(compoundKeys[j], ECompoundConstProperties::HEAT_OF_VAPORIZATION_AT_NORMAL_BOILING_POINT);
					if (heatOfVaporizationLiquid == heatOfVaporizationVapor)
					{
						liquidVaporIndices.push_back(std::make_pair(i, j));
					}
				}
			}
		}
	}
	if (liquidVaporIndices.size() == 0)
		RaiseError("No compounds for phase transition found. Please add vapor or liquid form of compound or check material database entries for heat of vaporiation.");
	else
	{
		for (int i = 0; i < liquidVaporIndices.size(); i++)
			ShowInfo("Compounds found for phase change: " + GetCompoundName(compoundKeys[liquidVaporIndices[i].first]) + "," + GetCompoundName(compoundKeys[liquidVaporIndices[i].second]));
		compoundsKeyIndexPhasePartnerIndex[liquidVaporIndices[0].first] = std::make_pair(EPhase::LIQUID, liquidVaporIndices[0].second);
		compoundsKeyIndexPhasePartnerIndex[liquidVaporIndices[0].second] = std::make_pair(EPhase::VAPOR, liquidVaporIndices[0].first);
		indicesOfVaporOfPhaseChangingCompound = liquidVaporIndices[0];

		for (int i = 0; i < compoundKeys.size(); i++) // Liquid form
			if (LiqudCompoundsDistribution[i] > 0 && compoundsKeyIndexPhasePartnerIndex[i].first == EPhase::UNDEFINED)
				compoundsKeyIndexPhasePartnerIndex[i] = std::make_pair(EPhase::LIQUID, -1);
	}
	if (liquidVaporIndices.size() > 1)
		RaiseWarning("Only the first pair will be considered for phase change calculations.");

	C_PWaterVapor = GetCompoundProperty(compoundKeys[liquidVaporIndices[0].second], ECompoundTPProperties::HEAT_CAPACITY_CP); // Water vapor heat capacity
	Delta_h0 = GetCompoundProperty(compoundKeys[liquidVaporIndices[0].first], ECompoundConstProperties::HEAT_OF_VAPORIZATION_AT_NORMAL_BOILING_POINT) / GetCompoundProperty(compoundKeys[liquidVaporIndices[0].first], MOLAR_MASS); // Latent heat (evaporation heat) at 0 degree
	molarMassPhaseChangingLiquid = CompoundMolarMasses[liquidVaporIndices[0].first];

	/// Propterties of vapors
	if (m_inGasStream->GetPhaseFraction(_time, EPhase::GAS) != 1)
		RaiseError("Port: InletFluidizationGas contains non gaseous phases.");
	
	std::vector<double> VaporCompoundsDistribution = m_inGasStream->GetPhase(EPhase::GAS)->GetCompoundsDistribution(_time);

	density tempRhoGas = 0;
	heatCapacity tempHeatCapacityGas = 0;
	thermalConductivity templambdaGas = 0;
	dynamicViscosity tempEtaGas = 0;
	molarMass tempMolarMassGas = 0;
	//temperature tempCriticaleqData.temperatures = 0;
	for (int i = 0; i < compoundKeys.size(); i++)
		if (VaporCompoundsDistribution[i] > 0 && compoundsKeyIndexPhasePartnerIndex[i].first == EPhase::UNDEFINED)
		{
			compoundsKeyIndexPhasePartnerIndex[i] = std::make_pair(EPhase::GAS, -1);
			tempRhoGas += VaporCompoundsDistribution[i] * CompoundDensities[i];
			tempHeatCapacityGas += VaporCompoundsDistribution[i] * CompoundHeatCapacities[i];
			templambdaGas += VaporCompoundsDistribution[i] * CompoundThermalConductivities[i];
			tempEtaGas += VaporCompoundsDistribution[i] * CompoundDynViscosities[i];
			tempMolarMassGas += VaporCompoundsDistribution[i] * CompoundMolarMasses[i];
			//tempCriticaleqData.temperatures += VaporCompoundsDistribution[i] * CompoundCriticaleqData.temperatures[i];
		}
	rhoGas = tempRhoGas;
	C_PGas = tempHeatCapacityGas;
	lambdaGas = templambdaGas;
	etaGas = tempEtaGas;
	molarMassGas = tempMolarMassGas;
	//T_critGas = tempCriticaleqData.temperatures;

	for (int i = 0; i < compoundsKeyIndexPhasePartnerIndex.size(); i++)
		if (compoundsKeyIndexPhasePartnerIndex[i].first == EPhase::UNDEFINED)
			RaiseError("Compound " + GetCompoundName(compoundKeys[i]) + "could not be fitted into solid, liquid or vapor/gas phase");

	CompoundsKeyIndexPhasePartnerIndex = compoundsKeyIndexPhasePartnerIndex;

	if (debugToggle){
		std::ostringstream  os;
		os << "Debug info\n"
			<< "\trhoGas: " << rhoGas << "\n"
			<< "\tetaGas: " << etaGas << "\n"
			<< "\tC_PGas: " << C_PGas << "\n"
			<< "\tlambdaGas: " << lambdaGas << "\n"
			<< "\tmolarMassGas: " << molarMassGas << "\n" << "\n"

			<< "\trhoWater: " << rhoWater << "\n"
			<< "\tC_PWaterLiquid: " << C_PWaterLiquid << "\n"
			<< "\tC_PWaterVapor: " << C_PWaterVapor << "\n"
			<< "\tDelta_h0: " << Delta_h0 << "\n"
			<< "\tlambdaWater: " << lambdaWater << "\n"
			<< "\tmolarMassPhaseChangingLiquid: " << molarMassPhaseChangingLiquid << "\n" << "\n"

			<< "\trhoParticle: " << rhoParticle << "\n"
			<< "\tbeta_GP: ";
		if (beta_GP < 0)
			os << "automatic calculation";
		else
			os << beta_GP;
		os  << "\n"
			<< "\tC_PParticle: " << C_PParticle << "\n";
		ShowInfo(os.str());
		os.str("");
	} // DEBUG
}

void CDryerBatch::CheckHeightDiscretizationLayers(double _time)
{
	N_particle = this->GetConstIntParameterValue("N_el");
	std::vector<double> sizeGrid = GetNumericGrid(DISTR_SIZE);
	double maxDistributionSize = sizeGrid.back();
	if (chamber.at(0).height / N_particle < maxDistributionSize)
	{
		this->RaiseWarning("Height of discretization layer smaller than biggest particles.");
		int64_t tempN_el = N_particle -  1;
		while (chamber.at(0).height / tempN_el < maxDistributionSize)
		{
			tempN_el--;
			if (tempN_el == 0)
				break;
		}
		N_particle = tempN_el;
		if (N_particle < 1)
			this->RaiseError("Size of larges particles " + std::to_string(maxDistributionSize) + " m is larger than total bed height " + std::to_string(chamber.at(0).height) + " m.");
		this->RaiseWarning("Adjusted number of discretization layers to: " + std::to_string(N_particle));
	}
	for (auto& element : chamber)
		if (element.layers < 1)
			element.layers = 1;
	chamber.at(0).layers = N_particle;
	double sectionsFilledWithBed = 0;
	if (particlesGlobal)
	{
		std::vector<double> inGasStreamTimePoints = m_inGasStream->GetAllTimePoints();
		double maxInGasFlow = 0;
		heighestFlowTimepoint = 0;
		for (auto& element : inGasStreamTimePoints)
			if (m_inGasStream->GetMassFlow(element) > maxInGasFlow)
			{
				maxInGasFlow = m_inGasStream->GetMassFlow(element);
				heighestFlowTimepoint = element;
			}
		sectionsFilledWithBed = DetermineSectionsFilledWithBed(heighestFlowTimepoint, m_holdup->GetTemperature(heighestFlowTimepoint));
		// sectionsFilledWithBed = CalculateBedHeightOrDetermineSectionsFilledWithBed(heighestFlowTimepoint, std::max(m_holdup->GetTemperature(0), m_inGasStream->GetTemperature(heighestFlowTimepoint)), false);
		for (int i = 1; i < std::ceil(sectionsFilledWithBed); i++)
		{
			chamber.at(i).layers = std::ceil(chamber.at(0).layers * chamber.at(i).height / chamber.at(0).height);
			N_particle += chamber.at(i).layers;
		}
	}
	N_total = 0;
	for (auto& section : chamber)
		N_total += section.layers;
}

bool CDryerBatch::CheckForSmallBiot(double _time) const // Needed for uniform temperature distribution (Diss Soeren, chapter 2.2.2, page 27)
{
	const double Bi = CalcBiotNumber(_time);
	if (Bi < SmallBiotNumber)
		return true;
	else
		return false;
}

//bool CDryerBatch::LoadDimensions(std::string path)
//{
//	return false;
//}

/// Functions to calculate/return moisture related properties ///
moistureContent CDryerBatch::GetGasSaturationMoistureContent(temperature temperatureGas, pressure pressureGas)
{
	if (Y_sat == 0) // Moisture content calculation under assumption of ideal gas law
	{
		pressure P_sat = GetCompoundProperty(compoundKeys[indicesOfVaporOfPhaseChangingCompound.first], ECompoundTPProperties::VAPOR_PRESSURE, temperatureGas, pressureGas); // Saturation partial pressure of phase changing compound at particle temperature [Pa]					// Equilibrium partial pressure [Pa]
		const double ratioMM = molarMassPhaseChangingLiquid / molarMassGas;
		bool pressureReductionWarning = false;
		int i = 0;
		while (pressureGas - P_sat < 0)
		{
			P_sat = P_sat * 0.999;
			i++;
		}
		if (pressureReductionWarning)
			RaiseWarning("Pressure difference was non positive.\n Reducing saturation partial pressure by " + std::to_string(0.1 * i) + " % .\n Temperature : " + std::to_string(temperatureGas));
		const moistureContent Y_sat = ratioMM * P_sat / (pressureGas - P_sat); // Moisture content in kg per kg dry gas
		return Y_sat;
	}
	else // same as user input
	{
		return Y_sat;
	}
}

moistureContent CDryerBatch::CalculateGasEquilibriumMoistureContent(temperature temperatureParticle, pressure pressureGas, double RH) const
{
	// Moisture content calculation
	// Assumption: ideal gas law
	if (RH == 0)
		return 0;

	pressure P_sat = GetCompoundProperty(compoundKeys[indicesOfVaporOfPhaseChangingCompound.first], ECompoundTPProperties::VAPOR_PRESSURE, temperatureParticle, pressureGas); // Saturation partial pressure of phase changing compound at particle temperature [Pa]					// Equilibrium partial pressure [Pa]

	const double ratioMM = molarMassPhaseChangingLiquid / molarMassGas;

	const pressure P_eq = P_sat * RH;

	const moistureContent Y_eq = ratioMM * P_eq / (pressureGas - P_eq); // Moisture content in kg per kg dry gas
	return Y_eq;
}

double CDryerBatch::GetRelativeHumidity(moistureContent Y, temperature temperature, pressure pressure) //const
{
	/*
	// To Do: make variables fit
	const double M_air = 1.;				// Reference mass of air [kg]
	const double M_wg = M_air * _y;	// Mass of gaseous water [kg]

	const double n_air = M_air / molarMassGas;	// Amount of substance of air [mol]
	const double n_wg = M_wg / molarMassPhaseChangingLiquid;		// Amount of substance of gaseous water [mol]

	const double P_wg = n_wg / (n_wg + n_air) * _p; // Partial pressure of gaseous water [Pa]
	const double P_sat = GetCompoundProperty(compoundKeys[indicesOfVaporOfPhaseChangingCompound.first], ECompoundTPProperties::VAPOR_PRESSURE, _t, _p); // Saturated vapor pressure of water [Pa]

	const double phi = P_wg / P_sat; // Relative humidity [-]
	*/
	const moistureContent y_sat = GetGasSaturationMoistureContent(temperature, pressure);
	double phi = Y / y_sat;
	if (phi > 1)
		phi = 1;
	if (phi < 0)
		phi = 0;
	return phi;
}

moistureContent CDryerBatch::CalcuateSolidEquilibriumMoistureContent(double _time, temperature temperature, double RH)
{
	std::vector<moistureContent> EquilibriumMoistures(compoundKeys.size());
	moistureContent particleEquilibriumMoistureContent = 0;

	std::vector<massFraction> SolidCompoundsDistribution = m_holdup->GetPhase(EPhase::SOLID)->GetCompoundsDistribution(_time);

	for (int i = 0; i < compoundKeys.size(); i++)
		if (compoundKeys.at(i) == eqData.compoundKey)
			EquilibriumMoistures[i] = GetParticleEquilibriumMoistureContent(temperature, RH);
		else
			EquilibriumMoistures[i] = GetCompoundProperty(compoundKeys[i], ECompoundTPProperties::EQUILIBRIUM_MOISTURE_CONTENT, temperature, RH);
		// Equilibrium moisture content depends on temperatur and humidity
	for (int i = 0; i < compoundKeys.size(); i++)
		if (EquilibriumMoistures.at(i) < 0)
		{
			std::stringstream os;
			os << "Moisture content of " << GetCompoundName(compoundKeys.at(i)) << "at " << _time << "s returned a negativ moisture content for " << temperature << "K at "<< RH*100 <<"%.";
			RaiseError(os.str());
			os.str() = "";
		}

	for (int i = 0; i < compoundKeys.size(); i++)
		particleEquilibriumMoistureContent += SolidCompoundsDistribution[i] * EquilibriumMoistures[i];

	particleEquilibriumMoistureContent *= moistureScaler;
	return std::max(particleEquilibriumMoistureContent, minMoistureContent);
}

//double CDryerBatch::CalcuateEquilibriumRelativeHumidity(double _time, temperature temperature, double X) const
//{
//	return 0;
//}

double CDryerBatch::CalculateNormalizedDryingCurve(moistureContent X, moistureContent Xeq)
{
	// Neglect case for which X is smaller than X_eq, particle would take moisture from gas
	if (X <= Xeq)
	{
		return 0;
	}

	switch (dryingCurveSetting) // calculate evaporation rate if REA, or relative drying rate if NCDC
	{
		case 0: // REA: calculate evaporation rate in form of inverted REA
			return 1. - REA(X - Xeq);
			break;
		case 1: // NCDC: calculate relative drying rate 
			if (X < X_cr) // 2nd drying period
			{
				double deta = (X - Xeq) / (X_cr - Xeq); // Normalized moisture content [-]
				double dnu = k_dc * deta / (1. + deta * (k_dc - 1.)); // Normalized Drying Curve: model of van Meel (1958)
				/*if (dnu != dnu)
					RaiseError("Normalized Drying Curve based model of van Meel resulted in NaN");*/
				// If the moisture content X is larger than critical moisture content X_cr, nu is limited to 1 (first drying period)
				return dnu;
			}
			else // 1st drying period
				return 1;
			break;
		case 2: // No curve: only 1st drying period
			return 1;
			break;
	}
}

/// Supporting functions ///
double CDryerBatch::GetAvgConstCompoundProperty(double _time, EPhase phase, ECompoundConstProperties  property) const
{
	std::vector<double> CompoundProperties(compoundKeys.size());
	for (int i = 0; i < compoundKeys.size(); i++)
		CompoundProperties[i] = GetCompoundProperty(compoundKeys[i], property);

	std::vector<double> PhaseCompoundsDistribution = m_holdup->GetPhase(phase)->GetCompoundsDistribution(_time);

	double avgProperty = 0;
	for (int i = 0; i < compoundKeys.size(); i++)
		avgProperty += PhaseCompoundsDistribution[i] * CompoundProperties[i];
	return avgProperty;
}

double CDryerBatch::GetAvgTPCompoundProperty(double _time, EPhase phase, ECompoundTPProperties  property, double temperature, double pressure) const
{
	std::vector<double> CompoundProperties(compoundKeys.size());
	for (int i = 0; i < compoundKeys.size(); i++)
		CompoundProperties[i] = GetCompoundProperty(compoundKeys[i], property, temperature, pressure);

	std::vector<double> PhaseCompoundsDistribution = m_holdup->GetPhase(phase)->GetCompoundsDistribution(_time);

	double avgProperty = 0;
	for (int i = 0; i < compoundKeys.size(); i++)
		avgProperty += PhaseCompoundsDistribution[i] * CompoundProperties[i];
	return avgProperty;
}

/// Heat transfer ///
double CDryerBatch::CalckAc(double alphaIn, double alphaOut, double L, std::vector<double> d/* Inner to outer diameter*/, std::vector<double> lambda) const
{
	double wall = 0;
	if (lambda.size() == d.size() - 1)
		for (int i = 0; i < lambda.size(); i++)
			wall += log(d[i + 1] / d[i]) / lambda[i];
	double aIn = 0;
	if ((alphaIn * d.front()) != 0)
		aIn = 1. / (alphaIn * d.front());
	double aOut = 0;
	if ((alphaOut * d.back()) != 0)
		aIn = 1. / (alphaOut * d.back());
	double sum = aIn + 1. / 2 * wall + aOut;
	const double kA = 1. / (1. / MATH_PI / L * sum);
	return kA;
	/*
		lambda of 1.4404 steel = 15 W/m/K at 20°C
							T°C rho	al	cp	lam
		S31603				-100	15.0		 [73, 217]
		Wst-Nr. 1.4404		0		16.0 466
		X2CrNiMoi17-12-2	20 7956	16.2 470 12.7
		ASTM/AISI 316L		100		16.7 486 13.8
							200		17.1 501 15.5
							400		18.1 518 18.6
							600		18.8 539 21.7
							800		19.3 557 24.8

		https://doi.org/10.1007/978-3-662-57572-7_3
	*/
}

double CDryerBatch::CalckAp(double alphaIn, double alphaOut, std::vector<double> A, std::vector<double> delta, std::vector<double> lambda) const
{
	double wall = 0;
	if (lambda.size() == A.size() - 2 && lambda.size() == delta.size())
		for (int i = 0; i < lambda.size(); i++)
			wall += delta[i] / (lambda[i] * A[i + 1]);
	double aIn = 0;
	if ((alphaIn * A.front()) != 0)
		aIn = 1. / (alphaIn * A.front());
	double aOut = 0;
	if ((alphaOut * A.back()) != 0)
		aIn = 1. / (alphaOut * A.back());
	double sum = aIn + wall + aOut;
	const double kA = 1. / sum;
	return kA;
}

double CDryerBatch::CalcAlphaOutside(double _time, const double h, const double D, const double Ts, EShape shape) const
{
	// VDI Wärmeatlas F2
	const double Tstar = 0.5 * (Ts + T_inf); // Temperature for properties
	const double C_PGas = GetAvgTPCompoundProperty(_time, EPhase::GAS, ECompoundTPProperties::HEAT_CAPACITY_CP, Tstar);
	const double etaGas = GetAvgTPCompoundProperty(_time, EPhase::GAS, ECompoundTPProperties::VISCOSITY, Tstar);
	const double lambdaGas = GetAvgTPCompoundProperty(_time, EPhase::GAS, ECompoundTPProperties::THERMAL_CONDUCTIVITY, Tstar);
	const double rhoGas = GetAvgTPCompoundProperty(_time, EPhase::GAS, ECompoundTPProperties::DENSITY, Tstar);

	const double nyGas = etaGas / rhoGas;
	const double Pr = C_PGas * etaGas / lambdaGas;
	const double beta = 1 / (Tstar);
	const double deltaT = abs(Ts - T_inf);
	const double Gr = beta * STANDARD_ACCELERATION_OF_GRAVITY * deltaT * pow(h, 3) / pow(nyGas, 2);
	const double Ra = Gr * Pr;

	// ToDo - Add Ra under angle
	// ToDo - Implement surface line instead of height

	const double f1 = pow(1. + pow(0.492 / Pr, 9. / 16), -16. / 9);
	const double NuP = pow(0.825 + 0.387 * pow(Ra * f1, 1. / 6), 2);
	const double Nu = shape==EShape::CYLINDRICAL? NuP + 0.435 * h / D:NuP;
	const double alpha = Nu * lambdaGas / h;
	// Expected range 2.5-25 https://www.sciencedirect.com/topics/engineering/convection-heat-transfer-coefficient
	return alpha;
}

double CDryerBatch::CalculateAlpha_PW(double _t_p, double _t_g, double _p, double _time) const
{
	/// Get unit parameters
	double eps_mf = eps_0;				// Bed porosity at minimum fluidization, equal to INITIAL bed porosity [-]
	double eps_fb = CalculateBedPorosity(_time);						// Bed porosity at operating conditions [-]
	double dp_p = CalculateHoldupSauter(_time);							// Sauter diameter of particles [m]

	// Calculate modified free path of the gas molecule
	double C_A = 2.8;
	double gamma = 1. / (pow(10., 0.6 - (1000. / _t_g + 1.) / C_A) + 1); // gamma, A.3.85

	double l = 2. * (2. / gamma - 1.) * sqrt(2. * MATH_PI * MOLAR_GAS_CONSTANT * _t_g / molarMassGas) * lambdaGas / _p / (2. * C_PGas - MOLAR_GAS_CONSTANT / molarMassGas); //modified free path of the gas molecule, A.3.84

	//Calculate parameter Z and N
	double C_k = 2.6; //C_k, A.3.82
	double Nu_pw_max = 4. * ((1. + 2. * l / dp_p) * log(1. + dp_p / 2. / l) - 1.); //Nu_PW,max, A.3.83
	double Z = 1. / 6. * rhoParticle * C_PParticle / lambdaGas * sqrt(STANDARD_ACCELERATION_OF_GRAVITY * pow(dp_p, 3.) * (eps_fb - eps_mf) / 5. / (1 - eps_mf) / (1 - eps_fb)); //Parameter Z, A.3.80
	double N = Nu_pw_max / C_k / Z; //Parameter N, A.3.81

	//Calculate Nusselt Number and heat transfer coefficient
	double Nu_pw = (1. - eps_fb) * Z * (1. - exp(-N));
	double alpha_pw = Nu_pw * lambdaGas / dp_p;

	return alpha_pw;
}

double CDryerBatch::CalculateAlpha_GW(double _time) const
{
	const double Pr = CalculatePrandtl(_time);
	const double Ar = CalculateArchimedes(_time);
	const double NuG = 0.009 * pow(Pr, 1. / 3) * sqrt(Ar);
	double d = CalculateHoldupSauter(_time);
	const double alphaGW = NuG * lambdaGas / d;
	return alphaGW;
}

double CDryerBatch::CalculateReynolds(double _time, size_t section) const
{
	double d = (chamber.at(section).dimensionsInternal.at(0).first + chamber.at(section).dimensionsInternal.at(0).second) / 2;
	return CalculateGasVel(_time,section)* d* rhoGas / etaGas;
}

double CDryerBatch::CalculateAlpha_GW(double _time, size_t section) const
{
	const double Re = CalculateReynolds(_time, section);
	const double Pr = CalculatePrandtl(_time);
	double di = (chamber.at(section).dimensionsInternal.at(0).first + chamber.at(section).dimensionsInternal.at(0).second) / 2;
	double Nu_m2 = 1.616 * std::cbrt(Re * Pr * di / chamber.at(section).height);
	double Nu2300 = std::cbrt(pow(3.66, 3) + pow(0.7, 3) + pow(Nu_m2 - 0.7, 3));
	double zeta = pow(1.8 * log10(Re) - 1.5, -2) / 8;
	double Nu1e4 = zeta * Re * Pr * (1 + pow(di / chamber.at(section).height, 2. / 3)) / (1. + 12.7 * sqrt(zeta) * (pow(Pr, 2. / 3) - 1));
	double Nu = 0;
	if (Re < 2300)
		Nu = Nu2300;
	else if (Re > 1e4)
		Nu = Nu1e4;
	else
	{
		double gamma = (Re - 2300) / (1e4 - 2300);
		Nu = (1 - gamma) * Nu2300 + gamma * Nu1e4;
	}
	const double alphaGW = Nu * lambdaGas / di;
	return alphaGW;
}

double CDryerBatch::CalcAlphaParticleGas(double _time) const
{
	const double Pr = CalculatePrandtl(_time);
	const double Re = CalculateReynolds(_time);
	const double eps = CalculateBedPorosity(_time); // Bed porosity
	const double Nu_lam = CalculateNusseltSherwoodLam(Re / eps, Pr);
	const double Nu_turb = CalculateNusseltSherwoodTurb(Re / eps, Pr);
	const double Nu_SP = CalculateNusseltSherwoodSingleParticle(Nu_lam, Nu_turb);
	const double d32 = CalculateHoldupSauter(_time); // Sauter diameter of holdup
	const double Nu_bed = CalculateNusseltSherwoodBed(eps, Nu_SP);
	const double alpha = Nu_bed * lambdaWater / d32; // Heat transfer coefficient
	return alpha;
}

/// Mass transfer ///
double CDryerBatch::CalculateDiffusionCoefficient(double _time, double avgGasTemperature, double filmTemperature, double pressure) const
{
	switch (DiffCoeff)
	{
		case 0: // Dosta (2010): https://doi.org/10.1016/j.powtec.2010.07.018
			return (23e-5) * pow(avgGasTemperature / T_ref, 1.81); 
			break;
	
		case 1: // Tsotsas
			return 2.252 / m_holdup->GetPressure(_time) * pow(avgGasTemperature / T_ref, 1.81);
			break;

		case 2: // correlation Poos & Varju (2020): https://doi.org/10.1016/j.ijheatmasstransfer.2020.119500
			const double T_critGas = GetAvgConstCompoundProperty(_time, EPhase::GAS, CRITICAL_TEMPERATURE); // critical temperature air
			const double T_critPcL = GetCompoundProperty(compoundKeys[indicesOfVaporOfPhaseChangingCompound.first], CRITICAL_TEMPERATURE); // critical temperature of liquid on particle (water)
			const double V_critGas = GetAvgConstCompoundProperty(_time, EPhase::GAS, MOLAR_MASS) / GetAvgConstCompoundProperty(_time, EPhase::GAS, CONST_PROP_USER_DEFINED_01) * 1e6;
			const double V_critPcL = GetCompoundProperty(compoundKeys[indicesOfVaporOfPhaseChangingCompound.first], MOLAR_MASS) / GetCompoundProperty(compoundKeys[indicesOfVaporOfPhaseChangingCompound.first], CONST_PROP_USER_DEFINED_01) * 1e6;
			return (1.498e-6 * pow(filmTemperature, 1.81) * pow(1 / (molarMassGas * 1000) + 1 / (molarMassPhaseChangingLiquid * 1000), 0.5))
				/ ((pressure / STANDARD_CONDITION_P) * pow(T_critGas * T_critPcL, 0.1405) * pow(pow(V_critGas, 0.4) + pow(V_critPcL, 0.4), 2));
			break;
	}
}

double CDryerBatch::CalculateBetaPA(double _time, double avgGasTemperature, double filmTemperature, double D) const
{
	if (avgGasTemperature != avgGasTemperature || filmTemperature != filmTemperature || avgGasTemperature < 0 || filmTemperature < 0)
		return -1;
	if (beta_GP >= 0)
		return beta_GP;

	//// VDI M5.2
	//const double D_Dosta = (23e-5) * pow(avgGasTemperature / T_ref, 1.81); // Diffustion coefficient water vapor in air
	//// Dosta 2010 A.3
	//const double D_Poos = CalculateDiffusionCoefficient(_time,filmTemperature);
	//// https://doi.org/10.1016/j.ijheatmasstransfer.2020.119500
	//const double D_Tsotsas = 2.252 / m_holdup->GetPressure(_time) * pow(avgGasTemperature / T_ref, 1.81);
	//double D = 1;
	//switch (DiffCoeff) {
	//case 1:
	//	D = D_Dosta;
	//	break;
	//case 2:
	//	D = D_Tsotsas;
	//	break;
	//default:
	//	D = D_Poos;
	//	break;
	//}
	const double Sc = etaGas / (D * rhoGas);
	const double Re = CalculateReynolds(_time);
	const double eps = CalculateBedPorosity(_time); // Bed porosity
	const double Sh_lam = CalculateNusseltSherwoodLam(Re / eps, Sc);
	const double Sh_turb = CalculateNusseltSherwoodTurb(Re / eps, Sc);
	const double Sh = CalculateNusseltSherwoodSingleParticle(Sh_lam, Sh_turb);
	const double d32 = CalculateHoldupSauter(_time); // Sauter diameter of holdup
	const double beta_FA = CalculateNusseltSherwoodBed(eps, Sh) * D / d32; // Mass transfer coefficient // ToDO - Look at this!
	// Dosta 2010 A.3

	return beta_FA;
}

/// Bed and particle properties ///
double CDryerBatch::CalculateHoldupSauter(double _time) const
{
	std::vector<double> sizeGrid = GetNumericGrid(DISTR_SIZE);
	std::vector<double> q3_holdup = m_holdup->GetPSD(_time, PSD_q3);
	double d32 = GetSauterDiameter(sizeGrid, q3_holdup);
	return d32;
}

area CDryerBatch::CalculateParticleSurfaceArea(double _time) const
{
	std::vector<double> Grid = GetNumericGrid(DISTR_SIZE);
	const double M_tot = m_holdup->GetPhaseMass(_time, EPhase::SOLID);
	std::vector<double> q_3 = m_holdup->GetPSD(_time, PSD_q3, EPSDGridType::DIAMETER);
	const double A = GetSpecificSurface(Grid, q_3) / rhoParticle * M_tot;
	return A;
}

double CDryerBatch::CalculateBedHeight(double _time, double particleTemperature)
{
	const double massSolid = m_holdup->GetPhaseMass(_time, EPhase::SOLID);
	const double densitySolid = GetAvgTPCompoundProperty(_time, EPhase::SOLID, ECompoundTPProperties::DENSITY, particleTemperature, m_holdup->GetPressure(_time));
	const double volumeSolid =  massSolid / densitySolid;
	const double eps = CalculateBedPorosity(_time);
	const double volumeBed = volumeSolid / (1. - eps);
	double volumeOfFilledSections = 0;
	size_t numberOfFilledSection = 0;
	double heightOfFilledSection = 0;
	while (volumeBed > (CalculateSectionVolume(numberOfFilledSection)+ volumeOfFilledSections))
	{
		volumeOfFilledSections += CalculateSectionVolume(numberOfFilledSection);
		heightOfFilledSection += chamber.at(numberOfFilledSection).height;
		numberOfFilledSection++;
		if (numberOfFilledSection >= chamber.size())
			RaiseError("Bed expands out of chamber.");
	}

	double r = chamber.at(numberOfFilledSection).dimensionsInternal.at(0).first / 2;
	double R = chamber.at(numberOfFilledSection).dimensionsInternal.at(0).second / 2;
	double h = chamber.at(numberOfFilledSection).height;
	const double tanWallAngle = (R - r) / h;
	const double heightBed = heightOfFilledSection -std::cbrt(-MATH_PI * pow(r, 3) * pow(tanWallAngle, 3) - 3 * pow(tanWallAngle, 4) * volumeBed) / (std::cbrt(MATH_PI) * pow(tanWallAngle, 2))-r / tanWallAngle;
	return heightBed;
}

double CDryerBatch::CalculateBedHeightOrDetermineSectionsFilledWithBed(double _time, double particleTemperature, bool outputHeight)
{
	const double massSolid = m_holdup->GetPhaseMass(_time, EPhase::SOLID);
	const double densitySolid = GetAvgTPCompoundProperty(_time, EPhase::SOLID, ECompoundTPProperties::DENSITY, particleTemperature, m_holdup->GetPressure(_time));
	const double volumeSolid = massSolid / densitySolid;
	const double eps = CalculateBedPorosity(_time);
	const double volumeBed = volumeSolid / (1. - eps);
	double volumeOfFilledSections = 0;
	size_t numberOfFilledSection = 0;
	double heightOfFilledSection = 0;
	while (volumeBed > (CalculateSectionVolume(numberOfFilledSection) + volumeOfFilledSections))
	{
		volumeOfFilledSections += CalculateSectionVolume(numberOfFilledSection);
		heightOfFilledSection += chamber.at(numberOfFilledSection).height;
		numberOfFilledSection++;
		if (numberOfFilledSection >= chamber.size())
			RaiseError("Bed expands out of chamber.");
	}

	double r = chamber.at(numberOfFilledSection).dimensionsInternal.at(0).first / 2;
	double R = chamber.at(numberOfFilledSection).dimensionsInternal.at(0).second / 2;
	double h = chamber.at(numberOfFilledSection).height;
	const double tanWallAngle = (R - r) / h;
	const double bedHeight = -std::cbrt(-MATH_PI * pow(r, 3) * pow(tanWallAngle, 3) - 3 * pow(tanWallAngle, 4) * volumeBed) / (std::cbrt(MATH_PI) * pow(tanWallAngle, 2)) - r / tanWallAngle;
	const double totalBedHeight = heightOfFilledSection + bedHeight;
	if (outputHeight)
		return totalBedHeight;
	else
		return numberOfFilledSection + bedHeight / h;
}

double CDryerBatch::CalculateBedPorosity(double _time, bool homogeniusFluidization) const
{
	const double Ar = CalculateArchimedes(_time);
	const double ReL = 42.9 * (1.-eps_0)*(sqrt(1.+pow(eps_0,3)*Ar/(3214*pow(1.-eps_0,2)))-1);
	const double ReA = homogeniusFluidization ? 18*pow(sqrt(1.+sqrt(Ar) / 9) - 1, 2) : sqrt(4 * Ar / 3);
	const double Re = CalculateReynolds(_time);
	const double n = log(ReL / ReA) / log(eps_0);
	const double eps = pow(Re / ReA, 1 / n);
	return eps;
}

void CDryerBatch::SetupChamber()
{
	if (chamber.size() == 0)
	{
		chamberSection section;
		section.name = "main";
		section.shape = EShape::CYLINDRICAL;
		section.dimensionsInternal = { std::make_pair(0.180, 0.250) };
		section.height = { 0.337 };
		section.wallThicknesses = { 0.005 };
		section.thermalConductivities = { 15 };
		section.layers = 1;
		chamber.push_back(section);

		section.name = "adapter";
		section.dimensionsInternal = { std::make_pair(0.250, 0.500) };
		section.height = { 0.345 };
		chamber.push_back(section);

		section.name = "expander";
		section.dimensionsInternal = { std::make_pair(0.500, 0.500) };
		section.height = { 0.950 }; // 0.500 from Glatt, 0.95 measured
		chamber.push_back(section);
	}
}

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

double CDryerBatch::CalculateOverallHeatTransferCoefficientCylinder(size_t section, double alphaInternal, double alphaExternal, double heightUsage)
{
	if (chamber.at(section).thermalConductivities.size() != chamber.at(section).wallThicknesses.size() || heightUsage > 1 || heightUsage <= -1)
	{
		std::ostringstream  os;
		os << "CalculateOverallHeatTransferCoefficientCylinder has encountered an error in section: ";
		os << section;
		os << ".\nCheck if numer of wall thicknesses is equal to thermal conductivities.";
		os << "\nOr heightUsage was out of bounds.";
		RaiseError(os.str());
	}
		const double tanWallAngle = (chamber.at(section).dimensionsInternal.at(0).second - chamber.at(section).dimensionsInternal.at(0).first) / (2 * chamber.at(section).height);
	
	double h, R, r, R_Plus_r;

	if (heightUsage > 0)
	{
		h = heightUsage * chamber.at(section).height;
		r = chamber.at(section).dimensionsInternal.at(0).first / 2;
		R = h * tanWallAngle + r;
		R_Plus_r = (R + r);
	}
	else
	{
		h = (1-heightUsage) * chamber.at(section).height;
		R = chamber.at(section).dimensionsInternal.at(0).second / 2;
		r = R - h * tanWallAngle;
		R_Plus_r = (R + r);
	}

	double wall = 0;
	for (int i = 0; i < chamber.at(section).thermalConductivities.size(); i++)
	{
		double sumWallThicknesses = 0;
		for (int j = 0; j < chamber.at(section).wallThicknesses.size() - 1; j++)
			sumWallThicknesses += chamber.at(section).wallThicknesses.at(j);
		wall += log((R_Plus_r + 2 * (chamber.at(section).wallThicknesses.at(i) + sumWallThicknesses)) / (R_Plus_r + sumWallThicknesses)) / chamber.at(section).thermalConductivities.at(i);
	}
	double sumWallThicknesses = 0;
	for (int i = 0; i < chamber.at(section).wallThicknesses.size(); i++)
		sumWallThicknesses += chamber.at(section).wallThicknesses.at(i);
	double aIn = 0;
	if (alphaInternal * R_Plus_r != 0)
		aIn = 1. / (alphaInternal * R_Plus_r);
	double aOut = 0;
	if (alphaExternal * (R_Plus_r + sumWallThicknesses) != 0)
		aIn = 1. / (alphaExternal * (R_Plus_r + sumWallThicknesses));
	double sum = aIn + 1. / 2 * wall + aOut;

	const double m = sqrt(pow(R - r, 2) + pow(h, 2));
	const double kA = 1. / (sum / (MATH_PI * m));
	return kA;
}

double CDryerBatch::CalculateOverallHeatTransferCoefficientTopPlate(double alphaInternal, double alphaExternal)
{
	if (chamber.back().thermalConductivities.size() != chamber.back().wallThicknesses.size())
		bool error = true;
	double A = 0;
	if (chamber.back().shape == EShape::CYLINDRICAL)
		A = pow(chamber.back().dimensionsInternal.at(0).second,2) * MATH_PI/4;
	else
		A = chamber.back().dimensionsInternal.at(0).second * chamber.back().dimensionsInternal.at(1).second;

	double wall = 0;
	for (int i = 0; i < chamber.back().thermalConductivities.size(); i++)
		wall += chamber.back().wallThicknesses.at(i) / (chamber.back().thermalConductivities.at(i) * A);
	double aIn = 0;
	if ((alphaInternal * A) != 0)
		aIn = 1. / (alphaInternal * A);
	double aOut = 0;
	if ((alphaExternal * A) != 0)
		aIn = 1. / (alphaExternal * A);
	double sum = aIn + wall + aOut;
	if (sum == 0)
		bool error = true;
	const double kA = 1. / sum;
	return kA;
}

/// Fluiddynamic ///
double CDryerBatch::CalculateMinFluidizeVel(double _time) const
{
	double d32 = CalculateHoldupSauter(_time);
	double u_mf = this->u_mf;
	if (u_mf == 0)
	{
		double Ar_mf = CalculateArchimedes(_time);
		double Re_mf = 33.7 * (sqrt(1. + 3.6e-5 * Ar_mf) - 1);
		u_mf = Re_mf * etaGas / (d32 * rhoGas);
	}
	return u_mf;
}

double CDryerBatch::CalculateGasVel(double _time, size_t section) const
{
	double mFlow_gasIn = m_inGasStream->GetMassFlow(_time);
	double VFlow_gasIn = mFlow_gasIn / rhoGas;
	double area_bed = MATH_PI * pow(chamber.at(0).dimensionsInternal.at(0).first, 2) / 4;
	double u_gasIn = VFlow_gasIn / area_bed;
	double u_mf = CalculateMinFluidizeVel(_time);
	double u_gasHoldup = (u_gasIn - u_mf) / 3 + u_mf; // porous plate distributor
	if (section > 0)
		return u_gasHoldup * (pow(chamber.at(0).dimensionsInternal.at(0).first, 2) / pow(chamber.at(section).dimensionsInternal.at(0).first, 2));
	else
		return u_gasHoldup;
}

double CDryerBatch::GetParticleEquilibriumMoistureContent(double temperature, double RH) const
{
	if (eqData.RHs.empty() || eqData.temperatures.empty() || eqData.equilibriumMoistureContents.empty() || temperature != temperature || RH != RH)
		return -1; // Error
	std::set<double>::iterator itRHmax = eqData.RHs.end();
	std::advance(itRHmax, -1);
	// RH bounds [0,1] check
	if (RH > 1)
		RH = 1;
	if (RH < 0)
		RH = 0;

	// Temperature bounds [0,] check
	if (temperature < 0)
		temperature = 0;

	std::set<double>::iterator itT, itRH, itlowT, itupT, itlowRH, itupRH;

	itlowT = eqData.temperatures.lower_bound(temperature);
	if (itlowT == eqData.temperatures.end() || *itlowT != temperature)
		std::advance(itlowT, -1);
	itupT = eqData.temperatures.upper_bound(temperature);
	if (itlowT == eqData.temperatures.end())
	{
		itlowT = itupT;
		std::advance(itupT, 1);
	}
	if (itupT == eqData.temperatures.end())
	{
		itupT = itlowT;
		std::advance(itlowT, -1);
	}

	itlowRH = eqData.RHs.lower_bound(RH);
	if (itlowRH == eqData.RHs.end() || *itlowRH != RH)
		std::advance(itlowRH, -1);
	itupRH = eqData.RHs.upper_bound(RH);
	if (itlowRH == eqData.RHs.end())
	{
		itlowRH = itupRH;
		std::advance(itupRH, 1);
	}
	if (itupRH == eqData.RHs.end())
	{
		itupRH = itlowRH;
		std::advance(itlowRH, -1);
	}

	double Xeq;
	if (itlowT == eqData.temperatures.end() || itupT == eqData.temperatures.end() || itlowRH == eqData.RHs.end() || itupRH == eqData.RHs.end())
	{
		
		if (itlowT == eqData.temperatures.end() || itupT == eqData.temperatures.end())
		{
			itT = itlowT != eqData.temperatures.end() ? itlowT : itupT;
			if (itlowRH == eqData.RHs.end() || itupRH == eqData.RHs.end())
			{
				itRH = itlowRH != eqData.RHs.end() ? itlowRH : itupRH;
				Xeq = eqData.equilibriumMoistureContents.at(std::make_pair(*itT, *itRH));
				return Xeq;
			}
			else
			{
				double XeqTRHl = eqData.equilibriumMoistureContents.at(std::make_pair(*itT, *itlowRH));
				double XeqTRHu = eqData.equilibriumMoistureContents.at(std::make_pair(*itT, *itupRH));
				Xeq = lerp(XeqTRHl, XeqTRHu, (RH - *itlowRH) / (*itupRH - *itlowRH));
				return Xeq;
			}

		}
		else
		{
			itRH = itlowRH != eqData.RHs.end() ? itlowRH : itupRH;
			double XeqTlRH = eqData.equilibriumMoistureContents.at(std::make_pair(*itlowT, *itRH));
			double XeqTuRH = eqData.equilibriumMoistureContents.at(std::make_pair(*itupT, *itRH));
			Xeq = lerp(XeqTlRH, XeqTuRH, (temperature - *itlowT) / (*itupT - *itlowT));
			return Xeq;
		}
		return -1;
	}
	else
	{
		double XeqTlRHl, XeqTuRHl, XeqTlRHu, XeqTuRHu, XeqTl, XeqTu;

		XeqTlRHl = eqData.equilibriumMoistureContents.at(std::make_pair(*itlowT, *itlowRH));
		XeqTuRHl = eqData.equilibriumMoistureContents.at(std::make_pair(*itupT, *itlowRH));

		/*if (eqData.equilibriumMoistureContents.contains(std::make_pair(*itlowT, *itupRH)))
			XeqTlRHu = eqData.equilibriumMoistureContents.at(std::make_pair(*itlowT, *itupRH));
		else
			XeqTlRHu = eqData.equilibriumMoistureContents.at(std::make_pair(*itlowT, *eqData.maxRH.at(*itlowT)));*/

		try
		{
			XeqTlRHu = eqData.equilibriumMoistureContents.at(std::make_pair(*itlowT, *itupRH));
		}
		catch (std::out_of_range const& e)
		{
			XeqTlRHu = eqData.equilibriumMoistureContents.at(std::make_pair(*itlowT, *eqData.maxRH.at(*itlowT)));
		}

		/*if (eqData.equilibriumMoistureContents.contains(std::make_pair(*itupT, *itupRH)))
			XeqTuRHu = eqData.equilibriumMoistureContents.at(std::make_pair(*itupT, *itupRH));
		else
			XeqTuRHu = eqData.equilibriumMoistureContents.at(std::make_pair(*itupT, *eqData.maxRH.at(*itupT)));*/

		try
		{
			XeqTuRHu = eqData.equilibriumMoistureContents.at(std::make_pair(*itupT, *itupRH));
		}
		catch (std::out_of_range const& e)
		{
			XeqTuRHu = eqData.equilibriumMoistureContents.at(std::make_pair(*itupT, *eqData.maxRH.at(*itupT)));
		}

		XeqTl = lerp(XeqTlRHl, XeqTlRHu, (RH - *itlowRH) / (*itupRH - *itlowRH));
		XeqTu = lerp(XeqTuRHl, XeqTuRHu, (RH - *itlowRH) / (*itupRH - *itlowRH));
		Xeq = lerp(XeqTl, XeqTu, (temperature - *itlowT) / (*itupT - *itlowT));
		return Xeq;
	}
	return -1;
}

double CDryerBatch::GetEquilibriumRelativeHumidity(double temperature, double X) const
{
	if (eqData.RHs.empty() || eqData.temperatures.empty() || eqData.equilibriumMoistureContents.empty() || temperature != temperature || X != X)
		return -1; // Error
	//X = X / moistureScaler;

	if (X < 0)
		return 0;

	std::set<double>::iterator itlowT, itupT, itT, itRH;
	double XlowT, XupT, Xpolated, XT;

	itRH = eqData.RHs.end();
	std::advance(itRH, -1);

	itlowT = eqData.temperatures.lower_bound(temperature);
	if (itlowT == eqData.temperatures.end() || *itlowT != temperature)
		std::advance(itlowT, -1);
	itupT = eqData.temperatures.upper_bound(temperature);
	if (itlowT == eqData.temperatures.end())
	{
		itlowT = itupT;
		std::advance(itupT, 1);
	}
	if (itupT == eqData.temperatures.end())
	{
		itupT = itlowT;
		std::advance(itlowT, -1);
	}

	if (itlowT == eqData.temperatures.end() || itupT == eqData.temperatures.end())
	{
		itT = itlowT != eqData.temperatures.end() ? itlowT : itupT;
		//if (*itRH > *eqData.maxRH.at(*itT))
		//	itRH = eqData.maxRH.at(*itT);
		XT = eqData.equilibriumMoistureContents.at(std::make_pair(*itT, *itRH));
		if (X > XT)
			return 1;
		while (X < XT)
		{
			std::advance(itRH, -1);
			if (itRH == eqData.RHs.end())
				return -1;
			XT = eqData.equilibriumMoistureContents.at(std::make_pair(*itT, *itRH));
			if (X > XT)
				return *itRH;
		}
	}
	else
	{
		XlowT = eqData.equilibriumMoistureContents.at(std::make_pair(*itlowT, *itRH));
		XupT = eqData.equilibriumMoistureContents.at(std::make_pair(*itupT, *itRH));
		Xpolated = lerp(XlowT, XupT, (temperature - *itlowT) / (*itupT - *itlowT));
		if (X > Xpolated)
			return 1;
		while (X < XlowT)
		{
			std::advance(itRH, -1);
			if (itRH == eqData.RHs.end())
				return -1;
			XlowT = eqData.equilibriumMoistureContents.at(std::make_pair(*itlowT, *itRH));
			XupT = eqData.equilibriumMoistureContents.at(std::make_pair(*itupT, *itRH));
			Xpolated = lerp(XlowT, XupT, (temperature - *itlowT) / (*itupT - *itlowT));
			if (X > Xpolated)
				return *itRH;
		}
	}
	return -1;
}

bool CDryerBatch::InitializeMoistureContentDatabase(std::string path)
{
	//CCorrelation::GetParameters();
	std::ifstream  file(path);
	if (!file.is_open())
		RaiseError("Moisture content data file has incorrect path.");
	std::string line;
	std::vector<std::vector<std::string>> parsedCsv;
	while (std::getline(file, line))
	{
		std::stringstream lineStream(line);
		std::string cell;
		std::vector<std::string> parsedRow;
		while (std::getline(lineStream, cell, ','))
			parsedRow.push_back(cell);

		parsedCsv.push_back(parsedRow);
	}

	for (auto & element : parsedCsv.at(0))
		if (element != "")
			eqData.RHs.insert(std::stod(element));

	for (auto & element : parsedCsv)
		if (element.at(0) != "")
			eqData.temperatures.insert(std::stod(element.at(0)));

	for (int i = 1; i < parsedCsv.size(); i++)
	{
		for (int j = 1; j < parsedCsv.at(i).size(); j++)
			eqData.equilibriumMoistureContents.insert(std::make_pair(std::make_pair(std::stod(parsedCsv.at(i).at(0)), std::stod(parsedCsv.at(0).at(j))), std::stod(parsedCsv.at(i).at(j))));
		std::set<double>::iterator itRHmax = eqData.RHs.end();
		int fromMaxRH = eqData.RHs.size() - parsedCsv.at(i).size() + 2;
		std::advance(itRHmax, -fromMaxRH);
		eqData.maxRH.insert(std::make_pair(std::stod(parsedCsv.at(i).at(0)), itRHmax));
	}
	return true;
}

double CDryerBatch::DetermineSectionsFilledWithBed(double _time, double particleTemperature)
{
	if (particlesGlobal == false)
		return 0;
	double bedHeight = CalculateBedHeight(_time, particleTemperature);
	int sectionFilled = 0;
	double chamberHeight = 0;
	while (chamberHeight < bedHeight)
	{
		chamberHeight += chamber.at(sectionFilled).height;
		sectionFilled++;
		if (sectionFilled > chamber.size())
			return -1; // Bed larger than chamber
	}
	double filledHeight = chamberHeight - chamber.at(sectionFilled).height;
	return sectionFilled -1 + (bedHeight - filledHeight) / chamber.at(sectionFilled).height;
}

//double CDryerBatch::GetGasMassOfLayer(double _time, size_t layer, double gasTemperature, double particleTemperature)
//{
//	size_t sectionContainingLayer = DetermineSectionForLayer(layer);
//	size_t layersInPrevSections = 0;
//	for (size_t i = 0; i < sectionContainingLayer; i++)
//		layersInPrevSections += chamber.at(i).layers;
//	double sectionsFilledWithBed = DetermineSectionsFilledWithBed(_time, particleTemperature);
//	// double sectionsFilledWithBed = CalculateBedHeightOrDetermineSectionsFilledWithBed(_time, particleTemperature, false);
//	double sectionVolume = CalculateSectionVolume(sectionContainingLayer);
//	double sectionGasMass;
//	double heightUsage, fullSections;
//	heightUsage = modf(sectionsFilledWithBed, &fullSections);
//	double sectionGasMassVoid = sectionVolume * GetAvgTPCompoundProperty(_time, EPhase::GAS, ECompoundTPProperties::DENSITY, gasTemperature, m_holdup->GetPressure(_time));
//	double eps = CalculateBedPorosity(_time);
//	double sectionGasMassEPS = sectionVolume * eps * GetAvgTPCompoundProperty(_time, EPhase::GAS, ECompoundTPProperties::DENSITY, gasTemperature, m_holdup->GetPressure(_time));
//	double particleVolume = m_holdup->GetPhaseMass(_time, EPhase::SOLID) / GetAvgTPCompoundProperty(_time, EPhase::SOLID, ECompoundTPProperties::DENSITY, particleTemperature, m_holdup->GetPressure(_time));
//	double eps2 = 1.-particleVolume / sectionVolume;
//	if (sectionContainingLayer > fullSections)
//		sectionGasMass = sectionGasMassVoid;
//	else if (sectionContainingLayer < fullSections)
//		sectionGasMass = sectionGasMassEPS;
//	else
//		sectionGasMass = heightUsage * sectionGasMassEPS + (1 - heightUsage) * sectionGasMassVoid;
//
//	if (chamber.at(sectionContainingLayer).shape != EShape::CYLINDRICAL)
//		RaiseError("Layer mass could not be determined for section with non cylindrical shape.\nMissing implementation.");
//
//	const double r = chamber.at(sectionContainingLayer).dimensionsInternal.at(0).first / 2;
//	const double R = chamber.at(sectionContainingLayer).dimensionsInternal.at(0).second / 2;
//	const double h = chamber.at(sectionContainingLayer).height;
//	const double tanWallAngle = (R - r) / h;
//	const double r1 = r + h * (layer - layersInPrevSections) / chamber.at(sectionContainingLayer).layers * tanWallAngle;
//	const double R1 = r + h * (layer - layersInPrevSections + 1) / chamber.at(sectionContainingLayer).layers * tanWallAngle;
//	const double layerVolume = chamber.at(sectionContainingLayer).height * MATH_PI * (R1*R1+R1*r1+r1*r1) / (3*chamber.at(sectionContainingLayer).layers);
//	const double layerGasMass = sectionGasMass * layerVolume / sectionVolume;
//	return layerGasMass;
//}

std::vector<double> CDryerBatch::GetGasMassOfLayers(double _time, double gasTemperature, double particleTemperature)
{
	std::vector<double> layerMasses(N_total, 0);
	double chamberVolume = 0;
	for (size_t section = 0; section < chamber.size(); section++)
		chamberVolume += CalculateSectionVolume(section);

	double rhoParticle = GetAvgTPCompoundProperty(_time, EPhase::SOLID, ECompoundTPProperties::DENSITY, particleTemperature, m_holdup->GetPressure(_time));
	double rhoGas = GetAvgTPCompoundProperty(_time, EPhase::GAS, ECompoundTPProperties::DENSITY, gasTemperature, m_holdup->GetPressure(_time));

	const double particlesVolume = m_holdup->GetPhaseMass(_time, EPhase::SOLID) / rhoParticle;
	const double chamberGasMass = (chamberVolume - particlesVolume) * rhoGas;
	double sectionsFilledWithBed = DetermineSectionsFilledWithBed(_time, particleTemperature);
	//	double sectionsFilledWithBed = CalculateBedHeightOrDetermineSectionsFilledWithBed(_time, particleTemperature,false);
	double heightUsage, fullSections;
	heightUsage = modf(sectionsFilledWithBed, &fullSections);
	double sectionVolume = 0;
	double sectionGasMass = 0;
	std::vector<double>::reverse_iterator ritLayerMasses = layerMasses.rbegin();
	std::vector<chamberSection>::reverse_iterator ritChamberSections = chamber.rbegin();
	double assosiatedGasMass = 0;
	ritChamberSections++;
	double test = chamber.rend() - ritChamberSections;
	
	// Empty section above bed
	for (size_t section = chamber.size() -1; section > fullSections; section--)
	{
		sectionVolume = CalculateSectionVolume(section);
		sectionGasMass = sectionVolume * rhoGas;
		std::vector<double> sectionLayerGasMasses(chamber.at(section).layers, 0);
		const double r = chamber.at(section).dimensionsInternal.at(0).first / 2;
		const double R = chamber.at(section).dimensionsInternal.at(0).second / 2;
		const double h = chamber.at(section).height;
		const double tanWallAngle = (R - r) / h;
		for (size_t layer = 0; layer < chamber.at(section).layers; layer++)
		{
			const double r1 = r + h * (layer) / chamber.at(section).layers * tanWallAngle;
			const double R1 = r + h * (layer + 1) / chamber.at(section).layers * tanWallAngle;
			const double layerVolume = CalculateLayerVolume(section, R1, r1);
			sectionLayerGasMasses.at(layer) = sectionGasMass * layerVolume / sectionVolume;
		}
		std::vector<double>::reverse_iterator ritSectionLayerMasses = sectionLayerGasMasses.rbegin();
		for (; ritSectionLayerMasses != sectionLayerGasMasses.rend(); ritSectionLayerMasses++)
		{
			*ritLayerMasses = *ritSectionLayerMasses;
			ritLayerMasses++;
			assosiatedGasMass += *ritSectionLayerMasses;
		}
		
	}
	// Should be working.


	size_t layersWithParticles = 0;
	

	// Last section filled with bed
	{
		size_t section = fullSections;
		size_t layersUsed = std::ceil(heightUsage * chamber.at(section).layers);
		layersWithParticles += layersUsed;
		for (size_t i = 0; i < fullSections; i++)
			layersWithParticles += chamber.at(i).layers;

		std::vector<double> sectionLayerGasMasses(chamber.at(section).layers, 0);
		const double r = chamber.at(section).dimensionsInternal.at(0).first / 2;
		const double R = chamber.at(section).dimensionsInternal.at(0).second / 2;
		const double h = chamber.at(section).height;
		const double tanWallAngle = (R - r) / h;

		for (size_t layer = layersUsed; layer < chamber.at(section).layers; layer++)
		{
			const double r1 = r + h * (layer) / chamber.at(section).layers * tanWallAngle;
			const double R1 = r + h * (layer + 1) / chamber.at(section).layers * tanWallAngle;
			const double layerVolume = CalculateLayerVolume(section, R1, r1);
			sectionLayerGasMasses.at(layer) = sectionGasMass * layerVolume / sectionVolume;
		}

		for (size_t layer = 0; layer < layersUsed; layer++)
		{
			const double r1 = r + h * (layer) / chamber.at(section).layers * tanWallAngle;
			const double R1 = r + h * (layer + 1) / chamber.at(section).layers * tanWallAngle;
			const double layerVolume = CalculateLayerVolume(section, R1, r1) - (particlesVolume / layersWithParticles);
			sectionLayerGasMasses.at(layer) = sectionGasMass * layerVolume / sectionVolume;
		}

		std::vector<double>::reverse_iterator ritSectionLayerMasses = sectionLayerGasMasses.rbegin();
		for (; ritSectionLayerMasses != sectionLayerGasMasses.rend(); ritSectionLayerMasses++)
		{
			*ritLayerMasses = *ritSectionLayerMasses;
			ritLayerMasses++;
			assosiatedGasMass += *ritSectionLayerMasses;
		}
	}

	// Sections filled completely with particles
	if (fullSections > 0)
		for (double section = fullSections - 1; section > -1; section--)
		{
			sectionVolume = CalculateSectionVolume(section);
			sectionGasMass = sectionVolume * rhoGas;
			std::vector<double> sectionLayerGasMasses(chamber.at(section).layers, 0);
			const double r = chamber.at(section).dimensionsInternal.at(0).first / 2;
			const double R = chamber.at(section).dimensionsInternal.at(0).second / 2;
			const double h = chamber.at(section).height;
			const double tanWallAngle = (R - r) / h;
			for (size_t layer = 0; layer < chamber.at(section).layers; layer++)
			{
				const double r1 = r + h * (layer) / chamber.at(section).layers * tanWallAngle;
				const double R1 = r + h * (layer + 1) / chamber.at(section).layers * tanWallAngle;
				const double layerVolume = CalculateLayerVolume(section,R1,r1) - (particlesVolume / layersWithParticles);
				sectionLayerGasMasses.at(layer) = sectionGasMass * layerVolume / sectionVolume;
			}
			std::vector<double>::reverse_iterator ritSectionLayerMasses = sectionLayerGasMasses.rbegin();
			for (; ritSectionLayerMasses != sectionLayerGasMasses.rend(); ritSectionLayerMasses++)
			{
				*ritLayerMasses = *ritSectionLayerMasses;
				ritLayerMasses++;
				assosiatedGasMass += *ritSectionLayerMasses;
			}

		}

	if (chamberGasMass - assosiatedGasMass > 1e-16)
		RaiseError("Not all gas mass was accounted for in combined layer gas mass.");

	return layerMasses;
}

//size_t CDryerBatch::DetermineSectionForLayer(size_t layer) const
//{
//	size_t section = 0;
//	size_t numOfPrevLayers = chamber.at(section).layers;;
//	while (layer >= numOfPrevLayers)
//	{
//		section++;
//		numOfPrevLayers += chamber.at(section).layers;
//	}
//	return section;
//}

///  Main calculation functions ///
void CUnitDAEModel::CalculateResiduals(double _time, double* _vars, double* _ders, double* _res, void* _unit)
{
	const auto* unit = static_cast<CDryerBatch*>(_unit);
	auto* unit2 = static_cast<CDryerBatch*>(_unit);

	// Codes for checking _vars, _ders and _res
	std::vector vars(_vars, _vars + GetVariablesNumber());
	std::vector ders(_ders, _ders + GetVariablesNumber());
	std::vector  res(_res, _res + GetVariablesNumber());

	bool NaNvalue = false;
	for (int i = 0; i < GetVariablesNumber(); i++)
		if (_vars[i] != _vars[i])
			NaNvalue = true;

if (debugToggle){
	if (_time <= 1)
		counter++;
	progressCounterTotal++;
} // _DEBUG
	if (/*_time > 7929|| _time > 6000 || _time > 5948 || _time > 300 ||*/ _time > 2797.52)
		bool breakint = true;

	/// Read input parameters ///
	const pressure mPressure = unit->m_holdup->GetPressure(_time); // Pressure holdup [Pa]
	const massFlow mFlowSuspension = unit->m_inLiquidStream->GetMassFlow(_time);
	const massFlow mFlowSuspensionLiquid = unit->m_inLiquidStream->GetPhaseMassFlow(_time, EPhase::LIQUID); // Suspension mass flow [kg/s]
	const massFlow mFlowSuspensionSolid = unit->m_inLiquidStream->GetPhaseMassFlow(_time, EPhase::SOLID);
	const massFlow mFlowSuspensionGas = unit->m_inNozzleAirStream->GetPhaseMassFlow(_time, EPhase::GAS);
	//const double mFlowGasSupplemental = unit->m_inSuspensionStream->GetPhaseMassFlow(_time, EPhase::GAS); [kg/s]
	const massFlow mFlowGas = unit->m_inGasStream->GetPhaseMassFlow(_time, EPhase::GAS); // Gas mass flow [kg/s]
	const mass mHoldupSolid = unit->m_holdup->GetPhaseMass(_time, EPhase::SOLID); // Solid mass holdup [kg]
	const mass mHoldupGas = unit->m_holdup->GetPhaseMass(_time, EPhase::GAS); // Gas mass holdup [kg]
	const mass mHoldupLiquid = unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID); // Liquid mass holdup [kg]
	//const length H_bed = unit->H_bed; // Bed height [m]
	//const length d_bed = unit->d_bed; // Bed diameter [m]
	const moistureContent Y_in = unit->Y_in; // Y in [kg/kg dry]
	const area A_P = unit->A_P; // Total particle surface [m^2]
	const area varA_p = A_P; // Total particle surface [m^2]
	const double derA_p = 0; // Change of total particle surface [m^2/s]
	/*	If changed to a granulator, these varaibles need to change:
		phi, d/dt particle temperature, A_P, d/dt A_P
	*/
	const length Delta_f = unit->Delta_f; // Max. liquid film thickness [m]
	const massFraction w_l = unit->m_inLiquidStream->GetPhaseFraction(_time, EPhase::LIQUID);// unit->w_l; // Liquid mass fraction suspension [kg/kg]

	const temperature mGasInTemperature = unit->m_inGasStream->GetTemperature(_time); // Temperature of inGasStream general [K]
	const temperature mSuspensionTemperature = unit->m_inLiquidStream->GetTemperature(_time);

	/// DAE system ///
	// Values in _vars: determined by the solver & should not be changed, set as const!

	const temperature varTempParticle2 = _vars[m_iTempParticle]; // Temperature of the particles [K]
	const temperature varTempFilm2 = _vars[m_iTempFilm]; // Temperature of the film [K]

	//const double varPhi2 = _vars[m_iPhi2]; // Liquid film surface coverage

	size_t sectionsWithParticles = unit->particlesGlobal? std::ceil(unit2->DetermineSectionsFilledWithBed(_time, varTempParticle2)):0;
	// size_t sectionsWithParticles = unit->particlesGlobal ? std::ceil(unit2->CalculateBedHeightOrDetermineSectionsFilledWithBed(_time, varTempParticle2, false)) : 0;
	size_t layersWithParticles = 0;
	for (size_t section = 0; section < sectionsWithParticles; section++)
		layersWithParticles += unit->chamber.at(section).layers;
	
	// Calculated constants

	const moistureContent X = mHoldupLiquid / mHoldupSolid; // Moisture content of the bed [kg water per kg dry particle] 

	// Calculate other coefficients
	const double alpha_GP = unit->CalcAlphaParticleGas(_time); // Heat transfer coefficient gas to particle
	const double alpha_GF = alpha_GP; // Heat transfer coefficient gas to film
	const double alpha_PF = alpha_GP * unit->f_alpha;//unit->CalculateAlphaPF(_time); // Heat transfer coefficient particle to film


	// Model according to Dosta 2010

	// Height averaged temperature of the gas phase [K]
	// Dosta 2010 eq. 14a

	const temperature varAvTempGas = CalculateAverage(_vars,m_iTempGas, unit->N_particle);

	// Height averaged moisture content of the gas phase [kg/kg]
	// Dosta 2010 eq. 14b
	const moistureContent Y_av = CalculateAverage(_vars, m_iY_g, unit->N_particle);

	double rhoGas = unit->rhoGas; // Average density of gas phase for average gas temperature and system pressure [kg/m^3]
	rhoGas = unit->GetAvgTPCompoundProperty(_time, EPhase::GAS, ECompoundTPProperties::DENSITY, varAvTempGas, mPressure);
	double rhoLiquid = unit->rhoWater; // Average density of gas phase for average gas temperature and system pressure [kg/m^3]
	rhoLiquid = unit->GetAvgTPCompoundProperty(_time, EPhase::LIQUID, ECompoundTPProperties::DENSITY, varTempFilm2, mPressure);
	double rhoSolid= unit->rhoParticle; // Average density of solid phase for average particle temperature and system pressure [kg/m^3]
	rhoSolid = unit->GetAvgTPCompoundProperty(_time, EPhase::SOLID, ECompoundTPProperties::DENSITY, varTempParticle2, mPressure);

	double cPLiquidGaseous = unit->C_PWaterVapor;
	double cPSolid = unit->C_PParticle;
	cPSolid = unit->GetAvgTPCompoundProperty(_time, EPhase::SOLID, ECompoundTPProperties::HEAT_CAPACITY_CP, varTempParticle2, mPressure);

	const double D = unit->CalculateDiffusionCoefficient(_time, varAvTempGas, varTempFilm2, mPressure);
	const double beta = !unit->particlesGlobal? 0: unit->CalculateBetaPA(_time, varAvTempGas, varTempParticle2, D); // Mass transfer coefficient for film to gas [m/s]
	if (beta == -1)
		unit2->RaiseError("Mass transfer coefficient calculation error. Provided temperature = nan or < 0 K");

	double Y_sat = unit2->GetGasSaturationMoistureContent(varAvTempGas, mPressure); // Saturation moisture content of the gas (average) [kg/kg]

	// Normalized dyring curve [-]
	// Implementation according to Lehmann 2021
	double normalizedDryingCurve = 1;
	const double avgRH = unit2->GetRelativeHumidity(Y_av, varAvTempGas, mPressure); // Average relativ humidity of gas phase [-]

	double Xeq = !unit->particlesGlobal ? 0 : unit2->CalcuateSolidEquilibriumMoistureContent(_time, varTempParticle2, avgRH); // Equilibrium moisture content of particle [kg liquid/ kg dry solid]
	normalizedDryingCurve = unit2->CalculateNormalizedDryingCurve(X, Xeq);
	//Transfer = Xeq;

	bool liquidFilm = !unit->particlesGlobal ? false : normalizedDryingCurve > unit->phiCuttOff; // Liquid film of particle surface
	if (liquidFilm)
		normalizedDryingCurve = 1;

	double varPhi2 = 0;
	const double prevTime = unit->m_holdup->GetPreviousTimePoint(_time);

	const double RHeq = !unit->particlesGlobal ? 0 : unit->GetEquilibriumRelativeHumidity(varTempParticle2,X);// CalcuateSolidEquilibriumMoistureContent^-1 from Temp+X to RH
	const double Y_eq = this->Y_eq?unit->CalculateGasEquilibriumMoistureContent(varTempParticle2, mPressure, RHeq): Y_sat;

	// Mass stream of evaporating liquid gas side [kg/s]
	// Dosta 2010 eq. 16
	// Replaced with not DEA variable as not necessary
	double MFlowVaporGasSide = normalizedDryingCurve * beta * varA_p * rhoGas * (Y_eq - Y_av);
	// varA_p could be change to also represent swelling (X) and liquid film presence (phi,deltaF)
	if (MFlowVaporGasSide < 0)
		MFlowVaporGasSide = 0;

	//Transfer = MFlowVaporGasSide;
	//double MFlowVaporLiquidSide = DBL_MAX; // Mass to be evaporated liquid side [kg]

	// ToDo - Fix as kg is compared with kg/s
	//unit2->workingHoldup->SetMass(_time, 0);
	//unit2->workingHoldup->SetPhaseMass(_time,EPhase::LIQUID, (X - Xeq) * mHoldupSolid);
	//unit2->workingHoldup->AddStream(prevTime, _time, unit->m_inSuspensionStream);
	//if (unit2->workingHoldup->GetPhaseMass(_time, EPhase::LIQUID) - MFlowVaporGasSide * (_time - prevTime) < 0)
	//{
	//	MFlowVaporLiquidSide = unit2->workingHoldup->GetPhaseMass(_time, EPhase::LIQUID) / (_time - prevTime);
	//	liquidSideLimitedGlobal = true;
	//}
	//unit2->workingHoldup->SetMass(_time, 0);

	double MFlowVaporLiquidSide = (X - Xeq) * mHoldupSolid / (_time - prevTime) + mFlowSuspensionLiquid;
	MFlowVaporLiquidSide = (X - Xeq) * mHoldupSolid;
	if (MFlowVaporLiquidSide < 0)
		MFlowVaporLiquidSide = 0;

	//MFlowVaporLiquid = MFlowVaporLiquidSide;
	//liquidSideLimitedGlobal = MFlowVaporLiquidSide < MFlowVaporGasSide;
	const double MFlowVapor = std::min(MFlowVaporGasSide, MFlowVaporLiquidSide);
	VaporFlowStorage = MFlowVapor;

	const double XCutt = unit->REAinv(Xeq, 1. - unit->phiCuttOff);

	//liquidFilm = liquidFilm ? ((unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID) - ((_time - prevTime) * MFlowVapor) - XCutt * mHoldupSolid) / (unit->rhoWater * Delta_f * unit->A_P) > 0) : false;

	if (liquidFilm)
	{
		unit->m_holdup->AddStream(prevTime, _time, unit->m_inLiquidStream);
		varPhi2 = (unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID) - XCutt * mHoldupSolid - ((_time - prevTime) * MFlowVapor)) / (rhoLiquid * Delta_f * unit->A_P);
		unit->m_holdup->SetPhaseMass(_time, EPhase::LIQUID, mHoldupLiquid);
		unit->m_holdup->SetPhaseMass(_time, EPhase::SOLID, mHoldupSolid);
	}
	if (varPhi2 < 0)
		liquidFilm = false;
	TransferBool = liquidFilm;

	std::vector<double> Q_AP(unit->N_particle, 0);// Heat transfered from air to particle [W]
	double Q_AP3 = 0; // Sum of heat transfered from air to particle [W]
	if (unit->particlesGlobal)
	{
		for (int i = 0; i < unit->N_particle; i++)
			Q_AP[i] = alpha_GP * A_P / unit->N_particle * (1 - varPhi2) * (_vars[m_iTempGas + i] - varTempParticle2);
		
		for (int i = 0; i < Q_AP.size(); i++)
			Q_AP3 += Q_AP[i];
	}
	// Dosta 2010 eq. 23 - changed to prevent cooling of higher layers by subtracting the same amount of energy
	// Original was A_P / 1

	std::vector<double> Q_AF(unit->N_particle, 0);// Heat transfered from air to film [W]
	double Q_AF2 = 0;// Sum of heat transfered from air to film
	if (unit->particlesGlobal)
	{
		for (int i = 0; i < Q_AF.size(); i++)
			Q_AF[i] = alpha_GF * A_P / unit->N_particle * varPhi2 * (_vars[m_iTempGas + i] - varTempFilm2);
		
		for (int i = 0; i < Q_AF.size(); i++)
			Q_AF2 += Q_AF[i];
	}
	// Dosta 2010 eq. 21 - changed to prevent cooling of higher layers by subtracting the same amount of energy
	// Original was A_P / 1

	const double Q_PF2 = !unit->particlesGlobal?0:alpha_PF * A_P * varPhi2 * (varTempParticle2 - varTempFilm2); // Heat transfered from particle to film [W]
	// Dosta 2010 eq. 22

	std::pair<double, std::vector<double>> resultHeatLoss = CalculateChamberHeatLoss(_time, _unit, _vars); // Overall heat trasfer throught wall to environment
	const double Q_PW = resultHeatLoss.first; // Heat transfered from particles to wall and environment [W]
	std::vector<double> Q_GW = resultHeatLoss.second; // Heat transfered from gas to wall  and environment [W]
	heatLossTransfer = Q_PW;
	for (int i = 0; i < Q_GW.size(); i++)
		heatLossTransfer += Q_GW.at(i);
	
	const double HFlowVapor2 = MFlowVapor * (unit->C_PWaterVapor * varAvTempGas + unit->Delta_h0); // Enthaly stream of evporating liquid [J/s]
	// Dosta 2010 eq. 15
	const double H_susp2 =!unit->particlesGlobal?0: mSuspensionTemperature * (mFlowSuspensionLiquid * unit->C_PWaterLiquid * mFlowSuspensionSolid * unit->C_PParticle); // Enthaly stream of entering suspension [J/s]
	// Dosta 2010 eq. 17
	//const double H_gasSupp = unit->m_inSuspensionStream->GetTemperature(_time) * mFlowGasSupplemental * unit->C_PGas; // Enthaly stream of entering gas supplement [J/s]

	if (w_l == 0)
		unit2->RaiseError("Liquid mass fraction in suspension is 0.");
	
	const double varPhiFormula2 = (mFlowSuspensionLiquid - MFlowVapor) / (unit->rhoWater * Delta_f * varA_p) - varPhi2 / varA_p * derA_p; // Change in the degree of wetness of the particle [-/s]
	// Dosta 2010 eq. 13
	//_res[m_iPhi2] = _ders[m_iPhi2] - varPhiFormula2;
	const double dersPhi2 = varPhiFormula2;
	
	// Total particle surface area
	//const double varAFormula = 0; // Change of total particle surface area due to ganulation and attrition 
	//_res[m_iA_p] = derA_p - varAFormula;

	std::vector<double> HFlowInGas2(unit->N_total, 0); // Storage vector for enthaly stream entering height layers
	std::vector<double> HFlowOutGas2(unit->N_total, 0); // Storage vector for enthaly stream exiting height layers
	std::vector<double> sumEnergyStreamsGas(unit->N_total, 0);

	ParallelFor(unit->N_total, [&](size_t i)
		{
			double HFlowInGasFormula = 0; // Enthaly stream entering height layer Dosta
			if (i == 0)
				HFlowInGasFormula = mFlowGas * (unit->C_PGas * mGasInTemperature + Y_in * (unit->C_PWaterVapor * mGasInTemperature + unit->Delta_h0));
				// Dosta 2010 eq. 18 - addjusted for first layer with gas input
			else
				HFlowInGasFormula = (mFlowGas + ((i > unit->suspLayer) ? mFlowSuspensionGas : 0)) * (unit->C_PGas * _vars[m_iTempGas + (i - 1)] + _vars[m_iY_g + (i - 1)] * (unit->C_PWaterVapor * _vars[m_iTempGas + (i - 1)] + unit->Delta_h0));
				// Dosta 2010 eq. 18

			if (i == unit->suspLayer)
				HFlowInGasFormula += mFlowSuspensionGas * (unit->C_PGas * mSuspensionTemperature + Y_in * (unit->C_PWaterVapor * mSuspensionTemperature + unit->Delta_h0));

			HFlowInGas2[i] = HFlowInGasFormula;

			const double HFlowOutGasFormula = (mFlowGas + ((i >= unit->suspLayer) ? mFlowSuspensionGas : 0)) * (unit->C_PGas * _vars[m_iTempGas + i] + _vars[m_iY_g + i] * (unit->C_PWaterVapor * _vars[m_iTempGas + i] + unit->Delta_h0)); // Enthaly stream exiting height layer
			// Dosta 2010 eq. 19
			HFlowOutGas2[i] = HFlowOutGasFormula;

			if (i < unit->N_particle)
				sumEnergyStreamsGas[i] = HFlowInGas2[i] - HFlowOutGas2[i] - Q_AP[i] - Q_AF[i] - Q_GW[i] + (HFlowVapor2/* + H_gasSupp */) / unit->N_particle;
			else
				sumEnergyStreamsGas[i] = HFlowInGas2[i] - HFlowOutGas2[i] - Q_GW[i];
		});

	//std::vector<double> layerGasMasses = unit2->GetGasMassOfLayers(unit->heighestFlowTimepoint, varAvTempGas, varTempParticle2);


	ParallelFor(unit->N_total, [&](size_t i)
		{
			std::vector vars(_vars, _vars + GetVariablesNumber());
			std::vector ders(_ders, _ders + GetVariablesNumber());
			std::vector  res(_res, _res + GetVariablesNumber());

			double varYOutFormula2 = 0; // Change of gas moisture content
			if (i == 0)
				varYOutFormula2 = -(mFlowGas + ((i >= unit->suspLayer) ? mFlowSuspensionGas : 0)) * unit->N_total / mHoldupGas * (_vars[m_iY_g + i] - Y_in) + MFlowVapor / mHoldupGas;
				// Dosta 2010 eq. 12 addjusted for first layer with gas input
			else if (i < unit->N_particle)
				varYOutFormula2 = -(mFlowGas + ((i >= unit->suspLayer) ? mFlowSuspensionGas : 0)) * unit->N_total / mHoldupGas * (_vars[m_iY_g + i] - _vars[m_iY_g + (i - 1)]) + MFlowVapor / mHoldupGas;
				// Dosta 2010 eq. 12
			else
				varYOutFormula2 = -(mFlowGas + ((i >= unit->suspLayer) ? mFlowSuspensionGas : 0)) * unit->N_total / mHoldupGas * (_vars[m_iY_g + i] - _vars[m_iY_g + (i - 1)]);
			if (i == unit->suspLayer)
				varYOutFormula2 += mFlowSuspensionGas * Y_in * unit->N_total / mHoldupGas;
			_res[m_iY_g + i] = _ders[m_iY_g + i] - varYOutFormula2;

			// Change of gas temperature
			const double derTempGasFormula2 = unit->N_total / (mHoldupGas * (unit->C_PGas + unit->C_PWaterVapor * _vars[m_iY_g + i])) * sumEnergyStreamsGas[i]
											- _ders[m_iY_g + i] * ((unit->C_PWaterVapor * _vars[m_iTempGas + i] + unit->Delta_h0) / (unit->C_PGas + unit->C_PWaterVapor * _vars[m_iY_g + i]));
			// Dosta 2010 eq. 10
			
			_res[m_iTempGas + i] = _ders[m_iTempGas + i] - derTempGasFormula2;
		}
	);
	
	const double H_granulation = 0; // Granulation enthaly
	// Currently unsed, preparation for later

	const double deltaHParticleGrowth = 0; // Sum of all variables concering particle growth
	// Pre. for granulation
	const double sumEnergyStreamsParticle = Q_AP3 - Q_PW - Q_PF2 - deltaHParticleGrowth + (H_susp2 - HFlowVapor2) * (1 - varPhi2); // Sum of energy streams to particle
	const double derTempParticleFormula2 = sumEnergyStreamsParticle / (mHoldupSolid * (unit->C_PParticle + (liquidFilm ? XCutt:X) * unit->C_PWaterLiquid))+(liquidFilm?0: -varTempParticle2 * unit->C_PWaterLiquid / (unit->C_PParticle + unit->C_PWaterLiquid * X) * ((mFlowSuspension - MFlowVapor) / mHoldupSolid));// Change of particel temperature
	// Dosta 2010 eq. 11 with addjustments: no sum necessary as all values are not height dependent, division not applicable as no summation occures
	if (unit->particlesGlobal)
		_res[m_iTempParticle] = _ders[m_iTempParticle] - derTempParticleFormula2;
	doubleTransfer = derTempParticleFormula2;

	const double sumEnergyStreamsFilm = Q_PF2 + Q_AF2 + (H_susp2 - HFlowVapor2) * varPhi2; // Sum of energy streams to film
	// Dosta 2010 eq. 9 - no sum necessary as all values are not height dependent, division not applicable as no summation occures

	// Change in film temperature
	double derTempFilmFormula2 = 0;
	if (liquidFilm)
	{
		derTempFilmFormula2 = 1.0 / (unit->C_PWaterLiquid * mHoldupSolid * (X - XCutt)) * sumEnergyStreamsFilm
			- varTempFilm2*unit->C_PWaterLiquid/(unit->C_PWaterLiquid* (X - XCutt)) * ((mFlowSuspensionLiquid - MFlowVapor)/mHoldupSolid);
	}
	else
		derTempFilmFormula2 = derTempParticleFormula2;
	if (unit->particlesGlobal)
		_res[m_iTempFilm] = _ders[m_iTempFilm] - derTempFilmFormula2;

	

	// Codes for checking _vars, _ders and _res
	std::vector vars2(_vars, _vars + GetVariablesNumber());
	std::vector ders2(_ders, _ders + GetVariablesNumber());
	std::vector  res2(_res, _res + GetVariablesNumber());
}

void CUnitDAEModel::ResultsHandler(double _time, double* _vars, double* _ders, void *_unit)
{
	auto* unit = static_cast<CDryerBatch*>(_unit); // crash after 16700

	std::ostringstream  os;
	if (debugToggle){
		if (_time <= 1) {
			unit->ShowInfo(std::to_string(counter));
			counter = 0;
		}
	} //_DEBUG
	if (_time / (100 * progressCounter) >= 1)
	{
		os << "SimTime: " << 100 * progressCounter << " s has passed.";
		unit->ShowInfo(os.str());
		os.str("");
		progressCounter++;
	}
		
	unit->m_holdup->AddTimePoint(_time);

	massFlow mFlowInGas = unit->m_inGasStream->GetMassFlow(_time);
	const massFlow mFlowInLiquid = unit->m_inLiquidStream->GetPhaseMassFlow(_time, EPhase::LIQUID);
	const mass mHoldup = unit->m_holdup->GetMass(0);
	const mass mHoldupSolid = unit->m_holdup->GetPhaseMass(_time, EPhase::SOLID);
	const mass mHoldupGas = unit->m_holdup->GetPhaseMass(_time, EPhase::GAS);
	const pressure mPressure = unit->m_holdup->GetPressure(_time);
	mass mHoldupLiquid = unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID);
	const area A_P = unit->A_P;
	const length Delta_f = unit->Delta_f;
	const moistureContent Y_in = unit->Y_in;

	const massFlow mFlowVapor = VaporFlowStorage;
	const temperature mTempParticle = _vars[m_iTempParticle];
	const moistureContent mY_gOut = _vars[m_iY_g + (unit->N_total - 1)];
	const double mTempGasOut = _vars[m_iTempGas + (unit->N_total - 1)];
	const temperature mTempFilm = _vars[m_iTempFilm];

	const double prevTime = unit->m_holdup->GetPreviousTimePoint(_time);
	const mass mWaterNozzle = (_time - prevTime) * mFlowInLiquid;
	const mass mWaterAirIn = (_time - prevTime) * mFlowInGas * Y_in;
	const mass mVapor = (_time - prevTime) * mFlowVapor;

	unit->m_holdup->AddStream(prevTime, _time, unit->m_inLiquidStream);
	//mHoldupLiquid = unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID);

	unit->m_holdup->RemoveTimePointsAfter(_time);

	unit->m_holdup->SetPhaseMass(_time, EPhase::LIQUID, unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID) - mVapor);
	mHoldupLiquid = unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID);
	unit->m_holdup->SetTemperature(_time, mTempParticle);

	//Set outlet gas stream: incl. dry gas & water vapor
	// !!! different from inlet gas stream which contains only dry gas
	unit->m_outExhaustGasStream->CopyFromStream(_time, unit->m_inGasStream);
	unit->m_VaporStream->SetMassFlow(_time, mFlowVapor);
	unit->m_outExhaustGasStream->AddStream(_time, unit->m_VaporStream);
	unit->m_outExhaustGasStream->AddStream(_time, unit->m_inNozzleAirStream);
	unit->m_outExhaustGasStream->SetTemperature(_time, mTempGasOut);
	unit->m_outExhaustGasStream->SetPhaseFraction(_time, EPhase::SOLID, 0);
	unit->m_outExhaustGasStream->SetPhaseFraction(_time, EPhase::LIQUID, 0);
	unit->m_outExhaustGasStream->SetPhaseFraction(_time, EPhase::VAPOR, 1);
	const massFlow mFLowGasOut = unit->m_outExhaustGasStream->GetMassFlow(_time);

	const moistureContent X = mHoldupLiquid / mHoldupSolid;

	// Height averaged temperature of the gas phase [K]
	// Dosta 2010 eq. 14a
	const double avVarTempGas = CalculateAverage(_vars, m_iTempGas, unit->N_total);
	unit->SetStateVariable("Avg. temperature gas", avVarTempGas - unit->T_ref, _time);
	// Output of average gas temperature in celsius

	// Height averaged moisture content of the gas phase [kg/kg]
	// Dosta 2010 eq. 14b
	const double avY = CalculateAverage(_vars, m_iY_g, unit->N_total);
	unit->SetStateVariable("Avg. Y", avY, _time);


	const double avRH = unit->GetRelativeHumidity(avY, avVarTempGas, mPressure);
	const double Xeq = !unit->particlesGlobal?0:unit->CalcuateSolidEquilibriumMoistureContent(_time, mTempParticle, avRH);
	const double XCutt = unit->REAinv(Xeq, 1. - unit->phiCuttOff);
	double mPhi = unit->particlesGlobal ? (mHoldupLiquid - XCutt * mHoldupSolid) / (unit->rhoWater * Delta_f * unit->A_P) : 0;

	if (mPhi < 0)
		mPhi = 0;
	if (mPhi > 1)
		unit->RaiseWarning("Attention: Overwetting in the dryer at " + std::to_string(_time) + " s");
		// Give warning if degree of wetness (phi) exceeds 1
		// Warning texts in GUI
		

	if (unit->particlesGlobal) {
		unit->SetStateVariable("X", X, _time);
		unit->SetStateVariable("w%", X/(1.0+X)*100, _time);
		unit->SetStateVariable("phi", mPhi, _time);
		unit->SetStateVariable("solid mass in holdup", mHoldupSolid, _time);
		unit->SetStateVariable("liquid mass in holdup", unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID), _time);
		unit->SetStateVariable("evaporation rate [kg/s]", mFlowVapor, _time);
		unit->SetStateVariable("Temperature film", mTempFilm - unit->T_ref, _time);
		unit->SetStateVariable("Temperature particle", mTempParticle - unit->T_ref, _time);
		unit->SetStateVariable("Moisture precentage", unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID) * 100 / (unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID) + unit->m_holdup->GetPhaseMass(_time, EPhase::SOLID)), _time);
		unit->SetStateVariable("RH particles", unit->GetRelativeHumidity(avY, avVarTempGas, mPressure) * 100, _time);
	}
	unit->SetStateVariable("T_out", mTempGasOut - unit->T_ref,_time);
	unit->SetStateVariable("Y_out", mY_gOut, _time);
	unit->SetStateVariable("RH_out [%]", unit->GetRelativeHumidity(mY_gOut, mTempGasOut, mPressure) * 100, _time);
	unit->SetStateVariable("gas mass in holdup", mHoldupGas, _time);
	size_t probeLayer = unit->DetermineLayersInSectionFilledWithBed(0, unit->heightOfChamberTemperatureProbe / unit->chamber.at(0).height)-1;
	unit->SetStateVariable("Temperature chamber", _vars[m_iTempGas + probeLayer] - unit->T_ref, _time);
	ParallelFor(unit->N_total, [&](size_t i)
		{
			unit->SetStateVariable("Temperature gas height layer " + std::to_string(i), _vars[m_iTempGas + i] - unit->T_ref, _time);
			unit->SetStateVariable("Y height layer " + std::to_string(i), _vars[m_iY_g + i], _time);
		});

if (debugToggle){
		unit->SetStateVariable("Mass balance", (_time - prevTime==0?0:( - unit->m_holdup->GetMass(_time) + unit->m_holdup->GetMass(prevTime)) / (_time - prevTime)) + ((unit->m_inGasStream->GetMassFlow(_time)+unit->m_inNozzleAirStream->GetMassFlow(_time))*(1+Y_in) - unit->m_outExhaustGasStream->GetMassFlow(_time) + unit->m_inLiquidStream->GetMassFlow(_time)), _time);
		//unit->SetStateVariable("Energy balance", unit->m_inGasStream->GetMassFlow(_time) - unit->m_outExhaustGasStream->GetMassFlow(_time),_time);
		//unit->SetStateVariable("Mass balance new", -unit->m_holdup->GetMass(_time) + unit->m_holdup->GetMass(prevTime) + (_time - prevTime) * (unit->m_inGasStream->GetMassFlow(_time)*(1+Y_in) - unit->m_outExhaustGasStream->GetMassFlow(_time) + unit->m_inSuspensionStream->GetMassFlow(_time)), _time);

		double EnergyLiquidPhase = unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID) * unit->C_PWaterLiquid * mTempFilm;
		double EnergySolidPhase = unit->m_holdup->GetPhaseMass(_time, EPhase::SOLID) * unit->C_PParticle * mTempParticle;
		double EnergyGasPhase = unit->m_holdup->GetPhaseMass(_time, EPhase::GAS) * ((unit->C_PGas + unit->C_PWaterVapor * (unit->YavgOld - avY)) * avVarTempGas/* + unit->Delta_h0 * (unit->YavgOld - Y_av)*/);

		double EnergyInputGas = (_time - prevTime) * unit->m_inGasStream->GetMassFlow(_time) * ((unit->C_PGas + unit->Y_in * unit->C_PWaterVapor) * unit->m_inGasStream->GetTemperature(_time)/* + unit->Y_in * unit->Delta_h0*/);
		double EnergyInputLiquid = (_time - prevTime) * (unit->m_inLiquidStream->GetPhaseMassFlow(_time,EPhase::LIQUID) * unit->C_PWaterLiquid + unit->m_inLiquidStream->GetPhaseMassFlow(_time, EPhase::SOLID) * unit->C_PParticle) * unit->m_inLiquidStream->GetTemperature(_time);
		double EnergyOutputGas = (_time - prevTime) * (unit->m_outExhaustGasStream->GetMassFlow(_time) * ((unit->C_PGas + mY_gOut * unit->C_PWaterVapor) * mTempGasOut/* + mY_out * unit->Delta_h0*/));

		double HeatLoss = heatLossTransfer;
		/*
		const double EnergyChangeLiquidPhase = (-unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID) * mTempFilm + unit->m_holdup->GetPhaseMass(prevTime, EPhase::LIQUID) * unit->TempLiquidOld)* unit->C_PWaterLiquid;
		const double EnergyChangeSolidPhase = unit->m_holdup->GetPhaseMass(_time, EPhase::SOLID) * unit->C_PParticle * (-mTempParticle + unit->TempSolidOld);
		const double EnergyChangeGasPhase = unit->m_holdup->GetPhaseMass(_time, EPhase::GAS) * (unit->C_PGas + unit->C_PWaterVapor * (unit->YavgOld-Y_av)) * (-varAvTempGas + unit->TempGasOld);
		const double EnergyInputGas = (_time - prevTime) * unit->m_inGasStream->GetMassFlow(_time) * ((unit->C_PGas + unit->Y_in * unit->C_PWaterVapor) * unit->m_inGasStream->GetTemperature(_time));
		const double EnergyInputLiquid = (_time - prevTime) * unit->m_inSuspensionStream->GetMassFlow(_time) * (unit->C_PWaterLiquid * unit->w_l + unit->C_PParticle * (1 - unit->w_l))* unit->m_inSuspensionStream->GetTemperature(_time);
		const double EnergyOutputGas = (_time - prevTime) * (-unit->m_outExhaustGasStream->GetMassFlow(_time) * (
			(unit->C_PGas + mY_out * unit->C_PWaterVapor) * T_out
			+ (mY_out - unit->Y_in) * unit->Delta_h0));
		const double EnergyBalance = EnergyChangeLiquidPhase + EnergyChangeSolidPhase + EnergyChangeGasPhase + EnergyInputGas + EnergyInputLiquid + EnergyOutputGas;
		unit->TempLiquidOld = mTempFilm;
		unit->TempSolidOld = mTempParticle;
		unit->TempGasOld = varAvTempGas;
		unit->YavgOld = Y_av;
		*/
		double DeltaEnergyLiquid = EnergyLiquidPhase - unit->EnergyLiquidPhaseOld;
		double DeltaEnergySolid = EnergySolidPhase - unit->EnergySolidPhaseOld;
		double DeltaEnergyGas = EnergyGasPhase - unit->EnergyGasPhaseOld;
		double DeltaGasStreams = EnergyInputGas - EnergyOutputGas;
		double DeltaEnergy = DeltaEnergyLiquid + DeltaEnergySolid + DeltaEnergyGas;
		double EnergyBalance = -DeltaEnergy + DeltaGasStreams + EnergyInputLiquid - unit->HeatLossOld;
		/*unit->EnergyLiquidPhaseOld = EnergyLiquidPhase;
		unit->EnergyGasPhaseOld = EnergyGasPhase;
		unit->EnergySolidPhaseOld = EnergySolidPhase;
		unit->HeatLossOld = HeatLoss;
		unit->SetStateVariable("Energy balance new", EnergyBalance, _time);*/
		EnergyGasPhase = unit->m_holdup->GetPhaseMass(_time, EPhase::GAS) * ((unit->C_PGas + unit->C_PWaterVapor * avY) * avVarTempGas + unit->Delta_h0 * avY);
		EnergyInputGas = unit->m_inGasStream->GetMassFlow(_time) * ((unit->C_PGas + unit->Y_in * unit->C_PWaterVapor) * unit->m_inGasStream->GetTemperature(_time)+ unit->Y_in * unit->Delta_h0);
		EnergyInputLiquid = (unit->m_inLiquidStream->GetPhaseMassFlow(_time, EPhase::LIQUID) * unit->C_PWaterLiquid + unit->m_inLiquidStream->GetPhaseMassFlow(_time, EPhase::SOLID) * unit->C_PParticle) * unit->m_inLiquidStream->GetTemperature(_time);
		EnergyOutputGas = (unit->m_outExhaustGasStream->GetMassFlow(_time) * ((unit->C_PGas + mY_gOut * unit->C_PWaterVapor) * mTempGasOut + mY_gOut * unit->Delta_h0));
		EnergyBalance = (_time - prevTime == 0 ? 0 : -(EnergyGasPhase + EnergyLiquidPhase + EnergySolidPhase) + (unit->EnergyGasPhaseOld + unit->EnergyLiquidPhaseOld + unit->EnergySolidPhaseOld) / (_time - prevTime)) + EnergyInputGas+EnergyInputLiquid-EnergyOutputGas-unit->HeatLossOld;
		unit->SetStateVariable("Energy balance", EnergyBalance, _time);
		unit->EnergyLiquidPhaseOld = EnergyLiquidPhase;
		unit->EnergyGasPhaseOld = EnergyGasPhase;
		unit->EnergySolidPhaseOld = EnergySolidPhase;

		//unit->SetStateVariable("debug0", EnergyBalance!=0?mPhi/(-EnergyBalance):0, _time);// ~mPhi ~X
		/*
		unit->SetStateVariable("debug1",unit->GetGasSaturationMoistureContent(varAvTempGas, mPressure)/unit->Y_sat,_time);
	

		//unit->SetStateVariable("debug15", (mPressure - CachePressureRange <= inputCache[1][1] && inputCache[1][1] <= mPressure + CachePressureRange)?1:0, _time);
		//unit->SetStateVariable("debug17",resultCache[0],_time);
		unit->SetStateVariable("debug18", unit->GetRelativeHumidity(Y_av + Y_in, varAvTempGas, mPressure),_time);
		double RH = unit->GetRelativeHumidity(Y_av, varAvTempGas, mPressure);
		unit->SetStateVariable("debug19", RH, _time);
		double Xeq = unit->CalcuateSolidEquilibriumMoistureContent(_time, mTempParticle, RH);
		double deltaX = unit->m_holdup->GetPhaseMass(_time, EPhase::LIQUID) / mHoldupSolid - Xeq;
		unit->SetStateVariable("debug20", deltaX, _time);
		double REA = unit->REA(deltaX);
		unit->SetStateVariable("debug21", REA, _time);
		double curve = 1 - REA;
		unit->SetStateVariable("debug22", curve, _time);
		unit->SetStateVariable("debug24", Xeq-XeqTransfer, _time);
		unit->SetStateVariable("debug25", X-curveTransfer, _time);
		if (_time >= 134)
		{
			std::vector vars(_vars, _vars + GetVariablesNumber());
			std::vector ders(_ders, _ders + GetVariablesNumber());
			//std::vector  res(_res, _res + GetVariablesNumber());
			double deltaTime = _time - unit->m_holdup->GetPreviousTimePoint(_time);
			std::vector<double> DebugVector;
			DebugVector.reserve(vars.size() + ders.size()+3); // preallocate memory
			DebugVector.insert(DebugVector.end(), vars.begin(), vars.end());
			DebugVector.insert(DebugVector.end(), ders.begin(), ders.end());
			DebugVector.push_back(deltaX);
			DebugVector.push_back(MFlowVaporLiquid);
			DebugVector.push_back(mFlowVapor);
			//DebugCache.push_back(DebugVector);
			if (_time >= 136)
				double test = 1;
		}
		*/
		/*
		std::vector ders(_ders, _ders + GetVariablesNumber());
		std::pair<double,size_t> output = std::make_pair(0,-1);
		for (size_t i = 0; i < ders.size();i++ )
			if (abs(ders[i])>output.first)
				output = std::make_pair(ders[i], i);
		std::ofstream myfile;
		std::string fileLocation = "E:\\Dyssol\\output.csv";
		myfile.open(fileLocation, std::fstream::app);
		myfile << _time << "," << output.first << "," << output.second << "\n";
		myfile.close();*/
		unit->SetStateVariable("debug0", progressCounterTotal, _time);
		unit->SetStateVariable("debug1", liquidSideLimitedGlobal, _time);
		unit->SetStateVariable("debug2", TransferBool, _time);
		unit->SetStateVariable("debug3", doubleTransfer, _time);
		double X_eq = !unit->particlesGlobal?0:unit->CalcuateSolidEquilibriumMoistureContent(_time, mTempParticle, unit->GetRelativeHumidity(avY, avVarTempGas));
		unit->SetStateVariable("debug4",X_eq , _time);
		unit->SetStateVariable("debug5", X-X_eq, _time);
		unit->SetStateVariable("debug6", HeatLoss, _time);
		unit->SetStateVariable("debug7", (mTempFilm - unit->T_ref) - (mTempFilm - unit->T_ref) * TransferBool, _time);
		if (_time >= 99)
			bool breakPoint = true;
		//unit->SetStateVariable("debug3", RHTransfer, _time);
		//unit->SetStateVariable("debug3", unit->GetRelativeHumidity(Y_av, varAvTempGas, mPressure), _time);
		std::vector vars(_vars, _vars + GetVariablesNumber());
		std::vector ders(_ders, _ders + GetVariablesNumber());
	} // DEBUG
		//Temp particle climbs while temps of gas and film fall
		/*
		std::pair<double, std::vector<double>> result = CalculateTotalHeatLoss(_time, _unit,X,_vars, true);
		unit->ShowInfo("Heat loss - particles");
		unit->ShowInfo(std::to_string(result.first));
		unit->ShowInfo(std::to_string(result.second[0]));
		result = CalculateTotalHeatLoss(_time, _unit, X, _vars, false);
		unit->ShowInfo("Heat loss - gas");
		unit->ShowInfo(std::to_string(result.first));
		unit->ShowInfo(std::to_string(result.second[0]));
		unit->ShowInfo("");
		unit->ShowInfo(std::to_string(XeqTransfer));
		unit->ShowInfo("");*/
		//unit->ShowInfo(std::to_string(CalculateAlpha_PW(mTempParticle, X, varAvTempGas, mPressure, _unit, _time)* unit->AreaGasParticleWallHeatLoss* (mTempParticle - unit->TemperatureWall)));
		//unit->ShowInfo(std::to_string(unit->CalckA(CalculateAlpha_PW(mTempParticle, X, varAvTempGas, mPressure, _unit, _time), 0, 0.35, { 0.225 ,0.230 }, { 15 }) * (mTempParticle - unit->T_inf)));
		//unit->ShowInfo(std::to_string(unit->CalcAlphaOutside(_time, h, d.front(), /*unit->TemperatureWall*/0.4 * (std::max(unit->m_inGasStream->GetTemperature(_time), unit->m_holdup->GetTemperature(_time)) + unit->T_inf), unit->T_inf)));
}

/// Supporting functions ///
double CUnitDAEModel::CalculateAverage(double* _vars, size_t variableKey, int64_t end, int64_t start) const
{
	double tempAvg = 0;
	for (int64_t i = start; i < end; i++)
		tempAvg += _vars[variableKey + i];
	if (end - start != 0)
		tempAvg = tempAvg / (end - start);
	return tempAvg;
}

/// Heat loss calculation functions ///
std::pair<double, std::vector<double>> CUnitDAEModel::CalculateChamberHeatLoss(double _time, void* _unit, double* _vars)
{
	auto* unit = static_cast<CDryerBatch*>(_unit);
	const double sectionsFilledWithBed = unit->DetermineSectionsFilledWithBed(_time, _vars[m_iTempParticle]);
	//	const double sectionsFilledWithBed = unit->CalculateBedHeightOrDetermineSectionsFilledWithBed(_time, _vars[m_iTempParticle],false);
	double heightUsage, fullSections;
	heightUsage = modf(sectionsFilledWithBed, &fullSections);

	double Q_PW = 0;
	std::vector<double> Q_GW;

	// Fully filled sections
	for (int section = 0; section < fullSections; section++)
	{
		std::vector<double> Q_GWsection(unit->chamber.at(section).layers, 0);
		Q_PW += CalculateSectionHeatLoss(_time, _unit, _vars, section, &Q_GWsection);
		Q_GW.insert(Q_GW.end(),Q_GWsection.begin(), Q_GWsection.end());
	}

	// Partilly filled section

	std::vector<double> Q_GWsection(unit->chamber.at(fullSections).layers, 0);
	Q_PW += CalculateSectionHeatLoss(_time, _unit, _vars,fullSections, &Q_GWsection,heightUsage);
	Q_GW.insert(Q_GW.end(), Q_GWsection.begin(), Q_GWsection.end());

	// Fully empty sections
	for (size_t section = std::max(std::ceil(sectionsFilledWithBed),1.); section < unit->chamber.size(); section++)
	{
		std::vector<double> Q_GWsection(unit->chamber.at(section).layers, 0);
		Q_PW += CalculateSectionHeatLoss(_time, _unit, _vars, section, &Q_GWsection, 0);
		Q_GW.insert(Q_GW.end(), Q_GWsection.begin(), Q_GWsection.end());
	}
		
	// Top plate
	Q_GW.back() += CalculateTopPlateHeatLoss(_time,_unit,_vars);
	return std::make_pair(Q_PW, Q_GW);
}

double CUnitDAEModel::CalculateSectionHeatLoss(double _time, void* _unit, double* _vars, size_t section, std::vector<double>* Q_GW, double heightUsage)
{
	auto* unit = static_cast<CDryerBatch*>(_unit);

	bool noTemperatureDelta = unit->particlesGlobal ? ( (_vars[m_iTempParticle] - unit->T_inf) != 0 ? false : true) : true;
	if (noTemperatureDelta)
		for (int i = 0; i < Q_GW->size(); i++)
			if (_vars[m_iTempGas + (i + section * Q_GW->size())] - unit->T_inf != 0)
			{
				noTemperatureDelta = false;
				break;
			}
	if (noTemperatureDelta)
		return 0;

	const double alphaGPipe = unit->CalculateAlpha_GW(_time, section);
	const double alphaGFB = unit->CalculateAlpha_GW(_time);
	size_t layersOfPrevSections = 0;
	for (size_t i = 0; i < section; i++)
	{
		layersOfPrevSections += unit->chamber.at(i).layers;
	}
	const double varAvTempGas = CalculateAverage(_vars, m_iTempGas, unit->chamber.at(section).layers+layersOfPrevSections, layersOfPrevSections);
	const double alphaPFB = !unit->particlesGlobal ? 0 : unit->CalculateAlpha_PW(_vars[m_iTempParticle], varAvTempGas, unit->m_holdup->GetPressure(_time), _time);

	double sumWallThicknesses = 0;
	for (int i = 0; i < unit->chamber.at(section).wallThicknesses.size(); i++)
		sumWallThicknesses += unit->chamber.at(section).wallThicknesses.at(i);
	double outerDiameter = sumWallThicknesses + (unit->chamber.at(section).dimensionsInternal.at(0).first + unit->chamber.at(section).dimensionsInternal.at(0).second) / 2;
	double Q_PW = 0;

	double alphaOut = 0;
	double kA_P = 0;
	double kA_G = 0;

	size_t usedLayers = unit->DetermineLayersInSectionFilledWithBed(section, heightUsage);
	double usedLayerPercentage = 1. / unit->chamber.at(section).layers * usedLayers;
	size_t layersOfPreveusSections = 0;
	for (size_t i = 0; i < section; i++)
		layersOfPreveusSections += unit->chamber.at(i).layers;

	double temperatureWall0 = 0;
	double temperatureWall = 0;

	if (usedLayerPercentage > 0)
	{
		temperatureWall0 = 0;
		temperatureWall = 0.5 * (std::max(unit->m_inGasStream->GetTemperature(_time), unit->m_holdup->GetTemperature(_time)) - unit->T_inf) + unit->T_inf;
		while (abs(temperatureWall - temperatureWall0) > 1/*temperatureWall * GetRTol() * 10 + GetATol(unit->N_total)*/)
		{
			alphaOut = unit->CalcAlphaOutside(_time, unit->chamber.at(section).height, outerDiameter, temperatureWall);

			kA_P = unit->CalculateOverallHeatTransferCoefficientCylinder(section, alphaPFB, alphaOut, usedLayerPercentage);
			Q_PW = unit->particlesGlobal ? kA_P * (_vars[m_iTempParticle] - unit->T_inf) : 0;

			kA_G = unit->CalculateOverallHeatTransferCoefficientCylinder(section, alphaGFB, alphaOut, usedLayerPercentage);
			for (int i = 0; i < usedLayers; i++)
				Q_GW->at(i) = kA_G / usedLayers * (_vars[m_iTempGas + (i + layersOfPreveusSections)] - unit->T_inf);

			double Q = Q_PW;
			for (int i = 0; i < usedLayers; i++)
				Q += Q_GW->at(i);

			temperatureWall0 = temperatureWall;
			double alphaAout = unit->chamber.at(section).height * outerDiameter * MATH_PI * alphaOut;
			temperatureWall = Q / alphaAout + unit->T_inf;
		}
	}
	temperatureWall0 = 0;
	temperatureWall = 0.5 * (std::max(unit->m_inGasStream->GetTemperature(_time), unit->m_holdup->GetTemperature(_time)) - unit->T_inf) + unit->T_inf;
	while (abs(temperatureWall - temperatureWall0) > /*temperatureWall * GetRTol() * 10 + GetATol(unit->N_total)*/ 1)
	{
		alphaOut = unit->CalcAlphaOutside(_time, unit->chamber.at(section).height, outerDiameter, temperatureWall);

		kA_G = unit->CalculateOverallHeatTransferCoefficientCylinder(section, alphaGPipe, alphaOut, -usedLayerPercentage);
		for (size_t i = usedLayers; i < unit->chamber.at(section).layers; i++)
			Q_GW->at(i) = kA_G / (unit->chamber.at(section).layers -usedLayers) * (_vars[m_iTempGas + (i + layersOfPreveusSections)] - unit->T_inf);

		double Q = 0;
		for (size_t i = usedLayers; i < unit->chamber.at(section).layers; i++)
			Q += Q_GW->at(i);

		temperatureWall0 = temperatureWall;
		double alphaAout = unit->chamber.at(section).height * outerDiameter * MATH_PI * alphaOut;
		temperatureWall = Q / alphaAout + unit->T_inf;
	}

	return Q_PW;
}

double CUnitDAEModel::CalculateTopPlateHeatLoss(double _time, void* _unit, double* _vars)
{
	auto* unit = static_cast<CDryerBatch*>(_unit);

	if (_vars[m_iTempGas + (unit->N_total - 1)] - unit->T_inf == 0)
		return 0;

	const double alphaG = unit->CalculateAlpha_GW(_time, unit->chamber.size() - 1);
	double temperatureWall0 = 0;
	double temperatureWall = 0.5 * (std::max(unit->m_inGasStream->GetTemperature(_time), unit->m_holdup->GetTemperature(_time)) - unit->T_inf) + unit->T_inf;
	double Q_GWTop = 0;
	double alphaOut = 0;
	double kA_G = 0;

	while (abs(temperatureWall - temperatureWall0) > 1 /*temperatureWall * GetRTol()*10 + GetATol(unit->N_total)*/)
	{
		alphaOut = unit->CalcAlphaOutside(_time, unit->chamber.back().dimensionsInternal.at(0).second, 0, temperatureWall, EShape::RECTANGULAR);

		kA_G = unit->CalculateOverallHeatTransferCoefficientTopPlate(alphaG, alphaOut);
		Q_GWTop = kA_G  * (_vars[m_iTempGas + (unit->N_total-1)] - unit->T_inf);

		temperatureWall0 = temperatureWall;
		double alphaAout = pow(unit->chamber.back().dimensionsInternal.at(0).second,2)/4 * MATH_PI * alphaOut;
		temperatureWall = Q_GWTop / alphaAout + unit->T_inf;
	}
	return Q_GWTop;
}