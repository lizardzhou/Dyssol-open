/* Copyright (c) 2020, Dyssol Development Team. All rights reserved. This file is part of Dyssol. See LICENSE file for license information. */

#define DLL_EXPORT
#include "DryerSimpleBatch.h"

extern "C" DECLDIR CBaseUnit* DYSSOL_CREATE_MODEL_FUN()
{
	return new CDryerSimpleBatch();
}

//////////////////////////////////////////////////////////////////////////
/// Batch dryer with a water inlet and a fluidization gas inlet, and an outlet for exhaust gas
/// - Base on Simplified Model proposed by Peglow in 202205
/// - Without heat transfer
/// - Mass transfer of water from particle to fluidization gas
/// - Mass transfer is treated as ideal plug flow reactor and heat exchanger with phase change -> consider the efficiency during mass transfer
/// - Water mass balance holds: in + holdup = out

void CDryerSimpleBatch::CreateBasicInfo()
{
	/// Set basic unit info ///
	SetUnitName("Dryer (batch without heat transfer)");
	SetAuthorName("Xiye Zhou");
	SetUniqueID("003B3DA4-702F-4647-BD90-7F9A5FDD5E1C");
}

void CDryerSimpleBatch::CreateStructure()
{
	/// Add ports ///
	AddPort("InletWater", EUnitPort::INPUT);
	AddPort("InletFluidizationGas", EUnitPort::INPUT);
	AddPort("OutletExhaustGas", EUnitPort::OUTPUT);

	/// Add unit parameters ///
	AddConstRealParameter("Delta_f", 40, "mum", "Thickness of the water film on particles", 1, 100);
	AddConstRealParameter("Y_in", 5, "g/kg dry gas", "Moisture content of the fluidization gas", 0, 40);
	AddConstRealParameter("Y_sat", 20, "g/kg dry gas", "Saturation moisture content of the fluidization gas", 10, 100);
	AddConstRealParameter("A_P", 4, "m2", "Surface area of all particles in the granulator", 1, 100);
	// Bed properties, use for development of further complexer models
	AddStringParameter("Bed properties", "", "Currently not in use, designed for further development");
	AddConstRealParameter("d_bed", 0.2, "m", "Bed diameter", 0.1, 1);
	AddConstRealParameter("H_bed", 0.106, "m", "Bed height without fluidization", 0.01, 1);
	AddConstRealParameter("eps_0", 0.4, "-", "Bed porosity without fluidization", 0, 1);
	// Tolerance for solver
	AddStringParameter("Tolerance for solver", "", "");
	AddConstRealParameter("Relative tolerance", 0.0, "-", "Solver relative tolerance. Set to 0 to use flowsheet-wide value", 0);
	AddConstRealParameter("Absolute tolerance", 0.0, "-", "Solver absolute tolerance. Set to 0 to use flowsheet-wide value", 0);


	/// Add holdups ///
	AddHoldup("Holdup");

	/// Set this unit as user data of model ///
	m_model.SetUserData(this);
}

void CDryerSimpleBatch::Initialize(double _time)
{
	/// Check flowsheet parameters ///
	if (!IsPhaseDefined(EPhase::SOLID))		RaiseError("Solid phase has not been defined.");
	if (!IsPhaseDefined(EPhase::LIQUID))	RaiseError("Liquid phase has not been defined.");
	if (!IsPhaseDefined(EPhase::VAPOR))		RaiseError("Gas phase has not been defined.");
	if (!IsDistributionDefined(DISTR_SIZE))	RaiseError("Size distribution has not been defined.");

	/// Get holdup ///
	m_holdup = GetHoldup("Holdup");

	/// Get pointers to streams ///
	m_inWaterStream = GetPortStream("InletWater");
	m_inGasStream = GetPortStream("InletFluidizationGas");
	m_outExhaustGasStream = GetPortStream("OutletExhaustGas");

	/// Read input parameters ///
	const double Y_in = GetConstRealParameterValue("Y_in");
	
	/// Information of PSD, currently not necessary ///
	/*
	/// Get number of classes for PSD ///
	m_classesNum = GetClassesNumber(DISTR_SIZE); //n
	/// Get grid of PSD ///
	m_sizeGrid = GetNumericGrid(DISTR_SIZE); //d_min
	m_avgDiam = GetClassesMeans(DISTR_SIZE); //d_m,i
	m_classSize = GetClassesSizes(DISTR_SIZE); //Delta d
	
	/// Get initial PSD ///
	std::vector<double> vPSD = m_holdup->GetPSD(_time, PSD_q3);
	
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

	/// Clear all state variables in model ///
	m_model.ClearVariables();

	/// Add state variables to a model ///
	m_model.m_iX = m_model.AddDAEVariable(true, 0, 0, 0);	// Moisture content in the bed in [kg/kg]
	m_model.m_iMFlowVapor = m_model.AddDAEVariable(false, 0, 0, 0);	// Mass flow of water vapor in the gas going to particles
	m_model.m_iEtaPFTR = m_model.AddDAEVariable(false, 0, 0, 0);	// Efficiency of water mass transfer in the gas
	m_model.m_iNTU = m_model.AddDAEVariable(false, 0, 0, 0);	// Number of transfer units
	m_model.m_iPhi = m_model.AddDAEVariable(false, 0, 0, 0);	// Degree of wetness of the particles
	m_model.m_iMWaterBed = m_model.AddDAEVariable(false, 0, 0, 0);	// Water mass in the bed
	m_model.m_iYOut = m_model.AddDAEVariable(false, Y_in, 0, 0); // Moisture content of exhaust gas in [kg/kg]

	AddStateVariable("X", 0); 
	AddStateVariable("dX/dt", 0); 
	AddStateVariable("MFlowVapor", 0); 
	AddStateVariable("phi", 0); 
	AddStateVariable("Y_out", Y_in); 
	AddStateVariable("Eta_PFTR", 0);	
	
	/// Set tolerances to model ///
	const auto rtol = GetConstRealParameterValue("Relative tolerance");
	const auto atol = GetConstRealParameterValue("Absolute tolerance");
	m_model.SetTolerance(rtol != 0.0 ? rtol : GetRelTolerance(), atol != 0.0 ? atol : GetAbsTolerance());

	/// Set model to a solver ///
	if (!m_solver.SetModel(&m_model))
		RaiseError(m_solver.GetError());
}

void CDryerSimpleBatch::SaveState()
{
	m_solver.SaveState();
}

void CDryerSimpleBatch::LoadState()
{
	m_solver.LoadState();
}

void CDryerSimpleBatch::Simulate(double _timeBeg, double _timeEnd)
{
	if (!m_solver.Calculate(_timeBeg, _timeEnd))
		RaiseError(m_solver.GetError());
}

void CUnitDAEModel::CalculateResiduals(double _time, double* _vars, double* _ders, double* _res, void* _unit)
{
	const auto* unit = static_cast<CDryerSimpleBatch*>(_unit);
	const CStream* inWaterStream = unit->GetPortStream("InletWater");
	const CStream* inGasStream = unit->GetPortStream("InletFluidizationGas");
	const CHoldup* holdup = unit->GetHoldup("Holdup");


	/// Read input parameters /// 
	const double mFlowWater = inWaterStream->GetMassFlow(_time);
	const double mFlowGas = inGasStream->GetMassFlow(_time);
	const double mHoldup = holdup->GetMass(_time);
	const double Y_in = unit->GetConstRealParameterValue("Y_in") / 1000; // convert in [kg/kg]
	const double Y_sat = unit->GetConstRealParameterValue("Y_sat") / 1000; // convert in [kg/kg]
	const double A_P = unit->GetConstRealParameterValue("A_P");
	const double Delta_f = unit->GetConstRealParameterValue("Delta_f") * 1e-6; // convert in [m]

	/// DAE system ///
	// Values in _vars: should not be changed, must set to const !!!
	const double varMFlowVapor = _vars[m_iMFlowVapor]; 
	const double varEtaPFTR = _vars[m_iEtaPFTR];
	const double varNTU = _vars[m_iNTU];
	const double varPhi = _vars[m_iPhi];
	const double varMWaterBed = _vars[m_iMWaterBed];
	const double varYOut = _vars[m_iYOut];
	// Moisture content in the bed
	const double derX = (mFlowWater - varMFlowVapor) / mHoldup;
	_res[m_iX] = _ders[m_iX] - derX;
	// Mass flow of water vapor in the gas going to particles
	const double varMFlowVaporFormula = varEtaPFTR * mFlowGas * (Y_sat - Y_in);
	_res[m_iMFlowVapor] = varMFlowVapor - varMFlowVaporFormula;
	// Efficiency of water mass transfer of the gas
	const double varEtaPFTRFormula = 1 - exp(-varNTU);
	_res[m_iEtaPFTR] = varEtaPFTR - varEtaPFTRFormula;
	// Number of transfer units of water mass transfer of the gas
	const double varNTUFormula = rhoGas * beta_GP * A_P * varPhi / mFlowGas;
	_res[m_iNTU] = varNTU - varNTUFormula;
	// Degree of wetness of the particle
	const double varPhiFormula = ((varMWaterBed / rhoWater) / Delta_f) / A_P;
	_res[m_iPhi] = varPhi - varPhiFormula;
	// Water mass in the bed
	const double varMWaterBedFormula = mHoldup * _vars[m_iX];
	_res[m_iMWaterBed] = varMWaterBed - varMWaterBedFormula;
	// Moisture content of exhaust gas [in g/kg dry gas]
	const double varYOutFormula = varMFlowVapor / mFlowGas + Y_in;
	_res[m_iYOut] = varYOut - varYOutFormula;

	// Codes for checking _vars, _ders and _res, only turn on if needed
	///*
	#ifdef _DEBUG
		std::vector vars(_vars, _vars + GetVariablesNumber());
		std::vector ders(_ders, _ders + GetVariablesNumber());
		std::vector  res(_res, _res + GetVariablesNumber());
	#endif
	//*/
}

void CUnitDAEModel::ResultsHandler(double _time, double* _vars, double* _ders, void *_unit)
{
	auto* unit = static_cast<CDryerSimpleBatch*>(_unit);

	unit->m_holdup->AddTimePoint(_time);
	double mHoldup = unit->m_holdup->GetMass(_time);

	// calculate and apply transformation matrix: currently not necessary
	/*
	//const std::vector<double> inDistr = unit->m_holdup->GetDistribution(_time, DISTR_SIZE); //Mass fraction
	//const std::vector<double> inDistr_q3 = unit->m_holdup->GetPSD(_time, PSD_q3);
	//const std::vector<double> outDistr = Convertq3ToMassFractions(unit->m_sizeGrid, temp_q3);
	//CTransformMatrix TM;
	//CMyGranulatorWithoutThermo::CalculateTM(DISTR_SIZE, inDistr, outDistr, TM); // here only use mass fraction for inDistr and outDistr
	//unit->m_holdup->ApplyTM(_time, TM);
	*/

	// HACK: handle holdup as material stream
	static_cast<CBaseStream*>(unit->m_holdup)->Add(_time, *static_cast<CBaseStream*>(unit->m_inWaterStream));	

	unit->m_holdup->RemoveTimePointsAfter(_time);
	unit->m_holdup->SetMass(_time, mHoldup);

	// Give warning if degree of wetness (phi) exceeds 1
	const double A_P = unit->GetConstRealParameterValue("A_P");
	const double Delta_f = unit->GetConstRealParameterValue("Delta_f") * 1e-6; // convert in [m]
	const double phi = _vars[m_iPhi];
	const double X = _vars[m_iX];
	const double MWater_bed = _vars[m_iMWaterBed];
	const double MWater_bed_max = rhoWater * A_P * Delta_f;
	const double X_max = MWater_bed_max / mHoldup;
	if (phi > 1 || X > X_max || MWater_bed > MWater_bed_max)
	{
		unit->RaiseWarning("Attention: Overwetting in the dryer at " + std::to_string(_time) + " s");
	}

	//Set outlet gas stream: incl. dry gas & water vapor
	// !!! different from inlet gas stream which contains only dry gas
	const double Y_in = unit->GetConstRealParameterValue("Y_in") / 1000; // convert in [kg/kg]
	const double mFlowGas = unit->m_inGasStream->GetMassFlow(_time);
	unit->m_outExhaustGasStream->SetMassFlow(_time, mFlowGas * (1 + Y_in));
	//Assumuption: all outlet streams are in gas phase
	unit->m_outExhaustGasStream->SetPhaseFraction(_time, EPhase::SOLID, 0);
	unit->m_outExhaustGasStream->SetPhaseFraction(_time, EPhase::LIQUID, 0);
	unit->m_outExhaustGasStream->SetPhaseFraction(_time, EPhase::VAPOR, 1);

	unit->SetStateVariable("X" , _vars[m_iX] , _time);
	//unit->SetStateVariable("dX/dt", _ders[m_iX], _time);
	//unit->SetStateVariable("MFlowVapor", _ders[m_iMFlowVapor], _time);
	unit->SetStateVariable("phi" , _vars[m_iPhi] , _time);
	unit->SetStateVariable("Y_out", _vars[m_iYOut], _time);
	unit->SetStateVariable("Eta_PFTR", _vars[m_iEtaPFTR], _time);
	
}
