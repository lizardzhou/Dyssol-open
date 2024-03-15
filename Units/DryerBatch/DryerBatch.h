/* Copyright (c) 2020, Dyssol Development Team. All rights reserved. This file is part of Dyssol. See LICENSE file for license information. */

#pragma once
#include <math.h>

#include "UnitDevelopmentDefines.h"
#include "unitDeclarations.h"

class CUnitDAEModel : public CDAEModel
{
private:
	//double RHTransfer = 0;
	//double doubleTransfer = 0;
	//double MFlowVaporLiquid = 0;
	//double curveTransfer = 0;
	int64_t progressCounter = 1;
	//int progressCounterTens = 1;
	//int progressCounterTotal = 0;
	//bool liquidSideLimitedGlobal = false;
	//bool TransferBool = false;
	//double heatLossTransfer = 0;
	//std::vector<double> TransferVector;

public:
	// Internal unit for property calculation ///
	void* m_unit{};
	//bool Y_eq = true;

	/// Indices of state variables for DAE solver: 5 DAE variables ///
	// gas phase
	size_t m_iYOutGas{}; //0 - outlet gas, no height discretization
	size_t m_iTempOutGas{}; //1 - outlet gas, no height discretization
	// particle (solid) phase
	size_t m_iTempParticle{}; //2
	size_t m_iPhi{}; //3 - particle wetness degree
	//size_t m_iX{}; //particle moisture content
	// liquid phase (water film)
	size_t m_iTempFilm{}; //4 - water film (on particle surface) temperature
	// water vapor
	//size_t m_iMFlowVapor{}; //vapor flow (evaporation) rate
	//size_t m_iHFlowVapor{}; //vapor flow enthalpy
	// heat transfer
	//size_t m_iQFlow_GF{}; //heat transfer from air to water film, == Q_AP
	//size_t m_iQFlow_GP{}; //heat transfer from air to particle, == Q_AF
	//size_t m_iQFlow_PF{}; //heat transfer from particle to water film
	//size_t m_iQFlow_WE{}; // heat transfer from wall to environment (heat loss to ambient)
	// transfer coefficients
	//size_t m_iPr{}; //Prandtl number
	//size_t m_iDa{}; //water diffusion coefficient from liquid to gas

	// Debug
	//std::vector<double> derFormulaStorage; // Storage for debug purposses
	//double VaporFlowStorage = 0;
	//int counter = 0;
	//std::vector<double> QAP3;
	//double CacheTemperatureRange = 5;
	//double CachePressureRange = 100;
	
public:
	// Basics for DAESolver
	void CalculateResiduals(double _time, double* _vars, double* _ders, double* _res, void* _unit) override;
	void ResultsHandler(double _time, double* _vars, double* _ders, void* _unit) override;

	/**	Initialize the variables for debugging
	 *	\return None*/
	//void InitCounterVariables() 
	//{ 
	//	progressCounter = 1; 
	//	progressCounterTotal = 0;
	//};

	/**	Calculates average for the variable with the variable key between the start and end index
	 *	\param variables
	 *  \param key for variable
	 *  \param end index for variable of key
	 *  \param start index for variable of key
	 *	\return average of variable in unit of variable*/
	//double CalculateAverage(double* _vars, size_t variableKey, int64_t end, int64_t start =0) const;	
	/**	 Calculates heat loss over the entire chamber split up by layer
	 *	\param 
	 *	\return pair of heat loss through particles and vector for heat loss through gas all in watts*/
	//std::pair<double, std::vector<double>> CalculateChamberHeatLoss(double _time, void* _unit, double* _vars);
	//	Calculates heat loss of each chamber section
	//double CalculateSectionHeatLoss(double _time, void* _unit, double* _vars, size_t section, std::vector<double>* Q_GW, double heightUsage = 1);
	//	Calculates heat loss of top section
	//double CalculateTopPlateHeatLoss(double _time, void* _unit, double* _vars);

	bool debugToggle = false;
};

class CDryerBatch : public CDynamicUnit
{
private:
	CUnitDAEModel m_model;
	CDAESolver m_solver;

private:
	//std::vector<std::pair< EPhase, int>> CompoundsKeyIndexPhasePartnerIndex; // Storage vector for phase and phase change partner, index same as compoundKeys variable
	//void PullCompoundDataFromDatabase(double _time); // Reads all material properties using matieral database values
	//void CheckHeightDiscretizationLayers(double _time); // Adjusts number of height discretization layers if layer height is lower than max particle size
	//double minMoistureContent = 0;
	//double moistureScaler = 1;
	//massTransferCoefficient beta;// = 0.02; // Water mass transfer coefficient from gas to particle in [m/s]
	//double alpha_GF;
	//double alpha_GP; // == alpha_GF
	const double f_alpha = 1; // ratio alpha_PF / alpha_AP
	
	//double alpha_PF;
	//massTransferCoefficient beta_AF;
	//massTransferCoefficient beta_PF;
	//equilibriumMoistureContentData eqData;
	//std::set<double> RHs;
	//std::set<double> temperatures;
	//std::map< std::pair<double, double>, double> equilibriumMoistureContents;

public:
	//bool debugToggle = false;
	const temperature T_ref = STANDARD_CONDITION_T - 25; // Ref. temperature for enthalpy [K] - default 273.15 K
	temperature T_inf;// = T_ref + 20.5; // Ambient temperature [K] - default: Standard condition
	//mass mTotHoldup; // mass of solid + liquid in holdup, == user input
	// Gas phase
		density rhoGas = 1.2; // Density gas [kg/m^3] - default: air
		density rhoVapor = 0.8; // Density water vapor [kg/m3]
		dynamicViscosity etaGas = 1.8e-5; // Dynamic viscosity gas [Pa*s] - default: air
		heatCapacity C_PGas = 1200; // Heat capacity gas [J/(kg*K)] - default: air
		thermalConductivity lambdaGas = 0.025; // Thermal conductivity gas [W/(m*K)] - default: air
		molarMass molarMassGas = 0.028949; // Molar mass of gas mixture [kg/mol] - default: air
	
	// Inlet fluidization gas
		//massFlow mFlowInGas;
		//massFlow mFlowInGasDry;
		//moistureContent Y_inGas; // Moisture content of input gas stream [kg/kg]
		//double RH_inGas;
		//specificLatentHeat h_inGas; // enthalpy for inlet gas: determined by user input, in [J/kg]
		//temperature theta_inGas; 
	// Inlet nozzle gas
		//massFlow mFlowInNozzleGas;
		//massFlow mFlowInNozzleGasDry;
		//moistureContent Y_nozzle;
		//specificLatentHeat h_inNozzle;
		//temperature thetaNozzleGas;
		//specificLatentHeat h_nozzleGas;
	// Gas in holdup (whole plant, incl. chamber & expansion)
		mass mGasHoldup = 0.62; // mass of DRY gas in the plant (chamber + expansion part) [kg]
		//moistureContent Y_sat; // = 0.020; // Saturation moisture content of gas [kg/kg]
		
	// Liquid phase
		density rhoWater = 1000; // Density liquid [kg/m^3] - default: water
		heatCapacity C_PWaterLiquid = 4200; // Heat capacity liquid phase change compound[J / (kg * K)] - default: water
		heatCapacity C_PWaterVapor = 2000; // Heat capacity vapor phase change compound [J/(kg*K)] - default: water
		specificLatentHeat Delta_h0 = 2500e3; // Specific latent heat (evaporation heat) phase change compound at 0 degree [J/kg] - default: water
		thermalConductivity lambdaWater = 0.6; // Thermal conductivity [W/(m*K)] - default: water
		molarMass molarMassPhaseChangingLiquid = 0.018; // Molar mass of phase changing liquid [kg/mol] - default: water
	
	//double ratioMM; // = molarMassPhaseChangingLiquid / molarMassGas;
	// Liquid in holdup
		//mass mLiquidHoldup;
	// Spray liquid
		//massFlow mFlowSprayLiquid;
		//massFraction x_wSusp; 
		//temperature thetaSprayLiquid;
		//specificLatentHeat h_susp;
	// Particle (solid) phase
		//std::vector<double> Grid; // d_min
		//std::vector<double> q_3;
		//std::vector<double> avgClassDiam; // d_m,i
		//std::vector<double> classSize; // Delta d
		density rhoParticle = 1500; // Particle density (Cellets: skeletal density)
		heatCapacity C_PParticle = 1000; // Heat capacity
		thermalConductivity lambdaParticle = 0.2; // https://doi.org/10.1016/j.ijpharm.2017.10.018 MCC relative density ~= 0.7
	
		// Particle in holdup
		//mass mSolidHoldup;
		//length d32; // Sauter diameter
		//area A_P; // = 4; // total surface area of particle mass [m^2]
		//length Delta_f; // = 40e-6; // Thickness of the water film on particles [m]
		//moistureContent initX = 0;
		//drying kinetic parameters, CURRENTLY NOT IN USE
		//double k_dc; // = 3.5; // k for normalized drying curve
		//moistureContent X_cr; // = 0.025; // Critical moisture content [kg/kg]
	// Bed
		//length heightOfBed;
		//length diamOfBed;
		//double eps_0; // = 0.4;
		//double u_mf; // = 0;
	// process chamber
		//std::vector<chamberSection> chamber;
		//double heightOfChamber;
		//double heightOfChamberTemperatureProbe;// = 0.070;
		//double heightOfNozzle;
		//double heighestFlowTimepoint = 0;

	// Settings
		//bool calcBeta = GetCheckboxParameterValue("calcBeta");
		//bool calcY_sat = GetCheckboxParameterValue("calcY_sat");
		//bool calcNdc = GetCheckboxParameterValue("calcNdc");
		//size_t dryingCurveSetting = GetComboParameterValue("DryingCurve");
		//size_t dryingCurveSetting{};
		//double SmallBiotNumber = 0.1;
		//double phiCuttOff = 0.999;
	// REA function parameters
		//double REA1 = 0.96;
		//double REA2 = -19.63;
		//double REA3 = 0.73;

	//size_t suspLayer = 0;

	bool particlesGlobal = true;

	//size_t N_particle = 1; // Number of hight discretization layers of all sections containing particles
	//size_t N_total = 1;// Total number of hight discretization layers

	CHoldup* m_holdupSolid{}; // Holdup solid (particle)
	CHoldup* m_holdupLiquid{}; // Holdup liquid (liquid film on particle)
	CHoldup* m_holdupGas{}; // Holdup gas (moist air)
	CMaterialStream* m_inLiquidStream{}; // Input of water stream
	CMaterialStream* m_inNozzleAirStream{}; // Input nozzle air
	CMaterialStream* m_inGasStream{}; // Input gas (fluidization air) stream
	CMaterialStream* m_outExhaustGasStream{};	// Output of exhaust gas

	//CStream* m_VaporStream{};
	//CHoldup* m_Expander{}; //ToDo - check usage
	//CHoldup* workingHoldup{}; //ToDo - check usage

	// String keys of compounds for material database
	//std::vector<std::string> compoundKeys;
	// Indices of liquid (first) and vapor (second) form of phase changing compound in compoundsKey
	//std::pair<int,int> indicesOfVaporOfPhaseChangingCompound=std::make_pair(-1,-1);
	//bool debugHasBeenShown = false;

	//const temperature TempLiquidOld = T_ref;
	//const temperature TempGasOld = T_ref;
	//const temperature TempSolidOld = T_ref;
	//const moistureContent YavgOld = Y_inGas;
	//size_t DiffCoeff; // index of correlation for calcualting diffusion coefficient in the dropdown list

	//double EnergyLiquidPhaseOld = 0;
	//double EnergySolidPhaseOld = 0;
	//double EnergyGasPhaseOld = 0;
	//double HeatLossOld = 0;

public:
	void CreateBasicInfo() override;
	void CreateStructure() override;
	void Initialize(double _time) override;
	void SaveState() override;
	void LoadState() override;
	void Simulate(double _timeBeg, double _timeEnd) override;

/// Get compound or phase properties ///
	/**	Calculates sauter diameter of the particles
	 *	\param time stamp
	 *	\return sauter diameter in meter*/
	massFraction ConvertMoistContentToMassFrac(moistureContent Y) const;
	moistureContent ConvertMassFracToMoistContent(massFraction y) const;
	length CalculateHoldupSauter(double _time) const;
	area CalculateParticleSurfaceArea(double _time) const;
	moistureContent CalculateGasSaturationMoistureContent(temperature T_Gas, pressure pressureGas = STANDARD_CONDITION_P) const;
	double CalculateDiffusionCoefficient(double _time, double avgGasTemperature, double filmTemperature, double pressure = STANDARD_CONDITION_P) const;
	double CalculateGasRelativeHumidity(moistureContent Y, temperature temperature, pressure pressure = STANDARD_CONDITION_P); //const;
	// Calculates the ratio (Delta E_v / Delta E_v,eq) in case of REA
	//double REA(double deltaX) const 
	//{ 
	//	if (deltaX > 0) 
	//		return REA1 * exp(REA2 * pow(deltaX, REA3)); 
	//	else 
	//		return 1;
	//};
	//// Calculates the inverted REA, namely the evaporation rate in case of REA
	//double REAinv(double Xeq, double y) const 
	//{ 
	//	if (REA1 == 0 || REA2 == 0) 
	//		return 0; 
	//	else
	//		return pow(log(y / REA1) / REA2, 1. / REA3) + Xeq;
	//};
	//double GetAvgConstCompoundProperty(double _time, EPhase phase, ECompoundConstProperties  property) const;	
	//double GetAvgTPCompoundProperty(double _time, EPhase phase, ECompoundTPProperties property, double temperature, double pressure = STANDARD_CONDITION_P) const;

/// Dimensionless numbers ///
	// Reynolds number without height discretization
	dimensionlessNumber CalculateReynolds(double _time, length d32) const;
	// Reynolds number with height discretization for each layer
	//dimensionlessNumber CalculateReynolds(double _time, size_t section) const;
	dimensionlessNumber CalculatePrandtl(temperature avgGasTemperature) const;
	dimensionlessNumber CalculateSchmidt(double D_a) const;
	dimensionlessNumber CalculateArchimedes(length d32) const;
	dimensionlessNumber CalculateNusseltSherwood(double Nu_Sh_lam, double Nu_Sh_turb) const
	{
		return 2. + sqrt(pow(Nu_Sh_lam, 2) + pow(Nu_Sh_turb, 2));
	};
	dimensionlessNumber CalculateNusseltSherwoodLam(double Re, double Pr_Sc) const
	{
		return 0.664 * pow(Pr_Sc, 1 / 3) * sqrt(Re);
	};
	dimensionlessNumber CalculateNusseltSherwoodTurb(double Re, double Pr_Sc) const
	{
		dimensionlessNumber Nu_Sc = (0.037 * pow(Re, 0.8) * Pr_Sc) / (1. + 2.443 * pow(Re, -0.1) * (pow(Pr_Sc, 2. / 3) - 1));
		if (Nu_Sc < 0)
		{
			Nu_Sc = 0;
		}
		return Nu_Sc;
		
	};
	//	Calculate Biot number
	//dimensionlessNumber CalcBiotNumber(double _time, temperature avgGasTemperature, length d32) const
	//{
	//	return (CalculateAlpha_GP(_time, avgGasTemperature, d32) * m_holdup->GetPhaseMass(_time, EPhase::SOLID) / (rhoParticle * CalculateParticleSurfaceArea(_time) * lambdaParticle));
	//};
	//bool CheckForSmallBiot(double _time) const;


/// Heat transfer coefficients ///
	double CalculateAlpha_GP(double _time, temperature avgGasTemperature, length d32) const; // == alpha_GF
	double CalculateAlpha_PF(temperature tempWater, pressure pressureHoldup, length d32) const;

/// Mass transfer coefficient ///
	massTransferCoefficient CalculateBeta(double _time, length d32, double D_a) const;

/// X_eq, CURRENTLY NOT IN USE ///
	// Returns praticle equilibirum moisture content from material database, return 0 if no entry in data base is found [kg liqudi per kg dry solid]
	//moistureContent CalcuateSolidEquilibriumMoistureContent(double _time, temperature temperature, double RH);
	// Returns normalized drying curve
	double CalculateRelativeDryingRate(moistureContent X) const;
	//moistureContent CalculateGasEquilibriumMoistureContent(temperature temperatureParticle, pressure pressureGas, double ratioMM, double RH=1) const;
	////double CalculateEquilibriumRelativeHumidity(double _time, temperature temperature, double X) const;
	//double GetEquilibriumRelativeHumidity(double temperature, double X) const;
	////	Initializes variables containing equilibirum moisture content date
	//bool InitializeMoistureContentDatabase(std::string path);
	//moistureContent GetParticleEquilibriumMoistureContent(double temperature, double RH) const;


///  Hydrodynamics ///
	double CalculateMinFluidizeVel(double _time, length d32) const;
	double CalculateGasVel(double _time, length d32) const;
	double CalculateBedPorosity(double _time, length d32, bool homogeniusFluidization = true) const;
	//double CalculateBedHeight(double _time, double particleTemperature);	
	//double CalculateGasVel(double _time, size_t section = 0) const;	


///	Setup variables containing the dimensions and information of the chamber ///
	//void SetupChamber();

	//Testing function, called during initialzation
	//void Testing();

/// Calculate heat transfer coeffients for heat loss to the environment, CURRENTLY NOT IN USE ///
	//double CalculateAlpha_PW(double _t_p, double _t_g, double _p, double _time) const;
	//double CalculateAlpha_GW(double _time, temperature avgGasTemperature) const;
	////double CalculateAlpha_GW(double _time, size_t section) const;
	//double CalcAlphaOutside(double _time, const double h, const double D, const double Ts, EShape shape = EShape::CYLINDRICAL) const;
	//double CalckAc(double alphaIn, double alphaOut, double L, std::vector<double> d/* Inner to outer diameter*/, std::vector<double> lambda) const;
	//double CalckAp(double alphaIn, double alphaOut, std::vector<double> A, std::vector<double> delta, std::vector<double> lambda) const;

/// Calculation involving height discretization, CURRENTLY NOT IN USE ///
	/**	Calculates overall heat transfer coefficient for cylindrical section.
	 *	\param section id of section in chamber object
	 *	\param alphaInternal heat transfer coefficient internal
	 *	\param alphaExternal heat transfer coefficient external
	 *	\param heightUsage use >0 for section partailly filled with particles, <0 for sections partially void of particles and =0 for sections void of particles range -1 to 1
	 *	\return overall heat transfer coefficient*/
	//double CalculateOverallHeatTransferCoefficientCylinder(size_t section, double alphaInternal, double alphaExternal, double heightUsage = 0);
	//Calculates overall heat transfer coefficient for the top plate
	//double CalculateOverallHeatTransferCoefficientTopPlate(double alphaInternal, double alphaExternal);	
	//double DetermineSectionsFilledWithBed(double _time, double particleTemperature);	
	//double GetGasMassOfLayer(double _time, size_t layer, double gasTemperature, double particleTemperature);
	//size_t DetermineLayersInSectionFilledWithBed(size_t section, double heigthUsage) const
	//{
	//	return std::ceil(chamber.at(section).layers * heigthUsage);
	//};
	//	Calculates volume of a section
	//double CalculateSectionVolume(size_t section);	
	/**	Determine into which section a layer falls
	 *	\param index of layer
	 *	\return index of section*/
	//size_t DetermineSectionForLayer(size_t layer) const;
	//std::vector<double> GetSectionGasMass(double _time, double gasTemperature, double particleTemperature);
	//double CalculateLayerVolume(size_t section, double R1, double r1)
	//{
	//	return chamber.at(section).height* MATH_PI* (R1 * R1 + R1 * r1 + r1 * r1) / (3 * chamber.at(section).layers);
	//}
	//double CalculateBedHeightOrDetermineSectionsFilledWithBed(double _time, double particleTemperature, bool outputHeight = true);

	/// Other math functions ///
	//double lerp(double lowNum, double highNum, double midNum) const
	//{
	//	return lowNum + midNum * (highNum - lowNum);
	//}
	//double ratio(double x, double y)
	//{
	//	if (y == 0)
	//	{
	//		RaiseError("Denominator is zero.");
	//	}
	//	else
	//	{
	//		return x / y;
	//	}
	//}
};