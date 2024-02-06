/* Copyright (c) 2020, Dyssol Development Team. All rights reserved. This file is part of Dyssol. See LICENSE file for license information. */

#pragma once
#include <math.h>

#include "UnitDevelopmentDefines.h"
#include "unitDeclarations.h"

class CUnitDAEModel : public CDAEModel
{
private:
	double RHTransfer = 0;
	double doubleTransfer = 0;
	double MFlowVaporLiquid = 0;
	double curveTransfer = 0;
	int progressCounter = 1;
	int progressCounterTens = 1;
	int progressCounterTotal = 0;
	bool liquidSideLimitedGlobal = false;
	bool TransferBool = false;
	double heatLossTransfer = 0;
	std::vector<double> TransferVector;

public:
	// Internal unit for property calculation ///
	void* m_unit{};
	bool Y_eq = true;

	// Dosta model DEA variables
	size_t m_iY_g{};
	size_t m_iTempGas{};
	size_t m_iTempFilm{};
	size_t m_iTempParticle{};

	// Debug
	std::vector<double> derFormulaStorage; // Storage for debug purposses
	double VaporFlowStorage = 0;
	int counter = 0;
	//std::vector<double> QAP3;
	//double CacheTemperatureRange = 5;
	//double CachePressureRange = 100;
	
public:
	// Basics for DAESolver
	void CalculateResiduals(double _time, double* _vars, double* _ders, double* _res, void* _unit) override;
	void ResultsHandler(double _time, double* _vars, double* _ders, void* _unit) override;

	/**	Initialize the variables for debugging
	 *	\return None*/
	void InitCounterVaraibles() 
	{ 
		progressCounter = 1; 
		progressCounterTotal = 0;
	};

	/**	Calculates average for the variable with the variable key between the start and end index
	 *	\param variables
	 *  \param key for variable
	 *  \param end index for variable of key
	 *  \param start index for variable of key
	 *	\return average of variable in unit of variable*/
	double CalculateAverage(double* _vars, size_t variableKey, int64_t end, int64_t start =0) const;
	
	/**	 Calculates heat loss over the entire chamber split up by layer
	 *	\param 
	 *	\return pair of heat loss through particles and vector for heat loss through gas all in watts*/
	std::pair<double, std::vector<double>> CalculateChamberHeatLoss(double _time, void* _unit, double* _vars);
	
	//	Calculates heat loss of each chamber section
	double CalculateSectionHeatLoss(double _time, void* _unit, double* _vars, size_t section, std::vector<double>* Q_GW, double heightUsage = 1);

	//	Calculates heat loss of top section
	double CalculateTopPlateHeatLoss(double _time, void* _unit, double* _vars);

	bool debugToggle = false;
};

class CDryerBatch : public CDynamicUnit
{
private:
	CUnitDAEModel m_model;
	CDAESolver m_solver;

private:
	std::vector<std::pair< EPhase, int>> CompoundsKeyIndexPhasePartnerIndex; // Storage vector for phase and phase change partner, index same as compoundKeys variable
	void PullCompoundDataFromDatabase(double _time); // Reads all material properties using matieral database values
	void CheckHeightDiscretizationLayers(double _time); // Adjusts number of height discretization layers if layer height is lower than max particle size
	double minMoistureContent = 0;
	double moistureScaler = 1;
	massTransferCoefficient beta_GP = 0.02; // Water mass transfer coefficient from gas to particle in [m/s]
	equilibriumMoistureContentData eqData;
	//std::set<double> RHs;
	//std::set<double> temperatures;
	//std::map< std::pair<double, double>, double> equilibriumMoistureContents;
public:
	bool debugToggle = false;
	const temperature T_ref = STANDARD_CONDITION_T - 25; // Ref. temperature for enthalpy [K] - default 273.15 K
	temperature T_inf;// = T_ref + 20.5; // Ambient temperature [K] - default: Standart condition
	// Gas
		density rhoGas = 1.2; // Density gas [kg/m^3] - default: air
		dynamicViscosity etaGas = 1.8e-5; // Dynamic viscosity gas [Pa*s] - default: air
		heatCapacity C_PGas = 1200; // Heat capacity gas [J/(kg*K)] - default: air
		thermalConductivity lambdaGas = 0.025; // Thermal conductivity gas [W/(m*K)] - default: air
		molarMass molarMassGas = 0.028949; // Molar mass of gas mixture [kg/mol] - default: air
		moistureContent Y_in; // Moisture content of input gas stream [kg/kg]
		moistureContent Y_sat = 0.020; // Saturation moisture content of gas [kg/kg]
	// Liquid
		density rhoWater = 1000; // Density liquid [kg/m^3] - default: water
		heatCapacity C_PWaterLiquid = 4200; // Heat capacity liquid phase change compound[J / (kg * K)] - default: water
		heatCapacity C_PWaterVapor = 2000; // Heat capacity vapor phase change compound [J/(kg*K)] - default: water
		specificLatentHeat Delta_h0 = 2500e3; // Specific latent heat (evaporation heat) phase change compound at 0 degree [J/kg] - default: water
		thermalConductivity lambdaWater = 0.6; // Thermal conductivity [W/(m*K)] - default: water
		massFraction w_l = 1; // Liquid mass fraction of suspenstion [kg/kg] - default: 1
		molarMass molarMassPhaseChangingLiquid = 0.018; // Molar mass of phase changing liquid [kg/mol] - default: water
	// Particle
		density rhoParticle = 1500; // Particle density (Cellets: skeletal density)
		heatCapacity C_PParticle = 1000; // Heat capacity
		double k_dc = 3.5; // k for normalized drying curve
		moistureContent X_cr = 0.025; // Critical moisture content [kg/kg]
		area A_P = 4; // total surface area of particle mass [m^2]
		length Delta_f = 40e-6; // Thickness of the water film on particles [m]
		thermalConductivity lambdaParticle = 0.58;
		moistureContent InitX = 0;
	// Bed
		double eps_0 = 0.4;
		double u_mf = 0;
		std::vector<chamberSection> chamber;
		double heightOfChamberTemperatureProbe = 0.070;
		double heighestFlowTimepoint = 0;
	// Settings
		bool calcBeta = GetCheckboxParameterValue("calcBeta");
		bool calcY_sat = GetCheckboxParameterValue("calcY_sat");
		bool calcNdc = GetCheckboxParameterValue("calcNdc");
		size_t dryingCurveSetting = GetComboParameterValue("DryingCurve");
		double SmallBiotNumber = 0.1;
		double f_alpha = 1; // ratio alpha_PF / alpha_AP
		double phiCuttOff = 0.999;
	// REA function parameters
		double REA1 = 0.96;
		double REA2 = -19.63;
		double REA3 = 0.73;

	size_t suspLayer = 0;

	bool particlesGlobal = true;

	size_t N_particle = 1; // Number of hight discretization layers of all sections containing particles
	size_t N_total = 1;// Total number of hight discretization layers

	CHoldup* m_holdup{}; // Holdup
	CMaterialStream* m_inLiquidStream{}; // Input of water stream
	CMaterialStream* m_inNozzleAirStream{}; // Input nozzle air
	CMaterialStream* m_inGasStream{}; // Input gas (fluidization air) stream
	CMaterialStream* m_outExhaustGasStream{};	// Output of exhaust gas

	CStream* m_VaporStream{};
	//CHoldup* m_Expander{}; //ToDo - check usage
	//CHoldup* workingHoldup{}; //ToDo - check usage
	

	// String keys of compounds for material database
	std::vector<std::string> compoundKeys;
	// Indices of liquid (first) and vapor (second) form of phase changing compound in compoundsKey
	std::pair<int,int> indicesOfVaporOfPhaseChangingCompound=std::make_pair(-1,-1);
	bool debugHasBeenShown = false;

	temperature TempLiquidOld = T_ref;
	temperature TempGasOld = T_ref;
	temperature TempSolidOld = T_ref;
	moistureContent YavgOld = Y_in;
	int DiffCoeff; // index of correlation for calcualting diffusion coefficient in the dropdown list

	double EnergyLiquidPhaseOld = 0;
	double EnergySolidPhaseOld = 0;
	double EnergyGasPhaseOld = 0;
	double HeatLossOld = 0;

public:
	void CreateBasicInfo() override;
	void CreateStructure() override;
	void Initialize(double _time) override;
	void SaveState() override;
	void LoadState() override;
	void Simulate(double _timeBeg, double _timeEnd) override;

	/// Get compound or phase properties///
	/**	Calculates sauter diameter of the particles
	 *	\param time stamp
	 *	\return sauter diameter in meter*/
	double CalculateHoldupSauter(double _time) const;
	
	// Returns gas moisture content for relativ humidity times saturation pressure [kg liquid per kg dry gas]
	moistureContent GetGasSaturationMoistureContent(temperature temperatureGas, pressure pressureGas = STANDARD_CONDITION_P);

	// Returns praticle equilibirum moisture content from material database, return 0 if no entry in data base is found [kg liqudi per kg dry solid]
	moistureContent CalcuateSolidEquilibriumMoistureContent(double _time, temperature temperature, double RH);

	// Returns normalized drying curve
	double CalculateNormalizedDryingCurve(moistureContent X, moistureContent Xeq);

	//	Calculates the diffusion coefficient following the correlation published in https://doi.org/10.1016/j.ijheatmasstransfer.2020.119500
	double CalculateDiffusionCoefficient(double _time, double avgGasTemperature, double filmTemperature, double pressure = STANDARD_CONDITION_P) const;

	// Calculates air relative humidity
	double GetRelativeHumidity(moistureContent Y, temperature temperature, pressure pressure = STANDARD_CONDITION_P); //const;

	// Calculates the ratio (Delta E_v / Delta E_v,eq) in case of REA
	double REA(double deltaX) const 
	{ 
		if (deltaX > 0) 
			return REA1 * exp(REA2 * pow(deltaX, REA3)); 
		else 
			return 1;
	};

	// Calculates the inverted REA, namely the evaporation rate in case of REA
	double REAinv(double Xeq, double y) const 
	{ 
		if (REA1 == 0 || REA2 == 0) 
			return 0; 
		else
			return pow(log(y / REA1) / REA2, 1. / REA3) + Xeq;
	};

	double GetAvgConstCompoundProperty(double _time, EPhase phase, ECompoundConstProperties  property) const;
	

	double GetAvgTPCompoundProperty(double _time, EPhase phase, ECompoundTPProperties property, double temperature, double pressure = STANDARD_CONDITION_P) const;
	
	double CalcAlphaOutside(double _time, const double h, const double D, const double Ts, EShape shape=EShape::CYLINDRICAL) const;

	double CalckAc(double alphaIn, double alphaOut, double L, std::vector<double> d/* Inner to outer diameter*/, std::vector<double> lambda) const;

	double CalckAp(double alphaIn, double alphaOut, std::vector<double> A, std::vector<double> delta, std::vector<double> lambda) const;
	
	area CalculateParticleSurfaceArea(double _time) const;

	moistureContent CalculateGasEquilibriumMoistureContent(temperature temperatureParticle, pressure pressureGas, double RH=1) const;

	double CalculateEquilibriumRelativeHumidity(double _time, temperature temperature, double X) const;

	double GetEquilibriumRelativeHumidity(double temperature, double X) const;

	double CalcAlphaParticleGas(double _time) const;

	bool CheckForSmallBiot(double _time) const;

	double CalculateMinFluidizeVel(double _time) const;

	double CalculateGasVel(double _time, size_t section = 0) const;
	
	double CalculateBedPorosity(double _time, bool homogeniusFluidization = true) const;

	double CalculateAlpha_PW(double _t_p, double _t_g, double _p, double _time) const;

	double CalculateAlpha_GW(double _time) const;

	double CalculateAlpha_GW(double _time, size_t section) const;

	double CalculateBetaPA(double _time, double avgGasTemperature, double filmTemperature, double D) const;

	double GetParticleEquilibriumMoistureContent(double temperature, double RH) const;
	
	//	Initializes variables containing equilibirum moisture content date
	bool InitializeMoistureContentDatabase(std::string path);

	double CalculateBedHeight(double _time, double particleTemperature);
	
	double CalculateBedHeightOrDetermineSectionsFilledWithBed(double _time, double particleTemperature, bool outputHeight = true);

	//double GetGasMassOfLayer(double _time, size_t layer, double gasTemperature, double particleTemperature);

	/// Dimensionless numbers ///
	// Reynolds number 
	double CalculateReynolds(double _time) const
	{
		return CalculateGasVel(_time) * CalculateHoldupSauter(_time) * rhoGas / etaGas;
	};
	
	// Reynolds number of each section of process chamber
	double CalculateReynolds(double _time, size_t section) const;
	
	//	Calculates the Prandtl number of the gas
	double CalculatePrandtl(double _time) const
	{
		return C_PGas * etaGas / lambdaGas;
	};

	//	Calculates Archimedes number for the gas particle combination using the sauter diameter
	double CalculateArchimedes(double _time) const
	{
		return STANDARD_ACCELERATION_OF_GRAVITY * pow(CalculateHoldupSauter(_time), 3) * (rhoParticle - rhoGas) * rhoGas / pow(etaGas, 2);
	};
	
	//	Converts Nusselt/Sherwood number of single particle to Nusselt/Sherwood number of fluidized bed
	double CalculateNusseltSherwoodBed(double porosity, double Nu_Sh_SingleParticle) const
	{
		return (1. + 1.5 * (1. - porosity)) * Nu_Sh_SingleParticle;
	};

	//	Calculates the Nusselt/Sherwood number of a single particle
	double CalculateNusseltSherwoodSingleParticle(double Nu_Sh_lam, double Nu_Sh_turb) const
	{
		return 2. + sqrt(pow(Nu_Sh_lam, 2) + pow(Nu_Sh_turb, 2));
	};

	//	Calculates the Nusselt/Sherwood number for the larminar part
	double CalculateNusseltSherwoodLam(double Re, double Pr_Sc) const
	{
		return 0.664 * pow(Pr_Sc, 1. / 3) * sqrt(Re);
	};

	//	Calculates the Nusselt/Sherwood number for the turbulent part
	double CalculateNusseltSherwoodTurb(double Re, double Pr_Sc) const
	{
		return (0.037 * pow(Re, 0.8) * Pr_Sc) / (1. + 2.443 * pow(Re, -0.1) * (pow(Pr_Sc, 2. / 3) - 1));
	};

	//	Calculate Biot number
	double CalcBiotNumber(double _time) const
	{
		return (CalcAlphaParticleGas(_time) * m_holdup->GetPhaseMass(_time, EPhase::SOLID) / (rhoParticle * CalculateParticleSurfaceArea(_time) * lambdaParticle));
	};

	//	Setup variables containing the dimensions and information of the chamber
	void SetupChamber();

	//Testing function, called during initialzation
	//void Testing();

	/**	Calculates overall heat transfer coefficient for cylindrical section.
	 *	\param section id of section in chamber object
	 *	\param alphaInternal heat transfer coefficient internal
	 *	\param alphaExternal heat transfer coefficient external
	 *	\param heightUsage use >0 for section partailly filled with particles, <0 for sections partially void of particles and =0 for sections void of particles range -1 to 1
	 *	\return overall heat transfer coefficient*/
	double CalculateOverallHeatTransferCoefficientCylinder(size_t section, double alphaInternal, double alphaExternal, double heightUsage = 0);

	//	Calculates overall heat transfer coefficient fo the top plate
	double CalculateOverallHeatTransferCoefficientTopPlate(double alphaInternal, double alphaExternal);
	
	double DetermineSectionsFilledWithBed(double _time, double particleTemperature);
	
	size_t DetermineLayersInSectionFilledWithBed(size_t section, double heigthUsage) const
	{
		return std::ceil(chamber.at(section).layers * heigthUsage);
	};

	//	Calculates volume of a section
	double CalculateSectionVolume(size_t section);
	
	/**	Determine into which section a layer falls
	 *	\param index of layer
	 *	\return index of section*/
	//size_t DetermineSectionForLayer(size_t layer) const;

	std::vector<double> GetGasMassOfLayers(double _time, double gasTemperature, double particleTemperature);

	double CalculateLayerVolume(size_t section, double R1, double r1)
	{
		return chamber.at(section).height* MATH_PI* (R1 * R1 + R1 * r1 + r1 * r1) / (3 * chamber.at(section).layers);
	};

	// linear interpolation: calculate the linear interpolation between lowNum and highNum
	double lerp(double lowNum, double highNum, double midNum) const
	{
		return lowNum + midNum * (highNum - lowNum);
	}
};