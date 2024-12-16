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
	size_t m_iYOutGas{}; // outlet gas absolute moisture, no height discretization
	size_t m_iTempOutGas{}; // outlet gas temperature, no height discretization == holdup gas temperature
	// particle (solid) phase
	size_t m_iTempParticle{}; // particle teemperature
	size_t m_iPhi{}; // particle wetness degree
	// liquid phase (water film)
	size_t m_iTempFilm{}; // water film (on particle surface) temperature

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
	
	//equilibriumMoistureContentData eqData;
	//std::set<double> RHs;
	//std::set<double> temperatures;
	//std::map< std::pair<double, double>, double> equilibriumMoistureContents;

public:
	const temperature T_ref = 273.15; // Ref. temperature for enthalpy [K] - default 273.15 K
	temperature T_env = 298.15; // Environment temperature [K] - default: Standard condition
	// Air
		std::string keyGas = "4e3a9D257E6A9E4F03D1"; // compound key for air
		density rhoGas = 1.2; // Density gas [kg/m^3] - default: air
		density rhoVapor = 0.8; // Density water vapor [kg/m3]
		dynamicViscosity etaGas = 1.8e-5; // Dynamic viscosity gas [Pa*s] - default: air
		const heatCapacity C_PGas = 1004; // Heat capacity gas [J/(kg*K)] - default: air
		thermalConductivity lambdaGas = 0.025; // Thermal conductivity gas [W/(m*K)] - default: air
		molarMass molarMassGas = 0.028949; // Molar mass of gas mixture [kg/mol] - default: air
	// Water liquid and vapor
		std::string keyLiquid = "4b3f8A1A71A315EFB4E5"; // compound key for liquid water
		std::string keyVapor = "Es8yKDAVn3QLJwWy0vNK"; // compound key for water vapor
		density rhoWater = 1000; // Density liquid [kg/m^3] - default: water
		const heatCapacity C_PWaterLiquid = 4190; // Heat capacity liquid phase change compound[J / (kg * K)] - default: water
		const heatCapacity C_PWaterVapor = 1890; // Heat capacity vapor phase change compound [J/(kg*K)] - default: water
		const specificLatentHeat Delta_h0 = 2500e3; // Specific latent heat (evaporation heat) phase change compound at 0 degree [J/kg] - default: water
		thermalConductivity lambdaWater = 0.6; // Thermal conductivity [W/(m*K)] - default: water
		molarMass molarMassPhaseChangingLiquid = 0.018; // Molar mass of phase changing liquid [kg/mol] - default: water
	// Particle (solid)
		std::string keySolid = "F3SLLny30KKziHtRbOE2"; // compound key for Cellets
		double wadellFactor = 0.95;
		density rhoParticle = 1500; // Particle density (Cellets: skeletal density)
		heatCapacity C_PParticle = 1250; // Heat capacity of MCC in [J/(kg*K)]
		thermalConductivity lambdaParticle = 0.2; // https://doi.org/10.1016/j.ijpharm.2017.10.018 MCC relative density ~= 0.7
	// Bed
		//length heightOfBed;
		//length diamOfBed;
		//double eps_0; // = 0.4;
		//double u_mf; // = 0;
	// process chamber
		thermalConductivity lambdaWall = 16; // thermal conductivity stainless steel (estimated), in [W/(m*K)]
		length wallThickness = 5e-3;
		length radiusChamber = 0.2;
		length lengthChamber = 0.35;
		length radiusExpanderLow = 0.4;
		length lengthExpanderLow = 0.4;
		length radiusExpanderHigh = 0.6;
		length lengthExpanderHigh = 0.95;
		length radiusExhaustAirPipe = 0.14;
		length lengthExhaustAirPipe = 0.1;
		//std::vector<chamberSection> chamber;
		//double heightOfChamber;
		//double heightOfChamberTemperatureProbe;// = 0.070;
		//double heightOfNozzle;
		//double heighestFlowTimepoint = 0;

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

	// String keys of compounds for material database
	//std::vector<std::string> compoundKeys;
	// Indices of liquid (first) and vapor (second) form of phase changing compound in compoundsKey
	//std::pair<int,int> indicesOfVaporOfPhaseChangingCompound=std::make_pair(-1,-1);
	//bool debugHasBeenShown = false;

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
	double CalculateDiffusionCoefficient(double avgGasTheta) const // Dosta(2010) [m2/s]
	{
		return 2.3e-5 * pow((avgGasTheta + 273.15) / 273.15, 1.81);
	};
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
	dimensionlessNumber CalculateReynoldsMF(double _time, length d32) const;
	dimensionlessNumber CalculateReynolds(double _time, length d32) const;
	// Reynolds number with height discretization for each layer
	//dimensionlessNumber CalculateReynolds(double _time, size_t section) const;
	dimensionlessNumber CalculatePrandtl(temperature avgGasTemperature) const;
	dimensionlessNumber CalculateSchmidt(double D_a) const;
	dimensionlessNumber CalculateArchimedes(length d32) const;
	dimensionlessNumber CalculateNusseltSherwood(dimensionlessNumber Nu_Sh_lam, dimensionlessNumber Nu_Sh_turb) const
	{
		return 2. + sqrt(pow(Nu_Sh_lam, 2.) + pow(Nu_Sh_turb, 2.));
	};
	dimensionlessNumber CalculateNusseltSherwoodLam(dimensionlessNumber Re, dimensionlessNumber Pr_Sc) const
	{
		return 0.664 * pow(Pr_Sc, 1. / 3.) * sqrt(Re);
	};
	dimensionlessNumber CalculateNusseltSherwoodTurb(dimensionlessNumber Re, dimensionlessNumber Pr_Sc) const
	{
		dimensionlessNumber Nu_Sc = (0.037 * pow(Re, 0.8) * Pr_Sc) / (1. + 2.443 * pow(Re, -0.1) * (pow(Pr_Sc, 2. / 3.) - 1));
		if (Nu_Sc < 0)
		{
			Nu_Sc = 0;
		}
		return Nu_Sc;
	};
	dimensionlessNumber CalculateNusseltSherwoodApp(dimensionlessNumber Nu_Sh, double eps) const
	{
		return Nu_Sh * (1. + 1.5 * (1. - eps));
	};
	dimensionlessNumber CalculateNusseltSherwoodModify(dimensionlessNumber Re, dimensionlessNumber Sc_Pr, dimensionlessNumber Sh_Nu_app,  dimensionlessNumber AtoF) const // Groenewold & Tsotsas -> see Rieck diss. and Soeren diss.
	{
		return (Re * Sc_Pr) * log(1. + (Sh_Nu_app * AtoF) / (Re * Sc_Pr)) / AtoF;
	};
	dimensionlessNumber CalculateNusseltFree(double _time, temperature T_surface) const;
	dimensionlessNumber CalculateGrashof(temperature T_env, temperature T_surface) const
	{
		return STANDARD_ACCELERATION_OF_GRAVITY * pow(lengthChamber, 3) * (T_surface - T_env) / (pow(etaGas / rhoGas, 2) * T_env);
	};

	//	Calculate Biot number
	//dimensionlessNumber CalcBiotNumber(double _time, temperature avgGasTemperature, length d32) const
	//{
	//	return (CalculateAlpha_GP(_time, avgGasTemperature, d32) * m_holdup->GetPhaseMass(_time, EPhase::SOLID) / (rhoParticle * CalculateParticleSurfaceArea(_time) * lambdaParticle));
	//};
	//bool CheckForSmallBiot(double _time) const;

/// Heat transfer coefficients ///
	double CalculateAlpha_GP(double _time, temperature avgGasTemperature, length d32) const; 
	double CalculateAlpha_GF(double _time, temperature avgGasTemperature, length d32) const;
	double CalculateAlpha_PF(/*temperature tempWater, pressure pressureHoldup, length d32*/ double alpha_GP) const;

/// Heat loss through the walls ///
	temperature GetNewTempSurface(double _time, temperature TempSurfaceOld, temperature thetaInside, length d32) ;
	temperature IterateSurfaceTemp(double _time, temperature T_gasHoldup, length d32) ;
	double CalculateHeatLossWall(double _time, length wallThickness, length height, length radiusInner, temperature thetaIn, temperature thetaSurface, double lambdaWall, length d32) const;

/// Mass transfer coefficient ///
	massTransferCoefficient CalculateBeta(double _time, length d32, double avgGasTheta, double D_a) const;

/// Drying kinetics ///
	//moistureContent CalculateSolidEquilibriumMoistureContent(double _time, temperature temperature, double RH); // X_eq
	double CalculateRelativeDryingRate(moistureContent X) const;
	pressure CalculateGasSaturationPressure(temperature theta_Gas, pressure pressureGas) const; // Y_eq
	moistureContent CalculateGasEquilibriumMoistureContent(pressure pressureGas, pressure P_sat, double RH) const; // Y_eq
	double CalculateGasEquilibriumRelativeHumidity(/*double _time, temperature temperature,*/ moistureContent X) const; // RH_eq

	//double GetEquilibriumRelativeHumidity(double temperature, double X) const;
	////	Initializes variables containing equilibirum moisture content date
	//bool InitializeMoistureContentDatabase(std::string path);
	//moistureContent GetParticleEquilibriumMoistureContent(double temperature, double RH) const;


///  Hydrodynamics ///
	double CalculateMinFluidizeVel(double _time, length d32) const;
	double CalculateGasVel(double _time, length d32) const;
	double CalculateBedPorosity(double _time, length d32, bool homogeneousFluidization = false) const;
	double CalculateBedPorosityMF(double wadellFactor) const
	{
		return pow(1. / (14. * wadellFactor), 1. / 3.);
	}
	length CalculateFluidizedBedHeight(/*double _time, double particleTemperature*/length H_fix, double eps) const;
	dimensionlessNumber CalculateAtoF(length H_fb, length d32, double eps) const
	{
		return 6 * (1. - eps) * H_fb / d32; // Soeren diss. eq. (C.13)
	};
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