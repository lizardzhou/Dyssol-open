/* Copyright (c) 2020, Dyssol Development Team. All rights reserved. This file is part of Dyssol. See LICENSE file for license information. */

#pragma once

#include "Types.h"

/// Base-model ///
class CNLModelBase : public CNLModel
{
private:
	// Internal unit for property calculation
	void* m_unit{};
	// Compound keys
	SCompoundKeys m_compoundKeys;
	// Unit constants
	SConstantParameters m_constants;
	// Port streams
	SPorts m_inPorts;
	// Discretization parameters
	SDiscretizationParameters m_discretization;
	// TD parameters
	STDParameters m_TDParameters;
	// Simulation settings
	SSettings m_settings;

public:
	/// Setup ///
	void Setup(void* _unit, const SCompoundKeys& _compoundKeys, const SConstantParameters& _constants, const SPorts& _inPort, const SDiscretizationParameters& _discretization, const SSettings& _settings);
	/// Get functions ///
	void* GetMaterialsModel() const { return m_unit; }
	SCompoundKeys GetCompoundKeys() { return m_compoundKeys; }

	size_t GetHeightGridPoints() const { return m_discretization.heightGridPoints; }
	size_t GetRTGridPoints() const { return m_discretization.RTGridPoints; }
	size_t GetSizeClasses() const { return m_discretization.sizeClasses; }
	size_t GetMoistureClasses() const { return m_discretization.moistureClasses; }
	size_t GetFDHeightGridPoints() const { return m_discretization.FDHeightGridPoints; }

	double GetInitialBedHeight() const { return m_constants.bed.h_0; }
	double GetInitialBedMass() const { return m_constants.bed.m_0; }
	double GetManualMinimumVelocity() const { return m_constants.bed.u_mf_man; }
	double GetInitialBedPorosity() const { return m_constants.bed.eps_0; }
	double GetVibrationIntensity() const { return m_constants.bed.Lambda; }
	double GetWallTemperature() const { return m_constants.bed.T_Wall; }

	std::string GetCSVPath() { return m_settings.CSVPath; }
	double GetSettingsQpp() const { return m_settings.C_Qpp; }
	double GetSettingsQw() const { return m_settings.C_Qw; }
	double GetSettingsNtaumax() const { return m_settings.C_Ntaumax; }
	EArea GetSettingsApAbed() const { return m_settings.C_ApAbed; }
	double GetSettingsKtis() const { return m_settings.C_Ktis; }
	bool GetSettingsSwitchHT() const { return m_settings.HT; }
	bool GetSettingsSwitchMT() const { return m_settings.MT; }

	double GetSettingsH_NTU() const { return m_settings.C_H_NTU; }
	double GetSettingsBubbleVolFractOverwrite() const { return m_settings.C_BubbleVolFraction; }

	double GetBedDiameter() const { return m_constants.geometry.diameter; }
	double GetBedArea() const { return m_constants.geometry.area; }
	size_t GetNumberOfOrifices() const { return m_constants.geometry.numOrifice; }
	double GetBedCharacteristicLength() const { return m_constants.geometry.shape == EShape::RECTANGULAR ? m_constants.geometry.length1 : m_constants.geometry.diameter; }
	double GetDryerCircumference() const { return m_constants.geometry.circumference; }

	double GetReferenceTemperature() const { return m_constants.unit.T_0; }
	double GetEvaporationEnthalpy() const { return m_constants.unit.hv_H2O; }
	double GetDryingCurveK() const { return m_constants.unit.k_dc; }
	double GetPhiEq() const { return m_constants.unit.phi_eq_p; }
	double GetXcr() const { return m_constants.unit.X_cr; }

	double GetUmf() const { return m_TDParameters.fluidDynamics.u_mf; }
	double GetRemf() const { return m_TDParameters.fluidDynamics.Re_mf; }
	double GetU0() const { return m_TDParameters.fluidDynamics.u_0; }
	double GetTheta() const { return m_TDParameters.fluidDynamics.theta; }
	SGeldart GetGeldart() const { return m_TDParameters.geldart; }


	/// Calculation functions ///
	STDParameters CalculateTDParameters(double _time, EGeldart _geldart);
	SBoundaryCondition CalculateBoundaryConditions(double _time);
	SBoundaryCondition CalculateInitialConditions(SBoundaryCondition _bc);
	double CalculateTheta(const SGeldart& _geldart) const;
	double CalculateParticleMeanMoisture(double _time, const std::vector<double>& _mfrac_f);

	/// Set functions ///
	void SetUnit(void* _unit) { m_unit = _unit; }
	void SetCompoundKeys(const SCompoundKeys& _compoundKeys) { m_compoundKeys = SCompoundKeys(_compoundKeys); }
	void SetConstantParameters(const SConstantParameters& _constantParameters) { m_constants = SConstantParameters(_constantParameters); }
	void SetInPorts(const SPorts& _ports) { m_inPorts = SPorts(_ports); }
	void SetDiscretization(const SDiscretizationParameters& _discretizationParams) { m_discretization = SDiscretizationParameters(_discretizationParams); }
	void SetTDParameters(const STDParameters& _tdParameters) { m_TDParameters = STDParameters(_tdParameters); }
	void SetSettings(const SSettings& _settings) { m_settings = SSettings(_settings); }

	/// Stream functions ///
	SMYTP GetInGasMYTP(double _time) { return SMYTP{ m_inPorts.gas, _time, EPhase::VAPOR, GetCompoundKeys().air, GetCompoundKeys().H2Og }; }
	SMYTP GetTDInGasMYTP() const { return m_TDParameters.inGas; }
	SMYTP GetInParticleMYTP(double _time) { return SMYTP{ m_inPorts.particle, _time, EPhase::SOLID, GetCompoundKeys().particles, GetCompoundKeys().H2Ol }; }
	std::vector<double> GetInParticlePSD(double _time) const { return m_inPorts.particle->GetPSD(_time, PSD_MassFrac); }
	std::vector<double> GetInParticleq3(double _time) const { return m_inPorts.particle->GetPSD(_time, PSD_q3); }
	SMYTP GetTDInParticleMYTP() const { return m_TDParameters.inParticle; }
	double GetInGasMoisture(double _time) const;
	double GetInGasEnthalpy(double _time);
	std::vector<std::vector<double>> GetInParticleMassfractions(double _time) const;
	double GetInParticleMeanMoisture(double _time);
	std::vector<double> GetInParticleEnthalpy(double _time);
	double GetInParticleDensity(double _time) const;
	double GetInParticleSauter(double _time) const;
	double GetInGasDensity(double _time) const;
	double GetInGasViscosity(double _time) const;
	/// Material property functions ///
	// Constant properties
	double GetGasMolarMass() const;
	// Gas properties
	double GetDiffusionCoefficient(double _t, double _p) const;
	double GetDryGasDensity(double _t, double _p) const;
	double GetMoistGasDensity(double _t, double _p, double _y) const;
	double GetDryGasViscosity(double _t, double _p) const;
	double GetDryGasHeatCapacity(double _t, double _p) const;
	double GetWaterHeatCapacity(double _t, double _p) const;
	double GetParticleHeatCapacity(double _t, double _p) const;
	double GetDryGasThermalConductivity(double _t, double _p) const;
	double GetGasEquilibriumMoistureContent(double _t_p, double _phi_eq, double _p) const;
	double GetMoistGasEnthalpy(double _t, double _p, double _y) const;
	double GetMoistGasTemperature(double _h, double _y, double _t_prop, double _p_prop) const;
	double GetRelativeHumidity(double _y, double _t, double _p) const;
	double GetParticleEquilibriumMoistureContent(double _phi, double _t) const;
	double GetWaterVaporEnthalpy(double _t, double _p) const;
	// Particle properties
	double GetDryParticleDensity(double _t, double _p) const;
	double GetMoistParticleDensity(double _x_p, double _t, double _p) const;
	double GetMoistParticleHeatCapacity(double _x_p, double _t, double _p) const;
	double GetMoistParticleEnthalpy(double _t, double _p, double _x) const;
	double GetMoistParticleTemperature(double _h, double _x, double _t_prop, double _p_prop) const;

	/// Unit functions ///
	std::vector<double> GetParticleSizeGrid() const;
	std::vector<double> GetMeanDiameterVector() const;
	std::vector<double> GetMeanMoistureVector() const;
};


/// Fluid dynamics model ///
class CNLModelFD : public CNLModelBase
{
private:
	// NLVariable indices
	std::vector<size_t> m_nd_v_z;	// Indices of bubble diameters
	size_t m_h_fb{};				// Index of total bed height

	// Results of FD NLModel
	SResultsNLModelFD m_resultsNLModelFD;

public:
	CNLModelFD() = default;
	~CNLModelFD() override;

	/// Declaration / Initialization ///
	void AddVariables();
	void Initialize(double _time, const STDParameters& _tdParameters);

	/// Get functions ///
	SResultsNLModelFD GetResultsNLModelFD() { return m_resultsNLModelFD; }
	double GetBedHeight() const { return m_resultsNLModelFD.h_fb; }
	std::vector<double> GetHeightVector() { return m_resultsNLModelFD.h_z; }
	std::vector<double> GetDimensionlessHeightVector() { return m_resultsNLModelFD.zeta_z; }
	std::vector<double> GetTDBubbleDiameterVector() { return m_resultsNLModelFD.d_v_z; }
	std::vector<double> GetTDGasSplitVector() { return m_resultsNLModelFD.nu_z; }

	/// Set functions ///
	void SetResultsNLModelFD(const SResultsNLModelFD& _resultsNLModelFD) { m_resultsNLModelFD = _resultsNLModelFD; }


	/// Solver functions ///
	void CalculateFunctions(double* _vars, double* _func, void* _unit) override;
	void ResultsHandler(double _time, double* _vars, void* _unit) override;

	/// Additional functions ///
	double Calculated_v_z_update(double _d_v_z_1, double _d_v_z, double _eps_b_z, double _u_b_z, double _delta_h_z, size_t _heighIndex, const SGeldart& _geldart, size_t _num_or, double _lambda);
	double CalculatePsi(double _z, double _l_char_fb, const SGeldart& _geldart) const;

	double CalculatedA_0_GeldartC() const;

	/// Clear functions ///
	void ClearModel();
};

class CNLModelTD : public CNLModelBase
{
private:
	/// NLVariable indices ///
	// Bubble
	std::vector<size_t> m_iY_b_z; // Yb(zeta)
	std::vector<size_t> m_ih_b_z; // hb(zeta)
	// Suspension gas
	std::vector<size_t> m_iY_s_z; // Ys(zeta)
	std::vector<size_t> m_ih_s_z; // hs(zeta)
	// Particles
	std::vector<std::vector<std::vector<size_t>>> m_iX_p_dfi;
	std::vector<std::vector<std::vector<size_t>>> m_ih_p_dfi;

	/// Discretization ///
	size_t m_RTSubGridPoints{};
	size_t m_HeightSubGridPoints{};

	/// State vectors ///
	//	Particle phase (p)
	std::vector<std::vector<std::vector<double>>> m_X_p_dfi;		// 3D-vector for storage of moisture of particles of size class d, moisture class f and residence class i
	std::vector<std::vector<std::vector<double>>> m_h_p_dfi;		// 3D-vector for storage of enthalpy of particles of size class d, moisture class f and residence class i
	std::vector<std::vector<std::vector<double>>> m_X_p_dfi_update; // 3D-vector for storage of updated moisture of particles of size class d, moisture class f and residence class i
	std::vector<std::vector<std::vector<double>>> m_h_p_dfi_update; // 3D-vector for storage of updated enthalpy of particles of size class d, moisture class f and residence class i
	//	Suspension phase (s)
	std::vector<double> m_Y_s_z;		// 1D-vector for storage of water loading of suspension phase at height class z
	std::vector<double> m_h_s_z;		// 1D-vector for storage of enthalpy of suspension phase at height class z
	std::vector<double> m_Y_s_z_update;	// 1D-vector for storage of updated water loading of suspension phase at height class z
	std::vector<double> m_h_s_z_update;	// 1D-vector for storage of updated enthalpy of suspension phase at height class z
	//	Bubble phase (b)
	std::vector<double> m_Y_b_z;		// 1D-vector for storage of water loading of bubble phase at height class z
	std::vector<double> m_h_b_z;		// 1D-vector for storage of enthalpy of bubble phase at height class z
	std::vector<double> m_Y_b_z_update;	// 1D-vector for storage of updated water loading of bubble phase at height class z
	std::vector<double> m_h_b_z_update;	// 1D-vector for storage of updated enthalpy of bubble phase at height class z

	/// Unit parameters ///
	SBoundaryCondition m_BC;
	STDModelParameters m_parameters;
	SResultsNLModelTD m_results;

public:
	CNLModelTD() = default;
	~CNLModelTD() override;

	/// Set functions ///
	void SetBoundaryConditions(const SBoundaryCondition& _bc) { m_BC = _bc; }
	void SetTDModelParameters(const STDModelParameters& _parameters) { m_parameters = _parameters; }
	void SetTDModelResults(const SResultsNLModelTD& _resultsNLModelTD) { m_results = _resultsNLModelTD; }

	/// Get functions ///
	SBoundaryCondition GetBoundaryConditions() const { return m_BC; }
	SResultsNLModelTD GetTDModelResults() const { return m_results; }

	size_t GetRTSubGridpoints() const { return m_RTSubGridPoints; }
	size_t GetHeightSubGridPoints() const { return m_HeightSubGridPoints; }
	std::vector<double> GetRatioApAbedVector() const { return m_parameters.ApAbed_d; }
	double GetInSauter() const { return m_parameters.dp_sauter; }
	double GetBedPorosity() const { return m_parameters.eps_fb; }
	double GetTauMax() const { return m_parameters.tau_max; }
	double GetBedHeight() const { return m_parameters.h_fb; }
	std::vector<std::vector<double>> GetParticleSurfaceVector() const { return m_parameters.Aps_df; }
	std::vector<double> GetTDMeanDiameterVector() const { return m_parameters.d_d; }
	std::vector<std::vector<double>> GetMassFractionVector() const { return m_parameters.mfrac_df; }
	std::vector<double> GetHeightVector() const { return m_parameters.h_z; }
	double GetTDGasSplitFactor() const { return m_parameters.nu_0; }
	std::vector<double> GetRTTauVector() const { return m_parameters.tau_i; }
	std::vector<double> GetRTNumberDensityVector() const { return m_parameters.n_i; }
	std::vector<double> GetDimensionlessHeightVector() const { return m_parameters.zeta_z; }

	/// Declaration / Initialization ///
	void AddVariables(double _time);
	void Initialize(double _time, const STDParameters& _tdParameters, const SResultsNLModelFD& _resultsNLModelFD);

	/// Solver functions ///
	void CalculateFunctions(double* _vars, double* _func, void* _unit) override;
	void ResultsHandler(double _time, double* _vars, void* _unit) override;

	/// Exchange flow functions ///
	void ExchangeParticleSuspension(const std::vector<double>& _zeta_z, const std::vector<double>& _tau_i, const std::vector<double>& _n_i, const std::vector<double>& _y_s_z, const std::vector<double>& _t_s_z, const std::vector<std::vector<std::vector<double>>>& _x_p_dfi, const std::vector<std::vector<std::vector<double>>>& _t_p_dfi, std::vector<double>& _m_ps_z, std::vector<double>& _h_ps_z, std::vector<double>& _q_ps_z, std::vector<std::vector<std::vector<double>>>& _m_ps_dfi, std::vector<std::vector<std::vector<double>>>& _h_ps_dfi, std::vector<std::vector<std::vector<double>>>& _q_ps_dfi);
	void ExchangeSuspensionBubble(const std::vector<double>& _y_s_z, const std::vector<double>& _t_s_z, const std::vector<double>& _y_b_z, const std::vector<double>& _t_b_z, std::vector<double>& _m_sb_z, std::vector<double>& _h_sb_z, std::vector<double>& _q_sb_z);
	void ExchangeParticleParticle(const std::vector<double>& _tau_i, const std::vector<double>& _n_i, const std::vector<std::vector<std::vector<double>>>& _x_p_dfi, const std::vector<std::vector<std::vector<double>>>& _t_p_dfi, std::vector<std::vector<std::vector<double>>>& _q_pp_dfi);
	void ExchangeParticleWall(const std::vector<double>& _tau_i, const std::vector<double>& _n_i, const std::vector<std::vector<std::vector<double>>>& _x_p_dfi, const std::vector<std::vector<std::vector<double>>>& _t_p_dfi, std::vector<std::vector<std::vector<double>>>& _q_pp_dfi);
	/// Balance functions ///
	void BalanceParticlePhase(const std::vector<double>& _vtau_i, const std::vector<std::vector<std::vector<double>>>& _x_p_dfi, const std::vector<std::vector<std::vector<double>>>& _h_p_dfi, const std::vector<std::vector<std::vector<double>>>& _t_p_dfi, std::vector<std::vector<std::vector<double>>>& _x_p_dfi_update, std::vector<std::vector<std::vector<double>>>& _h_p_dfi_update, const std::vector<std::vector<std::vector<double>>>& _m_ps_dfi, const std::vector<std::vector<std::vector<double>>>& _q_ps_dfi, const std::vector<std::vector<std::vector<double>>>& _h_ps_dfi, const std::vector<std::vector<std::vector<double>>>& _q_pp_dfi, const std::vector<std::vector<std::vector<double>>>& _q_pw_dfi);
	void BalanceSuspensionPhase(const std::vector<double>& _zeta_z, const std::vector<double>& _y_s_z, const std::vector<double>& _h_s_z, const std::vector<double>& _t_s_z, std::vector<double>& _y_s_z_update, std::vector<double>& _h_s_z_update, const std::vector<double>& _m_ps_z, const std::vector<double>& _h_ps_z, const std::vector<double>& _q_ps_z, const std::vector<double>& _m_sb_z, const std::vector<double>& _h_sb_z, const std::vector<double>& _q_sb_z) const;
	void BalanceBubblePhase(const std::vector<double>& _zeta_z, const std::vector<double>& _y_b_z, const std::vector<double>& _h_b_z, const std::vector<double>& _t_b_z, std::vector<double>& _y_b_z_update, std::vector<double>& _h_b_z_update, const std::vector<double>& _m_sb_z, const std::vector<double>& _h_sb_z, const std::vector<double>& _q_sb_z) const;

	/// Additional functions ///
	static double IntegrateTrapez(double _val1, double _val2, double _step);
	void CalculateParticleTemperatures(const std::vector<std::vector<std::vector<double>>>& _h_p_dfi, const std::vector<std::vector<std::vector<double>>>& _x_p_dfi, std::vector<std::vector<std::vector<double>>>& _t_p_dfi) const;
	void CalculateGasTemperatures(const std::vector<double>& _h_g_z, const std::vector<double>& _y_g_z, std::vector<double>& _t_g_z) const;
	double CalculateTauMax(double _n_tau_max) const;
	std::vector<double> CalculateRTTauVector(double _tau_max) const;
	std::vector<double> CalculateRTDensityVector(double _tau_max) const;
	double CaclulateBedPorosity(double _time, double _h_fb) const;
	double CalculateMoistureMeanValue(const std::vector<std::vector<std::vector<double>>>& _vals, size_t _rtSubGridPoints) const;
	double CalculateTemperatureMeanValue(const std::vector<std::vector<std::vector<double>>>& _vals, size_t _rtSubGridPoints) const;
	std::vector<double> CalculateRatioApAbed(double _time, double _h_fb, const std::vector<std::vector<double>>& _mfrac_df);
	std::vector<std::vector<double>> CalculateParticleSurfaceVector(double _time, const std::vector<double>& _d_d);
	std::vector<std::vector<double>> CalculateMassFractionVector(double _time) const;
	double CalculateAlpha_PS(double _t_ps, double _p_g, double _d_p, double _ratioApAbed_d) const;
	double CalculateBeta_PS(double _t_ps, double _y_s, double _p, double _d_p, double _ratioApAbed_d) const;
	double CalculateNormalizedDryingCurve(double _x, double _x_eq) const;
	double CalculateAlphaA_SB(double _t_sb, double _p, double _beta_sb_A_sb) const;
	double CalculateBetaA_SB(double _t_sb, double _p, double _h_fb, double _m_g_dot) const;
	double CalculateAlpha_PW(double _t_p, double _x_p, double _t_g, double _p) const;
	double CalculateAlpha_SW(double _t_s, double _y_s, double _p, double _h_fb, double _v_s) const;
	double CalculateA_PW() const;
	double CalculateAlpha_PP(double _t_p, double _x_p, double _t_g, double _p) const;
	double CalculateA_PP(double _t_p_mean, double _p) const;

	/// Clear functions ///
	void ClearModel();
};

/// FB-dryer unit ///
class CVibratedFluidizedBedDryer : public CSteadyStateUnit
{
private:
	/// Solvers and models ///
	CNLModelFD m_FDNLModel{};
	CNLSolver m_FDNLSolver{};
	CNLModelTD m_TDNLModel{};
	CNLSolver m_TDNLSolver{};

public:
	/// Port streams ///
	CMaterialStream* m_inGas{};			// Pointer to inlet gas stream
	CMaterialStream* m_inParticle{};	// Pointer to inlet particle stream
	CMaterialStream* m_outGas{};		// Pointer to outlet gas stream
	CMaterialStream* m_outParticle{};	// Pointer to outlet particle stream

	/// Compound keys ///
	SCompoundKeys m_compKeys{};

public:
	void CreateBasicInfo() override;
	void CreateStructure() override;
	void Initialize(double _time) override;
	void Simulate(double _time) override;
	void SaveState() override;
	void LoadState() override;

	/// Property functions ///
	// Constant material properties
	double GetMolarMass(const std::string& _compKey) const;
	// TPD Compound properties
	double GetCompoundHeatCapacity(double _t_eval, double _p_eval, const std::string& _compKey) const;
	double GetCompoundGasDensity(double _t_eval, double _p_eval, const std::string& _compKey) const;
	double GetCompoundViscosity(double _t_eval, double _p_eval, const std::string& _compKey) const;
	double GetCompoundThermalConductivity(double _t_eval, double _p_eval, const std::string& _compKey) const;
	// Gas properties
	double GetMoistGasDensity(double _t, double _p, double _y, const std::string& _keyAir) const;
	double GetMoistGasEnthalpy(double _t, double _p, double _y, double _t_0, double _hv_H2O, const std::string& _keyAir, const std::string& _keyH2Ol) const;
	double GetMoistGasTemperature(double _h, double _y, double _t_eval, double _p_eval, double _t_0, double _hv_H2O, const std::string& _keyAir, const std::string& _keyH2Ol) const;
	double GetRelativeHumidity(double _y, double _t, double _p, const std::string& _keyAir, const std::string& _keyH2Ol) const;
	double GetGasEquilibriumMoistureContent(double _t_p, double _phi_eq, double _p, const std::string& _keyAir, const std::string& _keyH2Ol, const std::string& _keyH2Og) const;
	double GetDiffusionCoefficientOld(double _t, double _p, const std::string& _keyAir, const std::string& _keyH2Og) const;
	double GetDiffusionCoefficient(double _t, double _p) const;
	// Water properties
	double GetWaterVaporEnthalpy(double _t, double _p, double _t_0, double _hv_H2O, const std::string& _keyH2Ol) const;

	// Particle properties
	double GetDryParticleDensity(double _t_eval, double _p_eval, const std::string& _keyParticle) const;
	double GetMoistParticleDensity(double _x_p, double _t_eval, double _p_eval, const std::string& _keyParticle) const;
	double GetMoistParticleHeatCapacity(double _x_p, double _t, double _p, const std::string& _keyParticle, const std::string& _keyH2Ol) const;
	double GetMoistParticleEnthalpy(double _t, double _p, double _x, double _t_0, const std::string& _keyParticle, const std::string& _keyH2Ol) const;
	double GetMoistParticleTemperature(double _h, double _x, double _t_eval, double _p_eval, double _t_0, const std::string& _keyParticle, const std::string& _keyH2Ol) const;
	double GetParticleEquilibriumMoistureContent(double _phi, double _t, const std::string& _keyParticle) const { return GetCompoundProperty(_keyParticle, ECompoundTPProperties::EQUILIBRIUM_MOISTURE_CONTENT, _t, _phi); }

	std::vector<double> GetMeanDiameterVector() const;
	std::vector<double> GetMeanMoistureVector() const;

	/// Auxiliary functions ///
	SGeometry ReadGeometryParameter();

	/// Unused functions ///
	void WriteCSV();
	void WriteParameters(const std::string& _resultsPath);
	void WriteStateVariables(const std::string& _resultsPath);
	void WritePlots(const std::string& _resultsPath) const;
};
