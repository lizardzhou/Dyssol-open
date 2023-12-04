/* Copyright (c) 2022, Dyssol Development Team. All rights reserved. This file is part of Dyssol. See LICENSE file for license information. */

#pragma once

#include "UnitDevelopmentDefines.h"

enum class EShape : size_t
{
	CYLINDRICAL = 0, RECTANGULAR = 1
};

enum class EGeldart : size_t
{
	AUTO = 0, A, B, C, D
};

enum class EArea : size_t
{
	SAUTER_DIAMETER, MEAN_AND_FRACTIONS, MEAN
};

struct SCompoundKeys
{
	std::string particles;	// Key of the particle compound in the material database.
	std::string H2Ol;		// Key of the H2O (l) compound in the material database.
	std::string H2Og;		// Key of the H2O (g) compound in the material database.
	std::string air;		// Key of the air compound in the material database.

	SCompoundKeys() = default;
	SCompoundKeys(std::string _keyParticle, std::string _keyH2Ol, std::string _keyH2Og, std::string _keyAir)
		: particles{ std::move(_keyParticle) }
		, H2Ol{ std::move(_keyH2Ol) }
		, H2Og{ std::move(_keyH2Og) }
		, air{ std::move(_keyAir) }
	{}
};

struct SPorts
{
	CMaterialStream* gas{};
	CMaterialStream* particle{};

	SPorts() = default;
	SPorts(CMaterialStream* _gas, CMaterialStream* _particle)
		: gas{ _gas }
		, particle{ _particle }
	{}
};

struct SGeometry
{
	EShape shape{ EShape::CYLINDRICAL };	// Shape of fluidized bed reactor
	double diameter{};
	double area{};
	double length1{};
	double length2{};
	size_t numOrifice{};
	double circumference{};

	SGeometry() = default;
	// Conical
	SGeometry(double _diameter, size_t _numOrifice)
		: shape{ EShape::CYLINDRICAL }
		, diameter{ _diameter }
		, area{ 0.25 * MATH_PI * _diameter * _diameter }
		, numOrifice{ _numOrifice }
		, circumference{ MATH_PI * _diameter }
	{}
	// Rectangular
	SGeometry(double _length1, double _length2, size_t _numOrifice)
		: shape{ EShape::RECTANGULAR }
		, diameter{ 4 * (_length1 * _length2) / (2 * _length1 + 2 * _length2) }
		, area{ _length1 * _length2 }
		, length1{ fmin(_length1, _length2) }
		, length2{ fmax(_length1, _length2) }
		, numOrifice{ _numOrifice }
		, circumference{ 2 * _length1 + 2 * _length2 }
	{}
};

struct SUnit
{
	double T_0{};		// Reference temperature [K]
	double hv_H2O{};	// Heat of evaporation of H2O [J/kg]
	double phi_eq_p{};
	double X_cr{};
	double k_dc{};

	SUnit() = default;
	SUnit(double _t_0, double _hv_H2O, double _phi_eq_p, double _x_cr, double _k_dc)
		: T_0{ _t_0 }
		, hv_H2O{ _hv_H2O }
		, phi_eq_p{ _phi_eq_p }
		, X_cr{ _x_cr }
		, k_dc{ _k_dc }
	{}
};

struct SBed
{
	double eps_0{};				// Initial Bed porosity [-]
	double h_0{};				// Initial bed height [m]
	double u_mf_man{};			// Manual value for minimum fluidization velocity [m/s]
	double m_0{};				// Initial bed mass (dry) [kg]
	double Lambda{};			// Vibration intensity [-]
	double T_Wall{ 293.15 };	// Temperature at inner wall [K]

	SBed() = default;
	SBed(double _eps_0, double _h_0, double _u_mf_man, double _rho_s, double _f, double _af, double _t_Wall, const SGeometry& _geometryParams)
		: eps_0{ _eps_0 }
		, h_0{ _h_0 }
		, u_mf_man{ _u_mf_man }
		, m_0{ _geometryParams.area * _h_0 * _rho_s * (1 - _eps_0) }
		, Lambda{ pow(2. * MATH_PI * _f, 2.) * _af / STANDARD_ACCELERATION_OF_GRAVITY }
		, T_Wall{ _t_Wall }
	{}
};

struct SConstantParameters
{
	SUnit unit{};			// General unit parameters
	SGeometry geometry{};	// Geometric parameters of apparatus
	SBed bed{};				// Bed parameters

	SConstantParameters() = default;
	SConstantParameters(const SUnit& _unit, const SGeometry& _geometry, const SBed& _bed)
		: unit{ _unit }
		, geometry{ _geometry }
		, bed{ _bed }
	{}
};

struct SMYTP
{
	double m{};
	double Y{};
	double T{};
	double P{};
	CMaterialStream* m_stream{};

	SMYTP() = default;
	SMYTP(double _m_dry, double _y, double _t, double _p)
		: m{ _m_dry }
		, Y{ _y }
		, T{ _t }
		, P{ _p }
	{}
	SMYTP(CMaterialStream* _stream, double _time, EPhase _phase, const std::string& _keyCompound, const std::string& _keyH2O)
		: m{ _stream->GetCompoundMassFlow(_time, _keyCompound, _phase) }
		, Y{ _stream->GetCompoundMassFlow(_time, _keyH2O, _phase) / m }
		, T{ _stream->GetTemperature(_time) }
		, P{ _stream->GetPressure(_time) }
		, m_stream{ _stream }
	{}
};

struct SFluidDynamics
{
	double Re_mf{};
	double u_mf{};
	double u_0{};
	double theta{};

	SFluidDynamics() = default;
	SFluidDynamics(double _re_mf, double _u_mf, double _u_0, double _theta)
		: Re_mf{ _re_mf }
		, u_mf{ _u_mf }
		, u_0{ _u_0 }
		, theta{ _theta }
	{}
};

struct SGeldart
{
	EGeldart type{ EGeldart::AUTO };

	SGeldart() = default;
	SGeldart(EGeldart _geldartType)
		: type{ _geldartType }
	{}

	SGeldart(double _d_p, double _rho_p, double _rho_g)
	{
		if (_rho_p - _rho_g < 0. || (_d_p <= 0. || _rho_p - _rho_g < 1e2 || _rho_p - _rho_g > 1e4))
			type = EGeldart::AUTO;
		else
		{
			const double yl = _rho_p - _rho_g;			// delta_rho = rho_p - rho_g
			const double xl = _d_p * 1e6;				// diameter in um
			const double yc = 3.3 * (log10(yl) - 2.);	// Rho_p-Rho_g, cartesian
			const double xc = 3.3 * (log10(xl) - 1.);	// diameter, cartesian
			//Geldart Type C
			if (xc < 2.3 && yc < -2.5652 * xc + 5.9)
				type = EGeldart::C;
			//Geldart Type A
			else if (xc < 4.4 && yc >= -2.5652 * xc + 5.9 && yc < -3. * xc + 13.2)
				type = EGeldart::A;
			else if (xc < 4.4 && yc >= -3. * xc + 13.2)
			{
				//Geldart Type B
				if (yc < 6.4 && yc < 9.7 - xc)
					type = EGeldart::B;
				//Geldart Type C
				else
					type = EGeldart::C;
			}
			else if (xc >= 4.4)
			{
				if (yc < 6.4 && xc < 7.3)
				{
					//Geldart Type B
					if (yc < 9.7 - xc)
						type = EGeldart::B;
					//Geldart Type D
					else
						type = EGeldart::D;
				}
				//Geldart Type D
				else
					type = EGeldart::D;
			}
			else
				type = EGeldart::AUTO;
		}
	}
};

struct STDParameters
{
	SMYTP inGas{};
	SMYTP inParticle{};
	SFluidDynamics fluidDynamics{};
	SGeldart geldart{};

	STDParameters() = default;
	STDParameters(const SMYTP& _inGas, const SMYTP& _inParticle, const SFluidDynamics& _fluidDynamics, const SGeldart& _geldart)
		: inGas{ _inGas }
		, inParticle{ _inParticle }
		, fluidDynamics{ _fluidDynamics }
		, geldart{ _geldart }
	{}
};

struct SDiscretizationParameters
{
	size_t heightGridPoints{};		// Number of grid points of height discretization. [-]
	size_t RTGridPoints{};			// Number of grid points of residence time discretization. [-]
	size_t sizeClasses{};			// Number of particle size classes. [-]
	size_t moistureClasses{};		// Number of particle moisture classes. [-]
	size_t FDHeightGridPoints{};	// Number of grid points of height discretization for FD Model. [-]

	SDiscretizationParameters() = default;
	SDiscretizationParameters(size_t _heightGridPoints, size_t _rtGridPoints, size_t _sizeClasses, size_t _moistureClasses, size_t _fdHeightGridPoints)
		: heightGridPoints{ _heightGridPoints }
		, RTGridPoints{ _rtGridPoints }
		, sizeClasses{ _sizeClasses }
		, moistureClasses{ _moistureClasses }
		, FDHeightGridPoints{ _fdHeightGridPoints }
	{}
};

struct SBoundaryCondition
{
	double T_g{ 0 };
	double h_g{ 0 };
	double Y_g{ 0 };
	double T_p{ 0 };
	std::vector<double> h_p;
	std::vector<double> X_p;

	SBoundaryCondition() = default;
	SBoundaryCondition(double _t_g, double _h_g, double _y_g, double _t_p, const std::vector<double> _h_p, const std::vector<double>& _x_p)
		: T_g{ _t_g }
		, h_g{ _h_g }
		, Y_g{ _y_g }
		, T_p{ _t_p }
		, h_p{ _h_p }
		, X_p{ _x_p }
	{}
};

struct SResultsNLModelFD
{
	bool valid{ false };
	double h_fb{ 0 };
	std::vector<double> h_z;
	std::vector<double> d_v_z;
	std::vector<double> nu_z;
	std::vector<double> zeta_z;

	SResultsNLModelFD() = default;
	SResultsNLModelFD(double _h_fb, const std::vector<double>& _h_z, std::vector<double> _d_v_z, std::vector<double> _nu_z)
		: valid{ true }
		, h_fb{ _h_fb }
		, h_z{ _h_z }
		, d_v_z{ std::move(_d_v_z) }
		, nu_z{ std::move(_nu_z) }
	{
		zeta_z.resize(_h_z.size());
		for (size_t z = 0; z < zeta_z.size(); ++z)
			zeta_z[z] = _h_z[z] / _h_fb;
	}
};

struct SResultsNLModelTD
{
	bool valid{ false };
	std::vector<std::vector<std::vector<double>>> X_p_dfi;
	std::vector<std::vector<std::vector<double>>> h_p_dfi;
	std::vector<std::vector<std::vector<double>>> T_p_dfi;
	double T_p_mean{};
	std::vector<double> Y_s_z;
	std::vector<double> h_s_z;
	std::vector<double> T_s_z;
	std::vector<double> Y_b_z;
	std::vector<double> h_b_z;
	std::vector<double> T_b_z;
	double T_g_out{};
	double X_p_out{};
	double T_p_out{};
	std::vector<std::vector<double>> alpha_ps_dz;
	double Q_pw{};
	double BalanceH2ODeviation{};
	double BalanceEnthalpyDeviation{};

	SResultsNLModelTD() = default;
	SResultsNLModelTD(std::vector<std::vector<std::vector<double>>> _x_p_dfi, std::vector<std::vector<std::vector<double>>> _h_p_dfi, std::vector<std::vector<std::vector<double>>> _t_p_dfi, double _t_p_mean, std::vector<double> _y_s_z, std::vector<double> _h_s_z, std::vector<double> _t_s_z, std::vector<double> _y_b_z, std::vector<double> _h_b_z, std::vector<double> _t_b_z, double _t_g_out, double _x_p_out, double _t_p_out, std::vector<std::vector<double>> _alpha_ps_dz, double _q_pw, double _balanceH2ODeviation, double _balanceEnthalpyDeviation)
		: valid{ true }
		, X_p_dfi{ std::move(_x_p_dfi) }
		, h_p_dfi{ std::move(_h_p_dfi) }
		, T_p_dfi{ std::move(_t_p_dfi) }
		, T_p_mean{ _t_p_mean }
		, Y_s_z{ std::move(_y_s_z) }
		, h_s_z{ std::move(_h_s_z) }
		, T_s_z{ std::move(_t_s_z) }
		, Y_b_z{ std::move(_y_b_z) }
		, h_b_z{ std::move(_h_b_z) }
		, T_b_z{ std::move(_t_b_z) }
		, T_g_out{ _t_g_out }
		, X_p_out{ _x_p_out }
		, T_p_out{ _t_p_out }
		, alpha_ps_dz{ std::move(_alpha_ps_dz) }
		, Q_pw{ _q_pw }
		, BalanceH2ODeviation{ _balanceH2ODeviation }
		, BalanceEnthalpyDeviation{ _balanceEnthalpyDeviation }
	{}
};

struct SSettings
{
	std::string CSVPath{};
	double C_Qpp{};
	double C_Qw{};
	double C_Ntaumax{ 0.99 };
	EArea C_ApAbed{};
	double C_Ktis{ 1.0 };
	bool HT{ true };
	bool MT{ true };
	double C_H_NTU{ 0.05 };
	double C_BubbleVolFraction{ 2 };

	SSettings() = default;
	SSettings(std::string _csvPath, double _c_Qpp, double _c_Qw, double _c_Ntaumax, EArea _c_ApAbed, double _c_Ktis, bool _isHT, bool _isMT, double _c_H_NTU, double _c_BubbleVolFraction)
		: CSVPath{ std::move(_csvPath) }
		, C_Qpp{ _c_Qpp }
		, C_Qw{ _c_Qw }
		, C_Ntaumax{ _c_Ntaumax }
		, C_ApAbed{ _c_ApAbed }
		, C_Ktis{ _c_Ktis }
		, HT{ _isHT }
		, MT{ _isMT }
		, C_H_NTU{ _c_H_NTU }
		, C_BubbleVolFraction{ _c_BubbleVolFraction }
	{}
};

/// Thermodynamics-model ///
struct STDModelParameters
{
	std::vector<double> ApAbed_d;
	double dp_sauter{};
	double eps_fb{};
	double tau_max{};
	double h_fb{};
	std::vector<double> h_z;
	std::vector<double> zeta_z;
	double nu_0{};
	std::vector<double> tau_i;
	std::vector<double> n_i;
	std::vector<std::vector<double>> mfrac_df;
	std::vector<double> d_d;
	std::vector<std::vector<double>> Aps_df;

	STDModelParameters() = default;
	STDModelParameters(std::vector<double> _apAbed, double _dp_sauter, double _eps_fb, double _tau_max, double _h_fb, const std::vector<double>& _h_z, const std::vector<double>& _zeta_z, const double _nu_0, const std::vector<double>& _tau_i, const std::vector<double>& _n_i, const std::vector<std::vector<double>>& _mfrac_df, const std::vector<double>& _d_d, const std::vector<std::vector<double>>& _aps_df)
		: ApAbed_d{ std::move(_apAbed) }
		, dp_sauter{ _dp_sauter }
		, eps_fb{ _eps_fb }
		, tau_max{ _tau_max }
		, h_fb{ _h_fb }
		, h_z{ _h_z }
		, zeta_z{ _zeta_z }
		, nu_0{ _nu_0 }
		, tau_i{ _tau_i }
		, n_i{ _n_i }
		, mfrac_df{ _mfrac_df }
		, d_d{ _d_d }
		, Aps_df{ _aps_df }
	{}
};
