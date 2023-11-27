/* Copyright (c) 2020, Dyssol Development Team. All rights reserved. This file is part of Dyssol. See LICENSE file for license information. */

#pragma once

#include "UnitDevelopmentDefines.h"

class CUnitDAEModel : public CDAEModel
{
public:
	/// Important constants: later, can be read from Dyssol materials databank ///
	const double rhoGas = 1;
	const double rhoWater = 1000;
	const double rhoParticle = 2500;
	const double beta_GP = 0.02; // Water mass transfer coefficient from gas to particle in [m/s]
	
	/// Indices of state variables for solver ///
	size_t m_iX{};	// Moisture content in the bed in [kg/kg]
	size_t m_iMFlowVapor{};	// Mass flow of water vapor in the gas going to particles
	size_t m_iEtaPFTR{};	// Efficiency of water mass transfer of the gas
	size_t m_iNTU{}; // Number of transfer units of water mass transfer from particle to gas
	size_t m_iPhi{}; // Degree of wetness of the particle
	size_t m_iMWaterBed{}; // Water mass in the bed
	size_t m_iYOut{}; // Moisture content of exhaust gas in [kg/kg]
	
public:
	void CalculateResiduals(double _time, double* _vars, double* _ders, double* _res, void* _unit) override;
	void ResultsHandler(double _time, double* _vars, double* _ders, void* _unit) override;
};

class CDryerSimpleBatch : public CDynamicUnit
{
private:
	CUnitDAEModel m_model;
	CDAESolver m_solver;

public:
	CHoldup* m_holdup{};				// Holdup
	CMaterialStream* m_inWaterStream{};	// Input of water stream
	CMaterialStream* m_inGasStream{};	// Input gas stream
	CMaterialStream* m_outExhaustGasStream{};	// Output of exhaust gas

	/*
	size_t m_classesNum{};				// Number of classes for PSD
	std::vector<double> m_sizeGrid;		// Size grid for PSD
	std::vector<double> m_avgDiam;		// Average values of size grid for PSD
	std::vector<double> m_classSize;	// Class sizes of size grid for PSD
	double m_initMass{};				// Initial mass in the Granulator
	std::vector<double> m_diamRatio;	// Vector: stores ratio of two adjacent diameter values, for calculating G
	std::vector<double> m_massFrac;		// Mass fraction of classes
	*/

public:
	void CreateBasicInfo() override;
	void CreateStructure() override;
	void Initialize(double _time) override;
	void SaveState() override;
	void LoadState() override;
	void Simulate(double _timeBeg, double _timeEnd) override;
};
