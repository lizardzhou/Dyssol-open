/* Copyright (c) 2020, Dyssol Development Team. All rights reserved. This file is part of Dyssol. See LICENSE file for license information. */

#pragma once

#include "UnitDevelopmentDefines.h"

class CHumidifier : public CSteadyStateUnit
{
	CConstRealUnitParameter* m_upYg{};	// Unit parameter - Yg.
	CCompoundUnitParameter* m_upAir{};	// Unit parameter - Air.
	CCompoundUnitParameter* m_upH2Og{}; // Unit parameter - H2Og.

public:
	void CreateBasicInfo() override;
	void CreateStructure() override;
	void Initialize(double _time) override;
	void Simulate(double _time) override;
};
