/* Copyright (c) 2020, Dyssol Development Team. All rights reserved. This file is part of Dyssol. See LICENSE file for license information. */

#define DLL_EXPORT
#include "Humidifier.h"

extern "C" DECLDIR CBaseUnit* DYSSOL_CREATE_MODEL_FUN()
{
	return new CHumidifier();
}

//////////////////////////////////////////////////////////////////////////
/// Unit

void CHumidifier::CreateBasicInfo()
{
	/// Basic unit's info ///
	SetUnitName  ("Humidifier");
	SetAuthorName("Zhengyu Lu");
	SetUniqueID  ("9EAF72A9228A4DED968B693B21E14F23");

}

void CHumidifier::CreateStructure()
{
	/// Add ports ///
	AddPort("InGas" , EUnitPort::INPUT);
	AddPort("OutGas", EUnitPort::OUTPUT);

	/// Add unit parameters ///
	m_upYg   = AddConstRealParameter("Yg", 0, "g/kg", "Absolute humidity, dry basis.\nGram H2O per kg dry air", 0);
	m_upAir  = AddCompoundParameter("Air", "Compound to treat as air");
	m_upH2Og = AddCompoundParameter("Water vapor", "Compound to treat as water vapor");
}

void CHumidifier::Initialize(double _time)
{
	if (!IsPhaseDefined(EPhase::VAPOR))
		RaiseWarning("Vapor phase has not been defined.");
}

void CHumidifier::Simulate(double _time)
{
	const CMaterialStream* inStream = GetPortStream("InGas");
	CMaterialStream* outStream = GetPortStream("OutGas");

	outStream->CopyFromStream(_time, inStream);

	const std::string keyAir   = m_upAir->GetCompound();	// Compound key of air
	const std::string keyH2O_g = m_upH2Og->GetCompound();	// Compound key of H2O_g

	const double targetYg = GetConstRealParameterValue("Yg");

	const double massAir   = inStream->GetCompoundMassFlow(_time, keyAir, EPhase::VAPOR);
	const double massH2O_g = massAir * targetYg * 1e-3;		// Transfer the unit of mass of vapor from g to kg
	const double massTotal = massAir + massH2O_g;

	const double massFracAir   = massAir / massTotal;
	const double massFracH2O_g = massH2O_g / massTotal;

	outStream->SetMassFlow(_time, massTotal);
	outStream->SetCompoundFraction(_time, keyAir  , EPhase::VAPOR, massFracAir);
	outStream->SetCompoundFraction(_time, keyH2O_g, EPhase::VAPOR, massFracH2O_g);
	for (const auto& key : GetAllCompounds())	// Set the mass fractions of the rest compounds to 0.
		if (key != keyAir && key != keyH2O_g)
			outStream->SetCompoundFraction(_time, key, EPhase::VAPOR, 0.0);
}
