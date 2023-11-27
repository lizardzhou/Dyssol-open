/* Copyright (c) 2022, Dyssol Development Team. All rights reserved. This file is part of Dyssol. See LICENSE file for license information. */
// Author Alexander Thomas Hanke

#pragma once
#include "UnitDevelopmentDefines.h"

// Physical type definitions to increase readabilty

typedef double length; // m meter

typedef double temperature; // K kelvin

typedef double mass; // kg kilogram

//typedef double time; // s // ToDo - find out where pref. declared as function

typedef double area; // m^2

typedef double volume; // m^3

typedef double density; // kg/m^3

typedef double pressure; // Pa

typedef double dynamicViscosity; // Pa*s

typedef double heatCapacity; // J/(kg*K)

typedef double thermalConductivity; // W/(m*K)

typedef double molarMass; // kg/mol

typedef double moistureContent; // kg (liquid) /kg (dry gas/solid)

typedef double massFraction; // kg/kg

typedef double specificLatentHeat; // J/kg

typedef double massTransferCoefficient; // mol/(s*m2)/(mol/m3)

typedef double dimensionlessNumber; // -

typedef double massFlow; // kg/s

typedef double energyStream; // W

enum class EShape : size_t
{
	CYLINDRICAL = 0, RECTANGULAR = 1
};

struct chamberSection
{
	std::string name;
	EShape shape;
	std::vector<std::pair<double,double>> dimensionsInternal;
	double height;
	std::vector<double> wallThicknesses;
	std::vector<double> thermalConductivities;
	size_t layers;
};
struct equilibriumMoistureContentData
{
	std::string compoundKey;
	std::set<double> RHs;
	std::set<double> temperatures;
	std::map< std::pair<double, double>, double> equilibriumMoistureContents;
	std::map<double, std::set<double>::iterator> maxRH;
};