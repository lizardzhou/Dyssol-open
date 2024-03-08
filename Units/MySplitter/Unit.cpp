/* Copyright (c) 2020, Dyssol Development Team. All rights reserved. This file is part of Dyssol. See LICENSE file for license information. */

#define DLL_EXPORT
#include "Unit.h"

extern "C" DECLDIR CBaseUnit* DYSSOL_CREATE_MODEL_FUN()
{
	return new CUnit();
}

//////////////////////////////////////////////////////////////////////////
/// Unit

void CUnit::CreateBasicInfo()
{
	/// Basic unit's info ///
	SetUnitName  ("MySplitter");
	SetAuthorName("Mads Jessen");
	SetUniqueID  ("00000000000000004000000020000100"); 
}

void CUnit::CreateStructure()
{
	/// Add ports ///
	AddPort("In" , EUnitPort::INPUT);
	AddPort("Out1", EUnitPort::OUTPUT);
	AddPort("Out2", EUnitPort::OUTPUT);
	AddPort("Out3", EUnitPort::OUTPUT);

	/// Add unit parameters ///
	AddTDParameter       ("ParamTD"    , 0, "kg"        , "Unit parameter description");
	AddConstRealParameter("ParamConst" , 0, "s"         , "Unit parameter description");
	AddStringParameter   ("ParamString", "Initial value", "Unit parameter description");


}

void CUnit::Initialize(double _time)
{
	

}

void CUnit::Simulate(double _time)
{
	CStream* inStream  = GetPortStream("In");
	CStream* outStream1 = GetPortStream("Out1");
	CStream* outStream2 = GetPortStream("Out2");
	CStream* outStream3 = GetPortStream("Out3");

	outStream1->CopyFromStream(_time, inStream);
	outStream2->CopyFromStream(_time, inStream);
	outStream3->CopyFromStream(_time, inStream);


}

void CUnit::Finalize()
{

}
