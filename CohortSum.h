/*
 * CohortSum.h
 * Created on: May 7 2018
 * Author: Mike Bithell
 *
 *   Based on Infectionsum.h
 * Repast for High Performance Computing (Repast HPC)
 *
 *   Copyright (c) 2010 Argonne National Laboratory
 *   All rights reserved.
 *  
 *   Redistribution and use in source and binary forms, with 
 *   or without modification, are permitted provided that the following 
 *   conditions are met:
 *  
 *  	 Redistributions of source code must retain the above copyright notice,
 *  	 this list of conditions and the following disclaimer.
 *  
 *  	 Redistributions in binary form must reproduce the above copyright notice,
 *  	 this list of conditions and the following disclaimer in the documentation
 *  	 and/or other materials provided with the distribution.
 *  
 *  	 Neither the name of the Argonne National Laboratory nor the names of its
 *     contributors may be used to endorse or promote products derived from
 *     this software without specific prior written permission.
 *  
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE TRUSTEES OR
 *   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 *   EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *  Created on: Sep 2, 2010
 *      Author: nick
 */

#ifndef COHORTSUM_H_
#define COHORTSUM_H_

#include "repast_hpc/TDataSource.h"
class MadModel;
//----------------------------------------------------------
class CohortSum: public repast::TDataSource<int> {

private:
	MadModel* obs;

public:
	CohortSum(MadModel* Mobs);
	virtual ~CohortSum();

	int getData();
};
//----------------------------------------------------------
class StockSum: public repast::TDataSource<int> {

private:
	MadModel* obs;

public:
	StockSum(MadModel* Mobs);
	virtual ~StockSum();

	int getData();
};
//----------------------------------------------------------
class CohortAbundanceSum: public repast::TDataSource<double> {

private:
	MadModel* obs;

public:
	CohortAbundanceSum(MadModel* Mobs);
	virtual ~CohortAbundanceSum();

	double getData();
};
//----------------------------------------------------------
class CohortBiomassSum: public repast::TDataSource<double> {

private:
	MadModel* obs;

public:
	CohortBiomassSum(MadModel* Mobs);
	virtual ~CohortBiomassSum();

	double getData();
};
//----------------------------------------------------------
class StockBiomassSum: public repast::TDataSource<double> {

private:
	MadModel* obs;

public:
	StockBiomassSum(MadModel* Mobs);
	virtual ~StockBiomassSum();

	double getData();
};
//----------------------------------------------------------
class CohortOrganicPool: public repast::TDataSource<double> {

private:
	MadModel* obs;

public:
	CohortOrganicPool(MadModel* Mobs);
	virtual ~CohortOrganicPool();

	double getData();
};
//----------------------------------------------------------
class CohortResp: public repast::TDataSource<double> {

private:
	MadModel* obs;

public:
	CohortResp(MadModel* Mobs);
	virtual ~CohortResp();

	double getData();
};
//----------------------------------------------------------
class DispersalSum: public repast::TDataSource<double> {
private:
	MadModel* obs;

public:
    DispersalSum(MadModel* Mobs);
	~DispersalSum();

	double getData();
};
//----------------------------------------------------------
class ExtinctionSum: public repast::TDataSource<double> {
private:
	MadModel* obs;

public:
    ExtinctionSum(MadModel* Mobs);
	~ExtinctionSum();

	double getData();
};
//----------------------------------------------------------
class ProductionSum: public repast::TDataSource<double> {
private:
	MadModel* obs;

public:
    ProductionSum(MadModel* Mobs);
	~ProductionSum();

	double getData();
};
//----------------------------------------------------------
class CombinationSum: public repast::TDataSource<double> {
private:
	MadModel* obs;

public:
    CombinationSum(MadModel* Mobs);
	~CombinationSum();

	double getData();
};
#endif /* COHORTSUM_H_ */
