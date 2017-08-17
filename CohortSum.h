/*
 *   Repast for High Performance Computing (Repast HPC)
 *
 *   Copyright (c) 2010 Argonne National Laboratory
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with
 *   or without modification, are permitted provided that the following
 *   conditions are met:
 *
 *     Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *
 *     Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *
 *     Neither the name of the Argonne National Laboratory nor the names of its
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
 *
 *
 * InfectionSum.h
 *
 *  Created on: Oct 13, 2010
 *      Author: nick
 */

#ifndef COHORTSUM_H_
#define COHORTSUM_H_

#include "repast_hpc/TDataSource.h"

class MadObserver;
//----------------------------------------------------------
class CohortSum: public repast::TDataSource<int> {

private:
	MadObserver* obs;

public:
	CohortSum(MadObserver* Mobs);
	virtual ~CohortSum();

	int getData();
};
//----------------------------------------------------------
class CohortAbundanceSum: public repast::TDataSource<double> {

private:
	MadObserver* obs;

public:
	CohortAbundanceSum(MadObserver* Mobs);
	virtual ~CohortAbundanceSum();

	double getData();
};
//----------------------------------------------------------
class CohortBiomassSum: public repast::TDataSource<double> {

private:
	MadObserver* obs;

public:
	CohortBiomassSum(MadObserver* Mobs);
	virtual ~CohortBiomassSum();

	double getData();
};
//----------------------------------------------------------
class StockBiomassSum: public repast::TDataSource<double> {

private:
	MadObserver* obs;

public:
	StockBiomassSum(MadObserver* Mobs);
	virtual ~StockBiomassSum();

	double getData();
};

#endif /* COHORTSUM_H_ */
