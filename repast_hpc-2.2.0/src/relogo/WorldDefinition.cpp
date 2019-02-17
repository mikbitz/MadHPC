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
 *  WorldDefinition.cpp
 *
 *  Created on: Aug 5, 2010
 *      Author: nick
 */

#include "WorldDefinition.h"

#include "RelogoDiscreteSpaceAdder.h"
#include "RelogoContinuousSpaceAdder.h"
#include "relogo.h"

#include <iostream>

namespace repast {

namespace relogo {

WorldDefinition::WorldDefinition(int minX, int minY, int maxX, int maxY, bool wrapped, int buffer) :
	_dims(GridDimensions(Point<double> ((double)minX, (double)minY), Point<double> ((double)(maxX - minX + 1) , (double)(maxY - minY + 1)))), _wrapped(wrapped),
			_buffer(buffer)  {
}

WorldDefinition::~WorldDefinition() {
}

void WorldDefinition::defineNetwork(std::string name, bool directed, RelogoLinkContentManager* rlcm) {
	repast::SharedNetwork<RelogoAgent, RelogoLink, RelogoLinkContent, RelogoLinkContentManager>* net =
			new repast::SharedNetwork<RelogoAgent, RelogoLink, RelogoLinkContent, RelogoLinkContentManager>(name, directed, rlcm);
	networks.push_back(net);
}

void WorldDefinition::defineNetwork(bool directed, RelogoLinkContentManager* rlcm) {
	std::string name = directed ? DEFAULT_DIR_NET : DEFAULT_UNDIR_NET;
	defineNetwork(name, directed, rlcm);
}

}
}
