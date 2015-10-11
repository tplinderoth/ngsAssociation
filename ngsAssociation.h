/*
* ngsAssociation.h
*/

#ifndef NGSASSOCIATION_H_
#define NGSASSOCIATION_H_

#include <fstream>
#include "parsePileup.h"
#include "ArgParser.h"

const char* version = "0.1.4";

// FUNCTION PROTOTYPES

int processPileup (ArgParser* arg, std::string analysis);
unsigned int numCovered (const Pileup* data, const unsigned int mincov);
int doAssoc(Pileup* pile, std::vector< std::pair<std::string,double> >* treatment, std::ostream& os, int v = -1);
int doSummary (Pileup* pile, std::ostream& os, bool indivData, int v = -1);
void exitMessage (int status);

#endif /* NGSASSOCIATION_H_ */
