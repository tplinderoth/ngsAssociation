/*
*
* ngsAssociation - Association mapping and SNP calling using pooled or unpooled NGS data
* Copyright (C) 2015 Tyler P. Linderoth

* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.

* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.

* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
* Tyler's contact information: tylerp.linderoth@gmail.com
*
*/

#include <cstdio>
#include <iomanip>
#include "ngsAssociation.h"
#include "generalUtils.h"
#include "SiteStat.h"

int main (int argc, char** argv)
{
	// define variables
	int rc=0;
	std::string runtype;
	ArgParser runpar;

	// parse user input
	if((rc=runpar.parseInput(argc, argv, version)) != 0)
	{
		if (rc < 0)
		{
			exitMessage(rc);
			return rc;
		}
		else
			return 0;
	}
	else
	{
		if((rc=runpar.setStreams(runpar.inpileup_name(), runpar.outfile_name())))
		{
			exitMessage(rc);
			return rc;
		}
		else
			runtype=argv[1];
	}

	// process pileup file
	if(processPileup(&runpar,runtype))
	{
		exitMessage(1);
		return 0;
	}

	// exit program
	exitMessage(0);

	return 0;
}


int processPileup (ArgParser* arg, std::string analysis)
{
	// define variables
	int status=0;
	int verbose=-1;
	int type;
	std::string(line);
	std::vector< std::pair<std::string,double> > treat;
	unsigned int minpool=arg->minpool();
	unsigned int mincov=arg->mincov();

	// initialize pileup object
	Pileup pile;
	pile.setQualCode(arg->offset());
	pile.setMinQ(arg->minQ());
	pile.setpoolsz(arg->poolsz());
	getline(arg->input(), line);
	pile.setn(line);

	// check that there are >= minpool individuals in dataset
	if (arg->minpool() > pile.nInd())
	{
		fprintf(stderr,"Fewer than %u individuals in dataset --> decrease -minpooln\n",arg->minpool());
		return 1;
	}

	// set variables for specified analysis
	if (analysis == "association")
	{
		type=0;
		status=ArgParser::parseTreatment(arg->treatmentf_name(), &treat, &pile, arg->lrmethod());
		if (status)
			return status;
	}
	else if (analysis == "summarize")
		type=1;
	else
		type=-1;

	// loop through the input
	while (!line.empty())
	{
		// store data in pileup object
		pile.getSeqDat(line);
		if (pile.fail())
		{
			fprintf(stderr,"Error parsing sequencing data --> check pileup\n");
			return 1;
		}
		// perform analysis if site meets coverage requirements
		if (numCovered(&pile, mincov) >= minpool)
		{
			if (type == 0)
				status=doAssoc(&pile, &treat, arg->output(), arg->lrmethod(), verbose);
			else if (type == 1)
				status=doSummary(&pile, arg->output(), arg->printInd(), verbose);
			else
			{
				fprintf(stderr,"Unrecognized command in processPileup: %s\n",analysis.c_str());
				return 1;
			}
		}

		// return control for analysis
		if (status)
		{
			switch(status)
			{
				case -1 :
					fprintf(stderr, "Skipping site ...\n");
					break;
				default :
					fprintf(stderr,"Unable to carry out %s analysis",analysis.c_str());
					return 1;
			}
		}

		// fetch next line
		getline(arg->input(),line);
	}

	return 0;
}

unsigned int numCovered (const Pileup* data, const unsigned int mincov)
{
	static std::vector<SiteData>::const_iterator indIter;
	unsigned int n = 0;
	for (indIter = data->seqdat.begin(); indIter != data->seqdat.end(); ++indIter)
	{
		if (indIter->cov() >= mincov)
			++n;
	}
	return n;
}


int doAssoc(Pileup* pile, std::vector<treatdat>* treatment, std::ostream& os, int lrmethod, int v)
{
	int rc=0;
	static double like [2]; /* [null, alternative] */
	static double lr;
	static std::vector<treatdat>::const_iterator tIter;
	int j=0;
	const static int prec=8;
	const static double thresh=-1e-6; /* round negative LR to zero if greater than this */
	int offset; /* used to set which MAFs to print */

	// find MAFs and calculate LR
	if (lrmethod == 1)
	{
		lr = SiteStat::assoclrt1(pile, treatment, like, &rc, v);
		offset = 1;
	}
	else if (lrmethod == 2)
	{
		lr = SiteStat::assoclrt2(pile, treatment, like, &rc, v);
		offset = 0;
	}
	else
	{
		fprintf(stderr, "Invalid LR calculation method %i passed to doAssoc subroutine\n", lrmethod);
		return -1;
	}
	if (rc)
		return rc;

	// print result
	if (lr < 0.0)
	{
		if (lr >= thresh)
			lr = 0.0; // round up to zero
	}
	os << pile->seqName() << "\t" << pile->position(); // chromosome and position
	for (j=0; j<2; ++j)
		os << "\t" << std::fixed << std::setprecision(prec) << like[j]; // likelihoods
	os << "\t" << std::fixed << std::setprecision(prec) << lr; // LR
	for (tIter=treatment->begin()+offset; tIter!=treatment->end(); ++tIter)
		os << "\t" << std::fixed << std::setprecision(prec) << tIter->second; // MAFs
	os << "\n";

	return rc;
}

int doSummary (Pileup* pile, std::ostream& os, bool indivData, int v)
{
	int rc=0;
	static double lr;
	static double maf;
	static std::vector<SiteData>::const_iterator poolIter;
	static std::vector<seqread>::const_iterator readIter;
	static unsigned int poolcov;
	const static int prec=8;
	const static int qprec=1;
	const static double thresh=-1e-6; /* round negative LR to zero if greater than this */

	// find MAF and calculate LR
	lr = SiteStat::snplr(pile, &maf, &rc, v);
	if (rc)
		return rc;

	// print summary information
	if (lr < 0.0)
	{
		if (lr >= thresh)
			lr = 0.0; // round up to zero
	}
	os << pile->seqName()
	<< "\t" << pile->position()
	<< "\t" << pile->refAllele()
	<< "\t" << pile->siteDepth() << ";" << pile->refcount() << ";" << pile->altcount()
	<< "\t" << pile->alleleCount('A') << ";" << pile->alleleCount('C') << ";" << pile->alleleCount('G') << ";" << pile->alleleCount('T') << ";" << pile->alleleCount('I')
	<< "\t" << std::fixed << std::setprecision(prec) << maf
	<< "\t" << std::fixed << std::setprecision(prec) << lr;

	// print individual pool reads and quality scores
	if (indivData)
	{
		// loop over all individuals
		for (poolIter = pile->seqdat.begin(); poolIter != pile->seqdat.end(); ++poolIter)
		{
			os << "\t" << poolIter->cov() << "\t";
			poolcov=poolIter->cov();
			if (poolcov > 0)
			{
				// loop over all reads for individual
				for (readIter = poolIter->rdat.begin(); readIter != poolIter->rdat.begin() + poolcov; ++readIter)
					os << readIter->first;
				// loop over all quality scores for individual
				os << "\t";
				for (readIter = poolIter->rdat.begin(); readIter != poolIter->rdat.begin() + poolcov; ++readIter)
				{
					os << std::fixed << std::setprecision(qprec) << readIter->second;
					if (static_cast <unsigned int> (std::distance(readIter, poolIter->rdat.begin() + poolcov)) > 1)
						os << ";";
				}
			}
			else
				os << "*\t*";
		}
	}
	os << "\n";

	return rc;
}

void exitMessage (int status)
{
	if (status)
		std::cerr << "--> exiting\n";
	else
		std::cerr << "Finished!\n";
}
