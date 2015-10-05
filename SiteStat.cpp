/*
* SiteStat.cpp
*/

#include <vector>
#include <math.h>
#include "SiteStat.h"
#include "randomVar.h"
#include "Matrix.h"
#include "generalUtils.h"
#include "bfgs.h"

// FUNCTION DEFINITIONS

double SiteStat::assoclrt (Pileup* pile, std::vector<treatdat>* treatment, double like [], int* status, int verb)
{
	/*
	 * returns the likelihood ratio that alleles are differentially associated among treatments
	 * seqdata: sequencing data
	 * treatment: different treatment identifiers
	 * like: stores -log likelihoods of [null, alternative]
	 */

	const int npars = 1; // number of parameters to optimize
	static double inval [npars]; // starting maf value
	static int nbounds [] = {2}; // signifies boundary conditions (see bfgs.h for description)
	static double lowbound [] = {0.0}; // lower bound on maf
	static double upbound [npars]; // upper bound on maf
	static std::vector<treatdat>::iterator tIter;
	double altlike, lr;

	// parameters for iterative optimization
	const double thresh=-1e-6; // cutoff for determining optimization failure
	static double nullstart [] = {0.0001, 0.25, 0.49}; // null model starting points
	static double altstart [] = {0.0, 0.0001, 0.5, 0.99}; // alternative model starting points
	const int nullpoints = sizeof(nullstart)/sizeof(double); // number of null optimization start points
	const int altpoints = sizeof(altstart)/sizeof(double); // number of alternative optimization start points

	// set major allele
	pile->setMajor(pile->empiricalMajor()); // using most common allele as major

	// calculate null MAF
	upbound[0] = 0.5;
	pile->assignTreatment("");
	inval[0]=mafguess(pile);
	like[0] = findmax_bfgs(npars, inval, pile, maflike, NULL, lowbound, upbound, nbounds, verb, status);

	if (*status)
	{
		*status=0;
		like[0]=multiOptim(pile, nullstart, nullpoints, inval, lowbound, upbound, nbounds, like[0], status, verb);
		if (*status)
		{
			optimErrorMsg(pile);
			return -1.0;
		}
	}

	for(tIter = treatment->begin(); tIter != treatment->end(); ++tIter)
	{
		if (tIter->first == "null")
		{
			tIter->second=altstart[0]=inval[0];
			break;
		}
	}

	// calculate alternative MAFs
	like[1]=0.0;
	upbound[0] = 1.0;
	for(tIter = treatment->begin(); tIter != treatment->end(); ++tIter)
	{
		if (tIter->first != "null")
		{
			pile->assignTreatment(tIter->first);
			inval[0]=altstart[0];
			altlike = findmax_bfgs(npars, inval, pile, maflike, NULL, lowbound, upbound, nbounds, verb, status);
			if (*status)
			{
				*status=0;
				altlike=multiOptim(pile, altstart, altpoints, inval, lowbound, upbound, nbounds, altlike, status, verb);
				if (*status)
				{
					optimErrorMsg(pile);
					return -1.0;
				}
			}
			tIter->second=inval[0];
			like[1] += altlike;
		}
	}

	// calculate LR
	lr = 2*like[0]-2*like[1];

	// check for optimization failure
	if (lr < thresh)
	{
		/* try multiple starting points */

		//null
		upbound[0] = 0.5;
		like[0]=multiOptim(pile, nullstart, nullpoints, inval, lowbound, upbound, nbounds, like[0], status, verb);
		for(tIter = treatment->begin(); tIter != treatment->end(); ++tIter)
		{
			if (tIter->first == "null")
			{
				tIter->second=altstart[0]=inval[0];
				break;
			}
		}
		// check if bfgs routine failed
		if (*status)
		{
			optimErrorMsg(pile);
			return -1.0;
		}

		// alternative
		like[1]=0.0;
		upbound[0] = 1.0;
		for(tIter = treatment->begin(); tIter != treatment->end(); ++tIter)
		{
			if (tIter->first != "null")
			{
				pile->assignTreatment(tIter->first);
				altlike=multiOptim(pile, altstart, altpoints, inval, lowbound, upbound, nbounds, 1.0/0.0, status, verb);
				tIter->second=inval[0];
				like[1] += altlike;
			}
		}
		// check if bfgs routine failed
		if (*status)
		{
			optimErrorMsg(pile);
			return -1.0;
		}

		// recalculate LR
		lr = 2*like[0]-2*like[1];
	}

	return lr;
}

double SiteStat::snplr (Pileup* pile, double* mlmaf, int* status, int verb)
{
	/*
	 * returns the likelihood ratio that a site is variable
	 * can return maximum likelihood estimate of MAF through mlmaf
	 * null model: maf = 0;
	 * alternative model: maf is a value in the range (0.0,0.5]
	 */

	const int npars = 1; // number of parameters to optimize
	static double inval [npars]; // starting maf value
	static int nbounds [] = {2}; // signifies boundary conditions (see bfgs.h for description)
	static double lowbound [] = {0.0}; // lower bound on maf
	static double upbound [] = {0.5}; // upper bound on maf
	double alt, null, lr;
	bool multiopt=false;

	// parameters for iterative optimization
	const double thresh=-1e-6; // cutoff for determining optimization failure
	static const double startval [] = {0.0001, 0.25, 0.49};
	const int npoints = sizeof(startval)/sizeof(double); // number of optimization start points

	// set major allele
	pile->setMajor(pile->empiricalMajor()); // using most common allele as major

	// find maximum likelihood estimate of maf

	inval[0] = mafguess(pile);
	alt = findmax_bfgs(npars, inval, pile, maflike, NULL, lowbound, upbound, nbounds, verb, status);
	if (*status)
	{
		*status=0;
		alt=multiOptim(pile, startval, npoints, inval, lowbound, upbound, nbounds, 1.0/0.0, status, verb);
		multiopt=true;
		if (*status)
		{
			optimErrorMsg(pile);
			return -1.0;
		}
	}
	*mlmaf=inval[0];

	// calculate null case
	inval[0]=0.0;
	null = maflike(inval, pile);

	// calculate LR
	lr = 2*null-2*alt;

	// check for optimization failure
	if (lr < thresh && !multiopt)
	{
		// optimize alternative maf at multiple start points
		alt=multiOptim(pile, startval, npoints, inval, lowbound, upbound, nbounds, 1.0/0.0, status, verb);
		if (*status)
		{
			optimErrorMsg(pile);
			return -1.0;
		}

		// recalculate LR
		lr = 2*null-2*alt;
	}

	return lr;
}

double SiteStat::multiOptim (Pileup* pile, const double* startval, int npoints, double* par, double* lb, double* ub, int* nbounds, double min, int* status, int verb)
{
	const int npars = 1;
	double like;
	double inval;
	int fail;
	double* inptr = &inval;

	*status=0;

	for (int i=0; i<npoints; ++i)
	{
		fail=0;
		inval = startval[i];
		like=findmax_bfgs(npars, inptr, pile, maflike, NULL, lb, ub, nbounds, verb, &fail);
		if (like < min)
		{
			min=like;
			*par=inval;
			*status=fail;
		}
	}
	return min;
}

double SiteStat::mafguess (Pileup* pile)
{
	char major;
	const static char a [] = {'A', 'C', 'G', 'T'};
	int i;
	double c = 0.0;
	double total = 0.0;
	double n;

	major = pile->majorid() ? pile->majorid() : pile->empiricalMajor();

	for (i=0; i<4; ++i)
	{
		n = static_cast<double>(pile->alleleCount(a[i]));
		if (a[i] != major)
			c += n;
		total += n;
	}

	return (c/total);
}

void SiteStat::optimErrorMsg (const Pileup* seqdata)
{
	fprintf(stderr,"Optimization failure for %s %u\n", seqdata->seqName().c_str(), seqdata->position());
}

double SiteStat::maflike (const double par [], const void* data)
{
/*
* adapted from Kim et al. 2010 in Genetic Epidemiology
* instead of using one error rate, the error for each read is considered
*
* seqdata has the reads and quality scores for the site
* maf =  alternate allele frequency
* poolsz = number of individuals each library represents (1 if a single indivdiual, > 1 if pooled)
* err = Prob(sequencing error); if not specified read quality scores from seqdata will be used
*/

	// define variables
	const Pileup* seqdata=static_cast<const Pileup*>(data);
	double maf=par[0];
	unsigned int s = seqdata->poolsz(); /* number of haplotypes used to construct each pool */
	static Array<double> ptmp(s+1, 0.0); /* array for storing P(X_r|k,err)*P(k|maf) */
	//double comp = 0.0; /* compensation for lost low-order bits in Kahan sum */
	static std::vector<SiteData>::const_iterator poolIter;
	static std::vector<seqread>::const_iterator readIter;
	unsigned int k = 0;
	double c = -1.0/0.0; /* scaling factor for underflow protection */
	double ppool = 0.0; /* probability of pool m, P(O_m|maf) */
	double like = 0.0; /* -log P(D_i|maf) */

	// identify major allele (most common read for site)
	char maj=seqdata->majorid();

	// loop over all pools (or individuals if s=2)
	for (poolIter = seqdata->seqdat.begin(); poolIter != seqdata->seqdat.end(); ++poolIter)
	{
		// check for pool missing data and check if pool is in the subset to be analyzed
		if (poolIter->cov() < 1 || (seqdata->treatment()!="" && poolIter->id() != seqdata->treatment()))
			continue;
		for (k = 0; k <= s; ++k)
			ptmp[k] = 0.0;
		c = -1.0/0.0;
        // sum over all possible true counts of minor allele in pool
        for (k = 0; k <= s; ++k)
        {
			// iterate over all reads and use read-specific error rate
			for (readIter = poolIter->rdat.begin(); readIter != poolIter->rdat.begin() + poolIter->cov(); ++readIter)
				ptmp[k] += log(readProb(readIter->first, maj, k, readIter->second, s, MAXQUAL));
				//kahanSum(log(readProb(readIter->first, maj, k, readIter->second, s, MAXQUAL)), &ptmp[k], &comp); /* caused rounding error in optimization */
			// apply the prior on the number of true minor alleles
			ptmp[k] += log(randomVar::binomProb(s,k,maf));
			//kahanSum(log(randomVar::binomProb(s,k,maf)), &ptmp[k], &comp); /* caused rounding error in optimization */
			if (ptmp[k] > c)
				c = ptmp[k];
        }
		// calculate probability of pool, P(O_m|maf), with underflow protection
		ppool = 0.0;
		for (k = 0; k <= s; ++k)
			ppool += exp(ptmp[k] - c);
		like += log(ppool) + c;
	}
	return -like;
}

char SiteStat::findMajor (const Pileup* sitedata)
{
	char maj='A';
	static const char a[] = {'A', 'C', 'G', 'T'};
	for (int i=0; i<4; ++i)
	{
		if (sitedata->alleleCount(a[i]) > sitedata->alleleCount(maj))
			maj=a[i];
	}
	return maj;
}

double SiteStat::readProb (const char obs, const char major, unsigned int nminor, unsigned int qscore, const unsigned int nchr, const unsigned int qmax)
{
	/* for probability matrix rows = # possible true minor alleles in pool, columns = quality score*/
	static Matrix<double> majtab(nchr+1, qmax+1, 0.0); /* stores all possible values for P(X_r|k,err), X_r = major allele */
	static Matrix<double> mintab(nchr+1, qmax+1, 0.0); /* stores all possible values for P(X_r|k,err), X_r = minor allele */
	static Matrix<double>* p = NULL;

	if (nminor <= nchr && qscore <= qmax)
	{
		p = obs == major ? &majtab : &mintab;
		return ((*p)[nminor][qscore] != 0.0 ? (*p)[nminor][qscore] : ((*p)[nminor][qscore]=pread(nminor, qscore, nchr, obs, major)));
	}
	else
		return pread(nminor, qscore, nchr, obs, major);
}

double SiteStat::pread (int k, double qscore, int s, const char obs, const char major)
{
	/*
	* k = number of true minor alleles
	* qscore = phred scaled base quality score
	* s = haploid pool sample size
	* obs = observed allele
	* major = major allele
	*/

	double err = pow(10, -qscore/10);

	if (obs == major)
		return (1.0-err)*(static_cast<double>(s-k)/s) + (err/3)*(static_cast<double>(k)/s);
	else
		return err*(static_cast<double>(s-k)/s) + (1.0 - err + (2*err)/3)*(static_cast<double>(k)/s);
}

void SiteStat::kahanSum(double summand, double* total, double* comp)
{
	double x = summand - *comp;
	double y = *total + x;
	*comp = (y - *total) - x;
	if (*comp != *comp) /* prevent numerical problems when *comp == NaN */
		*comp = 0.0;
	*total = y;
}

#undef MAXQUAL
