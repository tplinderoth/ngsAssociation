/*
* SiteStat.cpp
*/

#include <vector>
#include <math.h>
#include <limits>
#include "SiteStat.h"
#include "randomVar.h"
#include "Matrix.h"
#include "generalUtils.h"
#include "bfgs.h"

// FUNCTION DEFINITIONS

double SiteStat::assoclrt2 (Pileup* pile, std::vector<treatdat>* treatment, double like [], int* status, int verb)
{
	/*
	 * returns the likelihood ratio that alleles are differentially associated among treatments
	 * seqdata: sequencing data
	 * treatment: different treatment identifiers
	 * like: stores -log likelihoods of [null, alternative]
	 */

	//verb = 100; // debug

	const int npars = 1; // number of parameters to optimize
	static double inval [npars]; // starting maf value
	static int nbounds [] = {2}; // signifies boundary conditions (see bfgs.h for description)
	static double lowbound [] = {0.0}; // lower bound on maf
	static double upbound [npars]; // upper bound on maf
	static std::vector<treatdat>::iterator tIter;
	double lr;
	static Array<double> altlike(treatment->size()-1, 0.0);
	int i;
	bool weightcount = true;

	// parameters for iterative optimization
	const double thresh=-1e-6; // cutoff for determining optimization failure
	static double nullstart [] = {0.0001, 0.25, 0.49}; // null model starting points
	static double altstart [] = {0.0, 0.0001, 0.5, 0.99}; // alternative model starting points, first value is a dummy
	const int nullpoints = sizeof(nullstart)/sizeof(double); // number of null optimization start points
	const int altpoints = sizeof(altstart)/sizeof(double); // number of alternative optimization start points

	// set major and minor allele
	pile->setMajor(pile->empiricalMajor(weightcount)); // using most common allele as major
	pile->setMinor(pile->empiricalMinorFast(weightcount)); // using second most common allele

	// calculate null MAF
	upbound[0] = 0.5;
	pile->assignTreatment("");
	inval[0]=mafguess(pile, weightcount);
	like[0] = findmax_bfgs(npars, inval, pile, maflike, NULL, lowbound, upbound, nbounds, verb, status);

	if (*status)
	{
		*status=0;
		like[0]=multiOptim(pile, nullstart, nullpoints, inval, lowbound, upbound, nbounds, 1.0/0.0, status, verb);
		if (*status)
		{
			optimErrorMsg(pile);
			*status=-1.0;
			return std::numeric_limits<double>::quiet_NaN();
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
	i = 0;
	for(tIter = treatment->begin(); tIter != treatment->end(); ++tIter)
	{
		if (tIter->first != "null")
		{
			pile->assignTreatment(tIter->first);
			//inval[0]=altstart[0]; // this would be conservative
			inval[0]= treatfreqFast(pile, &tIter->first, pile->minorid());
			altlike[i] = findmax_bfgs(npars, inval, pile, maflike, NULL, lowbound, upbound, nbounds, verb, status);
			if (*status)
			{
				*status=0;
				altlike[i]=multiOptim(pile, altstart, altpoints, inval, lowbound, upbound, nbounds, 1.0/0.0, status, verb);
				if (*status)
				{
					optimErrorMsg(pile);
					*status=-1.0;
					return std::numeric_limits<double>::quiet_NaN();
				}
			}
			tIter->second=inval[0];
			like[1] += altlike[i];
			++i;
		}
	}

	// calculate LR
	lr = calclr(&like[0], &like[1]);

	// check for optimization failure
	if (lr < thresh)
	{
		/* try multiple starting points */

		//null
		upbound[0] = 0.5;
		pile->assignTreatment("");
		like[0]=multiOptim(pile, nullstart, nullpoints, inval, lowbound, upbound, nbounds, 1.0/0.0, status, verb);
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
			*status=-1.0;
			return std::numeric_limits<double>::quiet_NaN();
		}

		// alternative
		like[1]=0.0;
		upbound[0] = 1.0;
		i=0;
		for(tIter = treatment->begin(); tIter != treatment->end(); ++tIter)
		{
			if (tIter->first != "null")
			{
				pile->assignTreatment(tIter->first);
				altlike[i]=multiOptim(pile, altstart, altpoints, inval, lowbound, upbound, nbounds, 1.0/0.0, status, verb);
				tIter->second=inval[0];
				like[1] += altlike[i];
				++i;
			}
		}
		// check if bfgs routine failed
		if (*status)
		{
			optimErrorMsg(pile);
			*status=-1.0;
			return std::numeric_limits<double>::quiet_NaN();
		}

		// recalculate LR
		lr = calclr(&like[0], &like[1]);
		// check if optimization still failed
		if (lr < thresh)
		{
			optimErrorMsg(pile);
			*status=-1.0;
			return std::numeric_limits<double>::quiet_NaN();
		}
	}

	return lr;
}


double SiteStat::assoclrt1 (Pileup* pile, std::vector<treatdat>* treatment, double like [], int* status, int verb)
{
	/*
	 * returns the likelihood ratio that alleles are differentially associated among treatments
	 * seqdata: sequencing data
	 * treatment: different treatment identifiers
	 * like: stores -log likelihoods of [null, alternative]
	*/

	//verb = 100; // debug

	const int npars = 1; // number of parameters to optimize
	static double inval [npars]; // starting maf value
	static int nbounds [] = {2}; // signifies boundary conditions (see bfgs.h for description)
	static double lowbound [] = {0.0}; // lower bound on maf
	static double upbound [npars] = {1.0}; // upper bound on maf
	static std::vector<treatdat>::iterator tIter;
	double nullfreq, lr;
	bool weightcount = true;

	// parameters for iterative optimization
	const double thresh=-1e-6; // cutoff for determining optimization failure
	static double start [] = {0.1, 0.0001, 0.5, 0.99}; // optimization starting points, first value is a dummy
	const int startpoints = sizeof(start)/sizeof(double); // number of optimization start points

	// set major and minor allele
	pile->setMajor(pile->empiricalMajor(weightcount)); // using most common allele as major
	pile->setMinor(pile->empiricalMinorFast(weightcount)); // using second most common allele

	// treatment 1 MAF
	std::string t1id;
	like[0] = 0.0;
	tIter = treatment->begin();
	while(tIter->first == "null") ++tIter;
	t1id=tIter->first;
	pile->assignTreatment(t1id);
	inval[0]=treatfreqFast(pile, &t1id, pile->minorid());

	findmax_bfgs(npars, inval, pile, maflike, NULL, lowbound, upbound, nbounds, verb, status);
	if (*status)
	{
		*status=0;
		multiOptim(pile, start, startpoints, inval, lowbound, upbound, nbounds, 1.0/0.0, status, verb);
		if (*status)
		{
			optimErrorMsg(pile);
			*status=-1.0;
			return std::numeric_limits<double>::quiet_NaN();
		}
	}
	tIter->second=nullfreq=start[0]=inval[0];

	// treatment 2 MAF
	std::string t2id;
	like[1]=0.0;
	tIter = treatment->begin();
	while(tIter->first == "null" || tIter->first == t1id) ++tIter;
	t2id=tIter->first;
	pile->assignTreatment(t2id);
	inval[0]=treatfreqFast(pile, &t2id, pile->minorid());

	like[1] = findmax_bfgs(npars, inval, pile, maflike, NULL, lowbound, upbound, nbounds, verb, status);
	if (*status)
	{
		*status=0;
		like[1]=multiOptim(pile, start, startpoints, inval, lowbound, upbound, nbounds, 1.0/0.0, status, verb);
		if (*status)
		{
			optimErrorMsg(pile);
			*status=-1.0;
			return std::numeric_limits<double>::quiet_NaN();
		}
	}
	tIter->second=inval[0];

	// null likelihood
	inval[0] = nullfreq;
	like[0] = maflike(inval, pile);

	// calculate LR
	lr = calclr(&like[0], &like[1]);

	// check for optimization failure
	if (lr < thresh)
	{
		// try multiple starting points //

		// treatment 1
		tIter=treatment->begin();
		pile->assignTreatment(t1id);
		multiOptim(pile, start, startpoints, inval, lowbound, upbound, nbounds, 1.0/0.0, status, verb);
		while (tIter->first == "null" || tIter->first == t2id) ++tIter;
		tIter->second=nullfreq=start[0]=inval[0];

		// check if bfgs routine failed
		if (*status)
		{
			optimErrorMsg(pile);
			*status=-1.0;
			return std::numeric_limits<double>::quiet_NaN();
		}

		// treatment 2
		tIter=treatment->begin();
		pile->assignTreatment(t2id);
		like[1]=multiOptim(pile, start, startpoints, inval, lowbound, upbound, nbounds, 1.0/0.0, status, verb);
		while(tIter->first == "null" || tIter->first == t1id) tIter++;
		tIter->second=inval[0];

		// check if bfgs routine failed
		if (*status)
		{
			optimErrorMsg(pile);
			*status=-1.0;
			return std::numeric_limits<double>::quiet_NaN();
		}

		//null
		inval[0] = nullfreq;
		like[0] = maflike(inval, pile);

		// recalculate LR
		lr = calclr(&like[0], &like[1]);
		// check if optimization still failed
		if (lr < thresh)
		{
			optimErrorMsg(pile);
			*status=-1.0;
			return std::numeric_limits<double>::quiet_NaN();
		}
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

	//verb = 100; // debug

	const int npars = 1; // number of parameters to optimize
	static double inval [npars]; // starting maf value
	static int nbounds [] = {2}; // signifies boundary conditions (see bfgs.h for description)
	static double lowbound [] = {0.0}; // lower bound on maf
	static double upbound [] = {0.5}; // upper bound on maf
	double alt, null, lr;
	bool multiopt=false;
	bool weightcount = true;

	// parameters for iterative optimization
	const double thresh=-1e-6; // cutoff for determining optimization failure
	static const double startval [] = {0.0, 0.0001, 0.25, 0.49};
	const int npoints = sizeof(startval)/sizeof(double); // number of optimization start points

	// set major allele
	pile->setMajor(pile->empiricalMajor(weightcount)); // using most common allele as major
	pile->setMinor(pile->empiricalMinorFast(weightcount)); // using second most common allele

	// find maximum likelihood estimate of maf
	inval[0] = mafguess(pile, weightcount);
	alt = findmax_bfgs(npars, inval, pile, maflike, NULL, lowbound, upbound, nbounds, verb, status);
	if (*status)
	{
		*status=0;
		alt=multiOptim(pile, startval, npoints, inval, lowbound, upbound, nbounds, 1.0/0.0, status, verb);
		multiopt=true;
		if (*status)
		{
			optimErrorMsg(pile);
			*status=-1.0;
			return std::numeric_limits<double>::quiet_NaN();
		}
	}
	*mlmaf=inval[0];

	// calculate null case
	inval[0]=0.0;
	null = maflike(inval, pile);

	// calculate LR
	lr = calclr(&null, &alt);

	// check for optimization failure
	if (lr < thresh && !multiopt)
	{
		// optimize alternative maf at multiple start points
		alt=multiOptim(pile, startval, npoints, inval, lowbound, upbound, nbounds, 1.0/0.0, status, verb);
		if (*status)
		{
			optimErrorMsg(pile);
			*status=-1.0;
			return std::numeric_limits<double>::quiet_NaN();
		}

		// recalculate LR
		lr = calclr(&null, &alt);
		// check if optimization still failed
		if (lr < thresh)
		{
			optimErrorMsg(pile);
			*status=-1.0;
			return std::numeric_limits<double>::quiet_NaN();
		}
	}

	return lr;
}

double SiteStat::calclr (const double* null, const double* alt)
{
	/*
	 * null: -log likelihood of null hypothesis
	 * alt: -log likelihood of alternative hypothesis
	 */
	return 2*(*null-*alt);
}

double SiteStat::multiOptim (Pileup* pile, const double* startval, int npoints, double* par, double* lb, double* ub, int* nbounds, double min, int* status, int verb)
{
	const int npars = 1;
	double like;
	double inval;
	int fail;
	double* inptr = &inval;

	*status=-1.0;

	for (int i=0; i<npoints; ++i)
	{
		fail=0;
		inval = startval[i];
		like=findmax_bfgs(npars, inptr, pile, maflike, NULL, lb, ub, nbounds, verb, &fail);
		if (like <= min && !fail)
		{
			min=like;
			*par=inval;
			*status=fail;
		}
	}
	return min;
}


double SiteStat::mafguess (Pileup* pile, bool wt)
{
	char minor;
	const static char a [] = {'A', 'C', 'G', 'T'};
	int i;
	double c = 0.0;
	double total = 0.0;
	double n;

	minor = pile->minorid() ? pile->minorid() : pile->empiricalMinorFast(wt);

	for (i=0; i<4; ++i)
	{
		n = wt ? pile->wtalleleCount(a[i]) : static_cast<double>(pile->alleleCount(a[i]));
		if (a[i] == minor) c = n;
		total += n;
	}

	return (c/total);
}

double SiteStat::treatfreqFast (Pileup* pile, std::string* id, char allele)
{
	if (allele == 'N')
		return 0.0;
	return pile->treatcounts(id, allele)/pile->treatcounts(id);
}

double SiteStat::treatfreq (Pileup* pile, std::string* id, char allele)
{
	if (allele == 'N')
		return 0.0;
	unsigned int c = 0;
	unsigned int total = 0;
	static std::vector<SiteData>::const_iterator tIter;
	for(tIter = pile->seqdat.begin(); tIter != pile->seqdat.end(); ++tIter)
	{
		if (tIter->id() == *id)
		{
			c += tIter->cov(allele);
			total += tIter->cov();
		}
	}
	return(static_cast<double>(c)/total);
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
	char minor=seqdata->minorid();

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
			{
				if (readIter->first == maj || readIter->first == minor) // skip reads that are not major or minor
					ptmp[k] += log(readProb(readIter->first, minor, k, readIter->second, s, MAXQUAL));
				//kahanSum(log(readProb(readIter->first, maj, k, readIter->second, s, MAXQUAL)), &ptmp[k], &comp); /* caused rounding error in optimization */
			}
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

double SiteStat::readProb (const char obs, const char minor, unsigned int nminor, unsigned int qscore, const unsigned int nchr, const unsigned int qmax)
{
	// for probability matrix rows = # possible true minor alleles in pool, columns = quality score //
	static Matrix<double> mintab(nchr+1, qmax+1, 0.0); // stores all possible values for P(X_r|k,err), X_r = minor allele //
	static Matrix<double> majtab(nchr+1, qmax+1, 0.0); // stores all possible values for P(X_r|k,err), X_r = major allele //
	static Matrix<double>* p = NULL;

	if (nminor <= nchr && qscore <= qmax)
	{
		p = obs == minor ? &mintab : &majtab;
		return ((*p)[nminor][qscore] != 0.0 ? (*p)[nminor][qscore] : ((*p)[nminor][qscore]=pread(nminor, qscore, nchr, obs, minor)));
	}
	else
		return pread(nminor, qscore, nchr, obs, minor);
}

double SiteStat::pread (int k, double qscore, int s, const char obs, const char minor)
{
	/*
	* k = number of true minor alleles
	* qscore = phred scaled base quality score
	* s = haploid pool sample size
	* obs = observed allele
	* major = major allele
	*/

	double err = pow(10, -qscore/10.0);

	if (obs == minor)
		return (static_cast<double>(k)/s)*(1.0-err) + (static_cast<double>(s-k)/s)*(err/3.0);
	else
		return (static_cast<double>(k)/s)*err + (static_cast<double>(s-k)/s)*(1.0 - err + (2*err)/3);
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
