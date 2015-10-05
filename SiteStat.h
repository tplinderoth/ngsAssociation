/*
* SiteStat.h
*/

#ifndef _SITESTAT_H_
#define _SITESTAT_H_

#include "parsePileup.h"

#define MAXQUAL 41 // maximum quality score value

typedef std::pair<std::string,double> treatdat;

// FUNCTION DECLARATION

namespace SiteStat
{
	double assoclrt (Pileup* pile, std::vector<treatdat>* treatment, double like [], int* status, int verb = -1); /* returns likelihood ratio for a treatment/minor-allele association*/
	double snplr (Pileup* pile, double* mlmaf, int* status, int verb = -1); /* returns a likelihood ratio of a site being variable */
	double mafguess (Pileup* pile); /* returns the empirical, non-major, allele frequency */
	void optimErrorMsg (const Pileup* seqdata); /* optimization failure message */
	double maflike (const double par [], const void* data); /* returns -log P(D_i|f), f=maf*/
	char findMajor (const Pileup* sitedata); /* finds most common read type for a site */
	double pread (int k, double qscore, int s, const char obs, const char major); /* calculates P(X_r|X_r,true)*P(X_r,true|k,err) */
	double readProb (const char obs, const char major, unsigned int nminor,
		unsigned int qscore, const unsigned int nchr = 2, const unsigned int qmax = MAXQUAL); /* returns P(X_r|X_r,true)*P(X_r,true|k,err) from table */
	void kahanSum(double summand, double* total, double* comp); /* performs Kahan summation */
	double multiOptim (Pileup* pile, const double* startval, int npoints, double* par, double* lb, double* ub, int* nbounds, double min, int* status, int verb); /* performs optimization at multiple start points */
};

#endif /* _SITESTAT_H_ */

