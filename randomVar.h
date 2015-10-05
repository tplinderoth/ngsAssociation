/*
* randomVar.h
*/

#ifndef _RANDOMVAR_H_
#define _RANDOMVAR_H_

namespace randomVar
{
	double gammaln (const double z); /* returns ln(gamma(z)) */
	double factrl (const int n); /* returns n! */
	double logfactl (const int n); /* returns ln(n!) */
	double binomcoef (const int n, const int k); /* returns the binomial coefficient choose(n,k) */
	double binCoefTab (const int n, const int k); /* returns binomial coefficent choose(n,k) from table */
	double binomProb (const int n, const int k, const double p); /* Binomial(k;n,p) */
}

#endif /* _RANDOMVAR_H_ */
