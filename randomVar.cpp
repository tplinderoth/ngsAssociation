/*
* randomVar.cpp
*/

#include <math.h>
#include "randomVar.h"
#include "Matrix.h"

// FUNCTION DEFINITIONS

double randomVar::gammaln (const double z)
{
/*
* Numerical Recipes 2nd ed. pg 219
* returns ln(gamma(z)) for z > 0
*/

	if (z <= 0)
	{
		fprintf(stderr, "invalid value %f given to randomVar::gammaln\n", z);
		return -1;
	}

	int j;
	double x, y, tmp, ser;
	static const double cof[6]={76.18009172947146, -86.50532032941677,
		24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};

	y=x=z;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0; j<6; ++j)
		ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

double randomVar::factrl (const int n)
{
/*
* Numerical Recipes 2nd ed. pg 219
* returns n!
*/

	static int ntop=4;
	const static int overflow_thresh = 33;
	static double a[overflow_thresh]={1.0, 1.0, 2.0, 6.0, 24.0};
	int j;

	if (n < 0)
	{
		fprintf(stderr, "Negative value %i used in randomVar::factrl\n", n);
		return -1;
	}
	if (n > overflow_thresh - 1)
		return exp(gammaln(n+1.0));
	while (ntop<n)
	{
		j=ntop++;
		a[ntop]=a[j]*ntop;
	}
	return a[n];
}

double randomVar::logfactl (const int n)
{
/*
* Numerical Recipes 2nd ed. page 220
* returns ln(n!)
*/

	const static int tablesz=101;
	static double a[tablesz];

	if (n < 0)
	{
		fprintf(stderr, "Negative value %i used in randomVar::logfactrl\n", n);
		return -1;
	}
	else if (n <= 1.0)
		return 0.0;
	else if (n < tablesz)
		return (a[n] != 0.0 ? a[n] : (a[n]=gammaln(n+1.0)));
	else
		return gammaln(n+1.0);
}

double randomVar::binomcoef (const int n, const int k)
{
/*
* Numerical Recipes 2nd ed. page 220
* return the binomial coefficient choose(n,k)
*/

	if (n < 0.0 || k < 0.0)
	{
		fprintf(stderr, "Negative value used for randomVar::binomcoef: choose(%i,%i)\n",n,k);
		return -1;
	}
	else if (n < 1.0 || k > n)
		return 0.0;
	else
		return floor(0.5+exp(logfactl(n)-logfactl(k)-logfactl(n-k)));

}

double randomVar::binCoefTab (const int n, const int k)
{
/*
* returns the binomial coefficient choose(n,k)
* use when many different binomial coefficients need to be computed
* increased speed at the expense of memory usage for storing a table
*/

	static Matrix<double> coef(101, 101);

	if (n < 1.0 || k > n)
		return 0.0;

	static unsigned int nidx = static_cast<unsigned int> (n);
	static unsigned int kidx = static_cast<unsigned int> (k);

	if (nidx < coef.rown() && kidx < coef.coln())
		return(coef[nidx][kidx] != 0.0 ? coef[nidx][kidx] : (coef[nidx][kidx]=binomcoef(n,k)));
	else
		return binomcoef(n,k);
}

double randomVar::binomProb (const int n, const int k, const double p)
{
	if ( (k < 0.0 || k > n ) || p > 1.0 || p < 0.0)
	{
		fprintf(stderr, "invalid arguments to randomVar::binomProb n=%i, k=%i, p=%f\n",n,k,p);
		return -1;
	}
	else if (p == 0.0)
	{
		if (k > 0)
			return 0.0;
		else
			return 1.0;
	}
	else if (p == 1.0)
		if (k < n)
			return 0.0;
		else
			return 1.0;
	else
		return exp(logfactl(n) - logfactl(k) -logfactl(n-k) + static_cast<double>(k)*log(p) + static_cast<double>(n-k)*log(1.0-p));
}
