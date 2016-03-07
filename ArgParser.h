/*
 * ArgParser.h
 * used to handle user input and distribute it to other functions throughout the program
 */

#ifndef ARGPARSER_H_
#define ARGPARSER_H_

#include <fstream>
#include "parsePileup.h"

class ArgParser
{
public:
	/* FUNCTIONS */
	ArgParser ();
	int parseInput (const int c, char** v, const char* version); /* parse command line arguments to set data members*/
	int setStreams (std::string infile, std::string outfile); /* sets IO streams */
	static int parseTreatment (const std::string fname, std::vector< std::pair<std::string, double> >* treatment, Pileup* pile, int lrmethod); /* get treatment info from file */
	double minQ () const; /* return _minQ */
	double offset () const; /* return _Qoffset */
	unsigned int minpool () const; /* return _minind */
	unsigned int mincov () const; /* return _mincov */
	unsigned int poolsz () const; /* return _poosz */
	int lrmethod () const; /* return _lrstat */
	int printInd () const; /* returns _printInd */
	std::string inpileup_name() const; /* returns _infile */
	std::string outfile_name() const; /* returns _outfile */
	std::string treatmentf_name() const; /* returns _treatfile */
	std::istream& input (); /* returns &_is */
	std::ostream& output (); /* returns _os */
private:
	/* FUNCTIONS */
	int help (const char* arg, const char* version); /* calls subroutine help */
	void maininfo (const char* v); /* print main help */
	void summinfo (); /* print help information for 'summarize' subroutine */
	void associnfo (); /* print help information for 'association' subroutine */
	int commandCheck(const char* command); /* checks for valid commands */
	/* MEMBER VARIABLES */
	double _minQ; /* min base quality score */
	double _Qoffset; /* min possible ASCII decimal value used to encode quality scores */
	unsigned int _minpooln; /* min number of 'covered' pools */
	unsigned int _mincov; /* min per pool depth for a pool to be considered 'covered' */
	unsigned int _poolsz; /* haploid size of each pool */
	int _printIndiv; /* controls whether individual information should be printed */
	std::string _infile; /* name of input pileup file */
	std::string _outfile; /* name of output file */
	std::string _treatfile; /* name of input file containing pool treatments*/
	/* IO streams */
	std::fstream _fin;
	std::fstream _fout;
	std::istream& _is;
	std::ostream _os;
	int _lrstat; /* calculate association LR using 1 => TL method or 2 => SYK method */
	/* status flags */
	int _fail; /* flag to denote program errors */
};


#endif /* ARGPARSER_H_ */
