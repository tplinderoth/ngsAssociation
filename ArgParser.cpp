/*
 * ArgParser.cpp
 *
 */

#include "ArgParser.h"
#include "generalUtils.h"
#include <iomanip>
#include <cstring>

ArgParser::ArgParser ()
	: _minQ(20),
	  _Qoffset(33),
	  _minpooln(1),
	  _mincov(1),
	  _poolsz(2),
	  _printIndiv(0),
	  _is(std::cin),
	  _os(std::cout.rdbuf()),
	  _lrstat(1),
	  _fail(0)
{ }

int ArgParser::parseInput (const int c, char** v, const char* version)
{
	int argPos = 1;
	int increment = 0;

	if (c < 3)
	{
		if (c < 2)
		{
			maininfo(version);
			return 1;
		}
		else
			return (help(v[argPos],version));
	}

	const char* command = v[argPos];
	if(commandCheck(command))
		return -1;

	++argPos;
	while (argPos < c)
	{
		if (strcmp(v[argPos], "-lrstat") == 0)
		{
			_lrstat = atoi(v[argPos+1]);
			switch (_lrstat)
			{
				case 1 :
					break;
				case 2 :
					break;
				default :
					fprintf(stderr, "Invalid argument to -lrstat\n");
					return -1;
			}
		}
		else if (strcmp(v[argPos], "-infile") == 0)
		{
			_infile = (fexists(v[argPos + 1]) || strcmp(v[argPos+1], "-") == 0) ? v[argPos+1] : "";
			if (_infile.empty())
			{
				fprintf(stderr, "Couldn't open input file %s\n",v[argPos+1]);
				return -1;
			}
		}
		else if (strcmp(v[argPos], "-outfile") == 0)
		{
			_outfile = v[argPos+1];
		}
		else if (strcmp(v[argPos], "-treatments") == 0)
		{
			_treatfile = fexists(v[argPos+1]) ? v[argPos+1] : "";
			if (_treatfile.empty())
			{
				fprintf(stderr,"Couldn't open file of treatments %s\n",v[argPos+1]);
				return -1;
			}
		}
		else if (strcmp(v[argPos], "-Qoffset") == 0)
		{
			_Qoffset = atof(v[argPos+1]);
			if (_Qoffset < 0)
			{
				fprintf(stderr, "-Qoffset argument must be >= 0.0\n");
				return -1;
			}
		}
		else if (strcmp(v[argPos], "-minQ") == 0)
		{
			_minQ = atof(v[argPos+1]);
			if (_minQ < 0.0)
			{
				fprintf(stderr, "-minQ argument must be >= 0.0\n");
				return -1;
			}
		}
		else if (strcmp(v[argPos], "-minpooln") == 0)
		{
			_minpooln = atoi(v[argPos+1]);
			if (_minpooln < 0)
			{
				fprintf(stderr, "-minpooln argument must be >= 0\n");
				return -1;
			}
		}
		else if (strcmp(v[argPos], "-mincov") == 0)
		{
			_mincov = atoi(v[argPos+1]);
			if (_mincov < 0)
			{
				fprintf(stderr, "-mincov argument must be >= 0\n");
				return -1;
			}
		}
		else if (strcmp(v[argPos], "-poolsz") == 0)
		{
			_poolsz = atoi(v[argPos+1]);
			if (_poolsz < 1)
			{
				fprintf(stderr,"-poolsz argument must be >= 1\n");
				return -1;
			}
		}
		else if (strcmp(v[argPos], "-printIndiv") == 0)
		{
			_printIndiv = atoi(v[argPos+1]);
			if (_printIndiv < 0)
			{
				fprintf(stderr,"-printIndiv argument must be >= 0\n");
				return -1;
			}
		}
		else
		{
			fprintf(stderr, "Unknown option: %s\n", v[argPos]);
			return -1;
		}
		argPos += 2 + increment;
		increment = 0;
	}

	return 0;
}

int ArgParser::help (const char* arg, const char* version)
{
	if (strcmp(arg,"help") == 0 || strcmp(arg, "-help") == 0)
		maininfo(version);
	else if (strcmp(arg,"association") == 0)
		associnfo();
	else if (strcmp(arg,"summarize") == 0)
		summinfo();
	else
	{
		fprintf(stderr, "\nUnknown argument: %s\n",arg);
		return -1;
	}
	return 1;
}

int ArgParser::commandCheck(const char* command)
{
	if (strcmp(command,"association")==0)
		return 0;
	else if (strcmp(command,"summarize")==0)
		return 0;
	else
	{
		fprintf(stderr,"Invalid command: %s\n",command);
		return -1;
	}
}

void ArgParser::maininfo (const char* v)
{
	int w = 12;
	std::cerr << "\nngsAssociation\nversion " << v << "\n\nUsage: ngsAssociation [command] [arguments]\n"
	<< "\nCommands:\n"
	<< "\n" << std::setw(w) << std::left << "association" << "Test for an association between allele frequency and treatment"
	<< "\n" << std::setw(w) << std::left << "summarize" << "Sequencing and allele frequency information for sites"
	<< "\n\n";
}

void ArgParser::summinfo ()
{
	int w = 12;
	std::string(indDefault);
	indDefault = _printIndiv == true ? "SET" : "NOT SET";
	std::cerr << "\nUsage: ngsAssociation summarize [arguments]\n"
	<< "\nInput:\n"
	<< "\n" << std::setw(w) << std::left << "-infile" << std::setw(w) << "FILE|-" << "pileup format file of reads and quality scores; specifying '-' will read from STDIN [" << _infile << "]"
	<< "\n" << std::setw(w) << std::left << "-outfile" << std::setw(w) << "FILE" << "name of output file; if not provided results printed to STDOUT [" << _outfile << "]"
	<< "\n" << std::setw(w) << std::left << "-poolsz" << std::setw(w) << "INT" << "haploid sample size of each pool [" << _poolsz << "]"
	<< "\n" << std::setw(w) << std::left << "-Qoffset" << std::setw(w) << "FLOAT" << "minimum possible ASCII decimal value for base quality scores [" << _Qoffset << "]"
	<< "\n" << std::setw(w) << std::left << "-minQ" << std::setw(w) << "FLOAT" << "minimum base quality score to retain read [" << _minQ << "]"
	<< "\n" << std::setw(w) << std::left << "-minpooln" << std::setw(w) << "INT" << "minimum number of covered pools to retain site [" << _minpooln << "]"
	<< "\n" << std::setw(w) << std::left << "-mincov" << std::setw(w) << "INT" << "minimum number of reads for a pool to be considered 'covered' [" << _mincov << "]"
	<< "\n" << std::setw(w) << std::left << "-printIndiv" << std::setw(w) << "INT" << "output coverage and quality score information for each pool [" << _printIndiv << "]"
	<< "\n\nOutput by field:"
	<< "\n(1) sequence ID"
	<< "\n(2) position in sequence (1-base indexed)"
	<< "\n(3) reference allele"
	<< "\n(4) site coverage summary: total_site_depth;reference_allele_count;alternate_allele_count"
	<< "\n(5) count of each allele for site: A;C;G;T;INDEL"
	<< "\n(6) maximum likelihood MAF estimate"
	<< "\n(7) likelihood ratio that site is variable"
	<< "\n(8) pool coverage (if printIndiv=1, + per A;C;G;T coverage if printIndiv=2)"
	<< "\n(9) pool read bases (if printIndiv=3)"
	<< "\n(10) base quality scores for reads (if printIndiv=3)"
	<< "\n\nFields 8-10 are repeated for each pool depending on the value of -printIndiv."
	<< "\n\n";
}

void ArgParser::associnfo ()
{
	int w=12;
	std::cerr << "\nUsage: ngsAssociation association [arguments]\n"
	<< "\nInput:\n"
	<< "\n" << std::setw(w) << std::left << "-lrstat" << std::setw(w) << "1|2" << "method for calculating the likelihood ratio of association [" << _lrstat << "]"
	<< "\n" << std::setw(w) << std::left << "-infile" << std::setw(w) << "FILE|-" << "pileup format file of reads and quality scores; specifying '-' will read from STDIN [" << _infile << "]"
	<< "\n" << std::setw(w) << std::left << "-outfile" << std::setw(w) << "FILE" << "name of output file; if not provided results printed to STDOUT [" << _outfile << "]"
	<< "\n" << std::setw(w) << std::left << "-treatments" << std::setw(w) << "FILE" << "file of treatment IDs for pools [" << _treatfile << "]"
	<< "\n" << std::setw(w) << std::left << "-poolsz" << std::setw(w) << "INT" << "haploid sample size of each pool [" << _poolsz << "]"
	<< "\n" << std::setw(w) << std::left << "-Qoffset" << std::setw(w) << "FLOAT" << "minimum possible ASCII decimal value for base quality scores [" << _Qoffset << "]"
	<< "\n" << std::setw(w) << std::left << "-minQ" << std::setw(w) << "FLOAT" << "minimum base quality score to retain read [" << _minQ << "]"
	<< "\n" << std::setw(w) << std::left << "-minpooln" << std::setw(w) << "INT" << "minimum number of covered pools to retain site [" << _minpooln << "]"
	<< "\n" << std::setw(w) << std::left << "-mincov" << std::setw(w) << "INT" << "minimum number of reads for a pool to be considered 'covered' [" << _mincov << "]"
	<< "\n\n'treatments' file format: list pool treatment IDs in the order that they appear in the input pileup. Example:\ncontrol\ncontrol\ncase\ncase\ncontrol\ncase"
	<< "\n\nOutput by field:"
	<< "\n(1) sequence ID"
	<< "\n(2) position in sequence (1-base indexed)"
	<< "\n(3) -log null likelihood"
	<< "\n(4) -log alternative likelihood"
	<< "\n(5) likelihood ratio of allelic association"
	<< "\n(6) null MAF"
	<< "\n(7+) alternative MAF"
	<<"\n\nLR method 1: the null MAF is the frequency of the minor allele in the first treatment appearing in the -treatments file"
	<< "\nLR method 2: alternative (treatment) MAFs are listed in the order that treatments uniquely appear in the -treatments file."
	<< "\n\n";
}

int ArgParser::setStreams (const std::string infile, const std::string outfile)
{
	// open input stream
	if (!infile.empty())
	{
		if (infile != "-")
		{
			if (getFILE(_fin, infile.c_str(), "in")) // open the pileup file
			{
				_is.rdbuf(_fin.rdbuf()); // change input buffer from STDIN to fstream's buffer
				std::cerr << "Reading sequencing data from " << infile << "\n";
			}
			else
			{
				std::cerr << "Problem opening input file ...\n";
				_fail = 1;
				return (_fail);
			}
		}
		else
		std::cerr << "Reading sequencing data from standard input\n";
	}
	else
	{
		std::cerr << "Problem parsing input stream ... \n";
		_fail=1;
		return(_fail);
	}

	// open output steam
	if (!outfile.empty())
	{
		if (getFILE(_fout, outfile.c_str(), "out"))
		{
			_os.rdbuf(_fout.rdbuf()); // switch output stream's buffer to output file's buffer
			std::cerr << "Dumping results to " << outfile << "\n";
		}
		else
		{
			std::cerr << "Problem opening output file ...\n";
			_fail = 1;
			return(_fail);
		}
	}
	else
		std::cerr << "Dumping results to standard output\n";

	return 0;
}

int ArgParser::parseTreatment (const std::string fname, std::vector< std::pair<std::string, double> >* treatment, Pileup* pile, int lrmethod)
{
	std::fstream fin;
	const int n = 256;
	char s[n];
	unsigned int i = 0;
	std::vector< std::pair<std::string,double> >::const_iterator j;
	std::string id;
	int seen=0;

	if (pile->seqdat.empty())
	{
		fprintf(stderr,"Error: No pools to assign treatments to in call to ArgParser::parseTreatments\n");
		return -1;
	}

	if (!getFILE(fin, fname.c_str(), "in"))
		return -1;

	if (!treatment->empty())
		treatment->clear();
	treatment->reserve(3);
	treatment->push_back(std::make_pair("null", 0.0));

	while(fin.getline(s,n))
	{
		if (strcmp(s,"")==0)
			continue;
		if (i < pile->seqdat.size())
		{
			pile->seqdat[i].id() = s;
		}
		for (j=treatment->begin(); j!=treatment->end(); ++j)
		{
			seen=0;
			if (strcmp(s,j->first.c_str()) == 0)
			{
				seen=1;
				break;
			}
		}
		if (!seen)
		{
			treatment->push_back(std::make_pair(s, 0.0));
			id = s;
			pile->addtreatment(&id);
		}
		++i;
	}
	if (i != pile->seqdat.size())
	{
		fprintf(stderr, "Error: Number of pools in pileup and treatment ID's in %s differ\n",fname.c_str());
		return -1;
	}
	if (treatment->size() < 2)
	{
		fprintf(stderr, "Error: -treatments file has less than 2 different treatments\n");
		return -1;
	}
	if (lrmethod == 1 && treatment->size()-1 > 2)
	{
		fprintf(stderr, "Error: -lrstat 1 only works for 2 treatments\n");
		return -1;
	}

	return 0;
}

double ArgParser::minQ () const
{
	return _minQ;
}

double ArgParser::offset () const
{
	return _Qoffset;
}

unsigned int ArgParser::minpool () const
{
	return _minpooln;
}

unsigned int ArgParser::mincov () const
{
	return _mincov;
}

int ArgParser::printInd() const
{
	return _printIndiv;
}

std::string ArgParser::inpileup_name() const
{
	return _infile;
}

std::string ArgParser::outfile_name() const
{
	return _outfile;
}

std::string ArgParser::treatmentf_name() const
{
	return _treatfile;
}

std::istream& ArgParser::input ()
{
	return _is;
}

std::ostream& ArgParser::output ()
{
	return _os;
}

unsigned int ArgParser::poolsz () const
{
	return _poolsz;
}

int ArgParser::lrmethod () const
{
	return _lrstat;
}
