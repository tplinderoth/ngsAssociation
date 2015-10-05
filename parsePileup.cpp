// parsePileup.cpp

#include <cstdio>
#include <algorithm> // sort
#include <math.h> // pow
#include "parsePileup.h"
#include "generalUtils.h"


// SiteData structure constructor
SiteData::SiteData ()
	: depth(0)
{}

unsigned int SiteData::cov (char allele) const
{
	if (allele)
	{
		switch(allele)
		{
			case 'A' :
			case 'a' :
				return allecount[0];
			case 'C' :
			case 'c' :
				return allecount[1];
			case 'G' :
			case 'g' :
				return allecount[2];
			case 'T' :
			case 't' :
				return allecount[3];
			default :
				fprintf(stderr,"Unrecognized allele %c in call to SiteData::cov\n",allele);
				return 0;
		}
	}
	else
		return depth;
}

std::string SiteData::id() const
{
	if (!_id.empty())
		return _id;
	else
		return "";
}

std::string& SiteData::id()
{
	return _id;
}

// Pileup class constructor
Pileup::Pileup ()
	: _pos(0),
	  _encode(33.0),
	  _minQ(13.0),
	  _fail(0),
	  _nind(0),
	  _depthReserve(20),
	  _numalt(0),
	  _numref(0),
	  _alleles(),
	  _refallele('\0'),
	  _poolsz(2)
{
	for (int i=0; i<2; ++i)
		_majmin[i]='\0';
}

// Pileup::setQualCode assigns quality encoding
void Pileup::setQualCode (double code)
{
	if (code < 0.0)
	{
		fprintf(stderr, "Invalid quality score encoding '%f' in Pileup::getCode\n", code);
		_fail = 1;
	}
	_encode = code;
}

// Pileup::setMinQ assigns minimum quality score for reads to be considered
void Pileup::setMinQ (double q)
{
	if (q < 0.0)
	{
		fprintf(stderr, "Cannot set minimum quality score to < 0.0 in Pileup::setMinQ\n");
		_fail = 1;
	}
	_minQ = q;
}

int Pileup::getMinQ ()
{
 return _minQ;
}

int Pileup::getQualCode ()
{
	return _encode;
}

// Pileup::getSeqDat extracts information from pileup line
int Pileup::getSeqDat (const std::string& pile)
{
	int i;

	if (pile.empty())
	{
		fprintf(stderr, "No pileup line supplied to Pileup::getSeqDat\n");
		_fail = 1;
	}

	std::vector<std::string> pvec = split(pile, '\t');

	if (pvec.size() == 1)
	{
		fprintf(stderr, "Piluep line couldn't be split by Pileup::getSeqDat; check delimiter\n");
		_fail = 1;
	}

	// assign name
	if (!pvec[0].empty())
		_name = pvec[0];
	else
	{
		fprintf(stderr, "Attempt to assign empty container to name by Pileup::getSeqDat\n");
		_fail = 1;
	}

	// assign position
	if (!pvec[1].empty())
		_pos = atoi(pvec[1].c_str());
	else
	{
		fprintf(stderr, "Attempt to assign empty container to position by Pileup::getSeqDat\n");
		_fail = 1;
	}

	// get reference allele identity
	if (!pvec[2].empty())
		_refallele = toupper(pvec[2][0]);
	else
	{
		fprintf(stderr, "Attempt to assign empty container as reference allele by Pileup::getSeqDat\n");
		_fail = 1;
	}

	// get number of individuals
	if (seqdat.empty())
	{
		_nind = ExtractIndN(pile);
		if (_nind < 1)
		{
			std::cerr << "Number of individuals < 1\n";
			return -1;
		}
		initializeSeqdat(_nind);
	}

	// reset values
	_numref = 0;
	_numalt = 0;
	for (i = 0; i < 5; ++i)
		_alleles[i] = 0;
	for (i=0; i<2; ++i)
		_majmin[i] = '\0';

	// assign reads and quality scores
	size_t ind = 0;
	std::string reads;
	std::string qscores;
	std::vector<std::string>::const_iterator iter = pvec.begin() + 3;

	if ((*iter).empty())
	{
		fprintf(stderr, "Sequencing data required by Pileup::getSeqDat missing\n");
		_fail = 1;
	}

	for (iter = pvec.begin() + 3; iter != pvec.end(); ++iter)
	{
		++ind;
		seqdat[ind-1].depth = 0;
		if ( (*iter)[0] == '0' ) // no data for individual
		{
			while ((iter+1) != pvec.end() && !isdigit((*(iter+1))[0]))
				++iter;
			if ((iter+1) == pvec.end())
				break;
		}
		else
		{
			reads = *(++iter);
			while ( *((iter)->end()-1) == '^' ) // ascii "space" is mapping quality
			{
				reads += ' ';
				reads += *(++iter);
			}
			qscores = *(++iter);
			getReadDat(&reads, &qscores, ind-1);
		}
	}
	return 0;
}

// assigns read and quality scores for an individual
void Pileup::getReadDat (std::string* reads, std::string* qual, size_t ind)
{
	char r;
	unsigned int index = 0;
	unsigned int indsize = 0;
	size_t read_num = 0;
	int refidx = 0;
	std::string::iterator q = qual->begin();

	for (int i = 0; i < 4; ++i)
		seqdat[ind].allecount[i] = 0;

	switch (_refallele)
	{
                case 'A' :
                        refidx=0;
                        break;
                case 'C' :
                        refidx=1;
                        break;
                case 'G' :
                        refidx=2;
                        break;
                case 'T' :
                        refidx=3;
                        break;
		default :
			fprintf(stderr, "warning: unrecognized reference allele '%c' at %s position %d\n", _refallele, _name.c_str(), _pos);
			refidx = -1;
	}

	for(std::string::const_iterator r_iter = reads->begin(); r_iter != reads->end(); ++r_iter) // go over all reads
	{
		r = toupper(*r_iter);

		switch(r)
		{
			/* reference allele */
			case '.' :
			case ',' :
				recordRead(ind, q, read_num, index, refidx, _refallele);
				break;
			/* alternate allele */
			case 'A' :
				recordRead(ind, q, read_num, index, 0, r);
				break;
			case 'C' :
				recordRead(ind, q, read_num, index, 1, r);
				break;
			case 'G' :
				recordRead(ind, q, read_num, index, 2, r);
				break;
			case 'T' :
				recordRead(ind, q, read_num, index, 3, r);
				break;
			/* missing read or gap */
			case 'N' :
			case '*' :
			case '<' :
			case '>' :
				++q;
				++index;
				break;
			/* indel */
			case '+' :
			case '-' :
				indsize = indelSize(reads, (index+1));
                        	r_iter += indsize;
                        	index += indsize + 1;
                        	++_alleles[4];
				break;
			/* beginning of read */
			case '^' :
				++r_iter;
				index += 2;
				break;
			/* end of read */
			case '$' :
				++index;
				break;
			/* unrecognized symbol */
			default :
				fprintf(stderr, "Unrecognized symbol '%c' in %s position %u:\n%s\t%s",r,_name.c_str(),_pos,reads->c_str(),qual->c_str());
				_fail = 1;
		}

	}
	seqdat[ind].depth = read_num;
}

void Pileup::recordRead (const size_t ind, std::string::iterator& q, size_t& read_num, unsigned int& index, const int id, const char read)
{
	static const float factor = 0.5;
	bool add = true;
	double phredq = static_cast<double>(*q) - _encode;

	if (phredq < 0.0)
	{
		fprintf(stderr, "ERROR: Negative quality score '%f' at %s %u --> check quality score encoding\n",phredq,_name.c_str(),_pos);
		_fail = 1;
	}

	if (phredq >= _minQ)
	{
        	if (read_num >= seqdat[ind].rdat.size())
        	{
                	if (read_num < seqdat[ind].rdat.capacity())
                	{
                        	seqdat[ind].rdat.push_back(std::make_pair(read, phredq));
                        	add=false;
                	}
                	else
                        	seqdat[ind].rdat.resize( seqdat[ind].rdat.size() + factor * seqdat[ind].rdat.size() );
		}
		if (add)
		{
        		seqdat[ind].rdat[read_num].first = read;
        		seqdat[ind].rdat[read_num].second = phredq;
		}
		++read_num;
		if (read == _refallele)
			++_numref;
		else
			++_numalt;
		++seqdat[ind].allecount[id];
		++_alleles[id];
	}

	++q;
	++index;
}

// Pileup::indelSize gets the size of an indel + the number of digits comprising the size
unsigned int Pileup::indelSize (std::string* s, unsigned int start)
{
        unsigned int len = 0;
        std::string::const_iterator it = s->begin() + start;
        while ( isdigit(*it) && it != s->end())
        {
                ++len;
                ++it;
        }

        std::string indsize = s->substr (start,len);
        return ( atoi(indsize.c_str()) + indsize.length() );
}

// Pileup::fail returns value of _fail member
int Pileup::fail ()
{
	return _fail;
}

size_t Pileup::ExtractIndN (const std::string& line)
{
	size_t nind = 0;
	int inc = 0;
	if (!line.empty())
	{
		std::vector<std::string> ptoke = split(line, '\t');
		for (std::vector<std::string>::const_iterator iter = ptoke.begin() + 3; iter != ptoke.end(); ++iter)
		{
			++nind;
			inc = 0;
			if ( (*iter)[0] == '0' ) // no data for individual
			{
				while ((iter+1) != ptoke.end() && !isdigit((*(iter+1))[0]))
					++iter;
				if ((iter+1) == ptoke.end())
					break;
			}
            		else
            		{
            			while ( *((iter+(1+inc))->end()-1) == '^' ) // ascii "space" is mapping quality
            			{
            				++inc;
            			}
            			iter += 2 + inc;
            		}
		}
	}
	else
	{
		fprintf(stderr, "No sequencing data in string provided to Pileup::ExtractIndN\n");
		_fail = 1;
	}
	return nind;
}

void Pileup::setIndN (size_t n)
{
	if (n > 0)
		_nind = n;
	else
		fprintf(stderr, "Attempt to set nonpositive number of individuals in Pileup::setIndN\n");
}

void Pileup::initializeSeqdat (size_t n)
{
	if (n >= 0)
	{
		size_t start_size = _depthReserve;
		seqdat.resize(n);
		for (size_t i = 0; i < n; ++i)
		{
			seqdat[i].rdat.resize(start_size);
			seqdat[i].depth = start_size;
		}
		_nind=n;
	}
}

unsigned int Pileup::setn (std::string ins)
{
	unsigned int n=0;

	if (ins.empty())
	{
		fprintf(stderr,"Empty pileup line passed to Pileup::setn");
		return 0;
	}

	if (!seqdat.empty())
	{
		seqdat.clear();
	}
	else
	{
		n=ExtractIndN(ins);
		initializeSeqdat(n);
	}
	return n;
}

std::string Pileup::seqName () const
{
	return _name;
}

unsigned int Pileup::position () const
{
	return _pos;
}

void Pileup::setDepthReserve (size_t depth)
{
	_depthReserve = depth;
}

size_t Pileup::nInd () const
{
	return _nind;
}

size_t Pileup::altcount () const
{
	return _numalt;
}

size_t Pileup::refcount () const
{
	return _numref;
}

double Pileup::altfreq () const
{
	return  static_cast <double> (_numalt) / (_numref + _numalt);
}

size_t Pileup::siteDepth () const
{
	return _numref + _numalt;
}

unsigned int Pileup::alleleCount (const char allele) const
{
	unsigned int n = 0;
	switch (allele)
	{
		case 'A' :
		case 'a' :
			n = _alleles[0];
			break;
		case 'C' :
		case 'c' :
			n = _alleles[1];
			break;
		case 'G' :
		case 'g' :
			n = _alleles[2];
			break;
		case 'T' :
		case 't' :
			n = _alleles[3];
			break;
		case 'I' :
		case 'i' :
			n = _alleles[4];
			break;
		default :
			fprintf(stderr, "Unrecognized base '%c' in call to Pileup::alleleCount\n", allele);
	}
	return n;
}

char Pileup::refAllele () const
{
	return _refallele;
}

void Pileup::assignTreatment (std::string id)
{
	_treatment = id;
}

std::string Pileup::treatment () const
{
	return _treatment;
}

void Pileup::setpoolsz (unsigned int n)
{
	/* n is the haploid sample size of each pool
	 * for diploids, n=2
	 */
	_poolsz = n;
}

unsigned int Pileup::poolsz (int ploidy) const
{
	/*
	 * if ploidy argument is supplied, number of individuals comprising each pool returned
	 * otherwise, the haploid size of each pool is returned
	 */

	if (_poolsz % ploidy == 0)
		return _poolsz/ploidy;
	else
		fprintf(stderr,"Incorrect ploidy passed to Pileup::poolsz(): haploid size = %u, ploidy = %i\n",_poolsz,ploidy);
	return 0;
}

void Pileup::setMajor(char allele)
{
	allele = toupper(allele);
	switch (allele)
	{
		case 'A' :
		case 'C' :
		case 'G' :
		case 'T' :
		case 'N' :
			_majmin[0] = allele;
			break;
	default :
		fprintf(stderr,"Tried passing invalid allele '%c' to Pileup::setMajor\n",allele);
	}
}

void Pileup::setMinor(char allele)
{
	allele = toupper(allele);
	switch (allele)
	{
		case 'A' :
		case 'C' :
		case 'G' :
		case 'T' :
		case 'N' :
			_majmin[1] = allele;
			break;
	default :
		fprintf(stderr,"Tried passing invalid allele '%c' to Pileup::setMinor\n",allele);
	}
}

char Pileup::majorid () const
{
	return _majmin[0];
}

char Pileup::minorid () const
{
	return _majmin[1];
}

char Pileup::empiricalMajor ()
{
	char maj='A';
	static const char a[] = {'A', 'C', 'G', 'T'};
	for (int i=0; i<4; ++i)
	{
		if (_alleles[i] > alleleCount(maj))
			maj=a[i];
	}
	return maj;
}

bool Pileup::countCmp (std::pair<char,double> i, std::pair<char,double> j)
{
	return (i.second < j.second);
}


bool Pileup::baseCmp (std::pair<char,double> i, std::pair<char,double> j)
{
	return ((int)i.first < (int)j.first);
}

char Pileup::empiricalMinor ()
{
	static const seqread b [] = {std::make_pair('A',0.0), std::make_pair('C',0.0), std::make_pair('G',0.0), std::make_pair('T',0.0)};
	static std::vector<seqread> counts (b, b + sizeof(b)/sizeof(seqread));
	static std::vector<seqread>::iterator i;
	static std::vector<SiteData>::iterator j;

	// sort the elements of counts vector in A,C,G,T order and set counts to zero
	std::sort (counts.begin(), counts.end(), baseCmp);

	for (i=counts.begin(); i!=counts.end(); ++i)
		i->second=0.0;

	// count occurrence of each base at site weighted by the quality score
	for (j=seqdat.begin(); j!=seqdat.end(); ++j)
	{
		for (i=j->rdat.begin(); i!=j->rdat.end(); ++j)
		{
			switch (i->first)
			{
				case 'A' :
					counts[0].second += 1-error(i->second);
					break;
				case 'C' :
					counts[1].second += 1-error(i->second);
					break;
				case 'G' :
					counts[2].second += 1-error(i->second);
					break;
				case 'T' :
					counts[3].second += 1-error(i->second);
					break;
			}
		}
	}

	// sort the elements of counts vector in ascending count order and find second most common base
	std::sort (counts.begin(), counts.end(), countCmp);
	if (counts[2].second == 0.0)
		return 'N'; // site is fixed, i.e. no minor allele
	else
		return counts[2].first;
}

double Pileup::error (int q)
{
	static const int n = 71;
	static double err[n];

	if (q < 0)
	{
		fprintf(stderr,"Invalid quality score '%i' passed to Pileup:error\n",q);
		_fail = 1;
		return 1;
	}

	if (q < n)
		return (err[q] != 0.0 ? err[q] : (err[q]=scalePhred(q)));
	else
		return scalePhred(q);
}

double Pileup::scalePhred (double q)
{
	return pow(10, -q/10);
}
