#include <map>
#include <deque>
#include <string>
#include <vector>
#include <bitset>
#include <memory>
#include <cassert>
#include <sstream>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <unordered_map>

#include <tclap/CmdLine.h>
#include <tbb/parallel_sort.h>

#include <graphdump/graphdump.h>
#include <atomic>

#include <dnachar.h>
#include <streamfastaparser.h>
#include <junctionapi/junctionapi.h>

#include "bin_conv.hpp"


bool CompareJunctionsById(const TwoPaCo::JunctionPosition & a, const TwoPaCo::JunctionPosition & b)
{
	return a.GetId() < b.GetId();
}

bool CompareJunctionsByPos(const TwoPaCo::JunctionPosition & a, const TwoPaCo::JunctionPosition & b)
{
	return std::make_pair(a.GetChr(), a.GetPos()) < std::make_pair(b.GetChr(), b.GetPos());
}

struct EqClass
{
	int64_t label;
	std::vector<TwoPaCo::JunctionPosition> position;
};

int64_t Abs(int64_t x)
{
	return x > 0 ? x : -x;
}

const int64_t ID_POWER = 35;
int64_t reservedPath = int64_t(1) << (ID_POWER - 1);
const int64_t MAX_JUNCTION_ID = int64_t(1) << (ID_POWER - 4);
const int64_t MAX_SEGMENT_NUMBER = int64_t(1) << ID_POWER;
std::unordered_map<std::string, int64_t> ump;
int64_t pathId = 1;

class Segment
{
	public:
		Segment() {}
		Segment(int64_t segmentId) 
		{
			segmentId_ = segmentId;
		}
		Segment(TwoPaCo::JunctionPosition begin, TwoPaCo::JunctionPosition end, char posEdgeCh, char negEdgeCh)
		{
			bool uniquePath = false;
			int64_t absBeginId = Abs(begin.GetId());
			int64_t absEndId = Abs(end.GetId());
			if (absBeginId >= MAX_JUNCTION_ID || absEndId >= MAX_JUNCTION_ID)
			{
				throw std::runtime_error("A vertex id is too large, cannot generate GFA");
			}

			if (absBeginId < absEndId || (absBeginId == absEndId && absBeginId > 0))
			{
				uniquePath = posEdgeCh == 'N';
				segmentId_ = TwoPaCo::DnaChar::MakeUpChar(posEdgeCh);
				begin_ = begin;
				end_ = end;
			}
			else
			{
				uniquePath = negEdgeCh == 'N';
				segmentId_ = TwoPaCo::DnaChar::MakeUpChar(negEdgeCh);
				begin_ = TwoPaCo::JunctionPosition(begin.GetChr(), begin.GetPos(), -end.GetId());
				end_ = TwoPaCo::JunctionPosition(end.GetChr(), end.GetPos(), -begin.GetId());
			}

			if (!uniquePath)
			{
				if (begin_.GetId() < 0)
				{
					segmentId_ |= 1 << 2;
					segmentId_ |= Abs(begin_.GetId()) << 3;
				}
				else
				{
					segmentId_ |= begin_.GetId() << 3;
				}

				if (begin.GetId() != begin_.GetId())
				{
					segmentId_ = -segmentId_;
				}
			}
			else
			{
				segmentId_ = reservedPath++;
			}
		}

		int64_t GetSegmentId() const
		{
			return segmentId_;
		}	

		int64_t GetAbsSegmentId() const
		{
			return Abs(segmentId_);
		}

	private:
		int64_t segmentId_;
		TwoPaCo::JunctionPosition begin_;
		TwoPaCo::JunctionPosition end_;	
};

bool CompareJunctionClasses(const EqClass & a, const EqClass & b)
{
	return CompareJunctionsByPos(a.position[0], b.position[0]);
}

void GenerateGroupOutupt(const std::string & inputFileName)
{
	TwoPaCo::JunctionPosition pos;
	TwoPaCo::JunctionPositionReader reader(inputFileName.c_str());
	std::vector<EqClass> eqClass;
	std::vector<TwoPaCo::JunctionPosition> junction;
	while (reader.NextJunctionPosition(pos))
	{
		junction.push_back(pos);
	}

	std::sort(junction.begin(), junction.end(), CompareJunctionsById);
	for (size_t i = 0; i < junction.size();)
	{
		size_t j = i;
		for (; j < junction.size() && junction[i].GetId() == junction[j].GetId(); j++);
		std::sort(junction.begin() + i, junction.begin() + j, CompareJunctionsByPos);
		eqClass.push_back(EqClass());
		eqClass.back().label = junction[i].GetId();
		for (size_t k = i; k < j; k++)
		{
			eqClass.back().position.push_back(junction[k]);
		}

		i = j;
	}

	tbb::parallel_sort(eqClass.begin(), eqClass.end(), CompareJunctionClasses);
	/*for (auto junctionClass : eqClass)
	  {
	  for (auto j : junctionClass.position)
	  {
	  std::cout << j.GetChr() << ' ' << j.GetPos() << "; ";
	  }

	  std::cout << std::endl;
	  }*/

}

void GenerateOrdinaryOutput(const std::string & inputFileName)
{
	TwoPaCo::JunctionPosition pos;
	TwoPaCo::JunctionPositionReader reader(inputFileName.c_str());
	/*while (reader.NextJunctionPosition(pos))
	  {
	  std::cout << pos.GetChr() << ' ' << pos.GetPos() << ' ' << pos.GetId() << std::endl;
	  }*/
}

char Sign(int64_t arg)
{
	return arg >= 0 ? '+' : '-';
}

void ReadInputSequences(const std::vector<std::string> & genomes, std::vector<std::string> & chrSegmentId, std::vector<uint64_t> & chrSegmentLength, std::map<std::string, std::string> & fileName, bool noPrefix)
{
	size_t chrCount = 0;
	chrSegmentId.clear();	
	chrSegmentLength.clear();
	for (const std::string & chrFileName : genomes)
	{
		TwoPaCo::StreamFastaParser parser(chrFileName);
		while (parser.ReadRecord())
		{
			std::stringstream ssId;
			if (noPrefix)
			{
				ssId << parser.GetCurrentHeader();
			}
			else
			{
				ssId << "s" << chrCount << "_" << parser.GetCurrentHeader();
			}

			chrSegmentId.push_back(ssId.str());
			fileName[ssId.str()] = chrFileName;


			uint64_t size = 0;
			for (char ch; parser.GetChar(ch); ++size);
			chrSegmentLength.push_back(size);
		}
	}
}

std::ofstream output_bin;

class BinaryGenerator
{
	public:
		//static sdsl::int_vector<2> sLen;
		//static sdsl::int_vector<2> sSeq;
		static int sLenCount;
		// static sdsl::osfstream *sSeqOut;
		void Header(std::ostream & out) const
		{
			//out << "H\tVN:Z:1.0" << std::endl;
		}

		void ListInputSequences(const std::vector<std::string> & seq, std::map<std::string, std::string> & fileName, std::ostream & out) const
		{
			/*for (const auto & it : seq)
			  {
			  out << "S\t"
			  << it
			  << "\t*\tUR:Z:"
			  << fileName[it]
			  << std::endl;
			  }*/
		}

		void Segment(int64_t segmentId, uint64_t segmentSize, const std::string & body, std::ostream & out) const
		{
			sPack(Abs(segmentId), body, output_bin);
			/*out << "S\t"
			  << Abs(segmentId) << "\t"
			  << body << std::endl;
			 */
			/*
			//sdsl::int_vector<2> sSeq;
			//sdsl::int_vector<2> sLen;
			sPack(sSeq, sLen, body);
			sLenCount += sLen.size();
			// sdsl::osfstream out_tmp = *sSeqOut;
			//std::cout << sSeq.size() << std::endl;
			std::cout << sLen.size() << std::endl;
			for (int i = 0; i< sSeq.size(); i++)
			std::cout << sSeq[i];
			std::cout << std::endl;
			for (int i = 0; i< sLen.size(); i++)
			std::cout << sLen[i];
			std::cout << std::endl;
			//sdsl::append_to_file(sSeq, *sSeqOut);
			//sdsl::append_to_file(sLen, *sLenOut);
			 */
		}

		void Occurrence(int64_t segmentId, uint64_t segmentSize, const std::string & chrSegmentId, uint64_t chrSegmentSize, uint64_t begin, uint64_t end, uint64_t k, std::ostream & out) const
		{
			/*out << "C\t"
			  << Abs(segmentId) << '\t'
			  << Sign(segmentId) << '\t'
			  << chrSegmentId << "\t+\t"
			  << end << std::endl;
			 */
		}

		void Edge(int64_t prevSegmentId, uint64_t prevSegmentSize, int64_t segmentId, uint64_t segmentSize, uint64_t k, std::ostream & out) const
		{
			/*out << "L\t"
			  << Abs(prevSegmentId) << '\t'
			  << Sign(prevSegmentId) << '\t'
			  << Abs(segmentId) << '\t'
			  << Sign(segmentId) << '\t'
			  << k << 'M' << std::endl;
			 */
		}

		void FlushPath(std::vector<int64_t> & currentPath, const std::string & seqId, size_t k, std::ostream & out) const
		{

			if (currentPath.size() > 0)
			{
				/*out << "P\t" << seqId << '\t';
				  for (auto it = currentPath.begin(); it != currentPath.end() - 1; ++it)
				  {
				  out << Abs(*it) << Sign(*it) << ",";
				  }

				  out << Abs(currentPath.back()) << Sign(currentPath.back()) << "\t*" << std::endl;
				 */
				pPack(currentPath, output_bin);
				currentPath.clear();
			}
		}
};

class Gfa1Generator
{
	public:
		void Header(std::ostream & out) const
		{
			out << "H\tVN:Z:1.0" << std::endl;
		}

		void ListInputSequences(const std::vector<std::string> & seq, std::map<std::string, std::string> & fileName, std::ostream & out) const
		{
			for (const auto & it : seq)
			{
				out << "S\t"
					<< it
					<< "\t*\tUR:Z:"
					<< fileName[it]
					<< std::endl;
			}
		}

		void Segment(int64_t segmentId, uint64_t segmentSize, const std::string & body, std::ostream & out) const
		{
			out << "S\t"
				<< Abs(segmentId) << "\t" 
				<< body << std::endl;
		}

		void Occurrence(int64_t segmentId, uint64_t segmentSize, const std::string & chrSegmentId, uint64_t chrSegmentSize, uint64_t begin, uint64_t end, uint64_t k, std::ostream & out) const
		{
			out << "C\t" 
				<< Abs(segmentId) << '\t' 
				<< Sign(segmentId) << '\t'
				<< chrSegmentId << "\t+\t" 
				<< end << std::endl;
		}

		void Edge(int64_t prevSegmentId, uint64_t prevSegmentSize, int64_t segmentId, uint64_t segmentSize, uint64_t k, std::ostream & out) const
		{
			out << "L\t" 
				<< Abs(prevSegmentId) << '\t' 
				<< Sign(prevSegmentId) << '\t' 
				<< Abs(segmentId) << '\t' 
				<< Sign(segmentId) << '\t' 
				<< k << 'M' << std::endl;
		}

		void FlushPath(std::vector<int64_t> & currentPath, const std::string & seqId, size_t k, std::ostream & out) const
		{
			if (currentPath.size() > 0)
			{
				out << "P\t" << seqId << '\t';
				for (auto it = currentPath.begin(); it != currentPath.end() - 1; ++it)
				{
					out << Abs(*it) << Sign(*it) << ",";
				}

				out << Abs(currentPath.back()) << Sign(currentPath.back()) << "\t*" << std::endl;
				currentPath.clear();
			}
		}
};

std::string Gfa2Position(size_t pos, size_t length)
{
	std::stringstream ss;
	if (pos == length)
	{
		ss << pos << "$";
	}
	else
	{
		ss << pos;
	}

	return ss.str();
}

std::string Gfa2Segment(int64_t segment)
{
	std::stringstream ss;
	ss << Abs(segment) << Sign(segment);
	return ss.str();
}

class Gfa2Generator
{
	public:
		void Header(std::ostream & out) const
		{
			out << "H\tVN:Z:2.0" << std::endl;
		}

		void ListInputSequences(const std::vector<std::string> & seq, std::map<std::string, std::string> & fileName, std::ostream & out) const
		{

		}

		void Segment(int64_t segmentId, uint64_t segmentSize, const std::string & body, std::ostream & out) const
		{
			out << "S\t" 
				<< Abs(segmentId) << "\t" 
				<< segmentSize << "\t" 
				<< body << std::endl;
		}

		void Occurrence(int64_t segmentId, uint64_t segmentSize, const std::string & chrSegmentId, uint64_t chrSegmentSize, uint64_t begin, uint64_t end, uint64_t k, std::ostream & out) const
		{
			std::cout << "F\t"
				<< Abs(segmentId) << '\t'
				<< chrSegmentId << Sign(segmentId) << '\t'
				<< "0\t"
				<< segmentSize << "$" << "\t"
				<< Gfa2Position(begin, chrSegmentSize) << "\t"
				<< Gfa2Position(end + k, chrSegmentSize) << "\t"
				<< k << "M" << std::endl;
		}

		void Edge(int64_t prevSegmentId, uint64_t prevSegmentSize, int64_t segmentId, uint64_t segmentSize, uint64_t k, std::ostream & out) const
		{
			uint64_t prevSegmentStart;
			uint64_t prevSegmentEnd;
			uint64_t segmentStart;
			uint64_t segmentEnd;
			if (prevSegmentId > 0)
			{
				prevSegmentStart = prevSegmentSize - k;
				prevSegmentEnd = prevSegmentSize;
			}
			else
			{
				prevSegmentStart = 0;
				prevSegmentEnd = k;
			}

			if (segmentId > 0)
			{
				segmentStart = 0;
				segmentEnd = k;
			}
			else
			{
				segmentStart = segmentSize - k;
				segmentEnd = segmentSize;			
			}

			out << "E\t"
				<< Gfa2Segment(prevSegmentId)
				<< '\t' << Gfa2Segment(segmentId) << '\t'
				<< Gfa2Position(prevSegmentStart, prevSegmentSize) << '\t'
				<< Gfa2Position(prevSegmentEnd, prevSegmentSize) << '\t'
				<< Gfa2Position(segmentStart, segmentSize) << '\t'
				<< Gfa2Position(segmentEnd, segmentSize) << '\t'
				<< k << 'M' << std::endl;
		}

		void FlushPath(std::vector<int64_t> & currentPath, const std::string & seqId, size_t k, std::ostream & out) const
		{
			if (currentPath.size() > 0)
			{
				out << "O\t" << seqId << "p" << '\t';
				for (auto it = currentPath.begin(); it != currentPath.end() - 1; ++it)
				{
					out << Abs(*it) << Sign(*it) << " ";
				}

				out << Abs(currentPath.back()) << Sign(currentPath.back()) << std::endl;
				currentPath.clear();
			}
		}
};

	template<class G>
std::string getSegment(int64_t segmentId, std::string::iterator b, std::string::iterator e, G g, std::ostream *os, std::vector<bool> *seen)
{
	std::stringstream ss;
	if (segmentId > 0)
	{
		std::copy(b, e, std::ostream_iterator<char>(ss));
	}
	else
	{
		std::string buf = TwoPaCo::DnaChar::ReverseCompliment(std::string(b, e));
		std::copy(buf.begin(), buf.end(), std::ostream_iterator<char>(ss));
	}
	std::string seg = ss.str();
	if (ump.find(seg) == ump.end()) 
	{
		ump[seg] = segmentId > 0 ? pathId : -pathId;
		pathId++;
	}
	int64_t newSegmentId = ump[seg];
	if (!(*seen)[std::abs(newSegmentId)])
	{   
		g.Segment(newSegmentId, 0, seg, *os);
		(*seen)[std::abs(newSegmentId)] = true;
	}

	return seg;
}


bool isModifiable(int segmentId, std::vector<int64_t> *start_contigs, std::vector<int64_t> *end_contigs) {
	for (int segId : *start_contigs) {
		if (std::abs(segmentId) == std::abs(segId)) {
			return false;
		}
	}

	for (int segId : *start_contigs) {     
		if (std::abs(segmentId) == std::abs(segId)) {
			return false;
		}
	}

	return true;
}


	template<class G>
void GenerateGfaOutput(const std::vector<std::string> & genomes, size_t k, bool prefix, G g,
		std::string file_name, tbb::concurrent_queue<TwoPaCo::JunctionPosition> *queue, std::atomic<bool> * done)
{	
	std::vector<uint64_t> chrSegmentLength;
	std::vector<std::string> chrSegmentId;
	std::map<std::string, std::string> chrFileName;

	std::filebuf fb;
	fb.open(file_name, std::ios::out);
	std::ostream os(&fb);

	//std::cout << "H\tVN:Z:1.0" << std::endl;
	//g.Header(std::cout);
	g.Header(os);

	ReadInputSequences(genomes, chrSegmentId, chrSegmentLength, chrFileName, !prefix);
	g.ListInputSequences(chrSegmentId, chrFileName, os);

	std::vector<int64_t> currentPath;
	const int64_t NO_SEGMENT = 0;
	std::string chr;	
	int64_t seqId = NO_SEGMENT;
	int64_t prevSegmentId = NO_SEGMENT;
	int64_t prevSegmentSize = -1;
	TwoPaCo::JunctionPosition end;
	TwoPaCo::JunctionPosition begin;
	TwoPaCo::ChrReader chrReader(genomes);
	//TwoPaCo::JunctionPositionReader reader(inputFileName.c_str());
	std::vector<bool> seen(MAX_SEGMENT_NUMBER, 0);
	int64_t previousId = 0;

#ifdef _DEBUG
	std::map<int64_t, std::string> segmentBody;
#endif

	bool begin_set = (*queue).try_pop(begin);
	chrReader.NextChr(chr);



	std::vector<std::vector<TwoPaCo::JunctionPosition>> junc_pos_vectors;
	junc_pos_vectors.push_back({});
	std::vector<int64_t> start_contigs;
	std::vector<int64_t> end_contigs;
	bool first_contig = true;

	while (!begin_set) {
		begin_set = (*queue).try_pop(begin);
	}

	junc_pos_vectors.back().push_back(begin);


	// Hello
	while (!((*done).load(std::memory_order_relaxed)) || !(*queue).empty()) {
		if (!(*queue).empty()) {
			(*queue).try_pop(end);


			if (begin.GetChr() == end.GetChr())
			{     
				junc_pos_vectors.back().push_back(end);
				Segment nowSegment(begin, end, chr[begin.GetPos() + k], TwoPaCo::DnaChar::ReverseChar(chr[end.GetPos() - 1]));
				if (first_contig) {
					first_contig = false;
					start_contigs.push_back(nowSegment.GetSegmentId());	    
				}

				/*Segment nowSegment(begin, end, chr[begin.GetPos() + k], TwoPaCo::DnaChar::ReverseChar(chr[end.GetPos() - 1]));

				  int64_t segmentId = nowSegment.GetSegmentId();
				//currentPath.push_back(segmentId); 
				//uint64_t segmentSize = end.GetPos() + k - begin.GetPos();

				int bin = begin.GetIn(), bout = begin.GetOut(), ein = end.GetIn(), eout = end.GetOut();        
				int bid = begin.GetId(), eid = end.GetId();			
				if (bid < 0) std::swap(bin, bout);
				if (eid < 0) std::swap(ein, eout);
				std::cout << "Seg: " << begin.GetPos() + k << ' ' << chr[begin.GetPos() + k] << ' ' << end.GetPos()-1 << ' ' << TwoPaCo::DnaChar::ReverseChar(chr[end.GetPos() - 1]) << '\n';
				std::cout << "[" << bid << "] " << begin.GetPos() << '\t' << bin << ' ' << bout << '\n';
				std::cout << "[" << eid << "] " << end.GetPos() << '\t' << ein << ' ' << eout << '\n';;

				if (bin == 1 && bout == 1 && ein == 1 && eout > 1) 
				{
				//puts("1");
				seg = getSegment(segmentId, chr.begin() + begin.GetPos(), chr.begin() + end.GetPos() + k, g, &os, &seen);
				newSegmentId = ump[seg];
				currentPath.push_back(newSegmentId);
				}

				else if (bin == 1 && bout > 1 && ein > 1 && eout == 1) 
				{
				//puts("2");
				seg = getSegment(segmentId, chr.begin() + begin.GetPos() + 1, chr.begin() + end.GetPos() + k, g, &os, &seen);
				newSegmentId = ump[seg];
				currentPath.push_back(newSegmentId);

				}
				else if (bin > 1 && bout == 1 && ein == 1 && eout == 1) 
				{ 
				//puts("3");
				seg = getSegment(segmentId, chr.begin() + begin.GetPos() + 1, chr.begin() + end.GetPos() + k, g, &os, &seen);
				newSegmentId = ump[seg];
				currentPath.push_back(newSegmentId);
				}

				else if (bin == 1 && bout == 1 && ein > 1 && eout > 1) 
				{
				//puts("4");
				seg = getSegment(segmentId, chr.begin() + begin.GetPos(), chr.begin() + end.GetPos() + k - 1, g, &os, &seen);
				newSegmentId = ump[seg];
				currentPath.push_back(newSegmentId);
				// add complex node: k (simply k from start of end)
				seg2 = getSegment(42, chr.begin() + end.GetPos(), chr.begin() + end.GetPos() + k, g, &os, &seen);	// just a rnd +ve segment id
				newSegmentId2 = ump[seg2];
				currentPath.push_back(newSegmentId2);
				}

				else if (bin > 1 && bout > 1 && ein == 1 && eout == 1) 
				{
				//puts("5");
				seg = getSegment(segmentId, chr.begin() + begin.GetPos() + 1, chr.begin() + end.GetPos() + k, g, &os, &seen);
				newSegmentId = ump[seg];
				currentPath.push_back(newSegmentId);
				}
				else 
				{
				std::cout << "Shouldn't be printed\n";
				std::stringstream ss;       
				if (segmentId > 0)
				{
				std::copy(chr.begin() + begin.GetPos(), chr.begin() + end.GetPos() + k, std::ostream_iterator<char>(ss));
				}
				else
				{
				std::string buf = TwoPaCo::DnaChar::ReverseCompliment(std::string(chr.begin() + begin.GetPos(), chr.begin() + end.GetPos() + k));
				std::copy(buf.begin(), buf.end(), std::ostream_iterator<char>(ss));
				}
				//std::cout << ss.str() << " (" << bid << "): " << bin << ' ' << bout << " - (" << eid << "): " << ein << ' ' << eout << '\n';
			}*/
			/*
			   Segment nowSegment(begin, end, chr[begin.GetPos() + k], TwoPaCo::DnaChar::ReverseChar(chr[end.GetPos() - 1]));
			   int64_t segmentId = nowSegment.GetSegmentId();
			   currentPath.push_back(segmentId); 
			   uint64_t segmentSize = end.GetPos() + k - begin.GetPos();
			   if (!seen[Abs(segmentId)])
			   {
			//std::cout << "S\t" << Abs(segmentId) << "\t";
			std::stringstream ss;       
			if (segmentId > 0)
			{
			std::copy(chr.begin() + begin.GetPos(), chr.begin() + end.GetPos() + k, std::ostream_iterator<char>(ss));
			}
			else
			{
			std::string buf = TwoPaCo::DnaChar::ReverseCompliment(std::string(chr.begin() + begin.GetPos(), chr.begin() + end.GetPos() + k));
			std::copy(buf.begin(), buf.end(), std::ostream_iterator<char>(ss));
			}

			//std::cout << "Test junct " << begin.GetPos() << " : " << begin.GetIn() << ' ' << begin.GetOut() << ' ' << segmentId << std::endl;
			std::cout << begin.GetPos() << ' ' << end.GetPos() << '\n';
			std::cout << ss.str() << ' ' << begin.GetIn() << ' ' << begin.GetOut() << " - " << end.GetIn() << ' ' << end.GetOut() << '\n';
			g.Segment(segmentId, segmentSize, ss.str(), os);
			seen[Abs(segmentId)] = true;
			}

			g.Occurrence(segmentId, segmentSize, chrSegmentId[seqId], chrSegmentLength[seqId],
			begin.GetPos(), end.GetPos(), k, os);

			if (prevSegmentId != NO_SEGMENT)
			{
			g.Edge(prevSegmentId, prevSegmentSize, segmentId, segmentSize, k, os);
			}                             

			prevSegmentId = segmentId;
			prevSegmentSize = segmentSize;*/
			begin = end;
			previousId = nowSegment.GetSegmentId();
			}
			else
			{
				first_contig = true;
				//std::cout << "End position for transcript : " << begin.GetPos() << " - " << begin.GetChr() << "\n";
				end_contigs.push_back(previousId);
				junc_pos_vectors.push_back({});
				//g.FlushPath(currentPath, chrSegmentId[seqId], k, os);
				chrReader.NextChr(chr);
				//prevSegmentId = 0;
				begin = end;
				junc_pos_vectors.back().push_back(begin);

				if (begin.GetChr() != ++seqId)
				{
					throw std::runtime_error("The input is corrupted");
				}
			}

		}
	}

	end_contigs.push_back(previousId);
	//end_contigs.push_back(junc_pos_vectors.back().back());

	// Hello


	//g.FlushPath(currentPath, chrSegmentId[seqId], k, os);
	//fb.close();

	for (int i = 0; i < start_contigs.size(); i++) {
		std::cout << "start " << start_contigs[i] << " end " << end_contigs[i] << '\n';
	}

	TwoPaCo::JunctionPosition begin_n, end_n;
	bool first = true;
	seqId = NO_SEGMENT;
	std::map<std::string, std::vector<int64_t>> new_contig_map;

	TwoPaCo::ChrReader chrReader1(genomes);
	chrReader1.NextChr(chr);

        for (std::vector <TwoPaCo::JunctionPosition> path: junc_pos_vectors) {
            for (TwoPaCo::JunctionPosition jpos : path) {
                if (first) {
                    begin_n = jpos;
                    end_n = jpos;
                    first = false;
                    continue;
                }
                begin_n = end_n;
                end_n = jpos;
	        std::cout << begin_n.GetPos() << " - " << end_n.GetPos() << " ";

                Segment nowSegment(begin_n, end_n, chr[begin_n.GetPos() + k], TwoPaCo::DnaChar::ReverseChar(chr[end_n.GetPos() - 1]));
		int64_t segmentId = nowSegment.GetSegmentId();
		currentPath.push_back(segmentId);
		uint64_t segmentSize = end_n.GetPos() + k - begin_n.GetPos();

		if (!seen[Abs(segmentId)])
		    {
			//std::cout << "S\t" << Abs(segmentId) << "\t";
			std::stringstream ss;					
			if (segmentId > 0)
			    {
			        std::copy(chr.begin() + begin_n.GetPos(), chr.begin() + end_n.GetPos() + k, std::ostream_iterator<char>(ss));
                            }
			else
			    {
			        std::string buf = TwoPaCo::DnaChar::ReverseCompliment(std::string(chr.begin() + begin_n.GetPos(), chr.begin() + end_n.GetPos() + k));
				std::copy(buf.begin(), buf.end(), std::ostream_iterator<char>(ss));
			    }

			g.Segment(segmentId, segmentSize, ss.str(), os);
			seen[Abs(segmentId)] = true;
		   }



            }
            first = true;
            g.FlushPath(currentPath, chrSegmentId[seqId], k, os);
            chrReader1.NextChr(chr);
            ++seqId;
	    std::cout << "\n";
        }


	/* for (std::vector <TwoPaCo::JunctionPosition> path: junc_pos_vectors) {
		for (TwoPaCo::JunctionPosition jpos: path) {
			if (first) {
				begin_n = jpos;
				end_n = jpos;
				first = false;
				continue;
			}
			begin_n = end_n;
			end_n = jpos;

			Segment nowSegment(begin_n, end_n, chr[begin_n.GetPos() + k], TwoPaCo::DnaChar::ReverseChar(chr[end_n.GetPos() - 1]));
			int64_t segmentId = nowSegment.GetSegmentId();

			currentPath.push_back(segmentId); 
			uint64_t segmentSize = end_n.GetPos() + k - begin_n.GetPos();
			if (!seen[Abs(segmentId)])
			{
				//std::cout << "S\t" << Abs(segmentId) << "\t";
				std::stringstream ss;       
				if (segmentId > 0)
				{
					std::copy(chr.begin() + begin_n.GetPos(), chr.begin() + end_n.GetPos() + k, std::ostream_iterator<char>(ss));
				}
				else
				{
					std::string buf = TwoPaCo::DnaChar::ReverseCompliment(std::string(chr.begin() + begin_n.GetPos(), chr.begin() + end_n.GetPos() + k));
					std::copy(buf.begin(), buf.end(), std::ostream_iterator<char>(ss));
				}

				g.Segment(segmentId, segmentSize, ss.str(), os);
				seen[Abs(segmentId)] = true;
			}

			// check if it belongs in start, end (efficiently)
			// true: don't change
			// false:
			// begin = jpos, end = jpos+1
			// cin = bin, cout = eout
			// case: needs to be updated
			// check if already updated in global map (another path)
			// check if present in global set: do nothing
			// otherwise insert in set: 
			// set (key: string, val: segid) - to handle dups
			// currentPath update: (simple: cur, complex: cur, new)
		}
		first = true;
		g.FlushPath(currentPath, chrSegmentId[seqId], k, os);
		chrReader.NextChr(chr);
                //begin_n = end_n;
                //seqId++;
		if (begin.GetChr() != ++seqId) {
		    throw std::runtime_error("The input is corrupted");
	        } 
		// g.flush(currentPath)
		std::cout << '\n';
	}*/
	fb.close();
}

	template<class It>
void OutFastaBody(It begin, It end)
{
	int64_t count = 0;
	for (; begin != end; ++begin)
	{
		std::cout << *begin;
		if (++count % 80 == 0)
		{
			std::cout << std::endl;
		}
	}

	if (count % 80 != 0)
	{
		std::cout << std::endl;
	}
}

void GenerateFastaOutput(const std::string & inputFileName, const std::vector<std::string> & genomes, size_t k)
{
	std::vector<uint64_t> chrSegmentLength;
	std::vector<std::string> chrSegmentId;
	std::map<std::string, std::string> chrFileName;


	ReadInputSequences(genomes, chrSegmentId, chrSegmentLength, chrFileName, false);	

	std::vector<int64_t> currentPath;
	const int64_t NO_SEGMENT = 0;
	std::string chr;
	int64_t seqId = NO_SEGMENT;
	int64_t prevSegmentId = NO_SEGMENT;
	int64_t prevSegmentSize = -1;
	TwoPaCo::JunctionPosition end;
	TwoPaCo::JunctionPosition begin;
	TwoPaCo::ChrReader chrReader(genomes);
	TwoPaCo::JunctionPositionReader reader(inputFileName.c_str());
	std::vector<bool> seen(MAX_SEGMENT_NUMBER, 0);
	int64_t previousId = 0;

#ifdef _DEBUG
	std::map<int64_t, std::string> segmentBody;
#endif
	if (reader.NextJunctionPosition(begin))
	{
		chrReader.NextChr(chr);
		while (reader.NextJunctionPosition(end))
		{
			if (begin.GetChr() == end.GetChr())
			{
				Segment nowSegment(begin, end, chr[begin.GetPos() + k], TwoPaCo::DnaChar::ReverseChar(chr[end.GetPos() - 1]));
				int64_t segmentId = nowSegment.GetSegmentId();
				currentPath.push_back(segmentId);
				uint64_t segmentSize = end.GetPos() + k - begin.GetPos();
				if (!seen[Abs(segmentId)])
				{
					std::cout << ">" << Abs(segmentId) << std::endl;
					if (segmentId > 0)
					{
						OutFastaBody(chr.begin() + begin.GetPos(), chr.begin() + end.GetPos() + k);
					}
					else
					{
						std::string buf = TwoPaCo::DnaChar::ReverseCompliment(std::string(chr.begin() + begin.GetPos(), chr.begin() + end.GetPos() + k));
						OutFastaBody(buf.begin(), buf.end());
					}


					seen[Abs(segmentId)] = true;
				}

#ifdef _DEBUG
				int64_t absSegmentId = Abs(segmentId);
				std::string buf = segmentId > 0 ? std::string(chr.begin() + begin.GetPos(), chr.begin() + end.GetPos() + k) :
					TwoPaCo::DnaChar::ReverseCompliment(std::string(chr.begin() + begin.GetPos(), chr.begin() + end.GetPos() + k));
				if (segmentBody.count(absSegmentId) == 0)
				{
					segmentBody[absSegmentId] = buf;
				}
				else
				{
					assert(segmentBody[absSegmentId] == buf);
				}
#endif
				prevSegmentId = segmentId;
				prevSegmentSize = segmentSize;
				begin = end;
			}
			else
			{				
				chrReader.NextChr(chr);
				prevSegmentId = 0;
				begin = end;

				if (begin.GetChr() != ++seqId)
				{
					throw std::runtime_error("The input is corrupted");
				}
			}
		}
	}
}


void GenerateDotOutput(const std::string & inputFileName)
{
	TwoPaCo::JunctionPosition pos;
	TwoPaCo::JunctionPosition prevPos;
	TwoPaCo::JunctionPositionReader reader(inputFileName.c_str());
	std::cout << "digraph G\n{\n\trankdir = LR" << std::endl;

	while (reader.NextJunctionPosition(pos))
	{
		if (pos.GetChr() == prevPos.GetChr())
		{
			std::cout << '\t' << prevPos.GetId() << " -> " << pos.GetId() <<
				"[color=\"blue\", label=\"chr=" << prevPos.GetChr() << " pos=" << prevPos.GetPos() << "\"]" << std::endl;	
			std::cout << '\t' << -pos.GetId() << " -> " << -prevPos.GetId() <<
				"[color=\"red\", label=\"chr=" << prevPos.GetChr() << " pos=" << prevPos.GetPos() << "\"]" << std::endl;
		}

		prevPos = pos;
	}

	std::cout << "}" << std::endl;
}

int my_main(int argc, char * argv[], std::string file_name, tbb::concurrent_queue<TwoPaCo::JunctionPosition> *queue,
		std::atomic<bool> * done)
{
	std::vector<std::string> format;
	format.push_back("seq");
	format.push_back("group");
	format.push_back("dot");
	format.push_back("gfa1");
	format.push_back("gfa2");
	format.push_back("fasta");
	format.push_back("binary");
	std::stringstream formatString;
	std::copy(format.begin(), format.begin(), std::ostream_iterator<std::string>(formatString, "|"));
	try
	{
		TCLAP::CmdLine cmd("This utility converts the binary output of TwoPaCo to another format", ' ', "0.9.3");
		TCLAP::SwitchArg prefix("", "prefix", "Add a prefix to segments in GFA (in case if you have genomes with identical FASTA headers)", cmd, false);

		TCLAP::UnlabeledValueArg<std::string> inputFileName("infile",
				"input file name",
				false,
				"",
				"file name",
				cmd);

		TCLAP::ValuesConstraint<std::string> formatConstraint(format);
		TCLAP::ValueArg<std::string> outputFileFormat("f",
				"format",
				"Output format",
				true,
				format[0],
				&formatConstraint,
				cmd);

		TCLAP::MultiArg<std::string> seqFileName("s",
				"seqfile",
				"sequences file name",
				false,
				"",
				cmd);

		TCLAP::ValueArg<unsigned int> kvalue("k",
				"kvalue",
				"Value of k",
				true,
				25,
				"integer",
				cmd);

		cmd.parse(argc, argv);
		if (outputFileFormat.getValue() == format[0])
		{
			GenerateOrdinaryOutput(inputFileName.getValue());
		}
		else if (outputFileFormat.getValue() == format[1])
		{
			GenerateGroupOutupt(inputFileName.getValue());
		}
		else if (outputFileFormat.getValue() == format[2])
		{
			GenerateDotOutput(inputFileName.getValue());
		}
		else if (outputFileFormat.getValue() == format[3])
		{
			if (!seqFileName.isSet())
			{
				throw TCLAP::ArgParseException("Required argument missing\n", "seqfilename");
			}

			//GenerateGfaOutput(inputFileName.getValue(), seqFileName.getValue(), kvalue.getValue(), prefix.getValue(), Gfa1Generator());
			GenerateGfaOutput(seqFileName.getValue(), kvalue.getValue(), prefix.getValue(), Gfa1Generator(), file_name, queue, done);

		}
		else if (outputFileFormat.getValue() == format[4])
		{
			if (!seqFileName.isSet())
			{
				throw TCLAP::ArgParseException("Required argument missing\n", "seqfilename");
			}

			//GenerateGfaOutput(inputFileName.getValue(), seqFileName.getValue(), kvalue.getValue(), prefix.getValue(), Gfa2Generator());
			GenerateGfaOutput(seqFileName.getValue(), kvalue.getValue(), prefix.getValue(), Gfa2Generator(), file_name, queue, done);
		}
		else if (outputFileFormat.getValue() == format[5])
		{
			if (!seqFileName.isSet())
			{
				throw TCLAP::ArgParseException("Required argument missing\n", "seqfilename");
			}

			GenerateFastaOutput(inputFileName.getValue(), seqFileName.getValue(), kvalue.getValue());
		}
		else if (outputFileFormat.getValue() == format[6])
		{
			if (!seqFileName.isSet())
			{
				throw TCLAP::ArgParseException("Required argument missing\n", "seqfilename");
			}

			std::string file_bin = file_name;

			output_bin.open(file_bin, std::ios::out);
			//GenerateGfaOutput(inputFileName.getValue(), seqFileName.getValue(), kvalue.getValue(), prefix.getValue(), BinaryGenerator());
			GenerateGfaOutput(seqFileName.getValue(), kvalue.getValue(), prefix.getValue(), BinaryGenerator(), file_name, queue, done);
			output_bin.close();	

			/*
			   std::ifstream input;
			   input.open(file_bin, std::ios::in);
			//unPack(input);
			std::pair <bool, std::string> conv;
			while(isNext(input)) {
			conv = getNext(input);
			}
			input.close();
			 */	
		}
	}
	catch (TCLAP::ArgException &e)
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		return 1;
	}
	catch (std::runtime_error & e)
	{
		std::cerr << "error: " << e.what() << std::endl;
		return 1;
	}


	return 0;
}

void run_graph_dump(int argc, char * argv[], std::string file_name, tbb::concurrent_queue<TwoPaCo::JunctionPosition> *queue,
		std::atomic<bool> * done) {
	my_main(argc, argv, file_name, queue, done);
}
