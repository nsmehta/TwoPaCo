#ifndef _CANDIDATE_OCCURENCE_
#define _CANDIDATE_OCCURENCE_

#include "compressedstring.h"

namespace TwoPaCo
{
	template<size_t CAPACITY>
	class CandidateOccurence
	{
	public:
		static const size_t IS_PREV_N = 1;
		static const size_t IS_NEXT_N = 2;
		static const char TRUE_BIF = 'A';
		static const char FAKE_BIF = 'C';
		static const size_t ADDITIONAL_CHAR = 4 + 10;	// 4: IN_A, 5: IN_C, 6: IN_G, 7: IN_T, 8: IN_N, 9: OUT_A, 10: OUT_C, 11: OUT_G, 12: OUT_T, 13: OUT_N
		static const size_t MAX_SIZE = CAPACITY * 32;
		static const size_t NEXT_POS = MAX_SIZE - ADDITIONAL_CHAR;
		static const size_t PREV_POS = NEXT_POS + 1;
		static const size_t NMASK_POS = NEXT_POS + 2;
		static const size_t IS_BIF_POS = NEXT_POS + 3;
		static const size_t EDGE_POS = NEXT_POS + 4;
		static const size_t VERTEX_SIZE = MAX_SIZE - ADDITIONAL_CHAR;		

		CandidateOccurence(){}

		void SetInEdges(const char base, uint64_t posHash0, uint64_t negHash0, std::string::const_iterator pos, size_t vertexLength)
		{
			if (DnaChar::IsDefinite(base))
				body_.SetChar(EDGE_POS + DnaChar::MakeUpChar(base), TRUE_BIF);
			else
				body_.SetChar(EDGE_POS + 4, TRUE_BIF);
		}
		
		void SetOutEdges(const char base, uint64_t posHash0, uint64_t negHash0, std::string::const_iterator pos, size_t vertexLength)
		{
			if (DnaChar::IsDefinite(base))
				body_.SetChar(EDGE_POS + 5 + DnaChar::MakeUpChar(base), TRUE_BIF);
			else
				body_.SetChar(EDGE_POS + 5 + 4, TRUE_BIF);
		}

		int GetInEdges()
		{
			int ret = 0;
			for (int i = 0; i < 5; i++) {
				if (body_.GetChar(EDGE_POS + i) == TRUE_BIF) 
					ret++;
			}
			return ret;
		}
		
		int GetOutEdges()
		{
			int ret = 0;
			for (int i = 0; i < 5; i++) {
				if (body_.GetChar(EDGE_POS + 5 + i) == TRUE_BIF) 
					ret++;
			}
			return ret;
		}
		
		void Set(uint64_t posHash0,
			uint64_t negHash0,			
			std::string::const_iterator pos,
			size_t vertexLength,
			char posExtend,
			char posPrev,
			bool isBifurcation)
		{	
			if (posHash0 < negHash0 || (posHash0 == negHash0 && DnaChar::LessSelfReverseComplement(pos, vertexLength)))
			{
				body_.CopyFromString(pos, vertexLength);
				body_.SetChar(NEXT_POS, posExtend);
				body_.SetChar(PREV_POS, posPrev);
				body_.SetChar(NMASK_POS, EncodeNmask(posExtend, posPrev));
			}
			else
			{
				body_.CopyFromReverseString(pos, vertexLength);
				body_.SetChar(NEXT_POS, DnaChar::ReverseChar(posPrev));
				body_.SetChar(PREV_POS, DnaChar::ReverseChar(posExtend));
				body_.SetChar(NMASK_POS, EncodeNmask(posPrev, posExtend));
			}
			for (int i = 0; i < 10; i++) {
				body_.SetChar(NEXT_POS + 4 + i, FAKE_BIF);
			}
		}

		char Prev() const
		{
			char nmask = body_.RawChar(NMASK_POS);
			return nmask & IS_PREV_N ? 'N' : body_.GetChar(PREV_POS);
		}

		char Next() const
		{
			char nmask = body_.RawChar(NMASK_POS);
			return nmask & IS_NEXT_N ? 'N' : body_.GetChar(NEXT_POS);
		}

		bool IsBifurcation() const
		{
			return body_.GetChar(IS_BIF_POS) == TRUE_BIF;
		}

		void MakeBifurcation() const
		{
			body_.SetCharConcurrently(IS_BIF_POS, TRUE_BIF);
		}

		bool EqualBase(const CandidateOccurence & occurence) const
		{
			return CompressedString<CAPACITY>::EqualPrefix(VERTEX_SIZE, occurence.body_, body_);
		}

		uint64_t Hash() const
		{
			return body_.HashPrefix(VERTEX_SIZE);
		}

		CompressedString<CAPACITY> GetBase() const
		{
			CompressedString<CAPACITY> ret;
			ret.CopyPrefixFrom(body_, VERTEX_SIZE);
			return ret;
		}

		bool operator < (const CandidateOccurence & other) const
		{
			return CompressedString<CAPACITY>::LessPrefix(body_, other.body_, VERTEX_SIZE);
		}
		
	private:
		char EncodeNmask(char next, char prev)
		{
			return DnaChar::LITERAL[(next == 'N' ? IS_NEXT_N : 0) | (prev == 'N' ? IS_PREV_N : 0)];
		}

		CompressedString<CAPACITY> body_;
	};


	inline size_t CalculateNeededCapacity(size_t vertexLength)
	{
		size_t bufSize = vertexLength + CandidateOccurence<1>::ADDITIONAL_CHAR;
		return bufSize / UNIT_CAPACITY + (bufSize % UNIT_CAPACITY == 0 ? 0 : 1);
	}
}

#endif
