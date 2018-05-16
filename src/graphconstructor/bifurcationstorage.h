#ifndef _BIFURCATION_STORAGE_H_
#define _BIFURCATION_STORAGE_H_

#include "common.h"
#include "compressedstring.h"


namespace TwoPaCo
{
	template<size_t CAPACITY>
	class BifurcationStorage
	{
	public:
		typedef CompressedString<CAPACITY> DnaString;
		BifurcationStorage(){}

		uint64_t GetDistinctVerticesCount() const
		{
			return bifurcationKey_.size();
		}

		uint64_t GetTotalVerticesCount() const
		{
			return bifurcationKey_.size() * 2;
		}

		static bool LessBif(std::pair<DnaString, std::pair<int, int> > b1, std::pair<DnaString, std::pair<int, int> > b2)
		{
			return DnaString::Less(b1.first, b2.first);
		}
		
		void Init(std::istream & bifurcationTempRead, std::istream & junctionTempRead, uint64_t verticesCount, uint64_t vertexLength, size_t threads)
		{
			uint64_t bitsPower = 0;
			vertexLength_ = vertexLength;
			while (verticesCount * 8 >= (uint64_t(1) << bitsPower))
			{
				++bitsPower;
			}

			size_t hashFunctionNumber = 3;
			bitsPower = max(bitsPower, size_t(24));
			bifurcationFilter_.assign(uint64_t(1) << bitsPower, false);
			hashFunction_.resize(hashFunctionNumber);
			for (HashFunctionPtr & ptr : hashFunction_)
			{
				ptr.reset(new HashFunction(vertexLength, bitsPower));
			}

			DnaString buf;
			int inDeg, outDeg;
			char inDeg_c, outDeg_c;
			std::string stringBuf(vertexLength, ' ');
			for (size_t i = 0; i < verticesCount; i++)
			{
				buf.ReadFromFile(bifurcationTempRead);
				if (!bifurcationTempRead)
				{
					throw StreamFastaParser::Exception("Can't read from a temporary file");
				}

				buf.ToString(stringBuf, vertexLength);
				junctionTempRead >> inDeg_c;
				junctionTempRead >> outDeg_c;
				inDeg = inDeg_c - '0';
				outDeg = outDeg_c - '0';
				bifurcationKey_.push_back(std::make_pair(buf, std::make_pair(inDeg, outDeg)));
				std::cout << stringBuf << ' ' << inDeg << ' ' << outDeg << '\n';
				for (HashFunctionPtr & ptr : hashFunction_)
				{
					uint64_t hf = ptr->hash(stringBuf);
					bifurcationFilter_[hf] = true;
				}
			}

			tbb::task_scheduler_init init(threads);
			tbb::parallel_sort(bifurcationKey_.begin(), bifurcationKey_.end(), LessBif);
		}

		std::pair<int64_t, std::pair<int, int> > GetId2(std::string::const_iterator pos) const
		{
			return GetId2(pos, true, true);
		}

		int64_t GetId(std::string::const_iterator pos) const
		{
			return GetId(pos, true, true);
		}
		/*
		int64_t GetId(std::string::const_iterator pos, const std::vector<HashFunctionPtr> & posVertexHash, const std::vector<HashFunctionPtr> & negVertexHash) const
		{
			bool posFound = true;
			bool negFound = true;
			int64_t ret = INVALID_VERTEX;
			for (size_t i = 0; i < posVertexHash.size() && (posFound || negFound); i++)
			{
				if (!bifurcationFilter_[posVertexHash[i]->hashvalue])
				{
					posFound = false;
				}

				if (!bifurcationFilter_[negVertexHash[i]->hashvalue])
				{
					negFound = false;
				}
			}

			return GetId(pos, posFound, negFound);
		}
		*/

		const std::vector<HashFunctionPtr>& GetHashFunctions() const
		{
			return hashFunction_;
		}

        

	private:
		int64_t GetId(std::string::const_iterator pos, bool posFound, bool negFound) const
		{
			DnaString bitBuf;
			int64_t ret = INVALID_VERTEX;
			if (posFound)
			{
				posFound = false;
				bitBuf.Clear();
				bitBuf.CopyFromString(pos, vertexLength_);
				auto it = std::lower_bound(bifurcationKey_.begin(), bifurcationKey_.end(), std::make_pair(bitBuf, std::make_pair(-1, -1)), LessBif);
				if (it != bifurcationKey_.end() && (*it).first == bitBuf)
				{
					posFound = true;
					ret = it - bifurcationKey_.begin() + 1;
				}
			}

			if (negFound && !posFound)
			{
				negFound = false;
				bitBuf.Clear();
				bitBuf.CopyFromReverseString(pos, vertexLength_);
				auto it = std::lower_bound(bifurcationKey_.begin(), bifurcationKey_.end(), std::make_pair(bitBuf, std::make_pair(-1, -1)), LessBif);
				if (it != bifurcationKey_.end() && (*it).first == bitBuf)
				{
					negFound = true;
					ret = -(it - bifurcationKey_.begin() + 1);
				}
			}
#ifdef _DEBUG
			bool found = false;
			for (size_t strand = 0; strand < 2; ++strand)
			{
				bitBuf.Clear();
				if (strand == 0)
				{
					bitBuf.CopyFromString(pos, vertexLength_);
				}
				else
				{
					bitBuf.CopyFromReverseString(pos, vertexLength_);
				}

				auto it = std::lower_bound(bifurcationKey_.begin(), bifurcationKey_.end(), std::make_pair(bitBuf, std::make_pair(-1, -1)), LessBif);
				if (it != bifurcationKey_.end() && (*it).first == bitBuf)
				{
					found = true;
				}
			}

			assert(found == (posFound || negFound));
#endif
			return ret;
		}

		std::pair<int64_t, std::pair<int, int> > GetId2(std::string::const_iterator pos, bool posFound, bool negFound) const
		{
			DnaString bitBuf;
			std::pair<int64_t, std::pair<int, int> > ret = std::make_pair(INVALID_VERTEX, std::make_pair(-1, -1));
			if (posFound)
			{
				posFound = false;
				bitBuf.Clear();
				bitBuf.CopyFromString(pos, vertexLength_);
				auto it = std::lower_bound(bifurcationKey_.begin(), bifurcationKey_.end(), std::make_pair(bitBuf, std::make_pair(-1, -1)), LessBif);
				if (it != bifurcationKey_.end() && (*it).first == bitBuf)
				{
					posFound = true;
					ret.first = it - bifurcationKey_.begin() + 1;
					ret.second.first = (*it).second.first;
					ret.second.second = (*it).second.second;
				}
			}

			if (negFound && !posFound)
			{
				negFound = false;
				bitBuf.Clear();
				bitBuf.CopyFromReverseString(pos, vertexLength_);
				auto it = std::lower_bound(bifurcationKey_.begin(), bifurcationKey_.end(), std::make_pair(bitBuf, std::make_pair(-1, -1)), LessBif);
				if (it != bifurcationKey_.end() && (*it).first == bitBuf)
				{
					negFound = true;
					ret.first = -(it - bifurcationKey_.begin() + 1);
					ret.second.first = (*it).second.first;
					ret.second.second = (*it).second.second;
				}
			}
#ifdef _DEBUG
			bool found = false;
			for (size_t strand = 0; strand < 2; ++strand)
			{
				bitBuf.Clear();
				if (strand == 0)
				{
					bitBuf.CopyFromString(pos, vertexLength_);
				}
				else
				{
					bitBuf.CopyFromReverseString(pos, vertexLength_);
				}

				auto it = std::lower_bound(bifurcationKey_.begin(), bifurcationKey_.end(), std::make_pair(bitBuf, std::make_pair(-1, -1)), LessBif);
				if (it != bifurcationKey_.end() && (*it).first == bitBuf)
				{
					found = true;
				}
			}

			assert(found == (posFound || negFound));
#endif
			return ret;
		}


		DISALLOW_COPY_AND_ASSIGN(BifurcationStorage<CAPACITY>);

		size_t vertexLength_;
		std::vector<bool> bifurcationFilter_;
		std::vector<std::pair<DnaString, std::pair<int, int> > > bifurcationKey_;
		std::vector<HashFunctionPtr> hashFunction_;
	};
}

#endif
