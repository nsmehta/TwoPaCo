#include <deque>
#include <ctime>
#include <memory>
#include <numeric>
#include <cassert>
#include <fstream>
#include <iostream>
#include <algorithm>

#include <boost/ref.hpp>
#include <boost/lockfree/spsc_queue.hpp>

#include <tbb/blocked_range.h>
#include <tbb/parallel_sort.h>
#include <tbb/parallel_reduce.h>

#include "ngramhashing/cyclichash.h"

#include "vertexenumerator.h"

namespace Sibelia
{
	const size_t VertexEnumerator::INVALID_VERTEX = -1;

	namespace
	{
		typedef CyclicHash<uint64_t> HashFunction;
		typedef std::unique_ptr<HashFunction> HashFunctionPtr;

		template<class F>
		bool PutInBloomFilter(ConcurrentBitVector & filter, std::vector<HashFunctionPtr> & hf, F f, char farg, uint64_t hash0)
		{
			for (size_t i = 0; i < hf.size(); i++)
			{
				uint64_t hvalue = i == 0 ? hash0 : ((*hf[i]).*f)(farg);
				filter.SetConcurrently(hvalue);
			}

			return true;
		}

		template<class F>
		bool IsInBloomFilter(const ConcurrentBitVector & filter, std::vector<HashFunctionPtr> & hf, F f, char farg, uint64_t hash0)
		{
			for (size_t i = 0; i < hf.size(); i++)
			{
				uint64_t hvalue = i == 0 ? hash0 : ((*hf[i]).*f)(farg);
				if (!filter.Get(hvalue))
				{
					return false;
				}
			}

			return true;
		}

		class VertexLess
		{
		public:
			VertexLess(size_t vertexSize) : vertexSize_(vertexSize)
			{

			}

			bool operator () (const uint64_t & a, const uint64_t & b) const
			{
				DnaString stra(vertexSize_, a);
				DnaString strb(vertexSize_, b);
				return stra.GetBody() < strb.GetBody();
			}


		private:
			size_t vertexSize_;
		};

		size_t CharIndex(char ch)
		{
			return std::find(DnaString::LITERAL.begin(), DnaString::LITERAL.end(), ch) - DnaString::LITERAL.begin();
		}

		struct Task
		{
			bool isFinal;
			uint64_t start;
			std::string str;
#ifdef _DEBUG
			static const size_t TASK_SIZE = 36;
#else
			static const size_t TASK_SIZE = 1 << 18;
#endif			
			static const size_t GAME_OVER = SIZE_MAX;
			Task() {}
			Task(uint64_t start, bool isFinal, std::string && str) : start(start), isFinal(isFinal), str(std::move(str)) {}
		};

		DnaString kmer;

		const size_t QUEUE_CAPACITY = 16;
		typedef boost::lockfree::spsc_queue<Task> TaskQueue;
		typedef std::unique_ptr<TaskQueue> TaskQueuePtr;

		std::string TempFile(size_t record)
		{
			std::stringstream ss;
			ss << record << ".bin";
			return ss.str();
		}

		uint64_t NormHash(const std::vector<HashFunctionPtr> & posVertexHash, const std::vector<HashFunctionPtr> & negVertexHash)
		{
			return std::min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue);
		}

		bool Within(uint64_t hvalue, uint64_t low, uint64_t high)
		{
			return hvalue >= low && hvalue <= high;
		}

		DnaString Generate(std::string::const_iterator it, size_t size)
		{
			DnaString str;
			for (size_t i = 0; i < size; i++)
			{
				str.AppendBack(*(it + i));
			}

			return str;
		}

		std::string RevComp(const std::string & str)
		{
			std::string ret;
			for (auto it = str.rbegin(); it != str.rend(); ++it)
			{
				ret.push_back(DnaString::Reverse(*it));
			}

			return ret;
		}

		void InitializeHashFunctions(const std::vector<HashFunctionPtr> & seed,
			std::vector<HashFunctionPtr> & posEdgeHash,
			std::vector<HashFunctionPtr> & negEdgeHash,
			const std::string & fragment,
			size_t length,
			size_t offset = 0)
		{
			for (size_t i = 0; i < seed.size(); i++)
			{
				posEdgeHash[i] = HashFunctionPtr(new HashFunction(*seed[i]));
				negEdgeHash[i] = HashFunctionPtr(new HashFunction(*seed[i]));
				for (auto it = fragment.begin() + offset; it != fragment.begin() + length + offset; ++it)
				{
					posEdgeHash[i]->eat(*it);
				}

				assert(posEdgeHash[i]->hashvalue == posEdgeHash[i]->hash(fragment.substr(offset, length)));
				for (std::string::const_reverse_iterator it(fragment.begin() + length + offset); it != fragment.rend() - offset; ++it)
				{
					char c = DnaString::Reverse(*it);
					negEdgeHash[i]->eat(DnaString::Reverse(*it));
				}

				assert(negEdgeHash[i]->hashvalue == negEdgeHash[i]->hash(RevComp(fragment.substr(offset, length))));
			}
		}

		void CandidateCheckingWorker(uint64_t low,
			uint64_t high,
			const std::vector<HashFunctionPtr> & hashFunction,
			const ConcurrentBitVector & bitVector,
			size_t vertexLength,
			TaskQueue & taskQueue,
			FastSet & fastSet,
			boost::mutex & outMutex,
			size_t & total)
		{
			std::vector<uint64_t> output;
			while (true)
			{
				Task task;
				if (taskQueue.pop(task))
				{
					output.clear();
					if (task.start == Task::GAME_OVER)
					{
						break;
					}

					if (task.str.size() < vertexLength)
					{
						continue;
					}

					if (task.start == 0)
					{
						DnaString posVertex = Generate(task.str.begin(), vertexLength);
						DnaString negVertex = posVertex.RevComp();
						//	if (Within(NormHash(posVertexHash, negVertexHash), low, high))
						{
							fastSet.Insert(FastSet::Record(posVertex, negVertex, 'A', 'A'));
							fastSet.Insert(FastSet::Record(posVertex, negVertex, 'C', 'A'));
						}
					}

					if (task.isFinal)
					{
						DnaString posVertex = Generate(task.str.end() - vertexLength, vertexLength);
						DnaString negVertex = posVertex.RevComp();
						//	if (Within(NormHash(seed, posVertex, negVertex), low, high))
						{
							fastSet.Insert(FastSet::Record(posVertex, negVertex, 'A', 'A'));
							fastSet.Insert(FastSet::Record(posVertex, negVertex, 'C', 'A'));
						}
					}

					if (task.str.size() >= vertexLength + 2)
					{
						DnaString posVertex = Generate(task.str.begin() + 1, vertexLength);
						DnaString negVertex = posVertex.RevComp();
						std::vector<HashFunctionPtr> posVertexHash(hashFunction.size());
						std::vector<HashFunctionPtr> negVertexHash(hashFunction.size());
						InitializeHashFunctions(hashFunction, posVertexHash, negVertexHash, task.str, vertexLength, 1);
						for (size_t pos = 1;; ++pos)
						{
							char posPrev = task.str[pos - 1];
							char posExtend = task.str[pos + vertexLength];
							assert(posVertex.RevComp() == negVertex);
							//		if (Within(NormHash(posVertexHash, negVertexHash), low, high))
							{
								size_t inCount = 0;
								size_t outCount = 0;
								for (int i = 0; i < DnaString::LITERAL.size() && inCount < 2 && outCount < 2; i++)
								{
									char nextCh = DnaString::LITERAL[i];
									char revNextCh = DnaString::Reverse(nextCh);
									if (nextCh == posPrev)
									{
										++inCount;
									}
									else
									{
										uint64_t posHash0 = posVertexHash[0]->hash_prepend(nextCh);
										uint64_t negHash0 = negVertexHash[0]->hash_extend(revNextCh);
										assert(posHash0 == posVertexHash[0]->hash(std::string(1, nextCh) + task.str.substr(pos, vertexLength)));
										assert(negHash0 == negVertexHash[0]->hash(RevComp(std::string(1, nextCh) + task.str.substr(pos, vertexLength))));
										if ((posHash0 <= negHash0 && IsInBloomFilter(bitVector, posVertexHash, &HashFunction::hash_prepend, nextCh, posHash0)) ||
											(posHash0 > negHash0 && IsInBloomFilter(bitVector, negVertexHash, &HashFunction::hash_extend, revNextCh, negHash0)))
										{
											outCount++;
										}
									}

									if (nextCh == posExtend)
									{
										++outCount;
									}
									else
									{
										uint64_t posHash0 = posVertexHash[0]->hash_extend(nextCh);
										uint64_t negHash0 = negVertexHash[0]->hash_prepend(revNextCh);
										assert(posHash0 == posVertexHash[0]->hash(task.str.substr(pos, vertexLength) + nextCh));
										assert(negHash0 == negVertexHash[0]->hash(RevComp(task.str.substr(pos, vertexLength) + nextCh)));
										if ((posHash0 <= negHash0 && IsInBloomFilter(bitVector, posVertexHash, &HashFunction::hash_extend, nextCh, posHash0)) ||
											(posHash0 > negHash0 && IsInBloomFilter(bitVector, negVertexHash, &HashFunction::hash_prepend, revNextCh, negHash0)))
										{
											outCount++;
										}
									}
								}

								if (inCount > 1 || outCount > 1)
								{
									fastSet.Insert(FastSet::Record(posVertex, negVertex, posExtend, posPrev));
									if (posVertex == negVertex)
									{
										fastSet.Insert(FastSet::Record(posVertex, negVertex, DnaString::Reverse(posPrev), DnaString::Reverse(posExtend)));
									}
								}
							}

							if (pos + vertexLength + 1 < task.str.size())
							{
								char negExtend = DnaString::Reverse(posExtend);
								posVertex.AppendBack(posExtend);
								negVertex.AppendFront(negExtend);
								char posPrev = posVertex.PopFront();
								char negPrev = negVertex.PopBack();
								for (size_t i = 0; i < hashFunction.size(); i++)
								{
									posVertexHash[i]->update(posPrev, posExtend);
									negVertexHash[i]->reverse_update(negExtend, negPrev);
									assert(posVertexHash[i]->hashvalue == posVertexHash[i]->hash(posVertex.ToString()));
									assert(negVertexHash[i]->hashvalue == negVertexHash[i]->hash(negVertex.ToString()));
								}
							}
							else
							{
								break;
							}
						}
					}
				}
			}
		}

		void CountingWorkerOneRound(uint64_t low, uint64_t high, const std::vector<HashFunctionPtr> & hashFunction, ConcurrentBitVector & filter, size_t edgeLength, TaskQueue & taskQueue)
		{
			while (true)
			{
				Task task;
				if (taskQueue.pop(task))
				{
					if (task.start == Task::GAME_OVER)
					{
						break;
					}

					if (task.str.size() < edgeLength)
					{
						continue;
					}

					size_t vertexLength = edgeLength - 1;
					std::vector<HashFunctionPtr> posVertexHash(hashFunction.size());
					std::vector<HashFunctionPtr> negVertexHash(hashFunction.size());
					InitializeHashFunctions(hashFunction, posVertexHash, negVertexHash, task.str, vertexLength);
					for (size_t pos = 0; pos + edgeLength - 1 < task.str.size(); ++pos)
					{
						char nextCh = task.str[pos + edgeLength - 1];
						char revNextCh = DnaString::Reverse(nextCh);
						uint64_t posHash0 = posVertexHash[0]->hash_extend(nextCh);
						uint64_t negHash0 = negVertexHash[0]->hash_prepend(revNextCh);
						if (posHash0 < negHash0)
						{
							PutInBloomFilter(filter, posVertexHash, &HashFunction::hash_extend, nextCh, posHash0);
						}
						else
						{
							PutInBloomFilter(filter, negVertexHash, &HashFunction::hash_prepend, revNextCh, negHash0);
						}

						char prevCh = task.str[pos];
						for (size_t i = 0; i < hashFunction.size(); i++)
						{
							posVertexHash[i]->update(prevCh, nextCh);
							assert(posVertexHash[i]->hashvalue == posVertexHash[i]->hash(task.str.substr(pos + 1, vertexLength)));
							negVertexHash[i]->reverse_update(DnaString::Reverse(nextCh), DnaString::Reverse(prevCh));
							assert(negVertexHash[i]->hashvalue == negVertexHash[i]->hash(RevComp(task.str.substr(pos + 1, vertexLength))));
						}
					};
				}
			}
		}

		void DistributeTasks(const std::vector<std::string> & fileName, size_t overlapSize, std::vector<TaskQueuePtr> & taskQueue)
		{
			for (size_t file = 0; file < fileName.size(); file++)
			{
				size_t record = 0;
				const std::string & nowFileName = fileName[file];
				for (StreamFastaParser parser(nowFileName); parser.ReadRecord(); record++)
				{
					char ch;
					std::string buf;
					uint64_t prev = 0;
					uint64_t start = 0;
					bool over = false;
					do
					{
						over = !parser.GetChar(ch);
						if (!over)
						{
							start++;
							buf.push_back(ch);
						}

						if (buf.size() >= overlapSize && (buf.size() == Task::TASK_SIZE || over))
						{
							for (bool found = false; !found;)
							{
								for (TaskQueuePtr & q : taskQueue)
								{
									if (q->write_available() > 0)
									{
										std::string overlap;
										if (!over)
										{
											overlap.assign(buf.end() - overlapSize, buf.end());
										}

										q->push(Task(prev, over, std::move(buf)));
										prev = start - overlapSize + 1;
										buf.swap(overlap);
										found = true;
										break;
									}
								}
							}
						}

					} while (!over);
				}
			}

			for (size_t i = 0; i < taskQueue.size(); i++)
			{
				while (taskQueue[i]->write_available() == 0)
				{
					boost::this_thread::sleep_for(boost::chrono::nanoseconds(1000000));
				}

				taskQueue[i]->push(Task(Task::GAME_OVER, true, std::string()));
			}
		}
	}

	VertexEnumerator::VertexEnumerator(const std::vector<std::string> & fileName,
		size_t vertexLength,
		size_t filterSize,
		size_t hashFunctions,
		size_t rounds,
		size_t threads,
		size_t aggregationThreads,
		const std::string & tmpFileName) :
		vertexSize_(vertexLength)
	{
		uint64_t realSize = uint64_t(1) << filterSize;
		std::cout << "Threads = " << threads << std::endl;
		std::cout << "Aggregation threads = " << aggregationThreads << std::endl;
		std::cout << "Hash functions = " << hashFunctions << std::endl;
		std::cout << "Filter size = " << realSize << std::endl;
		std::cout << "Files: " << std::endl;
		for (const std::string & fn : fileName)
		{
			std::cout << fn << std::endl;
		}

		std::cout << std::string(80, '-') << std::endl;

		if (vertexLength > 30)
		{
			throw std::runtime_error("The vertex size is too large");
		}

		std::vector<HashFunctionPtr> hashFunction(hashFunctions);
		for (HashFunctionPtr & ptr : hashFunction)
		{
			ptr = HashFunctionPtr(new HashFunction(vertexLength, filterSize));
		}

		size_t edgeLength = vertexLength + 1;
		uint64_t low = 0;
		for (size_t round = 0; round < rounds; round++)
		{
			time_t mark = time(0);
			size_t totalRecords = 0;
			FastSet candidate(vertexLength, 1 << 4);
			uint64_t high = round == rounds - 1 ? UINT64_MAX : (UINT64_MAX / rounds) * (round + 1);
			{
				std::vector<TaskQueuePtr> taskQueue;
				ConcurrentBitVector bitVector(realSize);
				std::vector<boost::thread> workerThread(threads);
				std::cout << "Round " << round << ", " << low << ":" << high << std::endl;
				std::cout << "Counting\tEnumeration\tAggregation" << std::endl;
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					taskQueue.push_back(TaskQueuePtr(new TaskQueue(QUEUE_CAPACITY)));
					workerThread[i] = boost::thread(CountingWorkerOneRound,
						low,
						high,
						boost::cref(hashFunction),
						boost::ref(bitVector),
						edgeLength,
						boost::ref(*taskQueue[i]));
				}

				DistributeTasks(fileName, edgeLength, taskQueue);
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					workerThread[i].join();
				}

				std::cout << time(0) - mark << "\t";
				mark = time(0);
				boost::mutex tmpFileMutex;
				std::ofstream tmpFile(tmpFileName.c_str(), std::ios_base::binary);
				if (!tmpFile)
				{
					throw StreamFastaParser::Exception("Can't open the temporary file");
				}

				for (size_t i = 0; i < workerThread.size(); i++)
				{
					workerThread[i] = boost::thread(CandidateCheckingWorker,
						low,
						high,
						boost::cref(hashFunction),
						boost::cref(bitVector),
						vertexLength,
						boost::ref(*taskQueue[i]),
						boost::ref(candidate),
						boost::ref(tmpFileMutex),
						boost::ref(totalRecords));
				}

				DistributeTasks(fileName, vertexLength + 1, taskQueue);
				for (size_t i = 0; i < taskQueue.size(); i++)
				{
					workerThread[i].join();
				}
				

				std::cout << time(0) - mark << "\t";
			}

			mark = time(0);
			std::ifstream tmpFile(tmpFileName.c_str(), std::ios_base::binary);
			if (!tmpFile)
			{
				throw StreamFastaParser::Exception("Can't open the temporary file");
			}

			uint64_t falsePositives;
			candidate.DumpBifurcations(std::back_inserter(bifurcation_), falsePositives);
			std::cout << time(0) - mark << std::endl;
			std::cout << "Vertex count = " << bifurcation_.size() << std::endl;
			std::cout << "FP count = " << falsePositives << std::endl;
			std::cout << std::string(80, '-') << std::endl;
			low = high + 1;
		}

		tbb::parallel_sort(bifurcation_.begin(), bifurcation_.end());
	}

	size_t VertexEnumerator::GetVerticesCount() const
	{
		return bifurcation_.size();
	}

	size_t VertexEnumerator::GetId(const DnaString & vertex) const
	{
		DnaString check[2] = { vertex, vertex.RevComp() };
		for (DnaString str : check)
		{
			std::vector<uint64_t>::const_iterator it = std::lower_bound(bifurcation_.begin(), bifurcation_.end(), str.GetBody());
			if (it != bifurcation_.end() && *it == str.GetBody())
			{
				return it - bifurcation_.begin();
			}
		}

		return INVALID_VERTEX;
	}
}