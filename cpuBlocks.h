#pragma once
#include "block2d.h"
#include "haloManager.h"
#include "communication.h"

template <class T>
class CpuBlocks
{
private:
    std::vector<Communication<T>> comms;

public:
    std::vector<Block2d<T>> blocks;

    CpuBlocks(const std::vector<BoxRank> &boxRanks, std::array<bool, 2> isPeriodic, size_t overlap)
    {
        HaloManager haloManager{boxRanks, isPeriodic, overlap};

        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        for (auto &&boxRank : boxRanks)
            if (boxRank.rank == rank)
                blocks.emplace_back(boxRank.box, overlap);

        for (auto &&block : blocks)
        {
            auto neighbors = haloManager.FindNeighbours(block.GetGlobalOwnBox());
            comms.emplace_back(block, neighbors);
        }
    };

    void CommunicateAsync()
    {
        for (auto &&comm : comms)
            comm.CommunicateAsync();
    }

    void Await()
    {
        for (auto &&comm : comms)
            comm.Await();
    }

    void Communicate()
    {
        CommunicateAsync();
        Await();
    }
};
