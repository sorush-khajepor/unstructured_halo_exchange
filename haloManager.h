#include <mpi.h>
#include "block2d.h"
#include <vector>

namespace Dir
{
    enum Type
    {
        center,
        right,
        top,
        left,
        bottom,
        topright,
        topleft,
        bottomleft,
        bottomright
    };
}
std::array<std::array<int, 2>, 9> Dir2Array{{{0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1}, {1, 1}, {-1, 1}, {-1, -1}, {1, -1}}};
int IsXPeriodicDir[9]{0, 1, 0, 1, 0, 0, 0, 0, 0};
int IsYPeriodicDir[9]{0, 0, 1, 0, 1, 0, 0, 0, 0};
int IsXyPeriodicDir[9]{0, 1, 1, 1, 1, 1, 1, 1, 1};

struct BoxRank
{
    Box box;
    int rank;
};

struct Neighbor
{
    size_t rank;
    Box box;
    Dir::Type dir;
    MPI_Datatype inType, outType;
    int inTag, outTag;
};

class HaloManager
{

private:
public:
    size_t overlap;
    std::vector<Neighbor> neighbors;
    // this rank block boxes
    Box box;
    Box ghostBox;
    Box localGhostBox;
    Box domainBox;
    int isPeriodicDir[9];

    int rank, size;

    HaloManager(const std::vector<Box> &domainBoxes_, std::array<bool, 2> isPeriodic, size_t overlap_) : overlap{overlap_}
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        SetIsPeriodicDir(isPeriodic);
        std::vector<BoxRank> boxRank;
        for (size_t i = 0; i < domainBoxes_.size(); i++)
        {
            boxRank.push_back(BoxRank{.box = domainBoxes_[i], .rank = i});
        }

        // this rank domain box
        box = domainBoxes_[rank];
        ghostBox = ToGhostBox(box);
        localGhostBox = RelativeToGhostBox(ghostBox);
        SetDomainBox(boxRank);
        SetNeighbours(boxRank);
    }
    void SetIsPeriodicDir(std::array<bool, 2> periodic)
    {
        if (periodic[0] && !periodic[1])
            for (size_t i = 0; i < 9; i++)
                isPeriodicDir[i] = IsXPeriodicDir[i];
        else if (!periodic[0] && periodic[1])
            for (size_t i = 0; i < 9; i++)
                isPeriodicDir[i] = IsYPeriodicDir[i];
        else if (periodic[0] && periodic[1])
            for (size_t i = 0; i < 9; i++)
                isPeriodicDir[i] = IsXyPeriodicDir[i];
    }
    void SetDomainBox(std::vector<BoxRank> &boxRanks)
    {
        domainBox = boxRanks[0].box;
        for (size_t i = 1; i < boxRanks.size(); i++)
            domainBox.Merge(boxRanks[i].box);
    }

    Box GetPeriodicBox(Box box, Dir::Type dir, bool &isThere)
    {
        if (!isPeriodicDir[dir])
        {
            isThere = false;
            return box;
        }

        isThere = true;
        auto dirVec = Dir2Array[dir];
        auto [lx, ly] = domainBox.GetExtent();
        return box.Translate(dirVec[0] * lx, dirVec[1] * ly);
    }

    Box PeriodicBoxToInDomainBox(const Box &box)
    {
        auto [lx, ly] = domainBox.GetExtent();
        return Box{(box.ix0 + lx) % lx, (box.ix1 + lx) % lx,
                   (box.iy0 + ly) % ly, (box.iy1 + ly) % ly};
    }
    int HashBox(const Box &box)
    {
        return box.ix0 + box.ix1 + box.iy0 + box.iy1;
    }

    void SetNeighbours(const std::vector<BoxRank> &boxRanks)
    {
        for (auto &boxRank : boxRanks)
        {
            if (boxRank.rank != rank)
                AddNeighbour(boxRank);

            for (size_t dir = 0; dir < 9; dir++)
            {
                bool isThere = false;
                auto pbox = GetPeriodicBox(boxRank.box, (Dir::Type)dir, isThere);
                if (isThere)
                    AddNeighbour(BoxRank{pbox, boxRank.rank});
            }
        }
    }

    void AddNeighbour(const BoxRank &boxRank)
    {
        auto inbox = ghostBox.clone().Intersect(boxRank.box);
        auto outbox = box.clone().Intersect(ToGhostBox(boxRank.box));

        if (inbox.IsValid() && outbox.IsValid())
        {
            auto intag = HashBox(PeriodicBoxToInDomainBox(inbox));
            auto outtag = HashBox(PeriodicBoxToInDomainBox(outbox));
            auto localInbox = RelativeToGhostBox(inbox);
            auto localOutbox = RelativeToGhostBox(outbox);
            neighbors.push_back(Neighbor{
                .rank = boxRank.rank,
                .box = PeriodicBoxToInDomainBox(boxRank.box),
                .inType = MakeMpiType(localGhostBox, localInbox),
                .outType = MakeMpiType(localGhostBox, localOutbox),
                .inTag = intag,
                .outTag = outtag});
        }
    }

    MPI_Datatype MakeMpiType(const Box &box, const Box &subBox)
    {
        MPI_Datatype type;
        auto size = box.GetExtent();
        auto subSize = subBox.GetExtent();
        int subStart[2]{subBox.ix0, subBox.iy0};
        MPI_Type_create_subarray(2, size.data(), subSize.data(), subStart,
                                 MPI_ORDER_C, MPI_DOUBLE, &type);
        MPI_Type_commit(&type);
        return type;
    }

    Box ToGhostBox(const Box &box)
    {
        return box.clone().Expand(overlap);
    }

    Box RelativeToGhostBox(const Box &box)
    {
        return box.clone().Translate(-ghostBox.ix0, -ghostBox.iy0);
    }
};