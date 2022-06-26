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
    Box domainBox;
    std::vector<BoxRank> boxRanks;

    int isPeriodicDir[9];

    HaloManager(const std::vector<BoxRank> &boxRanks_, std::array<bool, 2> isPeriodic, size_t overlap_) : boxRanks{boxRanks_}, overlap{overlap_}
    {
        SetIsPeriodicDir(isPeriodic);
        SetDomainBox(boxRanks);
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

    auto FindNeighbours(const Box& fbox)
    {
        std::vector<Neighbor> neighbors;
        for (auto &boxRank : boxRanks)
        {
            if (!fbox.IsEqual(boxRank.box))
                AddNeighbour(fbox, boxRank, neighbors);

            for (size_t dir = 0; dir < 9; dir++)
            {
                bool isThere = false;
                auto pbox = GetPeriodicBox(boxRank.box, (Dir::Type)dir, isThere);
                if (isThere)
                    AddNeighbour(fbox, BoxRank{pbox, boxRank.rank}, neighbors);
            }
        }
        return neighbors;
    }

    void AddNeighbour(const Box& fbox, const BoxRank &boxRank, std::vector<Neighbor>& neighbors )
    {
        auto ghostBox = ToGhostBox(fbox);
        auto localGhostBox = RelativeToGhostBox(fbox, ghostBox);
        auto inbox = ghostBox.clone().Intersect(boxRank.box);
        auto outbox = fbox.clone().Intersect(ToGhostBox(boxRank.box));

        if (inbox.IsValid() && outbox.IsValid())
        {
            auto intag = HashBox(PeriodicBoxToInDomainBox(inbox));
            auto outtag = HashBox(PeriodicBoxToInDomainBox(outbox));
            auto localInbox = RelativeToGhostBox(fbox,inbox);
            auto localOutbox = RelativeToGhostBox(fbox,outbox);
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

    Box RelativeToGhostBox(const Box& fbox,const Box &box)
    {
        auto ghostBox = ToGhostBox(fbox);
        return box.clone().Translate(-ghostBox.ix0, -ghostBox.iy0);
    }
};