#include <mpi.h>
#include "block2d.h"
#include <vector>

enum Dir
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

struct BoxRank
{
    Box box;
    int rank;
};

struct Neighbor
{
    size_t rank;
    Box box;
    Dir dir;
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

    int rank, size;

    HaloManager(const std::vector<Box> &domainBoxes_, size_t overlap_) : overlap{overlap_}
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        bool periodic[2]{true, true};
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
        AddPeriodicBoxes(boxRank, periodic);
        SetNeighbours(boxRank);
    }
    void SetDomainBox(std::vector<BoxRank> &boxRanks) {
        domainBox = boxRanks[0].box;
        for (size_t i = 1; i < boxRanks.size(); i++)
            domainBox.Merge(boxRanks[i].box);
    }
    void AppendBoxRanks(std::vector<BoxRank> &boxRanks, const std::vector<BoxRank> &other)
    {
        boxRanks.insert(boxRanks.end(), other.begin(), other.end());
    }
    auto Translate(std::vector<BoxRank> boxRanks, int dx, int dy)
    {
        for (auto &&item : boxRanks)
        {
            item.box.Translate(dx, dy);
        }
        return boxRanks;
    }
    void AddPeriodicBoxes(std::vector<BoxRank>& boxRanks, bool periodic[2])
    {
        auto [dx, dy] = domainBox.GetExtent();
        auto right = Translate(boxRanks, dx, 0);
        auto left = Translate(boxRanks, -dx, 0);

        auto top = Translate(boxRanks, 0, dy);
        auto bottom = Translate(boxRanks, 0, -dy);

        auto topright = Translate(boxRanks, dx, dy);
        auto topleft = Translate(boxRanks, -dx, dy);
        auto bottomleft = Translate(boxRanks, -dx, -dy);
        auto bottomright = Translate(boxRanks, dx, -dy);

        if (periodic[0])
        {
            AppendBoxRanks(boxRanks, right);
            AppendBoxRanks(boxRanks, left);
        }
        if (periodic[1])
        {
            AppendBoxRanks(boxRanks, top);
            AppendBoxRanks(boxRanks, bottom);
        }

        if (periodic[0] && periodic[1])
        {
            AppendBoxRanks(boxRanks, topright);
            AppendBoxRanks(boxRanks, topleft);
            AppendBoxRanks(boxRanks, bottomleft);
            AppendBoxRanks(boxRanks, bottomright);
        }
    }
    Box PeriodicBoxToInDomainBox(const Box& box){
            auto [dx,dy] = domainBox.GetExtent();
            return Box{(box.ix0+dx)%dx, (box.ix1+dx)%dx , (box.iy0+dy)%dy, (box.iy1+dy)%dy};
    }
    int HashBox(const Box& box){
             return box.ix0 + box.ix1 + box.iy0 + box.iy1;
    }
    
    void SetNeighbours(const std::vector<BoxRank> &boxRanks)
    {
        for (size_t i = 0; i < boxRanks.size(); i++)
        {
            auto& boxRank=boxRanks[i];

            if (boxRank.rank == rank && i<size)
                continue;
            auto inbox = ghostBox.clone().Intersect(boxRank.box);
            auto outbox = box.clone().Intersect(ToGhostBox(boxRank.box));
            auto intag = HashBox(PeriodicBoxToInDomainBox(inbox));
            auto outtag = HashBox(PeriodicBoxToInDomainBox(outbox));

            if (inbox.IsValid() && outbox.IsValid())
            {
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
    }

    MPI_Datatype MakeMpiType(const Box &box, const Box &subBox)
    {
        MPI_Datatype type;
        auto size = box.GetExtent();
        auto subSize = subBox.GetExtent();
        int subStart[2]{subBox.ix0, subBox.iy0};
        MPI_Type_create_subarray(2, size.data(), subSize.data(), &subStart[0],
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