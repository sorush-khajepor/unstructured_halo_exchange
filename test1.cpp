#include "cpuBlocks.h"

void fillBlock(Block2d<double> &block, int start)
{
    auto box = block.GetGhostBox();
    for (int ix = box.ix0; ix <= box.ix1; ix++)
        for (int iy = box.iy0; iy <= box.iy1; iy++)
            block(ix, iy) = start++;
}

int main()
{

    MPI_Init(NULL, NULL);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /*
      size_t overlap = 2;
      Box b0{0, 5, 0, 10};
      Box b1{6, 10, 0, 10};
      Box b2{0, 5, 11, 15};
      Box b3{6, 10, 11, 15};

      std::vector<Box> boxes{b0, b1, b2, b3};
    */
    /*
      size_t overlap = 2;
      Box b0{0, 20, 0, 7};
      Box b1{0, 5,  8, 14};
      Box b2{15, 20, 8, 14};

      std::vector<Box> boxes{b0, b1, b2};
    */

    size_t overlap = 2;
    Box b0{0, 4, 0, 7};
    Box b1{5, 10, 0, 7};
    Box b2{11, 16, 0, 7};
    Box b3{17, 22, 0, 7};

    std::vector<BoxRank> boxRanks{
        BoxRank{.box = b0, .rank = 0},
        BoxRank{.box = b1, .rank = 1},
        BoxRank{.box = b2, .rank = 2},
        BoxRank{.box = b3, .rank = 2},
    };

    std::array<bool, 2> isPeriodic{true, true};
    
    CpuBlocks<double> cpuBlocks{boxRanks, isPeriodic, overlap};
    
    for (size_t i = 0; i < cpuBlocks.blocks.size(); i++)
    {
        fillBlock(cpuBlocks.blocks[i], 10000*rank+1000*i);
    }
    
    cpuBlocks.Communicate();
    std::cout << "Rank:" << rank << "\n";
    for (auto &&block : cpuBlocks.blocks)
        block.Print();
    std::cout << "\n";

    MPI_Finalize();

    return 0;
}