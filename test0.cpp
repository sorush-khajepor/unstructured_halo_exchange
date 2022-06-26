#include "block2d.h"
#include "haloManager.h"
#include "communication.h"

void fillBlock(Block2d<double> &block, int rank)
{
  auto box = block.GetGhostBox();
  int i = rank * 1000;
  for (int ix = box.ix0; ix <= box.ix1; ix++)
    for (int iy = box.iy0; iy <= box.iy1; iy++)
      block(ix, iy) = i++;
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
      BoxRank{.box = b2, .rank = 2} /*,
      BoxRank{.box = b3, .rank = 2}, */
  };

  std::array<bool, 2> isPeriodic{true, true};
  HaloManager haloManager{boxRanks, isPeriodic, overlap};

  Block2d<double> block(boxRanks[rank].box, overlap);

  fillBlock(block, rank);

  auto neighbors = haloManager.FindNeighbours(block.GetGlobalOwnBox());

  Communication<double> comm(block, neighbors);
  comm.Communicate();

  std::cout << "Rank:" << rank << "\n";

  block.Print();
  std::cout << "\n";

  MPI_Finalize();

  return 0;
}