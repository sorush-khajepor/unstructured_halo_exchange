
#pragma once
#include <vector>
#include <array>
#include "box2d.h"
#include <iostream>
#include <iomanip>

template <class T>
class Block2d
{
private:
    std::vector<T> nodes;

public:
    size_t sizeX, sizeY;
    int globalX0, globalY0;
    size_t overlap;

    Block2d() = delete;
    Block2d(const Box& box, size_t overlap_) : overlap{overlap_}
    {
        auto extent = box.GetExtent();
        sizeX = extent[0]+2*overlap;
        sizeY = extent[1]+2*overlap;
        globalX0 = box.ix0 - overlap;
        globalY0 = box.iy0 - overlap;
        nodes.assign(sizeX * sizeY, T{});
    }
    auto assign(std::vector<T> nodes_)
    {
        nodes = nodes_;
        if (sizeX * sizeY != nodes.size())
            throw;
    }
    auto cartToIndex(size_t ix, size_t iy)
    {
        return iy + ix * sizeY;
    }
    auto &operator()(size_t ix, size_t iy)
    {
        return nodes[cartToIndex(ix, iy)];
    }
    auto GetGhostBox() const
    {
        return Box{0, sizeX - 1, 0, sizeY - 1};
    }
    auto GetGlobalGhostBox() const
    {
        return GetGhostBox().Translate(globalX0,globalY0);
    }
    auto GetOwnBox() const { return GetGhostBox().Expand(-overlap); }
    auto GetGlobalOwnBox() const { return GetGlobalGhostBox().Expand(-overlap); }
    auto ToString(const Box& box)
    {
        std::ostringstream s;
        for (int iy = box.iy1; iy >= box.iy0   ; iy--)
        {
            s << "\n";
            for (int ix = box.ix0; ix <= box.ix1; ix++)
            {
                s << std::setw(5) << operator()(ix, iy) << " ";
            }
        }
        s << "\niy↑  ix-> \n";
        return s.str();
    }
    auto Print(const Box& box)
    {
        std::cout << ToString(box);
    }
    auto Print(){ Print(GetGhostBox());}
};