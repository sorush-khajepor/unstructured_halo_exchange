#include <array>
#include <iostream>

enum Axis
{
    x,
    y
};

struct Box
{
    int ix0 = 0, ix1 = 0, iy0 = 0, iy1 = 0;
    Box() = default;
    Box(int ix0_, int ix1_, int iy0_, int iy1_) : ix0{ix0_}, ix1{ix1_}, iy0{iy0_}, iy1{iy1_} {}
    std::array<int, 2> GetExtent() const
    {
        return std::array<int, 2>{ix1 - ix0 + 1, iy1 - iy0 + 1};
    }

    Box clone() const
    {
        return Box{*this};
    }

    bool IsEqual(const Box &box) const
    {
        return ix0 == box.ix0 && ix1 == box.ix1 &&
               iy0 == box.iy0 && iy1 == box.iy1;
    }

    Box &Translate(int dx, int dy)
    {
        ix0 += dx;
        ix1 += dx;
        iy0 += dy;
        iy1 += dy;
        return *this;
    }
    Box &Expand(int thickness, Axis axis)
    {
        ix0 -= thickness * (1 - axis);
        ix1 += thickness * (1 - axis);
        iy0 -= thickness * axis;
        iy1 += thickness * axis;
        return *this;
    }
    Box &Expand(int thickness)
    {
        return Expand(thickness, Axis::x).Expand(thickness, Axis::y);
    }

    Box &Merge(const Box &other)
    {
        ix0 = std::min(ix0, other.ix0);
        ix1 = std::max(ix1, other.ix1);
        iy0 = std::min(iy0, other.iy0);
        iy1 = std::max(iy1, other.iy1);
        return *this;
    }

    Box &Intersect(const Box &other)
    {
        ix0 = std::max(ix0, other.ix0);
        ix1 = std::min(ix1, other.ix1);
        iy0 = std::max(iy0, other.iy0);
        iy1 = std::min(iy1, other.iy1);
        return *this;
    }
    bool IsValid() const
    {
        return ix0 <= ix1 && iy0 <= iy1;
    }
    bool Contains(int ix, int iy) const
    {
        return ix0 <= ix && iy0 <= iy &&
               ix <= ix1 && iy <= iy1;
    }
    bool Contains(const Box &box) const
    {
        return Contains(box.ix0, box.iy0) && Contains(box.ix1, box.iy1);
    }

    void Print() const
    {
        std::cout << "ixmin=" << ix0 << ", iymin=" << iy0 << "\n";
        std::cout << "ixmax=" << ix1 << ", iymax=" << iy1 << "\n";
    }
};