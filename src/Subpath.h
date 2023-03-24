//
// Created by creeper on 23-3-24.
//
#pragma once
#include <vector>
#include "intersection.h"
struct Subpath
{
    std::vector<PathVertex> path;
    void append(const PathVertex& pv) { path.push_back(pv); }
    uint length() const { return path.size(); }
    PathVertex& operator[](int i) { return path[i]; }
    PathVertex operator[](int i) const { return path[i]; }
};

