//
//  Algorithm.cpp
//  vdbTest
//
//  Created by Hanjie Liu on 8/11/17.
//  Copyright © 2017 Hanjie Liu. All rights reserved.
//

#include "Algorithm.hpp"
#include <openvdb/tools/LevelSetFilter.h>

using namespace std;
using namespace openvdb;

bool PFIAir::algorithme::inflate(const float offset, FloatGrid::Ptr) {
    return true;
}

void PFIAir::algorithme::printValues(FloatGrid::Ptr p) {
//    for (FloatGrid::ValueAllIter iter = p -> beginValueAll() ; iter; ++iter) {
//        cout << iter.getValue() << endl;
//    }
    
    for (FloatGrid::ValueOnCIter iter = p -> cbeginValueOn() ; iter; ++iter) {
        cout << iter.getValue() << endl;
    }
}

void PFIAir::algorithme::changeActiveVoxelValues(FloatGrid::Ptr p, float offset) {
    for (FloatGrid::ValueOnIter iter = p -> beginValueOn() ; iter; ++iter) {
        float dist = iter.getValue();
        iter.setValue(dist - offset);
    }
}

//void skeletonization(FloatGrid& p) {
//    using namespace openvdb;
//
//    tools::LevelSetFilter<FloatGrid> tracker = tools::LevelSetFilter<FloatGrid>(p);
//    
//}


// write loops to check signed vs unsigned sizes
// interseting voxel or bounding box
// corner 
//'FloatGrid::Ptr' (aka 'std::__1::shared_ptr<openvdb::v4_0_2::Grid<openvdb::v4_0_2::tree::Tree<openvdb::v4_0_2::tree::RootNode<openvdb::v4_0_2::tree::InternalNode<openvdb::v4_0_2::tree::InternalNode<openvdb::v4_0_2::tree::LeafNode<float, 3>, 4>, 5> > > > >')
//
//'tools::LevelSetTracker<FloatGrid>' (aka 'LevelSetTracker<Grid<Tree<RootNode<InternalNode<InternalNode<LeafNode<float, 3U>, 4U>, 5U> > > > >')

