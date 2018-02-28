//
//  Algorithm.hpp
//  vdbTest
//
//  Created by Hanjie Liu on 8/11/17.
//  Copyright © 2017 Hanjie Liu. All rights reserved.
//

#ifndef Algorithm_hpp
#define Algorithm_hpp

#include <stdio.h>
#include <stdio.h>
#include <string>
#include <openvdb/openvdb.h>

namespace PFIAir {
    namespace algorithme {
        using namespace openvdb;
        
        bool inflate(const float offset, FloatGrid::Ptr);
        void changeActiveVoxelValues(FloatGrid::Ptr p, float offset);
        
        void skeletonization(FloatGrid::Ptr p);
        
        // debug
        void printValues(FloatGrid::Ptr);
    }
}

#endif /* Algorithm_hpp */
