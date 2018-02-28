//
//  Container.hpp
//  vdbTest
//
//  Created by Hanjie Liu on 7/28/17.
//  Copyright © 2017 Hanjie Liu. All rights reserved.
//

#ifndef IOManager_hpp
#define IOManager_hpp

#include <stdio.h>
#include <string>
#include <openvdb/openvdb.h>

namespace PFIAir {
    using namespace openvdb;
    using std::string;
    using std::vector;
    
    class Container{
    private:
        string _filename;
        
        vector<Vec3s> _points = vector<Vec3s>();
        vector<Vec3I> _indicesTri = vector<Vec3I>();
        vector<Vec4I> _indicesQuad = vector<Vec4I>();
        
        math::Transform _scale = math::Transform();
        
    protected:
        double getEuclideanFromIndices(int p1, int p2);
        
    public:
        Container();
        
        string getModelName() {return _filename;}
        
        // i/o
        void loadMeshModel(const string filename);
        static void exportModel(const string name, FloatGrid::Ptr model);
        
        // conversion
        FloatGrid::Ptr getWaterTightLevelSet();
        FloatGrid::Ptr getWaterTightLevelSetWithBandWidth(float w);
        FloatGrid::Ptr getUnsignedDistanceField(float bandwidth);
        
        // scale
        void setScale(Vec3d scale);
        float computeAverageEdgeLength();
        
        // center
        void computeMeshCenter();
        void OLDcomputeMeshCenter();

    };
}

#endif /* IOManager_hpp */
