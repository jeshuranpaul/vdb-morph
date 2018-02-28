#ifndef GRIDOPERATIONS_H
#define GRIDOPERATIONS_H

#include <openvdb/util/NullInterrupter.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <fstream>

#include "CommonOperations.h"

namespace GridOperations {
    
    struct MatMul {
        openvdb::math::Mat3s M;
        MatMul(const openvdb::math::Mat3s& mat): M(mat) {}
        inline void operator()(const openvdb::Vec3SGrid::ValueOnIter& iter) const {
            try {
                iter.setValue(M.transform(*iter));
            }
            catch(...) {
                
            }
        }
    };
    
    openvdb::GridBase::Ptr readFile(const std::string);
    void writeToFile(const std::string, openvdb::FloatGrid::Ptr);
    void writeOnlySurface(const openvdb::FloatGrid::Ptr);
    bool checkIfSurface(openvdb::FloatGrid::ValueOnIter, const openvdb::FloatGrid::Ptr);
    
    double measureGrid(const openvdb::FloatGrid::Ptr);
    double calcPlatonicScale(int, double);
    
    void writeActiveVoxelVals();
    void setDefaultBackgroundValue(openvdb::FloatGrid::Ptr&);
    
    
    openvdb::GridBase::Ptr readFile(const std::string file_name) {
        openvdb::io::File file(file_name);
        file.open();
        openvdb::GridBase::Ptr baseGrid;
        for (openvdb::io::File::NameIterator nameIter = file.beginName(); nameIter != file.endName(); ++nameIter)
        {
            baseGrid = file.readGrid(nameIter.gridName());
        }
        file.close();
        return baseGrid;
    }
    
    void writeToFile(const std::string file_name, openvdb::FloatGrid::Ptr grid_pointer) {
        grid_pointer->setGridClass(openvdb::GRID_LEVEL_SET);
        openvdb::io::File file(file_name);
        openvdb::GridPtrVec grids;
        grids.push_back(grid_pointer);
        file.write(grids);
        file.close();
    }
    
    void writeOnlySurface(const openvdb::FloatGrid::Ptr grid_pointer) {
        double bg = grid_pointer->background();
        int surface = 0, non_surface = 0;
        std::cout << grid_pointer->activeVoxelCount() << std::endl;
        for (openvdb::FloatGrid::ValueOnIter iter = grid_pointer->beginValueOn(); iter; ++iter) {
            if(iter.getValue() < 0) {
                if(!GridOperations::checkIfSurface(iter, grid_pointer)) {
                    ++non_surface;
                    iter.setValue(bg);
                }
                else ++surface;
            }
        }
        std::cout << surface << std::endl;
        std::cout << non_surface << std::endl;
        GridOperations::writeToFile("only_surface.vdb", grid_pointer);
    }
    
    bool checkIfSurface(openvdb::FloatGrid::ValueOnIter iterator, const openvdb::FloatGrid::Ptr grid_pointer) {
        bool found_positive = false;
        openvdb::Coord coord = iterator.getCoord();
        openvdb::Coord six_connected[6] = {
            openvdb::Coord(coord.x() - 1, coord.y(), coord.z()),
            openvdb::Coord(coord.x() + 1, coord.y(), coord.z()),
            openvdb::Coord(coord.x(), coord.y() - 1, coord.z()),
            openvdb::Coord(coord.x(), coord.y() + 1, coord.z()),
            openvdb::Coord(coord.x(), coord.y(), coord.z() - 1),
            openvdb::Coord(coord.x(), coord.y(), coord.z() + 1)
        };
        
        openvdb::FloatGrid::Accessor accessor = grid_pointer->getAccessor();
        
        for(int i = 0; i < 6; i++) {
            if(accessor.getValue(six_connected[i]) > 0) {
                found_positive = true;
                break;
            }
        }
        return found_positive;
    }
    
    double measureGrid(const openvdb::FloatGrid::Ptr grid_pointer) {
        openvdb::util::NullInterrupter* interrupt = nullptr;
        openvdb::tools::LevelSetMeasure <openvdb::FloatGrid, openvdb::util::NullInterrupter> ls_measure(*(grid_pointer), interrupt);
        
        double area = 0.0, volume = 0.0;
        ls_measure.measure(area, volume);
        
        return area;
    }
    
    double calcPlatonicScale(int face, double radius) {
        double _pi = 3.14159,
        sa_sphere = 4 * _pi * openvdb::math::Pow(radius, 2),
        ret_scale = 1.0;
        
        double openvdb_platonic_def_side_length_arr[5] = {1.632991, 1.0, 1.414214, 0.708876, 1.051463};
        
        switch(face) {
            case 4:
                ret_scale = openvdb::math::Sqrt(sa_sphere / (openvdb::math::Sqrt(3.0) * openvdb::math::Pow(openvdb_platonic_def_side_length_arr[0], 2)));
                break;
            case 6:
                ret_scale =  openvdb::math::Sqrt(sa_sphere / (6 * openvdb::math::Pow(openvdb_platonic_def_side_length_arr[1], 2)));
                break;
            case 8:
                ret_scale = openvdb::math::Sqrt(sa_sphere / (2 * openvdb::math::Sqrt(3.0) * openvdb::math::Pow(openvdb_platonic_def_side_length_arr[2], 2)));
                break;
            case 12:
                ret_scale = openvdb::math::Sqrt(sa_sphere / (20.6457288071 * openvdb::math::Pow(openvdb_platonic_def_side_length_arr[3], 2)));
                break;
            case 20:
                ret_scale = openvdb::math::Sqrt(sa_sphere / (5 * openvdb::math::Sqrt(3.0) * openvdb::math::Pow(openvdb_platonic_def_side_length_arr[4], 2)));
                break;
        }
        return ret_scale;
    }
    
    void writeActiveVoxelVals() {
        double sphere_radius = 5.0, scale = 5.0, voxel_size = 0.03, half_width = 30.0;
        openvdb::FloatGrid::Ptr grids[6] = {CommonOperations::getSphereVolume(sphere_radius, voxel_size, half_width),
            CommonOperations::getPlatonicVolume(4, scale, voxel_size, half_width),
            CommonOperations::getPlatonicVolume(6, scale, voxel_size, half_width),
            CommonOperations::getPlatonicVolume(8, scale, voxel_size, half_width),
            CommonOperations::getPlatonicVolume(12, scale, voxel_size, half_width),
            CommonOperations::getPlatonicVolume(20, scale, voxel_size, half_width)};
        
        std::string labels[6] = {"Sphere", "Tetra", "Cube", "Octa", "Dodec", "Icoso"};
        std::string filename = "";
        std::ofstream file;
        openvdb::Coord curr_coord;
        
        file.open ("active_vals.txt");
        file << "Sphere radius - " << sphere_radius << "\n";
        file << "Scale - " << scale << "\n";
        file << "Voxel size - " << voxel_size << "\n";
        file << "Half width - " << half_width << "\n\n";
        
        for(int i=0; i < 6; i++) {
            openvdb::FloatGrid::Ptr grid_pointer = grids[i];
            filename = "output/stock_vdbs/r" + std::to_string((int)sphere_radius) + "s" + std::to_string((int)scale) + "/" + labels[i] + ".vdb";
            GridOperations::writeToFile(filename, grid_pointer);
            
            file << "For level set " << labels[i] << "\n";
            for (openvdb::FloatGrid::ValueOnIter iter = grid_pointer->beginValueOn(); iter; ++iter) {
                curr_coord = iter.getCoord();
                if(curr_coord.z() == 0)
                    file << "At Coord: " << curr_coord << " = " << iter.getValue() << "\n";
            }
            file << "\n\n\n";
        }
        file.close();
    }
    
    void setDefaultBackgroundValue(openvdb::FloatGrid::Ptr& grid_pointer) {
        for (openvdb::FloatGrid::ValueOffIter iter = grid_pointer->beginValueOff(); iter; ++iter) {
            if(iter.getValue() < 0) iter.setValue(-1.0);
            else iter.setValue(1.0);
        }
    }
}

#endif