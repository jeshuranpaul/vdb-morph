#ifndef TESTS_H
#define TESTS_H

#include <fstream> 
#include "MorphOperations.h"
#include "CommonOperations.h"

namespace Tests {
    void testingOnOffVoxels(openvdb::FloatGrid::Ptr& grid_pointer) {
        double bg = grid_pointer->background();
        openvdb::Vec3d voxel_size = grid_pointer->voxelSize();
        
        int count = 0, less = 0, more = 0;
        
        count = 0, less = 0, more = 0;
        for (openvdb::FloatGrid::ValueOnIter iter = grid_pointer->beginValueOn(); iter; ++iter) {
            ++count;
            if(iter.getValue() < voxel_size[0] + 0.1) ++less;
            else ++more;
        }
        CommonOperations::display("count on", count);
        CommonOperations::display("less", less);
        CommonOperations::display("more", more);
        
        
        
        for (openvdb::FloatGrid::ValueAllIter iter = grid_pointer->beginValueAll(); iter; ++iter) {
            ++count;
            if(iter.getValue() < -0.3 ) {
                ++less;
                iter.setValue(bg);
            }
            else ++more;
        }
        CommonOperations::display("count all", count);
        CommonOperations::display("less", less);
        CommonOperations::display("more", more);
        
        
        
        count = 0, less = 0, more = 0;
        for (openvdb::FloatGrid::ValueOffIter iter = grid_pointer->beginValueOff(); iter; ++iter) {
            ++count;
            if(iter.getValue() < 0) ++less;
            else ++more;
        }
        CommonOperations::display("count off", count);
        CommonOperations::display("less", less);
        CommonOperations::display("more", more);
        
        count = 0, less = 0, more = 0;
        for (openvdb::FloatGrid::ValueOnIter iter = grid_pointer->beginValueOn(); iter; ++iter) {
            ++count;
            if(iter.getValue() < 0) ++less;
            else ++more;
        }
        CommonOperations::display("count on", count);
        CommonOperations::display("less", less);
        CommonOperations::display("more", more);
        
        GridOperations::writeToFile("settobackround.vdb", grid_pointer);
    }
    
    void setBackgroundTest() {
        openvdb::FloatGrid::Ptr source_grid1 = CommonOperations::getSphereVolume(1.0, 0.03, 3);
        openvdb::FloatGrid::Accessor acc = source_grid1->getAccessor();
        double value;
        double bg = source_grid1->background();
        
        openvdb::Coord xyz;
        
        int min_index = -100, max_index = 100;
        
        for(int i = -100; i < 100; i++)
            for(int j = -100; j < 0; j++)
                for(int k = -100; k < 100; k++) {
                    xyz = openvdb::Coord(i, j, k);
                    
                    acc.setValueOn(xyz, bg);
                }
        
        
        GridOperations::writeToFile("grid.vdb", source_grid1);
        
    }
    
    void runExperiments() { 
        double sphere_radius = 1.0, scale = 1.0, voxel_size = 0.03, half_width = 3.0;
        openvdb::FloatGrid::Ptr grid_pointer = CommonOperations::getPlatonicVolume(8, scale, voxel_size, half_width);
        std::ofstream file;
        file.open ("bg_vals.txt");
        
        for (openvdb::FloatGrid::ValueOffIter iter = grid_pointer->beginValueOff(); iter; ++iter) {
            file << iter.getValue() << "\n";
        }
        file.close();
        GridOperations::writeToFile("bg_changed.vdb", grid_pointer);
    }
    
    void old_script() {
        std::ofstream file;
        
        double dt_array[1] = {0.1};
        std::string dt_arr_str[1] = {"0.1"};
        
        double voxel_size_arr[1] = {0.025};
        std::string voxel_size_str_arr[1] = {"0.025"};
        
        int faces_arr[5] = {4, 6, 8, 12, 20};
        double sphere_radius = 1.0;
        double source_half_width = 3.0, target_half_width = 40.0;
        std::string table = "", path = "";
        
        for(int k = 0; k < sizeof(dt_array)/sizeof(dt_array[0]); k++) {                     //dt
            for(int j = 0; j < sizeof(voxel_size_arr)/sizeof(voxel_size_arr[0]); j++) {     //voxel size
                
                file.open ("log.txt");
                file << "dt - " << dt_arr_str[k] << "\n";
                file << "Voxel size - " << voxel_size_str_arr[j] << "\n";
                file << "Source hw - " << source_half_width << "\n";
                file << "Target hw - " << target_half_width << "\n";
                
                file << "Sphere radius - " << sphere_radius << "\n";
                
                file.close();
                
                openvdb::FloatGrid::Ptr source_grid1 = nullptr;
                openvdb::FloatGrid::Ptr target_grid1 = nullptr;
                openvdb::FloatGrid::Ptr source_grid2 = nullptr;
                openvdb::FloatGrid::Ptr target_grid2 = nullptr;
                
                table += "<tr><td colspan='4'>dt = " + dt_arr_str[k] + ", Voxel size = " + voxel_size_str_arr[j] + ", Radius = " + std::to_string(sphere_radius) + "</td></tr>";
                
                
                
                for(int i = 0; i < sizeof(faces_arr)/sizeof(faces_arr[0]); i++) {           //faces arr
                    file.open ("log.txt", std::ios_base::app);
                    file << "Platonic scale - " << GridOperations::calcPlatonicScale(faces_arr[i], sphere_radius) << "\n\n";
                    file.close();
                    
                    source_grid1 = CommonOperations::getSphereVolume(sphere_radius, voxel_size_arr[j], source_half_width);
                    target_grid1 = CommonOperations::getPlatonicVolume(faces_arr[i], GridOperations::calcPlatonicScale(faces_arr[i], sphere_radius), voxel_size_arr[j], target_half_width);
                    
                    //                setDefaultBackgroundValue(source_grid1);
                    //                setDefaultBackgroundValue(target_grid1);
                    
                    path = "output/r1sx/dt" + dt_arr_str[k] + "/" + voxel_size_str_arr[j] + "/s" + std::to_string(faces_arr[i]) + "";
                    double total_energy1 = MorphOperations::morphModels(source_grid1, target_grid1, dt_array[k], path, table);
                    
                    table += "<td></td></tr>";
                    
                    source_grid2 = CommonOperations::getPlatonicVolume(faces_arr[i], GridOperations::calcPlatonicScale(faces_arr[i], sphere_radius), voxel_size_arr[j], source_half_width);
                    target_grid2 = CommonOperations::getSphereVolume(sphere_radius, voxel_size_arr[j], target_half_width);
                    
                    //                setDefaultBackgroundValue(source_grid2);
                    //                setDefaultBackgroundValue(target_grid2);
                    
                    path = "output/r1sx/dt" + dt_arr_str[k] + "/" + voxel_size_str_arr[j] + "/" + std::to_string(faces_arr[i]) + "s";
                    double total_energy2 = MorphOperations::morphModels(source_grid2, target_grid2, dt_array[k], path, table);
                    
                    file.open ("log.txt", std::ios_base::app);
                    file << "Mean energy - " << (total_energy1 + total_energy2) / 2 << "\n\n\n\n";
                    file.close();
                    
                    table += "<td>" + std::to_string((total_energy1 + total_energy2) / 2) + "</td></tr>";
                    
                    source_grid1 = nullptr;
                    target_grid1 = nullptr;
                    source_grid2 = nullptr;
                    target_grid2 = nullptr;
                }
                file.open ("log.txt", std::ios_base::app);
                file << table + "\n\n\n\n";
            }
        }
    }
    
    void performMorphPermutation() {
        float source_half_width = 3.0,
        target_half_width = 40.0,
        sphere_radius = 1.0;
        
        float voxel_size = 0.025;
        std::string voxel_size_str = "0.03";
        
        float dt = 0.01;
        std::string dt_str = "0.01";
        
        std::ofstream file;
        std::string table = "";
        
        int grid_face[6] = {0, 4, 6, 8, 12, 20};
        std::string labels[6] = {"Sphere", "Tetrahedron", "Cube", "Octahedron", "Dodecahedron", "Icosohedron"};
        
        std::string path = "";
        //    int grid_face[2] = {8, 12};
        //    std::string labels[2] = {"Octahedron", "Dodecahedron"};
        
        openvdb::FloatGrid::Ptr source_grid1 = nullptr;
        openvdb::FloatGrid::Ptr target_grid1 = nullptr;
        openvdb::FloatGrid::Ptr source_grid2 = nullptr;
        openvdb::FloatGrid::Ptr target_grid2 = nullptr;
        
        file.open ("log.txt");
        file.close();
        
        for(int i = 0; i < sizeof(grid_face)/sizeof(grid_face[0]); i++) {
            table += "<tr><td colspan='6'>Permuting for " + std::to_string(grid_face[i]) + " faces</td></tr>";
            
            for(int j = 0; j < sizeof(grid_face)/sizeof(grid_face[0]); j++) {
                if(j <= i) continue;
                
                file.open ("log.txt", std::ios_base::app);
                file << labels[i] << " to " << labels[j] << "\n\n";
                file.close();
                
                if(grid_face[i] == 0) {
                    source_grid1 = CommonOperations::getSphereVolume(sphere_radius, voxel_size, source_half_width);
                    target_grid2 = CommonOperations::getSphereVolume(sphere_radius, voxel_size, target_half_width);
                }
                else {
                    source_grid1 = CommonOperations::getPlatonicVolume(grid_face[i], GridOperations::calcPlatonicScale(grid_face[i], sphere_radius), voxel_size, source_half_width);
                    target_grid2 = CommonOperations::getPlatonicVolume(grid_face[i], GridOperations::calcPlatonicScale(grid_face[i], sphere_radius), voxel_size, target_half_width);
                }
                
                if(grid_face[j] == 0) {
                    source_grid2 = CommonOperations::getSphereVolume(sphere_radius, voxel_size, source_half_width);
                    target_grid1 = CommonOperations::getSphereVolume(sphere_radius, voxel_size, target_half_width);
                }
                else {
                    source_grid2 = CommonOperations::getPlatonicVolume(grid_face[j], GridOperations::calcPlatonicScale(grid_face[j], sphere_radius), voxel_size, source_half_width);
                    target_grid1 = CommonOperations::getPlatonicVolume(grid_face[j], GridOperations::calcPlatonicScale(grid_face[j], sphere_radius), voxel_size, target_half_width);
                }
                
                table += "<tr><td colspan='6'>" + std::to_string(grid_face[i]) + " to " + std::to_string(grid_face[j]) + " faces</td></tr>";
                
                path = "output/r1sx/dt" + dt_str + "/" + voxel_size_str + "/" + std::to_string(grid_face[i]) + "_" + std::to_string(grid_face[j]) + "";
                double total_energy1 = MorphOperations::morphModels(source_grid1, target_grid1, dt, path, table);
                
                table += "<td></td></tr>";
                
                path = "output/r1sx/dt" + dt_str + "/" + voxel_size_str + "/" + std::to_string(grid_face[j]) + "_" + std::to_string(grid_face[i]) + "";
                double total_energy2 = MorphOperations::morphModels(source_grid2, target_grid2, dt, path, table);
                
                table += "<td>" + std::to_string((total_energy1 + total_energy2) / 2) + "</td></tr>";
                
                source_grid1 = nullptr;
                target_grid1 = nullptr;
                source_grid2 = nullptr;
                target_grid2 = nullptr;
                
                std::cout << table << std::endl;
                
                file.open ("log.txt", std::ios_base::app);
                file << "\n" << table << "\n\n\n\n";
                file.close();
            }
        }
    }
}

#endif