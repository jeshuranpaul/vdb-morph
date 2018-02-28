//
//  Container.cpp
//  vdbTest
//
//  Created by Hanjie Liu on 7/28/17.
//  Copyright © 2017 Hanjie Liu. All rights reserved.
//
// code partially from http://www.cplusplus.com/doc/tutorial/files/
// and https://stackoverflow.com/questions/236129/most-elegant-way-to-split-a-string

#include "Container.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <openvdb/openvdb.h>
#include <openvdb/points/PointConversion.h>
#include <openvdb/points/PointCount.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/MeshToVolume.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <tuple>
#include <algorithm>
#include <iterator>

#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>
#include <boost/graph/eccentricity.hpp>

namespace PFIAir {
    using namespace openvdb;
    using namespace std;
    
    Container::Container() {
        initialize();
    }
    
    void Container::setScale(Vec3d scale) {
        this -> _scale.preScale(scale);
    }
    
    /// import index and face defined model files like obj or smf file
    void Container::loadMeshModel(const string filename) {
        this -> _filename = filename;
        
        string line;
        ifstream myfile (filename);
        
        if (myfile.is_open())
        {
            while ( getline (myfile,line) )
            {
                
                istringstream iss(line);
                vector<string> tokens{istream_iterator<string>{iss},
                    istream_iterator<string>{}};
                
                if (tokens.size() > 0) {
                    // vertices
                    if (!tokens[0].compare("v")) {
                        assert(tokens.size() == 4);

                        this -> _points.push_back(Vec3s(stof(tokens[1]),stof(tokens[2]),stof(tokens[3])));
                        continue;
                    }
                    
                    // facets
                    if (!tokens[0].compare("f")) {
                        assert(tokens.size() == 4 || tokens.size() == 5 );
                        
                        if (tokens.size() == 4) {
                            // fix obj indexing issue by -1
                            _indicesTri.push_back(Vec3I(stoi(tokens[1]) - 1,stoi(tokens[2]) - 1,stoi(tokens[3]) - 1));
                        } else {
                            _indicesQuad.push_back(Vec4I(stoi(tokens[1]) - 1,stoi(tokens[2]) - 1,stoi(tokens[3]) - 1,stoi(tokens[4]) - 1));
                        }
                    }
                }
            }
            myfile.close();
        }
        
        else cout << "Unable to open file" << endl;
    }
    
    void Container::OLDcomputeMeshCenter() {
        using namespace boost;

        typedef std::pair<int, int> Edge;
        std::vector<Edge> edges;

        for (int i = 0; i < _indicesTri.size(); i++) {
            int point1 = _indicesTri[i].x();
            int point2 = _indicesTri[i].y();
            int point3 = _indicesTri[i].z();

            Edge edge1 = Edge(point1, point2);
            Edge edge2 = Edge(point1, point3);
            Edge edge3 = Edge(point2, point3);

            edges.push_back(edge1);
            edges.push_back(edge2);
            edges.push_back(edge3);
        }

        for (int i = 0; i < _indicesQuad.size(); i++) {
            int point1 = _indicesQuad[i].x();
            int point2 = _indicesQuad[i].y();
            int point3 = _indicesQuad[i].z();
            int point4 = _indicesQuad[i].w();

            Edge edge1 = Edge(point1, point2);
            Edge edge2 = Edge(point2, point3);
            Edge edge3 = Edge(point3, point4);
            Edge edge4 = Edge(point4, point1);

            edges.push_back(edge1);
            edges.push_back(edge2);
            edges.push_back(edge3);
            edges.push_back(edge4);
        }

        typedef adjacency_list<vecS, vecS, undirectedS
        > Graph;

        Graph g(edges.begin(), edges.end(), edges.size());


                shared_array_property_map<double, property_map<Graph, vertex_index_t>::const_type>
                centrality_map(num_vertices(g), get(boost::vertex_index, g));
        
                brandes_betweenness_centrality(g, centrality_map);


       // cout << num_vertices(g) << endl;

                double max = 0;
                int maxIndex = 0;
                for (int i = 0; i < num_vertices(g); i++) {
                    if (centrality_map[i] > max) {
                        max = centrality_map[i];
                        maxIndex = i;
                    }
                }
        
    }
    
    typedef std::vector<std::pair<double,Vec3s>> Result;
    
    bool comp (std::pair<double,Vec3s> left, std::pair<double,Vec3s> right) {
        return left.first > right.first;
    }
    void insert( Result &cont, std::pair<double,Vec3s> value ) {
        Result::iterator it = std::lower_bound( cont.begin(), cont.end(), value, comp);
        cont.insert( it, value );
    }
    
    bool edgeDup(std::vector<std::pair<int, int>> edges, int a, int b) {
        for (int i = 0; i < edges.size(); i++) {
            if ((edges[i].first == a && edges[i].second == b) || (edges[i].first == b && edges[i].second == a)) {
                return false;
            }
        }
        return true;
    }
    
    double Container::getEuclideanFromIndices(int p1, int p2){
        using namespace openvdb;
        Vec3s point1 = _points[p1];
        Vec3s point2 = _points[p2];
        
        Vec3s diff = point1.sub(point2, point1);
        
        return (double)diff.length();
    }
    
    void Container::computeMeshCenter() {
        using namespace boost;

        typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
//        typedef adjacency_list<vecS, vecS, undirectedS,boost::no_property, EdgeWeightProperty
//        > Graph;
        typedef adjacency_list<vecS, vecS, undirectedS, boost::no_property, EdgeWeightProperty> Graph;
        typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t;
        typedef boost::graph_traits<Graph>::edge_descriptor edge_t;

        std::map<int, vertex_t> m;
        std::vector<int> ver;
        
        Graph g;
        std::vector<std::pair<int, int>> edges;
        
        for (int i = 0; i < _indicesQuad.size(); i++) {
            int point1 = _indicesQuad[i].x();
            int point2 = _indicesQuad[i].y();
            int point3 = _indicesQuad[i].z();
            int point4 = _indicesQuad[i].w();
            
            if (m.find(point1) == m.end()) {
                m[point1] = (vertex_t)add_vertex(g);
                ver.push_back(point1);
            }
            if (m.find(point2) == m.end()) {
                m[point2] = (vertex_t)add_vertex(g);
                ver.push_back(point2);
            }
            if (m.find(point3) == m.end()) {
                m[point3] = (vertex_t)add_vertex(g);
                ver.push_back(point3);
            }
            
            if (m.find(point4) == m.end()) {
                m[point4] = (vertex_t)add_vertex(g);
                ver.push_back(point4);
            }
            
            vertex_t d1 = m.find(point1)->second;
            vertex_t d2 = m.find(point2)->second;
            vertex_t d3 = m.find(point3)->second;
            vertex_t d4 = m.find(point4)->second;
            
            if (edgeDup(edges, point1, point2)) {
                EdgeWeightProperty e = getEuclideanFromIndices(point1, point2);
                add_edge(d1, d2, e, g);
                edges.push_back(std::pair<int, int>(point1, point2));
            }
            if (edgeDup(edges, point2, point3)) {
                EdgeWeightProperty e = getEuclideanFromIndices(point2, point3);
                add_edge(d2, d3, e, g);
                edges.push_back(std::pair<int, int>(point2, point3));
            }
            if (edgeDup(edges, point3, point4)) {
                EdgeWeightProperty e = getEuclideanFromIndices(point3, point4);
                add_edge(d3, d4, e, g);
                edges.push_back(std::pair<int, int>(point3, point4));
            }
            if (edgeDup(edges, point4, point1)) {
                EdgeWeightProperty e = getEuclideanFromIndices(point4, point1);
                add_edge(d4, d1, e, g);
                edges.push_back(std::pair<int, int>(point4, point1));
            }
        }
        
        for (int i = 0; i < _indicesTri.size(); i++) {
            int point1 = _indicesTri[i].x();
            int point2 = _indicesTri[i].y();
            int point3 = _indicesTri[i].z();
            
            if (m.find(point1) == m.end()) {
                m[point1] = (vertex_t)add_vertex(g);
                ver.push_back(point1);
            }
            if (m.find(point2) == m.end()) {
                m[point2] = (vertex_t)add_vertex(g);
                ver.push_back(point2);
            }
            if (m.find(point3) == m.end()) {
                m[point3] = (vertex_t)add_vertex(g);
                ver.push_back(point3);
            }

            vertex_t d1 = m.find(point1)->second;
            vertex_t d2 = m.find(point2)->second;
            vertex_t d3 = m.find(point3)->second;

            //add_edge(d1,d2,EdgeWeightProperty(2),g);
            if (edgeDup(edges, point1, point2)) {
                EdgeWeightProperty e = getEuclideanFromIndices(point1, point2);
                add_edge(d1,d2,e,g);
                edges.push_back(std::pair<int, int>(point1, point2));
            }
            if (edgeDup(edges, point1, point3)) {
                EdgeWeightProperty e = getEuclideanFromIndices(point1, point3);
                add_edge(d1,d3,e,g);
                edges.push_back(std::pair<int, int>(point1, point3));
            }
            if (edgeDup(edges, point2, point3)) {
                EdgeWeightProperty e = getEuclideanFromIndices(point2, point3);
                add_edge(d2,d3,e,g);
                edges.push_back(std::pair<int, int>(point2, point3));
            }
        }
        
//        shared_array_property_map<double, property_map<Graph, vertex_index_t>::const_type>
//        centrality_map(num_vertices(g), get(boost::vertex_index, g));
//        
        property_map<Graph, edge_weight_t>::type w = get(edge_weight, g);
//
//        brandes_betweenness_centrality(g, boost::centrality_map(centrality_map).weight_map(w));
        
        typedef graph_traits<Graph>::edge_descriptor Edge;
        typedef exterior_vertex_property<Graph, double> DistanceProperty;
        typedef DistanceProperty::matrix_type DistanceMatrix;
        typedef DistanceProperty::matrix_map_type DistanceMatrixMap;
        typedef constant_property_map<Edge, double> ECCWeightMap;
        
        DistanceMatrix distances(num_vertices(g));
        DistanceMatrixMap dm(distances, g);
        ECCWeightMap wm(1);
        //floyd_warshall_all_pairs_shortest_paths(g, dm, w);
        floyd_warshall_all_pairs_shortest_paths(g, dm, weight_map(w));
        
        typedef boost::exterior_vertex_property<Graph, double> EccentricityProperty;
        typedef EccentricityProperty::container_type EccentricityContainer;
        typedef EccentricityProperty::map_type EccentricityMap;
        
        int r, d;
        EccentricityContainer eccs(num_vertices(g));
        EccentricityMap em(eccs, g);
        boost::tie(r, d) = all_eccentricities(g, dm, em);
        
        Result sorted_result;
        map<int, vertex_t>::iterator it;
        int c = 0;
        for ( it = m.begin(); it != m.end(); it++ ) {
            
            insert(sorted_result, std::pair<double, Vec3s>(em[c], _points[ver.at(c)]));
            
            c++;
        }
        
        return;
        

//        Result sorted_result;
//        map<int, vertex_t>::iterator it;
//        int c = 0;
//        for ( it = m.begin(); it != m.end(); it++ ) {
//
//            insert(sorted_result, std::pair<double, Vec3s>(centrality_map[c], _points[ver.at(c)]));
//
//            c++;
//        }
//        
//        for (int i = 0; i < sorted_result.size(); i++) {
//            cout << sorted_result[i].first << " " << sorted_result[i].second << endl;
//        }
//        
//        cout << "\nvertices: " << num_vertices(g) << endl;
//        cout << "edges: " << num_edges(g) << endl;
    }
    
    FloatGrid::Ptr Container::getWaterTightLevelSet() {
        return tools::meshToLevelSet<FloatGrid>(_scale, _points, _indicesTri, _indicesQuad);
    }
    
    FloatGrid::Ptr Container::getWaterTightLevelSetWithBandWidth(float w) {
        return tools::meshToLevelSet<FloatGrid>(_scale, _points, _indicesTri, _indicesQuad, w);
    }
    
    FloatGrid::Ptr Container::getUnsignedDistanceField(float bandwidth) {
        return tools::meshToUnsignedDistanceField<FloatGrid>(_scale, _points, _indicesTri, _indicesQuad, bandwidth);
    }
    
    void Container::exportModel(const string name, FloatGrid::Ptr model) {
        io::File file(name);
        GridPtrVec grids;
        grids.push_back(model);
        file.write(grids);
        file.close();
    }
    
    float Container::computeAverageEdgeLength() {
        vector<float> tri_avg = vector<float>();
        vector<float> quad_avg = vector<float>();

        // average length for triangles
        for (int i = 0; i < _indicesTri.size(); i++) {
            Vec3s point1 = _points[_indicesTri[i].x()];
            Vec3s point2 = _points[_indicesTri[i].y()];
            Vec3s point3 = _points[_indicesTri[i].z()];

            Vec3s edge1 = point1 - point2;
            Vec3s edge2 = point1 - point3;
            Vec3s edge3 = point2 - point3;
            
            float avg = (edge1.length() + edge2.length() + edge3.length()) / 3.0;
            tri_avg.push_back(avg);
        }
        
        // average length for quads
        for (int i = 0; i < _indicesQuad.size(); i++) {
            Vec3s point1 = _points[_indicesQuad[i].x()];
            Vec3s point2 = _points[_indicesQuad[i].y()];
            Vec3s point3 = _points[_indicesQuad[i].z()];
            Vec3s point4 = _points[_indicesQuad[i].w()];
            
            Vec3s edge1 = point1 - point2;
            Vec3s edge2 = point2 - point3;
            Vec3s edge3 = point3 - point4;
            Vec3s edge4 = point4 - point1;

            
            float avg = (edge1.length() + edge2.length() + edge3.length() + edge4.length()) / 4.0;
            quad_avg.push_back(avg);
        }
        
        // compute overall average
        float tri_sum = 0, quad_sum = 0;
        for (int i = 0; i < tri_avg.size(); i++) {
            tri_sum += tri_avg[i];
        }
        
        for (int i = 0; i < quad_avg.size(); i++) {
            quad_sum += quad_avg[i];
        }
        
        float weighted = (tri_sum + quad_sum) / (tri_avg.size() + quad_avg.size());
        
        return weighted;
    }
}
