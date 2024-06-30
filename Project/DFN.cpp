#include "DFN.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <set>

using namespace std;

double crossProduct2D(double x1, double y1, double x2, double y2) {
    return x1 * y2 - y1 * x2;
}

bool segmentsIntersect2D(const Vertex& p1, const Vertex& p2, const Vertex& q1, const Vertex& q2, Vertex& intersection) {
    double s1_x = p2.x - p1.x;
    double s1_y = p2.y - p1.y;
    double s2_x = q2.x - q1.x;
    double s2_y = q2.y - q1.y;

    double s = (-s1_y * (p1.x - q1.x) + s1_x * (p1.y - q1.y)) / (-s2_x * s1_y + s1_x * s2_y);
    double t = ( s2_x * (p1.y - q1.y) - s2_y * (p1.x - q1.x)) / (-s2_x * s1_y + s1_x * s2_y);

    if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
        intersection.x = p1.x + (t * s1_x);
        intersection.y = p1.y + (t * s1_y);
        intersection.z = p1.z + (t * (p2.z - p1.z));
        return true;
    }

    return false;
}

bool segmentsIntersect2D_vista_x(const Vertex& p1, const Vertex& p2, const Vertex& q1, const Vertex& q2, Vertex& intersection) {
    double s1_y = p2.y - p1.y;
    double s1_z = p2.z - p1.z;
    double s2_y = q2.y - q1.y;
    double s2_z = q2.z - q1.z;

    double s = (-s1_z * (p1.y - q1.y) + s1_y * (p1.z - q1.z)) / (-s2_y * s1_z + s1_y * s2_z);
    double t = ( s2_y * (p1.z - q1.z) - s2_z * (p1.y - q1.y)) / (-s2_y * s1_z + s1_y * s2_z);

    if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
        intersection.x = q1.x + (s * (q2.x - q1.x));
        intersection.y = p1.y + (t * s1_y);
        intersection.z = p1.z + (t * s1_z);
        return true;
    }

    return false;
}

set<Vertex> pointInFractures_x(const set<Vertex>& intersections, const Fracture& f) {
    set<Vertex> result;
    double x_min = 100, x_max = 0;
    for (size_t n = 0; n < f.vertices.size(); n++) {
        if (f.vertices[n].x > x_max) {
            x_max = f.vertices[n].x;
        }
        if (f.vertices[n].x < x_min) {
            x_min = f.vertices[n].x;
        }
    }
    for (const Vertex& point : intersections) {
        if ((point.x >= x_min && point.x <= x_max)) {
            result.insert(point);
        }
    }
    return result;
}

set<Vertex> pointInFractures_z(const set<Vertex>& intersections, const Fracture& f) {
    set<Vertex> result;
    double z_min = 100, z_max = 0;
    for (size_t n = 0; n < f.vertices.size(); n++) {
        if (f.vertices[n].z > z_max) {
            z_max = f.vertices[n].z;
        }
        if (f.vertices[n].z < z_min) {
            z_min = f.vertices[n].z;
        }
    }
    for (const Vertex& point : intersections) {
        if ((point.z >= z_min && point.z <= z_max)) {
            result.insert(point);
        }
    }
    return result;
}

vector<Vertex> calculateIntersectionPoints(const Fracture& f1, const Fracture& f2) {
    set<Vertex> uniqueIntersections;
    Vertex intersection;
    for (size_t i = 0; i < f1.vertices.size(); ++i) {
        size_t next_i = (i + 1) % f1.vertices.size();
        for (size_t j = 0; j < f2.vertices.size(); ++j) {
            size_t next_j = (j + 1) % f2.vertices.size();
            if (segmentsIntersect2D(f1.vertices[i], f1.vertices[next_i], f2.vertices[j], f2.vertices[next_j], intersection)) {
                uniqueIntersections.insert(intersection);
            }
        }
    }
    uniqueIntersections = pointInFractures_z(uniqueIntersections, f2);

    for (size_t i = 0; i < f1.vertices.size(); ++i) {
        size_t next_i = (i + 1) % f1.vertices.size();
        for (size_t j = 0; j < f2.vertices.size(); ++j) {
            size_t next_j = (j + 1) % f2.vertices.size();
            if (segmentsIntersect2D_vista_x(f1.vertices[i], f1.vertices[next_i], f2.vertices[j], f2.vertices[next_j], intersection)) {
                uniqueIntersections.insert(intersection);
            }
        }
    }
    uniqueIntersections = pointInFractures_x(uniqueIntersections, f1);

    return vector<Vertex>(uniqueIntersections.begin(), uniqueIntersections.end());
}

bool printIntersection(const Fracture& f1, const Fracture& f2) {
    vector<Vertex> intersectionPoints = calculateIntersectionPoints(f1, f2);
    if (!intersectionPoints.empty()) {
        cout << "Fracture " << f1.id << " intersects with Fracture " << f2.id << endl;
        cout << "Intersection Points:" << endl;
        for (const Vertex& point : intersectionPoints) {
            cout << "(" << point.x << ", " << point.y << ", " << point.z << ")" << endl;
        }
        return true;
    }
    cout << "Fracture " << f1.id << " does not intersect with Fracture " << f2.id << endl;
    return false;
}

bool pointOnSegment(const Vertex& point1, const Vertex& point2, const Vertex& point3) {
    double APx = point1.x - point2.x;
    double APy = point1.y - point2.y;
    double APz = point1.z - point2.z;

    double ABx = point3.x - point2.x;
    double ABy = point3.y - point2.y;
    double ABz = point3.z - point2.z;

    double crossProduct1 = APx * ABy - APy * ABx;
    double crossProduct2 = APz * ABx - APx * ABz;

    if (abs(crossProduct1) < 1e-8 && abs(crossProduct2) < 1e-8) {
        double dotProduct = APx * ABx + APy * ABy + APz * ABz;
        double AB_squared = ABx * ABx + ABy * ABy + ABz * ABz;
        double t = dotProduct / AB_squared;

        if (t >= 0 && t <= 1) {
            return true;
        }
    }

    return false;
}

bool calculateTips(const Trace& t, const Fracture& f){
    int count = 0;
    for(size_t i=0; i<t.intersectionPoints.size(); i++){
        for(size_t j=0; j<f.vertices.size(); j++){
            Vertex point1 = t.intersectionPoints[i];
            Vertex point2 = f.vertices[j];
            Vertex point3 = f.vertices[(j+1)% f.vertices.size()];
            if(pointOnSegment(point1, point2, point3)){
                count++;
            }
        }
    }
    if(count>=2){
        return true;
    }
    return false;
}

vector<Trace> calculateTraces(const vector<Fracture>& fractures) {
    int numFractures = fractures.size();
    vector<Trace> Traces;
    int ids = 0;
    for (int i = 0; i < numFractures - 1; i++) {
        for (int j = i + 1; j < numFractures; j++) {
            if (printIntersection(fractures[i], fractures[j])) {
                Trace Intersection;
                Intersection.id = ids;
                ids++;
                Intersection.fractureId1 = fractures[i].id;
                Intersection.fractureId2 = fractures[j].id;
                Intersection.intersectionPoints = calculateIntersectionPoints(fractures[i], fractures[j]);
                Intersection.Tips1 = calculateTips(Intersection, fractures[i]);
                Intersection.Tips2 = calculateTips(Intersection, fractures[j]);
                Traces.push_back(Intersection);
            }
        }
    }
    return Traces;
}

vector<Fracture> readDFNFile(const string& filename) {
    ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        cerr << "Error: Unable to open input file." << endl;
        return {};
    }

    string line;
    int numFractures = 0;

    while (getline(inputFile, line)) {
        if (line.empty() || line[0] == '#')
            continue;

        stringstream ss(line);
        ss >> numFractures;
        if (ss.fail()) {
            cerr << "Error: Unable to read number of fractures from file." << endl;
            return {};
        }
        break;
    }

    vector<Fracture> fractures;
    while (getline(inputFile, line)) {
        if (line.empty() || line[0] == '#')
            continue;

        stringstream ss(line);
        int id, numVertices;
        char delimiter;
        ss >> id >> delimiter >> numVertices;
        if (ss.fail() || delimiter != ';') {
            cerr << "Error: Unable to read fracture data from file." << endl;
            return {};
        }

        Fracture fracture;
        fracture.id = id;
        fracture.vertices.resize(numVertices);

        vector<double> xCoords(numVertices), yCoords(numVertices), zCoords(numVertices);

        while (getline(inputFile, line)) {
            if (line.empty() || line[0] == '#')
                continue;

            stringstream vertexStream(line);
            for (int i = 0; i < numVertices; ++i) {
                char dummy;
                vertexStream >> xCoords[i] >> dummy;
                if (dummy != ';' && i != numVertices - 1) {
                    cerr << "Error: Expected ';' delimiter in x coordinates." << endl;
                    return {};
                }
            }
            break;
        }

        while (getline(inputFile, line)) {
            if (line.empty() || line[0] == '#')
                continue;

            stringstream vertexStream(line);
            for (int i = 0; i < numVertices; ++i) {
                char dummy;
                vertexStream >> yCoords[i] >> dummy;
                if (dummy != ';' && i != numVertices - 1) {
                    cerr << "Error: Expected ';' delimiter in y coordinates." << endl;
                    return {};
                }
            }
            break;
        }

        while (getline(inputFile, line)) {
            if (line.empty() || line[0] == '#')
                continue;

            stringstream vertexStream(line);
            for (int i = 0; i < numVertices; ++i) {
                char dummy;
                vertexStream >> zCoords[i] >> dummy;
                if (dummy != ';' && i != numVertices - 1) {
                    cerr << "Error: Expected ';' delimiter in z coordinates." << endl;
                    return {};
                }
            }
            break;
        }

        for (int i = 0; i < numVertices; ++i) {
            fracture.vertices[i] = {xCoords[i], yCoords[i], zCoords[i]};
        }

        fractures.push_back(fracture);
    }
    inputFile.close();
    return fractures;
}
