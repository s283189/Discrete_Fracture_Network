#pragma once

#include <vector>
#include <string>
#include <set>

using namespace std;

struct Vertex {
    double x, y, z;

    bool operator<(const Vertex& other) const {
        if (x != other.x) return x < other.x;
        if (y != other.y) return y < other.y;
        return z < other.z;
    }
};

struct Fracture {
    int id;
    vector<Vertex> vertices;
};

struct Trace {
    int id;
    int fractureId1;
    int fractureId2;
    bool Tips1;
    bool Tips2;
    vector<Vertex> intersectionPoints;
};

bool segmentsIntersect2D(const Vertex& p1, const Vertex& p2, const Vertex& q1, const Vertex& q2, Vertex& intersection);
bool segmentsIntersect2D_vista_x(const Vertex& p1, const Vertex& p2, const Vertex& q1, const Vertex& q2, Vertex& intersection);
set<Vertex> pointInFractures_x(const set<Vertex>& intersections, const Fracture& f);
set<Vertex> pointInFractures_z(const set<Vertex>& intersections, const Fracture& f);
vector<Vertex> calculateIntersectionPoints(const Fracture& f1, const Fracture& f2);
bool printIntersection(const Fracture& f1, const Fracture& f2);
vector<Trace> calculateTraces(const vector<Fracture>& fractures);
bool pointOnSegment(const Vertex& point1, const Vertex& point2, const Vertex& point3);
bool calculateTips(const Trace& t, const Fracture& f);
vector<Fracture> readDFNFile(const string& filename);


