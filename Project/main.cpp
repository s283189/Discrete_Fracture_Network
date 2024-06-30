#include "DFN.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <gtest/gtest.h>

using namespace std;

// GOOGLE TEST

// Test per crossProduct2D
TEST(CrossProduct2DTest, BasicAssertions) {
    EXPECT_DOUBLE_EQ(crossProduct2D(1, 2, 3, 4), -2);
    EXPECT_DOUBLE_EQ(crossProduct2D(0, 0, 0, 0), 0);
    EXPECT_DOUBLE_EQ(crossProduct2D(1, 0, 0, 1), 1);
}

// Test per segmentsIntersect2D
TEST(SegmentsIntersect2DTest, IntersectingSegments) {
    Vertex p1 = {0, 0, 0}, p2 = {2, 2, 0};
    Vertex q1 = {0, 2, 0}, q2 = {2, 0, 0};
    Vertex intersection;

    EXPECT_TRUE(segmentsIntersect2D(p1, p2, q1, q2, intersection));
    EXPECT_DOUBLE_EQ(intersection.x, 1);
    EXPECT_DOUBLE_EQ(intersection.y, 1);
}

TEST(SegmentsIntersect2DTest, NonIntersectingSegments) {
    Vertex p1 = {0, 0, 0}, p2 = {1, 1, 0};
    Vertex q1 = {2, 2, 0}, q2 = {3, 3, 0};
    Vertex intersection;

    EXPECT_FALSE(segmentsIntersect2D(p1, p2, q1, q2, intersection));
}

// Test per pointOnSegment
TEST(PointOnSegmentTest, PointOnSegment) {
    Vertex point1 = {1, 1, 1};
    Vertex point2 = {0, 0, 0};
    Vertex point3 = {2, 2, 2};

    EXPECT_TRUE(pointOnSegment(point1, point2, point3));
}

TEST(PointOnSegmentTest, PointNotOnSegment) {
    Vertex point1 = {3, 3, 3};
    Vertex point2 = {0, 0, 0};
    Vertex point3 = {2, 2, 2};

    EXPECT_FALSE(pointOnSegment(point1, point2, point3));
}

// Test per calculateIntersectionPoints
TEST(CalculateIntersectionPointsTest, BasicIntersections) {
    Fracture f1 = {1, {{0, 0, 0}, {2, 2, 0}}};
    Fracture f2 = {2, {{0, 2, 0}, {2, 0, 0}}};
    auto intersections = calculateIntersectionPoints(f1, f2);

    ASSERT_EQ(intersections.size(), 1);
    EXPECT_DOUBLE_EQ(intersections[0].x, 1);
    EXPECT_DOUBLE_EQ(intersections[0].y, 1);
}

TEST(CalculateIntersectionPointsTest, NoIntersections) {
    Fracture f1 = {1, {{0, 0, 0}, {1, 1, 0}}};
    Fracture f2 = {2, {{2, 2, 0}, {3, 3, 0}}};
    auto intersections = calculateIntersectionPoints(f1, f2);

    EXPECT_EQ(intersections.size(), 0);
}

int main(int argc, char **argv) {
    string filenameInput = "C:\\Desktop\\PDF UNIVERSITA'\\PROGRAMMAZIONE E CALCOLO SCIENTIFICO\\PROGETTO\\Progetto_PCS_2024\\Project\\DFN\\FR3_data.txt";
    ofstream outputFile("output.txt");

    vector<Fracture> fractures = readDFNFile(filenameInput);
    vector<Trace> Traces;

    if (fractures.empty()) {
        cerr << "No fractures read from the file." << endl;
        return 1;
    }

    Traces = calculateTraces(fractures);
    outputFile << scientific << setprecision(16);

    // Scrivi il numero di tracce nel file
    outputFile << "# Number of Traces" << endl;
    outputFile << Traces.size() << endl;

    // Scrivi i dettagli delle tracce nel file
    outputFile << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;
    for (const Trace& trace : Traces) {
        for (size_t i = 0; i < trace.intersectionPoints.size(); i += 2) {
            const Vertex& p1 = trace.intersectionPoints[i];
            const Vertex& p2 = trace.intersectionPoints[i + 1];
            outputFile << trace.id << "; " << trace.fractureId1 << "; " << trace.fractureId2 << "; "
                       << p1.x << "; " << p1.y << "; " << p1.z << "; "
                       << p2.x << "; " << p2.y << "; " << p2.z << endl;
        }
    }

    outputFile << "#Passanti" << endl;
    outputFile << "#TraceId; FractureId" << endl;

    for (const Trace& trace : Traces) {
        vector<int> passingFractureIds; // per memorizzare gli ID delle fratture per cui la traccia è passante

        for (const Fracture& fracture : fractures) {
            if (calculateTips(trace, fracture)) {
                passingFractureIds.push_back(fracture.id);
            }
        }

        // Stampa l'ID della traccia e gli ID delle fratture per cui è passante
        if (!passingFractureIds.empty()) {
            for (int fractureId : passingFractureIds) {
                outputFile << trace.id << "; " << fractureId << endl;
            }
        }
    }

    // Ora affrontiamo le tracce non passanti
    outputFile << "#Non Passanti" << endl;
    outputFile << "#TraceId" << endl;

    for (const Trace& trace : Traces) {
        bool isPassing = false;

        for (const Fracture& fracture : fractures) {
            if (calculateTips(trace, fracture)) {
                isPassing = true;
                break;
            }
        }

        // Se la traccia non è passante per nessuna frattura, la stampiamo
        if (!isPassing) {
            // Stampiamo solo l'ID della traccia senza ID frattura
            outputFile << trace.id << endl;
        }
    }
        outputFile.close();

    // GOOGLE TEST
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();

    return 0;
}

