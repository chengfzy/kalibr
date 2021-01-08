#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <iostream>
#include "bsplines/BSpline.hpp"
#include "cc/Heading.hpp"

using namespace std;
using namespace Eigen;
using namespace fmt;
using namespace cc;
using namespace bsplines;

int main(int argc, char* argv[]) {
    google::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);

    cout << Section("", false) << Section("B-Spline Tets", false) << Section("", false);

    // define spline of degree = 3(order = 4)
    BSpline spline(4);
    // set knot and control points
    vector<double> knots = {0, 0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5, 5};
    MatrixXd ctrlPoints(1, 9);
    ctrlPoints << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9;
    cout << format("knots: {}", knots) << endl;
    cout << format("control points: {}", ctrlPoints) << endl;
    spline.setKnotsAndCoefficients(knots, ctrlPoints);

    // the number of variables
    cout << Section("Number of Variables");
    cout << format("numCoefficientsRequired = {}", spline.numCoefficientsRequired(8)) << endl;
    cout << format("numKnotsRequired = {}", spline.numKnotsRequired(8)) << endl;
    cout << format("minimumKnotsRequired = {}", spline.minimumKnotsRequired()) << endl;
    cout << format("numValidTimeSegments(12) = {}", spline.numValidTimeSegments(12)) << endl;
    cout << format("numValidTimeSegments = {}", spline.numValidTimeSegments()) << endl;
    cout << format("tMin = {}", spline.tMin()) << endl;
    cout << format("tMax = {}", spline.tMax()) << endl;

    // evaluate
    cout << Section("Evaluation");
    vector<double> inputU{0, 1, 2, 3, 4, 4.5, 5};
    for (auto v : inputU) {
        MatrixXd j1, j2;
        auto y0 = spline.eval(v);
        auto y1 = spline.evalDAndJacobian(v, 1, &j1, nullptr);
        auto y2 = spline.evalDAndJacobian(v, 2, &j2, nullptr);
        cout << format("f({0}) = {1}, f'({0}) = {2}, f''({0}) = {3}", v, y0[0], y1[0], y2[0]) << endl;
        cout << format("\tJ1 = {}, J2 = {}", j1, j2) << endl;
    }

    // segment quadratic integral
    cout << Section("Segment Quadratic Integral");
    VectorXd WDiag(1);
    WDiag[0] = 1;
    MatrixXd int1 = spline.segmentQuadraticIntegralDiag(WDiag, 1, 1);
    cout << format("integral 1 = \n{}", int1) << endl;
    MatrixXd int2 = spline.segmentQuadraticIntegralDiag(WDiag, 1, 2);
    cout << format("integral 2 = \n{}", int2) << endl;

    google::ShutdownGoogleLogging();
    google::ShutDownCommandLineFlags();
    return 0;
}