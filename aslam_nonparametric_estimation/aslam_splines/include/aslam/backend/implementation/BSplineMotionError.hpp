#include <stdio.h>
#include <aslam/backend/BSplineMotionError.hpp>
#include "sm/DebugInfo.h"

using namespace std;

namespace aslam {
namespace backend {

template <class SPLINE_T>
BSplineMotionError<SPLINE_T>::BSplineMotionError(spline_t* splineDV, Eigen::MatrixXd W) : _splineDV(splineDV), _W(W) {
    initialize(splineDV, W, 2);  // default: acceleration criterion
}

template <class SPLINE_T>
BSplineMotionError<SPLINE_T>::BSplineMotionError(spline_t* splineDV, Eigen::MatrixXd W, unsigned int errorTermOrder)
    : _splineDV(splineDV), _W(W) {
    initialize(splineDV, W, errorTermOrder);
}

template <class SPLINE_T>
BSplineMotionError<SPLINE_T>::~BSplineMotionError() {}

template <class SPLINE_T>
void BSplineMotionError<SPLINE_T>::initialize(spline_t* /* splineDV */, Eigen::MatrixXd /* W */,
                                              unsigned int errorTermOrder) {
    // check spline order:
    int splineOrder = _splineDV->spline().splineOrder();
    if (splineOrder <= (int)errorTermOrder && splineOrder >= 2) {
        errorTermOrder = splineOrder - 1;
        std::cout << "! Invalid ErrorTermOrder reduced to " << errorTermOrder << std::endl;
    }
    sbm_t Qsp = _splineDV->spline().curveQuadraticIntegralSparse(_W, errorTermOrder);

    // set spline design variables
    unsigned int numberOfSplineDesignVariables = _splineDV->spline().numVvCoefficients();
    _coefficientVectorLength = _splineDV->spline().coefficients().rows() * numberOfSplineDesignVariables;

    Qsp.cloneInto(_Q);
    //      std::cout << "next:" << std::endl;

    //      std::cout << Q.cols() << ":" << Q.rows() << std::endl;
    //      std::cout << _Q->cols() << ":" << _Q->rows() << std::endl;

    // Tell the super class about the design variables:
    // loop the design variables and add to vector:
    std::vector<aslam::backend::DesignVariable*> dvV;
    for (unsigned int i = 0; i < numberOfSplineDesignVariables; i++) {
        dvV.push_back(_splineDV->designVariable(i));
    }
    setDesignVariables(dvV);

    // print debug info
    printInfo("BSplineMotionError<SPLINE_T>::initialize()");
    cout << "errorTermOrder = " << errorTermOrder << endl;
    cout << "numberOfSplineDesignVariables = " << numberOfSplineDesignVariables << endl;
    cout << "numDesignVariables() = " << numDesignVariables() << endl;
    cout << "_coefficientVectorLength = " << _coefficientVectorLength << endl;
    cout << "W size = " << _W.rows() << " X " << _W.cols() << ", |W| = " << _W.norm() << endl;
    cout << "Q size = " << _Q.rows() << " X " << _Q.cols() << endl;
    cout << "|sqrtInvR()| = " << sqrtInvR().norm() << endl;
}

// NOTE by CC. There are 2 main propose in evaluateErrorImplementation(), (1) calculate the squared weight error of
// object function, e^T * W * e; (2) calculate the error of e and call ErrorTermFs<T>::setError(const
// Eigen::MatrixBase<DERIVED>& e) in parent class. But in the BSplineMotionError class, the (2) task is hard to set, and
// then the incremental optimization method based on Gaussian-Newton cannot finished. That's is why the
// evaluateJacobiansImplementation() function isn't implemented.

/// \brief evaluate the error term and return the weighted squared error e^T invR e
template <class SPLINE_T>
double BSplineMotionError<SPLINE_T>::evaluateErrorImplementation() {
    // the error is a scalar: c' Q c, with c the vector valued spline coefficients stacked

    const double* cMat = &((_splineDV->spline()).coefficients()(0, 0));
    Eigen::Map<const Eigen::VectorXd> c = Eigen::VectorXd::Map(cMat, _coefficientVectorLength);

    // Q*c :
    // create result container:
    Eigen::VectorXd Qc(_Q.rows());  // number of rows of Q:
    Qc.setZero();
    //  std::cout << Qc->rows() << ":" << Qc->cols() << std::endl;
    _Q.multiply(&Qc, c);

    return c.transpose() * (Qc);
}

/// \brief evaluate the jacobians
template <class SPLINE_T>
void BSplineMotionError<SPLINE_T>::evaluateJacobiansImplementation(
    aslam::backend::JacobianContainer& _jacobians) const {
    double errorImpl = sqrt((const_cast<BSplineMotionError<SPLINE_T>*>(this))->evaluateErrorImplementation());

#if 0
    // jacobian = Q * c + Q * dc ~= Q * c ?
    const double* cMat = &((_splineDV->spline()).coefficients()(0, 0));
    Eigen::Map<const Eigen::VectorXd> c = Eigen::VectorXd::Map(cMat, _coefficientVectorLength);
    Eigen::VectorXd jac(_Q.rows());
    jac.setZero();
    _Q.multiply(&jac, c);
    jac *= 2;

    printInfo("BSplineMotionError<SPLINE_T>::evaluateJacobiansImplementation()", DebugInfoType::Paragraph, true);
    cout << "numDesignVariables() = " << numDesignVariables()
         << ", v.minimalDimensions() = " << designVariable(0)->minimalDimensions() << endl;
    cout << "c size = " << c.size() << ", |c| = " << c.norm() << endl;
    cout << "Q size = " << _Q.rows() << " X " << _Q.cols() << ", |Q| = " << _Q.toDense().norm() << endl;
    cout << "Jac size = " << jac.size() << ", |Jac| = " << jac.norm() << endl;
    cout << "|sqrtInvR()| = " << sqrtInvR().norm() << ", |error| = " << error().norm() << ", errorImpl = " << errorImpl
         << endl;

    // resize output jacobian
    // _jacobians.reset(_splineDV->spline().coefficients().rows());

    // assign data
    for (size_t i = 0; i < numDesignVariables(); ++i) {
        const DesignVariable* v = designVariable(i);
        // cout << "[" << i << "], columnBase = " << v->columnBase() << ", blockIndex = " << v->blockIndex()
        //      << ", minDim = " << v->minimalDimensions() << endl;
        _jacobians.add((const_cast<BSplineMotionError<SPLINE_T>*>(this))->designVariable(i),
                       jac.segment(i * v->minimalDimensions(), v->minimalDimensions()).transpose());
    }
#else
    // -Q * C = - Jac * f(x), f(x)^2 = J(x) = 0.5 * c' * Q * c is scalar; => Jac = Q * C / sqrt(0.5 * c' * Q * c)
    const double* cMat = &((_splineDV->spline()).coefficients()(0, 0));
    Eigen::Map<const Eigen::VectorXd> c = Eigen::VectorXd::Map(cMat, _coefficientVectorLength);

    // calculate fx
    Eigen::VectorXd Qc(_Q.rows());
    Qc.setZero();
    _Q.multiply(&Qc, c);
    // double fx2 = 0.5 * c.transpose() * Qc; // f(x)^2 = 0.5 * c' * Q * c
    double fx2 = c.transpose() * Qc;  // f(x)^2 = c' * Q * c

    // calculate jac = 2 * Q * c / f(x)
    Eigen::VectorXd jac(_Q.rows());
    jac.setZero();
    // jac = Qc / 2.0 / sqrt(fx2);
    jac = Qc / sqrt(fx2);  // jac = Qc / fx

    printInfo("BSplineMotionError<SPLINE_T>::evaluateJacobiansImplementation()", DebugInfoType::Paragraph, true);
    cout << "numDesignVariables() = " << numDesignVariables()
         << ", v.minimalDimensions() = " << designVariable(0)->minimalDimensions() << endl;
    cout << "Q size = " << _Q.rows() << " X " << _Q.cols() << ", |Q| = " << _Q.toDense().norm() << endl;
    cout << "Qc size = " << Qc.rows() << " X " << Qc.cols() << ", |Qc| = " << Qc.norm() << endl;
    cout << "jac size = " << jac.rows() << " X " << jac.cols() << ", |jac| = " << jac.norm() << endl;
    cout << "|sqrtInvR()| = " << sqrtInvR().norm() << ", |error| = " << error().norm() << ", errorImpl = " << errorImpl
         << endl;

    // assign data
    for (size_t i = 0; i < numDesignVariables(); ++i) {
        const DesignVariable* v = designVariable(i);
        // cout << "[" << i << "], columnBase = " << v->columnBase() << ", blockIndex = " << v->blockIndex()
        //      << ", minDim = " << v->minimalDimensions() << endl;
        _jacobians.add((const_cast<BSplineMotionError<SPLINE_T>*>(this))->designVariable(i),
                       jac.segment(i * v->minimalDimensions(), v->minimalDimensions()).transpose());
    }
#endif

    // this is an error...
    // SM_THROW(Exception, "This is currently unsupported");
}

template <class SPLINE_T>
void BSplineMotionError<SPLINE_T>::buildHessianImplementation(SparseBlockMatrix& outHessian, Eigen::VectorXd& outRhs,
                                                              bool /* useMEstimator */) {
    // get the coefficients:
    Eigen::MatrixXd coeff = _splineDV->spline().coefficients();
    // create a column vector of spline coefficients
    int dim = coeff.rows();
    int seg = coeff.cols();
    // build a vector of coefficients:
    Eigen::VectorXd c(dim * seg);
    // rows are spline dimension
    for (int i = 0; i < seg; i++) {
        c.block(i * dim, 0, dim, 1) = coeff.block(0, i, dim, 1);
    }

    // right hand side:
    Eigen::VectorXd b_u(_Q.rows());  // number of rows of Q:

    b_u.setZero();
    /*      std::cout <<"b" << std::endl;
          for(int i = 0 ; i < b_u->rows(); i++)
              std::cout << (*b_u)(i) << std::endl;
                      std::cout <<"/b" << std::endl;  */

    _Q.multiply(&b_u, c);

#if false
    printInfo("BSplineMotionError<SPLINE_T>::buildHessianImplementation()", DebugInfoType::Paragraph, true);
    double errorImpl = sqrt((const_cast<BSplineMotionError<SPLINE_T>*>(this))->evaluateErrorImplementation());
    cout << "c size = " << c.size() << ", |c| = " << c.norm() << endl;
    cout << "Q size = " << _Q.rows() << " X " << _Q.cols() << ", |Q| = " << _Q.toDense().norm() << endl;
    cout << "bu size = " << b_u.size() << ", |b_u| = " << b_u.norm() << endl;
    cout << "|sqrtInvR()| = " << sqrtInvR().norm() << ", |error| = " << error().norm() << ", errorImpl = " << errorImpl
         << endl;
#endif

    // place the hessian elements in the correct place:

    // build hessian:
    for (size_t i = 0; i < numDesignVariables(); i++) {
        if (designVariable(i)->isActive()) {
            // get the block index
            int colBlockIndex = designVariable(i)->blockIndex();
            int rows = designVariable(i)->minimalDimensions();
            int rowBase = outHessian.colBaseOfBlock(colBlockIndex);

            // <- this is our column index
            //_numberOfSplineDesignVariables
            for (size_t j = 0; j <= i; j++)  // upper triangle should be sufficient
            {
                if (designVariable(j)->isActive()) {
                    int rowBlockIndex = designVariable(j)->blockIndex();

                    // select the corresponding block in _Q:
                    Eigen::MatrixXd* Qblock = _Q.block(i, j, false);  // get block and do NOT allocate.

                    if (Qblock) {  // check if block exists
                        // get the Hessian Block
                        const bool allocateIfMissing = true;
                        Eigen::MatrixXd* Hblock = outHessian.block(rowBlockIndex, colBlockIndex, allocateIfMissing);
                        *Hblock += *Qblock;  // insert!
                    }
                }
            }

            outRhs.segment(rowBase, rows) -= b_u.segment(i * rows, rows);
        }
    }

    // std::cout << "OutHessian" << outHessian.toDense() << std::endl;

    // show outRhs:
    //  for (int i = 0;  i < outRhs.rows(); i++)
    //     std::cout << outRhs(i) <<  " : " << (*b_u)(i) << std::endl;
}

template <class SPLINE_T>
Eigen::VectorXd BSplineMotionError<SPLINE_T>::rhs() {
    const double* cMat = &((_splineDV->spline()).coefficients()(0, 0));
    Eigen::Map<const Eigen::VectorXd> c = Eigen::VectorXd::Map(cMat, _coefficientVectorLength);

    // right hand side:
    Eigen::VectorXd b_u(_Q.rows());  // number of rows of Q:
    b_u.setZero();
    _Q.multiply(&b_u, -c);

    printInfo("BSplineMotionError<SPLINE_T>::rhs()", DebugInfoType::Paragraph, true);
    cout << "c size = " << c.size() << ", |c| = " << c.norm() << endl;
    cout << "Q size = " << _Q.rows() << " X " << _Q.cols() << ", |Q| = " << _Q.toDense().norm() << endl;
    cout << "bu size = " << b_u.size() << ", |b_u| = " << b_u.norm() << endl;

    return b_u;
}

}  // namespace backend
}  // namespace aslam
