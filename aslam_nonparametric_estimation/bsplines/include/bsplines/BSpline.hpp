/**
 * @file   BSpline.hpp
 * @author Paul Furgale <paul.furgale@utoronto.ca>
 * @date   Fri Feb 11 13:51:57 2011
 *
 * @brief  A class to facilitate state estimation for vehicles in 3D
 *         space using B-splines.
 *
 *
 */

#ifndef _BSPLINE_HPP
#define _BSPLINE_HPP
#include <sparse_block_matrix/sparse_block_matrix.h>
#include <Eigen/Core>
#include <sm/assert_macros.hpp>
#include <vector>

namespace bsplines {
class BiVector;
}

namespace Eigen {
namespace internal {
template <>
struct functor_traits<bsplines::BiVector> {
    enum { Cost = 1, PacketAccess = false, IsRepeatable = true };
};
}  // namespace internal
}  // namespace Eigen

namespace bsplines {

class BiVector {
    enum { Cost = 1, PacketAccess = false, IsRepeatable = true };

  private:
    const int startIndex_;
    const double endValue_;
    const Eigen::VectorXd localBi_;

  public:
    BiVector(int startIndex, const Eigen::VectorXd& localBi, const double& endValue)
        : startIndex_(startIndex), endValue_(endValue), localBi_(localBi){};

    double operator()(int i, int j = 0) const {
        // kill unused parameter warning
        static_cast<void>(j);
        i -= startIndex_;
        if (i < 0) {
            return endValue_;
        }
        if (i >= (localBi_.rows())) {
            return 0;
        }
        return i >= 0 ? localBi_(i) : 0;
    }
};

/**
 * @brief A class to facilitate state estimation for vehicles in 3D space using B-Splines
 *
 * NOTE by CC:
 * (1) The B-Spline definition in the class is some different to wiki and MTU course
 *  (a) Basic function N[i,0](u) = 1 in TUM, B[i,1](t) = 1 in this definition
 *  (b) Basic function N[i,p](u) is a degree p polynomial in u, but B[i,k](t) is k-1 polynomial in this definition
 *  (c) m = n + p + 1 for MTU, m = n + k for this. And control points(coefficient) number is n+1, knots number is m+1,
 * spline degree is p
 * (2) Some terminology in this class
 *  (a) Spline degree: k
 *  (b) Polynomial degree: k - 1
 *  (c) m = n + k for this definition. To distinguish the knots number for open and clamped curve, I denote m + 1 for
 * the knots number of open curve, and m1 + 1 for clamped one, they have the relation m1 = m + 2(k-1)
 *  (d) Valid time segment: m
 *  (e) Knots numbers: m1 + 1 = m + 2(k - 1) + 1 = m + 2k - 1
 *  (f) Coefficient matrix, the same mean of control points matrix
 *  (g) Coefficient number: n + 1
 * (3) When evaluating the B-Spline value at t, only part of control points(ie, coefficient matrix) will active in the
 * calculation, Vk = {V[i-k+1], V[i-k+2],..., V[i]}. The matrix of Vk in code is called local parameters.
 *
 * I am more familiar with MTU course definition. So, although I have add some comment to this class, but
 * maybe they are not quite correct.
 * Some function about integral I could understand yet.
 */
class BSpline {
  public:
    /**
     * @brief A base class for BSpline exceptions
     */
    SM_DEFINE_EXCEPTION(Exception, std::runtime_error);

    /**
     * @brief Create a spline of the specified order. The resulting B-spline will be a series of piecewise polynomials
     * of degree splineOrder - 1.
     * @param splineOrder   The order of the spline
     */
    explicit BSpline(int splineOrder);

    /**
     * A destructor.
     */
    ~BSpline() = default;

    /**
     * @brief Get the order of the spline
     * @return The order of the spline
     */
    inline int splineOrder() const { return splineOrder_; }

    /**
     * @brief Get the degree of polynomial used by the spline
     * @return The degree of polynomial used by the spline
     */
    inline int polynomialDegree() const { return splineOrder_ - 1; }

    /**
     * @brief Get the number of coefficients required for a specified number of valid time segments
     * @param numTimeSegments The number of time segments required
     * @return The number of coefficients required for a specified number of valid time segments
     */
    int numCoefficientsRequired(int numTimeSegments) const;

    /**
     * @brief Get the minimum number of knots required to have at least one valid time segment
     * @return The minimum number of knots required to have at least one valid time segment
     */
    int minimumKnotsRequired() const;

    /**
     * @brief Get the number of knots required for a specified number of valid time segments
     * @param numTimeSegments The number of time segments required
     * @return The number of knots required for a specified number of valid time segments
     */
    int numKnotsRequired(int numTimeSegments) const;

    /**
     * @brief Get the number of valid time segments for a given number of knots
     * @param numKnots  The number of knots.
     * @return The number of valid time segments
     */
    int numValidTimeSegments(int numKnots) const;

    /**
     * @brief Get the number of valid time segments for a given for the current knot sequence
     * @return The number of valid time segments for a given for the current knot sequence
     */
    int numValidTimeSegments() const;

    /**
     * @brief Get the basic matrix active on the i-th time segment
     * @param i The index of the time segment
     * @return The basic matrix active on the i-th time segment
     */
    const Eigen::MatrixXd& basisMatrix(int i) const;

    /**
     * @brief Get the minimum time that the spline is well-defined on
     * @return The minimum time that the spline is well-defined on
     */
    const double& tMin() const;

    /**
     * @brief Get the maximum time that the spline is well-defined on. Because B-spline are defined on half-open
     * intervals, the spline curve is well defined up to but not including this time.
     * @return The minimum time that the spline is well-defined on
     */
    const double& tMax() const;

    /**
     * @brief Get the time interval that the spline is well-defined on [tMin, tMax)
     * @return The time interval that the spline is well-defined on [tMin, tMax)
     */
    std::pair<double, double> timeInterval() const;

    /**
     * @brief Get the time interval of a single spline segment
     * @param i The index of the time segment
     * @return The time interval of the i-th spline segment
     */
    std::pair<double, double> timeInterval(int i) const;

    /**
     * @brief Set the knots and coefficients of the spline, each column of the coefficient matrix is interpreted as a
     * single, vector-valued spline coefficient
     * @param knots         A non-decreasing knot sequence
     * @param coefficients  A set of spline coefficients
     */
    void setKnotsAndCoefficients(const std::vector<double>& knots, const Eigen::MatrixXd& coefficients);

    /**
     * @brief Set the knots and coefficients of the spline, each column of the coefficient matrix is interpreted as a
     * single, vector-valued spline coefficient.
     * @param knots         A non-decreasing knot sequence
     * @param coefficients  A set of spline coefficients
     */
    void setKnotVectorAndCoefficients(const Eigen::VectorXd& knots, const Eigen::MatrixXd& coefficients);

    /**
     * @brief Set the coefficient matrix
     * @param coefficients Coefficient matrix
     */
    void setCoefficientMatrix(const Eigen::MatrixXd& coefficients);

    /**
     * @brief Sets the coefficient matrix from the stacked vector of coefficients
     * @param coefficients The stacked vector of coefficients
     */
    void setCoefficientVector(const Eigen::VectorXd& coefficients);

    /**
     * @brief Get the knot vector
     * @return Knot vector
     */
    inline const std::vector<double>& knots() const { return knots_; }

    /**
     * @brief Get the knot vector with a column matrix
     * @return Knot vector with a column matrix
     */
    Eigen::VectorXd knotVector() const;

    /**
     * @brief Get the coefficient matrix, each column of the coefficient matrix is interpreted as a single,
     * vector-valued spline coefficient.
     * @return Coefficient matrix
     */
    inline const Eigen::MatrixXd& coefficients() const { return coefficients_; }

    /**
     * @brief Get the stacked vector of coefficients
     * @return Stacked vector of coefficients
     */
    Eigen::VectorXd coefficientVector();

    /**
     * @brief Get the number of total coefficients the spline currently uses
     * @return The number of total coefficients the spline currently uses
     */
    int numCoefficients() const;

    int coefficientVectorLength() const;

    /**
     * @brief This is equivalent to spline.coefficients().cols()  = n
     * @return The number of vector-valued coefficient columns the spline currently uses
     */
    int numVvCoefficients() const;

    /**
     * @brief Get a map to a single coefficient column. This allows the user to pass around what is essentially a
     * pointer to a single column in the coefficient matrix.
     * @param i The column of the coefficient matrix to return. 0 < i < coefficients().cols() = n
     * @return A map to column i of the coefficient matrix.
     */
    Eigen::Map<const Eigen::VectorXd> vvCoefficientVector(int i) const;

    /**
     * @brief Get a map to a single coefficient column. This allows the user to pass around what is essentially a
     * pointer to a single column in the coefficient matrix.
     * @param i The column of the coefficient matrix to return. 0 < i < coefficients().cols() = n
     * @return A map to column i of the coefficient matrix.
     */
    Eigen::Map<Eigen::VectorXd> vvCoefficientVector(int i);

    /**
     * @brief Evaluate the spline curve at time t
     * @param t The time to evaluate the spline curve
     * @return The value of the spline curve at the time t
     */
    Eigen::VectorXd eval(double t) const;

    /**
     * @brief Evaluate the derivative of the spline curve at time t
     * @param t                 The time to evaluate the spline derivative
     * @param derivativeOrder   The order of the derivative. This must be >= 0
     * @return The value of the derivative of the spline curve evaluated at t
     */
    Eigen::VectorXd evalD(double t, int derivativeOrder) const;

    /**
     * @brief Evaluate the derivative of the spline curve at time t and retrieve the Jacobian of the value with respect
     * to small changed in the parameter vector(control points, coefficient matrix). The Jacobian only refers to the
     * local parameter vector (part of coefficient matrix). The indices of the local parameters with respect to the full
     * parameter vector can be retrieved using localCoefficientVectorIndices()
     * @param t                 The time to evaluate the spline derivative
     * @param derivativeOrder   The order of the derivative. This must be >= 0
     * @return The value of the derivative of the spline curve evaluated at t and the Jacobian
     */
    std::pair<Eigen::VectorXd, Eigen::MatrixXd> evalDAndJacobian(double t, int derivativeOrder) const;

    /**
     * @brief Evaluate the derivative of the spline curve at time t and retrieve the Jacobian of the value with respect
     * to small changed in the parameter vector(control points, coefficient matrix). The Jacobian only refers to the
     * local parameter vector (part of coefficient matrix)
     * @param t                     The time to evaluate the spline derivative
     * @param derivativeOrder       The order of the derivative. This must be >= 0
     * @param jacobian              A pointer to the Jacobian matrix to fill in
     * @param coefficientIndices    A pointer to an int vector that will be filled with the local coefficient indices
     * (which index is activated for calculating Jacobian matrix)
     * @return The value of the derivative of the spline curve evaluated at t
     */
    Eigen::VectorXd evalDAndJacobian(double t, int derivativeOrder, Eigen::MatrixXd* jacobian,
                                     Eigen::VectorXi* coefficientIndices) const;

    /**
     * @brief Get the local basic matrix evaluated at the time t. For vector-valued spline coefficients of dimension n1,
     * and a B-Spline of order k, this matrix will be n1 x (k*n1).
     *
     * Evaluating the B-spline at time t, eval(t,O) is equivalent to evaluating Phi(t,O) * localCoefficientVector(t)
     *
     * @param t                 The time to evaluate the local basis matrix.
     * @param derivativeOrder   The derivative order to return (0 is no derivative)
     * @return  The local basis matrix evaluated at time t
     */
    Eigen::MatrixXd Phi(double t, int derivativeOrder = 0) const;

    /**
     * @brief Get the local basic matrix evaluated at the time t. For vector-valued spline coefficients of dimension n1,
     * and a B-Spline of order k, this matrix will be n1 x (k*n1).
     *
     * Evaluating the B-spline at time t, eval(t,O) is equivalent to evaluating Phi(t,O) * localCoefficientVector(t)
     *
     * @param t                 The time to evaluate the local basis matrix.
     * @param derivativeOrder   The derivative order to return (0 is no derivative)
     * @return  The local basis matrix evaluated at time t
     */
    Eigen::MatrixXd localBasisMatrix(double t, int derivativeOrder = 0) const;

    /**
     * @brief Get the local coefficient matrix evaluated at the time t. For vector-valued spline coefficients of
     * dimension n1, and a B-Spline of order k, this matrix will be n1 x k
     * @param t The time being queried
     * @return The local coefficient matrix active at time t
     */
    Eigen::MatrixXd localCoefficientMatrix(double t) const;

    /**
     * @brief Get the local coefficient vector evaluated at the time t. For vector-valued spline coefficients of
     * dimension n1, and a B-Spline of order k, this vector will be k*n1 x 1. Evaluating the B-spline at time t,
     * eval(t,O) is equivalent to evaluating Phi(t,O) * localCoefficientVector(t)
     * @param t The time being queried
     * @return The local coefficient vector active at time t
     */
    Eigen::VectorXd localCoefficientVector(double t) const;

    /**
     * @brief Get the local coefficient vector for segment i
     * @param segmentIdx The segment index
     * @return The local coefficient vector active on time segment i
     */
    Eigen::VectorXd segmentCoefficientVector(int segmentIdx) const;

    /**
     * @brief Return a map to a single coefficient column.
     *
     * This allows the user to pass around what is essentially a pointer to a single column in the coefficient matrix
     * @tparam D    The dimension of vector-valued spline coefficients n1
     * @param i     The column of the coefficient matrix to return. 0 <= i < n = coefficients().cols()
     * @return  A map to column i of the coefficient matrix
     */
    template <int D>
    Eigen::Map<Eigen::Matrix<double, D, 1> > fixedSizeVvCoefficientVector(int i) {
        SM_ASSERT_EQ_DBG(Exception, D, coefficients_.rows(),
                         "Size mismatch between requested vector size and actual vector size");
        SM_ASSERT_GE_LT(Exception, i, 0, coefficients_.cols(), "Index out of range");
        return Eigen::Map<Eigen::Matrix<double, D, 1> >(&coefficients_(0, i), coefficients_.rows());
    }

    /**
     * @brief Return a map to a single coefficient column, const version.
     *
     * This allows the user to pass around what is essentially a pointer to a single column in the coefficient matrix
     * @tparam D    The dimension of vector-valued spline coefficients n1
     * @param i     The column of the coefficient matrix to return. 0 <= i < n = coefficients().cols()
     * @return  A map to column i of the coefficient matrix
     */
    template <int D>
    Eigen::Map<const Eigen::Matrix<double, D, 1> > fixedSizeVvCoefficientVector(int i) const {
        SM_ASSERT_EQ_DBG(Exception, D, coefficients_.rows(),
                         "Size mismatch between requested vector size and actual vector size");
        SM_ASSERT_GE_LT(Exception, i, 0, coefficients_.cols(), "Index out of range");
        return Eigen::Map<const Eigen::Matrix<double, D, 1> >(&coefficients_(0, i), coefficients_.rows());
    }

    /**
     * @brief Get the indices of the local coefficients active at time t. Only part of control points(coefficient
     * matrix) is activated in evaluating B-Spline at t, ie, Vk = { V[i-k+1], V[i-k+2], ... , V[i] }, the indices is the
     * index of every element of Vk with respect to the coefficient matrix, ie, [(i-k+1) * n1, (i+1) * n1)
     * @param t The time being queried
     * @return The indices of the local coefficients active at time t
     */
    Eigen::VectorXi localCoefficientVectorIndices(double t) const;

    /**
     * @brief Get the indices of the local coefficients active on segment i
     * @param segmentIdx The segment being queried
     * @return The indices of the local coefficients active on this segment
     */
    Eigen::VectorXi segmentCoefficientVectorIndices(int segmentIdx) const;

    /**
     * @brief Get the indices of the local vector-valued coefficients active at time t. Only part of control
     * points(coefficient matrix) is activated in evaluating B-Spline at t, ie, Vk = { V[i-k+1], V[i-k+2], ... , V[i] },
     * the indices is the row index of Vj with respect to the coefficient matrix. ie, [i-k+1, i+1)
     * @param t The time being queried
     * @return The indices of the local vector-valued coefficients active at time t
     */
    Eigen::VectorXi localVvCoefficientVectorIndices(double t) const;

    /**
     * @brief Get the indices of the local vector-valued coefficients active on segment i
     * @param segmentIdx The segment being queried
     * @return The indices of the local vector-valued coefficients active at time t
     */
    Eigen::VectorXi segmentVvCoefficientVectorIndices(int segmentIdx) const;

    /**
     * @brief Update the local coefficient vector
     * @param t The time used to select the local coefficients.
     * @param c The local coefficient vector
     */
    void setLocalCoefficientVector(double t, const Eigen::VectorXd& c);

    // NOTE by CC: I cannot understand below 4 init spline methods right now. Maybe it need B-Spline curve fitting and
    // interpolation background knowledge, and the conversion between B-Spline and Bezier. Furthermore, the lambda
    // meaning is also problem.

    /**
     * Initialize a spline from two times and two positions. The spline will be initialized to
     * have one valid time segment \f$[t_0, t_1)\f$ such that \f$\mathbf b(t_0) = \mathbf p_0\f$,
     * \f$\mathbf b(t_1) = \mathbf p_1\f$,
     * \f$\dot{\mathbf b}(t_0) = \frac{\mathbf{p_1} - \mathbf p_0}{t_1 - t_0}\f$, and
     * \f$\dot{\mathbf b}(t_1) = \frac{\mathbf{p_1} - \mathbf p_0}{t_1 - t_0}\f$.
     *
     * @param t0 The start of the time interval.
     * @param t1 The end of the time interval
     * @param p0 The position at the start of the time interval.
     * @param p1 The position at the end of the time interval.
     */
    void initSpline(double t0, double t1, const Eigen::VectorXd& p0, const Eigen::VectorXd& p1);

    void initSpline2(const Eigen::VectorXd& times, const Eigen::MatrixXd& interpolationPoints, int numSegments,
                     double lambda);

    void initSpline3(const Eigen::VectorXd& times, const Eigen::MatrixXd& interpolationPoints, int numSegments,
                     double lambda);

    void initSplineSparse(const Eigen::VectorXd& times, const Eigen::MatrixXd& interpolationPoints, int numSegments,
                          double lambda);

    void initSplineSparseKnots(const Eigen::VectorXd& times, const Eigen::MatrixXd& interpolationPoints,
                               const Eigen::VectorXd knots, double lambda);

    /**
     * Add a curve segment that interpolates the point p, ending at time t.
     *
     * If the new time corresponds with the first knot past the end of the curve,
     * the existing curve is perfectly preserved. Otherwise, the existing curve
     * will interpolate its current position at the current endpoint and the new
     * position at the new endpoint but will not necessarily match the last segment
     * exactly.
     *
     * @param t The time of the point to interpolate. This must be greater than t_max()
     * @param p The point to interpolate at time t.
     */
    void addCurveSegment(const double& t, const Eigen::VectorXd& p);

    /**
     * Add a curve segment that interpolates the point p, ending at time t.
     *
     * If the new time corresponds with the first knot past the end of the curve,
     * the existing curve is perfectly preserved. Otherwise, the existing curve
     * will interpolate its current position at the current endpoint and the new
     * position at the new endpoint but will not necessarily match the last segment
     * exactly.
     *
     * @param t The time of the point to interpolate. This must be greater than t_max()
     * @param p The point to interpolate at time t.
     * @param lambda a smoothness parameter. Higher for more smooth.
     */
    void addCurveSegment2(const double& t, const Eigen::VectorXd& p, const double& lambda);

    /**
     * @brief Removes a curve segment from the left by removing one knot and one coefficient vector.
     * After calling this function, the curve will have one fewer segment. The new minimum time will be
     * timeInterval(0).first
     */
    void removeCurveSegment();

    /**
     * Get the \f$ \mathbf V_i \f$ matrix associated with the integral over the segment.
     *
     * @param segmentIndex
     *
     * @return the \f$ \mathbf V_i \f$ matrix
     */
    Eigen::MatrixXd Vi(int segmentIndex) const;

    Eigen::VectorXd evalIntegral(double t1, double t2) const;
    inline Eigen::VectorXd evalI(double t1, double t2) const { return evalIntegral(t1, t2); }

    Eigen::MatrixXd Mi(int segmentIndex) const;
    Eigen::MatrixXd Bij(int segmentIndex, int columnIndex) const;
    Eigen::MatrixXd U(double t, int derivativeOrder) const;
    Eigen::VectorXd u(double t, int derivativeOrder) const;
    int segmentIndex(double t) const;
    Eigen::MatrixXd Dii(int segmentIndex) const;
    Eigen::MatrixXd Di(int segmentIndex) const;

    /**
     * Get the b_i(t) for i in localVvCoefficientVectorIndices (@see #localVvCoefficientVectorIndices).
     *
     * @param t The time being queried.
     *
     * @return [b_i(t) for i in localVvCoefficientVectorIndices].
     *
     */
    Eigen::VectorXd getLocalBiVector(double t) const;
    void getLocalBiInto(double t, Eigen::VectorXd& ret) const;

    /**
     * Get the cumulative (tilde) b_i(t) for i in localVvCoefficientVectorIndices (@see
     * #localVvCoefficientVectorIndices).
     *
     * @param t The time being queried.
     *
     * @return [tilde b_i(t) for i in localVvCoefficientVectorIndices].
     *
     */
    Eigen::VectorXd getLocalCumulativeBiVector(double t) const;

    Eigen::CwiseNullaryOp<BiVector, Eigen::VectorXd> getBiVector(double t) const {
        return Eigen::CwiseNullaryOp<BiVector, Eigen::VectorXd>(numValidTimeSegments(), 1,
                                                                BiVector(segmentIndex(t), getLocalBiVector(t), 0));
    }

    Eigen::CwiseNullaryOp<BiVector, Eigen::VectorXd> getCumulativeBiVector(double t) const {
        return Eigen::CwiseNullaryOp<BiVector, Eigen::VectorXd>(
            numValidTimeSegments(), 1, BiVector(segmentIndex(t), getLocalCumulativeBiVector(t), 1));
    }

    Eigen::MatrixXd segmentIntegral(int segmentIdx, const Eigen::MatrixXd& W, int derivativeOrder) const;

    Eigen::MatrixXd segmentQuadraticIntegral(const Eigen::MatrixXd& W, int segmentIdx, int derivativeOrder) const;
    Eigen::MatrixXd segmentQuadraticIntegralDiag(const Eigen::VectorXd& Wdiag, int segmentIdx,
                                                 int derivativeOrder) const;
    Eigen::MatrixXd curveQuadraticIntegral(const Eigen::MatrixXd& W, int derivativeOrder) const;
    Eigen::MatrixXd curveQuadraticIntegralDiag(const Eigen::VectorXd& Wdiag, int derivativeOrder) const;

    sparse_block_matrix::SparseBlockMatrix<Eigen::MatrixXd> curveQuadraticIntegralSparse(const Eigen::MatrixXd& W,
                                                                                         int derivativeOrder) const;
    sparse_block_matrix::SparseBlockMatrix<Eigen::MatrixXd> curveQuadraticIntegralDiagSparse(
        const Eigen::VectorXd& Wdiag, int derivativeOrder) const;

    void initConstantSpline(double t_min, double t_max, int numSegments, const Eigen::VectorXd& constant);

  private:
    /**
     * @brief Check the knot sequence is valid, it will throw an exception if it's not valid
     * @param knots The knot sequence to verify
     */
    void verifyKnotSequence(const std::vector<double>& knots);

    /**
     * @brief Initialize the basis matrices based on the current knot sequence. There is one basis matrix for each valid
     * time segment defined by the spline.
     *
     * Implemented using the recursive basis matrix algorithm from
     * Qin, Kaihuai, General matrix representations for B-splines, The Visual Computer (2000) 16:177–186
     */
    void initializeBasisMatrices();

    /**
     * @brief The recursive function to calculate basis matrix M.
     *           [ M_{k-1}(i) ]        [      0       ]
     *  M_k(i) = [            ] * A +  [              ] * B = M1 * A + M2 * B
     *           [     0      ]        [ [ M_{k-1}(i) ]
     *
     *  M_1(i) = [1]
     *  size of A and B: (k-1) X k
     *
     * Ref: Qin, Kaihuai, General matrix representations for B-splines, The Visual Computer (2000) 16:177–186.
     *
     * @param k The degree(order) of the matrix
     * @param i The time segment index
     * @return  Basic matrix
     */
    Eigen::MatrixXd M(int k, int i);

    /**
     * @brief A helper function to calculate d0 for producing th basic matrix M.
     *  d(0,j) = (t[i] - t[j]) / (t[j+k-1] - t[j]).
     * Defined in Qin, Kaihuai, General matrix representations for B-splines, The Visual Computer (2000) 16:177–186.
     *
     * @param k Spline degree(order)
     * @param i Time segment index
     * @param j Time segment index
     * @return  Value d0
     */
    double d0(int k, int i, int j);

    /**
     * @brief A helper function to calculate d1 for producing the M matrices.
     *  d(1,j) = (t[i+1] - t[i]) / (t[j+k-1] - t[j]).
     * Defined in Qin, Kaihuai, General matrix representations for B-splines, The Visual Computer (2000) 16:177–186.
     *
     * @param k Spline degree(order)
     * @param i Time segment index
     * @param j Time segment index
     * @return  Value d1
     */
    double d1(int k, int i, int j);

    /**
     * @brief An internal function to find the segment of the knot sequence that the time t falls in. The function
     * returns the value u = (t - t[i]) / (t[i+1] - t[i]) and the index i
     * @param t The time being queried
     * @return A pair with the first value u = (t - t[i]) / (t[i+1] - t[i]) and the second value the index i
     */
    std::pair<double, int> computeUAndTIndex(double t) const;

    /**
     * @brief An internal function to find the segment of the knot sequence that the time t falls in. The function
     * returns the width of the knot segment dt = t[i+1] - t[i] and the index i
     * @param t The time being queried
     * @return  The pair with the first value dt = t[i+1] -t[i] and the second value the index i
     */
    std::pair<double, int> computeTIndex(double t) const;

    /**
     * @brief Compute the vector U for a spline of order k, this is an k X 1 vector.
     *
     * At derivative order 0 (no derivative), this vector is U(t) = [1, u(t), u(t)^2, ... , u(t)^{k-1}]^T.
     *
     * For higher derivative order n, the vector is U^{(n)}(t) = d^n U(t) / dt^{n}
     * @param u                 The value u(t)
     * @param segmentIndex      The index of t
     * @param derivativeOrder   Derivative order
     * @return The vector of U
     */
    Eigen::VectorXd computeU(double u, int segmentIndex, int derivativeOrder) const;

    int basisMatrixIndexFromStartingKnotIndex(int startingKnotIndex) const;
    int startingKnotIndexFromBasisMatrixIndex(int basisMatrixIndex) const;
    const Eigen::MatrixXd& basisMatrixFromKnotIndex(int knotIndex) const;

  private:
    int splineOrder_;            // The order of the spline
    std::vector<double> knots_;  // The knot sequence used by the B-spline
    // The basis matrices for each time segment the B-spline is defined over
    std::vector<Eigen::MatrixXd> basisMatrices_;

    // The coefficient matrix(control points) used by the B-Spline. Each column can be seen as a single vector-valued
    // spline coefficient. This is stored explicitly in column major order to ensure that each column (i.e. a single
    // vector-valued spline coefficient) is stored in contiguous memory. This allows one to, for example, map a single
    // spline coefficient using the Eigen::Map type.
    // NOTE by CC: the control points is V = {V_0, V_1,.., V_n], consider Vj is size of n1 x 1, then the size of V is
    // n1 x n
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> coefficients_;
};

}  // namespace bsplines

#endif /* _BSPLINE_HPP */
