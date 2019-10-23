#pragma once
#include <fmt/format.h>
#include <glog/logging.h>
#include <sm/boost/null_deleter.hpp>
#include "aslam/backend/EuclideanPoint.hpp"
#include "aslam/backend/RotationQuaternion.hpp"
#include "aslam/backend/TransformationBasic.hpp"
#include "sm/kinematics/Transformation.hpp"

namespace cc {

class TransformationDesignVariable {
  public:
    explicit TransformationDesignVariable(const sm::kinematics::Transformation& transformation,
                                          bool rotationActive = true, bool translationActive = true)
        : initTransformation(transformation),
          q(transformation.q()),
          t(transformation.t()),
          basicDesignVar(q.toExpression(), t.toExpression()) {
        q.setActive(rotationActive);
        t.setActive(translationActive);
        expression = basicDesignVar.toExpression();
    }

    inline aslam::backend::TransformationExpression& toExpression() { return expression; }
    inline int numDesignVariables() const { return 2; }
    Eigen::Matrix4d T() const { return expression.toTransformationMatrix(); }

    aslam::backend::DesignVariable::Ptr designVariable(int i) {
        if (i == 0) {
            return aslam::backend::DesignVariable::Ptr(&q, sm::null_deleter());
        } else if (i == 1) {
            return aslam::backend::DesignVariable::Ptr(&t, sm::null_deleter());
        } else {
            LOG(ERROR) << fmt::format("index out of range: {} >= 2", i);
        }
    }

  private:
    sm::kinematics::Transformation initTransformation;
    aslam::backend::RotationQuaternion q;
    aslam::backend::EuclideanPoint t;
    aslam::backend::TransformationBasic basicDesignVar;
    aslam::backend::TransformationExpression expression;
};

}  // namespace cc