#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include <cmath>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  MatrixXd CalculateJacobian(const VectorXd& x_state);

  /**
  * A helper method to convert Polar coords to cartesian.
  */
  VectorXd PolarToCartesian(const VectorXd& polar);

  /**
  * A helper method to convert cartesian coords to Polar.
  */
  VectorXd CartesianToPolar(const VectorXd& cartesian);

  /**
  * A helper method to normalize an angle between -Pi and Pi
  */
  double normalizePhi(double phi);

};

#endif /* TOOLS_H_ */
