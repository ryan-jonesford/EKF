#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

/**
 * Name: Tools
 * Description: Tools class creates helper functions used for
 * 				implementing and evaluating Kalman filters
**/
Tools::Tools() {}

/**
 * Name: Tools
 * Description: Tools class Destructor
**/
Tools::~Tools() {}

/**
 * Name: CalculateRMSE
 * Return: Root mean squared values expressed as VectorXd
 * Description: Calculates the RMSE
**/
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  
 	int vector_size = estimations[0].size();
	VectorXd rmse = VectorXd::Zero(vector_size);

	// the estimation vector size should not be zero
	// the estimation vector size should equal ground truth vector size
	if( !estimations.size() ||
		estimations.size() != ground_truth.size())
 	{
		cout << "Tools::CalculateRMSE: Invalid estimation "
				"or ground_truth Vector" << endl;
		return rmse;
	}

	//accumulate squared residuals
	VectorXd residuals = VectorXd::Zero(vector_size);
	VectorXd resid(vector_size);
	for(int i=0; i < estimations.size(); ++i)
 	{
		resid = (estimations[i] - ground_truth[i]);
		resid = resid.array()*resid.array();
		residuals += resid;
	}

	//calculate the mean
	rmse = residuals.array()/vector_size;

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

/**
 * Name: CalculateJacobian
 * Return: Jacobian matrix expressed as MatrixXd
 * Description: Calculates a Jacobian matrix from a given state
**/
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

	MatrixXd Hj = MatrixXd::Zero(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);
	float px2py2 = px*px+py*py;
	float px2py2_sqrt = sqrt(px2py2);

	//check division by (or close to) zero
	if( fabs(px2py2) < 0.0001 )
  	{
		cout << "Tools::CalculateJacobian: Division by Zero detected!! Bailing out." << endl;
		return Hj;
	}
	
	//compute the Jacobian matrix
	Hj(0,0) = px/px2py2_sqrt;
	Hj(0,1) = py/px2py2_sqrt;
	Hj(1,0) = -py/px2py2;
	Hj(1,1) = px/px2py2;
	Hj(2,0) = (py*(vx*py-vy*px))/(pow(px2py2,(3.0/2.0)));
	Hj(2,1) = (px*(vy*px-vx*py))/(pow(px2py2,(3.0/2.0)));
	Hj(2,2) = px/px2py2_sqrt;
	Hj(2,3) = py/px2py2_sqrt;

	return Hj;
}

/**
 * Name: PolarToCartesian
 * Return: VectorXd in cartesian coordinates
 * Description: Convert Polar coords to cartesian
**/
VectorXd Tools::PolarToCartesian(const VectorXd& polar){
	float rho = polar(0);
	float phi = polar(1);
	// normalize angle between -pi and pi
	while (phi < -M_PI){
		phi += 2 * M_PI;
	}
	while ( phi > M_PI ){
		phi -= 2 * M_PI;
	}
	if( fabs(phi) > M_PI )
	{
		cout << "Tools::PolarToCartesian: Phi outside of -pi to pi: "
			 << phi << "\nBailing out." << endl;
		exit(1);
	}
	float px = cos(phi) + rho;
	float py = sin(phi) + rho;
    VectorXd cartesian = VectorXd(4);
	cartesian(0) = px;
	cartesian(1) = py;
	return cartesian;
}

/**
 * Name: CartesianToPolar
 * Return: VectorXd in Polar coordinates
 * Description: Convert cartesian coords to Polar
**/
VectorXd Tools::CartesianToPolar(const VectorXd& cartesian){
	// grab the cartesian coords
	float px = cartesian(0);
	float py = cartesian(1);
	float vx = cartesian(2);
	float vy = cartesian(3);

	// setup to convert to polar
	float px2py2 = px*px+py*py;
	//check division by (or close to) zero
	if( fabs(px2py2) < 0.0001 )
  	{
		cout << "Tools::CalculateJacobian: Division by Zero detected!! Bailing out." << endl;
		exit(1);
	}
	float sqrt_px2_py2 = sqrt(px2py2);
	float atan__py_div_px = atan2(py,px);

	// normalize angle between -pi and pi
	while (atan__py_div_px < -M_PI){
		atan__py_div_px += 2 * M_PI;
	}
	while ( atan__py_div_px > M_PI ){
		atan__py_div_px -= 2 * M_PI;
	}
	if( fabs(atan__py_div_px) > M_PI )
	{
		cout << "Tools::CalculateJacobian: Phi outside of -pi to pi: "
			 << atan__py_div_px << "\nBailing out." << endl;
		exit(1);
	}
	
	// populate the polar vector and return
	VectorXd polar = VectorXd(3);
	polar << 	sqrt_px2_py2,
				atan__py_div_px,
				( px*vx + py*vy ) / sqrt_px2_py2;

	return polar;
}