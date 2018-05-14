#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

/**
 * Name: Predict
 * Return: None
 * Description: Implements Kalman Filter Prediction
**/
void KalmanFilter::Predict() {
  x_ = (F_*x_);
  P_ = F_*P_*F_.transpose() + Q_;
}

/**
 * Name: Update
 * Return: None
 * Description: Updates the state by using standard Kalman Filter equations
**/
void KalmanFilter::Update(const VectorXd &z) {
  
  VectorXd z_pred = H_ * x_;
	VectorXd y      = z - z_pred;
	MatrixXd Ht     = H_.transpose();
	MatrixXd S      = H_ * P_ * Ht + R_;
	MatrixXd PHt    = P_ * Ht;
	MatrixXd K      = PHt * S.inverse();

	//new estimate
	x_ = x_ + (K * y);
	int x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

/**
 * Name: UpdateEKF
 * Return: None
 * Description: Updates the state by using Extended Kalman Filter equations
**/
void KalmanFilter::UpdateEKF(const VectorXd &z) {
  MatrixXd y  = z - tools.CartesianToPolar(x_);
  while (y(1) < -M_PI){
		y(1) += 2 * M_PI;
	}
	while ( y(1) > M_PI ){
		y(1) -= 2 * M_PI;
	}
	if( fabs(y(1)) > M_PI )
	{
		cout << "Tools::PolarToCartesian: y(1) outside of -pi to pi: "
			 << y(1) << "\nBailing out." << endl;
		exit(1);
	}
	MatrixXd Ht = H_.transpose();
  MatrixXd S  = H_ * P_ * Ht + R_;
  MatrixXd K  = P_ * Ht * S.inverse();

	//new estimate
	x_ = x_ + (K * y);
	int x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}
