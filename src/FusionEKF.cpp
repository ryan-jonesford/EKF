#include "FusionEKF.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_  = MatrixXd(2, 2);
  R_radar_  = MatrixXd(3, 3);
  H_laser_  = MatrixXd(2, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0,      0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0,      0,
              0,    0.0009, 0,
              0,    0,      0.09;
  
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // hard code noise for now
  noise_ax_ = 9, noise_ay_ = 9;

  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ <<  1, 0, 1, 0,
              0, 1, 0, 1,
              0, 0, 1, 0,
              0, 0, 0, 1;
  ekf_.P_ = MatrixXd::Zero(4,4);
  ekf_.Q_ = MatrixXd::Zero(4,4);

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}


/**
 * Name: ProcessMeasurement
 * Return: None
 * Description: Updates the state predictions from measurements
**/
void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    // first measurement
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    // set timestamp
    previous_timestamp_ = measurement_pack.timestamp_;
    
    // Use a switch statement for sensor type switching
    switch (measurement_pack.sensor_type_)
    {
      case (MeasurementPackage::RADAR):
      {
        VectorXd px_py = tools.PolarToCartesian(measurement_pack.raw_measurements_);
        ekf_.x_(0) =  px_py(0);
        ekf_.x_(1) =  px_py(1);

        ekf_.P_ <<  10, 0, 0, 0,
                    0, 10, 0, 0,
                    0, 0, 100, 0,
                    0, 0, 0, 100;
        break;
      }
      case (MeasurementPackage::LASER):
      {
        ekf_.x_(0) = measurement_pack.raw_measurements_(0);
        ekf_.x_(1) = measurement_pack.raw_measurements_(1);
        ekf_.P_ <<  1, 0, 0, 0,
                    0, 1, 0, 0,
                    0, 0, 1000, 0,
                    0, 0, 0, 1000;
        break; 
      }
      default:
        cerr << "FusionEKF::ProcessMeasurement: Invalid Measurment device" << endl;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  //Update time, dt_1 expressed in seconds
  float dt_1 = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  // If dt_1 is out of bounds, end program
  if( dt_1 < 0 || dt_1 > MAX_DT )
  {
    cerr << "FusionEKF::ProcessMeasurement: Something went wrong... dt = "
         << dt_1 << endl;
    exit(1);
  }

  // Update the state transition matrix F
  ekf_.F_(0,2) = dt_1;
  ekf_.F_(1,3) = dt_1;

  // Update the process noise covariance matrix
  float dt_2 = dt_1 * dt_1;
  float dt_3 = dt_2 * dt_1;
  float dt_4 = dt_3 * dt_1;

  ekf_.Q_(0,0) = ( dt_4 / 4.0 ) * noise_ax_;
  ekf_.Q_(0,2) = ( dt_3 / 2.0 ) * noise_ax_;
  ekf_.Q_(1,1) = ( dt_4 / 4.0 ) * noise_ay_;
  ekf_.Q_(1,3) = ( dt_3 / 2.0 ) * noise_ay_;
  ekf_.Q_(2,0) = ( dt_3 / 2.0 ) * noise_ax_;
  ekf_.Q_(2,2) = ( dt_2 / 1.0 ) * noise_ax_;
  ekf_.Q_(3,1) = ( dt_3 / 2.0 ) * noise_ay_;
  ekf_.Q_(3,3) = ( dt_2 / 1.0 ) * noise_ay_;

  ekf_.Predict();
  
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  // Use a switch statement for sensor type
  switch (measurement_pack.sensor_type_)
  {
    case (MeasurementPackage::RADAR):
    {
      ekf_.R_ = R_radar_;
      ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
      break;
    }
    case (MeasurementPackage::LASER):
    {
      ekf_.R_=R_laser_;
      ekf_.H_=H_laser_;
      ekf_.Update(measurement_pack.raw_measurements_);
      break; 
    }
    default:
      cerr << "FusionEKF::ProcessMeasurement: Invalid Measurment device" << endl;
  }

  // print the output
  // cout << "x_ = " << ekf_.x_ << endl;
  // cout << "P_ = " << ekf_.P_ << endl;
}
