#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
   * End DO NOT MODIFY section for measurement noise values
   */

   n_x_ = 5;
   n_aug_ = n_x_ + 2;
   lambda_ = 3 - n_aug_;
   Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
   weights_ = VectorXd(2 * n_aug_ + 1);
   time_us_ = 0;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

   if (!is_initialized_) {
     P_ << 1, 0, 0, 0, 0,
           0, 1, 0, 0, 0,
           0, 0, 0.3 * 0.3, 0, 0,
           0, 0, 0, 0.03 * 0.03, 0,
           0, 0, 0, 0, 1;
     if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
       double rho = meas_package.raw_measurements_(0);
       double phi = meas_package.raw_measurements_(1);
       double rho_d = meas_package.raw_measurements_(2);
       double x = rho * cos(phi);
       double y = rho * sin(phi);
       double v = rho_d;
       x_ << x, y, v, 0, 0;
     } else {
       x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0, 0, 0;
     }
     is_initialized_ = true;
     return;
   }
   double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
   time_us_ = meas_package.timestamp_;
   Prediction(delta_t);

   if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
     UpdateRadar(meas_package);
   }
   if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
     UpdateLidar(meas_package);
   }
}

void UKF::Prediction(double delta_t) {

   // Augmented state matrix
   VectorXd x_aug = VectorXd(n_aug_);
   x_aug.fill(0);
   x_aug.head(n_x_) = x_;

   // Augmented state covariance matrix
   MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
   P_aug.fill(0);
   MatrixXd Q_ = MatrixXd(2, 2);
   Q_ << std_a_ * std_a_, 0,
         0, std_yawdd_ * std_yawdd_;
   P_aug.topLeftCorner(n_x_, n_x_) = P_;
   P_aug.bottomRightCorner(2, 2) = Q_;

   // Generating sigma points
   MatrixXd Xsig = MatrixXd(n_aug_, 2 * n_aug_ + 1);
   Xsig.fill(0);
   Xsig.col(0) = x_aug;

   // Calculating square root of P_aug
   MatrixXd A = MatrixXd(n_aug_, n_aug_);
   A = P_aug.llt().matrixL();

   for (int i = 0; i < n_aug_; i++) {
     Xsig.col(i + 1)          = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
     Xsig.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
   }

   // Predicting Sigma points
   Xsig_pred_.fill(0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     double px = Xsig(0, i);
     double py = Xsig(1, i);
     double v = Xsig(2, i);
     double psi = Xsig(3, i);
     double psi_d = Xsig(4, i);

     double px_p, py_p;
     double v_p = v;
     double psi_p = psi + (psi_d * delta_t);
     double psi_d_p = psi_d;

     if (fabs(psi_d) > 0.001) {
       px_p = px + ((v / psi_d) * (sin(psi + (psi_d * delta_t)) - sin(psi)));
       py_p = py + ((v / psi_d) * (- cos(psi + (psi_d * delta_t)) + cos(psi)));
     } else {
       px_p = px + (v * cos(psi) * delta_t);
       py_p = py + (v * sin(psi) * delta_t);
     }
     px_p = px_p + (0.5 * delta_t * delta_t * cos(psi) * Xsig(5, i));
     py_p = py_p + (0.5 * delta_t * delta_t * sin(psi) * Xsig(5, i));
     v_p = v_p + (delta_t * Xsig(5, i));
     psi_p = psi_p + (0.5 * delta_t * delta_t * Xsig(6, i));
     psi_d_p = psi_d_p + (delta_t * Xsig(6, i));

     Xsig_pred_(0, i) = px_p;
     Xsig_pred_(1, i) = py_p;
     Xsig_pred_(2, i) = v_p;
     Xsig_pred_(3, i) = psi_p;
     Xsig_pred_(4, i) = psi_d_p;
   }

   // Calculating state x
   weights_(0) = lambda_ / (lambda_ + n_aug_);
   for (int i = 1; i < 2 * n_aug_ + 1; i++) {
     weights_(i) = 0.5 / (lambda_ + n_aug_);
   }
   x_.fill(0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     x_ = x_ + weights_(i) * Xsig_pred_.col(i);
   }

   // Calculating state covariance
   P_.fill(0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     VectorXd x_diff = Xsig_pred_.col(i) - x_;
     while (x_diff(3) > M_PI) x_diff(3) = x_diff(3) - 2 * M_PI;
     while (x_diff(3) < -M_PI) x_diff(3) = x_diff(3) + 2 * M_PI;
     P_ = P_ + (weights_(i) * x_diff * (x_diff.transpose()));
   }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {

   VectorXd z_ = meas_package.raw_measurements_;
   int n_z = 2;
   MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
   Zsig.fill(0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     Zsig(0, i) = Xsig_pred_(0, i);
     Zsig(1, i) = Xsig_pred_(1, i);
   }

   // Calculating mean
   VectorXd z_pred = VectorXd(n_z);
   z_pred.fill(0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     z_pred = z_pred + weights_(i) * Zsig.col(i);
   }

   // Calculating covariance
   MatrixXd S = MatrixXd(n_z, n_z);
   S.fill(0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     VectorXd z_diff = Zsig.col(i) - z_pred;
     S = S + (weights_(i) * z_diff * z_diff.transpose());
   }
   MatrixXd R = MatrixXd(n_z, n_z);
   R << std_laspx_ * std_laspx_, 0,
        0, std_laspy_ * std_laspy_;

   S = S + R;

   // Calculating cross correlation matrix
   MatrixXd T = MatrixXd(n_x_, n_z);
   T.fill(0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     VectorXd x_diff = Xsig_pred_.col(i) - x_;
     VectorXd z_diff = Zsig.col(i) - z_pred;
     T = T + (weights_(i) * x_diff * (z_diff.transpose()));
   }

   // Calculating Kalman Gain
   MatrixXd K = T * (S.inverse());

   // Calculating updated state mean
   VectorXd z_diff = z_ - z_pred;
   x_ = x_ + K * z_diff;

   // Calculating updated state covariance
   P_ = P_ - (K * S * (K.transpose()));
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {

   VectorXd z_ = meas_package.raw_measurements_;
   int n_z = 3;
   MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
   Zsig.fill(0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     double px = Xsig_pred_(0, i);
     double py = Xsig_pred_(1, i);
     double v = Xsig_pred_(2, i);
     double psi = Xsig_pred_(3, i);
     double psi_d = Xsig_pred_(4, i);

     double rho = sqrt(px * px + py * py);
     double phi = atan2(py, px);
     double rho_d = (px * cos(psi) * v + py * sin(psi) * v) / (sqrt(px * px + py * py));
     Zsig(0, i) = rho;
     Zsig(1, i) = phi;
     Zsig(2, i) = rho_d;
   }

   // Calculating mean
   VectorXd z_pred = VectorXd(n_z);
   z_pred.fill(0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     z_pred = z_pred + weights_(i) * Zsig.col(i);
   }

   // Calculating covariance
   MatrixXd S = MatrixXd(n_z, n_z);
   S.fill(0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     VectorXd z_diff = Zsig.col(i) - z_pred;
     while (z_diff(1) > M_PI) z_diff(1) = z_diff(1) - 2 * M_PI;
     while (z_diff(1) < -M_PI) z_diff(1) = z_diff(1) + 2 * M_PI;
     S = S + (weights_(i) * z_diff * (z_diff.transpose()));
   }
   MatrixXd R = MatrixXd(n_z, n_z);
   R << std_radr_ * std_radr_, 0, 0,
        0, std_radphi_ * std_radphi_, 0,
        0, 0, std_radrd_ * std_radrd_;

   S = S + R;

   // Calculating cross correlation matrix
   MatrixXd T = MatrixXd(n_x_, n_z);
   T.fill(0);

   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     VectorXd x_diff = Xsig_pred_.col(i) - x_;
     VectorXd z_diff = Zsig.col(i) - z_pred;
     while (x_diff(3) > M_PI) x_diff(3) = x_diff(3) - 2 * M_PI;
     while (x_diff(3) < -M_PI) x_diff(3) = x_diff(3) + 2 * M_PI;

     while (z_diff(1) > M_PI) z_diff(1) = z_diff(1) - 2 * M_PI;
     while (z_diff(1) < -M_PI) z_diff(1) = z_diff(1) + 2 * M_PI;

     T = T + (weights_(i) * x_diff * (z_diff.transpose()));
   }

   // Calculating Kalman Gain
   MatrixXd K = T * (S.inverse());

   // Calculating updated state mean
   VectorXd z_diff = z_ - z_pred;
   while (z_diff(1) > M_PI) z_diff(1) = z_diff(1) - 2 * M_PI;
   while (z_diff(1) < -M_PI) z_diff(1) = z_diff(1) + 2 * M_PI;
   x_ = x_ +  K * z_diff;

   // Calculating updated state covariance
   P_ = P_ - (K * S * (K.transpose()));
}
