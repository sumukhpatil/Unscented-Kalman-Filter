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
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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

  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
   if (!is_initialized_) {
     MeasurementPackage meas_package;
     x_ << meas_package.raw_measurements_[0],
           meas_package.raw_measurements_[1],
           0,
           meas_package.raw_measurements_[3],
           0;
     P_ << 1, 0, 0, 0, 0,
           0, 1, 0, 0, 0,
           0, 0, 1000, 0, 0,
           0, 0, 0, 1, 0,
           0, 0, 0, 0, 1000;

     long previous_time = meas_package.timestamp_;
     is_initialized_ = true;
   }
   int n_x_ = 5;
   int n_aug_ = n_x_ + 2;
   double lambda_ = 3 - n_aug_;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location.
   * Modify the state vector, x_. Predict sigma points, the state,
   * and the state covariance matrix.
   */

   // Augmented state matrix
   VectorXd x_aug = VectorXd(n_aug_);
   x_aug.fill(0);
   x_aug.head(n_x_) = x_;

   // Process covariance matrix
   MatrixXd Q_ = MatrixXd(2, 2);
   Q_ << std_a_ * std_a_, 0,
         0, std_yawdd_ * std_yawdd_;

   // Augmented covariance matrix
   MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
   P_aug.fill(0);
   P_aug.topLeftCorner(n_x_, n_x_) = P_;
   P_aug.bottomRightCorner(2, 2) = Q_;

   // Augmented sigma points
   MatrixXd X_sig = MatrixXd(n_aug_, 2 * n_aug_ + 1);
   X_sig.fill(0);
   X_sig.col(0) = x_aug;

   // square root matrix
   MatrixXd A = MatrixXd(n_aug_, n_aug_);
   A = P_aug.llt().matrixL();
   for (int i = 0; i < n_aug_; i++) {
     X_sig.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
     X_sig.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
   }

   // Predicting Xsig matrix
   MatrixXd Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
   Xsig_pred_.fill(0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     double px = X_sig(0, i);
     double py = X_sig(1, i);
     double v = X_sig(2, i);
     double psi = X_sig(3, i);
     double psi_d = X_sig(4, i);

     double px_p, py_p;
     double v_p = v;
     double psi_p = psi + psi_d * delta_t;
     double psi_d_p = psi_d;

     if (fabs(psi_d) > 0.001) {
       px_p = px + ((v / psi_d) * (sin(psi + (psi_d * delta_t)) - sin(psi)));
       py_p = py + ((v / psi_d) * (-cos(psi + (psi_d * delta_t)) + cos(psi)));
     } else {
       double px_p = px + (v * cos(psi * delta_t));
       double py_p = py + (v * sin(psi * delta_t));
     }

     px_p = px_p + (0.5 * delta_t * delta_t * cos(psi) * X_sig(5, i));
     py_p = py_p + (0.5 * delta_t * delta_t * sin(psi) * X_sig(5, i));
     v_p = v_p + (delta_t * X_sig(i, 5));
     psi_p = psi_p + (0.5 * delta_t * delta_t * X_sig(6, i));
     psi_d_p = psi_d_p + (delta_t * X_sig(6, i));

     Xsig_pred_(0, i) = px_p;
     Xsig_pred_(1, i) = py_p;
     Xsig_pred_(2, i) = v_p;
     Xsig_pred_(3, i) = psi_p;
     Xsig_pred_(4, i) = psi_d_p;
   }

   // calculating state x
   weights_ = VectorXd(2 * n_aug_ + 1);
   weights_(0) = lambda_ / (lambda_ + n_aug_);
   for (int i = 1; i < 2 * n_aug_ + 1; i++) {
     weights_(i) = 0.5 / (lambda_ + n_aug_);
   }
   double px_m, py_m, v_m, psi_m, psi_d_m = 0;
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     px_m = px_m + (weights_(i) * Xsig_pred_(0, i));
     py_m = py_m + (weights_(i) * Xsig_pred_(1, i));
     v_m = v_m + (weights_(i) * Xsig_pred_(2, i));
     psi_m = psi_m + (weights_(i) * Xsig_pred_(3, i));
     psi_d_m = psi_d_m + (weights_(i) * Xsig_pred_(4, i));
   }
   x_(0) = px_m;
   x_(1) = py_m;
   x_(2) = v_m;
   x_(3) = psi_m;
   x_(4) = psi_d_m;

   // calculating state covariance
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     VectorXd x_diff = VectorXd(n_x_);
     x_diff = Xsig_pred_.col(i) - x_;
     while (x_diff(3) > M_PI) x_diff(3) = x_diff(3) - 2 * M_PI;
     while (x_diff(3) < -M_PI) x_diff(3) = x_diff(3) + 2 * M_PI;
     P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
   }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
   double n_z_radar = 3;
   MatrixXd Z_Sig = MatrixXd(n_z_radar, 2 * n_aug_ + 1);
   VectorXd z_pred = VectorXd(n_z_radar);
   MatrixXd S = MatrixXd(n_z_radar, n_z_radar);
   S.fill(0);

   double _rho, phi, _rho_dot = 0;
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     double px = Xsig_pred_(0, i);
     double py = Xsig_pred_(1, i);
     double v = Xsig_pred_(2, i);
     double psi = Xsig_pred_(3, i);
     double psi_d = Xsig_pred_(4, i);

     _rho = sqrt(px * px + py * py);
     phi = atan2(py, px);
     _rho_dot = ((px * cos(psi) * v) + (py * sin(psi) * v)) / sqrt(px * px + py * py);
     Z_Sig(0, i) = _rho;
     Z_Sig(1, i) = phi;
     Z_Sig(2, i) = _rho_dot;
   }
   z_pred.fill(0);
   double rho_m, phi_m, rho_d_m = 0;
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        rho_m = rho_m + (weights_(i) * Z_Sig(0, i));
        phi_m = phi_m + (weights_(i) * Z_Sig(1, i));
        rho_d_m = rho_d_m + (weights_(i) * Z_Sig(2, i));
   }
   z_pred(0) = rho_m;
   z_pred(1) = phi_m;
   z_pred(2) = rho_d_m;

   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     VectorXd z_diff = VectorXd(n_z_radar);
     z_diff = Z_Sig.col(i) - z_pred;
     while (z_diff(1) > M_PI) z_diff(1) = z_diff(1) - 2 * M_PI;
     while (z_diff(1) < -M_PI) z_diff(1) = z_diff(1) + 2 * M_PI;
     S = S + (weights_(i) * z_diff * z_diff.transpose());
   }
   MatrixXd R = MatrixXd(n_z_radar, n_z_radar);
   R << std_radr_ * std_radr_, 0, 0,
        0, std_radphi_ * std_radphi_, 0,
        0, 0, std_radrd_ * std_radrd_;
   S = S + R;

   MatrixXd T = MatrixXd(n_x_, n_z_radar);
   T.fill(0);
   VectorXd Tx_diff = VectorXd(n_x_);
   VectorXd Tz_diff = VectorXd(n_z_radar);
   for (int i = 0; i < 2 * n_aug_ + 1; i++){
     Tx_diff = Xsig_pred_.col(i) - x_;
     Tz_diff = Z_Sig.col(i) - z_pred;
     while (Tx_diff(3) > M_PI) Tx_diff(3) = Tx_diff(3) - 2 * M_PI;
     while (Tx_diff(3) < -M_PI) Tx_diff(3) = Tx_diff(3) + 2 * M_PI;
     while (Tz_diff(1) > M_PI) Tz_diff(1) = Tz_diff(1) - 2 * M_PI;
     while (Tz_diff(1) < -M_PI) Tz_diff(1) = Tz_diff(1) + 2 * M_PI;
     T = T + (weights_(i) * Tx_diff * Tz_diff.transpose());
   }

   // Calculating Kalman Gain
   MatrixXd K = MatrixXd(n_x_, n_z_radar);
   K.fill(0);
   K = T * S.inverse();
   VectorXd diff_z = VectorXd(n_z_radar);
   diff_z.fill(0);
   diff_z = meas_package.raw_measurements_ - z_pred;
   while (diff_z(1) > M_PI) diff_z(1) = diff_z(1) - 2 * M_PI;
   while (diff_z(1) < -M_PI) diff_z(1) = diff_z(1) - 2 * M_PI;

   // Updated state
   x_ = x_ + K * diff_z;

   // Updated covariance
   P_ = P_ - (K * S * K.transpose());

}
