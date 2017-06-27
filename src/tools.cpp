#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  //**
  //TODO:
  //  * Calculate the RMSE here.
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if(estimations.size() != ground_truth.size()
  		|| estimations.size() ==0){

  	cout << "Invalid Estimations or Ground Truth Data.." << endl;
  	return rmse;
  }

  //compute and accumulate squared residuals
  for(unsigned int i=0; i<estimations.size(); ++i){

  	VectorXd residual = estimations[i] - ground_truth[i];

  	//coefficient-wise multiplication
  	residual = residual.array() * residual.array();
  	rmse += residual;
  }

  //Calculate the Mean
  rmse = rmse/estimations.size();

  //Calculate the Sq Root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
	MatrixXd Hj(3,4);
  Hj << 0,0,0,0,
      0,0,0,0,
      0,0,0,0;

	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

    float px_sq = px*px;
    float py_sq = py*py;
    float pxy_sq = px_sq + py_sq;
    float pxy_sqrt = sqrt(pxy_sq);
    float pxy_sq_ex = pow(pxy_sq, 3/2);
	
	//check division by zero
	//compute the Jacobian matrix
	if (pxy_sq!=0 && pxy_sqrt!=0 && pxy_sq_ex!=0) 
        Hj << px/pxy_sqrt, py/pxy_sqrt, 0, 0,
              -py/pxy_sq, px/pxy_sq, 0, 0,
               py*(vx*py - vy*px)/pxy_sq_ex, px*(vy*px - vx*py)/pxy_sq_ex, px/pxy_sqrt, py/pxy_sqrt;
    
	return Hj;  
}
