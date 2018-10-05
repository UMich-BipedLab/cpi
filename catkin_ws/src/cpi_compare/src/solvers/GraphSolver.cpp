/**
 * MIT License
 * Copyright (c) 2018 Kevin Eckenhoff
 * Copyright (c) 2018 Patrick Geneva
 * Copyright (c) 2018 Guoquan Huang
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


#include "GraphSolver.h"

using namespace std;
using namespace gtsam;



/**
 * This will get two pose measurements and will initialize the system to the true orientation and pose
 * The velocity is calculated through the difference between the two poses.
 */
void GraphSolver::addtrue_pose(double timestamp, Eigen::Vector4d q_GtoI, Eigen::Vector3d p_IinG) {

    // Request access to the imu measurements
    std::unique_lock<std::mutex> lock(truth_mutex);

    // Append to our pose vector
    true_times.push_back(timestamp);
    true_qGtoI.push_back(q_GtoI);
    true_pIinG.push_back(p_IinG);

}



/**
 * This function handles new IMU measurements.
 * We just append this to our IMU vectors which will be used to create a preintegrated measurement later
 */
void GraphSolver::addmeasurement_imu(double timestamp, Eigen::Vector3d linacc, Eigen::Vector3d angvel) {

    {
        // Request access to the imu measurements
        std::unique_lock<std::mutex> lock(imu_mutex);

        // Append this new measurement to the array
        imu_times.push_back(timestamp);
        imu_linaccs.push_back(linacc);
        imu_angvel.push_back(angvel);
    }

}


/**
 * This function takes in the calculated plane variables and their covariances.
 * The planes are seen in the LIDAR frame at the given timestep.
 * We first create a preintegrated measurement up to this time, which we can then insert into the graph.
 */
void GraphSolver::addmeasurement_uv(double timestamp, std::vector<uint> leftids, std::vector<Eigen::Vector2d> leftuv,
                                    std::vector<uint> rightids, std::vector<Eigen::Vector2d> rightuv) {

    // Return if we don't actually have any plane measurements
    if(leftids.empty() || rightids.empty())
        return;

    // Return if we don't actually have any IMU measurements
    if(imu_times.size() < 2)
        return;

    // We should try to initialize now
    // Or add the current a new IMU measurement and state!
    if(!systeminitalized) {

        trytoinitalize(timestamp);

        // Return if we have not initialized the system yet
        if(!systeminitalized)
            return;

    } else {

        // Store our IMU messages as they get deleted during propagation
        std::deque<double> imu_timestemp(imu_times);
        std::deque<Eigen::Vector3d> imu_linaccstemp(imu_linaccs);
        std::deque<Eigen::Vector3d> imu_angveltemp(imu_angvel);

        //==========================================================================
        // PREINTEGRATION IMU FACTORS
        //==========================================================================

        // Model 1 CPI
        ImuFactorCPIv1 imuFactorMODEL1 = createimufactor_cpi_v1(timestamp, values_initialMODEL1);
        imu_times = imu_timestemp;
        imu_linaccs = imu_linaccstemp;
        imu_angvel = imu_angveltemp;
        graph_newMODEL1->add(imuFactorMODEL1);
        graphMODEL1->add(imuFactorMODEL1);

        // Model 1 CPI
        ImuFactorCPIv2 imuFactorMODEL2 = createimufactor_cpi_v2(timestamp, values_initialMODEL2);
        imu_times = imu_timestemp;
        imu_linaccs = imu_linaccstemp;
        imu_angvel = imu_angveltemp;
        graph_newMODEL2->add(imuFactorMODEL2);
        graphMODEL2->add(imuFactorMODEL2);

        // Forster discrete preintegration
        //ImuFactorCPIv1 imuFactorFORSTER = createimufactor_discrete(timestamp, values_initialFORSTER);
        imu_times = imu_timestemp;
        imu_linaccs = imu_linaccstemp;
        imu_angvel = imu_angveltemp;
        //graph_newFORSTER->add(imuFactorFORSTER);
        //graphFORSTER->add(imuFactorFORSTER);
        graph_newFORSTER->add(imuFactorMODEL1);
        graphFORSTER->add(imuFactorMODEL1);

        // Forster2 discrete preintegration
        PreintegratedCombinedMeasurements* preint_gtsam;
        CombinedImuFactor imuFactorFORSTER2 = createimufactor_discrete2(timestamp, values_initialFORSTER2, preint_gtsam);
        graph_newFORSTER2->add(imuFactorFORSTER2);
        graphFORSTER2->add(imuFactorFORSTER2);


        // Debug print
        /*cout << "=======================================================" << endl;
        cout << "MODEL 1: " << endl << imuFactorMODEL1 << endl;
        cout << "=======================================================" << endl;
        cout << "MODEL 2: " << endl << imuFactorMODEL2 << endl;
        cout << "=======================================================" << endl;
        cout << "DISCRETE: " << endl << imuFactorFORSTER << endl;
        cout << "=======================================================" << endl;
        */

        //==========================================================================
        // NEW PREDICTED STATE
        //==========================================================================

        // Original models
        JPLNavState newstateMODEL1 = getpredictedstate_v1(imuFactorMODEL1, values_initialMODEL1);
        JPLNavState newstateMODEL2 = getpredictedstate_v2(imuFactorMODEL2, values_initialMODEL2);
        //JPLNavState newstateFORSTER = getpredictedstate_v1(imuFactorFORSTER, values_initialFORSTER);
        JPLNavState newstateFORSTER = getpredictedstate_v1(imuFactorMODEL1, values_initialFORSTER);
        NavState newstateFORSTER2 = getpredictedstate_discrete(preint_gtsam, values_initialFORSTER2);
        imuBias::ConstantBias newbiasFORSTER2 = values_initialFORSTER2.at<imuBias::ConstantBias>(B(ct_state));
	      gtsam::Pose3 newposeFORSTER2 = newstateFORSTER2.pose();


        //==========================================================================
        // Ground Truth Initialization
        //==========================================================================
       
        // If we have ground truth (like from simulation, then use that)
        // NOTE: We use ground truth poses to set initial values
        if(!true_times.empty() && config->useGroundTruthInitValues) {
          JPLQuaternion q_GtoI = true_qGtoI.at(true_qGtoI.size()-1);
          Vector3 p_IinG = true_pIinG.at(true_pIinG.size()-1);
          Bias3 ba = Bias3(0.002899, 0.000054, -0.000197);
          Bias3 bg = Bias3(-0.00000, 0.00002, -0.00002);

          newstateMODEL1.set_q(q_GtoI);
          newstateMODEL1.set_p(p_IinG);
          //newstateMODEL1.set_ba(ba);
          //newstateMODEL1.set_bg(bg);
          newstateMODEL2.set_q(q_GtoI);
          newstateMODEL2.set_p(p_IinG);
          //newstateMODEL2.set_ba(ba);
          //newstateMODEL2.set_bg(bg);
          newstateFORSTER.set_q(q_GtoI);
          newstateFORSTER.set_p(p_IinG);
          //newstateFORSTER.set_ba(ba);
          //newstateFORSTER.set_bg(bg);
          newposeFORSTER2 = gtsam::Pose3(gtsam::Quaternion(q_GtoI(3), q_GtoI(0), q_GtoI(1), q_GtoI(2)), gtsam::Point3(p_IinG));
        }

        // Debug print
        /*cout << "=======================================================" << endl;
        cout << "MODEL 1: " << endl << newstateMODEL1 << endl;
        cout << "=======================================================" << endl;
        cout << "MODEL 2: " << endl << newstateMODEL2 << endl;
        cout << "=======================================================" << endl;
        cout << "DISCRETE: " << endl << newstateFORSTER << endl;
        cout << "=======================================================" << endl;
        cout << "DISCRETE 2: " << endl << newstateFORSTER2 << endl;
        cout << "=======================================================" << endl;*/

       
        // Move node count forward in time
        ct_state++;

        // Append to our node vectors
        values_newMODEL1.insert(X(ct_state), newstateMODEL1);
        values_newMODEL2.insert(X(ct_state), newstateMODEL2);
        values_newFORSTER.insert(X(ct_state), newstateFORSTER);
        values_newFORSTER2.insert(Y(ct_state), newposeFORSTER2);
        values_newFORSTER2.insert(V(ct_state), newstateFORSTER2.v());
        values_newFORSTER2.insert(B(ct_state), newbiasFORSTER2);

        values_initialMODEL1.insert(X(ct_state), newstateMODEL1);
        values_initialMODEL2.insert(X(ct_state), newstateMODEL2);
        values_initialFORSTER.insert(X(ct_state), newstateFORSTER);
        values_initialFORSTER2.insert(Y(ct_state), newposeFORSTER2);
        values_initialFORSTER2.insert(V(ct_state), newstateFORSTER2.v());
        values_initialFORSTER2.insert(B(ct_state), newbiasFORSTER2);

        // Append to our fix lag smoother timestamps
        newTimestampsMODEL1[X(ct_state)] = timestamp;
        newTimestampsMODEL2[X(ct_state)] = timestamp;
        newTimestampsFORSTER[X(ct_state)] = timestamp;
        newTimestampsFORSTER2[Y(ct_state)] = timestamp;

    }

    // Assert our vectors are equal (note will need to remove top one eventually)
    assert(leftids.size() == leftuv.size());
    assert(rightids.size() == rightuv.size());

    // Debug printing
    cout << "[FEAT]: " << leftids.size() << " left features | " << rightids.size() << " right features" << endl;

    // Request access
    std::unique_lock<std::mutex> lock(features_mutex);

    // If we are using inverse depth, then lets call on it
    if(config->useInverseDepth) {
        process_feat_inverse(timestamp,leftids,leftuv,rightids,rightuv);
    } else {
        process_feat_normal(timestamp, leftids, leftuv, rightids, rightuv);
    }

}


/**
 * If we should optimize using the fix lag smoothers
 * We will send in new measurements and nodes, and get back the estimate of the current state
 */
void GraphSolver::optimize() {

    // Return if not initialized
    if(!systeminitalized && ct_state < 2)
        return;

    // Start our timer
    boost::posix_time::ptime t1(boost::posix_time::microsec_clock::local_time());

    // Perform smoothing update
    try {
        smootherBatchMODEL1->update(*graph_newMODEL1, values_newMODEL1, newTimestampsMODEL1);
        values_initialMODEL1 = smootherBatchMODEL1->calculateEstimate();
    } catch(gtsam::IndeterminantLinearSystemException &e) {
        ROS_ERROR("CPI MODEL 1 gtsam indeterminate linear system exception!");
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // Perform smoothing update
    try {
        smootherBatchMODEL2->update(*graph_newMODEL2, values_newMODEL2, newTimestampsMODEL2);
        values_initialMODEL2 = smootherBatchMODEL2->calculateEstimate();
    } catch(gtsam::IndeterminantLinearSystemException &e) {
        ROS_ERROR("CPI MODEL 2 gtsam indeterminate linear system exception!");
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // Perform smoothing update
    try {
        smootherBatchFORSTER->update(*graph_newFORSTER, values_newFORSTER, newTimestampsFORSTER);
        values_initialFORSTER = smootherBatchFORSTER->calculateEstimate();
    } catch(gtsam::IndeterminantLinearSystemException &e) {
        ROS_ERROR("FORSTER gtsam indeterminate linear system exception!");
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // Perform smoothing update
    try {
        smootherBatchFORSTER2->update(*graph_newFORSTER2, values_newFORSTER2, newTimestampsFORSTER2);
        values_initialFORSTER2 = smootherBatchFORSTER2->calculateEstimate();
    } catch(gtsam::IndeterminantLinearSystemException &e) {
        ROS_ERROR("FORSTER2 gtsam indeterminate linear system exception!");
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // Remove the used up nodes
    values_newMODEL1.clear();
    values_newMODEL2.clear();
    values_newFORSTER.clear();
    values_newFORSTER2.clear();

    // Clear used timestamps
    newTimestampsMODEL1.clear();
    newTimestampsMODEL2.clear();
    newTimestampsFORSTER.clear();
    newTimestampsFORSTER2.clear();

    // Remove the used up factors
    graph_newMODEL1->resize(0);
    graph_newMODEL2->resize(0);
    graph_newFORSTER->resize(0);
    graph_newFORSTER2->resize(0);

    // Debug print time
    boost::posix_time::ptime t2(boost::posix_time::microsec_clock::local_time());
    ROS_INFO("[GRAPH]: %d states | %d features | %d edges [%.5f seconds to optimize]",(int)ct_state+1,(int)ct_features,(int)graph_newMODEL1->size(),(t2-t1).total_microseconds()*1e-6);

}


void GraphSolver::optimizeLM() {

    // Return if not initialized
    if(!systeminitalized || ct_state < 587)
        return;

    // Start our timer
    boost::posix_time::ptime t1(boost::posix_time::microsec_clock::local_time());

    // Perform smoothing update
    try {
      cout << "initial error = " << graphMODEL1->error(values_initialMODEL1) << endl;
      LevenbergMarquardtOptimizer optimizer(*graphMODEL1, values_initialMODEL1);
      values_initialMODEL1 = optimizer.optimize();
      cout << "final error = " << graphMODEL1->error(values_initialMODEL1) << endl;
    } catch(gtsam::IndeterminantLinearSystemException &e) {
        ROS_ERROR("CPI MODEL 1 gtsam indeterminate linear system exception!");
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // Perform smoothing update
    try {
      cout << "initial error = " << graphMODEL2->error(values_initialMODEL2) << endl;
      LevenbergMarquardtOptimizer optimizer(*graphMODEL2, values_initialMODEL2);
      values_initialMODEL2 = optimizer.optimize();
      cout << "final error = " << graphMODEL2->error(values_initialMODEL2) << endl;
    } catch(gtsam::IndeterminantLinearSystemException &e) {
        ROS_ERROR("CPI MODEL 2 gtsam indeterminate linear system exception!");
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // Perform smoothing update
    try {
      cout << "initial error = " << graphFORSTER->error(values_initialFORSTER) << endl;
      LevenbergMarquardtOptimizer optimizer(*graphFORSTER, values_initialFORSTER);
      values_initialFORSTER = optimizer.optimize();
      cout << "final error = " << graphFORSTER->error(values_initialFORSTER) << endl;
    } catch(gtsam::IndeterminantLinearSystemException &e) {
        ROS_ERROR("FORSTER gtsam indeterminate linear system exception!");
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // Perform smoothing update
    try {
      cout << "initial error = " << graphFORSTER2->error(values_initialFORSTER2) << endl;
      LevenbergMarquardtOptimizer optimizer(*graphFORSTER2, values_initialFORSTER2);
      values_initialFORSTER2 = optimizer.optimize();
      cout << "final error = " << graphFORSTER2->error(values_initialFORSTER2) << endl;
    } catch(gtsam::IndeterminantLinearSystemException &e) {
        ROS_ERROR("FORSTER2 gtsam indeterminate linear system exception!");
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // Debug print time
    boost::posix_time::ptime t2(boost::posix_time::microsec_clock::local_time());
    ROS_INFO("[GRAPH]: %d states | %d features | %d edges [%.5f seconds to optimize]",(int)ct_state+1,(int)ct_features,(int)graphMODEL1->size(),(t2-t1).total_microseconds()*1e-6);
}


void GraphSolver::optimizeISAM2() {

    // Return if not initialized
    if(!systeminitalized && ct_state < 2)
        return;

    // Start our timer
    boost::posix_time::ptime t1(boost::posix_time::microsec_clock::local_time());

    // Perform smoothing update
    /*try {
        isam2MODEL1->update(*graph_newMODEL1, values_newMODEL1);
        isam2MODEL1->update();
        isam2MODEL1->update();
        values_initialMODEL1 = isam2MODEL1->calculateEstimate();
    } catch(gtsam::IndeterminantLinearSystemException &e) {
        ROS_ERROR("CPI MODEL 1 gtsam indeterminate linear system exception!");
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // Perform smoothing update
    try {
        isam2MODEL2->update(*graph_newMODEL2, values_newMODEL2);
        isam2MODEL2->update();
        isam2MODEL2->update();
        values_initialMODEL2 = isam2MODEL2->calculateEstimate();
    } catch(gtsam::IndeterminantLinearSystemException &e) {
        ROS_ERROR("CPI MODEL 2 gtsam indeterminate linear system exception!");
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // Perform smoothing update
    try {
        isam2FORSTER->update(*graph_newFORSTER, values_newFORSTER);
        isam2FORSTER->update();
        isam2FORSTER->update();
        values_initialFORSTER = isam2FORSTER->calculateEstimate();
    } catch(gtsam::IndeterminantLinearSystemException &e) {
        ROS_ERROR("FORSTER gtsam indeterminate linear system exception!");
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
    }*/
    
    // Perform smoothing update
    try {
      cout << "initial error = " << graphFORSTER2->error(values_initialFORSTER2) << endl;
      isam2FORSTER2->update(*graph_newFORSTER2, values_newFORSTER2);
        isam2FORSTER2->update();
        //isam2FORSTER2->update();
        values_initialFORSTER2 = isam2FORSTER2->calculateEstimate();
      cout << "final error = " << graphFORSTER2->error(values_initialFORSTER2) << endl;
    } catch(gtsam::IndeterminantLinearSystemException &e) {
        ROS_ERROR("FORSTER2 gtsam indeterminate linear system exception!");
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // Remove the used up nodes
    values_newMODEL1.clear();
    values_newMODEL2.clear();
    values_newFORSTER.clear();
    values_newFORSTER2.clear();

    // Clear used timestamps
    newTimestampsMODEL1.clear();
    newTimestampsMODEL2.clear();
    newTimestampsFORSTER.clear();
    newTimestampsFORSTER2.clear();

    // Remove the used up factors
    graph_newMODEL1->resize(0);
    graph_newMODEL2->resize(0);
    graph_newFORSTER->resize(0);
    graph_newFORSTER2->resize(0);

    // Debug print time
    boost::posix_time::ptime t2(boost::posix_time::microsec_clock::local_time());
    ROS_INFO("[GRAPH]: %d states | %d features | %d edges [%.5f seconds to optimize]",(int)ct_state+1,(int)ct_features,(int)graph_newMODEL1->size(),(t2-t1).total_microseconds()*1e-6);
}


/**
 * This will try to take the current IMU vector and initalize
 * If there are enough IMU, we should find the current orientation and biases
 * NOTE: This assumes that we are starting from rest!!
 */
void GraphSolver::trytoinitalize(double timestamp) {

    // If we have already initialized, then just return
    if(systeminitalized)
        return;

    // Wait for enough IMU readings
    if(imu_times.size() < (size_t)config->imuWait)
        return;

    //==========================================================================
    // INITIALIZE!
    //==========================================================================

    // Request access to the imu measurements
    std::unique_lock<std::mutex> lock1(imu_mutex);
    std::unique_lock<std::mutex> lock2(truth_mutex);

    // Sum up our current accelerations and velocities
    Eigen::Vector3d linsum = Eigen::Vector3d::Zero();
    Eigen::Vector3d angsum = Eigen::Vector3d::Zero();
    for(size_t i=0; i<imu_times.size(); i++) {
        linsum += imu_linaccs.at(i);
        angsum += imu_angvel.at(i);
    }

    // Calculate the mean of the linear acceleration and angular velocity
    Eigen::Vector3d linavg = Eigen::Vector3d::Zero();
    Eigen::Vector3d angavg = Eigen::Vector3d::Zero();
    linavg = linsum/imu_times.size();
    angavg = angsum/imu_times.size();

    // Get z axis, which alines with -g (z_in_G=0,0,1)
    Eigen::Vector3d z_axis = linavg/linavg.norm();

    // Create an x_axis
    Eigen::Vector3d e_1(1,0,0);

    // Make x_axis perpendicular to z
    Eigen::Vector3d x_axis = e_1-z_axis*z_axis.transpose()*e_1;
    x_axis= x_axis/x_axis.norm();

    // Get z from the cross product of these two
    Eigen::Matrix<double,3,1> y_axis = skew_x(z_axis)*x_axis;

    // From these axes get rotation
    Eigen::Matrix<double,3,3> Ro;
    Ro.block(0,0,3,1) = x_axis;
    Ro.block(0,1,3,1) = y_axis;
    Ro.block(0,2,3,1) = z_axis;

    // Create our state variables
    JPLQuaternion q_GtoI = rot_2_quat(Ro);
    Vector3 v_IinG(0,0,0);
    Vector3 p_IinG(0,0,0);

    // Set our biases equal to our noise (subtract our gravity from accelerometer bias)
    Bias3 bg = angavg;
    Bias3 ba = linavg - quat_2_Rot(q_GtoI)*config->gravity;

    // If we have ground truth (like from simulation, then use that)
    // NOTE: We still need to start from zero here
    if(!true_times.empty() && config->useGroundTruthInit) {
        q_GtoI = true_qGtoI.at(true_qGtoI.size()-1);
        p_IinG = true_pIinG.at(true_pIinG.size()-1);
        //ba = Bias3(0.002899, 0.000054, -0.000197);
        //bg = Bias3(0.00000, 0.00002, -0.00002);
        ba = Bias3(0.0, 0.0, 0.0);
        bg = Bias3(0.0, 0.0, 0.0);
    }

    // Since this is the true pose, lets just assign a small error
    // Note: this isn't ideal, but works fine since we are using true pose for simulation
    Eigen::Matrix<double,15,15> cov = 1e-8*Eigen::Matrix<double,15,15>::Identity();

    //==========================================================================
    // CREATE PRIOR FACTORS AND INITALIZE GRAPHS
    //==========================================================================

    // Create our prior factor and add it to our graph
    JPLNavStatePrior priorfactor(X(ct_state), cov, q_GtoI, bg, v_IinG, ba, p_IinG);
    graph_newMODEL1->add(priorfactor);
    graph_newMODEL2->add(priorfactor);
    graph_newFORSTER->add(priorfactor);
    graphMODEL1->add(priorfactor);
    graphMODEL2->add(priorfactor);
    graphFORSTER->add(priorfactor);

    // Create gtsam prior factor and add it to FORSTER2 graph
    gtsam::Pose3 prior_pose = gtsam::Pose3(gtsam::Quaternion(q_GtoI(3), q_GtoI(0), q_GtoI(1), q_GtoI(2)), gtsam::Point3(p_IinG));
    imuBias::ConstantBias prior_imu_bias;
    
    noiseModel::Diagonal::shared_ptr pose_noise_model = noiseModel::Diagonal::Sigmas((Vector(6) << 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1).finished());
    noiseModel::Diagonal::shared_ptr velocity_noise_model = noiseModel::Isotropic::Sigma(3, 1e-2);
    noiseModel::Diagonal::shared_ptr bias_noise_model = noiseModel::Isotropic::Sigma(6, config->sigma_wa);

    graph_newFORSTER2->add(PriorFactor<Pose3>(Y(ct_state), prior_pose, pose_noise_model));
    graph_newFORSTER2->add(PriorFactor<Vector3>(V(ct_state), v_IinG, velocity_noise_model));
    graph_newFORSTER2->add(PriorFactor<imuBias::ConstantBias>(B(ct_state), prior_imu_bias, bias_noise_model));
    graphFORSTER2->add(PriorFactor<Pose3>(Y(ct_state), prior_pose, pose_noise_model));
    graphFORSTER2->add(PriorFactor<Vector3>(V(ct_state), v_IinG, velocity_noise_model));
    graphFORSTER2->add(PriorFactor<imuBias::ConstantBias>(B(ct_state), prior_imu_bias, bias_noise_model));


    // Add our initial state
    values_newMODEL1.insert(X(ct_state), JPLNavState(q_GtoI, bg, v_IinG, ba, p_IinG));
    values_newMODEL2.insert(X(ct_state), JPLNavState(q_GtoI, bg, v_IinG, ba, p_IinG));
    values_newFORSTER.insert(X(ct_state), JPLNavState(q_GtoI, bg, v_IinG, ba, p_IinG));
    values_initialMODEL1.insert(X(ct_state), JPLNavState(q_GtoI, bg, v_IinG, ba, p_IinG));
    values_initialMODEL2.insert(X(ct_state), JPLNavState(q_GtoI, bg, v_IinG, ba, p_IinG));
    values_initialFORSTER.insert(X(ct_state), JPLNavState(q_GtoI, bg, v_IinG, ba, p_IinG));

    // Add initial state to FORSTER2 model
    values_newFORSTER2.insert(Y(ct_state), prior_pose);
    values_newFORSTER2.insert(V(ct_state), v_IinG);
    values_newFORSTER2.insert(B(ct_state), prior_imu_bias);
    
    values_initialFORSTER2.insert(Y(ct_state), prior_pose);
    values_initialFORSTER2.insert(V(ct_state), v_IinG);
    values_initialFORSTER2.insert(B(ct_state), prior_imu_bias);

    // Append to our fix lag smoother timestamps
    newTimestampsMODEL1[X(ct_state)] = timestamp;
    newTimestampsMODEL2[X(ct_state)] = timestamp;
    newTimestampsFORSTER[X(ct_state)] = timestamp;
    newTimestampsFORSTER2[Y(ct_state)] = timestamp;

    // Clear all old imu messages (keep the last two)
    imu_times.erase(imu_times.begin(), imu_times.end()-1);
    imu_linaccs.erase(imu_linaccs.begin(), imu_linaccs.end()-1);
    imu_angvel.erase(imu_angvel.begin(), imu_angvel.end()-1);

    // Set our initialized to true!
    systeminitalized = true;

    // Debug info
    ROS_INFO("\033[0;32m[INIT]: orientation = %.4f, %.4f, %.4f, %.4f\033[0m",q_GtoI(0),q_GtoI(1),q_GtoI(2),q_GtoI(3));
    ROS_INFO("\033[0;32m[INIT]: bias gyro = %.4f, %.4f, %.4f\033[0m",bg(0),bg(1),bg(2));
    ROS_INFO("\033[0;32m[INIT]: velocity = %.4f, %.4f, %.4f\033[0m",v_IinG(0),v_IinG(1),v_IinG(2));
    ROS_INFO("\033[0;32m[INIT]: bias accel = %.4f, %.4f, %.4f\033[0m",ba(0),ba(1),ba(2));
    ROS_INFO("\033[0;32m[INIT]: position = %.4f, %.4f, %.4f\033[0m",p_IinG(0),p_IinG(1),p_IinG(2));

}











