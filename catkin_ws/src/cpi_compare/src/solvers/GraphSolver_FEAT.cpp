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
 * OUT OF OR IN CONNECTION iWITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


#include "GraphSolver.h"




/**
 * This will add normal UV factors with the feature being represented as a 3D position
 * This might have some instabilities depending on the feature layout in the simulation
 */
void GraphSolver::process_feat_normal(double timestamp, std::vector<uint> leftids, std::vector<Eigen::Vector2d> leftuv,
                                      std::vector<uint> rightids, std::vector<Eigen::Vector2d> rightuv) {

    //==============================================================================
    // Loop through LEFT features
    for(size_t i=0; i<leftids.size(); i++) {
        // Check to see if it is already in the graph
        if(measurement_lookup.find(leftids.at(i)) != measurement_lookup.end() && (ct_state-measurement_state_lookup[leftids.at(i)]) < (size_t)config->uvWindowSize) {
            Eigen::Matrix<double,2,2> sqrtQ = config->sigma_camera_sq*Eigen::Matrix<double,2,2>::Identity();
            JPLImageUVFactor factor(X(ct_state),F(measurement_lookup[leftids.at(i)]),sqrtQ,leftuv.at(i),config->R_C0toI,config->p_IinC0);
            graphMODEL1->add(factor);
            graphMODEL2->add(factor);
            graphFORSTER->add(factor);
            graph_newMODEL1->add(factor);
            graph_newMODEL2->add(factor);
            graph_newFORSTER->add(factor);

            // Insert measurements to a smart factor
            //measurement_smart_lookup_left[leftids.at(i)]->add(Point2(leftuv.at(i)), Y(ct_state));
            continue;
        }
        // Next check to see if this feature is in our queue
        // We know this is in our queue, so lets add this new measurement to it
        if(measurement_queue.find(leftids.at(i)) != measurement_queue.end()) {
            measurement_queue[leftids.at(i)].leftids.push_back(leftids.at(i));
            measurement_queue[leftids.at(i)].leftuv.push_back(leftuv.at(i));
            measurement_queue[leftids.at(i)].leftstateids.push_back(ct_state);
            continue;
        }
        // Else this is a new feature, so lets just add it as a new feature
        // Create our new feature object
        feature featnew;
        featnew.leftids.push_back(leftids.at(i));
        featnew.leftuv.push_back(leftuv.at(i));
        featnew.leftstateids.push_back(ct_state);
        // Push into our queue
        measurement_queue[leftids.at(i)] = featnew;
    }

    //==============================================================================
    // Loop through RIGHT features
    for(size_t i=0; i<rightids.size(); i++) {
        // Check to see if it is already in the graph
        if(measurement_lookup.find(rightids.at(i)) != measurement_lookup.end() && (ct_state-measurement_state_lookup[rightids.at(i)]) < (size_t)config->uvWindowSize) {
            Eigen::Matrix<double,2,2> sqrtQ = config->sigma_camera_sq*Eigen::Matrix<double,2,2>::Identity();
            JPLImageUVFactor factor(X(ct_state),F(measurement_lookup[rightids.at(i)]),sqrtQ,rightuv.at(i),config->R_C1toI,config->p_IinC1);
            graphMODEL1->add(factor);
            graphMODEL2->add(factor);
            graphFORSTER->add(factor);
            graph_newMODEL1->add(factor);
            graph_newMODEL2->add(factor);
            graph_newFORSTER->add(factor);

            // Insert measurements to a smart factor
           //measurement_sma  rt_lookup_right[rightids.at(i)]->add(Point2(rightuv.at(i)), Y(ct_state));
            auto pos = std::find(leftids.begin(), leftids.end(), rightids.at(i));
            if (pos != leftids.end()) { // We know that the same landmark is seen by both left and right cameras
              Cal3_S2Stereo::shared_ptr K(new Cal3_S2Stereo(1, 1, 0, 0, 0, 0.1));
              measurement_smart_lookup_stereo[rightids.at(i)]->add(StereoPoint2(leftuv.at(pos-leftids.begin())[0], rightuv.at(i)[0], rightuv.at(i)[1]), Y(ct_state), K);
              std::cout << measurement_smart_lookup_stereo[rightids.at(i)]->isValid() << std::endl;
              measurement_smart_lookup_stereo[rightids.at(i)]->print();

            }
           continue;
        }
        // Check to see if this feature is in our queue
        // We know this is in our queue, so lets add this new measurement to it
        if(measurement_queue.find(rightids.at(i)) != measurement_queue.end()) {
            measurement_queue[rightids.at(i)].rightids.push_back(rightids.at(i));
            measurement_queue[rightids.at(i)].rightuv.push_back(rightuv.at(i));
            measurement_queue[rightids.at(i)].rightstateids.push_back(ct_state);
            continue;
        }
        // Else this is a new feature, so lets just add it as a new feature
        // Create our new feature object
        feature featnew;
        featnew.rightids.push_back(rightids.at(i));
        featnew.rightuv.push_back(rightuv.at(i));
        featnew.rightstateids.push_back(ct_state);
        // Push into our queue
        measurement_queue[rightids.at(i)] = featnew;
    }


    //==============================================================================
    // Track which features should be removed from the queue
    std::vector<int> toremove;
    // Create our initializer class
    FeatureInitializer initializerMODEL1(config, values_initialMODEL1, ct_state);
    FeatureInitializer initializerMODEL2(config, values_initialMODEL2, ct_state);
    FeatureInitializer initializerFORSTER(config, values_initialFORSTER, ct_state);
    FeatureInitializer initializerFORSTER2(config, values_initialFORSTER2, ct_state);
    
    int ct_successes = 0;
    int ct_failures = 0;
    // Lastly lets add all the features that have reached the
    for(auto& measurement: measurement_queue) {

        // Boolean logic statements, a correct queue feature does not have the min pose requirement
        // A queue feature should also have a non-empty state ID that is equal to the current state (i.e. is being actively tracked)
        bool leftminposes = measurement.second.leftuv.size() < (size_t)config->minPoseFeatureInit;
        bool rightminposes = measurement.second.rightuv.size() < (size_t)config->minPoseFeatureInit;
        bool leftlost = measurement.second.leftstateids.empty() || measurement.second.leftstateids.at(measurement.second.leftstateids.size()-1) != ct_state;
        bool rightlost = measurement.second.rightstateids.empty() || measurement.second.rightstateids.at(measurement.second.rightstateids.size()-1) != ct_state;

        // If the feature has lost track, then we should remove it if doesn't have the needed UV measurement size
        if(leftlost && rightlost && leftminposes && rightminposes) {
            toremove.push_back(measurement.first);
            //ROS_ERROR("Removing feature #%d as it has lost tracking (%d leftuv, %d rightuv)",measurement.first,(int)measurement.second.leftuv.size(),(int)measurement.second.rightuv.size());
            continue;
        }

        // Remove if less then half right features as there are left ones
        if(measurement.second.rightuv.size() < 0.5*measurement.second.leftuv.size()) {
            toremove.push_back(measurement.first);
            //ROS_ERROR("Removing feature #%d as not enough right features (%d leftuv, %d rightuv)",measurement.first,(int)measurement.second.leftuv.size(),(int)measurement.second.rightuv.size());
            continue;
        }

        // Skip if it has not reached max size, but is still being actively tracked
        if(leftminposes && rightminposes)
            continue;

        // We are good, either the left or right has enough poses to initialize
        // Store our key so that it will be removed from the queue
        toremove.push_back(measurement.first);

        // Check that all IDs are the same
        if(!measurement.second.leftids.empty()) assert(std::equal(measurement.second.leftids.begin()+1, measurement.second.leftids.end(), measurement.second.leftids.begin()));
        if(!measurement.second.rightids.empty()) assert(std::equal(measurement.second.rightids.begin()+1, measurement.second.rightids.end(), measurement.second.rightids.begin()));

        // Initialize the 3D position of the feature
        bool successMODEL1 = initializerMODEL1.initialize_feature(measurement.second);
        std::pair<int, feature> measurementMODEL1 = measurement;
        bool successMODEL2 = initializerMODEL2.initialize_feature(measurement.second);
        std::pair<int, feature> measurementMODEL2 = measurement;
        bool successFORSTER = initializerFORSTER.initialize_feature(measurement.second);
        std::pair<int, feature> measurementFORSTER = measurement;
        bool successFORSTER2 = initializerFORSTER2.initialize_feature(measurement.second);
        std::pair<int, feature> measurementFORSTER2 = measurement;

        // If not successful skip this feature
        if(!successMODEL1 || !successMODEL2 || !successFORSTER || !successFORSTER2) {
            ct_failures++;
            continue;
        }

        // Ensure we have at least one left feature (we pick this to be our anchor)
        if(measurementMODEL1.second.leftuv.empty() || measurementMODEL1.second.rightuv.empty())
            continue;

        // Move feature ID forward in time, and add to our lookup data structure
        // Note: we set our anchor to be the first left camera state
        ct_features++;
        measurement_lookup[measurement.first] = ct_features;
        measurement_state_lookup[measurement.first] = ct_state;

        // Lets add the new feature to graph
        // NOTE: normal 3d feature since we are NOT using inverse depth
        values_newMODEL1.insert(F(ct_features), gtsam::Point3(measurementMODEL1.second.pos_FinG));
        values_newMODEL2.insert(F(ct_features), gtsam::Point3(measurementMODEL2.second.pos_FinG));
        values_newFORSTER.insert(F(ct_features), gtsam::Point3(measurementFORSTER.second.pos_FinG));
        values_newFORSTER2.insert(F(ct_features), gtsam::Point3(measurementFORSTER2.second.pos_FinG));
        values_initialMODEL1.insert(F(ct_features), gtsam::Point3(measurementMODEL1.second.pos_FinG));
        values_initialMODEL2.insert(F(ct_features), gtsam::Point3(measurementMODEL2.second.pos_FinG));
        values_initialFORSTER.insert(F(ct_features), gtsam::Point3(measurementFORSTER.second.pos_FinG));
        
        // Append to our fix lag smoother timestamps
        newTimestampsMODEL1[F(ct_features)] = timestamp;
        newTimestampsMODEL2[F(ct_features)] = timestamp;
        newTimestampsFORSTER[F(ct_features)] = timestamp;
        newTimestampsFORSTER2[F(ct_features)] = timestamp;

        // create a smart stereo factor for the new feature
        noiseModel::Isotropic::shared_ptr measurementNoise = noiseModel::Isotropic::Sigma(3, 0.1);
        Cal3_S2Stereo::shared_ptr K(new Cal3_S2Stereo(1, 1, 0, 0, 0, 0.1));
        SmartStereoFactor::shared_ptr smartfactor(new SmartStereoFactor(measurementNoise));
        measurement_smart_lookup_stereo[measurement.first] = smartfactor;


        // Next lets add all LEFT factors (all graphs have the same measurements!!!)
        for(size_t j=0; j<measurementMODEL1.second.leftuv.size(); j++) {
            Eigen::Matrix<double,2,2> sqrtQ = config->sigma_camera_sq*Eigen::Matrix<double,2,2>::Identity();
            JPLImageUVFactor factor(X(measurementMODEL1.second.leftstateids.at(j)),F(ct_features),sqrtQ,measurementMODEL1.second.leftuv.at(j),config->R_C0toI,config->p_IinC0);
            graphMODEL1->add(factor);
            graphMODEL2->add(factor);
            graphFORSTER->add(factor);
            graph_newMODEL1->add(factor);
            graph_newMODEL2->add(factor);
            graph_newFORSTER->add(factor);
            //std::cout << "left: "<< measurementMODEL1.second.leftstateids.at(j) << " ";
            
            // Insert measurements to a smart factor
            //smartfactor_left->add(Point2(measurementMODEL1.second.leftuv.at(j)), Y(measurementMODEL1.second.leftstateids.at(j)));
        }
        //std::cout << std::endl;

        // Add smart factor to FORSTER2 model
        //graphFORSTER2->push_back(smartfactor_left);
        //graph_newFORSTER2->push_back(smartfactor_left);

        // create a smart factor for the new feature
        //Pose3 sensor_P_body = Pose3(Rot3(config->R_C1toI), Point3(config->p_IinC1));
        //SmartFactor::shared_ptr smartfactor_right(new SmartFactor(measurementNoise, K, sensor_P_body.inverse()));
        //measurement_smart_lookup_right[measurement.first] = smartfactor_right;
        
        // Next lets add all RIGHT factors (all graphs have the same measurements!!!)
        for(size_t j=0; j<measurementMODEL1.second.rightuv.size(); j++) {
            Eigen::Matrix<double,2,2> sqrtQ = config->sigma_camera_sq*Eigen::Matrix<double,2,2>::Identity();
            JPLImageUVFactor factor(X(measurementMODEL1.second.rightstateids.at(j)),F(ct_features),sqrtQ,measurementMODEL1.second.rightuv.at(j),config->R_C1toI,config->p_IinC1);
            graphMODEL1->add(factor);
            graphMODEL2->add(factor);
            graphFORSTER->add(factor);
            graph_newMODEL1->add(factor);
            graph_newMODEL2->add(factor);
            graph_newFORSTER->add(factor);
            //std::cout << "right: " << measurementMODEL1.second.rightstateids.at(j) << " ";

            // Insert measurements to a smart factor
            //smartfactor_right->add(Point2(measurementMODEL1.second.rightuv.at(j)), Y(measurementMODEL1.second.rightstateids.at(j)));
            auto pos = std::find(measurementMODEL1.second.leftstateids.begin(), measurementMODEL1.second.leftstateids.end(), measurementMODEL1.second.rightstateids.at(j));
            if (pos != measurementMODEL1.second.leftstateids.end()) {
              //smartfactor->add(StereoPoint2())
              //std::cout<< " left meas: " << measurementMODEL1.second.leftstateids.at(pos-measurementMODEL1.second.leftstateids.begin()) << " " << measurementMODEL1.second.leftuv.at(pos-measurementMODEL1.second.leftstateids.begin()) << std::endl;
              //std::cout<< "right meas: " << measurementMODEL1.second.rightstateids.at(j) << " " << measurementMODEL1.second.rightuv.at(j) << std::endl;
              smartfactor->add(StereoPoint2(measurementMODEL1.second.leftuv.at(pos-measurementMODEL1.second.leftstateids.begin())[0], measurementMODEL1.second.rightuv.at(j)[0], measurementMODEL1.second.rightuv.at(j)[1]), Y(measurementMODEL1.second.rightstateids.at(j)), K);
              //smartfactor->print();
              //std::cout << smartfactor->isValid() << std::endl;
            }
        
        }
        //std::cout << std::endl;

        // Add smart factor to FORSTER2 model
        //graphFORSTER2->push_back(smartfactor);
        //graph_newFORSTER2->push_back(smartfactor);

        // Record our success
        ct_successes++;

    }

    // Debug info
    if(ct_successes+ct_failures > 0) {
        ROS_WARN("[FEAT]: %d of %d features successfully initialized (non-inverse depth)", ct_successes, ct_successes + ct_failures);
    }


    //==============================================================================
    // Remove the "used up" features from the queue
    for (auto&& key : toremove) {
        measurement_queue.erase(key);
    }



}






/**
 * This will add all factors, such that we are using inverse depth for features!
 * Each inverse feature has a anchor node, and then all measurements will relate to this anchor
 */
void GraphSolver::process_feat_inverse(double timestamp, std::vector<uint> leftids, std::vector<Eigen::Vector2d> leftuv,
                                       std::vector<uint> rightids, std::vector<Eigen::Vector2d> rightuv) {


    //==============================================================================
    // Loop through LEFT features
    for(size_t i=0; i<leftids.size(); i++) {
        // Check to see if it is already in the graph
        if(measurement_lookup.find(leftids.at(i)) != measurement_lookup.end() && (ct_state-measurement_state_lookup[leftids.at(i)]) < (size_t)config->uvWindowSize) {
            Eigen::Matrix<double,2,2> sqrtQ = config->sigma_camera_sq*Eigen::Matrix<double,2,2>::Identity();
            InvUVFactor factor(X(ct_state),X(measurement_anchor_lookup[measurement_lookup[leftids.at(i)]]),F(measurement_lookup[leftids.at(i)]),sqrtQ,leftuv.at(i),config->R_C0toI,config->p_IinC0,config->R_C0toI,config->p_IinC0);
            graphMODEL1->add(factor);
            graphMODEL2->add(factor);
            graphFORSTER->add(factor);
            graph_newMODEL1->add(factor);
            graph_newMODEL2->add(factor);
            graph_newFORSTER->add(factor);
            continue;
        }
        // Next check to see if this feature is in our queue
        // We know this is in our queue, so lets add this new measurement to it
        if(measurement_queue.find(leftids.at(i)) != measurement_queue.end()) {
            measurement_queue[leftids.at(i)].leftids.push_back(leftids.at(i));
            measurement_queue[leftids.at(i)].leftuv.push_back(leftuv.at(i));
            measurement_queue[leftids.at(i)].leftstateids.push_back(ct_state);
            continue;
        }
        // Else this is a new feature, so lets just add it as a new feature
        // Create our new feature object
        feature featnew;
        featnew.leftids.push_back(leftids.at(i));
        featnew.leftuv.push_back(leftuv.at(i));
        featnew.leftstateids.push_back(ct_state);
        // Push into our queue
        measurement_queue[leftids.at(i)] = featnew;
    }

    //==============================================================================
    // Loop through RIGHT features
    for(size_t i=0; i<rightids.size(); i++) {
        // Check to see if it is already in the graph
        if(measurement_lookup.find(rightids.at(i)) != measurement_lookup.end() && (ct_state-measurement_state_lookup[rightids.at(i)]) < (size_t)config->uvWindowSize) {
            Eigen::Matrix<double,2,2> sqrtQ = config->sigma_camera_sq*Eigen::Matrix<double,2,2>::Identity();
            InvUVFactor factor(X(ct_state),X(measurement_anchor_lookup[measurement_lookup[rightids.at(i)]]),F(measurement_lookup[rightids.at(i)]),sqrtQ,rightuv.at(i),config->R_C1toI,config->p_IinC1,config->R_C0toI,config->p_IinC0);
            graphMODEL1->add(factor);
            graphMODEL2->add(factor);
            graphFORSTER->add(factor);
            graph_newMODEL1->add(factor);
            graph_newMODEL2->add(factor);
            graph_newFORSTER->add(factor);
            continue;
        }
        // Check to see if this feature is in our queue
        // We know this is in our queue, so lets add this new measurement to it
        if(measurement_queue.find(rightids.at(i)) != measurement_queue.end()) {
            measurement_queue[rightids.at(i)].rightids.push_back(rightids.at(i));
            measurement_queue[rightids.at(i)].rightuv.push_back(rightuv.at(i));
            measurement_queue[rightids.at(i)].rightstateids.push_back(ct_state);
            continue;
        }
        // Else this is a new feature, so lets just add it as a new feature
        // Create our new feature object
        feature featnew;
        featnew.rightids.push_back(rightids.at(i));
        featnew.rightuv.push_back(rightuv.at(i));
        featnew.rightstateids.push_back(ct_state);
        // Push into our queue
        measurement_queue[rightids.at(i)] = featnew;
    }


    //==============================================================================
    // Track which features should be removed from the queue
    std::vector<int> toremove;
    // Create our initializer class
    FeatureInitializer initializerMODEL1(config, values_initialMODEL1, ct_state);
    FeatureInitializer initializerMODEL2(config, values_initialMODEL2, ct_state);
    FeatureInitializer initializerFORSTER(config, values_initialFORSTER, ct_state);
    int ct_successes = 0;
    int ct_failures = 0;
    // Lastly lets add all the features that have reached the
    for(auto& measurement: measurement_queue) {

        // Boolean logic statements, a correct queue feature does not have the min pose requirement
        // A queue feature should also have a non-empty state ID that is equal to the current state (i.e. is being actively tracked)
        bool leftminposes = measurement.second.leftuv.size() < (size_t)config->minPoseFeatureInit;
        bool rightminposes = measurement.second.rightuv.size() < (size_t)config->minPoseFeatureInit;
        bool leftlost = measurement.second.leftstateids.empty() || measurement.second.leftstateids.at(measurement.second.leftstateids.size()-1) != ct_state;
        bool rightlost = measurement.second.rightstateids.empty() || measurement.second.rightstateids.at(measurement.second.rightstateids.size()-1) != ct_state;

        // If the feature has lost track, then we should remove it if doesn't have the needed UV measurement size
        if(leftlost && rightlost && leftminposes && rightminposes) {
            toremove.push_back(measurement.first);
            //ROS_ERROR("Removing feature #%d as it has lost tracking (%d leftuv, %d rightuv)",measurement.first,(int)measurement.second.leftuv.size(),(int)measurement.second.rightuv.size());
            continue;
        }

        // Remove if less then half right features as there are left ones
        if(measurement.second.rightuv.size() < 0.5*measurement.second.leftuv.size()) {
            toremove.push_back(measurement.first);
            //ROS_ERROR("Removing feature #%d as not enough right features (%d leftuv, %d rightuv)",measurement.first,(int)measurement.second.leftuv.size(),(int)measurement.second.rightuv.size());
            continue;
        }

        // Skip if it has not reached max size, but is still being actively tracked
        if(leftminposes && rightminposes)
            continue;

        // We are good, either the left or right has enough poses to initialize
        // Store our key so that it will be removed from the queue
        toremove.push_back(measurement.first);

        // Check that all IDs are the same
        if(!measurement.second.leftids.empty()) assert(std::equal(measurement.second.leftids.begin()+1, measurement.second.leftids.end(), measurement.second.leftids.begin()));
        if(!measurement.second.rightids.empty()) assert(std::equal(measurement.second.rightids.begin()+1, measurement.second.rightids.end(), measurement.second.rightids.begin()));

        // Initialize the 3D position of the feature
        bool successMODEL1 = initializerMODEL1.initialize_feature(measurement.second);
        std::pair<int, feature> measurementMODEL1 = measurement;
        bool successMODEL2 = initializerMODEL2.initialize_feature(measurement.second);
        std::pair<int, feature> measurementMODEL2 = measurement;
        bool successFORSTER = initializerFORSTER.initialize_feature(measurement.second);
        std::pair<int, feature> measurementFORSTER = measurement;

        // If not successful skip this feature
        if(!successMODEL1 || !successMODEL2 || !successFORSTER) {
          ct_failures++;
            continue;
        }

        // Ensure we have at least one left feature (we pick this to be our anchor)
        if(measurementMODEL1.second.leftuv.empty() || measurementMODEL1.second.rightuv.empty())
            continue;

        // Ensure we have the same anchor frame in the left and right!
        if(measurementMODEL1.second.leftstateids.at(0) != measurementMODEL1.second.rightstateids.at(0))
            continue;

        // Move feature ID forward in time, and add to our lookup data structure
        // Note: we set our anchor to be the first left camera state
        ct_features++;
        measurement_lookup[measurement.first] = ct_features;
        measurement_state_lookup[measurement.first] = ct_state;
        measurement_anchor_lookup[ct_features] = measurementMODEL1.second.leftstateids.at(0);

        // Lets add the new feature to graph
        // NOTE: we use the inverse depth value in the ANCHOR frame as the 3x1 vec
        values_newMODEL1.insert(F(ct_features), gtsam::Point3(measurementMODEL1.second.pos_FinA_inv));
        values_newMODEL2.insert(F(ct_features), gtsam::Point3(measurementMODEL2.second.pos_FinA_inv));
        values_newFORSTER.insert(F(ct_features), gtsam::Point3(measurementFORSTER.second.pos_FinA_inv));
        values_initialMODEL1.insert(F(ct_features), gtsam::Point3(measurementMODEL1.second.pos_FinA_inv));
        values_initialMODEL2.insert(F(ct_features), gtsam::Point3(measurementMODEL2.second.pos_FinA_inv));
        values_initialFORSTER.insert(F(ct_features), gtsam::Point3(measurementFORSTER.second.pos_FinA_inv));

        // Append to our fix lag smoother timestamps
        newTimestampsMODEL1[F(ct_features)] = timestamp;
        newTimestampsMODEL2[F(ct_features)] = timestamp;
        newTimestampsFORSTER[F(ct_features)] = timestamp;

        // Insert our anchor measurements, note they are only a function of the feature!
        Eigen::Matrix<double,2,2> sqrtQA = config->sigma_camera_sq*Eigen::Matrix<double,2,2>::Identity();
        InvAnchorFactor factorL(F(ct_features),sqrtQA,measurementMODEL1.second.leftuv.at(0),config->R_C0toI,config->p_IinC0,config->R_C0toI,config->p_IinC0);
        graphMODEL1->add(factorL);
        graphMODEL2->add(factorL);
        graphFORSTER->add(factorL);
        graph_newMODEL1->add(factorL);
        graph_newMODEL2->add(factorL);
        graph_newFORSTER->add(factorL);
        InvAnchorFactor factorR(F(ct_features),sqrtQA,measurementMODEL1.second.rightuv.at(0),config->R_C1toI,config->p_IinC1,config->R_C0toI,config->p_IinC0);
        graphMODEL1->add(factorR);
        graphMODEL2->add(factorR);
        graphFORSTER->add(factorR);
        graph_newMODEL1->add(factorR);
        graph_newMODEL2->add(factorR);
        graph_newFORSTER->add(factorR);

        // Next lets add all LEFT factors (all graphs have the same measurements!!!)
        for(size_t j=1; j<measurementMODEL1.second.leftuv.size(); j++) {
            Eigen::Matrix<double,2,2> sqrtQ = config->sigma_camera_sq*Eigen::Matrix<double,2,2>::Identity();
            InvUVFactor factor(X(measurementMODEL1.second.leftstateids.at(j)),X(measurement_anchor_lookup[ct_features]),F(ct_features),sqrtQ,measurementMODEL1.second.leftuv.at(j),config->R_C0toI,config->p_IinC0,config->R_C0toI,config->p_IinC0);
            graphMODEL1->add(factor);
            graphMODEL2->add(factor);
            graphFORSTER->add(factor);
            graph_newMODEL1->add(factor);
            graph_newMODEL2->add(factor);
            graph_newFORSTER->add(factor);
        }

        // Next lets add all RIGHT factors (all graphs have the same measurements!!!)
        for(size_t j=1; j<measurementMODEL1.second.rightuv.size(); j++) {
            Eigen::Matrix<double,2,2> sqrtQ = config->sigma_camera_sq*Eigen::Matrix<double,2,2>::Identity();
            InvUVFactor factor(X(measurementMODEL1.second.rightstateids.at(j)),X(measurement_anchor_lookup[ct_features]),F(ct_features),sqrtQ,measurementMODEL1.second.rightuv.at(j),config->R_C1toI,config->p_IinC1,config->R_C0toI,config->p_IinC0);
            graphMODEL1->add(factor);
            graphMODEL2->add(factor);
            graphFORSTER->add(factor);
            graph_newMODEL1->add(factor);
            graph_newMODEL2->add(factor);
            graph_newFORSTER->add(factor);
        }

        // Record our success
        ct_successes++;

    }

    // Debug info
    if(ct_successes+ct_failures > 0) {
        ROS_WARN("[FEAT]: %d of %d features successfully initialized (inverse depth)", ct_successes, ct_successes + ct_failures);
    }


    //==============================================================================
    // Remove the "used up" features from the queue
    for (auto&& key : toremove) {
        measurement_queue.erase(key);
    }


}
