/*
 * particle_filter.cpp
 *
 *  completed on: May 23, 2017
 *      Author: David Awad
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"

#define MIN               0.001
#define PI                3.14159
#define DEFAULT_THRESHOLD 40.0
#define NUM_PARTICLES     300

using namespace std;

// shared random engine
default_random_engine gen;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // Set the number of particles. Initialize all particles to first position (based on estimates of
  //   x, y, theta and their uncertainties from GPS) and all weights to 1.
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  // set number of particles some large number for testing
  num_particles = NUM_PARTICLES;

  // resize array of particle weights
  weights.resize(num_particles) ;

  // Standard deviations for x, y, and psi
  double std_x, std_y, std_theta;

  std_x     = std[0]; // meters
  std_y     = std[1]; // meters
  std_theta = std[2]; // radians

  // creates a Gaussian distribution for x, y, and psi.
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_psi(theta, std_theta);

  particles.resize(num_particles) ;

  // adding random gaussian noise to each particle
  for (int i = 0; i < num_particles; i++) {
    // initialize particles to positions
    particles[i].id     = i ;
    particles[i].x      = dist_x(gen) ;
    particles[i].y      = dist_y(gen) ;
    particles[i].theta  = dist_psi(gen) ;
    particles[i].weight = 1.0 ;
  }

  is_initialized = true;
  cout << "Initialization complete" << endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.

  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/

  for (int i = 0; i < num_particles; i++) {

    const double theta = particles[i].theta;

    if ( fabs(yaw_rate) > MIN) {
      const double velocity_yaw_rate = velocity / yaw_rate;
      const double yaw_rate_delta_t  = yaw_rate * delta_t;

      // calculate noise for changing yaw rate
      particles[i].x     += velocity_yaw_rate * (sin(theta + yaw_rate_delta_t) - sin(theta));
      particles[i].y     += velocity_yaw_rate * (cos(theta) - cos(theta + yaw_rate_delta_t));
      particles[i].theta += yaw_rate * delta_t ;

    } else { // the car is going essentially straight, use the alternative formulas
      const double velocity_delta_t = velocity * delta_t ;

      particles[i].x += velocity_delta_t * cos(theta) ;
      particles[i].y += velocity_delta_t * sin(theta) ;
    }

    // add random Gaussian noise to our computed particles
    normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
    normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
    normal_distribution<double> dist_psi(particles[i].theta, std_pos[2]);

    // set values of particles to our new values that include gaussian noise
    particles[i].x     = dist_x(gen) ;
    particles[i].y     = dist_y(gen) ;
    particles[i].theta = dist_psi(gen) ;
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
  //   observed measurement to this particular landmark.

  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
  //   implement this method and use it as a helper during the updateWeights phase.

  // iterate through all observed landmarks and then compare them to our predictions and to then associate them
  // NOTE: this function didn't need to be implemented so it has been left incomplete due to time :)
  //   - sorry, David.
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
    std::vector<LandmarkObs> observations, Map map_landmarks) {
  // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
  //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution

  // process for each particle:
  // - Rotate coordinate system from car's to global map where points are
  // - find the transform from the car to the landmark
  // - and then apply it from where a point is,
  // - then check how close our distance to the point/guess is
  // - relative to our actual measurement for the landmark

  // clear out weights of particles so far
  weights.clear();

  // iterate through all particles
  for (int i = 0; i < num_particles; i++) {

    long double weight_placeholder = 1.0 ;

    // iterate through all observations of that particle
    for (unsigned int j = 0; j < observations.size(); j++) {

      // find the specific angle, the the trig functions for that angle
      const double particle_theta = particles[i].theta ;
      const double cos_theta      = cos(particle_theta) ;
      const double sin_theta      = sin(particle_theta) ;

      // Calculate X from particle's perspective
      // p_x + (obs_x * cos(theta)) - (obs_y * sin(theta))
      double transformed_observation_x =  particles[i].x + (observations[j].x * cos_theta) - (observations[j].y * sin_theta);

      // Calculate Y from particle's perspective
      // p_y + (obs_x * sin(theta)) + (obs_y * cos(theta))
      double transformed_observation_y =  particles[i].y + (observations[j].x * sin_theta) + (observations[j].y * cos_theta);

      // reasonable threshold for an observation to be linked to a landmark
      double current_threshold = DEFAULT_THRESHOLD ;
      // set current id to a particular landmark, start with -1
      long landmark_index = -1;

      // find the closest landmark to this observation
      for (unsigned int k = 0;  k < map_landmarks.landmark_list.size(); k++) {

        // start by taking the distance from this observation to a given landmark
        // save the index of the landmark with the smallest distance
        double difference = dist(transformed_observation_x,
                                 transformed_observation_y,
                                 map_landmarks.landmark_list[k].x_f,
                                 map_landmarks.landmark_list[k].y_f );

        // check if threshold is less than previous
        if (difference < current_threshold) {
          // if it is, then we just found a closer landmark
          landmark_index = map_landmarks.landmark_list[k].id_i -1; // quick -1 for indexing
          // if the new distance is higher it won't update the landmark_index
          current_threshold = difference ;
        }
      }

      // we now have a landmark that is somewhat close to a landmark we've observed,
      // let's calculate our probability distribution using the gaussian formula
      // WARNING: Math incoming
      if (landmark_index >= 0) {

        // x,y values for the points prediction
        const long double x               = transformed_observation_x ;
        const long double y               = transformed_observation_y ;

        // values for the nearest landmark's measured position
        const long double u_x             = map_landmarks.landmark_list[landmark_index].x_f ;
        const long double u_y             = map_landmarks.landmark_list[landmark_index].y_f ;

        // (x - µ_x)^2
        const long double x_ux_squared    = (x - u_x) * (x - u_x);
        // (y - µ_y)^2
        const long double y_uy_squared    = (y - u_y) * (y - u_y);

        // standard error σ_x and σ_y
        const double sigma_x              = std_landmark[0];
        const double sigma_y              = std_landmark[1];

        // cache constant value to multiply with the exponential term
        // 1/ (2 * π * σ_x * σ_y)
        const long double common_denom    = 1 / (2 * PI * sigma_x * sigma_y);

        // cache terms for computing exponential
        // (x - µ_x)^2/(2σ_x^2) + (y - µ_y)^2/(2σy_x^2)
        const long double exp_terms       = (x_ux_squared / (sigma_x * sigma_x)) +
                                            (y_uy_squared / (sigma_y * sigma_y));

        // calculate the Multivariate-Gaussian Probability with this particular observation
        // previous placeholder * cached constants * e^(-1/2 * cached terms)
        const long double current_multivariate = common_denom * exp( (-1/2.) * (exp_terms));

        // update our weight placeholder by multiplying it with the new distribution
        weight_placeholder *= current_multivariate;
      }
    }

    // by taking the product of all of our variate probabilities
    // we then get a good understanding of how accurate this particle is,
    // which is a perfect usecase for it's relative weight
    particles[i].weight = weight_placeholder ;
    weights.push_back(particles[i].weight);
  }
}

void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight.
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  // use either the resampling wheel or some other method to read each variable
  discrete_distribution<> dist(weights.begin(), weights.end());

  // create resampled particles array and store our new particles there
  vector<Particle> resampled_particles ;
  resampled_particles.resize(num_particles) ;

  // use a discrete distribution with the input weights and
  // resample accordingly based on number of occurrences of similar weights
  // this way we will sample similar weights more frequently
  for(int n = 0; n < num_particles; n++) {
    int new_index = dist(gen);
    resampled_particles[n] = particles[new_index];
  }
  particles = resampled_particles ;
}

void ParticleFilter::write(std::string filename) {
  // You don't need to modify this file.
  std::ofstream dataFile;
  dataFile.open(filename, std::ios::app);
  for (int i = 0; i < num_particles; ++i) {
    dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
  }
  dataFile.close();
}
