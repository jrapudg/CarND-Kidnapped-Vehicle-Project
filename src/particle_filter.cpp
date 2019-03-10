/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using namespace std;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1.
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method
   *   (and others in this file).
  */
  if (is_initialized){
    return;
  }
  std::default_random_engine gen;

  num_particles = 100;  // TODO: Set the number of particles

  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);

  for(unsigned i = 0; i < num_particles; ++i){
    Particle particle;
    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0;

    particles.push_back(particle);
}

  // The filter is now initialized.
  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
   std::default_random_engine gen;

   std::normal_distribution<double> dist_x(0, std_pos[0]);
   std::normal_distribution<double> dist_y(0, std_pos[1]);
   std::normal_distribution<double> dist_theta(0, std_pos[2]);

  // Calculate new state.
  for (unsigned i = 0; i < num_particles; i++) {

  double theta = particles[i].theta;

    if ( fabs(yaw_rate) < 0.0001 ) { // To avoid division by 0
      particles[i].x += velocity * delta_t * cos( theta );
      particles[i].y += velocity * delta_t * sin( theta );

    } else {
      particles[i].x += velocity / yaw_rate * ( sin( theta + yaw_rate * delta_t ) - sin( theta ) );
      particles[i].y += velocity / yaw_rate * ( cos( theta ) - cos( theta + yaw_rate * delta_t ) );
      particles[i].theta += yaw_rate * delta_t;
    }

    // plus noise
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each
   *   observed measurement and assign the observed measurement to this
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will
   *   probably find it useful to implement this method and use it as a helper
   *   during the updateWeights phase.
   */

   unsigned n_obvs = observations.size();
   unsigned n_pred = predicted.size();
   //double sensor_range = 50;

   for (unsigned i = 0; i < n_obvs; ++i) { // For each observation

         // Maximum distance is square root of 2 times the range of sensor
         double min_distance = std::numeric_limits<double>::max();
         //Not possible to access id
         int nearest_landmark_Id = -1;

         for (unsigned j = 0; j < n_pred; ++j ) { // For each predition.

                 double delta_x = observations[i].x - predicted[j].x;
                 double delta_y = observations[i].y - predicted[j].y;

                 double distance = sqrt(delta_x * delta_x + delta_y * delta_y);

                 // If the "distance" is less than min, stored the id and update min.
                 if (distance < min_distance) {
                   min_distance = distance;
                   nearest_landmark_Id = predicted[j].id;
               }
             }
             // Update the observation identifier.
             observations[i].id = nearest_landmark_Id;
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian
   *   distribution. You can read more about this distribution here:
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system.
   *   Your particles are located according to the MAP'S coordinate system.
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
*/
   double weight_normalizer = 0.0;

   for (unsigned i = 0; i < num_particles; ++i){

         double x_part = particles[i].x;
         double y_part = particles[i].y;
         double theta_part = particles[i].theta;

         //(1) Transform observation from car's coordinates to map's coordinates

         vector<LandmarkObs> map_coordinates_observations;

         for(unsigned j = 0; j < observations.size(); ++j){
               double x_obs = observations[j].x;
               double y_obs = observations[j].y;

               LandmarkObs map_coordinates_obs;
               map_coordinates_obs.id = j;
               // transform to map y coordinate
               map_coordinates_obs.x = x_part + (cos(theta_part) * x_obs) - (sin(theta_part) * y_obs);
               map_coordinates_obs.y = y_part + (sin(theta_part) * x_obs) + (cos(theta_part) * y_obs);
               map_coordinates_observations.push_back(map_coordinates_obs);
         }
      // (2) Only keep observations from sensor Range
          vector<LandmarkObs> in_range_landmarks;
          for (unsigned k = 0; k < map_landmarks.landmark_list.size(); k++) {
            Map::single_landmark_s current_landmark = map_landmarks.landmark_list[k];
            if ((fabs((x_part - current_landmark.x_f)) <= sensor_range) && (fabs((y_part - current_landmark.y_f)) <= sensor_range)){
              in_range_landmarks.push_back(LandmarkObs {current_landmark.id_i, current_landmark.x_f, current_landmark.y_f});}
          }

          // (3) Association
          dataAssociation(in_range_landmarks, map_coordinates_observations);

          // (4) Weights calculations

          particles[i].weight = 1.0;

          double sigma_x = std_landmark[0];
          double sigma_y = std_landmark[1];
          double sigma_x_2 = pow(sigma_x, 2);
          double sigma_y_2 = pow(sigma_y, 2);
          double prob_normalizer = (1.0/(2.0 * M_PI * sigma_x * sigma_y));

          //Calculate the weight of particle based on the multivariate Gaussian probability function
          for (unsigned l = 0; l < map_coordinates_observations.size(); ++l) {
            double x_map = map_coordinates_observations[l].x;
            double y_map = map_coordinates_observations[l].y;
            double id_map = map_coordinates_observations[l].id;
            double prob = 1.0;

            for (unsigned n = 0; n < in_range_landmarks.size(); ++n) {
              double in_range_landmark_x = in_range_landmarks[n].x;
              double in_range_landmark_y = in_range_landmarks[n].y;
              double in_range_landmark_id = in_range_landmarks[n].id;

              if (id_map == in_range_landmark_id) {
                prob = prob_normalizer * exp(-1.0 * ((pow((x_map - in_range_landmark_x), 2)/(2.0 * sigma_x_2)) + (pow((y_map - in_range_landmark_y), 2)/(2.0 * sigma_y_2))));
                // To avoid multiplication by zero
                if (prob == 0) {
                  particles[i].weight *= 0.00001;
                } else {
                  particles[i].weight *= prob;
                }

              }
            }
        }

      //weight_normalizer += particles[i].weight;

      }

//(5) Normalize the weights of all particles since resmapling using probabilistic approach
/*for (unsigned i = 0; i < particles.size(); ++i) {
  //particles[i].weight /= weight_normalizer;
  weights[i] = particles[i].weight;
  }
  */
}

void ParticleFilter::resample() {
// TODO: Resample particles with replacement with probability proportional to their weight.
// NOTE: You may find std::discrete_distribution helpful here.
//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
default_random_engine gen;
vector<Particle> resampledParticles;

  // Get weights and max weight.
  vector<double> filtered_weights;
  double mw = numeric_limits<double>::min();
  for(int i = 0; i < num_particles; i++) {
    filtered_weights.push_back(particles[i].weight);
    if ( particles[i].weight > mw ) {
      mw = particles[i].weight;
    }
  }
  // Creating distributions.
  uniform_real_distribution<double> dist_continuous(0.0, mw);
  uniform_int_distribution<int> dist_discrete(0, num_particles - 1);

  // Generating index.
  int index = dist_discrete(gen);

  double beta = 0.0;

  // the wheel
  for(int i = 0; i < num_particles; i++) {
    beta += dist_continuous(gen) * 2.0;
    while(beta > filtered_weights[index]) {
      beta -= filtered_weights[index];
      index = (index + 1) % num_particles;
    }
    resampledParticles.push_back(particles[index]);
  }

  particles = resampledParticles;
}


void ParticleFilter::SetAssociations(Particle& particle,
                                     const vector<int>& associations,
                                     const vector<double>& sense_x,
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association,
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
