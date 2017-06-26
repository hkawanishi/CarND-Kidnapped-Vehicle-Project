/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 50;

	std::default_random_engine gen;

	std::normal_distribution<double> dist_x(x, std[0]);
	std::normal_distribution<double> dist_y(y, std[1]);
	std::normal_distribution<double> dist_theta(theta, std[2]);

	for (int i = 0; i < num_particles; i++)
	{
		Particle particle;
		particle.id = i;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);

		particle.weight = 1;

		particles.push_back(particle);
		weights.push_back(1);

		//std::cout << "i = " << i << "(" << particle.x << "," << particle.y << "," << particle.theta << ")\n";
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;

	for (int i = 0; i < num_particles; i++)
	{
		double new_x;
		double new_y;
		double new_theta;

		if (yaw_rate == 0){
			new_x = particles[i].x + velocity*delta_t*cos(particles[i].theta);
			new_y = particles[i].y + velocity*delta_t*sin(particles[i].theta);
			new_theta = particles[i].theta;
		} else {
			new_x = particles[i].x + velocity/yaw_rate * (sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta));
			new_y = particles[i].y + velocity/yaw_rate * (cos(particles[i].theta)-cos(particles[i].theta+yaw_rate*delta_t));
			new_theta = particles[i].theta + yaw_rate*delta_t;
		}

		normal_distribution<double> dist_x(new_x, std_pos[0]);
	  normal_distribution<double> dist_y(new_y, std_pos[1]);
	  normal_distribution<double> dist_theta(new_theta, std_pos[2]);


		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);

		//std::cout << "i = " << i << "(" << particles[i].x << "," << particles[i].y << "," << particles[i].theta << ")\n";

	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	for (int i = 0; i < num_particles; i++){

		vector<LandmarkObs> trans_observations;
		LandmarkObs obs;
		for (int j = 0; j < observations.size(); j++)
		{
			LandmarkObs trans_obs;
			obs = observations[j];
			//cout << "obs are " << obs.x << ", " << obs.y << "\n";
			//cout << "partcile coord is " << particles[i].x << ", " << particles[i].y << "\n";

			trans_obs.x = particles[i].x+(obs.x*cos(particles[i].theta)-obs.y*sin(particles[i].theta));
			trans_obs.y = particles[i].y+(obs.x*sin(particles[i].theta)+obs.y*cos(particles[i].theta));
			//cout << "trans_obs are " << trans_obs.x << ", " << trans_obs.y << "\n";
			trans_observations.push_back(trans_obs);
		}
		//cout << "total observed: " << trans_observations.size() << "\n";

		vector<LandmarkObs> predicted;
		LandmarkObs landmarks_in_range;

		// go through all landmarks and add only the ones in the sensor_range
    for (int j = 0; j < map_landmarks.landmark_list.size(); j++)
    {
    	float diff_x = fabs(map_landmarks.landmark_list[j].x_f -  particles[i].x);
    	float diff_y = fabs(map_landmarks.landmark_list[j].y_f -  particles[i].y);
    	double calc_dist = sqrt(diff_x*diff_x + diff_y*diff_y);
      //if ((diff_x < sensor_range) && (diff_y < sensor_range))
      if (calc_dist < sensor_range)
      {
      	landmarks_in_range.x = map_landmarks.landmark_list[j].x_f;
      	landmarks_in_range.y = map_landmarks.landmark_list[j].y_f;
      	landmarks_in_range.id = map_landmarks.landmark_list[j].id_i;

      	//cout << map_landmarks.landmark_list[j].x_f << ", " << map_landmarks.landmark_list[j].y_f << ", " << map_landmarks.landmark_list[j].id_i << "\n";
      	predicted.push_back(landmarks_in_range);
      }
    }
    //cout << "total predicted:" << predicted.size() << "\n";

    for (int j = 0; j < trans_observations.size(); j++)
		{
			float minimum_dist = 99999999999; // pick a large number.
			int min_id;
			for (int k = 0; k < predicted.size(); k++){
				float distance = dist(trans_observations[j].x,trans_observations[j].y,predicted[k].x,predicted[k].y);
				if (distance < minimum_dist) {
					minimum_dist = distance;
					min_id = predicted[k].id;
				}
			}
			trans_observations[j].id = min_id;
			//cout << "trans_obs are " << trans_observations[j].x << ", " << trans_observations[j].y << "\n";
			//cout << "j = " << j << ", " << trans_observations[j].id << "\n";
		}

		double sig_x = std_landmark[0];
		double sig_y = std_landmark[1];
		double first_term = 1./(2.*M_PI*sig_x*sig_y);
		particles[i].weight = 1.0;

		for (int j = 0; j < trans_observations.size(); j++)
		{
			double land_id = trans_observations[j].id;
			double predict_x, predict_y;
			for (int k = 0; k < predicted.size(); k++)
			{
				if (predicted[k].id == land_id){
					predict_x = predicted[k].x;
					predict_y = predicted[k].y;
					//cout << "j = " << j << " k = " << k << ", predict_x = " << predicted[k].x << ", predict_y = " << predicted[k].y << ", id = " << predicted[k].id << "\n";
					break;
				}
			}

			double weight;
			double x_term = (trans_observations[j].x-predict_x)*(trans_observations[j].x-predict_x)/(2.*sig_x*sig_x);
			double y_term = (trans_observations[j].y-predict_y)*(trans_observations[j].y-predict_y)/(2.*sig_y*sig_y);
			weight = first_term * exp(-(x_term+y_term));
			//cout << "j = " << j << ", weight = " << weight << "\n";
			if (weight > 0){
				particles[i].weight *= weight;
			}

		}
		weights[i] = particles[i].weight;
		//cout << "i = " << i << " particles weight = " << particles[i].weight << "\n";
  }

}	

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	default_random_engine gen;
	discrete_distribution<int> distribution(weights.begin(), weights.end());

	vector<Particle> resample_particles;


	for (int i = 0; i < num_particles; i++)
	{
		resample_particles.push_back(particles[distribution(gen)]);
	}
	
	particles = resample_particles;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
