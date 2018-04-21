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
	num_particles = 100;

	default_random_engine generator;
    normal_distribution<double> dist_x(0, std[0]);
    normal_distribution<double> dist_y(0, std[1]);
    normal_distribution<double> dist_theta(0, std[2]);

	for (int i = 0; i<num_particles; i++) {
		Particle oneParticle;
		oneParticle.id     = i;	
		oneParticle.x 	   = x + dist_x(generator);
		oneParticle.y 	   = y + dist_y(generator);
		oneParticle.theta  = theta + dist_theta(generator);
		oneParticle.weight = 1.f;

		particles.push_back(oneParticle);
		/* Check Initialization */
		// cout << "particles  " << particles[i].id << " x = " << particles[i].x << endl;
	}
	is_initialized = true;
}


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine generator;
    normal_distribution<double> dist_x(0, std_pos[0]);
    normal_distribution<double> dist_y(0, std_pos[1]);
    normal_distribution<double> dist_theta(0, std_pos[2]);

	// Update each particle //
	for (int i = 0; i<num_particles; i++) {
		/* Motion Model */
		particles[i].x += velocity/yaw_rate*(sin(particles[i].theta+yaw_rate)-sin(particles[i].theta)) + dist_x(generator);
		particles[i].y += velocity/yaw_rate*(-cos(particles[i].theta+yaw_rate)+cos(particles[i].theta)) + dist_y(generator);
		particles[i].theta += yaw_rate + dist_theta(generator);
		// cout << "particles  " << particles[i].id << " x = " << particles[i].x << endl;
	}

}


// #define MIN(a,b) ({ typeof(a) _a = a; \ // Macro to compute min for any type
//         typeof(b) _b = b; \
//         (((_a)<(_b))? (_a) : (_b)); })

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations, std::vector<LandmarkObs>& associations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	// Assign predicted landmark to observed land mark
	// loop over each predicted landmark 
	std::vector<LandmarkObs>copyObserv = observations;

	double minError;	
    for (uint32_t i = 0; i < predicted.size(); i++) {
    		// setup variable observation size
    	   	uint32_t prevPos = 0;
    	    for (uint32_t j = 0; j < copyObserv.size(); j++) {
    	    	  double eachError = dist(predicted[i].x,predicted[i].y,copyObserv[j].x,copyObserv[j].y);
    	    	  if (eachError < minError || j == 0) // initialize or find min Error 
    	    	  {
    	    	  	minError = eachError;  // set to min error
    	    	  	associations[i].id = copyObserv[j].id; // associate predicted landmark with observed landmark
    	    	  	associations[i].x = copyObserv[j].x;
    	    	  	associations[i].y = copyObserv[j].y;    	    	  	    	    	  	
    	    	  	copyObserv.erase(copyObserv.begin()+j-prevPos); // remove observation from vector for future data associations
    	    	  	prevPos = j+1;
    	    	  }
    	    }
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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
	std::vector<double> transform_row[3];
	std::vector<double> out;
	std::vector<double> observ; // single set of observations 

	double heading_part, x_part, y_part;
	
	std::vector<LandmarkObs> mapObs;

    for(int i = 0; i < map_landmarks.landmark_list.size(); i++) // set map_landmarks to map Obs vector
    {
		mapObs[i].id = map_landmarks.landmark_list[i].id_i;
		mapObs[i].x  = map_landmarks.landmark_list[i].x_f;
		mapObs[i].y  = map_landmarks.landmark_list[i].y_f;	
    }

	for (int i = 0; i < num_particles; i++) {
		heading_part = particles[i].theta; // particle heading
		x_part 		 = particles[i].x;
		y_part 		 = particles[i].y;
		// transform each land mark observation into the map frame w.r.t to the particle
		transform_row[0] = {cos(heading_part), -sin(heading_part), x_part}; 
		transform_row[1] = {sin(heading_part), cos(heading_part), y_part},
		transform_row[2] = {0.0, 0.0, 1.0};
		std::vector<LandmarkObs> particleObs;
		std::vector<LandmarkObs> associations;

        for(int i = 0; i < observations.size(); i++)
        {
 			observ[0] = observations[i].x;
			observ[1] = observations[i].y;
			observ[2] = 1.f;       	
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					out[j] += transform_row[j][k]*observ[k]; // output land mark observations in the map frame w.r.t particle 
				}
			}
			particleObs[i].id = observations[i].id; // set land mark observations for particle
			particleObs[i].x = out[0]; // transformed x
			particleObs[i].y = out[1]; // transformed y			
		}
		dataAssociation(particleObs, mapObs, associations); // associate each observed land mark for the particle to the map
		std::vector<int> association_ids;
		std::vector<double> sense_x;
		std::vector<double> sense_y;

		// extract struct into individiual vectors
		for(int i = 0; i < associations.size(); i++)
        {
        	association_ids[i] = associations[i].id;
        	sense_x[i] = associations[i].x;
        	sense_y[i] = associations[i].y;
        }

		SetAssociations(particles[i],association_ids,sense_x,sense_y); // set all vectors of associations of id, position x, y per particle
	
		// // compute the multi-variable gaussian distribution by taking the associated map land mark with the observed land mark
		// double distro
		// for (int i =; i < particleObs.size(); i++)
		// {
			// normal_distribution<double> dist_theta(error, error);

		// }
	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

    Particle particle_test;

    return particle_test;
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
