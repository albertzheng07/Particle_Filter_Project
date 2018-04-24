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
#include <map>

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
		oneParticle.weight = 1.0/num_particles;

		particles.push_back(oneParticle);
		weights.push_back(oneParticle.weight);	
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
		if (fabs(yaw_rate) < 0.0001)
		{	
			/* Motion Model */
			particles[i].x += velocity*delta_t*(cos(particles[i].theta)); // v*dt*cos(theta)
			particles[i].y += velocity*delta_t*(sin(particles[i].theta)); // v*dt*sin(theta)
		}
		else
		{
			particles[i].x += velocity/yaw_rate*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta)); // v/w*dt*(-sin(theta)+sin(theta+w*dt))
			particles[i].y += velocity/yaw_rate*(-cos(particles[i].theta+yaw_rate*delta_t)+cos(particles[i].theta)); // v/w*dt*(cos(theta)-cos(theta+w*dt))
			particles[i].theta += yaw_rate*delta_t ; // theta = w*dt		
		}
		// Add noise generators
		particles[i].x +=  dist_x(generator);
		particles[i].y +=  dist_y(generator);
		particles[i].theta += dist_theta(generator);
	}
}


LandmarkObs ParticleFilter::dataAssociation(LandmarkObs predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
   	LandmarkObs neighbor;

   	double minError = 0.0;
    for (uint32_t j = 0; j < observations.size(); j++) {
    	  double eachError = dist(predicted.x,predicted.y,observations[j].x,observations[j].y);
    	  if (eachError < minError || j == 0) // initialize or find min Error 
    	  {
    	  	minError = eachError;  // set to min error
    	  	neighbor.id = observations[j].id; // associate predicted landmark with observed landmark
    	  	neighbor.x = observations[j].x;
    	  	neighbor.y = observations[j].y; 
    	  }
    }

    return neighbor;
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
	std::vector<double> transform_row[2];
	double observ[3]; // single set of observations 
	double sum_w = 0.0; // sum of weights
	
	std::vector<LandmarkObs> mapObs;

    for(int landmarkInd = 0; landmarkInd < map_landmarks.landmark_list.size(); landmarkInd++) // set map_landmarks to map Obs vector
    {
    	LandmarkObs sLandMark;
		sLandMark.id = map_landmarks.landmark_list[landmarkInd].id_i;
		sLandMark.x  = map_landmarks.landmark_list[landmarkInd].x_f;
		sLandMark.y  = map_landmarks.landmark_list[landmarkInd].y_f;	
		mapObs.push_back(sLandMark);
    }

  const double denom = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);  
  const double x_denom = 2 * std_landmark[0] * std_landmark[0];
  const double y_denom = 2 * std_landmark[1] * std_landmark[1];

	//for (int i = 0; i < num_particles; i++) {
	for (int partInd = 0; partInd < num_particles; partInd++) 
	{    
		double heading_part  = particles[partInd].theta; // particle heading
		double x_part 		 = particles[partInd].x;
		double y_part 		 = particles[partInd].y;
		// transform each land mark observation into the map frame w.r.t to the particle
		transform_row[0] = {cos(heading_part), -sin(heading_part), x_part}; 
		transform_row[1] = {sin(heading_part), cos(heading_part), y_part};

		double currentWeight = 1.0; // initialize current weight for each particle

        for(int observInd = 0; observInd < observations.size(); observInd++)
        {
			std::vector<double> transform_out = {0.0, 0.0};
 			observ[0] = observations[observInd].x;
			observ[1] = observations[observInd].y;
			observ[2] = 1.0;
			// // Transformation from observation->map frame
			for (int j = 0; j < 2; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					transform_out[j] += transform_row[j][k]*observ[k]; // output land mark observations in the map frame w.r.t particle 
				}
			}
			LandmarkObs sParticleObs;
			sParticleObs.id = observations[observInd].id; // set land mark observations for particle
			sParticleObs.x = transform_out[0]; // transformed x
			sParticleObs.y = transform_out[1]; // transformed y

			LandmarkObs nearestNeighbor = dataAssociation(sParticleObs, mapObs); // find nearest neighbor map land mark to particle observation 

			// get multi-variate gaussian distribution between observed landmark by particle and nearest neighbor landmark from the map
			double xerr = transform_out[0]-nearestNeighbor.x;
			double yerr = transform_out[1]-nearestNeighbor.y;
    		double exponent = ((xerr * xerr) / x_denom) + ((yerr * yerr) / y_denom);
      		currentWeight *= denom * exp(-exponent); // multiply all observation pdfs		
		}
		particles[partInd].weight = currentWeight; 
		weights[partInd] = particles[partInd].weight; 	
		sum_w += weights[partInd];	
	}
	// Normalize weights of each particle
	for (int partInd = 0; partInd < num_particles; partInd++)
	{
		particles[partInd].weight = particles[partInd].weight/sum_w;
		weights[partInd] = particles[partInd].weight; 
	}		

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  	std::default_random_engine gen;
	std::discrete_distribution<double> distro(weights.begin(),weights.end());

    std::vector<Particle> new_particles(num_particles); // generate new particle vector for replacement
   
    for(int n=0; n<num_particles; n++) {
        new_particles[n] = particles[distro(gen)]; // randomly sample distribution and get index of particle for new particle list
    }

    particles = new_particles; // replace old particle list with new particles
}

void ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
