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
	num_particles = 50;

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
		/* Check Initialization */
		// cout << "particles  " << particles[i].id << " x = " << particles[i].x << endl;
		// cout << "particles  " << particles[i].id << " x = " << particles[i].weight << endl;		
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
			particles[i].x += velocity*(cos(particles[i].theta)) + dist_x(generator);
			particles[i].y += velocity*(sin(particles[i].theta)) + dist_y(generator);
			particles[i].theta += dist_theta(generator);
			// cout << "particles  " << particles[i].id << " x = " << particles[i].x << endl;
		}
		else
		{
			particles[i].x += velocity/yaw_rate*(sin(particles[i].theta+yaw_rate)-sin(particles[i].theta)) + dist_x(generator);
			particles[i].y += velocity/yaw_rate*(-cos(particles[i].theta+yaw_rate)+cos(particles[i].theta)) + dist_y(generator);
			particles[i].theta += yaw_rate + dist_theta(generator);		
		}
	}
}


void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations, std::vector<LandmarkObs>& associations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	// Assign predicted landmark to observed land mark
	// loop over each predicted landmark 
	std::vector<LandmarkObs>copyObserv = observations;

	double minError;	
	if (predicted.size() > copyObserv.size())
	{
		cout << "Prediction size is greater than observation size" << endl;
	}
    for (uint32_t i = 0; i < predicted.size(); i++) {
    	 //    for (uint32_t j = 0; j < copyObserv.size(); j++) {
    	 //    cout << "copyObserv " << j << " = "<< copyObserv[j].id << endl;
    	 //    cout << "copyObserv " << j << " = "<< copyObserv[j].x << endl;    	    
    		// }   
    	    // cout << "copyObserv size = "<< copyObserv.size() << endl;

    		// setup variable observation size
    	   	uint32_t prevPos = 0;
    	   	minError = 0.0;
    	    for (uint32_t j = 0; j < copyObserv.size(); j++) {
    	    	  double eachError = dist(predicted[i].x,predicted[i].y,copyObserv[j].x,copyObserv[j].y);
    	    	  if (eachError < minError || j == 0) // initialize or find min Error 
    	    	  {
    	    	  	minError = eachError;  // set to min error
    	    	  	LandmarkObs sAssoc;
    	    	  	sAssoc.id = copyObserv[j].id; // associate predicted landmark with observed landmark
    	    	  	sAssoc.x = copyObserv[j].x;
    	    	  	sAssoc.y = copyObserv[j].y; 
    	    	  	if (j == 0)
    	    	  	{
    	    	  		associations.push_back(sAssoc); 
    	    	  	}
    	    	  	else
    	    	  	{
    	    	  		associations.at(i) = sAssoc; // replace updated association at association position 	    	  	    	    	  	

    	    	  	} 
    	    	  	prevPos = j;
    	    	  }
    	    	  if (j == copyObserv.size()-1) // reached the end of the availabile observations
    	    	  {
    	    	  	// cout << "Erased copyObserv id = "<< copyObserv[prevPos].id << endl;

    	    	  	copyObserv.erase(copyObserv.begin()+prevPos); // remove observation from vector for future data associations
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
	std::vector<double> out = {0.0, 0.0, 0.0};
	double observ[3]; // single set of observations 
	double sum_w = 0.0; // sum of weights
	double heading_part, x_part, y_part;
	
	std::vector<LandmarkObs> mapObs;
    // cout << "test 2" << endl;

    for(int landmarkInd = 0; landmarkInd < map_landmarks.landmark_list.size(); landmarkInd++) // set map_landmarks to map Obs vector
    {
    	LandmarkObs sLandMark;
		sLandMark.id = map_landmarks.landmark_list[landmarkInd].id_i;
		sLandMark.x  = map_landmarks.landmark_list[landmarkInd].x_f;
		sLandMark.y  = map_landmarks.landmark_list[landmarkInd].y_f;	
		mapObs.push_back(sLandMark);
    }
    // cout << "test 3" << endl;

	//for (int i = 0; i < num_particles; i++) {
	for (int partInd = 0; partInd < num_particles; partInd++) 
	{    
		heading_part = particles[partInd].theta; // particle heading
		x_part 		 = particles[partInd].x;
		y_part 		 = particles[partInd].y;
		// transform each land mark observation into the map frame w.r.t to the particle
		transform_row[0] = {cos(heading_part), -sin(heading_part), x_part}; 
		transform_row[1] = {sin(heading_part), cos(heading_part), y_part},
		transform_row[2] = {0.0, 0.0, 1.0};
		std::vector<LandmarkObs> particleObs;
		std::vector<LandmarkObs> associations;
    	// cout << "test 4" << endl;

        for(int observInd = 0; observInd < observations.size(); observInd++)
        {
 			observ[0] = observations[observInd].x;
			observ[1] = observations[observInd].y;
			observ[2] = 1.0;     
    		// cout << "test 5" << endl;

			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					// cout << "transform_row[" << j << "]"<< "[" << k << "] =" << transform_row[j][k] << endl;
					// cout << "observ[" << k << "] = " << observ[k] << endl;
					out[j] += transform_row[j][k]*observ[k]; // output land mark observations in the map frame w.r.t particle 
					// cout << "out[" << j << "] = " << out[j] << endl;					
				}
			}
			LandmarkObs sParticleObs;
			sParticleObs.id = observations[observInd].id; // set land mark observations for particle
			sParticleObs.x = out[0]; // transformed x
			sParticleObs.y = out[1]; // transformed y
			particleObs.push_back(sParticleObs);			
		}
    	// cout << "test 5" << endl;		
		dataAssociation(particleObs, mapObs, associations); // associate each observed land mark for the particle to the map TODO check assocations make sense
		std::vector<int> association_ids;
		std::vector<double> sense_x;
		std::vector<double> sense_y;
    	// cout << "test 6" << endl;		
		// for(int observInd = 0; observInd < particleObs.size(); observInd++)
  //       {
  //   		cout << "particleObs id = " << particleObs[observInd].id << endl;
		// 	cout << "particleObs x" << particleObs[observInd].x << endl;
		// 	cout << "particleObs y" << particleObs[observInd].y << endl;
  //       }
		// cout << "particleObs size = " <<particleObs.size() << endl;    	
		// cout << "associations size = " <<associations.size() << endl;
		// extract struct into individiual vectors 
		for(int assocInd = 0; assocInd < associations.size(); assocInd++)
        {
        	association_ids.push_back(associations[assocInd].id);
   //  		cout << "associations id = " << associations[assocInd].id << endl;
			// cout << "associations x" << associations[assocInd].x << endl;
			// cout << "associations y" << associations[assocInd].y << endl;			
        	sense_x.push_back(associations[assocInd].x);
        	sense_y.push_back(associations[assocInd].y);
        }
    	// cout << "test 7" << endl;		

		SetAssociations(particles[partInd],association_ids,sense_x,sense_y); // set all vectors of associations of id, position x, y per particle
		// particles[partInd].sense_x = update_particle.sense_x;
		// particles[partInd].sense_y = update_particle.sense_y;
		// particles[partInd].associations = update_particle.association_ids;

		// for(int i = 0; i < associations.size(); i++)
  //       {
		// 	cout << "particleObs x " << i << particleObs[i].x << endl;
		// 	cout << "particleObs y " << i << particleObs[i].y << endl;        	
		// 	cout << "associationsObs x " << i << particles[partInd].sense_x[i] << endl;
		// 	cout << "associationsObs y " << i << particles[partInd].sense_y[i] << endl;
		// 	cout << "particle Obs " << i << " distance err = " << dist(particleObs[i].x,particleObs[i].y,particles[partInd].sense_x[i],particles[partInd].sense_y[i]) << endl;
  //       }
		// cout << "particleObs size = " <<particleObs.size() << endl;    	
		// cout << "weights size = " << weights.size() << endl;    	

		double currentWeight = 0.0;			

		// // compute the multi-variable gaussian distribution by taking the associated map land mark with the observed land mark
		for (int partObsInd = 0; partObsInd < particleObs.size(); partObsInd++)
		{
			// cout << "partObsInd"  << partObsInd << endl;
			// cout << "particleObs x" << particleObs[partInd].x << endl;
			// cout << "particleObs y" << particleObs[partInd].y << endl;
			// cout << "associationsObs x" << particles[partInd].sense_x[partObsInd] << endl;
			// cout << "associationsObs y" << particles[partInd].sense_y[partObsInd] << endl;	
			// check sensor range 
			double errorDist = dist(particleObs[partObsInd].x,particleObs[partObsInd].y, 
									particles[partInd].sense_x[partObsInd], particles[partInd].sense_y[partObsInd]);
			double xerr = particleObs[partObsInd].x-particles[partInd].sense_x[partObsInd];
			double yerr = particleObs[partObsInd].y-particles[partInd].sense_y[partObsInd];
			double exponent = xerr*xerr/(2.0*std_landmark[0])+yerr*yerr/(2.0*std_landmark[1]);
			double denom = M_PI*2.0*std_landmark[0]*std_landmark[1];
			if (errorDist < sensor_range) // check that land mark is within sensor range else add no weight
			{
				// cout << "test 8a" << endl;		
				currentWeight = denom*exp(-exponent);
			}
			else
			{
				currentWeight = 0.0;
			}	
			particles[partInd].weight += currentWeight; // update weight by multiplying all pdfs together
			// cout << "denom" << denom << endl;
			// cout << "exp(-exponent)" << exp(-exponent) << endl;
			// cout << "particle Obs " << partObsInd << " currentWeight = " << currentWeight << endl;

			// cout << "particle Obs " << partObsInd << " errorDist = " << errorDist << endl;																									
			// cout << "particle " << partInd << " weight = " << particles[partInd].weight << endl;
			// cout << "particle " << partInd << " weight = " << particles[partInd].weight << endl;									
		}
		// cout << "test 8" << endl;
		// cout << "weights" << partInd << " = " << weights[partInd] << endl;
		sum_w += particles[partInd].weight;
	}
	// Normalize weights of each particle
	for (int partInd = 0; partInd < num_particles; partInd++)
	{
		particles[partInd].weight = particles[partInd].weight/sum_w;
		weights[partInd] = particles[partInd].weight; 
		cout << "Normalize particle " << partInd << " weight = " << particles[partInd].weight << endl;
	}		

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	// particles[i].weight

  	std::default_random_engine gen;
	std::discrete_distribution<double> distro(weights.begin(),weights.end());

    std::map<int, int> m; // sample n particles using initial weights 
    
    for(int n=0; n<num_particles; n++) {
    	particles[n].weight = 0;
        ++m[distro(gen)];
    }

    int particle_index;
    for (auto p : m)
    {
    	particle_index = p.first;
    	particles[particle_index].weight = (double)p.second/num_particles;
    	cout << " particle " << p.first << " samples = " << p.second << "times" << endl;
    }
	for (int partInd = 0; partInd < num_particles; partInd++)
	{    
		weights[partInd] = particles[particle_index].weight;
    	cout << "Sample particle " << partInd << " weight = " << particles[partInd].weight << endl;
    }
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
