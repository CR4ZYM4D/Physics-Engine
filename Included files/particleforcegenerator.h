#pragma once
#ifndef PARTICLE_FORCE_GENERATOR
#define PARTICLE_FORCE_GENERATOR
#include "particles.h"
#include <vector> // included vector/list as we will use it to create a register for particles and bodies and add all force acting on it from here

class ForceGenerator { //creating an interface named ForceGenerator which we will use for the genralization of forces 

	public:
		virtual void updateForce(Particle* particle, real duration); //function that will be used to update  the force acting on the particle for given duration

};

class ParticleGravity : public ForceGenerator { //interfaces that inherit ForceGenerator for generalization of gravity & drag respectively

		Vector3D gravity;
	public:
		ParticleGravity(const Vector3D& gravity);
		virtual void updateForce(Particle* particle, real duration);

};

class ParticleDrag : public ForceGenerator {

	real dragcoefficient;
	real dragcoefficientsq; // square of drag coefficient
	public:
		ParticleDrag(real dragcoefficient , real dragcoefficientsq);
		virtual void updateForce(Particle* particle, real duration);

};

class ForceRegistry{ //interface that will be used to implement the force register to compute all force acting on a particle

	private:
		struct ParticleRegistry {

			Particle* particle;
			ForceGenerator* fg;

		};
		typedef std::vector<ParticleRegistry> Register; // defining a data type Register that is a vector that holds objects of type struct ParticleRegistry
		Register registrations; // creating a vector registrations of type Register

	public:
		void add(Particle* particle, ForceGenerator* fg); // will be used to add the given fg obejct to register vector
		void remove(Particle* particle, ForceGenerator* fg); // will be used to remove the fg object from regoster which has completed its duration
		void clear(); // will be used to clear all forces acting on particle from register
		void updateForces( real duration);
};

class ParticleSpring : public ForceGenerator { //creating interface for a spring

	Particle* other;

	real springConstant;

	real lknot;

	public :

		ParticleSpring(Particle* other, real springConstant, real lknot);

		virtual void updateForce(Particle* particle , real duration);

};

class AnchoredSpring : public ForceGenerator { //creating interface for a spring with one end fixed such as one hanging from a ceiling

	protected :
			
		Vector3D* anchor;

		real springConstant;

		real lknot;

	public :

		AnchoredSpring();

		AnchoredSpring(Vector3D* anchor, real springConst , real restLength); //create a new anchored spring

		const Vector3D* getAnchor() const { //retrieve anchor point
			return anchor; 
		}

		void init(Vector3D* anchor , real springConst , real lknot); //set the spring's properties

		virtual void updateForce(Particle* particle, real duration);
};

class FakeSpring : public ForceGenerator { //creating interface to compute force applied by a fully compressed/stiff spring

	Vector3D* anchor;
	real springConstant;
	real lknot;
	real damping;

public:

	FakeSpring(Vector3D* anchor, real SprongConstant, real lknot , real damping);

	virtual void updateForce(Particle* particle, real duration);

};

class Bungee : public ForceGenerator { //creating interface for a bungee/elastic cord

	Particle* other;

	real springConstant;

	real lknot;

	public:

		Bungee(Particle* other, real springCOnstant, real lknot);

		virtual void updateForce(Particle* particle, real duration);

};

class AnchoredBungee : public AnchoredSpring { //interface for bungee cord anchored like an anchored spring

	public:

		virtual void updateForce(Particle* particle , real duration);

};

class Buoyancy : public ForceGenerator { //interface to generate the bouyant force acting on a particle

	real volume; // stores volume of given body/ particle
	real height; // stores elevation of top-most point of particle
	real liquidLevel; // stores height till which liquid/ water is present (water bed is at y=0)
	real liquidDensity; // stores density

	public:

		Buoyancy(real volume, real height, real liquidLevel, real liquidDensity = 1000.0f); //setting default density as that of water

		virtual void updateForce(Particle* particle, real duration);

};
#endif