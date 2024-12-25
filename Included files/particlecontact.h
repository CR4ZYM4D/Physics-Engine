#include "Included Files/particles.h"
#pragma once

#ifndef PARTICLE_CONTACT
#define PARTICLE_CONTACT

class ParticleContactResolver;

class ParticleContacts {

	public:
		Particle* particles[2]; // stores pointers of 2 colliding particles

		real coeffrestitute; // stores coefficient of restitution

		Vector3D contactnormal ; // stores vector of normal of contact point w.r.t. first object 

		real penetration; // will be set to +ve if two free bodies are colliding or bodies connected via rope or rod will be at a distance gretaer than max length apart

		friend class ParticleContactResolver; //friend class can call all variables & methods(public , pvt & protected) of the base class in which it is declared as well

		Vector3D particleMovement[2];

	private:
		void resolveVelocity(real duration); // resolves the velocities of the objects post collision

		void resolvePenetration(real duration); // resolves the overlapping/intersection of two objects

	protected:
		void resolve(real duration);

		real calculateSeperationSpeed() const; // calculates speed with which the particles are moving away from each other before collision 

};

class ParticleContactResolver {

protected:
	unsigned iterations;  // stores the maximum number of iterations to allow   

	unsigned iterationsUsed; // stores the number of itertions done

public:
	ParticleContactResolver(unsigned iterations);

	void setIterations(unsigned iterations);

	void resolveContects(ParticleContacts* contactArray, unsigned numContacts, real duration);

};

class ParticleContactGenerator {

	public:
		virtual unsigned addContact(ParticleContacts* contact , unsigned limit)const = 0;
};

#endif 