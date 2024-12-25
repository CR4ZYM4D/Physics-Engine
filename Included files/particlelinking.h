#pragma once
#ifndef PARTICLE_LINKING
#define PARTICLE_LINKING

#include "particlecontact.h"
#include "particles.h"

class ParticleLink :public ParticleContactGenerator {

	public:
		Particle* particles[2]; // holds the pair of particles that are connected together

		virtual unsigned addContact(ParticleContacts* contact , unsigned limit)const = 0;
	protected:
		real currentLength()const;
};

class Cable : public ParticleLink{

	public : 
		real maxLength;

		real coeffrestitution;

		virtual unsigned addContact(ParticleContacts* contact, unsigned limit)const;
};

class Rod : public ParticleLink {

	public:
		real length;

		virtual unsigned addContact(ParticleContacts* contact, unsigned limit)const;

};

class ParticleConstraints : public ParticleContactGenerator { //will be used to connect particle to fixed anchor points

	public: 
		Particle* particle;

		Vector3D anchor;

		virtual unsigned addContact(ParticleContacts* contact, unsigned limit)const = 0;

	protected:
		real currentLength()const; // returns current distance b/w particle and anchor

};

class CableConstraint : public ParticleConstraints{ // particle is connected via cable

	public:
		real maxLength;
		
		real coeffrestitution;

		virtual unsigned addContact(ParticleContacts* contact, unsigned limit)const;
};

class RodConstraints : public ParticleConstraints {

	public:
		real length;
		
		virtual unsigned addContact(ParticleContacts* contact, unsigned limit)const;

};
#endif