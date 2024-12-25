#pragma once
#include "particlelinking.h"
#include "particleforcegenerator.h"

#ifndef PARTICLE_WORLD_H
#define PARTICLE_WORLD_H

class ParticleWorld {

	public :
		typedef std::vector<Particle*> ParticleList;
		typedef std::vector<ParticleContactGenerator*> ContactGeneratorList;

	protected :
		ParticleList plist; //list of particles

		bool calculateIterations;

		ForceRegistry registry;

		ParticleContactResolver resolver; //resollver for each contact

		ContactGeneratorList cgenlist; // list of contact generators 

		ParticleContacts* contacts; // list of contacts

		unsigned maxContacts; // stores the max no of contacts allowed 

	public:
		ParticleWorld(unsigned maxContacts, unsigned iterations = 0);

		~ParticleWorld();

		unsigned generateContacts();

		void integrate(real duration);

		void runPhysics(real duration);

		void startFrame(); //resets the force accumulators to run simulation for 

		ParticleList& getParticles(); // return particle list

		ContactGeneratorList& getContactGenerators(); // returns contact generators
		ForceRegistry& getForceRegistries(); //return force registry
};


class GroundContacts : ParticleContactGenerator { //special class for particles in conntact with the ground

	ParticleWorld::ParticleList* plist;

	public:
		void init(ParticleWorld::ParticleList* plist);

		virtual unsigned addContact(ParticleContacts* contact, unsigned limit)const;
};
#endif