#include "Included Files/particlecontact.h"

void ParticleContacts::resolve(real duration) { //calling it to resolve the velocty and overlapping

	resolveVelocity(duration);
	resolvePenetration(duration);

}

real ParticleContacts::calculateSeperationSpeed()const {

	Vector3D relativeSpeed = particles[0]->getVelocity();
	if (particles[1]) //checking if the other particle is not null i.e. an anchor
		relativeSpeed -= particles[1]->getVelocity();
	return relativeSpeed*contactnormal;

}

void ParticleContacts::resolveVelocity(real duration) {

	real separatingVelocity = calculateSeperationSpeed();

	if (separatingVelocity > 0) //check if particles are drifting away
		return;

	real newSeparatingSpeed = -separatingVelocity * coeffrestitute; //impulse speed applied by cable


	//calculating accelerations of each particle
	Vector3D netAcceleration = particles[0]->getAcceleration();

	if (particles[1])
		netAcceleration -= particles[1]->getAcceleration();

	real accelerationCausedSpeed = netAcceleration * contactnormal * duration; // calculates separation speed caused solely by the acceleration of the particles 

	if(accelerationCausedSpeed < 0){ //removing this speed if it was bringing the particles closer
	
		newSeparatingSpeed -= accelerationCausedSpeed * coeffrestitute;

		if (newSeparatingSpeed < 0) // making sure we dont accidentally reverse the speed direction
			newSeparatingSpeed = 0;

	}

	real deltaSpeed = newSeparatingSpeed - separatingVelocity;

	//applying velocity change inversely proportional to mass of each object
	
	real totalInverseMass = particles[0]->getInverseMass();

	if (particles[1])
		totalInverseMass += particles[1]->getInverseMass();

	if (totalInverseMass <= 0) //checking if all particles are not of infinite mass
		return;

	Vector3D impulsePerUnitMass = contactnormal * deltaSpeed; //since impulse would be applied in direction of contact normal and 1 unit of mass would have an impulse equal to the velocity

	Vector3D finalVelocity = particles[0]->getVelocity() + impulsePerUnitMass * (totalInverseMass / particles[0]->getInverseMass());

	particles[0]->setVelocity(finalVelocity);

	if (particles[1]) {

		finalVelocity = particles[1]->getVelocity() - impulsePerUnitMass * (totalInverseMass / particles[1]->getInverseMass());// - sign because it will rebound / bounce back in the opposite direction
		particles[1]->setVelocity(finalVelocity);

	}

}

void ParticleContacts :: resolvePenetration(real duration){

	if (penetration <= 0) //checking if particles/bodies are not overlapping at all
		return;

	//the distance that each body bounces off would be proprtional to the ration of their mass to total mass
	 
	real totalInverseMass = particles[0]->getInverseMass();

	if (particles[1])
		totalInverseMass += particles[1]->getInverseMass();

	if (totalInverseMass <= 0)
		return;

	Vector3D reboundPerUnitMass = contactnormal * (penetration * totalInverseMass); //since particles would bounce back in direction of contact normal

	particleMovement[0] = reboundPerUnitMass * (1 / particles[0]->getInverseMass());

	if (particles[1])
		particleMovement[1] = reboundPerUnitMass * (1 / particles[1]->getInverseMass());
	//else
	//	particleMovement[1].clear();

	//setting new positions of particles/bodies

	particles[0]->setPosition(particles[0]->getPosition() + particleMovement[0]);

	if (particles[1])
		particles[1]->setPosition(particles[1]->getPosition() + particleMovement[1]);

}

ParticleContactResolver :: ParticleContactResolver(unsigned iterations) : iterations(iterations){}

void ParticleContactResolver::setIterations(unsigned iterations) {

	ParticleContactResolver::iterations = iterations;

}

void ParticleContactResolver::resolveContects(ParticleContacts* contactArray, unsigned numContacts, real duration) {

	unsigned i;

	iterationsUsed = 0;

	while (iterationsUsed < iterations) 
	{

		real max = REAL_MAX;

		unsigned maxIndex = numContacts;

		//finding the particle/body with that is approaching with highest velocity

		for (i = 0; i < numContacts; i++) {

			real sepVel = contactArray[i].calculateSeperationSpeed();

			if (sepVel < max && (sepVel < 0 || contactArray[i].penetration>0)) {
				
				max = sepVel;
				maxIndex = i;
			
			}

		}

		//checking if any particle needs to be resolved

		if (maxIndex == numContacts)
			break;

		//resolving the contact
		contactArray[maxIndex].resolve(duration);

		//updating penetration for all other particles
		Vector3D* move = contactArray[maxIndex].particleMovement;

		for(i=0 ; i<numContacts ; i++)
		{
		
			if (contactArray[i].particles[0] == contactArray[maxIndex].particles[0]) {

				contactArray[i].penetration -= move[0] * contactArray[i].contactnormal;

			}

			else if (contactArray[i].particles[0] == contactArray[maxIndex].particles[1])
			{
				contactArray[i].penetration -= move[1] * contactArray[i].contactnormal;
			}
			
			if (contactArray[i].particles[1])
			{
				if (contactArray[i].particles[1] == contactArray[maxIndex].particles[0])
				{
					contactArray[i].penetration += move[0] * contactArray[i].contactnormal;
				}
				else if (contactArray[i].particles[1] == contactArray[maxIndex].particles[1])
				{
					contactArray[i].penetration += move[1] * contactArray[i].contactnormal;
				}
			}
		}
		iterationsUsed++;
	}
}