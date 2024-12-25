#include "Included files/particlelinking.h"

real ParticleLink::currentLength() const {	
	return (particles[0]->getPosition() - particles[1]->getPosition()).magnitude();
}

unsigned Cable::addContact(ParticleContacts* contact, unsigned int limit)const {

	real length = currentLength();

	if (length < maxLength) //checking if cable is over stretched
		return 0;

	//setting the two particles in contact of cable
	contact->particles[0] = particles[0];
	contact->particles[1] = particles[1];

	//calculating the contact normal
	Vector3D normal = particles[1]->getPosition() - particles[0]->getPosition();
	normal.normalization();
	contact->contactnormal = normal;

	//calculating how much the cable stretched 
	contact->penetration = length - maxLength;
	contact->coeffrestitute = coeffrestitution;
	return 1;

}

unsigned Rod::addContact(ParticleContacts* contact, unsigned int limit)const {

	real currlength = currentLength();

	if (currlength == length)
		return 0;

	//calcultaing contact normal vector
	Vector3D normal = particles[1]->getPosition() - particles[0]->getPosition();
	normal.normalization();

	contact->penetration = currlength - length;
	if (currlength > length)
		contact->contactnormal = normal;
	
	else {//we would need to send particles in opposite direction if the rod is compressed
	
		contact->contactnormal = normal * -1;
		contact->penetration = length - currlength;

	}

	contact->coeffrestitute = 0;
	return 1;

}

real ParticleConstraints::currentLength()const {
	return(particle->getPosition() - anchor).magnitude();
}

unsigned CableConstraint :: addContact(ParticleContacts* contact , unsigned int limit)const{

	real length = currentLength();

	if (length < maxLength)
		return 0;

	Vector3D normal = particle->getPosition() - anchor;
	normal.normalization();
	contact->contactnormal = normal;
	contact->penetration = length - maxLength;
	contact->coeffrestitute = coeffrestitution;
	return 1;

}

unsigned RodConstraints::addContact(ParticleContacts* contact, unsigned int limit)const {

	real currlength = currentLength();

	if (length == currlength)
		return 0;

	Vector3D normal = (particle->getPosition() - anchor)*-1;
	normal.normalization();

	if (currlength > length) {
		contact->penetration = currlength - length;
		contact->contactnormal = normal;
	}
	else {
		contact->contactnormal = normal * -1;
		contact->penetration = length - currlength;
	}
	contact->coeffrestitute = 0;
	return 1;

}