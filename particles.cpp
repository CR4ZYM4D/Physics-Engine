#include "Included files/particles.h"
#include <assert.h>
void Particle :: integrate(real duration){

	if (inverseMass <= 0.0f)
		return; // leave if particle is mass less

	assert(duration > 0.0); // assert function takes a condition parameter and executes only when given assertion is true else breaks

	//updating position , velocity & acceleration
	position.addScaledVector(velocity, duration); 

	Vector3D newAcceleration = acceleration;
	newAcceleration.addScaledVector(netForce, inverseMass);
	
	velocity.addScaledVector(newAcceleration, duration); 
	
	velocity *= drag; 

	clearAccumulator();
}

// defining all particles.h functions

void Particle::setMass(const real mass) {
	assert(mass != 0);
		Particle :: inverseMass = ((real)1.0) / mass; 
}

real Particle::getMass() const{
	if (Particle::inverseMass == 0.0f)
		return REAL_MAX;

		return ((real)1.0) / inverseMass;
	
}

void Particle :: setInverseMass(const real inverseMass){
	Particle :: inverseMass = inverseMass;
}


real Particle::getInverseMass() const {
	return inverseMass;
}

bool Particle::hasFiniteMass() const{
	return Particle::inverseMass > 0.0f;
}

void Particle::setDrag(const real Drag) {
	Particle::drag = drag;
}

real Particle::getDrag()const {
	return Particle::drag;
}

void Particle :: setPosition(const Vector3D &position){
	Particle::position = position;
}

void Particle::setPosition(const real  x, const real  y, const real z) {
	Particle::position.x = x;
	Particle::position.y = y;
	Particle::position.z = z;
}

Vector3D Particle::getPosition()const {
	return Particle::position;
}

void Particle::getPosition(Vector3D* position)const {
	*position = Particle::position;
}

void Particle::setVelocity(Vector3D& velocity) {
	Particle::velocity = velocity;
}

void Particle::setVelocity(const real  x, const real  y, const real z) {
	Particle::velocity.x = x;
	Particle::velocity.y = y;
	Particle::velocity.z = z;
}

Vector3D Particle::getVelocity ()const {
	return Particle::velocity;
}

void Particle::getVelocity(Vector3D *velocity)const{
	*velocity = Particle::velocity;
}

void Particle::setAcceleration(Vector3D& acceleration) {
	Particle::acceleration = acceleration;
}

void Particle::setAcceleration(const real  x, const real  y, const real z) {
	Particle::acceleration.x = x;
	Particle::acceleration.y = y;
	Particle::acceleration.z = z;
}

Vector3D Particle::getAcceleration()const {
	return Particle::acceleration;
}

void Particle::clearAccumulator() {
	netForce.clear();
}

void Particle :: addForce(const Vector3D& force){
	netForce += force;
}