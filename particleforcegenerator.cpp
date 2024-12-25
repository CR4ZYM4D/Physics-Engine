#include "Included files/particleforcegenerator.h"

void ForceRegistry::updateForces( real duration) { // creating the updateForce function of registry class

	Register::iterator i = registrations.begin(); // iterator of type Register to iterate thorugh the vector
	for (; i < registrations.end(); i++) {
		i->fg->updateForce(i->particle, duration); //calling the ForceGenerator object of i to update the force on the particle of i for given duration
	}

}

void ForceRegistry :: add(Particle* particle , ForceGenerator* fg){

	ForceRegistry::ParticleRegistry registration; //creating an object registration to  add to the registrations vector
	registration.particle = particle;
	registration.fg = fg;
	registrations.push_back(registration); // inserting the registration object at the end of the vector

}
ParticleGravity :: ParticleGravity(const Vector3D& gravity):gravity(gravity)
{	//setting value of const gravity to be used here as constructor parameter
}

void ParticleGravity::updateForce(Particle* particle, real duration) {
	
	if (!particle->hasFiniteMass())
		return;

	particle->addForce(gravity * particle ->getMass());
}

ParticleDrag :: ParticleDrag( real dragcoefficient , real dragcoeficientsq ) : dragcoefficient(dragcoefficient) , dragcoefficientsq(dragcoeficientsq){
}

void ParticleDrag :: updateForce(Particle* particle , real duration){

	Vector3D force;
	particle->getVelocity(& force);

	real dragval = force.magnitude()*(dragcoefficient + dragcoefficientsq*dragcoefficientsq);
	force.normalization();
	particle->addForce( force* (-dragval));

}

ParticleSpring::ParticleSpring(Particle* other, real springConst, real lknot):other(other) , springConstant(springConst) , lknot(lknot){
}

void ParticleSpring::updateForce(Particle* particle, real duration) {

	Vector3D force;
	particle->getPosition(&force); // getting position of point at which force is applied i.e. one end of spring 
	force -= other->getPosition(); // subtracting position of other point to get the current extended / compressed length vector of spring

	real magnitude = real_abs(force.magnitude() - lknot); //subtracting current length from rest length to obtain net change in length and then taking its absolute value
	force.normalization(); //normailzing current force to get its unit vector
	magnitude *= springConstant; //calcultaing magnitude of net force applied be spring
	force *= -magnitude; // force applies in opposite direction thus multiplying it's unit vector with -1 times the magnitude of force applied
	particle->addForce(force);

}

AnchoredSpring :: AnchoredSpring(){}

AnchoredSpring :: AnchoredSpring(Vector3D* a , real sc , real rl) : anchor(a) , springConstant(sc) , lknot(rl){}

void AnchoredSpring::init(Vector3D* anchor, real springConst, real restLength) { //initializing the values of new anchored spring object

	AnchoredSpring::anchor = anchor;
	AnchoredSpring::springConstant = springConst;
	AnchoredSpring::lknot = restLength;

}

void AnchoredSpring::updateForce(Particle* particle, real duration) {

	Vector3D force;
	particle->getPosition(&force);
	force -= *anchor;

	real magnitude = force.magnitude();
	magnitude = springConstant * (magnitude - lknot);
	force.normalization();
	force *= -magnitude;
	particle->addForce(force);

}

FakeSpring :: FakeSpring(Vector3D* a , real sc , real rl , real d):anchor(a) , springConstant(sc) , lknot(rl) , damping(d)
{}

void FakeSpring::updateForce(Particle* particle, real duration) {

	if (!particle->hasFiniteMass())
		return;

	Vector3D position;
	particle->getPosition(&position);
	position -= *anchor;

	real gamma = 0.5f * real_sqrt(4*springConstant - damping*damping);
	
	if (gamma == 0.0f)
		return;

	Vector3D c = position * (damping / (2.0f * gamma)) + particle->getVelocity() * (1.0f / gamma);

	Vector3D target = position * real_cos(gamma * duration) + c * real_sin(gamma * duration);
	target *= real_exp(0.5f * duration * damping);
	Vector3D acceleration = (target - position) * ((real)1.0 / duration * duration) - particle->getVelocity()* (real)(1.0 / duration);
	particle->addForce(acceleration * particle->getMass());

}

Bungee::Bungee(Particle* other , real springConst , real restLength):other(other) , springConstant(springConst) , lknot(restLength)
{}

void Bungee::updateForce(Particle* particle, real duration) {

	Vector3D force;
	particle->getPosition(&force); //getting position of one end of bungee
	force -= other->getPosition(); // getting position of other end and evaluating the length vector by vector subtraction

	real magnitude = force.magnitude(); //calculating current length vector magnitude
	
	if (magnitude <= lknot) //checking if bungee is currently compressed or at rest coz if it is, no force would be acting on it
		return;

	//computes the force if bungee is stretched
	magnitude = springConstant * (magnitude - lknot);
	force.normalization();
	force *= -magnitude;
	particle->addForce(force);

}

void AnchoredBungee::updateForce(Particle* particle, real duration) {

	Vector3D force;
	particle->getPosition(&force);
	force -= *anchor;
	real magnitude = force.magnitude();

	if (magnitude <= lknot)
		return;

	force.normalization();
	magnitude = springConstant*(magnitude-lknot);
	force *= -magnitude;
	particle->addForce(force);

}

Buoyancy :: Buoyancy (real v , real h , real ll , real ld): volume(v) , height(h) , liquidLevel(ll) , liquidDensity(ld)
{}

void Buoyancy :: updateForce(Particle* particle , real duration){

	Vector3D force(0.0f ,0.0f , 0.0f);
	real elevation = particle->getPosition().y; //obtaining y coordinate of particle 

	if (elevation >= liquidLevel + height) //checking if particle is completely above the water surface
		return;

	if (elevation <= liquidLevel - height) //checking if particle is completely submerged
		force.y = volume*liquidDensity;

	else {
		force.y = volume*liquidDensity* (real)((elevation - height - liquidLevel)/2*height);
	}

	particle->addForce(force);

}