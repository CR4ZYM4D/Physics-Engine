#include "Included Files/body.h"
#include <assert.h>
#include <memory.h>

// helper function to find the inverseInertiaTensorGlobal (iitG) matrix of a body using 
// the inverseInertiaTensorLocal (iitL) matrix,  a 3x4 rotation matrix rM and a orientation quaternion 
static inline void transformInertiaTensor(const Matrix3by3& iitL ,const Quaternion& orientation,const Matrix3by4& rM , Matrix3by3 iitG){
	
	//It is an auto generated code by code generator
	//I would be lying if I said I understood the maths applied here


}

// helper function to find the 3x4 transformation matrix of a body using its position vector and rotation quaternion 
static inline void calculateTransformMatrix(Matrix3by4& transformMatrix , const Vector3D& position , const Quaternion& rotation){
	transformMatrix.setOrientationAndPosition(rotation, position);
}

void RigidBody::calculateDerivedData() {
	
	//normalising the orientation quaternion
	orientation.normalise();

	//calculating the transform matrix
	calculateTransformMatrix(transformMatrix , position , orientation);

	//calculating the inverseInertiaTensorGlobal
	transformInertiaTensor(inverseInertiaTensorLocal , orientation , transformMatrix , inverseInertiaTensorGlobal);

}

void RigidBody::integrate(real duration) {

	// adding the linear acceleration caused by new forces
	previousFrameAcceleration = acceleration;
	previousFrameAcceleration.addScaledVector(forceAccumulated, inverseMass);

	//adding the angular acceleration caused by torques
	Vector3D angularAcceleration = inverseInertiaTensorGlobal.transformByVector(torqueAccumulated);

	//updating linear and angular velocities
	linearVelocity.addScaledVector(previousFrameAcceleration, duration);
	angularVelocity.addScaledVector(angularAcceleration, duration);

	//applying drag to linear and angular velocities
	linearVelocity *= real_pow(linearDamping, duration);
	angularVelocity *= real_pow(angularDamping, duration);

	//finding updated position and rotation/angle
	position.addScaledVector(linearVelocity, duration);
	orientation.addScaledVector(angularVelocity, duration);

	//recalculating the transform and inertia matrices
	calculateDerivedData();

	//clearing all force generators
	clearAccumulators();

}

void RigidBody::setMass(real mass) {

	if (mass != 0)
		inverseMass = (real)1.0f / mass;

}

real RigidBody :: getMass() const{

	return (real)1.0f / inverseMass;

}

void setInverseMass(real inverseMass) {

	RigidBody : inverseMass = inverseMass;

}

real RigidBody :: getInverseMass() const{
	
	return inverseMass;

}

bool RigidBody::hasFiniteMass() const{

	return inverseMass > 0.0f;

}

void RigidBody::setInertiaTensor(const Matrix3by3& m) {

	inverseInertiaTensorLocal.setInverse(m);

}

void RigidBody::getInertiaTensor(Matrix3by3* inertiaTensorL) const{

	(*inertiaTensorL).setInverse(inverseInertiaTensorLocal);

}

Matrix3by3 RigidBody::getInertiaTensor() const {

	Matrix3by3 itl;
	getInertiaTensor(&itl);
	return itl;

}

void RigidBody :: getInertiaTensorGlobal(Matrix3by3* itg)const{

	(*itg).setInverse(inverseInertiaTensorGlobal);

}

Matrix3by3 RigidBody :: getInertiaTensorGlobal()const{

	Matrix3by3 itg;
	getInertiaTensorGlobal(&itg);
	return itg;

}

void RigidBody :: setInverseInertiaTensorLocal(const Matrix3by3& m){

	RigidBody: inverseInertiaTensorLocal = m;

}

void RigidBody :: getInverseInertiaTensorLocal (Matrix3by3* iitl)const{

	*iitl = inverseInertiaTensorLocal;

}

Matrix3by3 RigidBody::getInverseInertiaTensorLocal()const {

	return inverseInertiaTensorLocal;

}

