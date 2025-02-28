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

void RigidBody :: setInverseMass(const real inverseMass) {

	RigidBody :: inverseMass = inverseMass;

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

	RigidBody:: inverseInertiaTensorLocal = m;

}

void RigidBody :: getInverseInertiaTensorLocal (Matrix3by3* iitl)const{

	*iitl = inverseInertiaTensorLocal;

}

Matrix3by3 RigidBody::getInverseInertiaTensorLocal()const {

	return inverseInertiaTensorLocal;

}

void RigidBody::setInverseInertiaTensorGlobal(const Matrix3by3& iitG){

	RigidBody :: inverseInertiaTensorGlobal = iitG;

}

void RigidBody::getInverseInertiaTensorGlobal(Matrix3by3* inverseInertiaTensorG) const
{

	*inverseInertiaTensorG = inverseInertiaTensorGlobal;

}

Matrix3by3 RigidBody::getInverseInertiaTensorGlobal() const {

	return inverseInertiaTensorGlobal;

}

void RigidBody::setDamping(const real linearDamping, const real angularDamping) {

	RigidBody::linearDamping = linearDamping;
	RigidBody::angularDamping = angularDamping;

}

void RigidBody::setLinearDamping(const real linearDamping) {

	RigidBody::linearDamping = linearDamping;

}

real RigidBody::getLinearDamping()const {

	return linearDamping;

}

void RigidBody::setAngularDamping(const real angularDamping) {

	RigidBody::angularDamping = angularDamping;

}

real RigidBody::getAngularDamping()const {

	return angularDamping;

}

void RigidBody::setPosition(const Vector3D& position) {

	RigidBody::position = position;

}

void RigidBody::setPosition(const real x, const real y, const real z) {

	RigidBody::position = Vector3D(x, y, z);

}

Vector3D RigidBody::getPosition()const {

	return position; 

}

void RigidBody::getPosition(Vector3D* p)const {

	*p = position;

}

void RigidBody::setOrientation(const Quaternion& orientation) {

	RigidBody::orientation = orientation;
	RigidBody :: orientation.normalise();

}

void RigidBody::setOrientation(const real r , const real i, const real j, const real k) {

	RigidBody::orientation = Quaternion(r,i,j,k);
	orientation.normalise();

}

Quaternion RigidBody::getOrientation()const {

	return orientation; 

}

void RigidBody::getOrientation(Quaternion* q)const {

	*q = orientation;

}

void RigidBody :: getOrientation(Matrix3by3* om)const{

	om->setOrientation(orientation);

}

void RigidBody :: getOrientation (real matrix[9])const{

	matrix[0] = 1 - 2 * (orientation.j * orientation.j + orientation.k * orientation.k);
	matrix[1] = 2 * (orientation.i * orientation.j + orientation.k * orientation.r);
	matrix[2] = 2 * (orientation.i * orientation.k - orientation.j * orientation.r);
	matrix[3] = 2 * (orientation.i * orientation.j - orientation.k * orientation.r);
	matrix[4] = 1 - 2 * (orientation.i * orientation.i + orientation.k * orientation.k);
	matrix[5] = 2 * (orientation.j * orientation.k + orientation.i * orientation.r);
	matrix[6] = 2 * (orientation.i * orientation.k + orientation.j * orientation.r);
	matrix[7] = 2 * (orientation.j * orientation.k - orientation.i * orientation.r);
	matrix[8] = 1 - 2 * (orientation.j * orientation.j + orientation.i * orientation.i);

}

void RigidBody::getTransform(Matrix3by4* transform)const {

	*transform = transformMatrix;

}

void RigidBody::getTransform(real matrix[16])const {

	memcpy(matrix , transformMatrix.data , sizeof(real)*12);
	matrix[13] = matrix[14] = matrix[12] = 0;
	matrix[15] = 1;

}

Matrix3by4 RigidBody :: getTransform()const{

	return transformMatrix;

}

void RigidBody :: getGLTransform(float matrix[16])const{
 
	 transformMatrix.fillGLArray(matrix);
 
 }

void RigidBody::setLinearVelocity(const Vector3D& velocity) {

	linearVelocity = velocity;

}

void RigidBody::setLinearVelocity(const real x  , const real y , const real z) {

	linearVelocity = Vector3D(x,y,z);

}

void RigidBody::getLinearVelocity(Vector3D* velocity) const{

	*velocity = linearVelocity;

}

Vector3D RigidBody::getLinearVelocity()const {

	return linearVelocity;

}

void RigidBody::addLinearVelocity(const Vector3D& velocity) {

	linearVelocity += velocity;

}

void RigidBody::setAngularVelocity(const Vector3D& velocity) {

	angularVelocity = velocity;

}

void RigidBody::setAngularVelocity(const real x , const real y , const real z) {

	angularVelocity = Vector3D(x,y,z);

}

void RigidBody::getAngularVelocity(Vector3D* velocity)const {

	*velocity = angularVelocity;

}

Vector3D RigidBody::getAngularVelocity()const {

	return angularVelocity;

}

void RigidBody::addAngularVelocity(const Vector3D& velocity) {

	angularVelocity += velocity;

}

void RigidBody::setAcceleration(const Vector3D& a) {

	acceleration = a;

}

void RigidBody::setAcceleration(const real x  ,const real y , const real z ) {

	acceleration = Vector3D(x,y,z);

}

void RigidBody::getAcceleration(Vector3D* a)const {

	*a = acceleration;

}

Vector3D RigidBody::getAcceleration()const {

	return acceleration;

}

void RigidBody :: getPreviousFrameAcceleration(Vector3D* a)const{

	*a = previousFrameAcceleration;

}

Vector3D RigidBody :: getPreviousFrameAcceleration() const {

	return previousFrameAcceleration;

}

void RigidBody::clearAccumulators() {

	forceAccumulated.clear();
	torqueAccumulated.clear();

}

void RigidBody :: addForce(const Vector3D& force){

	forceAccumulated += force;

}

void RigidBody::addForceAtGlobalPoint(const Vector3D& force, const Vector3D& p) {

	forceAccumulated += force;
	Vector3D deltatorque = p;
	deltatorque.addScaledVector(position, -1);
	torqueAccumulated += deltatorque.X(force);

}

void RigidBody::addForceAtLocalPoint(const Vector3D& force, const Vector3D& p) {

	forceAccumulated += force;
	Vector3D deltaTorque = p;
	torqueAccumulated += deltaTorque.X(force);

}

void RigidBody::addTorque(const Vector3D& t) {

	torqueAccumulated += t;

}

Vector3D RigidBody::getPointInGlobalSpace(const Vector3D& location)const {

	Vector3D ll = location;
	const Vector3D p = position;
	ll.addScaledVector(p,1);
	return ll;

}

Vector3D RigidBody::getPointInLocalSpace(const Vector3D& location)const {

	Vector3D ll = location;
	const Vector3D p = position;
	ll.addScaledVector(p,-1);
	return ll;

}

Vector3D RigidBody::getDirectionInGlobalSpace(const Vector3D& location)const {

	Vector3D dg = getPointInGlobalSpace(location);
	dg.normalization();
	return dg;

}

Vector3D RigidBody::getDirectionInLocalSpace(const Vector3D& location)const {

	Vector3D dl = getPointInLocalSpace(location);
	dl.normalization();
	return dl;

}

