#pragma once

#ifndef Particle_H //#if not defined function (#if function blocks are executed by preprocessor before compiler can look at them)
#define Particle_H
#include "fundamentals.h"
class Particle {
	
	protected :

		Vector3D position ; 
		Vector3D velocity ;
		Vector3D acceleration ; 
		Vector3D netForce; 
		real inverseMass ; //holds 1/mass of particle
		real drag ; // adding this to simulate the force of friction or drag which would prevent objects from moving infinitely 
		// drag is made real data type to use it as a multiplicative coefficient & prevent things from becoming  too complex
		// range of drag would be from 0 to 1 with 1 meaning frictionless or dragless

	public:

		void integrate(real duration);

		void setMass(const real Mass);

		real getMass() const;

		void setInverseMass(const real  inverseMass);

		real getInverseMass() const; // return inverseMass

		bool hasFiniteMass()const; // checkks if projectile is massless

		void setDrag(const real drag);

		real getDrag() const;

		void setPosition(const Vector3D& position);

		void setPosition(const real x, const real y, const real z);

		Vector3D getPosition() const;

		void getPosition(Vector3D* position)const;

		void setVelocity(Vector3D& velocity);

		void setVelocity(const real x, const real y, const real z);

		Vector3D getVelocity() const;

		void getVelocity(Vector3D* velocity)const;

		void setAcceleration(Vector3D& acceleration); // set acceleration by giving vector

		void setAcceleration(const real x, const real y, const real z); // by giving seperate axes values

		Vector3D getAcceleration() const; // return acceleration of particle

		void clearAccumulator(); // set total forces acting on particle to 0 

		void addForce(const Vector3D& force); // add a force acting on particle

};
#endif //used to end #if block