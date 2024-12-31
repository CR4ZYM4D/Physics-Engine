#pragma once

#ifndef BODY_H
#define BODY_H
#include "fundamentals.h"

class RigidBody { //class to store the methods and values of a body

	protected:
		
		real inverseMass; // holds 1/mass of the body to use directly to find acceleration

		//holds the inverse of the inertia tensor w.r.t the body's local axes to calculate 
		//the angular acceleration. The inertia Tensor is 3x3 matrix 
		// containing the moment of inertia of the body about each axis 
		Matrix3by3 inverseInertiaTensorLocal;
		
		real linearDamping;

		real angularDamping; // holds the coefficient or constant of damping/drag forces 
							//to the rotational forces acting on the body

		Vector3D position;

		Quaternion orientation; // holds the rotation/orientation of the body in 3 dimensions

		Vector3D linearVelocity;

		Vector3D angularVelocity;

		Matrix3by3 inverseInertiaTensorGlobal; //holds the inverse inertia tensor of the body w.r.t. 
											  // the global axes system

		Matrix3by4 transformMatrix; //transform matrix that will be used to interchange in between 
									//the global and local coordinates of some point in the body
		
		Vector3D forceAccumulated; // stores the total force accumulated on the body

		Vector3D torqueAccumulated; // stores the total torque accumulated on the body

		Vector3D acceleration; 

		Vector3D previousFrameAcceleration; // stores the linear acceleration of the body
										   //from the previous frame
		
		public:

			//function to calculate the transfom and inertia tensor 
			//matrices after the end of frame to get their updated values
			void calculateDerivedData();

			// linear function to integrate te acceleration and 
			//force acting on the body going into the next frame i.e.
			// it integrates the body forces and vectors forward in time using Newton's laws
			void integrate(real duration);

			void setMass(const real mass); // set the mass of the current body

			real getMass()const; //return the mass of the current body

			void setInverseMass(const real inverseMass);// set the inverse mass of the body

			real getInverseMass()const; //return the inverse mass of the current body

			bool hasFiniteMass()const; //checks if the current body does not have infinite mass

			void setInertiaTensor(const Matrix3by3& inertiaTensor); // sets the inertia tensor of the 
																	// current body as the parameter matrix
			
			// copies the inverse of the inverseInertiaTensorLocal matrix of the current body 
			// into the matrix whose pointer is passed as parameter 
			void getInertiaTensor(Matrix3by3* inertiaTensorLocal)const;

			Matrix3by3 getInertiaTensor() const; //return a new 3x3 matrix of the body's local inertia tensor

			//copies the global inertia tensor matrix of the body into the matrix  
			// whose pointer is passed as parameter
			void getInertiaTensorGlobal(Matrix3by3* inertiaTensorGlobal) const;

			Matrix3by3 getInertiaTensorGlobal() const; //return a new 3x3 matrix with the global 
			//inertia tensor of the body

			//sets the inverseInertiaTensorLocal matrix of the current body
			void setInverseInertiaTensorLocal(const Matrix3by3& inverseInertiaTensorLocal);

			void getInverseInertiaTensorLocal(Matrix3by3* inverseInertiaTensorL) const;

			Matrix3by3 getInverseInertiaTensorLocal()const; //return the inverseInertiaTensorLocal 
														   //of the body in a new matrix
			
			//sets the inverseInertiaTensorGlobal matrix of the current body
			void setInverseInertiaTensorGlobal(const Matrix3by3& inverseInertiaTensorGlobal);

			void getInverseInertiaTensorGlobal(Matrix3by3* inverseInertiaTensorG) const;

			Matrix3by3 getInverseInertiaTensorGlobal()const; //return the inverseInertiaTensorGlobal 
														   //of the body in a new matrix

			//similar as above functions to initialize or return the linearDamping and angularDamping coefficients
			void setDamping(const real linearDamping, const real angularDamping);

			void setLinearDamping(const real linearDamping);

			real getLinearDamping() const;
			
			void setAngularDamping(const real angularDamping);

			real getAngularDamping() const;

			//similar as above functions to initialize or return the body's postion
			void setPosition(const Vector3D& position);

			void setPosition(const real x, const real y ,const real z);

			Vector3D getPosition() const;
			
			void getPosition(Vector3D* position)const;

			//similar as above functions to initialize or return the body's orientation/rotation
			void setOrientation(const Quaternion& orientation);

			void setOrientation(const real r , const real i, const real j ,const real k);

			Quaternion getOrientation() const;
			
			void getOrientation(Quaternion* orientation)const;

			//fills the parameter matrix by the orientation of the body by transforming the orientation quaternion
			void getOrientation(Matrix3by3* orientationMatrix) const;

			//fills the parameter array by the orientation of the body by transforming the orientation quaternion
			void getOrientation(real matrix[9])const;

			void getTransform(Matrix3by4* transform)const; //fills the parameter matrix with the orientation 
														  //transform and position of the body

			void getTransform(real matrix[16])const;//fills the array with the orientation transform and position of the body

			Matrix3by4 getTransform()const;

			void getGLTransform(float matrix[16])const; //fills the matrix with the transpose of the matrix obtained by the 
														//getTransform method to make it suitable for use by OpenGL

			Vector3D getPointInLocalSpace(const Vector3D& location) const; //converts the coordinates of a point from the 
																		  //global coordinates to the local coordinates w.r.t. the body's center of mass 

			Vector3D getPointInGlobalSpace(const Vector3D& location)const; //converts the coordinates of a point w.r.t the body's local coordinates 
																		  //to the global coordinates

			Vector3D getDirectionInLocalSpace(const Vector3D& direction) const; //converts the direction of a point w.r.t global coordinates 
																				//into the body's local coordinates

			Vector3D getDirectionInGlobalSpace(const Vector3D& direction) const; //converts the direction of a point w.r.t body's local coordinates 
																				//into the global coordinates

			//functions to initialize/return or add the linear and angular velocities of the body
			void setLinearVelocity(const Vector3D& linearVelocity);

			void setLinearVelocity(const real x, const real y, const real z);

			void getLinearVelocity(Vector3D* linearVelocity)const;

			Vector3D getLinearVelocity()const;

			void addLinearVelocity(const Vector3D& deltaLinearVelocity);



			void setAngularVelocity(const Vector3D& angularVelocity);

			void setAngularVelocity(const real x, const real y, const real z);

			void getAngularVelocity(Vector3D* angularVelocity)const;

			Vector3D getAngularVelocity()const;

			void addAngularVelocity(const Vector3D& deltaAngularVelocity);

			// fills the given vector with the total linear acceleration on the body in the previous frame
			// the acceleration is set w.r.t global coordinates 
			void getPreviousFrameAcceleration(Vector3D* acceleration)const;

			Vector3D getPreviousFrameAccleration()const; //return the linear acceleration of the body in the 
														 //end of the previous frame in a new vector

			void clearAccumulators();//clear the force and torque accumulators

			void addForce(const Vector3D& force); // add a force acting in the center of mass of the body
			//the force vector is in global coordinates

			//add the force which is acting at some point on the body , not necessarily the center of mass. 
			// Thus contributing torque as well. Both the force and position vectors are w.r.t. global coordinates
			void addForceAtGlobalPoint(const Vector3D& force, const Vector3D& position);

			//add the force acting at some point on the body. the force vector is w.r.t. global coordinates while 
			// the position vector is w.r.t. local coordinates of the body
			void addForceAtLocalPoint(const Vector3D& force, const Vector3D& position);

			void addTorque(const Vector3D& torque);//adds the given torque acting on the body

			//functions to initialize/return the acceleration
			void setAcceleration(const Vector3D& acceleration);

			void setAcceleration(const real x, const real y, const real z);

			void getAcceleration(Vector3D* acceleration) const;

			Vector3D getAcceleration()const;
};

#endif