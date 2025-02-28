//#pragma once
//#ifndef BODY_FORCE_GENERATOR
//#define BODY_FORCE_GENERATOR
//#include "body.h"
//#include "particleforcegenerator.h"
//#include <vector>
//
//class BodyForceGenerator { //interface for the force generator of the body similar to that of the particle
//
//	public:
//		virtual void updateForce(RigidBody* body, real duration);
//
//};
//
//class Gravity : BodyForceGenerator {
//
//	Vector3D gravity;
//	public:
//		Gravity(const Vector3D& gravity);
//		virtual void updateForce(RigidBody* body, real duration);
//
//};
//
//class Spring : BodyForceGenerator {
//
//	Vector3D connection_point;//stores the global coordinates of the point at which 
//	//the spring is connected to the body 
//	Vector3D other_connection_point;//stores the connection point of the other body
//	RigidBody* otherBody;//stores the pointer of the other body
//	real spring_constant;//stores the spring constant
//	real rest_length;//stores the rest length of the spring
//
//	public:
//		Spring(Vector3D& local_connection_point , Vector3D& other_connection_point,
//			RigidBody* other, real spring_constant, real rest_length);
//
//		virtual void updateForce(RigidBody* body , real duration);
//};
//
//class Buoyancy : BodyForceGenerator {
//
//	real volume; //stores the volume of the body
//	real body_height; //stores the height of the body to calculte the max buoyant force
//	Vector3D center_of_buoyancy;//stores the center of buoyncy for the body
//	real liquid_elevation;//stores the global height till which the liquid is present 
//	real liquid_density;//stores the density of the liquid 
//
//	public:
//		Buoyancy(real v , real b_h , Vector3D& c_o_b , real l_e , real l_d = 1000.0f);
//
//		virtual void updateForce(RigidBody* body ,real duration);
//
//};
//
//class ForceRegistry {
//
//	private:
//		struct BodyRegistry {
//
//			RigidBody* body;
//			BodyForceGenerator* fg;
//
//		};
//		typedef std::vector<BodyRegistry> BodyRegister;
//		BodyRegister body_registrations;
//
//	public:
//
//		void add(RigidBody* body, BodyForceGenerator* fg);
//		void remove(RigidBody* body, BodyForceGenerator* fg);
//		void clear();
//		void updateForces(real duration);
//};
//
//class Explosion : public BodyForceGenerator , public ForceGenerator{
//
//	private:
//		
//		real min_radius;
//
//		real max_radius;
//
//		Vector3D detonation_point;
//
//		real explosion_duration;
//
//		real time_passed;
//
//		real explosion_force;
//
//		real tower_height;
//
//		real shockwave_speed;
//
//		real shockwave_thickness;
//
//		real concussion_force;
//
//		real concussion_duration;
//
//		real tower_radius;
//
//		real max_convection_force;
//
//		real convection_duration;
//
//	public:
//
//		Explosion();
//
//		virtual void updateForce(RigidBody* body, BodyForceGenerator* bfg);
//
//		virtual void updateForce(Particle* particle, ForceGenerator* fg);
//
//};
//
//
//
//#endif