//#include "Included files/bodyforcegenerator.h"
//
////function to update the forces acting on each body in the registry 
//void ForceRegistry::updateForces(real duration) {
//
//	BodyRegister::iterator i = body_registrations.begin();
//	for (; i < body_registrations.end(); i++) {
//		i->fg->updateForce(i->body, duration);
//	}
//
//}
//
////function to add a new body and its force .
//
//void ForceRegistry::add(RigidBody* body, BodyForceGenerator* bf) {
//
//	ForceRegistry::BodyRegistry registration;
//	registration.body = body;
//	registration.fg = bf;
//	body_registrations.push_back(registration);
//
//}
//
//
//
