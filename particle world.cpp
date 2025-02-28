//#include <cstddef>
//#include "Included files/particle world.h"
//
//ParticleWorld::ParticleWorld(unsigned maxContacts, unsigned iterations) : resolver(iterations), maxContacts(maxContacts) {
//
//	contacts = new ParticleContacts[maxContacts];
//	calculateIterations = (iterations == 0);
//
//}
//
//ParticleWorld :: ~ParticleWorld() {
//
//	delete[] contacts;
//
//}
//
//void ParticleWorld::startFrame() {
//
//	for (ParticleList::iterator p = plist.begin(); p != plist.end(); p++) {
//		(*p)->clearAccumulator();
//	}
//
//}
//
//unsigned ParticleWorld::generateContacts() {
//
//	unsigned limit = maxContacts; 
//	ParticleContacts* nextContact = contacts;
//
//	for (ContactGeneratorList::iterator g = cgenlist.begin(); g != cgenlist.end(); g++) {
//		unsigned used = (*g)->addContact(nextContact, limit);
//		limit -= used;
//		nextContact += used;
//
//		if (limit <= 0)
//			break;
//	}
//	return maxContacts - limit;
//}
//
//void ParticleWorld::integrate(real duration) {
//
//	for (ParticleList::iterator p = plist.begin(); p != plist.end(); p++){
//		(*p)->integrate(duration);
//	}
//
//}
//
//void ParticleWorld::runPhysics(real duration) {
//
//	//first apply all force generators
//	registry.updateForces(duration);
//
//	integrate(duration);
//
//	unsigned usedContacts = generateContacts();
//
//	if (usedContacts) {
//		if (calculateIterations)
//			resolver.setIterations(usedContacts * 2);
//		resolver.resolveContects(contacts, usedContacts, duration);
//	}
//
//}
//
//ParticleWorld::ParticleList& ParticleWorld ::getParticles() {
//
//	return plist;
//
//}
//
//ParticleWorld::ContactGeneratorList& ParticleWorld::getContactGenerators() {
//
//	return cgenlist;
//
//}
//
//ForceRegistry& ParticleWorld::getForceRegistries() {
//
//	return registry;
//
//}
//
//void GroundContacts::init(ParticleWorld::ParticleList* plist) {
//
//	GroundContacts::plist = plist;
//
//}
//
//unsigned GroundContacts::addContact(ParticleContacts* contacts, unsigned limit)const {
//
//	unsigned count = 0;
//	for (ParticleWorld::ParticleList::iterator p = plist->begin(); p != plist->end(); p++) {
//
//		real y = (*p)->getPosition().y;
//		if(y<0.0f){
//			contacts->contactnormal = Vector3D(0.0f, 1.0f, 0.0f);
//			contacts->particles[0] = *p;
//			contacts->particles[1] = NULL;
//			contacts->penetration = -y;
//			contacts->coeffrestitute = 0.2f;
//			contacts++;
//			count++;
//		}
//		if (count >= limit)
//			return count;
//	}
//	return count;
//}