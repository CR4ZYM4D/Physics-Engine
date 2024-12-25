#include "Included files/particles.h"
class Ballistic {
	
	enum shotType { //enum is a special data tyoe which can onky have values whch are defined inside of it
		
		Unused = 0, 
		Bullet, //these are all available type of projectiles
		Shell,
		Ray,
		Fireball
	
	};

	struct AmmoRound { //structure defining a single ammo round
		
		Particle particle;
		shotType type;
		unsigned startTime;

		void render(){} //function to render the single ammo round

	};

	const static unsigned totalRounds = 20; // total ammo capacity
	
	AmmoRound ammo[totalRounds]; 

	shotType currType;

	void fire();

	public :
		
		Ballistic(); //create new ballistic object
		
		virtual const char* getTitle();

		virtual void update(); // update particle position

		virtual void display();

		virtual void mouse(int button, int state, int x, int y); //to handle mouse click

		virtual void key(unsigned char key); // to handle key press

};

Ballistic::Ballistic()  { 
	
	for (AmmoRound* shot = ammo; shot < ammo + totalRounds; shot++) { //making entire mag unused
	
		shot->type = Unused;
	
	}

}

const char* Ballistic :: getTitle(){
	
	return "Projectile demo" ;

}

void Ballistic :: fire() {

	AmmoRound* shot;

	for (shot = ammo; shot < ammo + totalRounds; shot++) {
		if (shot->type == Unused)	
			break;
	}

	if (shot >= ammo + totalRounds)
		return;
	
	switch ( currType ){
		
		case Bullet :

			shot->particle.setMass(2.0f);
			shot->particle.setVelocity(25, 0, 0);
			shot->particle.setAcceleration(0,-0.5,0);
			shot->particle.setDrag(0.99f);
			break;

		case Shell : 

			shot->particle.setMass(200);
			shot->particle.setVelocity(40,30,0);
			shot->particle.setAcceleration(0 , -10 , 0);
			shot->particle.setDrag(0.99f);
			break;

		case Fireball :

			shot->particle.setMass(1);
			shot->particle.setVelocity(10, 0, 0);
			shot->particle.setAcceleration(0, 0.5, 0);
			shot->particle.setDrag(0.9f);
			break;

		case Ray :

			shot->particle.setMass(0.1f);
			shot->particle.setVelocity(100, 0, 0);
			shot->particle.setAcceleration(0, 0, 0);
			shot->particle.setDrag(0.99f);
			break;

	}

	shot->particle.setPosition(0, 1.6f, 0);
	shot->type = currType;
	shot->startTime = Timings::get().lastFrameTimeStamp;
	shot->particle.clearAccumulator();
}

void Ballistic::update() {
	
	for (AmmoRound* shot = ammo; shot <= ammo + totalRounds; shot++) {
		
		if (shot->type != Unused) {
		
			shot->particle.integrate(duration);
			if (shot->particle.getPosition().y < 0.0f || shot->particle.getPosition().x >= 400.0f) {
			
				shot->type = Unused;
			
			}
	
		}
		
		break;
	
	}

}

void Ballistic :: display()

void Ballistic::mouse(int button, int state, int x, int y)

void Ballistic :: key(unsigned char key){
	
	switch (key) {
	
		case 1: currType = Bullet;
		break;
	
		case 2: currType = Shell;
			break;
	
		case 3: currType = Fireball;
			break;
	
		case 4: currType = Ray;
			break;
	}

}