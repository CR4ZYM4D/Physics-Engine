#include "Included Files/body.h"
#include <assert.h>
#include <memory.h>

// helper function to find the inverseInertiaTensorGlobal (iitG) matrix of a body using 
// the inverseInertiaTensorLocal (iitL) matrix,  a 3x4 rotation matrix rM and a orientation quaternion 
static inline void transformInertiaTensor(const Matrix3by3& iitL ,const Quaternion& orientation,const Matrix3by4& rM , Matrix3by3 iitG){



}

// helper function to find the 3x4 transformation matrix of a body using its position vector and rotation quaternion 
static inline void calculateTransformMatrix(Matrix3by4& transformMatrix , const Vector3D& position , const Quaternion& rotation){
	transformMatrix.setOrientationAndPosition(rotation, position);
}



