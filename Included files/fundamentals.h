#pragma once
#include "accuracy.h"
#ifndef FUNDAMENTALS_H
#define FUNDAMENTALS_H

class Vector3D { // class to create vector object
public:
	real x;
	real y;
	real z;

	Vector3D() { // default constructor to make zero vector
		x = 0;
		y = 0;
		z = 0;
	}

	Vector3D(const real ax, const real ay, const real az) { // parameterized constructor
		x = ax;
		y = ay;
		z = az;
	}

	void invert() { //inverses a vector
		x *= -1;
		y *= -1;
		z *= -1;
	}

	real magnitude() const {
		return real_sqrt(x * x + y * y + z * z);
	}

	real magnitudeSquare() const {
		return (x * x + y * y + z * z);
	}

	void normalization() {
		(*this) *= (real)(1 / magnitude());
	}

	void operator*= (const real value) { /*overloading the *= operator so that we constantly dont need to multiply each value of
										 vector individually but just simply completely multiply the vector at once instead*/
		x *= value;
		y *= value;
		z *= value;
	}

	Vector3D operator*(const real value) const { //overloading the * operator as well for scaling returning a new scaled vector
		return Vector3D(x * value, y * value, z * value);
	}

	// similarly, we overload the + , - , += , -= operators as well for vector addition and subtraction

	void operator+= (const Vector3D& v) { //overloading the += operator
		x += v.x;
		y += v.y;
		z += v.z;
	}

	Vector3D operator+(Vector3D& v) const { //overloading the + operator 
		return Vector3D(x + v.x, y + v.y, z + v.z);
	}

	Vector3D operator+(Vector3D v) { //overloading the + operator 
		return Vector3D(x + v.x, y + v.y, z + v.z);
	}

	void operator-= (const Vector3D& v) { //overloading the -= operator
		x -= v.x;
		y -= v.y;
		z -= v.z;
	}

	Vector3D operator -=(Vector3D& v) {
		return Vector3D(x - v.x, y - v.y, z - v.z);
	}

	Vector3D operator-(Vector3D& v) { //overloading the - operator 
		return Vector3D(x - v.x, y - v.y, z - v.z);
	}

	Vector3D operator-(const Vector3D& v) { //overloading the - operator 
		return Vector3D(x - v.x, y - v.y, z - v.z);
	}

	Vector3D componentMultiply(Vector3D& v) { //function for component multiplication of two vectors returned in third vector
		//if there are two vectors v1 and v2 then the result of their component multiplication would be v1x*v2x , v1y*v2y , v1z*v2z 
		return Vector3D(x * v.x, y * v.y, z * v.z);
	}

	void componentMultiplier(Vector3D& v) {//component multiplication and stored in first vector
		x *= v.x;
		y *= v.y;
		z *= v.z;
	}

	real operator *=(const Vector3D& v)const { //overloading *= & * operator to perform scalar product upon giving vector parameter
		return x * v.x + y * v.y + z * v.z;
	}

	real operator *(const Vector3D& v) const{
		return x * v.x + y * v.y + z * v.z;
	}

	Vector3D X(Vector3D& v) { //function for vector product
		return Vector3D(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
	}

	Vector3D X(const Vector3D& v) { //function for vector product
		return Vector3D(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
	}

	void addScaledVector(Vector3D& v, real scale) {
		x += v.x * scale;
		y += v.y * scale;
		z += v.z * scale;
	}
	
	void addScaledVector(const Vector3D& v,const real scale) {
		x += v.x * scale;
		y += v.y * scale;
		z += v.z * scale;
	}

	void clear() {
		x = y = z = 0;
	}
};


class Quaternion { //will create a quaternion that will be used for rotations both partial and complete

	//our Quaternion will be of the form of a 4x1 matrix
	
	 //type of data structure in C/C++ in which all particles have same memory location
		public:

		union{
			struct {
				real r; //stores real part of the quaternion

				//the rest 3 will store the complex part of the quaternion
				real i;

				real j;

				real k;
			};
			real matrix[4]; //will store the quaternion as an array
		};
	

	Quaternion() : r(1) , i(0) , j(0) , k(0){} //default constructor which will create a quaternion with zero rotation

	//parameterized constructor given below

	Quaternion(real r , real i , real j , real k) : r(r) , i(i) , j(j) , k(k){}

	void normalise() { //used to normalise a quaternion to make its magnitude 1, similar to vectors

		real d = r*r + i*i + j*j + k*k;

		d = (real)(real(1.0) / real_sqrt(d));

		r *= d;
		i *= d;
		j *= d;
		k *= d;
	}

	void quaternionMultiply(const Quaternion &second){//will be used to multiply 2 quaternions which as a result will produce some rotation or transformation in the object
	
		Quaternion first = *this;
		r = first.r * second.r - first.i - second.i - first.j * second.j - first.k * second.k;
	
		i = first.r * second.i + second.r * first.i + first.j * second.k - first.k * second.j;

		j = first.r * second.j + first.j * second.r + first.i * second.k - first.k * second.i;

		k = first.r * second.k + first.k * second.r + first.i * second.j - first.j * second.i;
	
	}

	void addScaledVector(const Vector3D& v, real scale) {

		Quaternion q(0, v.x*scale , v.y*scale , v.z*scale); //creating a quaternion of given scaled vector

		q.quaternionMultiply(*this);
	
		r += q.r * (real(0.5));  //the 0.5 is done because of the theory of quaternion with the components representing sin and cos of half the total angle of rotation
	
		i += q.i * (real(0.5));

		j += q.j * (real(0.5));

		k += q.k * (real(0.5));
	}

	void rotateByVector(const Vector3D& v) {//to rotate a quaternion by a vector we simply need to multiply the two

		Quaternion q(0, v.x, v.y, v.z);

		(*this).quaternionMultiply(q);

	}
};


//creating a 3x4 matrix for rotations and transformations we will also use it as a 4x4 matrix sometimes by adding a 4th row as [0 0 0 1] 

//in each matrix the columns represent where the unit axes would point after a rotation or a transformation with the first, second and third columns representing the X, Y and Z axis respectively

class Matrix3by4{
	//in our 3x4 matrix, we assume that the last column of elements to be the positional offset part 
	public :
		real data[12]; //we create a 1D array so that things remain simpler to manage

		Matrix3by4() {//default constructor used to initialize first 3 rows of a 4x4 identity matrix
			
			data[0] = data[5] = data[10] = 1;
			data[1] = data[2] = data[3] = data[4] = data[6] = data[7] = data[8] = data[9] = data[11] = 0;

		}
		
		Matrix3by4 multiplyMatrices(const Matrix3by4& o) const { //multiplying two 3x4 matrices and storing it in a third matrix

			Matrix3by4 result;
			result.data[0] = data[0] * o.data[0] + data[1] * o.data[4] + data[2] * o.data[8];
			result.data[4] = data[4] * o.data[0] + data[5] * o.data[4] + data[6] * o.data[8];
			result.data[8] = data[8] * o.data[0] + data[9] * o.data[4] + data[10] * o.data[8];

			result.data[1] = data[0] * o.data[1] + data[1] * o.data[5] + data[2] * o.data[9];
			result.data[5] = data[4] * o.data[1] + data[5] * o.data[5] + data[6] * o.data[9];
			result.data[9] = data[8] * o.data[1] + data[9] * o.data[5] + data[10] * o.data[9];

			result.data[2] = data[0] * o.data[2] + data[1] * o.data[6] + data[2] * o.data[10];
			result.data[6] = data[4] * o.data[2] + data[5] * o.data[6] + data[6] * o.data[10];
			result.data[10] = data[8] * o.data[2] + data[9] * o.data[6] + data[10] * o.data[10];

			result.data[3] = data[0] * o.data[3] + data[1] * o.data[7] + data[2] * o.data[11] + data[3];
			result.data[7] = data[4] * o.data[3] + data[5] * o.data[7] + data[6] * o.data[11] + data[7];
			result.data[11] = data[8] * o.data[3] + data[9] * o.data[7] + data[10] * o.data[11] + data[11];
		
			return result;
		}

		Vector3D multiplyByVector(const Vector3D& v)const { //multiplying the matrix by given vector and returning result as vector

			return Vector3D(

				data[0] * v.x + data[1] * v.y + data[2] * v.z +data[3] ,
				
				data[4] * v.x + data[5] * v.y + data[6] * v.z  + data[7],

				data[8] * v.x + data[9] * v.y + data[10] * v.z +data[11]

			);

		}

		Vector3D transformByVector(const Vector3D& v)const { //transforming the given vector by matrix and returning a 
			return (*this).multiplyByVector(v);
		}

		real getDeterminant() const { //calculates the determinant of the given matrix
			return data[0]*data[5]*data[10]-
				data[0] * data[6] * data[9]-
				data[1] * data[4] * data[10]+
				data[1] * data[6] * data[8]+
				data[2] * data[4] * data[9]-
				data[2] * data[5] * data[8];
		}

		void setInverse(const Matrix3by4& m) { //sets the current matrix as the inverse of the parameter matrix
			
			real delta = m.getDeterminant();
			if (delta == 0)
				return;
			delta = real(1.0) / delta;

			data[0] = (m.data[5] * m.data[10] - m.data[6] * m.data[9]) * delta;
			data[4] = (m.data[6] * m.data[8] - m.data[4] * m.data[10]) * delta;
			data[8] = (m.data[4] * m.data[9] - m.data[5] * m.data[8]) * delta;

			data[1] = (m.data[2] * m.data[9] - m.data[1] * m.data[10]) * delta;
			data[5] = (m.data[0] * m.data[10] - m.data[2] * m.data[8]) * delta;
			data[9] = (m.data[1] * m.data[8] - m.data[0] * m.data[9]) * delta;

			data[2] = (m.data[1] * m.data[6] - m.data[2] * m.data[5]) * delta;
			data[6] = (m.data[2] * m.data[5] - m.data[0] * m.data[6]) * delta;
			data[10] = (m.data[0] * m.data[5] - m.data[1] * m.data[4]) * delta;

			data[3] = (m.data[1] * m.data[7] * m.data[10] + m.data[3] * m.data[6] * m.data[9] 
				+m.data[2]*m.data[5]*m.data[11] - m.data[1]*m.data[6]*m.data[11]
				- m.data[3]*m.data[5]*m.data[10] - m.data[2]*m.data[9]*m.data[7]) * delta;

			data[7] = (m.data[0] * m.data[6] * m.data[11] - m.data[0] * m.data[7] * m.data[10] 
				-m.data[2]*m.data[4]*m.data[11] + m.data[2]*m.data[7]*m.data[8]
				+ m.data[3]*m.data[4]*m.data[10] - m.data[3]*m.data[6]*m.data[8]) * delta;		
			
			data[11] = (m.data[0] * m.data[7] * m.data[9] + m.data[1] * m.data[4] * m.data[11] 
				+m.data[3]*m.data[5]*m.data[8] - m.data[0]*m.data[5]*m.data[11]
				- m.data[1]*m.data[7]*m.data[8] - m.data[3]*m.data[9]*m.data[4]) * delta;


		}

		Matrix3by4 inverse ()const { //function that would return the inverse of a 3x4 matrix in a new 3x4 matrix
			Matrix3by4 result;
			result.setInverse(*this);
			return result;
		}

		void invert() { //function that would inverse the original 3x4 matrix and store it in the same matrix
			setInverse(*this);
		}

		Vector3D transformByMatrix(Vector3D& v) { //function to transform the given direction vector by the 3x4 matrix a direction vector would be written as a 4x1 vector of [v.x v.y v.z 0]

			return Vector3D(
				data[0] * v.x + data[1] * v.y + data[2] * v.z ,

				data[4] * v.x + data[5] * v.y + data[6] * v.z ,

				data[8] * v.x + data[9] * v.y + data[10] * v.z 
			);

		}

		Vector3D transformUnitByInverseMatrix(Vector3D &v){ //transform the given direction vector by the inverse of the matrix assuming that the matrix is purely a rotation matrix which means that its inverse would be the same as its transform
		

			return Vector3D(

				v.x * data[0] + v.y * data[4] + v.z * data[8],
				v.x * data[1] + v.y * data[3] + v.z * data[7],
				v.x * data[2] + v.y * data[6] + v.z * data[10]

			);
		}

		Vector3D transformByInverseMatrix(Vector3D& v) const{ //transform the given vector in the same way as above but this time we first convert it into a directional vector by removing the offset part by subtracting the last column of the matrix from the vector components

			Vector3D tmp = v;
			tmp.x -= data[3];
			tmp.y -= data[7];
			tmp.z -= data[11];
		
			return Vector3D(

				tmp.x * data[0] + tmp.y * data[4] + tmp.z * data[8],
				tmp.x * data[1] + tmp.y * data[3] + tmp.z * data[7],
				tmp.x * data[2] + tmp.y * data[6] + tmp.z * data[10]

			);

		}

		Vector3D getAxisVector(int i) { //returns the vector along which the given axis is currently pointing at

			return Vector3D(data[i], data[i + 4], data[i + 8]);

		}

		void setOrientationAndPosition(const Quaternion& q, const Vector3D& v) { //used to convert the data of a rotation/orientation quaternion and position vector into a 3x4 matrix

			//assuming that the quaternion is a 4x1 matrix as [r x y z]
			//and vector is as [x y z] converting it into a 3x4 matrix by multiplying the quaternion 
			//with its transpose for the first 3 columns of the matrix and putting the 
			//	vector as the last column for the position
		
			data[0] = 1 - (2 * q.j * q.j + 2 * q.k * q.k);
			data[1] = 2 * q.i * q.j + 2 * q.k * q.r;
			data[2] = 2 * q.i * q.k - 2 * q.j * q.r;
			data[3] = v.x;

			data[4] = 2 * q.i * q.j - 2 * q.k * q.r;
			data[5] = 1 - (2 * q.i * q.i + 2 * q.k * q.k);
			data[6] = 2 * q.j * q.k + 2 * q.i * q.r;
			data[7] = v.y;

			data[8] = 2 * q.i * q.k + 2 * q.j * q.r;
			data[9] = 2 * q.j * q.k - 2 * q.i * q.r;
			data[10] = 1 - (2 * q.i * q.i + 2 * q.j * q.j);
			data[11] = v.z;

		}

		void fillGLArray(float array[16])const { //to fill the data of the given matrix into a 4x4 array for rendering by openGL

			//keep in mind that openGL uses a column major format while we were operating in a row major format,
			// so we will have to fill the array as the transpose of our original one 
			// the last column will be set as [0 0 0 1] as mentioned above to convert the 3x4 matrix to 4x4 

			array[0] = (float)data[0];
			array[4] = (float)data[1];
			array[8] = (float)data[2];
			array[12] = (float)data[3];

			array[1] = (float)data[4];
			array[5] = (float)data[5];
			array[9] = (float)data[6];
			array[13] = (float)data[7];

			array[2] = (float)data[8];
			array[6] = (float)data[9];
			array[10] = (float)data[10];
			array[14] = (float)data[11];

			array[3] = (float)0;
			array[7] = (float)0;
			array[11] = (float)0;
			array[15] = (float)1;

		}

};

class Matrix3by3 { //will be used to hold the inertia tensor of an object/body and thus used to calculate rotations
	//angular momentum  and angular velocity

	public:
		real data[9];
		
		Matrix3by3(){ //default constructor creates a nul matrix

			data[0] = data[1] = data[2] = data[3] = data[4] = data[5] = data[6] = data[7] = data[8] = 0;

		}

		Matrix3by3(const Vector3D& v1 , const Vector3D& v2 , const Vector3D& v3) {
			//parametrized constructor with each vector representing one column of the matrix

			data[0] = v1.x;
			data[3] = v1.y;
			data[6] = v1.z;

			data[1] = v2.x;
			data[4] = v2.y;
			data[7] = v2.z;

			data[2] = v3.x;
			data[5] = v3.y;
			data[8] = v3.z;

		}

		Matrix3by3(real c0, real c1, real c2, real c3, real c4, real c5, real c6, real c7, real c8) {
			//parameterized constructor with each element given individually 

			data[0] = c0; data[1] = c1; data[2] = c2;
			data[3] = c3; data[4] = c4; data[5] = c5;
			data[6] = c6; data[7] = c7; data[8] = c8;

		}

		//to create a diagonal matrix with the elements of the major diagonal given to you 

		void setDiagonal(real a, real b, real c) {
			setInertiaTensorCoeffs(a, b, c);
		}

		//to create the inertia tensor matrix for rotation with given values, inertia tensors are represented using diagonal matrices


		void setInertiaTensorCoeffs(real ix, real iy, real iz, real ixy = 0, real ixz = 0, real iyz = 0) {

			data[0] = ix;
			data[4] = iy;
			data[8] = iz;

			data[1] = data[3] = -ixy;
			data[2] = data[6] = -ixz;
			data[5] = data[7] = -iyz;

		}


		//set the value of the matrix as the inertia tensor of a cuboidal block with the vector
		//half denoting the half length of each side i.e. each component of the vector half denotes 
		//the distance of the face of cuboid lying perpendicular to that particular from the centre of the cuboid

		void setBlockInertiaTensor(Vector3D& half, real mass) {

			Vector3D sq = half.componentMultiply(half);
			setInertiaTensorCoeffs(0.3f * mass * (sq.y + sq.z), 0.3f * mass * (sq.x + sq.z), 0.3f * mass * (sq.y + sq.x));

		}


		//converts a vector into  a skew symmetric matrix,
		//multiplying a skew matrix with a second vector will give the same
		//result as doing the cross product of the two vectors

		void setSkewMatrix(const Vector3D v) {

			data[0] = data[4] = data[8] = 0;
			data[1] = -v.z;
			data[2] = v.y;
			data[3] = v.z;
			data[5] = -v.x;
			data[6] = -v.y;
			data[7] = v.x;

		}


		//sets the components of the matrix with 3 parameter 
		//vectors, each vector representing one column
		
		void setComponents(Vector3D& v1, Vector3D& v2, Vector3D& v3) {

			data[0] = v1.x;
			data[1] = v2.x;
			data[2] = v3.x;

			data[3] = v1.y;
			data[4] = v2.y;
			data[5] = v3.y;

			data[6] = v1.z;
			data[7] = v2.z;
			data[8] = v3.z;

		}


		//same functionality as multiplyByVector for the 3x4 matrix but for a 3x3 matrix instead

		Vector3D multiplyByVector(const Vector3D& v) const {

			return Vector3D(data[0] * v.x + data[1]*v.y + data[2]*v.z , 
				data[3] * v.x + data[4] * v.y + data[5] * v.z ,
				data[6] * v.x + data[7] * v.y + data[8] * v.z);

		}

		//same functionality as transformByVector for the 3x4 matrix but for a 3x3 matrix instead

		Vector3D transformByVector(const Vector3D& v) const {
			return (*this).multiplyByVector(v);
		}

		Vector3D transformByTranspose(const Vector3D& v)const { //returns a vector after transforming the given vector by the transpose of the matrix

			return Vector3D(data[0] * v.x + data[3] * v.y + data[6] * v.z,
				data[1] * v.x + data[4] * v.y + data[7] * v.z,
				data[2] * v.x + data[5] * v.y + data[8] * v.z);

		}

		Vector3D getRowVector(int i) const { //returns the i-1th row of the matrix 

			return Vector3D(data[i * 3], data[3 * i + 1], data[3 * i + 2]);
		
		}

		Vector3D getAxisVector(int i)const { //same functionality as that of the 4x4 matrix function
			return Vector3D(data[i], data[i + 3], data[i + 6]);
		}

		void setInverse(const Matrix3by3& m) { //inverts the matrix m and stores the inverted matrix in the current matrix

			real determinant = m.data[0] * (m.data[4] * m.data[8] - m.data[5] * m.data[7]) -
				m.data[1] * (m.data[3] * m.data[8] - m.data[5] * m.data[6]) +
				m.data[2] * (m.data[3] * m.data[7] - m.data[4] * m.data[6]);

			if (determinant == 0)
				return;

			real d = (real)(1.0) / determinant;

			data[0] = d * (m.data[4] * m.data[8] - m.data[5] * m.data[7]);
			data[1] = d * (m.data[2] * m.data[7] - m.data[1] * m.data[8]);
			data[2] = d * (m.data[1] * m.data[5] - m.data[2] * m.data[4]);
			data[3] = d * (m.data[5] * m.data[6] - m.data[3] * m.data[8]);
			data[4] = d * (m.data[0] * m.data[8] - m.data[2] * m.data[6]);
			data[5] = d * (m.data[2] * m.data[3] - m.data[0] * m.data[5]);
			data[6] = d * (m.data[3] * m.data[7] - m.data[4] * m.data[6]);
			data[7] = d * (m.data[1] * m.data[6] - m.data[0] * m.data[7]);
			data[8] = d * (m.data[0] * m.data[4] - m.data[1] * m.data[3]);

		}

		Matrix3by3 inverse() const { //returns a new matrix that stores the inverse of the current matrix

			Matrix3by3 res;
			res.setInverse(*this);
			return res;

		}

		void invert() { //inverse the current matrix and store it in this one itself
			setInverse(*this);
		}

		void setTranspose(const Matrix3by3& m) { //stores the transpose of the parameter matrix in the current one

			data[0] = m.data[0];
			data[4] = m.data[4];
			data[8] = m.data[8];
			data[1] = m.data[3];
			data[2] = m.data[6];
			data[3] = m.data[1];
			data[5] = m.data[7];
			data[6] = m.data[2];
			data[7] = m.data[5];

		}

		Matrix3by3 transpose() const { //return transpose of current matrix in a new matrix

			Matrix3by3 res;
			res.setTranspose(*this);
			return res;

		}

		Matrix3by3 multiplyMatrices(const Matrix3by3& o) { //multiplies given matrix with the current one and returns it in a new matrix 
			
				return Matrix3by3(
					data[0] * o.data[0] + data[1] * o.data[3] + data[2] * o.data[6],
					data[0] * o.data[1] + data[1] * o.data[4] + data[2] * o.data[7],
					data[0] * o.data[2] + data[1] * o.data[5] + data[2] * o.data[8],

					data[3] * o.data[0] + data[4] * o.data[3] + data[5] * o.data[6],
					data[3] * o.data[1] + data[4] * o.data[4] + data[5] * o.data[7],
					data[3] * o.data[2] + data[4] * o.data[5] + data[5] * o.data[8],

					data[6] * o.data[0] + data[7] * o.data[3] + data[8] * o.data[6],
					data[6] * o.data[1] + data[7] * o.data[4] + data[8] * o.data[7],
					data[6] * o.data[2] + data[7] * o.data[5] + data[8] * o.data[8]
				);
			
		}

		void multiplyMatrixInSelf(const Matrix3by3& o) { //multiplies current matrix with given matrix and 
														//stores the product in the current matrix

			real t1;
			real t2;
			real t3;

			t1 = data[0] * o.data[0] + data[1] * o.data[3] + data[2] * o.data[6];
			t2 = data[0] * o.data[1] + data[1] * o.data[4] + data[2] * o.data[7];
			t3 = data[0] * o.data[2] + data[1] * o.data[5] + data[2] * o.data[8];
			data[0] = t1;
			data[1] = t2;
			data[2] = t3;

			t1 = data[3] * o.data[0] + data[4] * o.data[3] + data[5] * o.data[6];
			t2 = data[3] * o.data[1] + data[4] * o.data[4] + data[5] * o.data[7];
			t3 = data[3] * o.data[2] + data[4] * o.data[5] + data[5] * o.data[8];
			data[3] = t1;
			data[4] = t2;
			data[5] = t3;

			t1 = data[6] * o.data[0] + data[7] * o.data[3] + data[8] * o.data[6];
			t2 = data[6] * o.data[1] + data[7] * o.data[4] + data[8] * o.data[7];
			t3 = data[6] * o.data[2] + data[7] * o.data[5] + data[8] * o.data[8];
			data[6] = t1;
			data[7] = t2;
			data[8] = t3;

		}

		void scale(const real scalar) { //scales or multiplies the entire matrix by the given number

			data[0] *= scalar; data[1] *= scalar; data[2] *= scalar;
			data[3] *= scalar; data[4] *= scalar; data[5] *= scalar;
			data[6] *= scalar; data[7] *= scalar; data[8] *= scalar;

		}
		void addInSelf(const Matrix3by3& o) { //similar to mulitplyInSelf but adds instead of multiplying

			data[0] += o.data[0]; data[1] += o.data[1]; data[2] += o.data[2];
			data[3] += o.data[3]; data[4] += o.data[4]; data[5] += o.data[5];
			data[6] += o.data[6]; data[7] += o.data[7]; data[8] += o.data[8];

		}

		void setOrientation(const Quaternion& q) { //used to convert the data of a rotation quaternion into the matrix  

			data[0] = 1 - 2 * (q.j * q.j + q.k * q.k);
			data[1] = 2 * (q.i * q.j + q.k * q.r);
			data[2] = 2 * (q.i * q.k - q.j * q.r);
			data[3] = 2 * (q.i * q.j - q.k * q.r);
			data[4] = 1 - 2 * (q.i * q.i + q.k * q.k);
			data[5] = 2 * (q.j * q.k + q.i * q.r);
			data[6] = 2 * (q.i * q.k + q.j * q.r);
			data[7] = 2 * (q.j * q.k - q.i * q.r);
			data[8] = 1 - 2 * (q.i * q.i + q.j * q.j);

		}
		
		static Matrix3by3 linearInterpolate(const Matrix3by3& m1, const Matrix3by3& m2, real prop){ //interpolates a pair of matrices

			Matrix3by3 result;
			for (int i = 0; i < 9; i++)
				result.data[i] = m1.data[i] * (1 - prop) + m2.data[i] * prop;
			return result;

		}
};

#endif