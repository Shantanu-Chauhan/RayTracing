#pragma once
#include <Eigen/StdVector> // For vectors, matrices (2d,3d,4d) and quaternions in f and d precision.
#include <Eigen/Geometry>
typedef Eigen::AlignedBox<float, 3> Bbox; // The BV type provided by Eigen
using namespace Eigen;

class Obj;
class Intersection {
public:
	Intersection() {
		t = std::numeric_limits<float>::max();
		objectHit = nullptr;
	}
	float t;
	Obj* objectHit;
	Vector3f P;		//Point of intersection in world coordinates
	Vector3f N;		//Normal of the surface at intersection point (world coordinates)
	Vector2f UV;    //Texture coordinates 2d or 3d
};
class Ray {
public:
	Ray() {
		Q = Vector3f(0.0f, 0.0f, 0.0f);
		D = Vector3f(0.0f, 0.0f, 0.0f);
	}
	Ray(Vector3f q, Vector3f d);
	Vector3f Q;
	Vector3f D;
	Vector3f eval(float t);
};


class Shape {
public:
	Shape() {
		parent = nullptr;
	}
	virtual bool Intersect(Ray* ray, Intersection&) { return true; }
	Obj* parent;
	Bbox hell;
};

class Sphere :public Shape {
public:
	Sphere(Vector3f center, float radius);
	bool Intersect(Ray* ray, Intersection&);
	Vector3f center;
	float radius;
};

class Slab {
public:
	Vector3f Normal; //Both the normals for the slabs are same
	float d0, d1;	 //Distance from the origin?
};

class Box :public Shape {
public:
	Box(Vector3f corner, Vector3f diagonal);
	bool Intersect(Ray* ray, Intersection& intersection);
	Slab slabb[3];
};

class Cylinder :public Shape {
public:
	Cylinder(Vector3f _base, Vector3f _axis, float _radius);
	bool Intersect(Ray* ray, Intersection& intersection);
	Vector3f axis;
	float radius;
	Vector3f base;
};
struct MeshData;
class Triangle :public Shape {
public:
	Triangle(MeshData* meshdata, unsigned int i0, unsigned int i1, unsigned int i2);
	bool Intersect(Ray* ray, Intersection& intersection);
	Vector3f v0, v1, v2, n0, n1, n2;
	Vector2f t0, t1, t2;
};

class Interval {
public:
	float T0, T1;	 //Begining and ending point along the ray
	Vector3f N0, N1; //Surface normals at t0 and t1 respectively
	Interval() // Initialize t0 to 0 and t1 to inifinity
	{
		T0 = 0.0f;
		T1 = std::numeric_limits<float>::max();
	}
	//Constructor(t0, t1, N0, N1) : Reorders t0, t1(and N0, N1) so t0 <= t1, and stores all 4.
	Interval(float t0, float t1, Vector3f n0, Vector3f n1);
	void Empty();
	//void Intersect(void* other);
	bool Intersect(Ray* ray, Slab slab);
};