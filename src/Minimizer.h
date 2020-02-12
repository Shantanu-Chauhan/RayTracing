//#pragma once
//#include<vector>
//#include <Eigen_unsupported/Eigen/BVH>
//#include <Eigen/StdVector> // For vectors, matrices (2d,3d,4d) and quaternions in f and d precision.
//#include <Eigen/Geometry>
//const float INF = std::numeric_limits<float>::max();
//typedef Eigen::AlignedBox<float, 3> Bbox; // The BV type provided by Eigen
//
//using namespace Eigen;
//
//class Ray;
//class Shape;
//class Intersection;
//class Minimizer {
//public:
//	Minimizer(std::vector<Shape*>::iterator start, std::vector<Shape*>::iterator end);
//	typedef float Scalar; // KdBVH needs Minimizer::Scalar defined
//	KdBVH<float, 3, Shape*> Tree;
//	// Stuff to track the minimal t and its intersection info
//	// Constructor
//	//Minimizer(const Ray& r) : ray(r) { ray = r; }
//	// Called by BVMinimize to intersect the ray with a Shape. This
//	// should return the intersection t, but should also track
//	// the minimum t and it's corresponding intersection info.
//	// Return INF to indicate no intersection.
//	float minimumOnObject(Shape* obj) {
//		// t = obj->intersect(ray); // or whatever
//		//… // Keep track the minimal intersection and object
//		//	return t;
//	}
//	// Called by BVMinimize to intersect the ray with a Bbox and
//	// returns the t value. This should be similar to the already
//	// written box (3 slab) intersection. (The difference being: a
//	// distance of zero should be returned if the ray starts within the Bbox.)
//	// Return INF to indicate no intersection.
//	float minimumOnVolume(const Bbox& box)
//	{
//		Vector3f L = box.min(); // Box corner
//		Vector3f U = box.max(); // Box corner
//	}
//};