#include"RaytracerHelper.h"
#include "raytrace.h"
Ray::Ray(Vector3f q, Vector3f d)
{
	Q = q;
	D = d;
}

Vector3f Ray::eval(float t)
{
	Vector3f result = Q + t * D;
	return result;
}

Sphere::Sphere(Vector3f _center, float _radius)
{
	center = _center;
	radius = _radius;
	boundingBox = Bbox(center - Vector3f(radius, radius, radius), center + Vector3f(radius, radius, radius));
}

bool Sphere::Intersect(Ray* ray, Intersection& intersection)
{
	Vector3f Qbar = ray->Q - center;
	float QbarDotD = Qbar.dot(ray->D);
	float QbarDotQbar = Qbar.dot(Qbar);
	float QbarDotDSquare = QbarDotD * QbarDotD;
	float RadiusSquare = radius * radius;
	float Discriminant = QbarDotDSquare - QbarDotQbar + RadiusSquare;

	if (Discriminant == 0.0f)
		return false;

	float T1 = -QbarDotD + std::sqrt(Discriminant);
	float T2 = -QbarDotD - std::sqrt(Discriminant);

	if (T1 < 0.0f && T2 < 0.0f)
		return false;
	float resultT = std::numeric_limits<float>::max();
	if (T1 < T2 && !signbit(T1))
		resultT = T1;
	else if (T2 < T1 && !signbit(T2))
		resultT = T2;

	if (resultT < intersection.t)
	{
		intersection.t = resultT;
		intersection.P = ray->eval(resultT);
		Vector3f normal(intersection.P - center);
		intersection.N = normal.normalized();
		float theta = atan2(normal.dot(Vector3f(0, 1, 0)), normal.dot(Vector3f(1, 0, 0)));
		float phi = acos(normal.dot(Vector3f(0, 0, 1)));
		intersection.UV = Vector2f(theta / 2.0f * 3.14f, phi / 3.14f);
		intersection.objectHit = this->parent;
		return true;
	}
	return false;
}

Box::Box(Vector3f corner, Vector3f diagonal)
{
	slabb[0].Normal = Vector3f(1, 0, 0);
	slabb[0].d0 = -corner.x();
	slabb[0].d1 = -corner.x() - diagonal.x();

	slabb[1].Normal = Vector3f(0, 1, 0);
	slabb[1].d0 = -corner.y();
	slabb[1].d1 = -corner.y() - diagonal.y();

	slabb[2].Normal = Vector3f(0, 0, 1);
	slabb[2].d0 = -corner.z();
	slabb[2].d1 = -corner.z() - diagonal.z();

	boundingBox = Bbox(corner, corner + diagonal);
}

bool Box::Intersect(Ray* ray, Intersection& intersection)
{
	Interval test;
	Interval solution;
	for (int i = 0; i < 3; i++)
	{
		if (test.Intersect(ray, slabb[i]))
		{
			if (test.T0 > solution.T0)
			{
				solution.N0 = test.N0;
				solution.T0 = test.T0;
			}
			if (test.T1 < solution.T1)
			{
				solution.T1 = test.T1;
				solution.N1 = test.N1;
			}
		}
	}
	if (solution.T0 > solution.T1)//no intersection
		return false;
	else
	{
		if (solution.T0 <= solution.T1 && !signbit(solution.T0))
		{
			if (solution.T0 < intersection.t)
			{
				intersection.P = ray->eval(solution.T0);
				intersection.N = solution.N0.normalized();
				intersection.t = solution.T0;
				intersection.objectHit = this->parent;
				return true;
			}
		}
		else if (solution.T1 < solution.T0 && !signbit(solution.T1))
		{
			if (solution.T1 < intersection.t)
			{
				intersection.P = ray->eval(solution.T1);
				intersection.N = solution.N1.normalized();
				intersection.t = solution.T1;
				intersection.objectHit = this->parent;
				return true;
			}
		}
		else
			return false;
	}
}

Cylinder::Cylinder(Vector3f _base, Vector3f _axis, float _radius)
{
	base = _base;
	axis = _axis;
	radius = _radius;
	Vector3f BASE = Vector3f(base + Vector3f(radius, radius, radius));
	Vector3f BASE2 = Vector3f(base - Vector3f(radius, radius, radius));
	Vector3f BASE3 = Vector3f(axis + base - Vector3f(radius, radius, radius));
	Vector3f BASE4 = Vector3f(axis + base - Vector3f(radius, radius, radius));
	float minx, maxx, miny, maxy, minz, maxz;
	minx = std::min(std::min(BASE.x(), BASE2.x()), std::min(BASE3.x(), BASE4.x()));
	miny = std::min(std::min(BASE.y(), BASE2.y()), std::min(BASE3.y(), BASE4.y()));
	minz = std::min(std::min(BASE.z(), BASE2.z()), std::min(BASE3.z(), BASE4.z()));

	maxx = std::max(std::max(BASE.x(), BASE2.x()), std::max(BASE3.x(), BASE4.x()));
	maxy = std::max(std::max(BASE.y(), BASE2.y()), std::max(BASE3.y(), BASE4.y()));
	maxz = std::max(std::max(BASE.z(), BASE2.z()), std::max(BASE3.z(), BASE4.z()));

	boundingBox = Bbox(Vector3f(minx, miny, minz), Vector3f(maxx, maxy, maxz));

}

bool Cylinder::Intersect(Ray* ray, Intersection& intersection)
{
	Ray* test = new Ray();
	Quaternionf q = Quaternionf::FromTwoVectors(axis, Vector3f::UnitZ());
	test->Q = q._transformVector(ray->Q - base);
	test->D = q._transformVector(ray->D);


	Interval a0;
	Interval b0;
	Interval c0;
	Interval solution;
	Slab basee;
	basee.Normal = Vector3f(0.0f, 0.0f, 1.0f);
	basee.d0 = 0.0f;
	basee.d1 = -(axis.norm());
	b0.Intersect(test, basee);
	float b = 2 * (test->D.x() * test->Q.x() + test->D.y() * test->Q.y());

	float bSq = b * b;

	float a = test->D.x() * test->D.x() + test->D.y() * test->D.y();

	float c = test->Q.x() * test->Q.x() + test->Q.y() * test->Q.y() - (radius * radius);

	float BSQMinus4AC = bSq - (4 * a * c);
	if (BSQMinus4AC < 0)
		return false;
	float t0, t1;


	t0 = (-b - std::sqrtf(BSQMinus4AC)) / (2 * a);
	t1 = (-b + std::sqrtf(BSQMinus4AC)) / (2 * a);
	//Calculate normal as well?
	if (t0 > t1)
	{
		float temp = t1;
		t1 = t0;
		t0 = temp;
	}
	c0.T0 = t0;
	c0.T1 = t1;

	Vector3f M0 = test->eval(t0);
	Vector3f M1 = test->eval(t1);
	c0.N0 = Vector3f(M0.x(), M0.y(), 0.0f);
	c0.N1 = Vector3f(M1.x(), M1.y(), 0.0f);

	if (a0.T0 >= b0.T0 && a0.T0 >= c0.T0)
	{
		solution.T0 = a0.T0;
		solution.N0 = a0.N0;
	}
	if (b0.T0 >= a0.T0 && b0.T0 >= c0.T0)
	{
		solution.T0 = b0.T0;
		solution.N0 = q.conjugate()._transformVector(b0.N0);
	}
	if (c0.T0 >= a0.T0 && c0.T0 >= b0.T0)
	{
		solution.T0 = c0.T0;
		solution.N0 = q.conjugate()._transformVector(c0.N0);
	}
	//-------------------T1
	if (a0.T1 <= b0.T1 && a0.T0 <= c0.T1)
	{
		solution.T1 = a0.T1;
		solution.N1 = a0.N1;
	}
	if (b0.T1 <= a0.T1 && b0.T1 <= c0.T1)
	{
		solution.T1 = b0.T1;
		solution.N1 = q.conjugate()._transformVector(b0.N1);
	}
	if (c0.T1 <= a0.T1 && c0.T1 <= b0.T1)
	{
		solution.T1 = c0.T1;
		solution.N1 = q.conjugate()._transformVector(c0.N1);
	}

	if (solution.T0 > solution.T1)
		return false;
	else
	{
		if (solution.T0 <= solution.T1 && !signbit(solution.T0))
		{
			if (solution.T0 < intersection.t)
			{
				intersection.P = ray->eval(solution.T0);
				intersection.N = solution.N0.normalized();
				intersection.t = solution.T0;
				intersection.objectHit = this->parent;
				return true;
			}
		}
		else if (solution.T1 < solution.T0 && !signbit(solution.T1))
		{
			if (solution.T1 < intersection.t)
			{
				intersection.P = ray->eval(solution.T1);
				intersection.N = solution.N1.normalized();
				intersection.t = solution.T1;
				intersection.objectHit = this->parent;
				return true;
			}
		}
	}
	return false;
}

Triangle::Triangle(MeshData* meshdata, unsigned int i0, unsigned int i1, unsigned int i2)
{
	v0 = meshdata->vertices[i0].pnt;
	n0 = meshdata->vertices[i0].nrm;
	t0 = meshdata->vertices[i0].tex;

	v1 = meshdata->vertices[i1].pnt;
	n1 = meshdata->vertices[i1].nrm;
	t1 = meshdata->vertices[i1].tex;

	v1 = meshdata->vertices[i2].pnt;
	n1 = meshdata->vertices[i2].nrm;
	t1 = meshdata->vertices[i2].tex;
}

bool Triangle::Intersect(Ray* ray, Intersection& intersection)
{
	Vector3f E1, E2;
	E1 = Vector3f(v1 - v0);
	E2 = Vector3f(v2 - v0);
	Vector3f p = ray->D.cross(E2);
	float d = p.dot(E1);
	if (d == 0.0f)
		//return No Intersection , Ray is parallel to triangle

		return false;
	Vector3f S = ray->Q - v0;
	float u = p.dot(S) / d;
	if (u < 0 || u>1)
		//return No Intersection,Ray intersects the plane but outside of E2 edge

		return false;
	Vector3f q = S.cross(E1);
	float v = ray->D.dot(q) / d;
	if (v < 0 || (u + v) >1)
		//return No Intersection,Ray intersects the plane but outside the other edges

		return false;
	float t = E2.dot(q) / d;
	if (t < 0)
		//return NoIntersection Ray's negative half intersects triangle 

		return false;
	//Intersection detected
	if (t < intersection.t)
	{
		intersection.P = ray->eval(t);
		//intersection.N = (1 - u - v) * n0[i] + u * n1[i] + v * n2[i];
		intersection.N = E2.cross(E1);
		intersection.N.normalize();
		intersection.UV = (1 - u - v) * t0 + u * t1 + v * t2;
		intersection.objectHit = this->parent;
		intersection.t = t;
		return true;
	}
	return false;
}

void Interval::Empty()
{
	T0 = 0.0f;
	T1 = -1.0f;
}

bool Interval::Intersect(Ray* ray, Slab slab)
{
	float SlabNormaldotRayD = slab.Normal.dot(ray->D);
	if (SlabNormaldotRayD != 0.0f)
	{
		//plane 0
		T0 = -(slab.d0 + slab.Normal.dot(ray->Q)) / SlabNormaldotRayD;
		N0 = -(slab.Normal);
		//plane 1
		T1 = -(slab.d1 + slab.Normal.dot(ray->Q)) / SlabNormaldotRayD;

		N0 = -slab.Normal;
		N1 = slab.Normal;

		if (T0 > T1)
		{
			float temp = T1;
			T1 = T0;
			T0 = temp;
			Vector3f tempN = N1;
			N1 = N0;
			N0 = tempN;
		}
		return true;
	}
	else
	{
		float SlabNormalDotRayQ = slab.Normal.dot(ray->Q);
		float s0, s1;
		s0 = SlabNormalDotRayQ + slab.d0;
		s1 = SlabNormalDotRayQ + slab.d1;
		if (s0 < 0.0f && s1 >= 0.0f || s1 < 0.0f && s0 >= 0.0f)
		{
			T0 = 0.0f;
			T1 = std::numeric_limits<float>::max();
			return true;
		}
		else
		{
			this->Empty();
			return true;
		}
	}
	//ray intersects both the slabs
}

Interval::Interval(float t0, float t1, Vector3f n0, Vector3f n1)
{
	if (t0 <= t1)
	{
		T0 = t0;
		T1 = t1;
		N0 = n0;
		N1 = n1;
	}
	else {
		T0 = t1;
		T1 = t0;
		N0 = n1;
		N1 = n0;
	}

}

