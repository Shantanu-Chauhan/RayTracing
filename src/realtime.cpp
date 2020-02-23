////////////////////////////////////////////////////////////////////////////////
// Temporary code.  Remove this from your raytracer.  This displays
// the contents of a scene file in realtime in a GLUT window.
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <fstream>
#include <vector>

#include "geom.h"
#include "raytrace.h"
#include "realtime.h"
#include"RaytracerHelper.h"
#include <random>
#include"Minimizer.h"

#include  <iterator>
std::mt19937_64 RNGen;
std::uniform_real_distribution<> myrandom(0.0, 1.0);

// Stupid C++ needs callbacks to be static functions.
static Realtime* globalRealtime = nullptr;
void CBDrawScene() { globalRealtime->DrawScene(); }
void CBRayTracerDrawScene() { globalRealtime->RayTracerDrawScene(); }
void CBReshapeWindow(int w, int h) { globalRealtime->ReshapeWindow(w, h); }
void CBKeyboardDown(unsigned char key, int x, int y) { globalRealtime->KeyboardDown(key, x, y); }
void CBKeyboardUp(unsigned char key, int x, int y) { globalRealtime->KeyboardUp(key, x, y); }
void CBMouseButton(int button, int state, int x, int y) { globalRealtime->MouseButton(button, state, x, y); }
void CBMouseMotion(int x, int y) { globalRealtime->MouseMotion(x, y); }
void CBAnimate(int value)
{
	glutTimerFunc(30, CBAnimate, 1);
	// atime = 360.0*glutGet(GLUT_ELAPSED_TIME)/12000;
	glutPostRedisplay();
}

unsigned int MakeVAO(MeshData* meshdata)
{
	unsigned int vao;
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);


	if (meshdata && meshdata->vertices.size() == 0) {
		std::cerr << "Missing meshdata->vertices in MakeVAO\n";
		exit(-1);
	}
	if (meshdata && meshdata->triangles.size() == 0) {
		std::cerr << "Missing meshdata->triangles in MakeVAO\n";
		exit(-1);
	}

	GLuint Pbuff;
	glGenBuffers(1, &Pbuff);
	glBindBuffer(GL_ARRAY_BUFFER, Pbuff);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 11 * meshdata->vertices.size(),
		&(meshdata->vertices[0].pnt.x()), GL_STATIC_DRAW);

	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 11 * sizeof(float), (void*)0);

	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 11 * sizeof(float), (void*)(3 * sizeof(float)));

	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 11 * sizeof(float), (void*)(6 * sizeof(float)));

	glEnableVertexAttribArray(3);
	glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, 11 * sizeof(float), (void*)(8 * sizeof(float)));

	glBindBuffer(GL_ARRAY_BUFFER, 0);


	GLuint Ibuff;
	glGenBuffers(1, &Ibuff);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, Ibuff);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * 3 * meshdata->triangles.size(),
		&(meshdata->triangles[0].x()), GL_STATIC_DRAW);

	glBindVertexArray(0);

	return vao;
}

MeshData* SphMesh()
{
	MeshData* meshdata = new MeshData();
	unsigned int n = 20;
	float d = 2.0f * PI / float(n * 2);
	for (unsigned int i = 0; i <= n * 2; i++) {
		float s = i * 2.0f * PI / float(n * 2);
		for (unsigned int j = 0; j <= n; j++) {
			float t = j * PI / float(n);
			float x = cos(s) * sin(t);
			float y = sin(s) * sin(t);
			float z = cos(t);
			meshdata->vertices.push_back(VertexData(Vector3f(x, y, z),
				Vector3f(x, y, z),
				Vector2f(s / (2 * PI), t / PI),
				Vector3f(sin(s), cos(s), 0.0)));
			if (i > 0 && j > 0) {
				meshdata->triangles.push_back(TriData((i - 1) * (n + 1) + (j - 1),
					(i - 1) * (n + 1) + (j),
					(i) * (n + 1) + (j)));
				meshdata->triangles.push_back(TriData((i - 1) * (n + 1) + (j - 1),
					(i) * (n + 1) + (j),
					(i) * (n + 1) + (j - 1)));
			}
		}
	}
	return meshdata;
}

MeshData* BoxMesh()
{
	Matrix4f face[6] = {
		Matrix4f::Identity(),
		rotate(180.0f * Radians, Vector3f(1.0f, 0.0f, 0.0f)),
		rotate(90.0f * Radians, Vector3f(1.0f, 0.0f, 0.0f)),
		rotate(-90.0f * Radians, Vector3f(1.0f, 0.0f, 0.0f)),
		rotate(90.0f * Radians, Vector3f(0.0f, 1.0f, 0.0f)),
		rotate(-90.0f * Radians, Vector3f(0.0f, 1.0f, 0.0f)) };

	Matrix4f half = translate(Vector3f(0.5f, 0.5f, 0.5f)) * scale(Vector3f(0.5f, 0.5f, 0.5f));
	MeshData* meshdata = new MeshData();
	for (unsigned int f = 0; f < 6; f++) {
		Matrix4f m4 = half * face[f];
		Matrix3f m3 = m4.block<3, 3>(0, 0);
		for (unsigned int i = 0; i < 2; i++) {
			for (unsigned int j = 0; j < 2; j++) {
				Vector4f p = m4 * Vector4f(float(2 * i) - 1.0f, float(2 * j) - 1.0f, 1.0f, 1.0f);
				Vector3f tnrm = m3 * Vector3f(0.0f, 0.0f, 1.0f);
				Vector3f ttan = m3 * Vector3f(1.0, 0.0, 0.0);
				meshdata->vertices.push_back(VertexData(Vector3f(p[0], p[1], p[2]),
					Vector3f(tnrm[0], tnrm[1], tnrm[2]),
					Vector2f(float(i), float(j)),
					Vector3f(ttan[0], ttan[1], ttan[2])));
				meshdata->triangles.push_back(TriData(4 * f + 0, 4 * f + 1, 4 * f + 3));
				meshdata->triangles.push_back(TriData(4 * f + 0, 4 * f + 3, 4 * f + 2));
			}
		}
	}
	return meshdata;
}

MeshData* CylMesh()
{
	MeshData* meshdata = new MeshData();
	unsigned int n = 20;
	float d = 2.0f * PI / float(n * 2);
	for (unsigned int i = 0; i <= n; i++) {
		float s = i * 2.0f * PI / float(n);
		float x = cos(s);
		float y = sin(s);

		meshdata->vertices.push_back(VertexData(Vector3f(x, y, 0.0f),
			Vector3f(x, y, 0.0f),
			Vector2f(s / (2 * PI), 0.0f),
			Vector3f(-sin(s), cos(s), 0.0f)));

		meshdata->vertices.push_back(VertexData(Vector3f(x, y, 1.0f),
			Vector3f(x, y, 0.0f),
			Vector2f(s / (2 * PI), 0.0f),
			Vector3f(-sin(s), cos(s), 0.0f)));

		if (i > 0) {
			meshdata->triangles.push_back(TriData((i - 1) * 2 + 1, (i - 1) * 2, (i) * 2));
			meshdata->triangles.push_back(TriData((i - 1) * 2 + 1, (i) * 2, (i) * 2 + 1));
		}
	}
	return meshdata;
}

////////////////////////////////////////////////////////////////////////
// Shader programming class;  Encapsulates a OpenGL Shader.
////////////////////////////////////////////////////////////////////////
void ShaderProgram::CreateShader(const std::string fname, const GLenum type)
{
	// Read a file into a string
	std::ifstream f;
	f.open(fname, std::ios_base::binary); // Open
	f.seekg(0, std::ios_base::end);       // Position at end
	int length = f.tellg();               // to get the length

	char* src = new char[length + 1];  // Create buffer of needed length
	f.seekg(0, std::ios_base::beg);      // Position at beginning
	f.read(src, length);             //   to read complete file
	f.close();                            // Close

	src[length] = char(0);            // Finish with a NULL

	// Create a shader, attach, send it the source, and compile it.
	int shader = glCreateShader(type);
	glAttachShader(program, shader);
	glShaderSource(shader, 1, &src, nullptr);
	glCompileShader(shader);

	// Get the compilation status
	int status;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &status);

	// If compilation status is not OK, get and print the log message.
	if (status != 1) {
		int length;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &length);
		char* buffer = new char[length];
		glGetShaderInfoLog(shader, length, nullptr, buffer);
		printf("Compile log(%s):\n%s\n", type == GL_VERTEX_SHADER ? "Vertex" : "Fragment", buffer);
		delete buffer;
	}
}

void ShaderProgram::LinkProgram()
{
	// Link program and check the status
	glLinkProgram(program);
	int status;
	glGetProgramiv(program, GL_LINK_STATUS, &status);

	// If link failed, get and print log
	if (status != 1) {
		int length;
		glGetProgramiv(program, GL_INFO_LOG_LENGTH, &length);
		char* buffer = new char[length];
		glGetProgramInfoLog(program, length, nullptr, buffer);
		printf("Link log:\n%s\n", buffer);
		delete buffer;
	}
}


void applyMaterial(Material* mat, const unsigned int program)
{
	int loc = glGetUniformLocation(program, "Kd");
	glUniform3fv(loc, 1, &mat->Kd[0]);

	loc = glGetUniformLocation(program, "Ks");
	glUniform3fv(loc, 1, &mat->Ks[0]);

	loc = glGetUniformLocation(program, "alpha");
	glUniform1f(loc, mat->alpha);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, mat->texid);
	loc = glGetUniformLocation(program, "tex");
	glUniform1i(loc, 0);

	loc = glGetUniformLocation(program, "emitter");
	glUniform1i(loc, 0);
}

////////////////////////////////////////////////////////////////////////
// Light: encapsulates a light and communiction with a shader.
////////////////////////////////////////////////////////////////////////
void applyLight(Material* mat, const unsigned int program)
{
	Vector3f Z;

	int loc = glGetUniformLocation(program, "Kd");
	glUniform3fv(loc, 1, &mat->Kd[0]);

	loc = glGetUniformLocation(program, "emitter");
	glUniform1i(loc, 1);
}

////////////////////////////////////////////////////////////////////////
// Obj: encapsulates objects to be drawn; uses OpenGL's VAOs
////////////////////////////////////////////////////////////////////////
Obj::Obj(MeshData* m, const Matrix4f& tr, Material* b)
	: meshdata(m), modelTR(tr), material(b)
{
	Vector4f sum(0, 0, 0, 0);
	//for (int i=0;  i<meshdata->vertices.size();  i++)
	//    sum += modelTR*Vector4f(v.pnt[0], v.pnt[1], v.pnt[2], 1.0);
	for (auto v : meshdata->vertices) // C++11 for loop
		sum += modelTR * Vector4f(v.pnt[0], v.pnt[1], v.pnt[2], 1.0);

	center = (sum / meshdata->vertices.size()).block<3, 1>(0, 0);

	std::cout << "center: " << center.transpose() << std::endl;

	vao = MakeVAO(meshdata);
}

void Obj::draw()
{
	glBindVertexArray(vao);
	glDrawElements(GL_TRIANGLES, 3 * meshdata->triangles.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
}


////////////////////////////////////////////////////////////////////////
// Realtime handles all realtime drawing/interaction
////////////////////////////////////////////////////////////////////////

// Constructor for Realtime.  Initializes OpenGL, GLUT,as well as the
// data elements of the class.

Realtime::Realtime()
{
	// Initialize the OpenGL bindings
	glbinding::Binding::initialize(false);

	globalRealtime = this;
	// Initialize GLUT
	int argc = 0;
	char* argv;

	glutInit(&argc, &argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitContextVersion(3, 3);
	glutInitContextProfile(GLUT_COMPATIBILITY_PROFILE);

	glutInitWindowSize(200, 200);
	glutCreateWindow("Class Framework");
	glutSetOption((GLenum)GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);

	glutIgnoreKeyRepeat(true);
	glutDisplayFunc(&CBDrawScene);
	glutReshapeFunc(&CBReshapeWindow);
	glutKeyboardFunc(&CBKeyboardDown);
	glutKeyboardUpFunc(&CBKeyboardUp);
	glutMouseFunc(&CBMouseButton);
	glutMotionFunc(&CBMouseMotion);
	glutTimerFunc(30, CBAnimate, 1);

	printf("OpenGL Version: %s\n", glGetString(GL_VERSION));
	printf("GLSL Version: %s\n", glGetString(GL_SHADING_LANGUAGE_VERSION));
	printf("Rendered by: %s\n", glGetString(GL_RENDERER));
	fflush(stdout);

	// Create the shader program
	lighting.CreateProgram();
	lighting.CreateShader("realtime.vert", GL_VERTEX_SHADER);
	lighting.CreateShader("realtime.frag", GL_FRAGMENT_SHADER);

	glBindAttribLocation(lighting.program, 0, "vertex");
	glBindAttribLocation(lighting.program, 1, "vertexNormal");
	glBindAttribLocation(lighting.program, 2, "vertexTexture");
	glBindAttribLocation(lighting.program, 3, "vertexTangent");
	lighting.LinkProgram();

	// Create the shader program
	rayTracer.CreateProgram();
	rayTracer.CreateShader("raytracer.vert", GL_VERTEX_SHADER);
	rayTracer.CreateShader("raytracer.frag", GL_FRAGMENT_SHADER);

	glBindAttribLocation(rayTracer.program, 0, "vertex");

	rayTracer.LinkProgram();

	glGenBuffers(1, &blockID); // Generates block
	unsigned int bindPoint = 1;

	glBindBuffer(GL_UNIFORM_BUFFER, blockID);
	glBufferData(GL_UNIFORM_BUFFER, sizeof(float) * width * height * 3, NULL, GL_STATIC_DRAW);

	// Several generic meshes which can be transfofrmed to *any* sphere, box, or cylinder.
	sphMesh = SphMesh();
	boxMesh = BoxMesh();
	cylMesh = CylMesh();

	// Initialize various member attributes
	//materials.push_back(new Material());

	nav = false;
	spin = 0.0f;
	tilt = 90.0f;
	speed = 0.05f;
	front = 0.1f;
	back = 1000.0f;

	shifted = false;
	leftDown = false;
	middleDown = false;
	rightDown = false;
	motionkey = 0;

}

// This function enters the event loop.
void Realtime::run(Color* image)
{
	ImagePointer = image;
	Tree.init(shapes.begin(), shapes.end());
	cDist = eye.norm();
	glutReshapeWindow(width, height);
	glutMainLoop();
}
void Realtime::DrawFSQ()
{
	static unsigned int quadVAO = 0;
	static unsigned int quadVBO;
	if (quadVAO == 0)
	{
		float quadVert[] = {
			-1.0f, 1.0f, 0.0f, 0.0f, 1.0f,
			-1.0f, -1.0f, 0.0f, 0.0f, 0.0f,
			1.0f, 1.0f, 0.0f, 1.0f, 1.0f,
			1.0f, -1.0f, 0.0f, 1.0f, 0.0f
		};

		// Setup plane vao
		glGenVertexArrays(1, &quadVAO);
		glGenBuffers(1, &quadVBO);
		glBindVertexArray(quadVAO);
		glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
		glBufferData(GL_ARRAY_BUFFER, sizeof(quadVert), &quadVert, GL_STATIC_DRAW);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
	}
	glBindVertexArray(quadVAO);
	glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
	glBindVertexArray(0);
}
// Called when the scene needs to be redrawn.
void Realtime::DrawScene()
{
	Vector3f viewDir = ViewDirection();
	Vector2f dir2 = Vector2f(viewDir.x(), viewDir.y()).normalized();
	if (motionkey == 'w')
		eye += speed * Vector3f(dir2.x(), dir2.y(), 0.0);
	if (motionkey == 's')
		eye -= speed * Vector3f(dir2.x(), dir2.y(), 0.0);
	if (motionkey == 'd')
		eye += speed * Vector3f(dir2.y(), -dir2.x(), 0.0);
	if (motionkey == 'a')
		eye -= speed * Vector3f(dir2.y(), -dir2.x(), 0.0);
	if (motionkey == 'e')
		eye -= speed * Vector3f(0.0f, 0.0f, -1.0f);
	if (motionkey == 'c')
		eye -= speed * Vector3f(0.0f, 0.0f, 1.0f);

	int loc;

	glClearColor(0.3, 0.3, 0.3, 1.0);
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	Matrix4f WorldView;
	Matrix4f R = toMat4(ViewQuaternion().conjugate());
	WorldView = R * translate(-eye);

	float rx = (ry * width) / height;
	Matrix4f WorldProj = frustum(-front * rx, front * rx, -front * ry, front * ry, front, back);

	lighting.Use();

	loc = glGetUniformLocation(lighting.program, "WorldProj");
	glUniformMatrix4fv(loc, 1, GL_FALSE, WorldProj.data());

	loc = glGetUniformLocation(lighting.program, "WorldView");
	glUniformMatrix4fv(loc, 1, GL_FALSE, WorldView.data());

	loc = glGetUniformLocation(lighting.program, "ambient");
	glUniform3fv(loc, 1, &ambient[0]);

	loc = glGetUniformLocation(lighting.program, "eyePos");
	glUniform3fv(loc, 1, &eye[0]);

	Vector3f lightEmit[8];
	Vector3f lightPosn[8];
	int  lightNum = lights.size();
	for (int i = 0; i < lightNum; i++) {
		lightPosn[i] = lights[i]->Center();
		lightEmit[i] = lights[i]->material->Kd;
	}

	loc = glGetUniformLocation(lighting.program, "lightNum");
	glUniform1i(loc, lightNum);

	loc = glGetUniformLocation(lighting.program, "lightPosn");
	glUniform3fv(loc, 8, &lightPosn[0][0]);

	loc = glGetUniformLocation(lighting.program, "lightEmit");
	glUniform3fv(loc, 8, &lightEmit[0][0]);

	// for (unsigned int i=0;  i<objs.size();  i++) {
	//     Material* material = objs[i]->material;
	//     Matrix4f& modelTR = objs[i]->modelTR;
	//     ...
	//     objs[i]->draw();
	for (auto obj : objs) {     // C++11 loop
		Material* material = obj->material;
		Matrix4f& modelTR = obj->modelTR;
		Matrix3f normalTR = modelTR.inverse().block<3, 3>(0, 0);

		loc = glGetUniformLocation(lighting.program, "ModelTr");
		glUniformMatrix4fv(loc, 1, GL_FALSE, modelTR.data());

		loc = glGetUniformLocation(lighting.program, "normalTR");
		glUniformMatrix3fv(loc, 1, GL_FALSE, normalTR.data());

		if (!material) {
			std::cerr << "No material associated with object\n";
			exit(-1);
		}
		if (material->isLight())
			applyLight(material, lighting.program);
		else
			applyMaterial(material, lighting.program);
		obj->draw();
	}

	lighting.Unuse();
	glutSwapBuffers();
}

// Write the image as a HDR(RGBE) image.  
#include "rgbe.h"
void WriteHdrImage(const std::string outName, const int width, const int height, Color* image)
{
	// Turn image from a 2D-bottom-up array of Vector3D to an top-down-array of floats
	float* data = new float[width * height * 3];
	float* dp = data;
	for (int y = height - 1; y >= 0; --y) {
		for (int x = 0; x < width; ++x) {
			Color pixel = image[y * width + x];
			*dp++ = pixel[0];
			*dp++ = pixel[1];
			*dp++ = pixel[2];
		}
	}

	// Write image to file in HDR (a.k.a RADIANCE) format
	rgbe_header_info info;
	char errbuf[100] = { 0 };

	FILE* fp = fopen(outName.c_str(), "wb");
	info.valid = false;
	int r = RGBE_WriteHeader(fp, width, height, &info, errbuf);
	if (r != RGBE_RETURN_SUCCESS)
		printf("error: %s\n", errbuf);

	r = RGBE_WritePixels_RLE(fp, data, width, height, errbuf);
	if (r != RGBE_RETURN_SUCCESS)
		printf("error: %s\n", errbuf);
	fclose(fp);

	delete data;
}


namespace Eigen
{
	Bbox bounding_box(Shape* obj)
	{
		return obj->BoundingBox;
	}
}

const float RussianRoulette = 0.8f;


Vector3f EvalRadiance(Obj* O)
{
	return O->material->Kd;
}
Vector3f SampleLobe(Vector3f N, float c, float phi)
{
	float s = sqrt(1.0f - (c * c));
	Vector3f K = Vector3f(s * cos(phi), s * sin(phi), c);// Vector centered around Z-axis
	Quaternionf q = Quaternionf::FromTwoVectors(Vector3f::UnitZ(), N);// q rotates Z to N
	return q._transformVector(K);// K rotated to N's frame
}
Vector3f EvalScattering(Vector3f N, Vector3f wi, Vector3f Kd)
{
	float NdotWiAbs = std::max(N.dot(wi), 0.0f);
	return NdotWiAbs * Kd / PI;
}

float PdfBRDF(Vector3f N, Vector3f wi)
{
	float NdotWiAbs = std::max(N.dot(wi), 0.0f);
	return NdotWiAbs / PI;
}

Vector3f SampleBRDF(Vector3f N)
{
	float sai1, sai2;
	sai1 = myrandom(RNGen);
	sai2 = myrandom(RNGen);
	return SampleLobe(N, sqrt(sai1), 2 * PI * sai2);
}


template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
	std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
	std::advance(start, dis(g));
	return start;
}

template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
	static std::random_device rd;
	static std::mt19937 gen(rd());
	return select_randomly(start, end, gen);
}

Intersection Realtime::SampleLight()
{
	Intersection test;
	Obj* L = *select_randomly(lights.begin(), lights.end());
	//SampleSphere
	Vector3f center = L->Center();
	float radius = static_cast<Sphere*>(L->shape)->radius;
	float sai1, sai2;
	sai1 = myrandom(RNGen);
	sai2 = myrandom(RNGen);
	float z = 2 * sai1 - 1.0f;
	float r = sqrt(1 - z * z);
	float a = 2 * PI * sai2;
	test.N = Vector3f(r * cos(a), r * sin(a), z);
	test.P = center + radius * test.N;
	test.objectHit = L;
	return test;
}

float GeometryFactor(Intersection A, Intersection B)
{
	Vector3f D = A.P - B.P;
	float AnDotD = A.N.dot(D);
	float BnDotD = B.N.dot(D);
	float DdotDSq = pow(D.dot(D),2);
	return std::max((AnDotD * BnDotD) / DdotDSq, 0.0f);
}

float PdfLight(Obj* O,int NumberOfLights)
{
	return 1 / (O->shape->Area() * NumberOfLights);
}
Vector3f Realtime::TracePath(Ray* ray)
{
	Vector3f C(0.0f, 0.0f, 0.0f);	//Accumulated light
	Vector3f W(1.0f, 1.0f, 1.0f);	//Accumulated weight
	Intersection P;// = new Intersection();
	Intersection Q;// = new Intersection();
	Minimizer miniP(ray, &P);
	float minDist = BVMinimize(Tree, miniP);
	Vector3f N = P.N;
	if (P.objectHit == nullptr)
		return C;
	if (P.objectHit->material->isLight())
	{
		return EvalRadiance(P.objectHit);
	}
	while (myrandom(RNGen) <= 0.8f)
	{
		Vector3f wi;
		Vector3f f;
		float p;
		//Explicit Light Correction
		bool check = true;
		{
			Intersection L = SampleLight();
			p = PdfLight(L.objectHit, lights.size()) / GeometryFactor(P, L);
			wi = L.P - P.P;
			Intersection I;
			Ray rayP;
			rayP.D = wi;
			rayP.Q = P.P;
			Minimizer miniI(&rayP, &I);
			BVMinimize(Tree, miniI);
			if (p > 0 && I.objectHit != nullptr && I.P == L.P)
			{
				f = EvalScattering(N, wi, L.objectHit->material->Kd);
				C += W.cwiseProduct(EvalRadiance(L.objectHit));
			}
		}
		wi = SampleBRDF(N);
		//wi.normalize();
		Ray NewRay;
		NewRay.D = wi;
		NewRay.Q = P.P;
		Minimizer miniQ(&NewRay, &Q);
		BVMinimize(Tree, miniQ);
		if (Q.objectHit == nullptr)
			break;
		f = EvalScattering(N, wi, P.objectHit->material->Kd);
		p = PdfBRDF(N, wi) * RussianRoulette;
		if (p < 0.000001f)
			break;
		W = W.cwiseProduct(f) / p;
		if (Q.objectHit->material->isLight())
		{
			C += W.cwiseProduct(EvalRadiance(Q.objectHit));
			break;
		}
		P.N = Q.N;
		P.objectHit = Q.objectHit;
		P.P = Q.P;
		P.t = Q.t;
		P.UV = Q.UV;
	}
	return C;
}

void Realtime::RayTracerDrawScene()
{
	//Tree.init(shapes.begin(), shapes.end());
	//Minimizer mini(shapes.begin(), shapes.end());
	float rx = (ry * width) / height;
	Vector3f X = rx * ViewQuaternion()._transformVector(Vector3f::UnitX());
	Vector3f Y = ry * ViewQuaternion()._transformVector(Vector3f::UnitY());
	Vector3f Z = -1 * ViewQuaternion()._transformVector(Vector3f::UnitZ());
	//Minimizer mini(shapes.begin(),shapes.end());
	int loc;
	for (int i = 1; i < 128; i++)
	{
#pragma omp parallel for schedule(dynamic, 1) // Magic: Multi-thread y loop
		for (int y = 0; y < height; y++) {

			fprintf(stderr, "Pass : %4d Rendering %4d\r",i, y);
			for (int x = 0; x < width; x++) {
				float dx, dy;
				//dx = 2.0f * (x + 0.5f) / width - 1.0f;
				//dy = 2.0f * (y + 0.5f) / height - 1.0f;
				dx = 2.0f * (x + myrandom(RNGen)) / width - 1.0f;
				dy = 2.0f * (y + myrandom(RNGen)) / height - 1.0f;
				Vector3f direction(X * dx + dy * Y + Z);
				direction.normalize();
				Ray ray;
				ray.Q = eye;
				ray.D = direction;
				//Intersection* frontMost = new Intersection();
				//Minimizer mini(ray, frontMost);
				//float minDist = BVMinimize(Tree, mini);
				/*for (int i = 0; i < shapes.size(); i++)
				{
					if(!shapes[i]->parent->material->isLight())
					shapes[i]->Intersect(ray, *frontMost);
				}*/
				Vector3f color;
				color = TracePath(&ray);

				//if (frontMost->objectHit == nullptr)
				//	color = Color(0.0, 0.0, 0.0);
				//else
				//{
				//	Vector3f Ia = Vector3f(0.2f, 0.2f, 0.2f);
				//	Vector3f Ii = Vector3f(1.1f, 1.1f, 1.1f);
				//	Vector3f L = lights[0]->center - frontMost->P;
				//	L.normalize();
				//	Vector3f N = frontMost->N;
				//	N.normalize();
				//	Vector3f V = ViewDirection();
				//	V.normalize();
				//	Vector3f H = L + V;
				//	H.normalize();
				//	float NL = std::max(N.dot(L), 0.0f);
				//	float NH = std::max(N.dot(H), 0.0f);
				//	NH = pow(NH, frontMost->objectHit->material->alpha);
				//	//color = frontMost->N;
				//	Vector3f C = frontMost->objectHit->material->Kd.cwiseProduct(Ia) + frontMost->objectHit->material->Kd.cwiseProduct(Ii * NL) + frontMost->objectHit->material->Ks.cwiseProduct(Ii * NH);
				//	color = C;
				//	//color = frontMost->objectHit->material->Kd;
				//	//color = Vector3f((frontMost->t-5.0f)/4.0f, (frontMost->t - 5.0f) / 4.0f, (frontMost->t - 5.0f) / 4.0f);
				//	//color = frontMost->P;
				//}
				Vector3f test = ImagePointer[y * width + x];
				test += color;
				ImagePointer[y * width + x] = test;
			}
		}
		//if(i==2|| i == 8 || i == 64 || i == 512 || i == 256 || i == 128 ){
		{
			glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			rayTracer.Use();
			glGenTextures(1, &texture);
			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, texture);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, (GLint)GL_REPEAT);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, (GLint)GL_REPEAT);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, (GLint)GL_LINEAR_MIPMAP_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, (GLint)GL_LINEAR);

			glTexImage2D(GL_TEXTURE_2D, 0, (GLint)GL_RGB32F, width, height, 0, GL_RGB, GL_FLOAT, &ImagePointer[0]);
			glGenerateMipmap(GL_TEXTURE_2D);

			
			loc = glGetUniformLocation(rayTracer.program, "pass");
			glUniform1i(loc, i);
			loc = glGetUniformLocation(rayTracer.program, "Image");
			glUniform1i(loc, 1);
			DrawFSQ();
			glutSwapBuffers();
		}
	}
	fprintf(stderr, "Image written\n");
	std::string inName = "testscene.scn";
	std::string hdrName = inName;

	hdrName.replace(hdrName.size() - 3, hdrName.size(), "hdr");
	WriteHdrImage(hdrName, width, height, ImagePointer);
}

// Called by GLUT when the window size is changed.
void Realtime::ReshapeWindow(int w, int h)
{
	if (w && h)
		glViewport(0, 0, w, h);
	width = w;
	height = h;
	// Force a redraw
	glutPostRedisplay();
}

// Called by GLUT for keyboard actions.
void Realtime::KeyboardDown(unsigned char key, int x, int y)
{
	printf("key: %c\n", key);
	switch (key) {
	case 9:
		nav = !nav;
		break;

	case 'v': {
		Quaternionf q = ViewQuaternion();
		printf("camera  %g %g %g   %g   q %g %g %g %g\n",
			eye[0], eye[1], eye[2], ry, q.w(), q.x(), q.y(), q.z());
		printf("screen %d %d\n", width, height);
		fflush(stdout); }
			  break;

	case 'w': case 's': case 'a': case 'd': case 'e': case 'c':
		motionkey = key;
		break;

	case 'r':
		glutDisplayFunc(&CBRayTracerDrawScene);
		break;
	case 27:                    // Escape key
	case 'q':
		glutLeaveMainLoop();
		break;
	}
}

void Realtime::KeyboardUp(unsigned char key, int x, int y)
{
	motionkey = 0;
	fflush(stdout);
}

// Called by GLut when a mouse button changes state.
void Realtime::MouseButton(int button, int state, int x, int y)
{
	// Record the position of the mouse click.
	mouseX = x;
	mouseY = y;

	// Test if the SHIFT keey was down for this mouse click
	shifted = glutGetModifiers() && GLUT_ACTIVE_SHIFT;

	// Ignore high order bits, set by some (stupid) GLUT implementation.
	button = button % 8;

	if (button == GLUT_LEFT_BUTTON) {
		leftDown = (state == GLUT_DOWN);
	}

	else if (button == GLUT_MIDDLE_BUTTON) {
		middleDown = (state == GLUT_DOWN);
	}

	else if (button == GLUT_RIGHT_BUTTON) {
		rightDown = (state == GLUT_DOWN);
	}

	else if (button == 3) {
		Vector3f C = eye + cDist * ViewDirection();
		cDist = pow(cDist, 1.0f / 1.02f);
		eye = C - cDist * ViewDirection();
	}

	else if (button == 4) {
		Vector3f C = eye + cDist * ViewDirection();
		cDist = pow(cDist, 1.02f);
		eye = C - cDist * ViewDirection();
	}

	// Force a redraw
	glutPostRedisplay();
	fflush(stdout);
}

void Realtime::MouseMotion(int x, int y)
{
	// Calculate the change in the mouse position
	int dx = x - mouseX;
	int dy = y - mouseY;

	if (leftDown) {        // Rotate light position
		if (nav) {
			spin += dx / 2.0f;
			tilt += dy / 2.0f;
		}
		else {
			Vector3f C = eye + cDist * ViewDirection();
			spin += dx / 2.0f;
			tilt += dy / 2.0f;
			eye = C - cDist * ViewDirection();
		}
	}

	else if (middleDown) {}

	// Record this position
	mouseX = x;
	mouseY = y;

	// Force a redraw
	glutPostRedisplay();
}


void Realtime::sphere(const Vector3f center, const float r, Material* mat)
{
	Matrix4f m = translate(center) * scale(Vector3f(r, r, r));
	Vector3f rrr(r, r, r);
	Obj* obj = new Obj(sphMesh, m, mat);
	obj->shape = new Sphere(center, r);
	obj->shape->parent = obj;
	objs.push_back(obj);
	shapes.push_back(obj->shape);
	if (mat->isLight())
		lights.push_back(obj);
}

void Realtime::box(const Vector3f base, const Vector3f diag, Material* mat)
{
	Matrix4f m = translate(base) * scale(Vector3f(diag[0], diag[1], diag[2]));
	Obj* obj = new Obj(boxMesh, m, mat);
	obj->shape = new Box(base, diag);
	obj->shape->parent = obj;
	objs.push_back(obj);
	shapes.push_back(obj->shape);
	if (mat->isLight())
		lights.push_back(obj);
}


void Realtime::cylinder(const Vector3f base, const Vector3f axis, const float radius, Material* mat)
{
	Vector3f Z(0.0f, 0.0f, 1.0f);
	Vector3f C = axis.normalized();
	Vector3f B = C.cross(Z);
	if (B.norm() < 1e-8)
		B = Vector3f(0, 1, 0);
	else
		B = B.normalized();
	Vector3f A = B.cross(C).normalized();
	Matrix4f R;
	R << A(0), B(0), C(0), 0.0f,    // Row-wise text (although init is col-wise)
		A(1), B(1), C(1), 0.0f,
		A(2), B(2), C(2), 0.0f,
		0.0f, 0.0f, 0.0f, 1.0;

	Matrix4f m = translate(base) * R * scale(Vector3f(radius, radius, axis.norm()));
	Vector3f rrr(radius, radius, radius);
	Obj* obj = new Obj(cylMesh, m, mat);
	obj->shape = new Cylinder(base, axis, radius);
	obj->shape->parent = obj;
	objs.push_back(obj);
	shapes.push_back(obj->shape);
	if (mat->isLight())
		lights.push_back(obj);
}

void Realtime::triangleMesh(MeshData* meshdata)
{
	Obj* obj = new Obj(meshdata, Matrix4f::Identity(), meshdata->mat);
	for (int i = 0; i < meshdata->triangles.size(); i++)
	{
		Shape* mpShape = new Triangle(meshdata, meshdata->triangles[i].x(), meshdata->triangles[i].y(), meshdata->triangles[i].z());
		mpShape->parent = obj;
		shapes.push_back(mpShape);
	}
	objs.push_back(obj);
	if (meshdata->mat->isLight())
		lights.push_back(obj);
}

