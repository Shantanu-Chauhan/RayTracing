////////////////////////////////////////////////////////////////////////////////
// Temporary code.  Remove this from your raytracer.  This displays
// the contents of a scene file in realtime in a GLUT window.
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <fstream>
#include <vector>

#include <glbinding/gl/gl.h>
#include <glbinding/Binding.h>
using namespace gl;
#include <freeglut.h>


const float Radians = PI/180.0f;    // Convert degrees to radians

////////////////////////////////////////////////////////////////////////
// Shader programming class;  Encapsulates a OpenGL Shader.
////////////////////////////////////////////////////////////////////////
class ShaderProgram
{
public:
    int program;
    
    void CreateProgram() { program = glCreateProgram(); }
    void Use() { glUseProgram(program); }
    void Unuse() { glUseProgram(0); }
    void CreateShader(const std::string fname, const GLenum type);
    void LinkProgram();

};

////////////////////////////////////////////////////////////////////////
// Obj: encapsulates objects to be drawn; uses OpenGL's VAOs
////////////////////////////////////////////////////////////////////////
class Obj
{
public:
    MeshData* meshdata;
    Matrix4f modelTR;
    Material* material;
    Vector3f center;
    unsigned int vao;
    Obj(MeshData* m, const Matrix4f& tr, Material* b);
	Shape* shape = nullptr;
    void draw();
    Vector3f Center() { return center; }
};

////////////////////////////////////////////////////////////////////////
// Realtime handles all realtime drawing/interaction
////////////////////////////////////////////////////////////////////////
/// 
class Ray;
class Intersection;
class Realtime
{
public:
	unsigned int texture;
	bool writeRayTracedImage;
    bool nav;
    char motionkey;
    float speed;
	Color* ImagePointer;

	std::vector<float> SkyDome;
	int SkyDomeWidth;
	int SkyDomeHeight;
	float SkyDomeRadius;
	float* pBuffer;
	float* pUDist;
    // Camera/viewing parameters
    Vector3f ambient;
    Vector3f eye;      // Position of eye for viewing scene
    Quaternionf orient;   // Represents rotation of -Z to view direction
	float D;// Distance from the eye to the plane of focus(in world coordinates)
	float W; // Size of the circle of confusion(in fraction of screen size e.g. 0.1)
	bool DOF;
    float ry;
    float front, back;
    float spin, tilt;
    float cDist;              // Distance from eye to center of scene
    //float lightSpin, lightTilt, lightDist;

    int mouseX, mouseY;
    bool shifted;
    bool leftDown;
    bool middleDown;
    bool rightDown;

    MeshData* sphMesh;
    MeshData* boxMesh;
    MeshData* cylMesh;

    ShaderProgram lighting;
	ShaderProgram rayTracer;

	KdBVH<float, 3, Shape*> Tree;

	unsigned int ColorBlockID;
	unsigned int blockID;


	void DrawFSQ();

    int width, height;
    void setScreen(const int _width, const int _height) { width=_width;  height=_height; }
    void setCamera(const Vector3f& _eye, const Quaternionf& _o, const float _ry,float _d,float _w,int _DOF)
	{
		eye = _eye; orient = _o; ry = _ry; D = _d; W = _w; DOF = _DOF;
	}
    void setAmbient(const Vector3f& _a) { ambient = _a; }
    int setTexture(const int width, const int height, unsigned char* image);
    
    std::vector<Obj*> objs;
    std::vector<Obj*> lights;
	std::vector<Shape*> shapes;
	std::vector<Shape*> BB;

    Quaternionf ViewQuaternion() {
        Quaternionf q = angleAxis((tilt-90.0f)*Radians, Vector3f(1,0,0))
                       *orient.conjugate()
                       *angleAxis(spin*Radians, Vector3f(0,0,1));
        return q.conjugate();
    }

    Vector3f ViewDirection() {
        return ViewQuaternion().toRotationMatrix() * Vector3f(0.0f, 0.0f, -1.0f);
    }
	Vector3f Realtime::TracePath(Ray ray);
    void DrawScene();
	void RayTracerDrawScene();
    void ReshapeWindow(int w, int h);
    void KeyboardUp(unsigned char key, int x, int y);
    void KeyboardDown(unsigned char key, int x, int y);
    void MouseButton(int button, int state, int x, int y);
    void MouseMotion(int x, int y);
    
    void sphere(const Vector3f center, const float r, Material* mat);
    void box(const Vector3f base, const Vector3f diag, Material* mat);
    void cylinder(const Vector3f base, const Vector3f axis, const float radius, Material* mat);

	void triangleMesh(MeshData* meshdata);
	Intersection SampleLight();
    Realtime();
    void run(Color* image);
};


