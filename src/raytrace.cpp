//////////////////////////////////////////////////////////////////////
// Provides the framework for a raytracer.
////////////////////////////////////////////////////////////////////////

#include <vector>

#ifdef _WIN32
    // Includes for Windows
    #include <windows.h>
    #include <cstdlib>
    #include <limits>
    #include <crtdbg.h>
#else
    // Includes for Linux
#endif

#include "geom.h"
#include "raytrace.h"
#include "realtime.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// A good quality *thread-safe* Mersenne Twister random number generator.
//#include <random>
//std::mt19937_64 RNGen;
//std::uniform_real_distribution<> myrandom(0.0, 1.0);
// Call myrandom(RNGen) to get a uniformly distributed random number in [0,1].
#include"rgbe.h"
void Scene::ReadHdrImage(const std::string inName, std::vector<float>& image,
	int& width, int& height)
{
	rgbe_header_info info;
	char errbuf[100] = { 0 };

	// Open file and read width and height from the header
	FILE* fp = fopen(inName.c_str(), "rb");
	if (!fp) {
		printf("Can't open file: %s\n", inName.c_str());
		return;
	}
	int rc = RGBE_ReadHeader(fp, &width, &height, &info, errbuf);
	if (rc != RGBE_RETURN_SUCCESS) {
		printf("RGBE read error: %s\n", errbuf);
		return;
	}

	// Allocate enough memory
	image.resize(3 * width * height);

	// Read the pixel data and close the file
	rc = RGBE_ReadPixels_RLE(fp, &image[0], width, height, errbuf);
	if (rc != RGBE_RETURN_SUCCESS) {
		printf("RGBE read error: %s\n", errbuf);
		return;
	}
	fclose(fp);

	printf("Read %s (%dX%d)\n", inName.c_str(), width, height);
}

Scene::Scene()
{ 
    realtime = new Realtime(); 
}

void Scene::Finit()
{
}

void Scene::triangleMesh(MeshData* mesh) 
{ 
    realtime->triangleMesh(mesh); 
}

Quaternionf Orientation(int i, 
                        const std::vector<std::string>& strings,
                        const std::vector<float>& f)
{
    Quaternionf q(1,0,0,0); // Unit quaternion
    while (i<strings.size()) {
        std::string c = strings[i++];
        if (c == "x")  
            q *= angleAxis(f[i++]*Radians, Vector3f::UnitX());
        else if (c == "y")  
            q *= angleAxis(f[i++]*Radians, Vector3f::UnitY());
        else if (c == "z")  
            q *= angleAxis(f[i++]*Radians, Vector3f::UnitZ());
        else if (c == "q")  {
            q *= Quaternionf(f[i+0], f[i+1], f[i+2], f[i+3]);
            i+=4; }
        else if (c == "a")  {
            q *= angleAxis(f[i+0]*Radians, Vector3f(f[i+1], f[i+2], f[i+3]).normalized());
            i+=4; } }
    return q;
}

////////////////////////////////////////////////////////////////////////
// Material: encapsulates surface properties
void Material::setTexture(const std::string path)
{
    int width, height, n;
    stbi_set_flip_vertically_on_load(true);
    unsigned char* image = stbi_load(path.c_str(), &width, &height, &n, 0);

    // Realtime code below:  This sends the texture in *image to the graphics card.
    // The raytracer will not use this code (nor any features of OpenGL nor the graphics card).
    glGenTextures(1, &texid);
    glBindTexture(GL_TEXTURE_2D, texid);
    glTexImage2D(GL_TEXTURE_2D, 0, n, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, image);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 100);
    glGenerateMipmap(GL_TEXTURE_2D);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, (int)GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, (int)GL_LINEAR_MIPMAP_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

    stbi_image_free(image);
}

void PreProcess(float** pBuffer,float** pUDist,int width,int height,std::vector<float>& image)
{
	// Pre-processing step: Marginal and conditional CDF
	*pBuffer = new float[width * (height + 1)];
	*pUDist = &(*pBuffer)[width * height];
	float* pSinTheta = new float[height];
	float angleFrac = PI / float(height);
	float theta = angleFrac * 0.5f;
	for (unsigned int i = 0; i < height; i++, theta += angleFrac)
		pSinTheta[i] = sin(theta);
	for (unsigned int i = 0, m = 0; i < width; i++, m += height) {
		float* pVDist = &(*pBuffer)[m];
		unsigned int k = i * 3;
		pVDist[0] = 0.2126f * image[k + 0] + 0.7152f * image[k + 1] + 0.0722f * image[k + 2];
		pVDist[0] *= pSinTheta[0];
		for (unsigned int j = 1, k = (width + i) * 3; j < height; j++, k += width * 3) {
			float lum = 0.2126 * image[k + 0] + 0.7152 * image[k + 1] + 0.0722 * image[k + 2];
			pVDist[j] = pVDist[j - 1] + lum * pSinTheta[j];
		}
		if (i == 0)
			(*pUDist)[i] = pVDist[height - 1];
		else
			(*pUDist)[i] = (*pUDist)[i - 1] + pVDist[height - 1];
	}
}
void Scene::Command(const std::vector<std::string>& strings,
                     std::vector<float>& f)
{
    if (strings.size() == 0) return;
    std::string c = strings[0];
    
    if (c == "screen") {
        // syntax: screen width height
        realtime->setScreen(int(f[1]),int(f[2]));
        width = int(f[1]);
        height = int(f[2]); }

	else if (c == "skydome") {
		
		// syntax: screen width height
		std::string file = "skyDome.hdr";
		ReadHdrImage(file, realtime->SkyDome, realtime->SkyDomeWidth, realtime->SkyDomeHeight);
		PreProcess(&(realtime->pBuffer), &(realtime->pUDist), realtime->SkyDomeWidth, realtime->SkyDomeHeight, realtime->SkyDome);
		realtime->SkyDomeRadius = f[4];
		realtime->sphere(Vector3f(f[1], f[2], f[3]), f[4], currentMat);
	}

    else if (c == "camera") {
        // syntax: camera x y z   ry   <orientation spec>
        // Eye position (x,y,z),  view orientation (qw qx qy qz),  frustum height ratio ry
        realtime->setCamera(Vector3f(f[1],f[2],f[3]), Orientation(5,strings,f), f[4],f[10],f[11],f[12]); 
	/*	f[5] = 0.0;
		f[6] = 0.0;*/
	}


    else if (c == "ambient") {
        // syntax: ambient r g b
        // Sets the ambient color.  Note: This parameter is temporary.
        // It will be ignored once your raytracer becomes capable of
        // accurately *calculating* the true ambient light.
        realtime->setAmbient(Vector3f(f[1], f[2], f[3])); }

    else if (c == "brdf")  {
        // syntax: brdf  r g b   r g b  alpha
        // later:  brdf  r g b   r g b  alpha  r g b ior
        // First rgb is Diffuse reflection, second is specular reflection.
        // third is beer's law transmission followed by index of refraction.
        // Creates a Material instance to be picked up by successive shapes
        currentMat = new Material(Vector3f(f[1], f[2], f[3]), Vector3f(f[4], f[5], f[6]), f[7],Vector3f(f[8], f[9], f[10]),f[11]); }

    else if (c == "light") {
        // syntax: light  r g b   
        // The rgb is the emission of the light
        // Creates a Material instance to be picked up by successive shapes
        currentMat = new Light(Vector3f(f[1], f[2], f[3])); }
   
    else if (c == "sphere") {
        // syntax: sphere x y z   r
        // Creates a Shape instance for a sphere defined by a center and radius
        realtime->sphere(Vector3f(f[1], f[2], f[3]), f[4], currentMat); }

    else if (c == "box") {
        // syntax: box bx by bz   dx dy dz
        // Creates a Shape instance for a box defined by a corner point and diagonal vector
        realtime->box(Vector3f(f[1], f[2], f[3]), Vector3f(f[4], f[5], f[6]), currentMat); }

    else if (c == "cylinder") {
        // syntax: cylinder bx by bz   ax ay az  r
        // Creates a Shape instance for a cylinder defined by a base point, axis vector, and radius
        realtime->cylinder(Vector3f(f[1], f[2], f[3]), Vector3f(f[4], f[5], f[6]), f[7], currentMat); }


    else if (c == "mesh") {
        // syntax: mesh   filename   tx ty tz   s   <orientation>
        // Creates many Shape instances (one per triangle) by reading
        // model(s) from filename. All triangles are rotated by a
        // quaternion (qw qx qy qz), uniformly scaled by s, and
        // translated by (tx ty tz) .
        Matrix4f modelTr = translate(Vector3f(f[2],f[3],f[4]))
                          *scale(Vector3f(f[5],f[5],f[5]))
                          *toMat4(Orientation(6,strings,f));
        ReadAssimpFile(strings[1], modelTr);  }

    
    else {
        fprintf(stderr, "\n*********************************************\n");
        fprintf(stderr, "* Unknown command: %s\n", c.c_str());
        fprintf(stderr, "*********************************************\n\n");
    }
}

void Scene::TraceImage(Color* image, const int pass)
{
    realtime->run(image);                          // Remove this (realtime stuff)

//#pragma omp parallel for schedule(dynamic, 1) // Magic: Multi-thread y loop
//    for (int y=0;  y<height;  y++) {
//
//        fprintf(stderr, "Rendering %4d\r", y);
//        for (int x=0;  x<width;  x++) {
//            Color color;
//            if ((x-width/2)*(x-width/2)+(y-height/2)*(y-height/2) < 100*100)
//                color = Color(myrandom(RNGen), myrandom(RNGen), myrandom(RNGen));
//            else if (abs(x-width/2)<4 || abs(y-height/2)<4)
//                color = Color(0.0, 0.0, 0.0);
//            else 
//                color = Color(1.0, 1.0, 1.0);
//            image[y*width + x] = color;
//        }
//    }

    fprintf(stderr, "\n");
}
