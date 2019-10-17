#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>



// VEC3 CLASS
template<typename T>
class Vec3
{
public:
   // INSTANCE FIELDS
   T x, y, z;
   // CONSTRUCTORS
   Vec3() : x(T(0)), y(T(0)), z(T(0)) {}              // default constructor
   Vec3(T xx) : x(xx), y(xx), z(xx) {}                // constructor w/ 1 arg
   Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}    // constructor w/ 3 args
   // OTHER FUNCTIONS
   Vec3& normalize()
   {
      // divides components of vector by magnitude
      T nor2 = length2();
      if (nor2 > 0) {
         T invNor = 1 / sqrt(nor2);
         x *= invNor, y *= invNor, z *= invNor;
      }
      return *this;
   }
   T dot(const Vec3<T> &v) const{ return x * v.x + y * v.y + z * v.z; }
   T length2() const { return x * x + y * y + z * z; }
   T length() const { return sqrt(length2()); }
   // OPERATOR FUNCTIONS
   Vec3<T> operator * (const T &f) const { return Vec3<T>(x * f, y * f, z * f); }
   Vec3<T> operator * (const Vec3<T> &v) const { return Vec3<T>(x * v.x, y * v.y, z * v.z); }
   Vec3<T> operator - (const Vec3<T> &v) const { return Vec3<T>(x - v.x, y - v.y, z - v.z); }
   Vec3<T> operator + (const Vec3<T> &v) const { return Vec3<T>(x + v.x, y + v.y, z + v.z); }
   Vec3<T>& operator += (const Vec3<T> &v) { x += v.x, y += v.y, z += v.z; return *this; }
   Vec3<T>& operator *= (const Vec3<T> &v) { x *= v.x, y *= v.y, z *= v.z; return *this; }
   Vec3<T> operator - () const { return Vec3<T>(-x, -y, -z); }
};
typedef Vec3<float> Vec3f;





// SPHERE CLASS
// used for physical spheres + light sources
class Sphere
{
public:
   // INSTANCE FIELDS
   Vec3f center;                           // position of the sphere
   float radius, radius2;                  // sphere radius and radius^2
   Vec3f surfaceColor, emissionColor;      // surface color and emission (light)
   float transparency, reflection;         // surface transparency and reflectivity
   // CONSTRUCTOR
   Sphere(
        const Vec3f &c,
        const float &r,
        const Vec3f &sc,
        const float &refl = 0,
        const float &transp = 0,
        const Vec3f &ec = 0) :
        center(c), radius(r), radius2(r * r), surfaceColor(sc), emissionColor(ec),
        transparency(transp), reflection(refl)
   { /* empty */ }
   // INTERSECTION FUNCTION
   // compute ray-sphere intersection (scratchapixel geometric solution)
   bool intersect(const Vec3f &rayorig, const Vec3f &raydir, float &t0, float &t1) const
    {
      Vec3f l = center - rayorig;         // compute length L by subtracting eye from center
      float tca = l.dot(raydir);          // compute tca by projecting L on D vector (dot product property)
      if (tca < 0) return false;             // return if t0 is behind eye (no intersection)
      
      float d2 = l.dot(l) - tca * tca;    // calculate d^2 using pythagorean thm
      if (d2 > radius2) return false;        // return if d greater than radius (no intersection)
      float thc = sqrt(radius2 - d2);     // calculate thc using pythagorean thm

      t0 = tca - thc;                     // compute 1st intersection point (distance along ray)
      t1 = tca + thc;                     // compute 2nd intersection point (distance along ray)
        
      return true;
   }
};





// TRACE FUNCTION
// 1. takes ray as arg, checks if it intersects w/ any geometry
// 2. if any intersection
//       a. compute intersection point, normal @ intersection, and shade point
//       b. return color for point
//    else return background color
#define MAX_RAY_DEPTH 5  // max recursion depth

// helper function for trace to mix color values
float mix(const float &a, const float &b, const float &mix)
{
   return b * mix + a * (1 - mix);
}
// main trace function
Vec3f trace(
   const Vec3f &rayorig,
   const Vec3f &raydir,
   const std::vector<Sphere> &spheres,
   const int &depth)
{
   float tnear = INFINITY;       // distance to nearest sphere, updates if sphere found
   const Sphere* sphere = NULL;  // pointer to nearest sphere, updates if sphere found
   // iterate through sphere list to find intersection
   for (unsigned i = 0; i < spheres.size(); ++i) {
      float t0 = INFINITY, t1 = INFINITY;    // "intersection points" (t in parametric formula)
      if (spheres[i].intersect(rayorig, raydir, t0, t1)) {
         if (t0 < 0)       // if t0 behind eye, replace it with t1
            t0 = t1;
         if (t0 < tnear) {          // if new point closer than previous tnear
            tnear = t0;                // update min distance
            sphere = &spheres[i];      // update nearest sphere to current sphere
         }
      } 
   }
   // if there's no intersection return black or background color
   //if (!sphere) return Vec3f(2);
   if (!sphere) return Vec3f(0.50, 0.23, 0.65);    // purple

   // else there is a sphere
   // initiate variables for drawing spheres
   Vec3f surfaceColor = 0;                   // color of the ray/surfaceof the object intersected by the ray
   Vec3f phit = rayorig + raydir * tnear;    // point of intersection
   Vec3f nhit = phit - sphere->center;       // normal at the intersection point
   nhit.normalize();                         // normalize normal direction
   float bias = 1e-4;                        // add some bias to the point from which we will be tracing
   bool inside = false;                      // boolean for determining if inside sphere

   // if normal and view direction are not opposite, we are inside the sphere
   if (raydir.dot(nhit) > 0) {
      nhit = -nhit;              // reverse normal direction
      inside = true;             // set inside to true
   }

   // REFLECTIVE/TRANSPARENT OBJECTS
   if ((sphere->transparency > 0 || sphere->reflection > 0) && depth < MAX_RAY_DEPTH)
   {
      // fresnel test
      float facingratio = -raydir.dot(nhit);
      float fresnelEffect = mix(pow(1 - facingratio, 3), 1, 0.8);

      // reflective
      Vec3f reflectVector = raydir - nhit * 2 * raydir.dot(nhit); // compute reflection direction (primary ref. over normal)
      reflectVector.normalize();                                  // normalize reflection vector
      Vec3f reflectionColor = trace(phit + nhit * bias, reflectVector, spheres, depth + 1); // recursively trace from hit point
      
      // refraction
      Vec3f refractionColor = 0;
      if (sphere->transparency) {
         // COPY
         float ior = 1.1, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface?
         float cosi = -nhit.dot(raydir);
         float k = 1 - eta * eta * (1 - cosi * cosi);

         Vec3f refractVector = raydir * eta + nhit * (eta *  cosi - sqrt(k));
         refractVector.normalize();
         refractionColor = trace(phit - nhit * bias, refractVector, spheres, depth + 1);
         // END COPY
      }

      // add apply values to surface color
      //surfaceColor = (reflectionColor * sphere->surfaceColor);
      surfaceColor = reflectionColor * fresnelEffect + refractionColor * (1 - fresnelEffect) * sphere->transparency * sphere->surfaceColor;
   }

   // MAT OBJECTS (PHONG SHADING)
   else {
      // iterate through list of spheres
      for (unsigned i = 0; i < spheres.size(); ++i) {
         // if current sphere is a light source
         if (spheres[i].emissionColor.x > 0) {
            Vec3f transmission = 1;                            // "power" of light on current pixel (affects occluded)
            Vec3f lightDirection = spheres[i].center - phit;   // get direction of light to point hit
            lightDirection.normalize();                        // normalize this "ray"

            // check all other spheres
            for (unsigned j = 0; j < spheres.size(); ++j) {
               // skip if same sphere
               if (i != j) {
                  float t0, t1;
                  // cast a ray from point to light source
                  // if intersect returns true (we hit another sphere on the way), this point is in shadow
                  // decrease transmission
                  if (spheres[j].intersect(phit + nhit * bias, lightDirection, t0, t1)) {
                     transmission = 0.3;     // 0 for full dark shadows; floats for partial
                     break;
                  }
               }
            }
            // diffuse contribution
            float diffuse = std::max(float(0), nhit.dot(lightDirection));

            // specular contribution
            Vec3f reflectVector = nhit * 2 * nhit.dot(lightDirection) - lightDirection;               // compute reflect vector
            float specularExponent = 25;                                                              // reflectfactor (tbd)
            float specular = powf(std::max(float(0), -reflectVector.dot(raydir)), specularExponent);  // specular formula

            // apply values to surface color
            surfaceColor += sphere->surfaceColor * transmission * diffuse * spheres[i].emissionColor; // diffuse
            surfaceColor += specular;                                                                 // specular
         }
      }
   }
   return surfaceColor + sphere->emissionColor;
}





// RENDER FUNCTION
// 1. compute a ray for each pixel of the image
// 2. trace it, return color
//   a. if ray hits a sphere, return color @ point
//   b. else return background color
void render(const std::vector<Sphere> &spheres)
{
   // pre-trace computations
   unsigned width = 640, height = 480;    // specify image size

   Vec3f *image = new Vec3f[width * height]; // array of Vec3's (color information per pixel)
   Vec3f *pixel = image;                     // initiate pointer for current pixel (start of array)

   float invWidth = 1 / float(width);        // inverse width for trace computations
   float invHeight = 1 / float(height);      // inverse height for trace computations

   float fov = 30;                                 // field of view (tbd)
   float aspectRatio = width / float(height);      // compute aspect ratio
   float angle = tan(M_PI * 0.5 * fov / 180.0);    // angle (tbd)
   

   // iterate through every pixel
   for (unsigned y = 0; y < height; ++y) {
      for (unsigned x = 0; x < width; ++x, ++pixel)
      {
         // compute adjusted x/y for image-plane conversion
         float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectRatio;
         float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
         // ray points to current pixel on x/y, -1 (into the screen)
         Vec3f raydir(xx, yy, -1);
         raydir.normalize();
         // trace ray, save returned color info into where *pixel points
         *pixel = trace(Vec3f(0), raydir, spheres, 0);
      }
   }

   // save results to PPM file
   std::ofstream ofs("./untitled.ppm", std::ios::out | std::ios::binary);
   ofs << "P6\n" << width << " " << height << "\n255\n";                   // PPM header
   // iterate through image array: for each pixel
   for (unsigned i = 0; i < width * height; ++i) {
      ofs << (unsigned char)(std::min(float(1), image[i].x) * 255) <<      // append R value 
             (unsigned char)(std::min(float(1), image[i].y) * 255) <<      // append G value
             (unsigned char)(std::min(float(1), image[i].z) * 255);        // append B value
    }
    ofs.close();
    delete [] image;
}





// MAIN FUNCTION
// 1. create scene (5 spheres + 1 light)
// 2. render scene
int main(int argc, char **argv)
{
   // initiate list of spheres
   std::vector<Sphere> spheres;

   // build scene
   // position, radius, surface color, reflectivity, transparency, emmission color
   
   // floor
   spheres.push_back(Sphere(Vec3f( 0.0, -10004, -20), 10000, Vec3f(0.33, 0.58, 0.0), 0, 0.0));

   // pumpkins
   spheres.push_back(Sphere(Vec3f(  0,  0.0, -20),       3, Vec3f(1.0, 0.46, 0.0), 0, 0.0));
   spheres.push_back(Sphere(Vec3f( -5.0, -0.4, -26),     2.6, Vec3f(1.0, 0.46, 0.0), 0, 0.0));

   // moon
   spheres.push_back(Sphere(Vec3f(  4.5,  5, -30),       3, Vec3f(0.6, 0.6, 0.6), 1, 0.2));

   // ghosts
   spheres.push_back(Sphere(Vec3f(  2.2,  -1, -10),       1.5, Vec3f(0.0, 0.9, 1), 0, 1));
   spheres.push_back(Sphere(Vec3f(  -2,  2.6, -10),       1.5, Vec3f(0.0, 0.9, 1), 0, 1));

   // lights
   spheres.push_back(Sphere(Vec3f( 5.0,     15, 0),     3, Vec3f(0.00, 0.00, 0.00), 0, 0.0, Vec3f(1.5)));
   spheres.push_back(Sphere(Vec3f(-25.0,     20, 0),     3, Vec3f(0.00, 0.00, 0.00), 0, 0.0, Vec3f(0.4)));

   // render
   render(spheres);
   return 0;
}