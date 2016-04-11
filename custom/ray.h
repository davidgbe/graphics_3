#ifndef RAY_H
#define RAY_H

#include <vector>
#include <light.h>
#include <sphere.h>
#include <vertex.h>
#include <triangle.h>

class Ray {
  public:
    Ray(double start_x, double start_y, double start_z, double vec_x, double vec_y, double vec_z);
    Ray(double vec_x, double vec_y, double vec_z);
    Ray(std::vector<double>& start, std::vector<double>& vec);
    Ray(std::vector<double>& vec);
    ~Ray();
    bool intersection(Sphere& s, double& t);
    bool intersection(Triangle& tri, double& t);

  private:
    void set_start_and_vector(double start_x, double start_y, double start_z, double vec_x, double vec_y, double vec_z);
    bool intersect_with_plane(double vec1[], double vec2[], double point[], double& t);
    bool intersect_with_plane(double normal[], double point[], double& t);

    std::vector<double> start;
    std::vector<double> vec;
};

#endif
