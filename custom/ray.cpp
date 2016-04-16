#include "ray.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <utilities.h>

//creates ray given a starting point and a vector
Ray::Ray(double start_x, double start_y, double start_z, double vec_x, double vec_y, double vec_z) {
  set_start_and_vector(start_x, start_y, start_z, vec_x, vec_y, vec_z);
}

Ray::Ray(double vec_x, double vec_y, double vec_z) {
  set_start_and_vector(0.0f, 0.0f, 0.0f, vec_x, vec_y, vec_z);
}

Ray::Ray(std::vector<double>& start, std::vector<double>& vec) {
  set_start_and_vector(start[0], start[1], start[2], vec[0], vec[1], vec[2]);
}

Ray::Ray(std::vector<double>& vec) {
  set_start_and_vector(0.0f, 0.0f, 0.0f, vec[0], vec[1], vec[2]);
}

Ray::~Ray() {}

//calculates a ray's intersection with a sphere, if it exists
bool Ray::intersection(Sphere& s, double intersection_point[3], double& t) {
  double a, b, c = 0.0;
  double center_start_diff;
  for(int i = 0; i < 3; ++i) {
    center_start_diff = start[i] - s.position[i];

    a += pow(vec[i], 2.0);
    b += (center_start_diff * vec[i]);
    c += pow(center_start_diff, 2.0);
  }
  b *= 2;
  c -= pow(s.radius, 2.0);

  //Ray starts inside the sphere
  if(c < 0.0) return false;
  //Ray starts on surface
  if(c == 0.0) {
    t = 0.0;
    return true;
  }

  double determinant = pow(b, 2.0) - 4.0 * a * c;

  if(determinant < 0.0) return false;
  if(determinant != 0.0) {
    determinant = sqrt(determinant);
  } 

  double twice_a = 2.0 * a;
  double neg_b = -1.0 * b;

  double t0 = (neg_b + determinant) / twice_a;
  double t1 = (neg_b - determinant) / twice_a;

  if(t0 < 0.0 && t1 < 0.0) return false;
  if(t0 > 0.0 && t1 > 0.0) {
    t = std::min(t0, t1);
  } else {
    t = std::max(t0, t1);
  }
  for(int i = 0; i < 3; ++i) {
    intersection_point[i] = start[i] + vec[i] * t;
  }
  return true;
}

//calculates a ray's intersection with a triangle, if it exists
bool Ray::intersection(Triangle& tri, double intersection_point[3], double& t) {
  double vec1[3];
  double vec2[3];

  for(int i = 0; i < 3; ++i) {
    vec1[i] = tri.v[1].position[i] - tri.v[0].position[i];
    vec2[i] = tri.v[2].position[i] - tri.v[0].position[i];
  }

  double interect_with_plane_t;
  //determine if ray intersects plane triangle lies in
  bool exists_intersect = intersect_with_plane(vec1, vec2, tri.v[0].position, interect_with_plane_t);
  if(!exists_intersect || interect_with_plane_t < 0.0) return false;
  //determine a plane to project to
  int plane = Utilities::non_ortho_plane(vec1, vec2);
  if(plane == -1) return false;
  double ray_intersect[3];
  for(int i = 0; i < 3; ++i) {
    ray_intersect[i] = start[i] + vec[i] * interect_with_plane_t;
  }

  double ray_intersect_2d[2];
  double vec1_2d[2];
  double vec2_2d[2];
  double tri_base_point[2];

  //project to 2D
  int to_set = 0;
  for(int i = 0; i < 3; ++i) {
    if(i != plane) {
      ray_intersect_2d[to_set] = ray_intersect[i];
      vec1_2d[to_set] = vec1[i];
      vec2_2d[to_set] = vec2[i];
      tri_base_point[to_set] = tri.v[0].position[i];
      ++to_set;
    }
  }
  
  //calculate u, v
  double u, v;

  if(vec2_2d[0] == 0.0) {
    if(vec2_2d[1] == 0.0) {
      return false;
    } else if (vec1_2d[0] == 0.0) {
      return false;
    } else {
      v = (ray_intersect_2d[0] - tri_base_point[0]) / vec1_2d[0];
      u = (ray_intersect_2d[1] - tri_base_point[1] - vec1_2d[1]*v) / vec2_2d[1];
    }

  } else if(vec1_2d[1] == 0.0) {
    if(vec2_2d[1] == 0.0) {
      return false;
    } else if (vec1_2d[0] == 0.0) {
      return false;
    } else {
      u = (ray_intersect_2d[1] - tri_base_point[1]) / vec2_2d[1];
      v = (ray_intersect_2d[0] - tri_base_point[0] - vec2_2d[0]) / vec1_2d[0];
    }
  } else {
    v = vec2_2d[0] * (ray_intersect_2d[1] - tri_base_point[1]);
    v -= (vec2_2d[1] * (ray_intersect_2d[0] - tri_base_point[0]));
    v /= (vec2_2d[0] * vec1_2d[1] - vec2_2d[1] * vec1_2d[0]);

    u = ray_intersect_2d[0] - tri_base_point[0] - vec1_2d[0] * v;
    u /= vec2_2d[0];
  }

  if(u < 0.0 || v < 0.0) return false;
  if(u + v > 1.0) return false; 
  t = interect_with_plane_t;
  if(t < 0.0) {
    return false;
  }

  for(int i = 0; i < 3; ++i) {
    intersection_point[i] = start[i] + vec[i] * t;
  }
  return true;
}

void Ray::set_start_and_vector(double start_x, double start_y, double start_z, double vec_x, double vec_y, double vec_z) {
  start.push_back(start_x);
  start.push_back(start_y);
  start.push_back(start_z);
  vec.push_back(vec_x);
  vec.push_back(vec_y);
  vec.push_back(vec_z);
}

bool Ray::intersect_with_plane(double vec1[], double vec2[], double point[], double& t) {
  double res[3];
  Utilities::normal(vec2, vec1, res);
  return intersect_with_plane(res, point, t);
}

//determine if ray intersect a plane
bool Ray::intersect_with_plane(double normal[], double point[], double& t) {
  double d = Utilities::dot_product(normal, point);
  double normal_dot_ray_vec = Utilities::dot_product(normal, vec);

  if(normal_dot_ray_vec == 0.0) return false;
  t = (d - Utilities::dot_product(normal, start)) / normal_dot_ray_vec;
  return true;
}
