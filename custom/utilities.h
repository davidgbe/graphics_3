#ifndef UTILITIES_H
#define UTILITIES_H

#include <math.h>
#include <iostream>
#include <vector>

class Utilities {
  public:
    static double calc_triangle_area(double vec1[], double vec2[]) {
      double res[3];
      Utilities::cross_product(vec1, vec2, res);
      return Utilities::calc_triangle_area(res);
    }

    static double calc_triangle_area(double unnormalized_normal[]) {
      double sum = 0.0;
      for(int i = 0; i < 3; ++i) {
        sum += pow(unnormalized_normal[i], 2.0);
      }
      return sum;
    }

    static void normalize(double vec[])
    {
      double mag = 0;
      for(int i = 0; i < 3; ++i) {
        mag += pow(vec[i], 2.0f);
      }
      if(mag == 0.0) {
        return;
      }
      mag = sqrt(mag);
      for(int i = 0; i < 3; ++i) {
        vec[i] = vec[i] / mag;
      }
    }

    static void normal(double vec1[], double vec2[], double res[]) {
      Utilities::cross_product(vec1, vec2, res);
      Utilities::normalize(res);
    }

    static void cross_product(double vec1[], double vec2[], double res[]) {
      res[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
      res[1] = -1.0f * (vec1[0] * vec2[2] - vec1[2] * vec2[0]);
      res[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
    }

    static double dot_product(std::vector<double>& vec1, double vec2[]) {
      Utilities::dot_product(vec2, vec1);
    }

    static double dot_product(double vec1[], std::vector<double>& vec2) {
      double sub_vec[3];
      for(int i = 0; i < 3; ++i) {
        sub_vec[i] = vec2[i];
      }
      return Utilities::dot_product(vec1, sub_vec);
    }

    static double dot_product(std::vector<double>& vec1, std::vector<double>& vec2) {
      double sub_vec[3];
      double sub_vec2[3];
      for(int i = 0; i < 3; ++i) {
        sub_vec[i] = vec1[i];
        sub_vec2[i] = vec2[i];
      }
      return Utilities::dot_product(sub_vec, sub_vec2);
    }

    static double dot_product(double vec1[], double vec2[]) {
      double sum = 0.0;
      for(int i = 0; i < 3; ++i) {
        sum += (vec1[i] * vec2[i]);
      }
      return sum;
    }

    static double magnitude(double vec[]) {
      return sqrt(Utilities::dot_product(vec, vec));
    }

    static double triangle_area(double vec1[], double vec2[]) {
      double total_area_vec[3];
      Utilities::cross_product(vec1, vec2, total_area_vec);
      return 0.5 * Utilities::magnitude(total_area_vec);
    }

    static int non_ortho_plane(double tri_vec_1[], double tri_vec_2[]) {
      double trial_vec[3];
      trial_vec[0] = 1.0;
      trial_vec[1] = 0.0;
      trial_vec[2] = 0.0;
      if(Utilities::dot_product(trial_vec, tri_vec_1) == 0.0) return 0;
      trial_vec[0] = 0.0;
      trial_vec[1] = 1.0;
      if(Utilities::dot_product(trial_vec, tri_vec_1) == 0.0) return 1;
      trial_vec[1] = 0.0;
      trial_vec[2] = 1.0;
      if(Utilities::dot_product(trial_vec, tri_vec_1) == 0.0) return 2;
      return -1;
    }

    static void mult_by_scalar(double vec[], double c) {
      for(int i = 0; i < 3; ++i) {
        vec[i] = c * vec[i];
      }
    }

    static void duplicate(double to_dup[], double dup[]) {
      for(int i = 0; i < 3; ++i) {
        to_dup[i] = dup[i];
      }
    }
};

#endif
