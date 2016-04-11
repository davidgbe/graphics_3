#ifndef VERTEX_H
#define VERTEX_H

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

#endif
