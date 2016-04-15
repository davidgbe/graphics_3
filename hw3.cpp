/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: <Your name here>
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <iostream>
#include <math.h>

#include <utilities.h>

#include <light.h>
#include <sphere.h>
#include <vertex.h>
#include <triangle.h>
#include <imageIO.h>

#include <ray.h>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

// #define WIDTH 40
// #define HEIGHT 30

double aspect_ratio = double (WIDTH) / HEIGHT;

//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

double a = 1.0;
double b = 1.0;
double c = 1.0;

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

bool find_closest_intersection(Ray& r, double& min_t, double min_coorinates[3], std::string& type, Triangle* min_tri, Sphere* min_sphere) {
  min_t = 1000000.0;
  type = "tri";
  double t;
  double temp_intersection[3];
  bool exists_intersect = false;

  for(int idx = 0; idx < num_triangles; ++idx) {
    bool intersect = r.intersection(triangles[idx], temp_intersection, t);
    if(intersect) {
      exists_intersect = true;
      if(t < min_t) {
        min_t = t;
        Utilities::duplicate(temp_intersection, min_coorinates);
        min_tri = &(triangles[idx]);
      }
    }
  }

  for(int idx = 0; idx < num_spheres; ++idx) {
    bool intersect = r.intersection(spheres[idx], temp_intersection, t);
    if(intersect) {
      exists_intersect = true;
      if(t < min_t) {
        std::cout << t << std::endl;
        type = "sphere";
        min_t = t;
        Utilities::duplicate(temp_intersection, min_coorinates);
        min_sphere = &(spheres[idx]);
      }
    }
  }

  return exists_intersect;
}

bool can_see_light(double point[3], Light& l)
{
  double t;
  double min_coorinates[3];
  std::string type;
  Triangle* min_tri;
  Sphere* min_sphere;


  Ray shadow_ray(point[0], point[1], point[2], l.position[0] - point[0], l.position[1] - point[1], l.position[2] - point[2]);

  bool exists_intersect = find_closest_intersection(shadow_ray, t, min_coorinates, type, min_tri, min_sphere);
  if(exists_intersect && t < 1.0) return false;
  return true;
}

void local_phong(double color_diffuse[3], double color_specular[3], double shininess, double point[3], double surface_normal[3], Light& light, double divide_mag, double colors[])
{
  double normal_to_light[3];
  for(int i = 0; i < 3; ++i) {
    normal_to_light[i] = light.position[i] - point[i];
  }
  Utilities::normalize(normal_to_light);
  double reflection_normal[3];
  double constant = 2.0 * Utilities::dot_product(normal_to_light, surface_normal);
  for(int i = 0; i < 3; ++i) {
    reflection_normal[i] = constant * surface_normal[i] - normal_to_light[i];
  }
  double viewer_normal[3];
  Utilities::duplicate(point, viewer_normal);
  Utilities::normalize(viewer_normal);

  double light_dot_normal = Utilities::dot_product(normal_to_light, surface_normal);
  double reflection_dot_view_norm = Utilities::dot_product(reflection_normal, viewer_normal);

  for(int i = 0; i < 3; ++i) {
    colors[i] = color_diffuse[i] * light_dot_normal;
    colors[i] += pow(color_specular[i] * reflection_dot_view_norm, shininess);
    colors[i] /= divide_mag;
    colors[i] += ambient_light[i];
  }
}

void local_phong_sphere(Sphere& s, double point[3], double colors[])
{
  double q = Utilities::magnitude(point);
  double divide_mag = 1.0 / ( a * pow(a, 2.0) + b + q + c);

  for(int i = 0; i < 3; ++i) {
    colors[i] = 0.0;
  }

  double surface_normal[3];
  for(int i = 0; i < 3; ++i) {
    surface_normal[i] = point[i] - s.position[i];
  }

  Utilities::normalize(surface_normal);

  for(int l = 0; l < num_lights; ++l) {
    if(!can_see_light(point, lights[l])) {
      continue;
    }

    double colors_from_light[3];

    local_phong(s.color_diffuse, s.color_specular, s.shininess, point, surface_normal, lights[l], divide_mag, colors_from_light);
    for(int i = 0; i < 3; ++i) {
      colors[i] += colors_from_light[i];
    }
  }
}

void get_barycentric(Triangle& t, double point[3], double barycentric[3])
{
  std::cout << "HERE" << std::endl;
  std::cout << t.v[0].position[2] << std::endl;
  std::cout << "HERE" << std::endl;
  for(int i = 0; i < 3; ++i) {
    std::cout << point[i] << std::endl;
  }
  double vec01[3];
  double vec02[3];
  double vec12[3];
  double vec0point[3];
  double vec1point[3];
  for(int i = 0; i < 3; ++i) {
    vec01[i] = t.v[1].position[i] - t.v[0].position[i];
    vec02[i] = t.v[2].position[i] - t.v[0].position[i];
    vec12[i] = t.v[2].position[i] - t.v[1].position[i];
    vec0point[3] = point[i] = t.v[0].position[i];
    vec1point[3] = point[i] = t.v[1].position[i];
  }
  double total_area = Utilities::triangle_area(vec01, vec02);
  std::cout << "AREA: " << total_area << std::endl;
  barycentric[0] = Utilities::triangle_area(vec01, vec0point) / total_area;
  barycentric[1] = Utilities::triangle_area(vec12, vec1point) / total_area;
  barycentric[2] = 1.0 - barycentric[1] - barycentric[0];
}

void interpolated_normal(double barycentric[3], Vertex v[3], double res[])
{
  for(int i = 0; i < 3; ++i) {
    res[i] = barycentric[0] * v[0].normal[i] + barycentric[1] * v[1].normal[i] + barycentric[2] * v[2].normal[i];
  }
}

void interpolated_color_diffuse(double barycentric[3], Vertex v[3], double res[])
{
  for(int i = 0; i < 3; ++i) {
    res[i] = barycentric[0] * v[0].color_diffuse[i] + barycentric[1] * v[1].color_diffuse[i] + barycentric[2] * v[2].color_diffuse[i];
  }
}

void interpolated_color_specular(double barycentric[3], Vertex v[3], double res[])
{
  for(int i = 0; i < 3; ++i) {
    res[i] = barycentric[0] * v[0].color_specular[i] + barycentric[1] * v[1].color_specular[i] + barycentric[2] * v[2].color_specular[i];
  }
}

double interpolated_shininess(double barycentric[3], Vertex v[3])
{
  return barycentric[0] * v[0].shininess + barycentric[1] * v[1].shininess + barycentric[2] * v[2].shininess;
}

void local_phong_triangle(Triangle& t, double point[3], double colors[])
{
  double q = Utilities::magnitude(point);
  double divide_mag = 1.0 / ( a * pow(a, 2.0) + b + q + c);

  for(int i = 0; i < 3; ++i) {
    colors[i] = 0.0;
  }

  double barycentric[3];
  get_barycentric(t, point, barycentric);

  double surface_normal[3];
  interpolated_normal(barycentric, t.v, surface_normal);
  Utilities::normalize(surface_normal);

  double surface_color_diffuse[3];
  interpolated_color_specular(barycentric, t.v, surface_color_diffuse);

  double surface_color_specular[3];
  interpolated_color_specular(barycentric, t.v, surface_color_specular);

  double surface_shininess = interpolated_shininess(barycentric, t.v);

  for(int l = 0; l < num_lights; ++l) {
    if(!can_see_light(point, lights[l])) {
      continue;
    }
    double colors_from_light[3];
    local_phong(surface_color_diffuse, surface_color_specular, surface_shininess, point, surface_normal, lights[l], divide_mag, colors_from_light);
    for(int i = 0; i < 3; ++i) {
      colors[i] += colors_from_light[i];
    }
  }
}

void generate_color_for_ray(Ray& r, double colors[3])
{
  double t;
  double min_coorinates[3];
  std::string type;
  Triangle* min_tri;
  Sphere* min_sphere;

  bool exists_intersect = find_closest_intersection(r, t, min_coorinates, type, min_tri, min_sphere);
  // std::cout << "INTERSECT" << std::endl;
  // std::cout << exists_intersect << std::endl;
  // std::cout << type << std::endl;
  // for(int i = 0; i < 3; ++i) {
  //   std::cout << min_coorinates[i] << std::endl;
  // }
  if(!exists_intersect) {
    for(int i = 0; i < 3; ++i) {
      colors[i] = ambient_light[i];
    }
  } else if(type == "tri") {
    std::cout << type << std::endl;
    local_phong_triangle(*min_tri, min_coorinates, colors);
  } else if(type == "sphere") {
    std::cout << type << std::endl;
    local_phong_sphere(*min_sphere, min_coorinates, colors);
  } 
}

Ray** cast_rays()
{
  double tan_fov_over_2 = tan(fov/2.0);
  Ray*** rays = new Ray**[WIDTH];

  double x_coord_base = -1.0 * aspect_ratio * tan_fov_over_2;

  for(unsigned int x = 0; x < WIDTH; ++x)
  {
    rays[x] = new Ray*[HEIGHT];
    double x_coord = x_coord_base + (float(x) / WIDTH) * -2.0 * x_coord_base;
    double y_coord_base = -1.0 * tan_fov_over_2;
    for(unsigned int y = 0; y < HEIGHT; ++y)
    {
      double y_coord = y_coord_base + (float(y) / HEIGHT) * -2.0 * y_coord_base;
      rays[x][y] = new Ray(x_coord, y_coord, -1.0);
      double colors[3];
      generate_color_for_ray(*(rays[x][y]), colors);

    }
  }
}

//MODIFY THIS FUNCTION
void draw_scene()
{
  cast_rays();
  //a simple test output
  for(unsigned int x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(unsigned int y=0; y<HEIGHT; y++)
    {
      plot_pixel(x, y, x % 256, y % 256, (x+y) % 256);
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else 
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

