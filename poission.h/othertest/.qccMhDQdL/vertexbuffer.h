#ifndef BASILISK_HEADER_17
#define BASILISK_HEADER_17
#line 1 "/home/dahuanghhc/basilisk/src/vertexbuffer.h"
/**
# Vertex buffers

They are used to store the vertex coordinates, normals and colors
computed by the the OpenGL commands typically used by [draw.h]().

These vertex buffers are the minimal information required to render
the objects. 

In combination with the ["dumb GL" implementation](gl/fb_dumb.c) this
allows to generate geometries using OpenGL commands but without any
OpenGL library. */

struct {
  // public
  Array * position, * normal, * color, * index;
  float modelview[16];
  int type;  // type = 0 -> lines, type = 1 -> mesh
  int dim;
  int vertex, nvertex;
  bool visible; // visible = true > only traverse visible vertices
  
  // private
  int line_loop, lines, line_strip ;
  int quads, polygon, fan;
  int state;
} VertexBuffer = {
  .visible = false, // traverse all vertices by default
  .modelview = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 }
};

static void vertex_buffer_push_index (unsigned int i)
{
  i -= VertexBuffer.vertex;
  array_append (VertexBuffer.index, &i, sizeof(unsigned int));
}

void vertex_buffer_setup()
{
  VertexBuffer.nvertex = 0;
  VertexBuffer.type = -1;
  VertexBuffer.dim = -1;
  VertexBuffer.position = array_new();
  VertexBuffer.normal = array_new();
  VertexBuffer.color = array_new();
  VertexBuffer.index = array_new();
}

void vertex_buffer_free()
{
  array_free (VertexBuffer.position);
  VertexBuffer.position = NULL;
  array_free (VertexBuffer.normal);
  VertexBuffer.normal = NULL;
  array_free (VertexBuffer.color);
  VertexBuffer.color = NULL;
  array_free (VertexBuffer.index);
  VertexBuffer.index = NULL;
}

static void vertex_buffer_glBegin (int state)
{
  if (VertexBuffer.index) {

    glGetFloatv (GL_MODELVIEW_MATRIX, VertexBuffer.modelview);

    bview * view = get_view();

    float q[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,
		    - view->tx, - view->ty, 3, 1 };
    matrix_multiply (q, VertexBuffer.modelview);
    for (int i = 0; i < 16; i++)
      VertexBuffer.modelview[i] = q[i];

    gl_build_rotmatrix ((float (*)[4])q, view->quat);
    swap (float, q[1], q[4]);
    swap (float, q[2], q[8]);
    swap (float, q[6], q[9]);
    matrix_multiply (q, VertexBuffer.modelview);
    for (int i = 0; i < 16; i++)
      VertexBuffer.modelview[i] = q[i];

    VertexBuffer.state = state;
    switch (state) {
    case GL_LINE_LOOP:
      VertexBuffer.line_loop = VertexBuffer.nvertex;
      break;
    case GL_LINES:
      VertexBuffer.lines = VertexBuffer.nvertex;
      break;
    case GL_LINE_STRIP:
      VertexBuffer.line_strip = VertexBuffer.nvertex;
      break;    
    case GL_QUADS:
      VertexBuffer.quads = VertexBuffer.nvertex;
      break;
    case GL_POLYGON:
      VertexBuffer.polygon = VertexBuffer.nvertex;
      break;
    case GL_TRIANGLE_FAN:
      VertexBuffer.fan = VertexBuffer.nvertex;
      break;
    default:
      fprintf (stderr, "glBegin (%d) not implemented yet\n", state);
      break;
    }
  }
  glBegin (state);
}

static void vertex_buffer_glEnd()
{
  glEnd();
  if (VertexBuffer.index) {
    int type = -1;
    switch (VertexBuffer.state) {

    case GL_LINE_LOOP:
      for (int i = VertexBuffer.line_loop; i < VertexBuffer.nvertex - 1; i++) {
	vertex_buffer_push_index (i);
	vertex_buffer_push_index (i + 1);
      }
      vertex_buffer_push_index (VertexBuffer.nvertex - 1);
      vertex_buffer_push_index (VertexBuffer.line_loop);
      type = 0;
      break;
    
    case GL_LINES:
      for (int i = VertexBuffer.lines; i < VertexBuffer.nvertex; i += 2) {
	vertex_buffer_push_index (i);
	vertex_buffer_push_index (i + 1);
      }
      type = 0;
      break;
    
    case GL_LINE_STRIP:
      for (int i = VertexBuffer.line_strip; i < VertexBuffer.nvertex - 1; i++) {
	vertex_buffer_push_index (i);
	vertex_buffer_push_index (i + 1);
      }
      type = 0;
      break;
    
    case GL_QUADS:
      for (int i = VertexBuffer.quads; i < VertexBuffer.nvertex; i += 4)
	for (int j = 1; j <= 2; j++) {
	  vertex_buffer_push_index (i);
	  vertex_buffer_push_index (i + j);
	  vertex_buffer_push_index (i + j + 1);
	}
      type = 1;
      break;
    
    case GL_POLYGON:
      for (int j = 1; j <= VertexBuffer.nvertex - VertexBuffer.polygon - 2;
	   j++) {
	vertex_buffer_push_index (VertexBuffer.polygon);
	vertex_buffer_push_index (VertexBuffer.polygon + j);
	vertex_buffer_push_index (VertexBuffer.polygon + j + 1);
      }
      type = 1;
      break;
    
    case GL_TRIANGLE_FAN:
      for (int i = VertexBuffer.fan + 1; i < VertexBuffer.nvertex - 1; i++) {
	vertex_buffer_push_index (VertexBuffer.fan);
	vertex_buffer_push_index (i);
	vertex_buffer_push_index (i + 1);
      }
      type = 1;
      break;
    
    default:
      break;
    }
    VertexBuffer.state = 0;
    if (VertexBuffer.type >= 0 && type >= 0) {
       // cannot mix lines and surfaces primitives
      assert (VertexBuffer.type == type);
    }
    else
      VertexBuffer.type = type;
  }
}

static void vertex_buffer_glColor3f (float r, float g, float b)
{
  glColor3f (r, g, b);
  if (VertexBuffer.color) {
    struct { float x, y, z; } color = {r, g, b}; // fixme: use r,g,b directly
    array_append (VertexBuffer.color, &color, 3*sizeof(float));
  }
}

static void vertex_buffer_glNormal3d (double nx, double ny, double nz)
{
  glNormal3d (nx, ny, nz);
  if (VertexBuffer.normal) {
    struct { float x, y, z; } normal = {nx, ny, nz};
    array_append (VertexBuffer.normal, &normal, 3*sizeof(float));
  }
}

static void vertex_buffer_glVertex3d (double x, double y, double z)
{
  glVertex3d (x, y, z);

  if (VertexBuffer.position) {
    if (VertexBuffer.dim < 3)
      VertexBuffer.dim = 3;
    float v[4] = {x, y, z, 1.};
    vector_multiply (v, VertexBuffer.modelview);
    array_append (VertexBuffer.position, v, 3*sizeof(float));
    VertexBuffer.nvertex++;
  }
}

static void vertex_buffer_glVertex2d (double x, double y)
{
  glVertex3d (x, y, 0.);
  
  if (VertexBuffer.position) {
    if (VertexBuffer.dim < 2)
      VertexBuffer.dim = 2;
    float v[4] = {x, y, 0, 1.};
    vector_multiply (v, VertexBuffer.modelview);
    array_append (VertexBuffer.position, v, 3*sizeof(float));
    VertexBuffer.nvertex++;
  }  
}

/**
Here we overload the default OpenGL commands, in order to call the
corresponding vertex buffer operations defined above. */

#define glBegin     vertex_buffer_glBegin
#define glEnd       vertex_buffer_glEnd
#define glVertex2d  vertex_buffer_glVertex2d
#define glVertex2f  vertex_buffer_glVertex2d
#define glVertex3d  vertex_buffer_glVertex3d
#define glColor3f   vertex_buffer_glColor3f
#define glNormal3d  vertex_buffer_glNormal3d

#endif
