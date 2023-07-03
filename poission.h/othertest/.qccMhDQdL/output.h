#ifndef BASILISK_HEADER_19
#define BASILISK_HEADER_19
#line 1 "/home/dahuanghhc/basilisk/src/output.h"
/**
# Output functions

## *output_field()*: Multiple fields interpolated on a regular grid (text format)

This function interpolates a *list* of fields on a *n+1 x n+1* regular
grid. The resulting data are written in text format in the file
pointed to by *fp*. The correspondance between column numbers and
variables is summarised in the first line of the file. The data are
written row-by-row and each row is separated from the next by a blank
line. This format is compatible with the *splot* command of *gnuplot*
i.e. one could use something like

~~~bash
gnuplot> set pm3d map
gnuplot> splot 'fields' u 1:2:4
~~~

The arguments and their default values are:

*list*
: list of fields to output. Default is *all*.

*fp*
: file pointer. Default is *stdout*.

*n*
: number of points along each dimension. Default is *N*.

*linear*
: use first-order (default) or bilinear interpolation. 

*box*
: the lower-left and upper-right coordinates of the domain to consider.
 Default is the entire domain. */

struct OutputField {
  scalar * list;
  FILE * fp;
  int n;
  bool linear;
  double box[2][2];
};

trace
void output_field (struct OutputField p)
{
  if (!p.list) p.list = all;
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = stdout;
  p.n++;
  if (p.box[0][0] == 0. && p.box[0][1] == 0. && 
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0;      p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }

  boundary (p.list);
  int len = list_len(p.list);
  double Delta = 0.999999*(p.box[1][0] - p.box[0][0])/(p.n - 1);
  int ny = (p.box[1][1] - p.box[0][1])/Delta + 1;
  double ** field = (double **) matrix_new (p.n, ny, len*sizeof(double));
  for (int i = 0; i < p.n; i++) {
    double x = Delta*i + p.box[0][0];
    for (int j = 0; j < ny; j++) {
      double y = Delta*j + p.box[0][1];
      if (p.linear) {
	int k = 0;
	for (scalar s in p.list)
	  field[i][len*j + k++] = interpolate (s, x, y);
      }
      else {
	Point point = locate (x, y);
	int k = 0;
	for (scalar s in p.list)
	  field[i][len*j + k++] = point.level >= 0 ? s[] : nodata;
      }
    }
  }

  if (pid() == 0) { // master
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], len*p.n*ny, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
@endif
    fprintf (p.fp, "# 1:x 2:y");
    int i = 3;
    for (scalar s in p.list)
      fprintf (p.fp, " %d:%s", i++, s.name);
    fputc('\n', p.fp);
    for (int i = 0; i < p.n; i++) {
      double x = Delta*i + p.box[0][0];
      for (int j = 0; j < ny; j++) {
	double y = Delta*j + p.box[0][1];
	//	map (x, y);
	fprintf (p.fp, "%g %g", x, y);
	int k = 0;
	for (scalar s in p.list)
	  fprintf (p.fp, " %g", field[i][len*j + k++]);
	fputc ('\n', p.fp);
      }
      fputc ('\n', p.fp);
    }
    fflush (p.fp);
  }
@if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, len*p.n*ny, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
@endif

  matrix_free (field);
}

/**
## *output_matrix()*: Single field interpolated on a regular grid (binary format)

This function writes a binary representation of a single field
interpolated on a regular *n x n* grid. The format is compatible with
the binary matrix format of gnuplot i.e. one could use

~~~bash
gnuplot> set pm3d map
gnuplot> splot 'matrix' binary u 2:1:3
~~~

The arguments and their default values are:

*f*
: a scalar field (compulsory).

*fp*
: file pointer. Default is *stdout*.

*n*
: number of points along each dimension. Default is *N*.

*linear*
: use first-order (default) or bilinear interpolation. */

struct OutputMatrix {
  scalar f;
  FILE * fp;
  int n;
  bool linear;
};

trace
void output_matrix (struct OutputMatrix p)
{
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = stdout;
  if (p.linear) {
    scalar f = p.f;
    boundary ({f});
  }
  float fn = p.n;
  float Delta = (float) L0/fn;
  fwrite (&fn, sizeof(float), 1, p.fp);
  for (int j = 0; j < p.n; j++) {
    float yp = (float) (Delta*j + X0 + Delta/2.);
    fwrite (&yp, sizeof(float), 1, p.fp);
  }
  for (int i = 0; i < p.n; i++) {
    float xp = (float) (Delta*i + X0 + Delta/2.);
    fwrite (&xp, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      float yp = (float)(Delta*j + Y0 + Delta/2.), v;
      if (p.linear)
	v = interpolate (p.f, xp, yp);
      else {
	Point point = locate (xp, yp);
	assert (point.level >= 0);
	v = val(p.f);
      }
      fwrite (&v, sizeof(float), 1, p.fp);
    }
  }
  fflush (p.fp);
}

/**
## Colormaps

Colormaps are arrays of (127) red, green, blue triplets. */

#define NCMAP 127

typedef void (* colormap) (double cmap[NCMAP][3]);

void jet (double cmap[NCMAP][3])
{
  for (int i = 0; i < NCMAP; i++) {
    cmap[i][0] = 
      i <= 46 ? 0. : 
      i >= 111 ? -0.03125*(i - 111) + 1. :
      i >= 78 ? 1. : 
      0.03125*(i - 46);
    cmap[i][1] = 
      i <= 14 || i >= 111 ? 0. : 
      i >= 79 ? -0.03125*(i - 111) : 
      i <= 46 ? 0.03125*(i - 14) : 
      1.;
    cmap[i][2] =
      i >= 79 ? 0. :
      i >= 47 ? -0.03125*(i - 79) :
      i <= 14 ? 0.03125*(i - 14) + 1.:
      1.;
  }
}

void cool_warm (double cmap[NCMAP][3])
{
  /* diverging cool-warm from:
   *  http://www.sandia.gov/~kmorel/documents/ColorMaps/CoolWarmFloat33.csv
   * see also:
   *  Diverging Color Maps for Scientific Visualization (Expanded)
   *  Kenneth Moreland
   */
  static double basemap[33][3] = {
    {0.2298057,   0.298717966, 0.753683153},
    {0.26623388,  0.353094838, 0.801466763},
    {0.30386891,  0.406535296, 0.84495867},
    {0.342804478, 0.458757618, 0.883725899},
    {0.38301334,  0.50941904,  0.917387822},
    {0.424369608, 0.558148092, 0.945619588},
    {0.46666708,  0.604562568, 0.968154911},
    {0.509635204, 0.648280772, 0.98478814},
    {0.552953156, 0.688929332, 0.995375608},
    {0.596262162, 0.726149107, 0.999836203},
    {0.639176211, 0.759599947, 0.998151185},
    {0.681291281, 0.788964712, 0.990363227},
    {0.722193294, 0.813952739, 0.976574709},
    {0.761464949, 0.834302879, 0.956945269},
    {0.798691636, 0.849786142, 0.931688648},
    {0.833466556, 0.860207984, 0.901068838},
    {0.865395197, 0.86541021,  0.865395561},
    {0.897787179, 0.848937047, 0.820880546},
    {0.924127593, 0.827384882, 0.774508472},
    {0.944468518, 0.800927443, 0.726736146},
    {0.958852946, 0.769767752, 0.678007945},
    {0.96732803,  0.734132809, 0.628751763},
    {0.969954137, 0.694266682, 0.579375448},
    {0.966811177, 0.650421156, 0.530263762},
    {0.958003065, 0.602842431, 0.481775914},
    {0.943660866, 0.551750968, 0.434243684},
    {0.923944917, 0.49730856,  0.387970225},
    {0.89904617,  0.439559467, 0.343229596},
    {0.869186849, 0.378313092, 0.300267182},
    {0.834620542, 0.312874446, 0.259301199},
    {0.795631745, 0.24128379,  0.220525627},
    {0.752534934, 0.157246067, 0.184115123},
    {0.705673158, 0.01555616,  0.150232812}	
  };
  
  for (int i = 0; i < NCMAP; i++) {
    double x = i*(32 - 1e-10)/(NCMAP - 1);
    int j = x; x -= j;
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (1. - x)*basemap[j][k] + x*basemap[j+1][k];
  }
}

void gray (double cmap[NCMAP][3])
{
  for (int i = 0; i < NCMAP; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = i/(NCMAP - 1.);
}

void randomap (double cmap[NCMAP][3])
{
  srand(0);
  for (int i = 0; i < NCMAP; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (noise() + 1.)/2.;
}

void blue_white_red (double cmap[NCMAP][3])
{
  for (int i = 0; i < (NCMAP + 1)/2; i++) {
    cmap[i][0] = i/((NCMAP - 1)/2.);
    cmap[i][1] = i/((NCMAP - 1)/2.);
    cmap[i][2] = 1.;
  }
  for (int i = 0; i < (NCMAP - 1)/2; i++) {
    cmap[i + (NCMAP + 1)/2][0] = 1.;
    cmap[i + (NCMAP + 1)/2][1] = cmap[(NCMAP - 3)/2 - i][1];
    cmap[i + (NCMAP + 1)/2][2] = cmap[(NCMAP - 3)/2 - i][1];
  }
}

/**
Given a colormap and a minimum and maximum value, this function
returns the red/green/blue triplet corresponding to *val*. */

typedef struct {
  unsigned char r, g, b;
} color;

color colormap_color (double cmap[NCMAP][3], 
		      double val, double min, double max)
{
  color c;
  if (val == nodata) {
    c.r = c.g = c.b = 0; // nodata is black
    return c;
  }
  int i;
  double coef;
  if (max != min)
    val = (val - min)/(max - min);
  else
    val = 0.;
  if (val <= 0.) i = 0, coef = 0.;
  else if (val >= 1.) i = NCMAP - 2, coef = 1.;
  else {
    i = val*(NCMAP - 1);
    coef = val*(NCMAP - 1) - i;
  }
  assert (i < NCMAP - 1);
  unsigned char * c1 = (unsigned char *) &c;
  for (int j = 0; j < 3; j++)
    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);
  return c;
}

/**
## Image/animation conversion

The open_image()/close_image() functions use pipes to convert PPM
images to other formats, including `.mp4`, `.ogv` and `.gif`
animations.

The functions check whether the 'ffmpeg' or 'convert' executables are
accessible, if they are not the conversion is disabled and the raw PPM
images are saved. An extra ".ppm" extension is added to the file name
to indicate that this happened. */

static const char * extension (const char * file, const char * ext) {
  int len = strlen(file);
  return len > 4 && !strcmp (file + len - 4, ext) ? file + len - 4 : NULL;
}

static const char * is_animation (const char * file) {
  const char * ext;
  if ((ext = extension (file, ".mp4")) ||
      (ext = extension (file, ".ogv")) ||
      (ext = extension (file, ".gif")))
    return ext;
  return NULL;
}

static struct {
  FILE ** fp;
  char ** names;
  int n;
} open_image_data = {NULL, NULL, 0};

static void open_image_cleanup()
{
  for (int i = 0; i < open_image_data.n; i++) {
    pclose (open_image_data.fp[i]);
    free (open_image_data.names[i]);
  }
  free (open_image_data.fp);
  free (open_image_data.names);
  open_image_data.fp = NULL;
  open_image_data.names = NULL;
  open_image_data.n = 0;
}

static FILE * open_image_lookup (const char * file)
{
  for (int i = 0; i < open_image_data.n; i++)
    if (!strcmp (file, open_image_data.names[i]))
      return open_image_data.fp[i];
  return NULL;
}

static bool which (const char * command)
{
  char * s = getenv ("PATH");
  if (!s)
    return false;
  char path[strlen(s) + 1];
  strcpy (path, s);
  s = strtok (path, ":");
  while (s) {
    char f[strlen(s) + strlen(command) + 2];
    strcpy (f, s);
    strcat (f, "/");
    strcat (f, command);
    FILE * fp = fopen (f, "r");
    if (fp) {
      fclose (fp);
      return true;
    }
    s = strtok (NULL, ":");
  }
  return false;
}

static FILE * ppm_fallback (const char * file, const char * mode)
{
  char filename[strlen(file) + 5];
  strcpy (filename, file);
  strcat (filename, ".ppm");
  FILE * fp = fopen (filename, mode);
  if (!fp) {
    perror (file);
#if _MPI
    MPI_Abort (MPI_COMM_WORLD, 1);
#endif
    exit (1);
  }
  return fp;
}

FILE * open_image (const char * file, const char * options)
{
  assert (pid() == 0);
  const char * ext;
  if ((ext = is_animation (file))) {
    FILE * fp = open_image_lookup (file);
    if (fp)
      return fp;

    int len = strlen ("ppm2???    ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "ppm2"); strcat (command, ext + 1);

    static int has_ffmpeg = -1;
    if (has_ffmpeg < 0) {
      if (which (command) && (which ("ffmpeg") || which ("avconv")))
	has_ffmpeg = true;
      else {
	fprintf (ferr,
		 "open_image(): cannot find '%s' or 'ffmpeg'/'avconv'\n"
		 "  falling back to raw PPM outputs\n", command);
	has_ffmpeg = false;
      }
    }
    if (!has_ffmpeg)
      return ppm_fallback (file, "a");

    static bool added = false;
    if (!added) {
      free_solver_func_add (open_image_cleanup);
      added = true;
    }      
    open_image_data.n++;
    qrealloc (open_image_data.names, open_image_data.n, char *);
    open_image_data.names[open_image_data.n - 1] = strdup (file);

    if (options) {
      strcat (command, " ");
      strcat (command, options);
    }
    strcat (command, !strcmp (ext, ".mp4") ? " " : " > ");
    strcat (command, file);
    qrealloc (open_image_data.fp, open_image_data.n, FILE *);
    return open_image_data.fp[open_image_data.n - 1] = popen (command, "w");
  }
  else { // !animation
    static int has_convert = -1;
    if (has_convert < 0) {
      if (which ("convert"))
	has_convert = true;
      else {
	fprintf (ferr,
		 "open_image(): cannot find 'convert'\n"
		 "  falling back to raw PPM outputs\n");
	has_convert = false;
      }
    }
    if (!has_convert)
      return ppm_fallback (file, "w");
    
    int len = strlen ("convert ppm:-   ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "convert ppm:- ");
    if (options) {
      strcat (command, options);
      strcat (command, " ");
    }
    strcat (command, file);
    return popen (command, "w");
  }
}

void close_image (const char * file, FILE * fp)
{
  assert (pid() == 0);
  if (is_animation (file)) {
    if (!open_image_lookup (file))
      fclose (fp);
  }
  else if (which ("convert"))
    pclose (fp);
  else
    fclose (fp);
}

/**
## *output_ppm()*: Portable PixMap (PPM) image output

Given a field, this function outputs a colormaped representation as a
[Portable PixMap](http://en.wikipedia.org/wiki/Netpbm_format) image.

If [ImageMagick](http://www.imagemagick.org/) is installed on the
system, this image can optionally be converted to any image format
supported by ImageMagick.

The arguments and their default values are:

*f*
: a scalar field (compulsory).

*fp*
: a file pointer. Default is stdout.

*n*
: number of pixels. Default is *N*.

*file*
: sets the name of the file used as output for
ImageMagick. This allows outputs in all formats supported by
ImageMagick. For example, one could use

~~~c
output_ppm (f, file = "f.png");
~~~

to get a [PNG](http://en.wikipedia.org/wiki/Portable_Network_Graphics)
image.

*min, max*
: minimum and maximum values used to define the
colorscale. By default these are set automatically using the *spread*
parameter. 

*spread*
: if not specified explicitly, *min* and *max* are set to the average
of the field minus (resp. plus) *spread* times the standard deviation.
By default *spread* is five. If negative, the minimum and maximum
values of the field are used.

*linear*
: whether to use bilinear or first-order interpolation. Default is 
first-order.

*box*
: the lower-left and upper-right coordinates of the domain to consider.
 Default is the entire domain.

*mask*
: if set, this field will be used to mask out (in black), the regions 
of the domain for which *mask* is negative. 

*map*
: the colormap: *jet*, *cool_warm* or *gray*. Default is *jet*.

*opt*
: options to pass to 'convert' or to the 'ppm2???' scripts (used
with *file*).
*/

struct OutputPPM {
  scalar f;
  FILE * fp;
  int n;
  char * file;
  double min, max, spread, z;
  bool linear;
  double box[2][2];
  scalar mask;
  colormap map;
  char * opt;
};

trace
void output_ppm (struct OutputPPM p)
{
  // default values
  if (p.n == 0) p.n = N;
  if (p.min == 0 && p.max == 0) {
    stats s = statsf (p.f);
    if (p.spread < 0.)
      p.min = s.min, p.max = s.max;
    else {
      double avg = s.sum/s.volume, spread = (p.spread ? p.spread : 5.)*s.stddev;
      p.min = avg - spread; p.max = avg + spread;
    }
  }
  if (p.box[0][0] == 0. && p.box[0][1] == 0. && 
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0;      p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }
  if (!p.map)
    p.map = jet;
  if (p.linear) {
    scalar f = p.f, mask = p.mask;
    if (mask.i)
      boundary ({f, mask});
    else
      boundary ({f});
  }
  
  double fn = p.n;
  double Delta = (p.box[1][0] - p.box[0][0])/fn;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;
  if (ny % 2) ny++;
  
  color ** ppm = (color **) matrix_new (ny, p.n, sizeof(color));
  double cmap[NCMAP][3];
  p.map (cmap);
  OMP_PARALLEL() {
    OMP(omp for schedule(static))
      for (int j = 0; j < ny; j++) {
	double yp = Delta*j + p.box[0][1] + Delta/2.;
	for (int i = 0; i < p.n; i++) {
	  double xp = Delta*i + p.box[0][0] + Delta/2., v;
	  if (p.mask.i) { // masking
	    if (p.linear) {
	      double m = interpolate (p.mask, xp, yp, p.z);
	      if (m < 0.)
		v = nodata;
	      else
		v = interpolate (p.f, xp, yp, p.z);
	    }
	    else {
	      Point point = locate (xp, yp, p.z);
	      if (point.level < 0 || val(p.mask) < 0.)
		v = nodata;
	      else
		v = val(p.f);
	    }
	  }
	  else if (p.linear)
	    v = interpolate (p.f, xp, yp, p.z);
	  else {
	    Point point = locate (xp, yp, p.z);
	    v = point.level >= 0 ? val(p.f) : nodata;
	  }
	  ppm[ny - 1 - j][i] = colormap_color (cmap, v, p.min, p.max);
	}
      }
  }
  
  if (pid() == 0) { // master
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, ppm[0], 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
		MPI_COMM_WORLD);
@endif
    if (!p.fp) p.fp = stdout;
    if (p.file)
      p.fp = open_image (p.file, p.opt);
    
    fprintf (p.fp, "P6\n%u %u 255\n", p.n, ny);
    fwrite (((void **) ppm)[0], sizeof(color), ny*p.n, p.fp);
    
    if (p.file)
      close_image (p.file, p.fp);
    else
      fflush (p.fp);
  }
@if _MPI
  else // slave
    MPI_Reduce (ppm[0], NULL, 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
		MPI_COMM_WORLD);
@endif
    
  matrix_free (ppm);
}

/**
## *output_grd()*: ESRI ASCII Grid format

The [ESRI GRD format](http://en.wikipedia.org/wiki/Esri_grid) is a
standard format for importing raster data into [GIS
systems](http://en.wikipedia.org/wiki/Geographic_information_system).

The arguments and their default values are:

*f*
: a scalar field (compulsory).

*fp*
: a file pointer. Default is stdout.

$\Delta$
: size of a grid element. Default is 1/N.

*linear*
: whether to use bilinear or first-order interpolation. Default is 
first-order.

*box*
: the lower-left and upper-right coordinates of the domain to consider.
 Default is the entire domain.

*mask*
: if set, this field will be used to mask out, the regions 
of the domain for which *mask* is negative. */

struct OutputGRD {
  scalar f;
  FILE * fp;
  double Delta;
  bool linear;
  double box[2][2];
  scalar mask;
};

trace
void output_grd (struct OutputGRD p)
{
  // default values
  if (!p.fp) p.fp = stdout;
  if (p.box[0][0] == 0. && p.box[0][1] == 0. && 
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0;      p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
    if (p.Delta == 0) p.Delta = L0/N;
  }
  if (p.linear) {
    scalar f = p.f, mask = p.mask;
    if (mask.i)
      boundary ({f, mask});
    else
      boundary ({f});
  }

  double Delta = p.Delta;
  int nx = (p.box[1][0] - p.box[0][0])/Delta;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;

  // header
  fprintf (p.fp, "ncols          %d\n", nx);
  fprintf (p.fp, "nrows          %d\n", ny);
  fprintf (p.fp, "xllcorner      %g\n", p.box[0][0]);
  fprintf (p.fp, "yllcorner      %g\n", p.box[0][1]);
  fprintf (p.fp, "cellsize       %g\n", Delta);
  fprintf (p.fp, "nodata_value   -9999\n");
  
  // data
  for (int j = ny-1; j >= 0; j--) {
    double yp = Delta*j + p.box[0][1] + Delta/2.;
    for (int i = 0; i < nx; i++) {
      double xp = Delta*i + p.box[0][0] + Delta/2., v;
      if (p.mask.i) { // masking
	if (p.linear) {
	  double m = interpolate (p.mask, xp, yp);
	  if (m < 0.)
	    v = nodata;
	  else
	    v = interpolate (p.f, xp, yp);
	}
	else {
	  Point point = locate (xp, yp);
	  if (point.level < 0 || val(p.mask) < 0.)
	    v = nodata;
	  else
	    v = val(p.f);
	}
      }
      else if (p.linear)
	v = interpolate (p.f, xp, yp);
      else {
	Point point = locate (xp, yp);
	v = point.level >= 0 ? val(p.f) : nodata;
      }
      if (v == nodata)
	fprintf (p.fp, "-9999 ");
      else
	fprintf (p.fp, "%f ", v);
    }
    fprintf (p.fp, "\n");
  }

  fflush (p.fp);
}

#if MULTIGRID

/**
## *output_gfs()*: Gerris simulation format

The function writes simulation data in the format used in
[Gerris](http://gfs.sf.net) simulation files. These files can be read
with GfsView.

The arguments and their default values are:

*fp*
: a file pointer. Default is *name* or stdout.

*list*
: a list of scalar fields to write. Default is *all*. 

*file*
: the name of the file to write to (mutually exclusive with *fp*).

*translate*
: whether to replace "well-known" Basilisk variables with their Gerris
equivalents.
*/

struct OutputGfs {
  FILE * fp;
  scalar * list;
  double t; // fixme: obsolete
  char * file;
  bool translate;
};

static char * replace (const char * input, int target, int with,
		       bool translate)
{
  if (translate) {
    if (!strcmp (input, "u.x"))
      return strdup ("U");
    if (!strcmp (input, "u.y"))
      return strdup ("V");
    if (!strcmp (input, "u.z"))
      return strdup ("W");
  }
  char * name = strdup (input), * i = name;
  while (*i != '\0') {
    if (*i == target)
      *i = with;
    i++;
  }
  return name;
}

trace
void output_gfs (struct OutputGfs p)
{
  char * fname = p.file;
  
@if _MPI
#if MULTIGRID_MPI
  not_mpi_compatible();
#endif // !MULTIGRID_MPI
  FILE * fp = p.fp;
  if (p.file == NULL) {
    long pid = getpid();
    MPI_Bcast (&pid, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    fname = qmalloc (80, char);
    snprintf (fname, 80, ".output-%ld", pid);
    p.fp = NULL;
  }
@endif // _MPI
  
  bool opened = false;
  if (p.fp == NULL) {
    if (fname == NULL)
      p.fp = stdout;
    else if (!(p.fp = fopen (fname, "w"))) {
      perror (fname);
      exit (1);
    }
    else
      opened = true;
  }
  
  scalar * list = p.list ? p.list : list_copy (all);

  restriction (list);
  fprintf (p.fp, 
	   "1 0 GfsSimulation GfsBox GfsGEdge { binary = 1"
	   " x = %g y = %g ",
	   0.5 + X0/L0, 0.5 + Y0/L0);
#if dimension == 3
  fprintf (p.fp, "z = %g ", 0.5 + Z0/L0);
#endif

  if (list != NULL && list[0].i != -1) {
    scalar s = list[0];
    char * name = replace (s.name, '.', '_', p.translate);
    fprintf (p.fp, "variables = %s", name);
    free (name);
    for (int i = 1; i < list_len(list); i++) {
      scalar s = list[i];
      if (s.name) {
	char * name = replace (s.name, '.', '_', p.translate);
	fprintf (p.fp, ",%s", name);
	free (name);
      }
    }
    fprintf (p.fp, " ");
  }
  fprintf (p.fp, "} {\n");
  fprintf (p.fp, "  Time { t = %g }\n", t);
  if (L0 != 1.)
    fprintf (p.fp, "  PhysicalParams { L = %g }\n", L0);
  fprintf (p.fp, "  VariableTracerVOF f\n");
  fprintf (p.fp, "}\nGfsBox { x = 0 y = 0 z = 0 } {\n");

@if _MPI
  long header;
  if ((header = ftell (p.fp)) < 0) {
    perror ("output_gfs(): error in header");
    exit (1);
  }
  int cell_size = sizeof(unsigned) + sizeof(double);
  for (scalar s in list)
    if (s.name)
      cell_size += sizeof(double);
  scalar index = new scalar;
  size_t total_size = header + (z_indexing (index, false) + 1)*cell_size;
@endif
  
  // see gerris/ftt.c:ftt_cell_write()
  //     gerris/domain.c:gfs_cell_write()
  foreach_cell() {
@if _MPI // fixme: this won't work when combining MPI and mask()
    if (is_local(cell))
@endif
    {
@if _MPI
      if (fseek (p.fp, header + index[]*cell_size, SEEK_SET) < 0) {
	perror ("output_gfs(): error while seeking");
	exit (1);
      }
@endif
      unsigned flags = 
	level == 0 ? 0 :
#if dimension == 1
	child.x == 1;
#elif dimension == 2
      child.x == -1 && child.y == -1 ? 0 :
	child.x == -1 && child.y ==  1 ? 1 :
	child.x ==  1 && child.y == -1 ? 2 : 
	3;
#else // dimension == 3
      child.x == -1 && child.y == -1 && child.z == -1  ? 0 :
	child.x == -1 && child.y == -1 && child.z ==  1  ? 1 :
	child.x == -1 && child.y ==  1 && child.z == -1  ? 2 : 
	child.x == -1 && child.y ==  1 && child.z ==  1  ? 3 : 
	child.x ==  1 && child.y == -1 && child.z == -1 ? 4 :
	child.x ==  1 && child.y == -1 && child.z ==  1 ? 5 :
	child.x ==  1 && child.y ==  1 && child.z == -1 ? 6 : 
	7;
#endif
      if (is_leaf(cell))
	flags |= (1 << 4);
      fwrite (&flags, sizeof (unsigned), 1, p.fp);
      double a = -1;
      fwrite (&a, sizeof (double), 1, p.fp);
      for (scalar s in list)
	if (s.name) {
	  if (s.v.x.i >= 0) {
	    // this is a vector component, we need to rotate from
	    // N-ordering (Basilisk) to Z-ordering (Gerris)
	    // fixme: this does not work for tensors
#if dimension >= 2
	    if (s.v.x.i == s.i) {
	      s = s.v.y;
	      a = is_local(cell) && s[] != nodata ? s[] : (double) DBL_MAX;
	    }
	    else if (s.v.y.i == s.i) {
	      s = s.v.x;
	      a = is_local(cell) && s[] != nodata ? - s[] : (double) DBL_MAX;
	    }
#endif
#if dimension >= 3
	    else
	      a = is_local(cell) && s[] != nodata ? s[] : (double) DBL_MAX;
#endif
	  }
	  else
	    a = is_local(cell) && s[] != nodata ? s[] : (double) DBL_MAX;
	  fwrite (&a, sizeof (double), 1, p.fp);
	}
    }
    if (is_leaf(cell))
      continue;
  }
  
@if _MPI
  delete ({index});
  if (!pid() && fseek (p.fp, total_size, SEEK_SET) < 0) {
    perror ("output_gfs(): error while finishing");
    exit (1);
  }
  if (!pid())
@endif  
    fputs ("}\n", p.fp);
  fflush (p.fp);

  if (!p.list)
    free (list);
  if (opened)
    fclose (p.fp);

@if _MPI
  if (p.file == NULL) {
    MPI_Barrier (MPI_COMM_WORLD);
    if (pid() == 0) {
      if (fp == NULL)
	fp = stdout;
      p.fp = fopen (fname, "r");
      size_t l;
      unsigned char buffer[8192];
      while ((l = fread (buffer, 1, 8192, p.fp)) > 0)
	fwrite (buffer, 1, l, fp);
      fflush (fp);
      remove (fname);
    }
    free (fname);
  }
@endif // _MPI
}

/**
## *dump()*: Basilisk snapshots

This function (together with *restore()*) can be used to dump/restore
entire simulations.

The arguments and their default values are:

*file*
: the name of the file to write to (mutually exclusive with *fp*). The
default is "dump".

*list*
: a list of scalar fields to write. Default is *all*. 

*fp*
: a file pointer. Default is stdout.

*unbuffered*
: whether to use a file buffer. Default is false.
*/

struct Dump {
  char * file;
  scalar * list;
  FILE * fp;
  bool unbuffered;
};

struct DumpHeader {
  double t;
  long len;
  int i, depth, npe, version;
  coord n;
};

static const int dump_version =
  // 161020
  170901;

static scalar * dump_list (scalar * lista)
{
  scalar * list = is_constant(cm) ? NULL : list_concat ({cm}, NULL);
  for (scalar s in lista)
    if (!s.face && !s.nodump && s.i != cm.i)
      list = list_add (list, s);
  return list;
}

static void dump_header (FILE * fp, struct DumpHeader * header, scalar * list)
{
  if (fwrite (header, sizeof(struct DumpHeader), 1, fp) < 1) {
    perror ("dump(): error while writing header");
    exit (1);
  }
  for (scalar s in list) {
    unsigned len = strlen(s.name);
    if (fwrite (&len, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing len");
      exit (1);
    }
    if (fwrite (s.name, sizeof(char), len, fp) < len) {
      perror ("dump(): error while writing s.name");
      exit (1);
    }
  }
  double o[4] = {X0,Y0,Z0,L0};
  if (fwrite (o, sizeof(double), 4, fp) < 4) {
    perror ("dump(): error while writing coordinates");
    exit (1);
  }
}

@if !_MPI
trace
void dump (struct Dump p)
{
  FILE * fp = p.fp;
  char def[] = "dump", * file = p.file ? p.file : p.fp ? NULL : def;

  char * name = NULL;
  if (file) {
    name = (char *) malloc (strlen(file) + 2);
    strcpy (name, file);
    if (!p.unbuffered)
      strcat (name, "~");
    if ((fp = fopen (name, "w")) == NULL) {
      perror (name);
      exit (1);
    }
  }
  assert (fp);
  
  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar size[];
  scalar * list = list_concat ({size}, dlist); free (dlist);
  struct DumpHeader header = { t, list_len(list), iter, depth(), npe(),
			       dump_version };
  dump_header (fp, &header, list);
  
  subtree_size (size, false);
  
  foreach_cell() {
    unsigned flags = is_leaf(cell) ? leaf : 0;
    if (fwrite (&flags, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing flags");
      exit (1);
    }
    for (scalar s in list)
      if (fwrite (&s[], sizeof(double), 1, fp) < 1) {
	perror ("dump(): error while writing scalars");
	exit (1);
      }
    if (is_leaf(cell))
      continue;
  }
  
  free (list);
  if (file) {
    fclose (fp);
    if (!p.unbuffered)
      rename (name, file);
    free (name);
  }
}
@else // _MPI
trace
void dump (struct Dump p)
{
  FILE * fp = p.fp;
  char def[] = "dump", * file = p.file ? p.file : p.fp ? NULL : def;

  if (fp != NULL || file == NULL) {
    fprintf (ferr, "dump(): must specify a file name when using MPI\n");
    exit(1);
  }

  char name[strlen(file) + 2];
  strcpy (name, file);
  if (!p.unbuffered)
    strcat (name, "~");
  FILE * fh = fopen (name, "w");
  if (fh == NULL) {
    perror (name);
    exit (1);    
  }

  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar size[];
  scalar * list = list_concat ({size}, dlist); free (dlist);
  struct DumpHeader header = { t, list_len(list), iter, depth(), npe(),
			       dump_version };

#if MULTIGRID_MPI
  for (int i = 0; i < dimension; i++)
    (&header.n.x)[i] = mpi_dims[i];
  MPI_Barrier (MPI_COMM_WORLD);
#endif

  if (pid() == 0)
    dump_header (fh, &header, list);
  
  scalar index = {-1};
  
  index = new scalar;
  z_indexing (index, false);
  int cell_size = sizeof(unsigned) + header.len*sizeof(double);
  int sizeofheader = sizeof(header) + 4*sizeof(double);
  for (scalar s in list)
    sizeofheader += sizeof(unsigned) + sizeof(char)*strlen(s.name);
  long pos = pid() ? 0 : sizeofheader;
  
  subtree_size (size, false);
  
  foreach_cell() {
    // fixme: this won't work when combining MPI and mask()
    if (is_local(cell)) {
      long offset = sizeofheader + index[]*cell_size;
      if (pos != offset) {
	fseek (fh, offset, SEEK_SET);
	pos = offset;
      }
      unsigned flags = is_leaf(cell) ? leaf : 0;
      fwrite (&flags, 1, sizeof(unsigned), fh);
      for (scalar s in list)
	fwrite (&s[], 1, sizeof(double), fh);
      pos += cell_size;
    }
    if (is_leaf(cell))
      continue;
  }

  delete ({index});
  
  free (list);
  fclose (fh);
  if (!p.unbuffered && pid() == 0)
    rename (name, file);
}
@endif // _MPI

trace
bool restore (struct Dump p)
{
  FILE * fp = p.fp;
  char * file = p.file;
  if (file && (fp = fopen (file, "r")) == NULL)
    return false;
  assert (fp);

  struct DumpHeader header;  
  if (fread (&header, sizeof(header), 1, fp) < 1) {
    fprintf (ferr, "restore(): error: expecting header\n");
    exit (1);
  }

#if TREE
  init_grid (1);
  foreach_cell() {
    cell.pid = pid();
    cell.flags |= active;
  }
  tree->dirty = true;
#else // multigrid
#if MULTIGRID_MPI
  if (header.npe != npe()) {
    fprintf (ferr,
	     "restore(): error: the number of processes don't match:"
	     " %d != %d\n",
	     header.npe, npe());
    exit (1);
  }
  dimensions (header.n.x, header.n.y, header.n.z);
  double n = header.n.x;
  int depth = header.depth;
  while (n > 1)
    depth++, n /= 2;
  init_grid (1 << depth);
#else // !MULTIGRID_MPI
  init_grid (1 << header.depth);
#endif
#endif // multigrid

  bool restore_all = (p.list == all);
  scalar * list = dump_list (p.list ? p.list : all);
  if (header.version == 161020) {
    if (header.len - 1 != list_len (list)) {
      fprintf (ferr,
	       "restore(): error: the list lengths don't match: "
	       "%ld (file) != %d (code)\n",
	       header.len - 1, list_len (list));
      exit (1);
    }
  }
  else { // header.version != 161020
    if (header.version != dump_version) {
      fprintf (ferr,
	       "restore(): error: file version mismatch: "
	       "%d (file) != %d (code)\n",
	       header.version, dump_version);
      exit (1);
    }
    
    scalar * input = NULL;
    for (int i = 0; i < header.len; i++) {
      unsigned len;
      if (fread (&len, sizeof(unsigned), 1, fp) < 1) {
	fprintf (ferr, "restore(): error: expecting len\n");
	exit (1);
      }
      char name[len + 1];
      if (fread (name, sizeof(char), len, fp) < 1) {
	fprintf (ferr, "restore(): error: expecting s.name\n");
	exit (1);
      }
      name[len] = '\0';

      if (i > 0) { // skip subtree size
	bool found = false;
	for (scalar s in list)
	  if (!strcmp (s.name, name)) {
	    input = list_append (input, s);
	    found = true; break;
	  }
	if (!found) {
	  if (restore_all) {
	    scalar s = new scalar;
	    free (s.name);
	    s.name = strdup (name);
	    input = list_append (input, s);
	  }
	  else
	    input = list_append (input, (scalar){INT_MAX});
	}
      }
    }
    free (list);
    list = input;

    double o[4];
    if (fread (o, sizeof(double), 4, fp) < 4) {
      fprintf (ferr, "restore(): error: expecting coordinates\n");
      exit (1);
    }
    origin (o[0], o[1], o[2]);
    size (o[3]);
  }

#if MULTIGRID_MPI
  long cell_size = sizeof(unsigned) + header.len*sizeof(double);
  long offset = pid()*((1 << dimension*(header.depth + 1)) - 1)/
    ((1 << dimension) - 1)*cell_size;
  if (fseek (fp, offset, SEEK_CUR) < 0) {
    perror ("restore(): error while seeking");
    exit (1);
  }
#endif // MULTIGRID_MPI
  
  scalar * listm = is_constant(cm) ? NULL : (scalar *){fm};
#if TREE && _MPI
  restore_mpi (fp, list);
#else
  foreach_cell() {
    unsigned flags;
    if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
      fprintf (ferr, "restore(): error: expecting 'flags'\n");
      exit (1);
    }
    // skip subtree size
    fseek (fp, sizeof(double), SEEK_CUR);
    for (scalar s in list) {
      double val;
      if (fread (&val, sizeof(double), 1, fp) != 1) {
	fprintf (ferr, "restore(): error: expecting a scalar\n");
	exit (1);
      }
      if (s.i != INT_MAX)
	s[] = val;
    }
    if (!(flags & leaf) && is_leaf(cell))
      refine_cell (point, listm, 0, NULL);
    if (is_leaf(cell))
      continue;
  }
  for (scalar s in all)
    s.dirty = true;
#endif
  
  scalar * other = NULL;
  for (scalar s in all)
    if (!list_lookup (list, s) && !list_lookup (listm, s))
      other = list_append (other, s);
  reset (other, 0.);
  free (other);
  
  free (list);
  if (file)
    fclose (fp);

  // the events are advanced to catch up with the time  
  while (iter < header.i && events (false))
    iter = inext;
  events (false);
  while (t < header.t && events (false))
    t = tnext;
  t = header.t;
  events (false);
  
  return true;
}

#endif // MULTIGRID

#endif
