#ifndef BASILISK_HEADER_4
#define BASILISK_HEADER_4
#line 1 "/home/dahuanghhc/basilisk/src/display.h"
/**
# Server-side display

This file implements the server-side of the interactive display of a
running Basilisk code. The client-side is typically done by the
[javascript implementation](jview/README).

This works using two main components:

* [The wsServer WebSocket library](wsServer/README.md)
* [Vertex buffers](vertexbuffer.h)

The initial state of the display can be controlled using the DISPLAY
macro, following these rules:

* `-DDISPLAY=1, -DDISPLAY`: play controls, start running immediately.
* `-DDISPLAY=-1`: play controls, initially paused.
* `-DDISPLAY=0`: no play controls, start running immediately.
* `#include "display.h"`: play controls, initially paused.

This can also be changed by setting `display_play = 0;` in `main()`. 

## Display port

Multiple servers can run simultaneously on the same machine but must
use different ports. When connecting, the server looks for a free port
in the default range (typically 7100 to 7200, see the DISPLAY_RANGE
macro below). The corresponding connection information (including the
URL of the [javascript interface](/src/jview/README)) is written in
the `display.html` file, which can be used to open the connection on
the client side (i.e. the web browser).

## Display control

The values of integer or double variables can be controlled
interactively using the `display_control()` macro, for example:

~~~literatec
display_control (value, 0, 100, "Short Name", "Long tooltip");
~~~

where the first argument is the name of the variable and the numbers
indicate the (optional) range (min and max values). The (optional)
name and tooltip are used to generate the [Basilisk
View](jview/README) interface. The only compulsory argument is the
variable name. 

## Default display

Default [drawing function calls](draw.h), which will be sent to an empty
[javascript client](jview/README), can be setup using something like:

~~~literatec
display ("squares (color = 'u.x', spread = -1);");
~~~

where the single quotes will be expanded to doubles quotes when
interpreted. Do not forget the trailing semi-column.

Note that these function calls will be appended to existing ones. To
overwrite existing calls, use:

~~~literatec
display ("squares (color = 'u.x', spread = -1);", true);
~~~
*/

#ifndef DISPLAY_JS
# define DISPLAY_JS "http://basilisk.fr/three.js/editor/index.html"
#endif

#ifndef DISPLAY_HOST
# define DISPLAY_HOST "localhost"
#endif

#ifndef DISPLAY_RANGE
# define DISPLAY_RANGE "7100:7200"
#endif

#if 1
# define debug(...)
#else
# define debug(...) fprintf (qerr, __VA_ARGS__), fflush(qerr)
#endif

#include <netdb.h>
#include <wsServer/include/ws.h>
#pragma autolink -L$BASILISK/wsServer -lws

#include "view.h"
#include "khash.h"

typedef struct {
  int fd;
  int iter;
} DisplayClient;

KHASH_MAP_INIT_STR(strhash, DisplayClient *)

static struct {
  khash_t(strhash) * objects;
  int sock, port;
  char * error;
  Array * controls;
} Display = { .sock = -1 };

static void display_display()
{
  debug ("**************************\n");
  for (khiter_t k = kh_begin (Display.objects); k != kh_end (Display.objects);
       ++k)
    if (kh_exist (Display.objects, k)) {
      debug ("%s", kh_key (Display.objects, k));
      DisplayClient * client = kh_value (Display.objects, k);
      while (client->fd >= 0) {
	debug (" %d/%d", client->fd, client->iter);
	client++;
      }
      debug ("\n");
    }
  debug ("--------------------------\n");
}

static char * read_file_into_buffer (FILE * fp)
{
  if (fseek (fp, 0L, SEEK_END) < 0)
    return NULL;
  long bufsize = ftell (fp);
  if (bufsize <= 0)
    return NULL;
  char * buf = malloc (sizeof(char)*(bufsize + 1));
  if (fseek (fp, 0L, SEEK_SET) < 0) {
    free (buf);
    return NULL;
  }  
  size_t newLen = fread (buf, sizeof(char), bufsize, fp);
  buf[newLen] = '\0'; /* Just to be safe. */
  return buf;
}

static void display_command (const char * command)
{
  debug ("display_command (%s)\n", command);
  
  vertex_buffer_setup();
  VertexBuffer.vertex = 0; // fixme

  // Temporarily redirect stderr to catch errors
  int bak = -1;
  if (pid() == 0) {
    fflush (stderr);
    bak = dup (2);
    FILE * fp = tmpfile();
    dup2 (fileno (fp), 2);
    fclose (fp);
  }

  char * line = strdup (command);
  bool status = process_line (line);
  free (line);
  
  free (Display.error);
  Display.error = NULL;
  if (status) {
    if (VertexBuffer.type < 0)
      VertexBuffer.type = 0;
  }
  else { // An error occured
    if (pid() == 0) {
      fflush (stderr);
      FILE * fp = fdopen (2, "r");
      Display.error = read_file_into_buffer (fp);
      int len = Display.error ? strlen (Display.error) : 0;
      if (len > 0 && Display.error[len - 1] == '\n')
	Display.error[len - 1] = '\0';
      fclose(fp);
    }
    else
      Display.error = strdup ("error (on slave process)");
    VertexBuffer.type = - 1;
  }

  if (pid() == 0) {
    // Restore stderr to its previous value
    dup2 (bak, 2);
    close (bak);
  }
    
  if (VertexBuffer.type < 0)
    debug ("error: '%s'\n", Display.error);
  else {
    if (pid() > 0) {
      if (VertexBuffer.normal->len < VertexBuffer.position->len)
	VertexBuffer.normal->len = 0;
      if (VertexBuffer.color->len < VertexBuffer.position->len)
	VertexBuffer.color->len = 0;
    }
    debug ("position: %ld, normal: %ld, color: %ld, index: %ld, type: %d\n", 
	   VertexBuffer.position->len,
	   VertexBuffer.normal->len,
	   VertexBuffer.color->len,
	   VertexBuffer.index->len, VertexBuffer.type);
  }
}

static int ws_send_array (int fd, Array * a, int status, int type,
			  unsigned int * shift)
{
@if _MPI
  if (pid() == 0) {
    void * p = NULL;
    long len;
    for (int pe = 0; pe < npe(); pe++) {
      if (pe == 0)
	p = a->p, len = a->len;
      else {
	MPI_Status status;
	MPI_Recv (&len, 1, MPI_LONG, pe, 22, MPI_COMM_WORLD, &status);
	if (len > 0) {
	  p = malloc (len);
	  MPI_Recv (p, len, MPI_BYTE, pe, 23, MPI_COMM_WORLD, &status);
	}
	else
	  p = NULL;
      }
      if (type == 0) // position
	shift[pe] = (pe > 0 ? shift[pe - 1] : 0) + len/(3*sizeof(float));
      else if (type == 1 && pe > 0) // index
	for (unsigned int i = 0; i < len/sizeof(unsigned int); i++)
	  ((unsigned int *) p)[i] += shift[pe - 1];
      if (status >= 0 && len > 0 && ws_send (fd, p, len) < len)
	status = -1;
      if (pe > 0)
	free (p);
    }
  }
  else { // pid() > 0
    MPI_Send (&a->len, 1, MPI_LONG, 0, 22, MPI_COMM_WORLD);
    if (a->len > 0)
      MPI_Send (a->p, a->len, MPI_BYTE, 0, 23, MPI_COMM_WORLD);
  }
@else
  if (status >= 0 && a->len > 0 && ws_send (fd, a->p, a->len) < a->len)
    status = -1;
@endif
  return status;
}

static int display_send (const char * command, int fd)
{
  int status = 0;

  debug ("sending '%s' to %d\n", command, fd);
  
  unsigned int commandlen = strlen (command);
  unsigned int errorlen = Display.error ? strlen (Display.error) : 0;
  
  int paddedlen = 4*ceil(commandlen/4.);
  size_t len = 2*sizeof(unsigned int) + paddedlen;

  unsigned int lens[] = {VertexBuffer.position->len,
			 VertexBuffer.normal->len,
			 VertexBuffer.color->len,
			 VertexBuffer.index->len}, glens[4];
  int type = VertexBuffer.type, gtype;
@if _MPI
  MPI_Reduce (lens, glens, 4, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Allreduce (&type, &gtype, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
@else
  for (int i = 0; i < 4; i++) glens[i] = lens[i];
  gtype = type;
@endif
  
  if (gtype < 0)
    len += errorlen;
  else
    len += 2*sizeof(int) +
      4*sizeof (unsigned int) + glens[0] + glens[1] + glens[2] + glens[3];

  if (pid() == 0) {
    if (ws_sendframe_init (fd, len, false, WS_FR_OP_BIN) < 0 ||
	ws_send (fd, &commandlen, sizeof(unsigned int)) < sizeof(unsigned int) ||
	ws_send (fd, &errorlen, sizeof(unsigned int)) < sizeof(unsigned int) ||
	ws_send (fd, command, commandlen) < commandlen)
      status = -1;
    
    // padding to a multiple of four
    for (int i = 0; i < paddedlen - commandlen && status >= 0; i++) {
      char c = '\0';
      if (ws_send (fd, &c, 1) < 1)
	status = -1;
    }
  }

  if (gtype < 0) {
    if (pid() == 0 && status >= 0 &&
	ws_send (fd, Display.error, errorlen) < errorlen)
      status = -1;
  }
  else {
    if (pid() == 0 && status >= 0 &&
	(ws_send (fd, &VertexBuffer.dim, sizeof(int)) < sizeof(int) ||
	 ws_send (fd, &gtype, sizeof(int)) < sizeof(int) ||
	 ws_send (fd, glens, 4*sizeof (unsigned int)) < 4*sizeof (unsigned int)))
      status = -1;
    unsigned int * shift = malloc (sizeof(unsigned int)*npe());
    status = ws_send_array (fd, VertexBuffer.position, status, 0, shift);
    status = ws_send_array (fd, VertexBuffer.normal, status, -1, shift);
    status = ws_send_array (fd, VertexBuffer.color, status, -1, shift);
    status = ws_send_array (fd, VertexBuffer.index, status, 1, shift);
    free (shift);
  }
  
@if _MPI
  MPI_Bcast (&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
@endif
  
  return status;
}

static void display_add (const char * command, int fd)
{
  debug ("adding '%s'\n", command);
  khiter_t k = kh_get (strhash, Display.objects, command);
  if (k == kh_end (Display.objects)) {
    int ret;
    k = kh_put (strhash, Display.objects, strdup (command), &ret);
    DisplayClient * client = malloc (sizeof(DisplayClient));
    client->fd = -1;
    kh_value (Display.objects, k) = client;
  }
  DisplayClient * clients = kh_value (Display.objects, k), * client = clients;
  int len = 0;
  while (client->fd >= 0) {
    if (client->fd == fd)
      fd = -1;
    client++, len++;
  }
  if (fd >= 0) { // not already in the array
    kh_value (Display.objects, k) = clients =
      realloc (clients, (len + 2)*sizeof (DisplayClient));
    clients[len].fd = fd;
    clients[len].iter = -1;
    clients[len + 1].fd = -1;
  }
  display_display();
}

#define JSON_BUILD(...) len += snprintf (build + len, 4096, __VA_ARGS__),     \
    build = realloc (build, len + 4096)

static void array_remove (khiter_t k, int fd)
{
  DisplayClient * clients = kh_value (Display.objects, k), * client = clients;
  int i = -1, len = 0;
  while (client->fd >= 0) {
    if (client->fd == fd) {
      if (i != -1)
	debug ("array_remove(): error! found multiple %d in '%s'\n",
	       fd, kh_key (Display.objects, k));
      i = len;
    }
    client++, len++;
  }
  if (i < 0)
    debug ("array_remove(): error! could not find %d in '%s'\n",
	   fd, kh_key (Display.objects, k));
  else if (len == 1) {
    free ((void *) kh_key (Display.objects, k));
    free ((void *) kh_value (Display.objects, k));
    kh_del (strhash, Display.objects, k);
  }
  else
    for (int j = i; j < len; j++)
      clients[j] = clients[j + 1];
}

static void display_remove (const char * command, int fd)
{
  debug ("removing '%s'\n", command);
  khiter_t k = kh_get (strhash, Display.objects, command);
  if (k == kh_end (Display.objects))
    debug ("display_remove(): error! could not find '%s' (%d)\n",
	   command, fd);
  else
    array_remove (k, fd);
  display_display();
}

typedef struct {
  char * name, * tooltip;
  void * ptr;
  double min, max;
  int size;
} DisplayControl;

static char * display_control_json()
{
  char * build = malloc (4096);
  int len = 0;
  JSON_BUILD ("#{");
  char sep = ' ';
  DisplayControl * d = Display.controls->p;
  for (int i = 0; i < Display.controls->len/sizeof(DisplayControl); i++, d++) {
    JSON_BUILD ("%c\n  \"%s\": { ", sep, d->name); sep = ',';
    if (d->tooltip)
      JSON_BUILD ("\"tooltip\": \"%s\", ", d->tooltip);
    switch (d->size) {
    case 4:
      JSON_BUILD ("\"type\": \"int\", \"value\": %d, \"min\": %g, \"max\": %g",
		  *((int *)d->ptr), d->min, d->max);
      break;
    case 8:
      JSON_BUILD ("\"type\": \"double\", \"value\": %g, "
		  "\"min\": %g, \"max\": %g",
		  *((double *)d->ptr), d->min, d->max);
      break;
    default:
      assert (false);
    }
    JSON_BUILD (" }");
  }
  JSON_BUILD ("}");
  return build;
}

static DisplayControl * display_control_lookup (const char * name)
{
  DisplayControl * d = Display.controls->p;
  for (int i = 0; i < Display.controls->len/sizeof(DisplayControl); i++, d++)
      if (!strcmp (d->name, name))
	return d;
  return NULL;
}

struct _DisplayControl {
  void * ptr;
  double min, max;
  char * name, * tooltip, * ptr_name;
  int size;
};

void display_control_internal (struct _DisplayControl p)
{
  DisplayControl d;
  if (!p.name)
    p.name = p.ptr_name;

  if (display_control_lookup (p.name))
    return;
    
  d.name = strdup (p.name);
  d.tooltip = p.tooltip ? strdup (p.tooltip) : NULL;
  d.ptr = p.ptr;
  d.size = p.size;
  d.min = p.min, d.max = p.max;
  array_append (Display.controls, &d, sizeof (DisplayControl));

  if (pid() == 0) {
    char * controls = display_control_json();
    ws_sendframe_txt (0, controls, true);
    free (controls);
  }
}

#undef display_control
#define display_control(val, ...) display_control_internal \
  (&(val), __VA_ARGS__, size = sizeof(val), ptr_name = #val)

static void display_control_update (const char * command, int fd)
{
  char * s = strdup (command), * s1 = strchr (command, ':');
  *s1++ = '\0';
  DisplayControl * d = display_control_lookup (command);
  if (d == NULL)
    debug ("display_control_update(): error! could not find '%s' (%d)\n",
	   command, fd);
  else {
    debug ("display_control_update (%s) = %s\n", command, s1);
    double val = atof(s1);
    if (d->max > d->min)
      val = clamp (val, d->min, d->max);
    switch (d->size) {
    case 4: *((int *)d->ptr) = val; break;
    case 8: *((double *)d->ptr) = val; break;
    default: assert (false);
    }

    if (pid() == 0) {
      char * controls = display_control_json();
      ws_sendframe_txt (- fd, controls, true);
      free (controls);
    }
  }
  free (s);
}

static char * bview_interface_json()
{
  char * build = malloc (4096);
  int len = 0;

  JSON_BUILD ("{\n");
  
  char p[4096] = {0};
  int i = 0;
  while (bview_interface[i].json) {
    JSON_BUILD ("%s", i ? ",\n" : "");
    len += bview_interface[i].json (p, build + len, 4096);
    build = realloc (build, len + 4096);
    JSON_BUILD ("\n");
    i++;
  }
  JSON_BUILD ("}");
  return build;
}

void display_onclose (int fd)
{  
  debug ("closing %d\n", fd);
  for (khiter_t k = kh_begin (Display.objects); k != kh_end (Display.objects);
       ++k)
    if (kh_exist (Display.objects, k))
      array_remove (k, fd);
  display_display();
}

void display_onmessage (int fd, const char * msg, size_t size, int type)
{
  if (type == WS_FR_OP_TXT) {
    if (!msg)
      fprintf (stderr, "error receiving data on websocket\n");
    else switch (msg[0]) {
      case '+': display_add (msg + 1, fd); break;
      case '-': display_remove (msg + 1, fd); break;
      case '#': display_control_update (msg + 1, fd); break;
      default: fprintf (stderr,
			"display_onmessage: error: unknown message type '%s'\n",
			msg);
	break;
      }
  }
  else
    fprintf (stderr, "display_onmessage: error: unexpected message type '%d'\n",
	     type);
}

void display_onopen (int fd)
{
  char * interface = bview_interface_json();
  char * controls = display_control_json();
  int status = 0;

  if (pid() == 0)
    if (ws_sendframe_txt (fd, interface, false) < 0 ||
	ws_sendframe_txt (fd, controls, false) < 0 ||
	(display_defaults && ws_sendframe_txt (fd, display_defaults, false) < 0))
      status = -1;
  free (interface);
  free (controls);
  
@if _MPI
  MPI_Bcast (&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
@endif

  debug ("open %d status %d\n", fd, status);
  
  if (status < 0) {
    display_onclose (fd);
    if (pid() == 0)
      close (fd);
  }
}

static void display_update (int i)
{
  for (khiter_t k = kh_begin (Display.objects); k != kh_end (Display.objects);
       ++k)
    if (kh_exist (Display.objects, k)) {
      DisplayClient * client = kh_value (Display.objects, k);
      while (client->fd >= 0) {
	if (client->iter < i)
	  break;
	client++;
      }
      if (client->fd >= 0) { // at least one client needs update
	const char * command = kh_key (Display.objects, k);
	display_command (command);
	client = kh_value (Display.objects, k);
	while (client->fd >= 0) {
	  if (client->iter < i) {
	      client->iter = i;
	      if (display_send (command, client->fd) < 0) {
		debug ("error sending '%s' to '%d'\n", command, client->fd);
		if (pid() == 0)
		  close (client->fd);
		display_onclose (client->fd);
		if (!kh_exist (Display.objects, k))
		  break;
	      }
	      else
		client++;
	  }
	  else
	    client++;
	}
	vertex_buffer_free();
	if (Display.error && kh_exist (Display.objects, k)) {
	  free ((void *) kh_key (Display.objects, k));
	  free ((void *) kh_value (Display.objects, k));
	  kh_del (strhash, Display.objects, k);
	}
      }
    }
}

@if _MPI
static Array * pack_messages (struct ws_message * messages)
{
  struct ws_message * msg = messages;
  Array * packed = array_new();
  while (msg && msg->fd >= 0) {
    array_append (packed, msg, sizeof(struct ws_message));
    array_append (packed, msg->msg, msg->size);
    msg++;
  }
  return packed;
}

static struct ws_message * unpack_messages (Array * packed)
{
  Array * array = array_new();
  char * p = packed->p;
  while (p - (char *) packed->p < packed->len) {
    struct ws_message * msg = (struct ws_message *) p;
    msg->msg = sysmalloc (msg->size + 1);
    p += sizeof(struct ws_message);
    memcpy (msg->msg, p, msg->size);
    msg->msg[msg->size] = '\0';
    array_append (array, msg, sizeof(struct ws_message));
    p += msg->size;
  }
  struct ws_message msg = {-1};
  array_append (array, &msg, sizeof(struct ws_message));
  struct ws_message * messages = array->p;
  free (array);
  return messages;
}
@endif
		 
int display_poll (int timeout)
{
  struct ws_message * messages = NULL, * msg;
  int nmsg = 0;
  if (pid() == 0) {
    msg = messages = ws_socket_poll (Display.sock, timeout);
    while (msg && msg->fd >= 0) msg++, nmsg++;
  }

@if _MPI
  Array * packed;
  if (pid() == 0)
    packed = pack_messages (messages);
  else
    packed = calloc (1, sizeof (Array));
  MPI_Bcast (&packed->len, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  if (packed->len > 0) {
    if (pid() != 0)
      packed->p = malloc (packed->len);
    MPI_Bcast (packed->p, packed->len, MPI_BYTE, 0, MPI_COMM_WORLD);
    if (pid() != 0)
      messages = unpack_messages (packed);
  }
  array_free (packed);
@endif

  msg = messages;
  nmsg = 0;
  while (msg && msg->fd >= 0) {
    switch (msg->type) {

    case WS_FR_OP_OPEN:
      display_onopen (msg->fd); break;

    case WS_FR_OP_TXT: case WS_FR_OP_BIN:
      display_onmessage (msg->fd, msg->msg, msg->size, msg->type);
      break;

    case WS_FR_OP_CLSE:
      display_onclose (msg->fd); break;
      
    default:
      assert (false);
      
    }
    sysfree (msg->msg);
    msg++, nmsg++;
  }
  sysfree (messages);
  return nmsg;
}

void display_url (FILE * fp)
{
  char hostname[1024];
  hostname[1023] = '\0';
  gethostname (hostname, 1023);
  struct hostent * h = gethostbyname (hostname);
  fprintf (fp, DISPLAY_JS "?ws://%s:%d", h->h_name, Display.port);
}

int display_usage = 20; // use 20% of runtime, maximum

#ifdef DISPLAY
 // negative: pause, zero: play, positive: step
int display_play = DISPLAY < 0 ? -1 : 0;
#else
int display_play = -1;
#endif

static void display_destroy()
{
  for (khiter_t k = kh_begin (Display.objects); k != kh_end (Display.objects);
       ++k)
    if (kh_exist (Display.objects, k)) {
      free ((void *) kh_key (Display.objects, k));
      free ((void *) kh_value (Display.objects, k));
    }
  kh_destroy (strhash, Display.objects);

  DisplayControl * d = Display.controls->p;
  for (int i = 0; i < Display.controls->len/sizeof(DisplayControl); i++, d++) {
    free (d->name);
    free (d->tooltip);
  }
  array_free (Display.controls);
  
  if (pid() == 0) {
    remove ("display.html");
    
    // fixme: close connection cleanly
    close (Display.sock);
  }
}

@init_solver
void display_init()
{
  if (pid() == 0) {
    const char * port = DISPLAY_RANGE;
    if (!strchr (port, ':'))
      Display.sock = ws_socket_open (atoi (port));
    else {
      char * s = strdup (port);
      char * s1 = strchr (s, ':');
      *s1++ = '\0';
      int pmax = atoi(s1);
      Display.port = atoi(s);
      while ((Display.sock = ws_socket_open (Display.port)) < 0 &&
	     Display.port < pmax)
	Display.port++;
      free (s);
    }
    if (Display.sock < 0) {
      char s[80];
      sprintf (s, "display(): could not open port '%s'", port);
      perror (s);
      exit (1);
    }

    FILE * fp = fopen ("display.html", "w");
    fputs ("<head><meta http-equiv=\"refresh\" content=\"0;URL=", fp);
    display_url (fp);
    fputs ("\"></head>\n", fp);
    fclose (fp);
  }

  Display.objects = kh_init (strhash);
  Display.controls = array_new();
  
  free_solver_func_add (display_destroy);

#ifndef DISPLAY
  display_control (display_play, -1, 1, "Run/Pause");
#elif DISPLAY != 0
  display_control (display_play, -1, 1, "Run/Pause");
#endif
#ifndef DISPLAY_NO_CONTROLS    
  display_control (display_usage, 0, 50, "Display %", 
		   "maximum % of runtime used by display");
#endif
}

event refresh_display (i++, last)
{
  do {
    if (display_play)
      display_update (i);
    if (display_poll (display_play ? - 1 : 0))
      display_update (i);
  } while (display_play < 0);

  static timer global_timer = {0};
  static double poll_elapsed = 0.;
  int refresh = (poll_elapsed <=
		 display_usage/100.*timer_elapsed (global_timer));
@if _MPI
  MPI_Bcast (&refresh, 1, MPI_INT, 0, MPI_COMM_WORLD);
@endif
  if (refresh) {
    global_timer = timer_start();
    display_update (i);
    poll_elapsed = timer_elapsed (global_timer);
  }
}

#endif
