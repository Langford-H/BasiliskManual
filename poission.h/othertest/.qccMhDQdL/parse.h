#ifndef BASILISK_HEADER_8
#define BASILISK_HEADER_8
#line 1 "/home/dahuanghhc/basilisk/src/parse.h"
/**
# Parser for interactive bview
*/

enum ParamsType { pstring, pint, punsigned, pbool, pfloat, pdouble, pcolormap };

typedef struct {
  char * key;
  enum ParamsType type;
  void * val;
  int n;
} Params;

static bool atobool (char * s)
{
  if (!strcmp (s, "true"))
    return true;
  if (!strcmp (s, "false"))
    return false;
  return atoi (s) != 0;
}
  
static bool args (Params * p, char * val)
{
  static char * name[] = { "string", "int", "unsigned",
			   "bool", "float", "double", "colormap" };  
  switch (p->type) {

  case pstring:
    if (val[0] != '"') {
      fprintf (stderr, "expecting a string for '%s' got '%s'\n", p->key, val);
      return false;
    }
    if (val[strlen(val) - 1] != '"') {
      fprintf (stderr, "unterminated quoted string '%s'\n", val);
      return false;
    }
    val[strlen(val) - 1] = '\0';
    char * s = &val[1];
    int nc = 0; // number of non-blank characters
    while (*s != '\0') {
      if (!strchr (" \t\n\r", *s))
	nc++;
      s++;
    }
    *((char **)p->val) = nc > 0 ? &val[1] : NULL;
    break;

  case pcolormap:
    if (!strcmp (val, "jet"))
      *((colormap *)p->val) = jet;
    else if (!strcmp (val, "cool_warm"))
      *((colormap *)p->val) = cool_warm;
    else if (!strcmp (val, "gray"))
      *((colormap *)p->val) = gray;
    else if (!strcmp (val, "randomap"))
      *((colormap *)p->val) = randomap;
    else {
      fprintf (stderr, "unknown colormap '%s'\n", val);
      return false;
    }
    break;

  case pint: case punsigned: case pbool: case pdouble: case pfloat:
    if (val[0] == '"') {
      fprintf (stderr, "expecting a %s for '%s' got %s\n",
	       name[p->type], p->key, val);
      return false;
    }
    if (!p->n) {
      switch (p->type) {
      case pint: *((int *)p->val) = atoi(val); break;
      case punsigned: *((unsigned *)p->val) = atoi(val); break;
      case pbool: *((bool *)p->val) = atobool(val); break;
      case pfloat: *((float *)p->val) = atof(val); break;
      case pdouble: *((double *)p->val) = atof(val); break;
      default: assert (false);
      }
    }
    else {
      if (val[0] != '{') {
	fprintf (stderr, "expecting an array for '%s' got %s\n", p->key, val);
	return false;
      }
      val++;
      int i = 0;
      char c = ',';
      while (i < p->n && c != '}') {
	char * s = strchr (val, ',');
	if (!s)
	  s = strchr (val, '}');
	if (!s) {
	  fprintf (stderr, "expecting an array for '%s' got %s\n", p->key, val);
	  return false;
	}
	c = *s;
	*s++ = '\0';
	switch (p->type) {
	case pint: ((int *)p->val)[i++] = atoi (val); break;
	case punsigned: ((unsigned *)p->val)[i++] = atoi (val); break;
	case pbool: ((bool *)p->val)[i++] = atobool (val); break;
	case pfloat: ((float *)p->val)[i++] = atof (val); break;
	case pdouble: ((double *)p->val)[i++] = atof (val); break;
	default: assert (false);
	}
	val = s;
      }
      if (c != '}') {
	fprintf (stderr, "expecting '}' for '%s' got %s\n", p->key, val);
	return false;
      }
    }
    break;

  default:
    assert (false);
  }
  return true;
}

static char * find_comma (char * s)
{
  int par = 0;
  while (*s != '\0') {
    if (*s == ',' && par == 0) {
      *s = '\0';
      return s + 1;
    }
    if (*s == '{')
      par++;
    else if (*s == '}')
      par--;
    s++;
  }
  return NULL;
}

static char * mystrtok (char * str, const char * delim)
{
  static char * s = NULL;
  char * start = str ? str : s;
  bool string = false;
  s = start;
  while (*s != '\0') {
    if (*s == '"')
      string = !string;
    if (!string && strchr(delim, *s))
      break;
    s++;
  }
  if (*s != '\0')
    *s++ = '\0';
  return start;
}

int parse_params (Params * params)
{
  char * s;
  int i = 0, n = 0;
  Params * p = params;
  while (p->key) p++, n++;
  if (!(s = mystrtok (NULL, ");")) || s[0] == '\n')
    return false;
  while (s && *s != '\0') {
    char * next = find_comma (s), * key = s;
    if ((s = strchr (key, '='))) {
      s[0] = '\0', s++;
      i = -1;
      Params * p = params;
      while (p->key && strcmp(p->key, key)) p++;
      if (!p->key) {
	fprintf (stderr, "unknown key '%s'\n", key);
	return false;
      }
      if (!args (p, s))
	return false;
    }
    else {
      if (i < 0) {
	fprintf (stderr, "anonymous value '%s' after keys\n", key);
	return false;
      }
      if (i >= n) {
	fprintf (stderr, "too many parameters: '%s' %d %d\n", key, i, n);
	return false;
      }
      if (!args (&params[i], key))
	return false;
      i++;
    }
    s = next;
  }
  return true;
}

#endif
