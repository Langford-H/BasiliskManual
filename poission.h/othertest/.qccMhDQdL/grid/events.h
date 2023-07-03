#ifndef BASILISK_HEADER_28
#define BASILISK_HEADER_28
#line 1 "/home/dahuanghhc/basilisk/src/grid/events.h"
#define INIT ev->expr[0]
#define COND ev->expr[1]
#define INC  ev->expr[2]
#define END_EVENT 1234567890

#define TEPS 1e-9

static void event_error (Event * ev, const char * s)
{
  fprintf (stderr, "%s:%d: error: %s\n", ev->file, ev->line, s);
  exit (1);
}

static void init_event (Event * ev)
{
  if (ev->arrayi || ev->arrayt) {
    ev->i = ev->t = -1;
    if (ev->arrayi)
      ev->i = ev->arrayi[0];
    else 
      ev->t = ev->arrayt[0];
    ev->a = 1;
    ev->expr[1] = NULL;
  }
  else {
    if (ev->nexpr > 0) {
      Expr init = NULL, cond = NULL, inc = NULL;
      for (int j = 0; j < ev->nexpr; j++) {
	int i = -123456; double t = i;
	(* ev->expr[j]) (&i, &t, ev);
	if (i == -123456 && t == -123456) {
	  /* nothing done to i and t: this must be the condition */
	  if (cond)
	    event_error (ev, "events can only use a single condition");
	  cond = ev->expr[j];
	}
	else {
	  /* this is either an initialisation or an increment */
	  int i1 = i; double t1 = t;
	  (* ev->expr[j]) (&i1, &t1, ev);
	  if (i1 == i && t1 == t) {
	    /* applying twice does not change anything: this is an
	       initialisation */
	    if (init)
	      event_error (ev, "events can only use a single initialisation");
	    init = ev->expr[j];
	  }
	  else {
	    /* this is the increment */
	    if (inc)
	      event_error (ev, "events can only use a single increment");
	    inc = ev->expr[j];
	  }
	}
      }
      INIT = init;
      COND = cond;
      INC  = inc;
      ev->nexpr = 0;
    }
    ev->i = ev->t = -1;
    if (INIT) {
      (* INIT) (&ev->i, &ev->t, ev);
      if (ev->i == END_EVENT || ev->t == END_EVENT) {
	ev->i = END_EVENT; ev->t = -1;
      }
    }
    else if (INC) {
      (* INC) (&ev->i, &ev->t, ev);
      if (ev->i != -1)
	ev->i = 0;
      if (ev->t != -1)
	ev->t = 0;
    }
  }
}

enum { event_done, event_alive, event_stop };

static int event_finished (Event * ev)
{
  ev->t = ev->i = -1;
  return event_done;
}

void event_register (Event event) {
  assert (Events);
  assert (!event.last);
  int n = 0, parent = -1;
  for (Event * ev = Events; !ev->last; ev++) {
    if (!strcmp (event.name, ev->name)) {
      assert (parent < 0);
      parent = n;
    }
    n++;
  }
  if (parent < 0) {
    qrealloc (Events, n + 2, Event);
    Events[n] = event;
    Events[n].next = NULL;
    Events[n + 1].last = true;
    init_event (&Events[n]);
  }
  else {
    Event * ev = qcalloc (1, Event);
    *ev = Events[parent];
    Events[parent] = event;
    Events[parent].next = ev;
    init_event (&Events[parent]);
  }
}

static int event_cond (Event * ev, int i, double t)
{
  if (!COND)
    return true;
  return (* COND) (&i, &t, ev);
}

#if DEBUG_EVENTS
static void event_print (Event * ev, FILE * fp)
{
  char * root = strstr (ev->file, BASILISK);
  fprintf (fp, "  %-25s %s%s:%d\n", ev->name, 
	   root ? "src" : "",
	   root ? &ev->file[strlen(BASILISK)] : ev->file, 
	   ev->line);
}
#endif

static int event_do (Event * ev, bool action)
{
  if ((iter > ev->i && t > ev->t) || !event_cond (ev, iter, t))
    return event_finished (ev);
  if (iter == ev->i || fabs (t - ev->t) <= TEPS) {
    if (action) {
      bool finished = false;
      for (Event * e = ev; e; e = e->next) {
#if DEBUG_EVENTS
	event_print (e, stderr);
#endif
	if ((* e->action) (iter, t, e))
	  finished = true;
      }
      if (finished) {
	event_finished (ev);
	return event_stop;
      }
    }
    if (ev->arrayi) { /* i = {...} */
      ev->i = ev->arrayi[ev->a++];
      if (ev->i < 0)
	return event_finished (ev);
    }
    if (ev->arrayt) { /* t = {...} */
      ev->t = ev->arrayt[ev->a++];
      if (ev->t < 0)
	return event_finished (ev);
    }
    else if (INC) {
      int i0 = ev->i;
      (* INC) (&ev->i, &ev->t, ev);
      if (i0 == -1 && ev->i != i0)
	ev->i += iter + 1;
      if (!event_cond (ev, iter + 1, ev->t))
	return event_finished (ev);
    }
    else if (INIT && !COND)
      return event_finished (ev);
  }
  return event_alive;
}

static void end_event_do (bool action)
{
#if DEBUG_EVENTS
  if (action)
    fprintf (stderr, "\nend events (i = %d, t = %g)\n", iter, t);
#endif
  for (Event * ev = Events; !ev->last; ev++)
    if (ev->i == END_EVENT && action)
      for (Event * e = ev; e; e = e->next) {
#if DEBUG_EVENTS
	event_print (e, stderr);
#endif
	e->action (iter, t, e);
      }
}

int events (bool action)
{
#if DEBUG_EVENTS
  if (action)
    fprintf (stderr, "\nevents (i = %d, t = %g)\n", iter, t);
#endif

  if (iter == 0)
    for (Event * ev = Events; !ev->last; ev++)    
      init_event (ev);

  int cond = 0, cond1 = 0;
  inext = END_EVENT; tnext = HUGE;
  for (Event * ev = Events; !ev->last && !cond; ev++)
    if (ev->i != END_EVENT && 
	(COND || (INIT && !COND && !INC) || ev->arrayi || ev->arrayt))
      cond = 1;
  for (Event * ev = Events; !ev->last; ev++) {
    int status = event_do (ev, action);
    if (status == event_stop) {
      end_event_do (action);
      return 0;
    }
    if (status == event_alive && ev->i != END_EVENT &&
	(COND || (INIT && !COND && !INC) || ev->arrayi || ev->arrayt))
      cond1 = 1;
    if (ev->t > t && ev->t < tnext)
      tnext = ev->t;
    if (ev->i > iter && ev->i < inext)
      inext = ev->i;
  }
  if ((!cond || cond1) && (tnext != HUGE || inext != END_EVENT)) {
    inext = iter + 1;
    return 1;
  }
  end_event_do (action);
  return 0;
}

void event (const char * name)
{
  for (Event * ev = Events; !ev->last; ev++)
    if (!strcmp (ev->name, name))
      for (Event * e = ev; e; e = e->next) {
#if DEBUG_EVENTS
	event_print (e, stderr);
#endif     
	(* e->action) (0, 0, e);
      }
}

double dtnext (double dt)
{
  if (tnext != HUGE && tnext > t) {
    unsigned int n = (tnext - t)/dt;
    assert (n < INT_MAX); // check that dt is not too small
    if (n == 0)
      dt = tnext - t;
    else {
      double dt1 = (tnext - t)/n;
      if (dt1 > dt + TEPS)
	dt = (tnext - t)/(n + 1);
      else if (dt1 < dt)
	dt = dt1;
      tnext = t + dt;
    }
  }
  else
    tnext = t + dt;
  return dt;
}

#endif
