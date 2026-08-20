#include <string.h>
#include <jansson.h>
#include <neko.h>

/*
 * Helper to compare a callback's scheme name, which is not null-terminated
 */
static int scheme_name_is(const char *scheme_name, int scheme_name_len,
                          const char *name) {
  return ((scheme_name_len == (int) strlen(name)) &&
          (strncmp(scheme_name, name, scheme_name_len) == 0));
}

/* Define initial conditions */
static void initial(const char *scheme_name, int scheme_name_len) {

  /*
   * In this case the check is redundant, but if both fluid and scalar
   * initial conditions should be defined, `scheme_name` can be used
   * to identify which solver has called the callback
   */
  if (scheme_name_is(scheme_name, scheme_name_len, "fluid")) {

    /*
     * To retrive a pointer to the data in a field from Neko's field
     * registry use neko_field()
     */
    neko_real *u = neko_field("u");
    int n = neko_field_size("u");

    for (int i = 0; i < n; i++) {
      u[i] = 1.0;
    }
  }
}

/* Define inflow conditions */
static void inflow(int *msk, int msk_size, neko_real t, int tstep) {
  int idx = 1; /* Note: Fortran indices */

  /*
   * In this case this check is redundant, but if different user provided
   * boundary conditions share the same callback, one could use the content
   * of the callback's field list to identify which condition should be
   * applied, for example a velocity condition passes the fields (u, v, w)
   */
  if (neko_cb_field_name_at_index(&idx, "u")) {

    /*
     * Inside a callback, the fields passed to it are retrived with
     * neko_cb_field_by_name() rather than from the field registry
     */
    neko_real *u = neko_cb_field_by_name("u");
    neko_real *v = neko_cb_field_by_name("v");
    neko_real *w = neko_cb_field_by_name("w");

    for (int i = 0; i < msk_size; i++) {
      int j = msk[i] - 1; /* Neko's mask is based on Fortran indices */
      u[j] = 1.0;
      v[j] = 0.0;
      w[j] = 0.0;
    }
  }
}

int main(int argc, char **argv) {

  /* Initialise Neko */
  neko_init();
  neko_job_info();

  /* Create a Neko case from a JSON file */
  json_error_t json_error;
  const char *case_json = json_dumps(json_load_file("cylinder.case", 0,
                                                    &json_error), JSON_COMPACT);
  int *neko_case = NULL;

  /*
   * To provide user callbacks, the case must be allocated before it is
   * initialised, such that the callbacks are registered before the case
   * resolves the "user" entries in the JSON object. Unused callbacks are
   * left undefined by passing a null pointer.
   */
  neko_case_allocate(&neko_case);
  neko_user_setup(&neko_case, initial, NULL, NULL, inflow, NULL, NULL);

  neko_case_init(&case_json, strlen(case_json), &neko_case);

  /* To solve the entire case we can call neko_solve() */
  /* neko_solve(&neko_case); */

  /* To manually step forward in time, call neko_step() */
  while (neko_case_time(&neko_case) < neko_case_end_time(&neko_case)) {
    neko_step(&neko_case);
  }

  /*
   * To retrive a pointer to the data in a field from Neko's field
   * registry use neko_field()
   */
  neko_real *u = neko_field("u");

  /* Cleanup */
  neko_case_free(&neko_case);
  neko_finalize();
}
