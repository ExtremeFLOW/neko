#include <string.h>
#include <jansson.h>
#include <neko.h>

int main(int argc, char **argv) {

  /* Initialise Neko */
  neko_init();
  neko_job_info();

  /* Create a Neko case from a JSON file */
  json_error_t json_error;
  const char *case_json = json_dumps(json_load_file("cylinder.case", 0,
                                                    &json_error), JSON_COMPACT);
  int *neko_case = NULL;
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
