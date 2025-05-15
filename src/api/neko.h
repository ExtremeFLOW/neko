#ifndef __NEKO_H

/* Initialise Neko */
void neko_init();

/* Finalize Neko */
void neko_finalize();

/* Display job information */
void neko_job_info();

/* Initialise a Neko case */
void neko_case_init(const char **case_json, int case_len, int *case_iptr);

/* Destroy a Neko case */
void neko_case_free(int *case_iptr);

/* Retrive the current time of a case */
double neko_case_time(int *case_iptr);

/* Retrive the end time of a case */
double neko_case_end_time(int *case_iptr);

/* Retrive the current time-step of a case */
int neko_case_tstep(int *case_iptr);

/* Solve a Neko case */
void neko_solve(int *case_iptr);

/* Compute a time-step for a neko case */
void neko_step(int *case_iptr);

#endif
