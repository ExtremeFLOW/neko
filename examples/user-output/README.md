## Custonized user output to a csv file
In this case we simulate a lid-driven cavity with smoothened belt velocity to fulfill the continuity equation in the corners. This case is ran as a pseudo-2D example (with a 2D mesh that then generates a case with one element in the spanwise direction).

The user file implements a customized user output in the user\_check function using the csv\_file\_t type.

The cavity flow is stable (steady) up to a Reynolds number of about 8000. The mesh size may be changed in the box file.
