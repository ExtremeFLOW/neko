rm interpolated_fields.csv
cp ../../../../unstructured_meshing/post_process_timedependent_fields/coordinates.csv .
./post_process_timedependent_fields
cp interpolated_fields.csv ../../../../unstructured_meshing/post_process_timedependent_fields/
