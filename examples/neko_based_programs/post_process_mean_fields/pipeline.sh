rm interpolated_fields.csv
cp ../../../../unstructured_meshing/coordinates.csv .
./post_process_mean_fields
cp interpolated_fields.csv ../../../../unstructured_meshing/
