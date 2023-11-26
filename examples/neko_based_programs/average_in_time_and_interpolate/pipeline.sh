rm interpolated_fields.csv
cp ../../../../unstructured_meshing/coordinates.csv .
./average_in_time_and_interpolate
cp interpolated_fields.csv ../../../../unstructured_meshing/
