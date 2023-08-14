# Make sure nek5000 tools are installed.

cp placeholderrea cylinder.rea

nekpath="/home/adalberto/Nek5000/"
tools=$nekpath"tools/"
bin=$nekpath"bin/"

nekopath="/home/adalberto/neko/"
contrib=$nekopath"contrib/rea2nbin/"

export PATH=$tools:$PATH
export PATH=$bin:$PATH
export PATH=$contrib:$PATH

# Generate file from gmsh
gmsh pipeMesh.geo -2 -order 2

# Convert to re2
gmsh2nek < gmsh2nek_input

# Convert to rea
re2torea < re2torea_input

# Extrude to 3d
n2to3 < n2to3_input

reatore2 < reatore2_input

# Convert to neko mesh
rea2nbin cylinder_3d.re2 cylinder.nmsh
