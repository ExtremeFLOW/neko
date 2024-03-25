#include <adios2.h>
#include <string>
#include <iostream>
#include <ctime>

// Global Adios2 variables
adios2::ADIOS adios;
adios2::IO io_asynchronous;
adios2::Engine writer_st;
adios2::Engine reader_st;
adios2::Variable<double> f2py_field;
adios2::Variable<double> py2f_field;
// Global C variables
int rank, size;
int reader_start;
int reader_count;

extern "C" void adios2_initialize_(
    const int *lxyz,
    const int *nelv,
    const int *offset_el,
    const int *glb_nelv,
    const int *gdim,
    const int *comm_int
){
    MPI_Comm comm = MPI_Comm_f2c(*comm_int);
    adios = adios2::ADIOS(comm);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    // Asynchronous IO.
    io_asynchronous = adios.DeclareIO("streamIO");
    io_asynchronous.SetEngine("SST");

    // Number of elements in my rank.
    unsigned int nel = static_cast<unsigned int>((*nelv));
    // Determine where my rank writes in the global array according to number of element in previous ranks
    unsigned int start = static_cast<unsigned int>(*offset_el);
    start *= static_cast<unsigned int>(*lxyz);
    // n is count, i.e number of entries in the array in my rank
    unsigned int n = static_cast<unsigned int> (*lxyz) * nel;
    // gn is the total size of the arrays, not per io rank 
    unsigned int gn = static_cast<unsigned int>((*glb_nelv)*(*lxyz));
 
    // Assign to global variables
    reader_start = start;
    reader_count = n;

    // If the process is asynchronous, define the relevant variables for writer_st
    f2py_field = io_asynchronous.DefineVariable<double>("f2py_field", {gn}, {start}, {n});
    
    // If asyncrhonous execution, open the global array
    std::cout << "create global array" << std::endl;
    writer_st = io_asynchronous.Open("globalArray_f2py", adios2::Mode::Write);
    reader_st = io_asynchronous.Open("globalArray_py2f", adios2::Mode::Read);

    // Put necesary information in a header stream
    writer_st.BeginStep();
    adios2::Variable<int> hdr_elems = io_asynchronous.DefineVariable<int>("global_elements");
    adios2::Variable<int> hdr_lxyz = io_asynchronous.DefineVariable<int>("points_per_element");
    adios2::Variable<int> hdr_gdim = io_asynchronous.DefineVariable<int>("problem_dimension");
    if( rank == 0 )
    {
       writer_st.Put(hdr_elems, static_cast<int> (*glb_nelv));
       writer_st.Put(hdr_lxyz,  static_cast<int> (*lxyz));
       writer_st.Put(hdr_gdim,  static_cast<int> (*gdim));
    }
    writer_st.EndStep();
}

extern "C" void adios2_finalize_(){
    std::cout << "Close global arrays" << std::endl;
    writer_st.Close();
    reader_st.Close();

}

extern "C" void adios2_stream_(
    const double *field
){
    writer_st.BeginStep();
    writer_st.Put<double>(f2py_field, field);
    writer_st.EndStep();
}

extern "C" void adios2_recieve_(
    double *field
){
    reader_st.BeginStep();
    py2f_field = io_asynchronous.InquireVariable<double>("py2f_field");
    py2f_field.SetSelection({{reader_start}, {reader_count}});
    reader_st.Get<double>(py2f_field, field);
    reader_st.EndStep();
}
