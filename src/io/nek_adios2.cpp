#include <adios2.h>
#include <mpi.h>
#include <thread>
#include <chrono>
#include <string>
#include <ctime>

// Global Adios2 variables
adios2::ADIOS adios;
adios2::IO io_stream;
adios2::Engine writer_st;
adios2::Engine reader_st;
adios2::Variable<double> f2py_field;
adios2::Variable<double> py2f_field;
unsigned int reader_start;
unsigned int reader_count;

static void init_wait()
{
    std::this_thread::sleep_for(std::chrono::seconds(2));
}

extern "C" void adios2_initialize_(
    const int *lxyz,
    const int *nelv,
    const int *offset_el,
    const int *glb_nelv,
    const int *gdim,
    const int *comm_int,
    const int *sync_comm_int,
    const int timeout_seconds
){
    int rank;
    MPI_Comm comm = MPI_Comm_f2c(*comm_int);
    MPI_Comm sync_comm = MPI_Comm_f2c(*sync_comm_int);
    adios = adios2::ADIOS(comm);
    MPI_Comm_rank(comm, &rank);
    // ADIOS2 IO.
    io_stream = adios.DeclareIO("streamIO");
    io_stream.SetEngine("SST");
    io_stream.SetParameters(
        {
	  {"OpenTimeoutSecs", std::to_string(timeout_seconds)}, 
	  {"RendezvousReaderCount", "1"},
	  {"QueueLimit", "1"},
	  {"QueueFullPolicy", "Block"}
	}
    );

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
    f2py_field = io_stream.DefineVariable<double>("f2py_field", {gn}, {start}, {n});

    //MPI_Barrier(sync_comm);
    //init_wait();
    writer_st = io_stream.Open("globalArray_f2py", adios2::Mode::Write);
    //init_wait();

    //MPI_Barrier(sync_comm);
    //init_wait();
    reader_st = io_stream.Open("globalArray_py2f", adios2::Mode::Read);
    //init_wait();

    //MPI_Barrier(sync_comm);
    //init_wait();

    // Put necesary information in a header stream
    writer_st.BeginStep();
    adios2::Variable<int> hdr_elems =
        io_stream.DefineVariable<int>("global_elements");
    adios2::Variable<int> hdr_lxyz =
        io_stream.DefineVariable<int>("points_per_element");
    adios2::Variable<int> hdr_gdim =
        io_stream.DefineVariable<int>("problem_dimension");
    if( rank == 0 )
    {
       writer_st.Put(hdr_elems, static_cast<int> (*glb_nelv));
       writer_st.Put(hdr_lxyz,  static_cast<int> (*lxyz));
       writer_st.Put(hdr_gdim,  static_cast<int> (*gdim));
    }
    writer_st.EndStep();
    //init_wait();

    //MPI_Barrier(sync_comm);
    //init_wait();
    //MPI_Barrier(sync_comm);
}

extern "C" void adios2_finalize_(){
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
    py2f_field = io_stream.InquireVariable<double>("py2f_field");
    py2f_field.SetSelection({{reader_start}, {reader_count}});
    reader_st.Get<double>(py2f_field, field);
    reader_st.EndStep();
}
