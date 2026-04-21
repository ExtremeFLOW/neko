#include <adios2.h>
#include <mpi.h>
#include <thread>
#include <chrono>
#include <string>
#include <iostream>
#include <ctime>

// Global Adios2 variables
adios2::ADIOS adios;
adios2::IO io_reader;
adios2::IO io_writer;
adios2::Engine writer_st;
adios2::Engine reader_st;
adios2::Variable<double> f2py_field;
adios2::Variable<double> py2f_field;
// Global C variables
int rank, size;
int sync_rank = -1;
int sync_size = -1;
unsigned int reader_start;
unsigned int reader_count;

static void dbg_log(const std::string &msg)
{
    std::cerr
        << "[NEKO-ADIOS]"
        << "[rank=" << rank << "/" << size << "]";
    if (sync_rank >= 0 && sync_size > 0)
    {
        std::cerr
            << "[sync_rank=" << sync_rank << "/" << sync_size << "]";
    }
    std::cerr << " " << msg << std::endl;
}

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
    const int *timeout_seconds
){
    MPI_Comm comm = MPI_Comm_f2c(*comm_int);
    MPI_Comm sync_comm = MPI_Comm_f2c(*sync_comm_int);
    adios = adios2::ADIOS(comm);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(sync_comm, &sync_rank);
    MPI_Comm_size(sync_comm, &sync_size);
    dbg_log(
        "initialize start"
        " timeout=" + std::to_string(*timeout_seconds) +
        " lxyz=" + std::to_string(*lxyz) +
        " nelv=" + std::to_string(*nelv) +
        " offset_el=" + std::to_string(*offset_el) +
        " glb_nelv=" + std::to_string(*glb_nelv) +
        " gdim=" + std::to_string(*gdim)
    );
    // ADIOS2 IOs.
    io_reader = adios.DeclareIO("streamReaderIO");
    io_reader.SetEngine("SST");
    io_reader.SetParameters(
        {{"OpenTimeoutSecs", std::to_string(*timeout_seconds)}}
    );

    io_writer = adios.DeclareIO("streamWriterIO");
    io_writer.SetEngine("SST");
    io_writer.SetParameters(
        {{"OpenTimeoutSecs", std::to_string(*timeout_seconds)}}
    );
    dbg_log("declared streamReaderIO and streamWriterIO");

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
    dbg_log(
        "local layout"
        " start=" + std::to_string(reader_start) +
        " count=" + std::to_string(reader_count) +
        " total=" + std::to_string(gn)
    );

    // If the process is asynchronous, define the relevant variables for writer_st
    f2py_field = io_writer.DefineVariable<double>(
        "f2py_field", {gn}, {start}, {n}
    );
    dbg_log("defined f2py_field");
    
    dbg_log("barrier before writer open");
    MPI_Barrier(sync_comm);
    dbg_log("barrier before writer open complete");
    init_wait();

    dbg_log("opening writer globalArray_f2py");
    writer_st = io_writer.Open("globalArray_f2py", adios2::Mode::Write);
    dbg_log("writer globalArray_f2py open complete");
    init_wait();

    dbg_log("barrier before reader open");
    MPI_Barrier(sync_comm);
    dbg_log("barrier before reader open complete");
    init_wait();
    dbg_log("opening reader globalArray_py2f");
    reader_st = io_reader.Open("globalArray_py2f", adios2::Mode::Read);
    dbg_log("reader globalArray_py2f open complete");
    init_wait();

    dbg_log("barrier after stream opens");
    MPI_Barrier(sync_comm);
    dbg_log("barrier after stream opens complete");
    init_wait();

    // Put necesary information in a header stream
    dbg_log("writer BeginStep for header");
    writer_st.BeginStep();
    adios2::Variable<int> hdr_elems =
        io_writer.DefineVariable<int>("global_elements");
    adios2::Variable<int> hdr_lxyz =
        io_writer.DefineVariable<int>("points_per_element");
    adios2::Variable<int> hdr_gdim =
        io_writer.DefineVariable<int>("problem_dimension");
    if( rank == 0 )
    {
       writer_st.Put(hdr_elems, static_cast<int> (*glb_nelv));
       writer_st.Put(hdr_lxyz,  static_cast<int> (*lxyz));
       writer_st.Put(hdr_gdim,  static_cast<int> (*gdim));
    }
    writer_st.EndStep();
    dbg_log("writer EndStep complete for header");
    init_wait();

    dbg_log("barrier after header write");
    MPI_Barrier(sync_comm);
    dbg_log("barrier after header write complete");
    init_wait();
    dbg_log("final init barrier");
    MPI_Barrier(sync_comm);
    dbg_log("initialize complete");
}

extern "C" void adios2_finalize_(){
    dbg_log("closing global arrays");
    writer_st.Close();
    reader_st.Close();
    dbg_log("global arrays closed");

}

extern "C" void adios2_stream_(
    const double *field
){
    dbg_log(
        "writer BeginStep for f2py_field"
        " count=" + std::to_string(reader_count)
    );
    writer_st.BeginStep();
    writer_st.Put<double>(f2py_field, field);
    writer_st.EndStep();
    dbg_log("writer EndStep complete for f2py_field");
}

extern "C" void adios2_recieve_(
    double *field
){
    dbg_log(
        "reader BeginStep for py2f_field"
        " start=" + std::to_string(reader_start) +
        " count=" + std::to_string(reader_count)
    );
    reader_st.BeginStep();
    py2f_field = io_reader.InquireVariable<double>("py2f_field");
    dbg_log("InquireVariable(py2f_field) completed");
    py2f_field.SetSelection({{reader_start}, {reader_count}});
    reader_st.Get<double>(py2f_field, field);
    reader_st.EndStep();
    dbg_log("reader EndStep complete for py2f_field");
}
