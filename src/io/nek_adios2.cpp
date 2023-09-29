#include <adios2.h>
#include <string>
#include <iostream>
#include <ctime>

adios2::ADIOS adios;
adios2::IO io;
adios2::IO ior;
adios2::IO io_head;
adios2::IO io_asynchronous;
adios2::Engine writer;
adios2::Engine writer_head;
adios2::Engine readr; 
adios2::Engine writer_st;
adios2::Variable<double> x;
adios2::Variable<double> y;
adios2::Variable<double> z;
adios2::Variable<int> lglelw;
adios2::Variable<int> lglelw_st;
adios2::Variable<double> p;
adios2::Variable<double> vx;
adios2::Variable<double> vy;
adios2::Variable<double> vz;
adios2::Variable<double> bm1;
adios2::Variable<double> p_st;
adios2::Variable<double> vx_st;
adios2::Variable<double> vy_st;
adios2::Variable<double> vz_st;
adios2::Variable<double> bm1_st;
adios2::Variable<double> vxr;
adios2::Variable<double> vyr;
adios2::Variable<double> vzr;
adios2::Variable<double> prr;
adios2::Variable<int> lglelr;
std::vector<double> vVXr;
std::vector<double> vVYr;
std::vector<double> vVZr;
std::vector<double> vVPrr;
std::vector<int> vVLGLELr;
// adios2::Variable<double> t;
double dataTime = 0.0;
std::clock_t startT;
std::clock_t elapsedT;
int rank, size;
int ifile;
int ifilew;
int ifstream;
int decide_stream_global;


extern "C" void adios2_setup_(
    const int *nval,
    const int *nelvin,
    const int *nelb,
    const int *nelgv,
    const int *nelgt,
    const double *xml,
    const double *yml,
    const double *zml,
    const int *if_asynchronous,
    const int *comm_int
){
    std::string configFile="adios2_config/config.xml";
    MPI_Comm comm = MPI_Comm_f2c(*comm_int);
    adios = adios2::ADIOS(configFile, comm);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    // Compressor writer.
    io = adios.DeclareIO("writer");
    // Mesh writer.
    io_head = adios.DeclareIO("writer0");
    // Asynchronous writer.
    io_asynchronous = adios.DeclareIO("writerISMPI");
    // Compressor reader
    ior = adios.DeclareIO("inputReader");

    // Determine if asyncrhonous operation will be needed for this set up
    unsigned int decide_stream = static_cast<unsigned int>((*if_asynchronous));
    decide_stream_global = decide_stream;
    
    // Number of elements in my rank.
    unsigned int nelv = static_cast<unsigned int>((*nelvin));

    // Determine where my rank writes in the global array according to number of element in previous ranks
    unsigned int start = static_cast<unsigned int>(*nelb);
    start *= static_cast<unsigned int>(*nval);

    // n is count, i.e number of entries in the array in my rank
    unsigned int n = static_cast<unsigned int> (*nval) * nelv;
    // gn is the total size of the arrays, not per io rank 
    unsigned int gn = static_cast<unsigned int>((*nelgv)*(*nval));
    std::cout << rank << ": " << gn << ", " << start << "," << n << std::endl;

    // Create the adios2 variables for writer that depend on the current start and n
    p = io.DefineVariable<double>("P_OUT", {gn}, {start}, {n});
    vx = io.DefineVariable<double>("VX_OUT", {gn}, {start}, {n});
    vy = io.DefineVariable<double>("VY_OUT", {gn}, {start}, {n});
    vz = io.DefineVariable<double>("VZ_OUT", {gn}, {start}, {n});
    bm1 = io.DefineVariable<double>("BM1_OUT", {gn}, {start}, {n});

    // Create the adios2 variables for writer0
    x = io_head.DefineVariable<double>("X", {gn}, {start}, {n});
    y = io_head.DefineVariable<double>("Y", {gn}, {start}, {n});
    z = io_head.DefineVariable<double>("Z", {gn}, {start}, {n});

    // If the process is asynchronous, define the relevant variables for writer_st
    if (decide_stream == 1){
	    p_st = io_asynchronous.DefineVariable<double>("P", {gn}, {start}, {n});
	    vx_st = io_asynchronous.DefineVariable<double>("VX", {gn}, {start}, {n});
	    vy_st = io_asynchronous.DefineVariable<double>("VY", {gn}, {start}, {n});
	    vz_st = io_asynchronous.DefineVariable<double>("VZ", {gn}, {start}, {n});
	    bm1_st = io_asynchronous.DefineVariable<double>("BM1", {gn}, {start}, {n});
    }


    // Do everything again for the global indices 
    nelv = static_cast<unsigned int>((*nelvin));
    start = static_cast<unsigned int>(*nelb);
    n = static_cast<unsigned int> (nelv);
    gn = static_cast<unsigned int>((*nelgv));
    // Define variable for compression writer
    lglelw = io.DefineVariable<int>("LGLEL_OUT", {gn}, {start}, {n});
    // Define variable for asyncrhonous writet
    if (decide_stream == 1){
    	lglelw_st = io_asynchronous.DefineVariable<int>("LGLEL", {gn}, {start}, {n});
    }

    // Write the mesh information only once (Currently commented out).
    //writer_head = io_head.Open("geo.bp", adios2::Mode::Write);
    //writer_head.Put<double>(x, xml);
    //writer_head.Put<double>(y, yml);
    //writer_head.Put<double>(z, yml);
    //writer_head.Close();
    //if(!rank)
	//std::cout << "geo.bp written" << std::endl;
    
    // If asyncrhonous execution, open the global array
    if (decide_stream == 1){
	std::cout << "create global array" << std::endl;
    	writer_st = io_asynchronous.Open("globalArray", adios2::Mode::Write);
    }

    // Initialize global variables for writing. This could be done in global definition    
    ifile = 0 ;
    ifilew = 0 ;
}

extern "C" void adios2_update_(
    const int *lglel,
    const double *pr,
    const double *u,
    const double *v,
    const double *w,
    const double *temp
){
    startT = std::clock();
    ifilew=ifilew+1;
    std:: string fileName= "out.f0000"+ std::to_string(ifilew) +".bp";
    writer = io.Open(fileName, adios2::Mode::Write);
    // Begin a step of the writer
    writer.BeginStep();
    writer.Put<double>(p, pr);
    writer.Put<double>(vx, u);
    writer.Put<double>(vy, v);
    writer.Put<double>(vz, w);
    writer.Put<int>(lglelw, lglel);
    writer.EndStep();
    writer.Close();
    dataTime += (std::clock() - startT) / (double) CLOCKS_PER_SEC;
}

extern "C" void adios2_stream_(
    const int *lglel,
    const double *pr,
    const double *u,
    const double *v,
    const double *w,
    const double *mass1,
    const double *temp
){
    startT = std::clock();
    // Begin a step of the writer
    writer_st.BeginStep();
    writer_st.Put<double>(p_st, pr);
    writer_st.Put<double>(vx_st, u);
    writer_st.Put<double>(vy_st, v);
    writer_st.Put<double>(vz_st, w);
    writer_st.Put<double>(bm1_st, mass1);
    writer_st.Put<int>(lglelw_st, lglel);
    writer_st.EndStep();
    dataTime += (std::clock() - startT) / (double) CLOCKS_PER_SEC;
}
extern "C" void adios2_read_(
    int *lglelrr,
    double *pr,
    double *v,
    double *u,
    double *w,
    const int *nval,
    const int *nelvin,
    const int *nelb,
    const int *nelgv,
    const int *nelgt,
    const int *comm_int,
    char *fname
){
    startT = std::clock();
    std::string configFile="adios2_config/config.xml";
    MPI_Comm comm = MPI_Comm_f2c(*comm_int);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    // See how much of the file my rank needs to read
    unsigned int nelv = static_cast<unsigned int>((*nelvin));
    unsigned int start = static_cast<unsigned int>(*nelb);
    start *= static_cast<unsigned int>(*nval);
    unsigned int n = static_cast<unsigned int> (*nval) * nelv;
    unsigned int gn = static_cast<unsigned int>((*nelgv)*(*nval));

    // See how much of the file my rank needs to read for the global indices
    unsigned int nelv3 = static_cast<unsigned int>((*nelvin));
    unsigned int start3 = static_cast<unsigned int>(*nelb);
    unsigned int n3 = static_cast<unsigned int> (nelv);
    unsigned int gn3 = static_cast<unsigned int>((*nelgv));
   
    // Write some debuffing info 
    if (rank==0){
    std::cout << "what adios is getting:"<< std::endl;
    std::cout << "nvals:"<< static_cast<unsigned int> (*nval) << std::endl;
    std::cout << "nelgv:"<< static_cast<unsigned int> (*nelgv) << std::endl;
    std::cout << "nelgt:"<< static_cast<unsigned int> (*nelgt) << std::endl;

    std::cout << "what adios is calculating:"<< std::endl;
    std::cout << "nelv:"<< nelv << std::endl;
    std::cout << "start:"<< start << std::endl;
    std::cout << "count:"<< n << std::endl;
    std::cout << "total number of entries:"<< gn << std::endl;
    }
    

    // Read the file with factory name
    ifile=ifile+1;
    // Make a new name every time the reader is called
    std:: string fileName= "out.f0000"+ std::to_string(ifile) +".bp";
    startT = std::clock();
    readr = ior.Open(fileName,adios2::Mode::Read); 

    int step = 1;
    bool firstStep = true;
    std::cout << "CAREFULL, YOU ARE SAYING EVERY STEP IS THE FIRST, SO IT ALLOCATES EVERY TIME" << std::endl;

    // Check if the variables that we are requesting are in the file
        vxr = ior.InquireVariable<double>("VX_OUT");
        if (!vxr)
        {
            std::cout << "Error: NO variable VX_OUT found. Unable to proceed. " << std::endl;
        }

        vyr = ior.InquireVariable<double>("VY_OUT");
        if (!vyr)
        {
            std::cout << "Error: NO variable VY_OUT found. Unable to proceed. " << std::endl;
        }

        vzr = ior.InquireVariable<double>("VZ_OUT");
        if (!vzr)
        {
            std::cout << "Error: NO variable VZ_OUT found. Unable to proceed. " << std::endl;
        }

        prr = ior.InquireVariable<double>("P_OUT");
        if (!prr)
        {
            std::cout << "Error: NO variable P_OUT found. Unable to proceed. " << std::endl;
        }
        lglelr = ior.InquireVariable<int>("LGLEL_OUT");
        if (!lglelr)
        {
            std::cout << "Error: NO variable LGLEL_OUT found. Unable to proceed. " << std::endl;
        }

    elapsedT = (std::clock() - startT) / (double) CLOCKS_PER_SEC;
    std::cout <<  "rank: " << rank << " open and inquire: " << elapsedT << "s." << std::endl;
    startT = std::clock();
        // Allocate variables if it is the first step
        if (firstStep)
        {
                readr.LockReaderSelections();
		vVXr.resize(n);
		vVYr.resize(n);
		vVZr.resize(n);
		vVPrr.resize(n);
		vVLGLELr.resize(n3);
		firstStep=false;
        }

    elapsedT = (std::clock() - startT) / (double) CLOCKS_PER_SEC;
    std::cout <<  "rank: " << rank << " resize vector: " << elapsedT << "s." << std::endl;
    startT = std::clock();
 
        // Select in the adios file the sections to read
    	vxr.SetSelection({{start}, {n}});
        vyr.SetSelection({{start}, {n}});
        vzr.SetSelection({{start}, {n}});
        prr.SetSelection({{start}, {n}});
        lglelr.SetSelection({{start3}, {n3}});

    elapsedT = (std::clock() - startT) / (double) CLOCKS_PER_SEC;
    std::cout <<  "rank: " << rank << " set selection: " << elapsedT << "s." << std::endl;
    startT = std::clock();
 
        // get the data that has been selected
        readr.Get<double>(vxr,vVXr.data());
        readr.Get<double>(vyr,vVYr.data());
        readr.Get<double>(vzr,vVZr.data());
        readr.Get<double>(prr,vVPrr.data());
        readr.Get<int>(lglelr,vVLGLELr.data());
    
    elapsedT = (std::clock() - startT) / (double) CLOCKS_PER_SEC;
    std::cout <<  "rank: " << rank << " get data: " << elapsedT << "s." << std::endl;
    startT = std::clock();

    // Close the reader
    readr.Close();	
    
    elapsedT = (std::clock() - startT) / (double) CLOCKS_PER_SEC;
    std::cout <<  "rank: " << rank << " Close: " << elapsedT << "s." << std::endl;
    startT = std::clock();

    //copy the data from the buffer to the array that is used in nek/neko
    if (rank==0){
    std::cout <<"Copying data into nek vector " << std::endl;
    }
    
    for (int i=0; i<vVZr.size(); i++){
    pr[i]=vVPrr[i];
    v[i]=vVXr[i];
    u[i]=vVYr[i];
    w[i]=vVZr[i];
    }
    for (int i=0; i<vVLGLELr.size(); i++){
    lglelrr[i]=vVLGLELr[i];
    }
    
    elapsedT = (std::clock() - startT) / (double) CLOCKS_PER_SEC;
    std::cout <<  "rank: " << rank << " copy into nek: " << elapsedT << "s." << std::endl;
    startT = std::clock();

    dataTime += (std::clock() - startT) / (double) CLOCKS_PER_SEC;

}

extern "C" void adios2_finalize_(){
    if (decide_stream_global == 1){
	std::cout << "Close global array" << std::endl;
    	writer_st.Close();
    	std::cout <<  "rank: " << rank << " in-situ time: " << dataTime << "s." << std::endl;
    }

}
