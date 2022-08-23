#include <adios2.h>
#include <string>
#include <iostream>
#include <ctime>

adios2::ADIOS adios;
adios2::IO io;
adios2::IO ior;
adios2::IO io_head;
adios2::Engine writer;
adios2::Engine writer_head;
adios2::Engine readr; 
adios2::Variable<double> x;
adios2::Variable<double> y;
adios2::Variable<double> z;
adios2::Variable<int> lglelw;
adios2::Variable<double> p;
adios2::Variable<double> vx;
adios2::Variable<double> vy;
adios2::Variable<double> vz;
adios2::Variable<double> bm1;
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


extern "C" void adios2_setup_(
    const int *nval,
    const int *nelvin,
    const int *nelb,
    const int *nelgv,
    const int *nelgt,
    const double *xml,
    const double *yml,
    const double *zml,
    const int *comm_int
){
    std::string configFile="config/config.xml";
    MPI_Comm comm = MPI_Comm_f2c(*comm_int);
    adios = adios2::ADIOS(configFile, comm);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    // In config, writer is opening the insituMPI engine. 
    // I guess this is due to the fact that this is just a write 
    // to let the compressor executable access the data.
    io = adios.DeclareIO("writer");
    // In config, writer 0 uses the BPfile engine, wich i guess stands for
    // bzip2
    io_head = adios.DeclareIO("writer0");
    if (!io.InConfigFile())
    {
        // if not defined by user, we can change the default settings
        // BPFile is the default engine
        io.SetEngine("BPFile");
        io.SetParameters({{"num_threads", "1"}});

        // ISO-POSIX file output is the default transport (called "File")
        // Passing parameters to the transport
    }
    ior = adios.DeclareIO("inputReader");
    // In config, writer 0 uses the BPfile engine, wich i guess stands for
    // bzip2
    if (!ior.InConfigFile())
    {
        // if not defined by user, we can change the default settings
        // BPFile is the default engine
        ior.SetEngine("BPFile");
        ior.SetParameters({{"num_threads", "1"}});

        // ISO-POSIX file output is the default transport (called "File")
        // Passing parameters to the transport
    }
    // Understand how the elements are divided among ranks.
    // set the pointers to the starting points in vectors of for each rank.
    //unsigned int nelv = static_cast<unsigned int>((*nelgv) / size);
    unsigned int nelv = static_cast<unsigned int>((*nelvin));
    unsigned int start = static_cast<unsigned int>(*nelb);
    //if((*nelgv)%size != 0){
    //    if(rank < ((*nelgv) % size)){
    //        ++nelv;
    //        start += static_cast<unsigned int>(rank);
    //    }else{
    //        start += static_cast<unsigned int>((*nelgv)%size);
    //    }
    //}
    start *= static_cast<unsigned int>(*nval);
    // n is count, i.e number of entries in the array in my rank
    unsigned int n = static_cast<unsigned int> (*nval) * nelv;
    // gn is the total size of the arrays, not per io rank 
    unsigned int gn = static_cast<unsigned int>((*nelgv)*(*nval));
    std::cout << rank << ": " << gn << ", " << start << "," << n << std::endl;
    // As I understand, here I am creating a big vector and I am saying
    // where in the vector my rank will write information to. This is
    // done to share information between executables I guess.
    x = io_head.DefineVariable<double>("X", {gn}, {start}, {n});
    y = io_head.DefineVariable<double>("Y", {gn}, {start}, {n});
    z = io_head.DefineVariable<double>("Z", {gn}, {start}, {n});
    p = io.DefineVariable<double>("P_OUT", {gn}, {start}, {n});
    vx = io.DefineVariable<double>("VX_OUT", {gn}, {start}, {n});
    vy = io.DefineVariable<double>("VY_OUT", {gn}, {start}, {n});
    vz = io.DefineVariable<double>("VZ_OUT", {gn}, {start}, {n});
    bm1 = io.DefineVariable<double>("BM1", {gn}, {start}, {n});

    // Do everything again for the temperature (I do not know why)
    //unsigned int nelt = static_cast<unsigned int>((*nelgt) / size);
    unsigned int nelt = static_cast<unsigned int>((*nelvin));
    start = static_cast<unsigned int>(*nelb);
    //if((*nelgt)%size != 0){
    //    if(rank < ((*nelgt) % size)){
    //        ++nelt;
    //        start += static_cast<unsigned int>(rank);
    //    }else{
    //        start += static_cast<unsigned int>((*nelgt)%size);
    //    }
    //}
    ifile = 0 ;
    ifilew = 0 ;
    start = start * static_cast<unsigned int>(*nval);
    n = static_cast<unsigned int> (*nval) * nelt;
    gn = static_cast<unsigned int>((*nelgt)*(*nval));
    // t = io.DefineVariable<double>("T", {gn}, {start}, {n});


    // Do everything again for the global indices 
    //unsigned int nelt = static_cast<unsigned int>((*nelvin));
    nelt = static_cast<unsigned int>((*nelvin));
    start = static_cast<unsigned int>(*nelb);
    //if((*nelgt)%size != 0){
    //    if(rank < ((*nelgt) % size)){
    //        ++nelt;
    //        start += static_cast<unsigned int>(rank);
    //    }else{
    //        start += static_cast<unsigned int>((*nelgt)%size);
    //    }
    //}
    ifile = 0 ;
    ifilew = 0 ;
    start = start;
    n = static_cast<unsigned int> (nelt);
    gn = static_cast<unsigned int>((*nelgt));
    lglelw = io.DefineVariable<int>("LGLEL_OUT", {gn}, {start}, {n});

    // Write the mesh information only once.
    writer_head = io_head.Open("geo.bp", adios2::Mode::Write);
    writer_head.Put<double>(x, xml);
    writer_head.Put<double>(y, yml);
    writer_head.Put<double>(z, yml);
    writer_head.Close();
    if(!rank)
	std::cout << "geo.bp written" << std::endl;
    // Open a global array that will be used by the compressor. 
    // As I see it, this is a buffer, where the compressor will read
    // and compress.
    //writer = io.Open("globalArray", adios2::Mode::Write);
}

extern "C" void adios2_setup_2d_(
    const int *nval,
    const int *nelgv,
    const int *nelgt,
    const double *xml,
    const double *yml,
    const int *comm_int
){
    // For comments here, just check the 3D version on top
    std::string configFile="config/config.xml";
    MPI_Comm comm = MPI_Comm_f2c(*comm_int);
    adios = adios2::ADIOS(configFile, comm);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    io = adios.DeclareIO("writer");
    io_head = adios.DeclareIO("writer0");
    if (!io.InConfigFile())
    {
        // if not defined by user, we can change the default settings
        // BPFile is the default engine
        io.SetEngine("BPFile");
        io.SetParameters({{"num_threads", "1"}});

        // ISO-POSIX file output is the default transport (called "File")
        // Passing parameters to the transport
    }
    unsigned int nelv = static_cast<unsigned int>((*nelgv) / size);
    unsigned int start = nelv * static_cast<unsigned int>(rank);
    if((*nelgv)%size != 0){
        if(rank < ((*nelgv) % size)){
            ++nelv;
            start += static_cast<unsigned int>(rank);
        }else{
            start += static_cast<unsigned int>((*nelgv)%size);
        }
    }
    start *= static_cast<unsigned int>(*nval);
    unsigned int n = static_cast<unsigned int> (*nval) * nelv;
    unsigned int gn = static_cast<unsigned int>((*nelgv)*(*nval));
    //std::cout << rank << ": " << gn << ", " << start << "," << n << std::endl;
    x = io_head.DefineVariable<double>("X", {gn}, {start}, {n});
    y = io_head.DefineVariable<double>("Y", {gn}, {start}, {n});
    // p = io.DefineVariable<double>("P", {gn}, {start}, {n});
    vx = io.DefineVariable<double>("VX", {gn}, {start}, {n});
    // vy = io.DefineVariable<double>("VY", {gn}, {start}, {n});
    unsigned int nelt = static_cast<unsigned int>((*nelgt) / size);
    start = nelt * static_cast<unsigned int>(rank);
    if((*nelgt)%size != 0){
        if(rank < ((*nelgt) % size)){
            ++nelt;
            start += static_cast<unsigned int>(rank);
        }else{
            start += static_cast<unsigned int>((*nelgt)%size);
        }
    }
    start *= static_cast<unsigned int>(*nval);
    n = static_cast<unsigned int> (*nval) * nelt;
    gn = static_cast<unsigned int>((*nelgt)*(*nval));
    // t = io.DefineVariable<double>("T", {gn}, {start}, {n});

    writer_head = io_head.Open("geo.bp", adios2::Mode::Write);
    writer_head.Put<double>(x, xml);
    writer_head.Put<double>(y, yml);
    writer_head.Close();
    //writer = io.Open("globalArray", adios2::Mode::Write);

}

extern "C" void adios2_update_2d_(
    const double *pr,
    const double *v,
    const double *u,
    const double *temp
){  
    // For comments, check the 3D version below
    startT = std::clock();
    writer.BeginStep();
    // writer.Put<double>(p, pr);
    writer.Put<double>(vx, v);
    // writer.Put<double>(vy, u);
    // writer.Put<double>(t, temp);
    writer.EndStep();
    dataTime += (std::clock() - startT) / (double) CLOCKS_PER_SEC;
}

extern "C" void adios2_update_(
    const int *lglel,
    const double *pr,
    const double *v,
    const double *u,
    const double *w,
    const double *temp
){
    startT = std::clock();
    ifilew=ifilew+1;
    std:: string fileName= "out.f0000"+ std::to_string(ifilew) +".bp";
    writer = io.Open(fileName, adios2::Mode::Write);
    // Begin a step of the writer. Remember that this will write to
    // the variable globalarray that will be read by the compressor.
    writer.BeginStep();
    writer.Put<double>(p, pr);
    writer.Put<double>(vx, v);
    writer.Put<double>(vy, u);
    writer.Put<double>(vz, w);
    writer.Put<int>(lglelw, lglel);
    //writer.Put<double>(t, temp);
    writer.EndStep();
    writer.Close();
    //++stepw;
    dataTime += (std::clock() - startT) / (double) CLOCKS_PER_SEC;
}

extern "C" void adios2_stream_(
    const int *lglel,
    const double *pr,
    const double *v,
    const double *u,
    const double *w,
    const double *mass1,
    const double *temp
){
    startT = std::clock();
    // Begin a step of the writer. Remember that this will write to
    // the variable globalarray that will be read by the compressor.
    writer.BeginStep();
    writer.Put<double>(p, pr);
    writer.Put<double>(vx, v);
    writer.Put<double>(vy, u);
    writer.Put<double>(vz, w);
    writer.Put<double>(bm1, mass1);
    writer.Put<int>(lglelw, lglel);
    //writer.Put<double>(t, temp);
    writer.EndStep();
    dataTime += (std::clock() - startT) / (double) CLOCKS_PER_SEC;
}
extern "C" void adios2_read_(
    int *lglelrr,
    double *pr,
    double *v,
    double *u,
    double *w,
    //const double *temp,
    const int *nval,
    const int *nelvin,
    const int *nelb,
    const int *nelgv,
    const int *nelgt,
    const int *comm_int,
    char *fname
){
    startT = std::clock();
    std::string configFile="config/config.xml";
    MPI_Comm comm = MPI_Comm_f2c(*comm_int);
    //adios = adios2::ADIOS(configFile, comm);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    // In config, writer is opening the insituMPI engine. 
    // I guess this is due to the fact that this is just a write 
    // to let the compressor executable access the data.
//    Coment all concerning ior since it was declared in the setup
//    ior = adios.DeclareIO("inputReader");
//    // In config, writer 0 uses the BPfile engine, wich i guess stands for//
//    // bzip2
//    if (!ior.InConfigFile())
//    {
//        // if not defined by user, we can change the default settings
//        // BPFile is the default engine
//        ior.SetEngine("BPFile");
//        ior.SetParameters({{"num_threads", "1"}});
//
//        // ISO-POSIX file output is the default transport (called "File")
//        // Passing parameters to the transport
//    }
    // Understand how the elements are divided among ranks.
    // set the pointers to the starting points in vectors of for each rank.
    //unsigned int nelv = static_cast<unsigned int>((*nelgv) / size);
    unsigned int nelv = static_cast<unsigned int>((*nelvin));
    unsigned int start = static_cast<unsigned int>(*nelb);
    //if((*nelgv)%size != 0){
    //    if(rank < ((*nelgv) % size)){
    //        ++nelv;
    //        start += static_cast<unsigned int>(rank);
    //    }else{
    //        start += static_cast<unsigned int>((*nelgv)%size);
    //    }
    //}
    start *= static_cast<unsigned int>(*nval);
    // n is count, i.e number of entries in the array in my rank
    unsigned int n = static_cast<unsigned int> (*nval) * nelv;
    // gn is the total size of the arrays, not per io rank 
    unsigned int gn = static_cast<unsigned int>((*nelgv)*(*nval));

    // These are the indicators for the global indices
    unsigned int nelv3 = static_cast<unsigned int>((*nelvin));
    unsigned int start3 = static_cast<unsigned int>(*nelb);
    unsigned int n3 = static_cast<unsigned int> (nelv);
    unsigned int gn3 = static_cast<unsigned int>((*nelgv));
    
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
    //std::cout << "rank"<< rank << ": " << gn << ", " << start << "," << n << std::endl;
    
    //std::cout << "Reading the file " << std::endl;

    // This file is here just to correctly advance the reader 
    // I have an error if I send the name from fortram. I would need to trim the string.
    ifile=ifile+1;
    // Make a new name every time the reader is called
    std:: string fileName= "out.f0000"+ std::to_string(ifile) +".bp";
    startT = std::clock();
    // Read the new name
    //adios2::Engine readr = ior.Open(fileName,adios2::Mode::Read,comm); 
    readr = ior.Open(fileName,adios2::Mode::Read); 
    //readr.BeginStep();

    //adios2::Engine readr = ior.Open("out.f00001.bp",adios2::Mode::Read,comm); 
    //adios2::Engine readr = ior.Open(fname,adios2::Mode::Read,comm); 
    int step = 1;
    bool firstStep = true;

    // I am commenting the while loop and all the breaks since they
    // appear to create problems when reading multiple files in succession.
    // the problem is somewhere in the "steps"

    //while (true)
    //{
     	//adios2::StepStatus status =
        //readr.BeginStep(adios2::StepMode::Read);
    //    if (status != adios2::StepStatus::OK)
    //    {
    //        break;
    //    }

        vxr = ior.InquireVariable<double>("VX_OUT");
        if (!vxr)
        {
            std::cout << "Error: NO variable VX_OUT found. Unable to proceed. " << std::endl;
    //        break;
        }

        vyr = ior.InquireVariable<double>("VY_OUT");
        if (!vyr)
        {
            std::cout << "Error: NO variable VY_OUT found. Unable to proceed. " << std::endl;
    //        break;
        }

        vzr = ior.InquireVariable<double>("VZ_OUT");
        if (!vzr)
        {
            std::cout << "Error: NO variable VZ_OUT found. Unable to proceed. " << std::endl;
    //        break;
        }

        prr = ior.InquireVariable<double>("P_OUT");
        if (!prr)
        {
            std::cout << "Error: NO variable P_OUT found. Unable to proceed. " << std::endl;
    //        break;
        }
        lglelr = ior.InquireVariable<int>("LGLEL_OUT");
        if (!lglelr)
        {
            std::cout << "Error: NO variable LGLEL_OUT found. Unable to proceed. " << std::endl;
    //        break;
        }

    elapsedT = (std::clock() - startT) / (double) CLOCKS_PER_SEC;
    std::cout <<  "rank: " << rank << " open and inquire: " << elapsedT << "s." << std::endl;
    startT = std::clock();
        if (firstStep)
        {
        	// Promise that we are not going to change the variable sizes
                // nor add new variables
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
 
        // std::cout << "Select how much is yours " 
           //               << std::endl;
        vxr.SetSelection({{start}, {n}});
        vyr.SetSelection({{start}, {n}});
        vzr.SetSelection({{start}, {n}});
        prr.SetSelection({{start}, {n}});
        lglelr.SetSelection({{start3}, {n3}});

    elapsedT = (std::clock() - startT) / (double) CLOCKS_PER_SEC;
    std::cout <<  "rank: " << rank << " set selection: " << elapsedT << "s." << std::endl;
    startT = std::clock();
 
        // std::cout << "get the data "
           //              << std::endl;
        readr.Get<double>(vxr,vVXr.data());
        readr.Get<double>(vyr,vVYr.data());
        readr.Get<double>(vzr,vVZr.data());
        readr.Get<double>(prr,vVPrr.data());
        readr.Get<int>(lglelr,vVLGLELr.data());
        //readr.Get<double>(vzr,w);

        //readr.EndStep();
        //++step;
    //}
    elapsedT = (std::clock() - startT) / (double) CLOCKS_PER_SEC;
    std::cout <<  "rank: " << rank << " get data: " << elapsedT << "s." << std::endl;
    startT = std::clock();

    readr.Close();	
    //++step;
    elapsedT = (std::clock() - startT) / (double) CLOCKS_PER_SEC;
    std::cout <<  "rank: " << rank << " Close: " << elapsedT << "s." << std::endl;
    startT = std::clock();

    //copy the data into w
    if (rank==0){
    std::cout <<"Copying data into nek vector " << std::endl;
    }
    //if (rank==0){
    for (int i=0; i<vVZr.size(); i++){
    //    std::cout <<"entry "<< i <<"is " << vVZr[i] << " ";
    //    std::cout <<"size of read vector " <<"is " << vvzr.size() << std::endl;
    pr[i]=vVPrr[i];
    v[i]=vVXr[i];
    u[i]=vVYr[i];
    w[i]=vVZr[i];
    }
    for (int i=0; i<vVLGLELr.size(); i++){
    //    std::cout <<"entry "<< i <<"is " << vVZr[i] << " ";
    //    std::cout <<"size of read vector " <<"is " << vvzr.size() << std::endl;
    lglelrr[i]=vVLGLELr[i];
    }
    //}
    elapsedT = (std::clock() - startT) / (double) CLOCKS_PER_SEC;
    std::cout <<  "rank: " << rank << " copy into nek: " << elapsedT << "s." << std::endl;
    startT = std::clock();

    dataTime += (std::clock() - startT) / (double) CLOCKS_PER_SEC;
    //}

}

extern "C" void adios2_finalize_(){
    writer.Close();
    std::cout <<  "rank: " << rank << " in-situ time: " << dataTime << "s." << std::endl;
}
