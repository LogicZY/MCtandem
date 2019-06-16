# 2019.6.15 updated.

# MCtandem is a parallel version of X!Tandem. We mainly modified the mscore.cpp, tandem.cpp and Makefile files. 

Instead of using .rar files, we have added the individual files to the repository, organizing things using folders. 
1. Modify the comments in the source code
2. Add script command file.
3. Comment the location of some optimization strategies added to the code.

#Please wait for us to issue the updates. You will be getting them in at most 14 days.

# 2018.7.6 updated.
#MCtandem

#An MIC cluster-enabled spectrum dot product algorithm for large-scale peptide identification

MCtandem is a parallel peptide database search sequencing tool based on Many Integrated Core (MIC) cluster. It significantly increases the efficiency of large-scale database search sequencing and thus is significant in deriving peptide sequences from mass spectrometry data.

#Installation and Requirements

Hardware requirements: MIC cluster. Each Intel MIC workstation besides containing an Intel Xeon E5-2640 v2 processor with 8 physical cores clocked at 2.0 GHz and 8G RAM, and an Intel Xeon Phi Coprocessor 7120P series with 61 cores at 1.33 GHz and 16G RAM. The running OS is NeoKylin 3.2, 64 bits.

A proper process/thread/memory affinity is the basis for optimal performance. Therefore, some default setting needs to be modified. The details of the configuration parameters are as follows:

Module load craype-hunepages2M

Export MKL_FAST_MEMORY_LIMIT = 0

Export OMP_PROC_BIND = TRUE

Export OMP_PLACES = threads

Export OMP_STACKSIZE = 512m

Export OMP_NUM_THREADS = 16

#Usage:

Single

Requirements:

ssh mic0;

mkdir tandem-mic

mkdir raw;

logout;

cd /home/lichuang/tandem-mic;

scp -r bin fasta mic0:/tmp/tandem-mic;

cd /home/lichaung/raw;

scp JNA-and-JGP-CAS-27Aug03_HUPO.mzXML Adult_Heart_Gel_Elite_54_f04.mgf HUMAN.fasta mic0:/tmp/raw;

cd /home/lichuang/expat/lib;

scp libexpat.so libexpat.so.1 libexpat.so.1.6.0 mic0:/tmp/tandem-mic;

eg:/MCtandem.exe input-path output-path or<param.xml>

Cluster

Requirements:

To build MCtandem, the following additional packages are required:

MPI.
Boost Serialization library.
eg:mpiexec -n N -machinefile ~./machinefile.txt ./MCtandem.exe input-path output-path
