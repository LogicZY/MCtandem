ssh mic0;
mkdir tandem-mic;
mkdir raw;
logout;
cd /home/lichuang/tandem-mic;
scp -r bin fasta mic0:/tmp/tandem-mic;
cd /home/lichaung/raw;
scp JNA-and-JGP-CAS-27Aug03_HUPO.mzXML Adult_Heart_Gel_Elite_54_f04.mgf HUMAN.fasta mic0:/tmp/raw;
cd /home/lichuang/expat/lib;
scp libexpat.so libexpat.so.1 libexpat.so.1.6.0 mic0:/tmp/tandem-mic;

