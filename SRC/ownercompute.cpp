/*
 Copyright (C) 2007 Robert D Bjornson, all rights reserved
 X!! tandem 
 This software is a component of the X! proteomics software
 development project

Use of this software governed by the Artistic license, as reproduced here:

The Artistic License for all X! software, binaries and documentation

Preamble
The intent of this document is to state the conditions under which a
Package may be copied, such that the Copyright Holder maintains some 
semblance of artistic control over the development of the package, 
while giving the users of the package the right to use and distribute 
the Package in a more-or-less customary fashion, plus the right to 
make reasonable modifications. 

Definitions
"Package" refers to the collection of files distributed by the Copyright 
	Holder, and derivatives of that collection of files created through 
	textual modification. 

"Standard Version" refers to such a Package if it has not been modified, 
	or has been modified in accordance with the wishes of the Copyright 
	Holder as specified below. 

"Copyright Holder" is whoever is named in the copyright or copyrights 
	for the package. 

"You" is you, if you're thinking about copying or distributing this Package. 

"Reasonable copying fee" is whatever you can justify on the basis of 
	media cost, duplication charges, time of people involved, and so on. 
	(You will not be required to justify it to the Copyright Holder, but 
	only to the computing community at large as a market that must bear 
	the fee.) 

"Freely Available" means that no fee is charged for the item itself, 
	though there may be fees involved in handling the item. It also means 
	that recipients of the item may redistribute it under the same
	conditions they received it. 

1. You may make and give away verbatim copies of the source form of the 
Standard Version of this Package without restriction, provided that 
you duplicate all of the original copyright notices and associated 
disclaimers. 

2. You may apply bug fixes, portability fixes and other modifications 
derived from the Public Domain or from the Copyright Holder. A 
Package modified in such a way shall still be considered the Standard 
Version. 

3. You may otherwise modify your copy of this Package in any way, provided 
that you insert a prominent notice in each changed file stating how and 
when you changed that file, and provided that you do at least ONE of the 
following: 

a.	place your modifications in the Public Domain or otherwise make them 
	Freely Available, such as by posting said modifications to Usenet 
	or an equivalent medium, or placing the modifications on a major 
	archive site such as uunet.uu.net, or by allowing the Copyright Holder 
	to include your modifications in the Standard Version of the Package. 
b.	use the modified Package only within your corporation or organization. 
c.	rename any non-standard executables so the names do not conflict 
	with standard executables, which must also be provided, and provide 
	a separate manual page for each non-standard executable that clearly 
	documents how it differs from the Standard Version. 
d.	make other distribution arrangements with the Copyright Holder. 

4. You may distribute the programs of this Package in object code or 
executable form, provided that you do at least ONE of the following: 

a.	distribute a Standard Version of the executables and library files, 
	together with instructions (in the manual page or equivalent) on 
	where to get the Standard Version. 
b.	accompany the distribution with the machine-readable source of the 
	Package with your modifications. 
c.	give non-standard executables non-standard names, and clearly 
	document the differences in manual pages (or equivalent), together 
	with instructions on where to get the Standard Version. 
d.	make other distribution arrangements with the Copyright Holder. 

5. You may charge a reasonable copying fee for any distribution of 
this Package. You may charge any fee you choose for support of 
this Package. You may not charge a fee for this Package itself. 
However, you may distribute this Package in aggregate with other 
(possibly commercial) programs as part of a larger (possibly 
commercial) software distribution provided that you do not a
dvertise this Package as a product of your own. You may embed this 
Package's interpreter within an executable of yours (by linking); 
this shall be construed as a mere form of aggregation, provided that 
the complete Standard Version of the interpreter is so embedded. 

6. The scripts and library files supplied as input to or produced as 
output from the programs of this Package do not automatically fall 
under the copyright of this Package, but belong to whomever generated 
them, and may be sold commercially, and may be aggregated with this 
Package. If such scripts or library files are aggregated with this 
Package via the so-called "undump" or "unexec" methods of producing 
a binary executable image, then distribution of such an image shall 
neither be construed as a distribution of this Package nor shall it 
fall under the restrictions of Paragraphs 3 and 4, provided that you 
do not represent such an executable image as a Standard Version of 
this Package. 

7. C subroutines (or comparably compiled subroutines in other languages) 
supplied by you and linked into this Package in order to emulate 
subroutines and variables of the language defined by this Package 
shall not be considered part of this Package, but are the equivalent 
of input as in Paragraph 6, provided these subroutines do not change 
the language in any way that would cause it to fail the regression 
tests for the language. 

8. Aggregation of this Package with a commercial distribution is always 
permitted provided that the use of this Package is embedded; that is, 
when no overt attempt is made to make this Package's interfaces visible 
to the end user of the commercial distribution. Such use shall not be 
construed as a distribution of this Package. 

9. The name of the Copyright Holder may not be used to endorse or promote 
products derived from this software without specific prior written permission. 

10. THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED 
WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF 
MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE. 

The End 
*/

#include <iostream>

#include "stdafx.h"
#include <algorithm>
#include <set>
#include "msequence.h"
#include "msequencecollection.h"
#include "msequenceserver.h"
#include "msequtilities.h"
#include "mspectrum.h"
#include "loadmspectrum.h"
#include "xmlparameter.h"
#include "mscore.h"
#include "mprocess.h"
#include "saxbiomlhandler.h"
#include "mbiomlreport.h"
#include "mrefine.h"
#include "serialize.h"
#include "timer.h"
#include "ownercompute.h"

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <string>
#include <sstream>

/* This wass necessary for some versions of openmpi because of a conflict between C++ headers and mpi headers
   See the comment in mpicxx.h
*/
#if 0
#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END
#endif

#include "mpi.h"

#define GATHER 1
#define SCATTER 2

/* These global functions are useful for code in other files to use */
string gnumprocs() {
    std::stringstream ss;
    ss << MPI::COMM_WORLD.Get_size();
    return ss.str();
}

int gmyproc() {
    return MPI::COMM_WORLD.Get_rank();
}

void gmyabort(int errorcode) {
    MPI::COMM_WORLD.Abort(errorcode);
}

ownercompute::ownercompute(int &argc, char **&argv) {
    int error_code;
    struct rlimit rlim;
    rlim.rlim_cur = rlim.rlim_max = RLIM_INFINITY;

    MPI::Init(argc,argv);
    num_processes = MPI::COMM_WORLD.Get_size();
    this_process = MPI::COMM_WORLD.Get_rank();
   
    if (this_process == 0) 
    t = new timer("master");
  }

  int ownercompute::me() {
    return this_process;
  }

  string ownercompute::numprocs() {
    std::stringstream ss;
    ss << num_processes;
    return ss.str();
  } 

  /* This will be called by num_processes processes, and each array of mprocess will have num_processes elements.  We want to
     merge all elements into process[0]'s array */
  void ownercompute::gather(mprocess **pProcess) {
    char *buf;
    int len, source;
    MPI::Status status;
    std::string s;

    if (me() > 0) { /* workers */
      /* serialize my data */
      s = save_process(*pProcess[me()]);
      buf = (char *) s.c_str();
      len=s.length()+1;
      /* send to process 0 */
      MPI::COMM_WORLD.Send(buf, len, MPI::CHAR, 0, GATHER);
    } else { /* master */
      t->addsplit("gather");
      for (int i=1; i<num_processes; ++i) {
        t->addsplit("probe");
        MPI::COMM_WORLD.Probe(MPI::ANY_SOURCE, GATHER, status);
        source = status.Get_source();
        len = status.Get_elements(MPI::CHAR);
        buf = new char[len];
        t->addsplit("get data");
        MPI::COMM_WORLD.Recv(buf, len, MPI::CHAR, source, GATHER);
        t->addsplit("got data");
        load_process(buf, *pProcess[source]);
        delete[] buf;
      }
    }
  }

/* master pushes each mprocess out to appropriate worker */
  void ownercompute::scatter(mprocess **pProcess) {
    char *buf;
    int len, i;
    MPI::Status status;
    std::string s;
    if (me() > 0) { /* workers */
        MPI::COMM_WORLD.Probe(0, SCATTER, status);
        len = status.Get_elements(MPI::CHAR);
        buf = new char[len];
        MPI::COMM_WORLD.Recv(buf, len, MPI::CHAR, 0, SCATTER);
        load_process(buf, *pProcess[me()]);
        delete[] buf;
    } else { /* master */
      t->addsplit("scatter");
      for (i=1; i<num_processes; ++i) {
        /* serialize worker data */
        s = save_process(*pProcess[i]);
        buf = (char *) s.c_str();
        /* send to process i */
        len=s.length()+1;
        MPI::COMM_WORLD.Send(buf, len, MPI::CHAR, i, SCATTER);
      }
      t->addsplit("send data");
    }
  }

void ownercompute::finalize() {
  MPI::Finalize();
  if (me() > 0)
    exit(0);
  }

void ownercompute::report() {
  t->addsplit("all done");
  t->print();
  }

