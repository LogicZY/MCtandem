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

#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <sys/timeb.h>
#include <ctime>

#include "boost.h"
#include "stdafx.h"
#include "msequence.h"
#include "msequencecollection.h"
#include "msequenceserver.h"
#include "msequtilities.h"
#include "mspectrum.h"
#include "mhistogram.h"
#include "xmlparameter.h"
#include "mscore.h"
#include "mprocess.h"
#include "serialize.h"

namespace boost {
  namespace serialization {

  template<class Archive>
  void serialize(Archive &ar, mprocess &p, const unsigned int version)
    {
      // omitted m_prcLog
      // omitted m_xmlPerformance
      // omitted m_xmlValues
      ar & p.m_vSpectra;
	ar & p.m_mapSequences; // a map containing all of the protein sequences discovered, indexed by their m_tUid value
	ar & p.m_vseqBest; // a vector of msequences used in the model refinement process
	ar & p.m_vstrModifications; //a vector containing the strings defining fixed modifications for a protein
	
	ar & p.m_tRefineModels; // total number of models generated by refinement
	ar & p.m_tRefineInput; // total number of sequences included in a refinement session
	ar & p.m_tRefinePartial; // the number of models discovered to have partial cleavage
	ar & p.m_tRefineUnanticipated; // the number of models discovered to have unanticpated cleavage
	ar & p.m_tRefineNterminal; // the number of models discovered to have modified N-terminii
	ar & p.m_tRefineCterminal; // the number of models discovered to have modified C-terminii
	ar & p.m_tRefinePam; // the number of models discovered to have point mutations
	ar & p.m_dRefineTime; // the time required to perform a refinement
	ar & p.m_tActive;	// total number of models remaining after each refinement step
	ar & p.m_bRefineCterm;  //true if processing 'refine, potential C-terminus modifications'. Set in mrefine::refine and 

	ar & p.m_viQuality; // contains the data quality scoring vector
	ar & p.m_bReversedOnly;
	ar & p.m_bSaps;
	ar & p.m_bAnnotation;

 	ar & p.m_bMinimalAnnotation;
  	ar & p.m_bSerialize;
 	ar & p.m_bCheckNg;
	ar & p.m_dNt;
 	ar & p.m_dNtAve;
 	ar & p.m_dNg;
 	ar & p.m_dNgAve;

	ar & p.m_strLastMods;
	ar & p.m_iCurrentRound;
	ar & p.m_bPermute;
	ar & p.m_bPermuteHigh;
	ar & p.m_bCrcCheck;
	ar & p.m_bQuickAcetyl;
	ar & p.m_bQuickPyro;
	ar & p.m_dEsum;

	// RDB TODO CAN I OMIT THIS?omitted ar & p.m_setRound;
	ar & p.m_vstrSaps;
	ar & p.m_vstrMods;
	ar & p.m_mapAnnotation;

	// omitted msemistate m_semiState; // maintains the state of the semi-enzymatic cleavage state machine
	// omitted mpyrostate m_pyroState; // maintains the state of the pyrolidone carboxylic acid detection state machine
	// omitted m_errValues;
	ar & p.m_dSearchTime; // total time elapsed during a protein modeling session process
	ar & p.m_lIonCount; // minimum sum of detected ions that are significant enough to store a sequence 
	ar & p.m_lThread; // thread number of this object
	ar & p.m_lThreads; // the total number of threads current active
	ar & p.m_lReversed; // the total number of peptides found where the reversed sequence was better than the forward sequence
	ar & p.m_dThreshold; // the current expectation value threshold
	ar & p.m_tContrasted; // the number of spectra subtracted using contrast angle redundancy detection
	ar & p.m_lStartMax; // set the maximum distance from N-terminus for a peptide
					  // normally set at an impossibly large value = 100,000,000
					  // for ragged N-terminii with potential modifications, set at a low but plausible value = 50
	ar & p.m_lCStartMax;
	// omitted char *m_pSeq; // a character pointer, used for temporary sequence information
	ar & p.m_bUn; // if true, cleave at all residues. if false, use cleavage specification in data input.
	ar & p.m_bUseHomologManagement; // set to true to use homologue management 
	ar & p.m_tMinResidues; // the minimum peptide length that will be scored
	ar & p.m_tMissedCleaves; // the maximum number of cleavage sites that can be missed
	ar & p.m_tPeptideCount; // the total number of peptide sequences generated during a process
	ar & p.m_tPeptideScoredCount; // the total number of peptide sequences scored during a process
	ar & p.m_tProteinCount; // the total number of protein sequences considered during a process
	ar & p.m_tSpectra; // the total number of spectra being modeled
	ar & p.m_tSpectraTotal; // the total number of spectra in the input file
	ar & p.m_tValid; // the number of valid peptide models
	ar & p.m_tTotalResidues; // the number of residues read
	ar & p.m_tSeqSize; // current length of the m_pSeq character array
	ar & p.m_tUnique; // the number of unique peptides found in a result
	string m_strOutputPath; // the path name of the XML output file
	// omitted mcleave m_Cleave; // the specification for a cleavage peptide bond
	// omitted msequence m_seqCurrent; // the msequence object that is currently being scored
#ifdef X_P3
	// omitted p3msequenceServer m_svrSequences; // the msequenceServer object that provides msequences to msequenceCollection
#else
	// omitted msequenceServer m_svrSequences; // the msequenceServer object that provides msequences to msequenceCollection
#endif
	// omitted mspectrumcondition m_specCondition; // the mspectrumcondition object that cleans up and normalized
										// spectra for further processing
	// omitted mscore* m_pScore; // the object that is used to score sequences and spectra
	// omitted mrefine* m_pRefine; // the object that is used to refine models

    }

  template<class Archive>
  void serialize(Archive &ar, maa &m, const unsigned int version)
    {
	ar & m.m_lPos; // the sequence position of the residue (N-terminal = 0)
	ar & m.m_dMod; // mass of the modification
	ar & m.m_cRes; // single letter abbreviation for the amino acid
	ar & m.m_cMut; // single letter abbreviation for a discovered point mutation
	ar & m.m_strId; // character string representing an external accession number for a mutation/modification
        ar & m.m_dPrompt; // prompt loss from modification mass
    }
  template<class Archive>
  void serialize(Archive &ar, mspectrum &s, const unsigned int version)
    {
	ar & s.m_tId; // an identification number
	ar & s.m_tCurrentSequence; // an identifier for the current sequence (used in mprocess)
	ar & s.m_uiType; 
	ar & s.m_fScore; // the convolution score
	ar & s.m_fHyper; // the hyper score
	ar & s.m_fScoreNext; // next best convolution score
	ar & s.m_fHyperNext; // next best hyper score
	ar & s.m_dRatio; // 
	ar & s.m_dExpect; // the expectation value
	ar & s.m_dProteinExpect; // the expectation value for the associated protein
	ar & s.m_dMH; // the parent ion mass + a proton
	ar & s.m_fI; // the parent ion intensity (if available)
	ar & s.m_fZ; // the parent ion charge
	ar & s.m_bRepeat; // a flag indicating that a better match for an individual peptide has already been found
	ar & s.m_bActive; // a flag indicating that a spectrum is available for scoring
	ar & s.m_vMI; // a vector containing the m/z - intensity information for fragment ions 
	ar & s.m_vMINeutral; // a vector containing the m/z - intensity information for fragment ions 
	ar & s.m_vseqBest; // a vector containing the highest scoring msequence objects
	ar & s.m_vdStats;
	ar & s.m_strDescription;
 	ar & s.m_strRt;

	ar & s.m_hHyper; // the histogram of hyper scores
	ar & s.m_hConvolute; // the histogram of convolution scores
	ar & s.m_chBCount; // the histogram of b-ion counts
	ar & s.m_chYCount; // the histogram of y-ion counts
    }
  template<class Archive>
  void serialize(Archive &ar, mi &m, const unsigned int version)
    {
  	ar & m.m_fM; // the m/z value
        ar & m.m_fI; // the intensity value
    }
  template<class Archive>
  void serialize(Archive &ar, mdomain &d, const unsigned int version)
    {
	ar & d.m_lS; // the start position of the peptide in the protein sequence (N-terminus = 0)
	ar & d.m_lE; // the end position of the peptide in the protein sequence
	ar & d.m_lMissedCleaves; // missed cleavages
	ar & d.m_fScore; // the convolution score for the peptide
        ar & d.m_fHyper; // the hyper score for the peptide 
	ar & d.m_dMH; // the mass of the peptide + a proton
	// double m_dDelta replaces float m_fDelta, starting with version 2006.02.01
	// because of an issue with the accuracy of this value
	ar & d.m_dDelta; // the mass difference between the mass of the peptide and the measured mass
	ar & d.m_bUn;
        ar & d.m_mapCount; // a map of the number of ions detected for each ion type
	ar & d.m_mapScore; // a map of the convolution scores for each ion type
	ar & d.m_vAa; // vector of modified amino acids
    }
}
}

std::string save_process(mprocess &p){
  std::ostringstream oss;
#if 0
  std::ofstream out_file("dump.out");
#endif

  boost::archive::text_oarchive oa(oss);
  oa << p;
  std::string s = oss.str();
#if 0
  const char *dummy=s.c_str();
  int i = strlen(dummy);
#endif
  return s;
}

void load_process(char *data, mprocess &p){
  std::string s(data);
  std::istringstream iss(s);
  boost::archive::text_iarchive ia(iss);
  ia >> p;
}

#if 0
void dummy_load_process(mprocess &p){
  std::ifstream ifs("dump.out");
  boost::archive::text_iarchive ia(ifs);
  ia >> p;
}
#endif

