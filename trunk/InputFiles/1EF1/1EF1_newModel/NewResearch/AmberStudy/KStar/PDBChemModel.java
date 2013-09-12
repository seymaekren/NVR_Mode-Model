////////////////////////////////////////////////////////////////////////////////////////////
// PDBChemModel.java
//
//  Version:           0.1
//
//
// authors:
//    initials    name            organization                email
//   ---------   --------------  ------------------------    ------------------------------
//     RHL        Ryan Lilien     Dartmouth College           ryan.lilien@dartmouth.edu
//
////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Written by Ryan Lilien (2000-2004)
 *
 * Reads a pdb file and creates a molecule object
 *
 * Rewritten by Ryan Lilien based on code by Neill White
 * Many functions have been added, others removed, most have had 
 *  at least some parts rewritten. Code rewrites have often been
 *  major to fix bugs or add functionality.
 * 
 * Based on software copyrighted, 1999, by Neill White. 
 *  The author hereby grants permission to use, copy, modify, and re-distribute
 *  this software and its documentation for any purpose, provided
 *  that existing copyright notices are retained in all copies and that this
 *  notice is included verbatim in any distributions. No written agreement,
 *  license, or royalty fee is required for any of the authorized uses.
 *  Modifications to this software may be distributed provided that
 *  the nature of the modifications are clearly indicated.
 *
 */

/*
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

 Contact Info:
   Bruce Donald
   HB 6211
   Computer Science Department
   Dartmouth College
   Hanover, NH 03755
   brd@cs.dartmouth.edu

 If you use or publish any results derived from the use of this program please cite:
   Ryan H. Lilien, Brian W. Stevens, Amy C. Anderson,
   Bruce R. Donald. "A Novel Ensemble-Based Scoring and
   Search Algorithm for Protein Redesign, and its
   Application to Modify the Substrate Specificity of
   the Gramicidin Synthetase A Phenylalanine Adenylation
   Enzyme." Proc. The Eighth Annual International Conference
   on Research in Computational Molecular Biology (RECOMB),
   San Diego, pp 46-57 (2004). 

 Copyright (C) 2004  Ryan H. Lilien and Bruce R. Donald

 <signature of Bruce Donald>, 16 May, 2004
 Bruce Donald, Professor of Computer Science
*/

import java.io.InputStream;
import java.io.BufferedReader;
import java.io.InputStreamReader;

class PDBChemModel {

  // Reads a pdb file and creates a molecule object
	PDBChemModel (Molecule m, InputStream is) throws Exception {
		
		int modelAtomNumber = 0;
		int residueNumber = 0;
		int strandNumber = 0;
		String atomName = "ZZZ", residueName = "ZZZ", strandName = "0"; 
		String lastResidueName = "ZZZ", fullResidueName = "ZZZZZZZZZ";
		String elementType = "  ";
		String curLine = null;
		String tmpStg = null;
		int tmpInt;
		char[] tmpChr = new char[15];
		float	x = 0f, y = 0f, z = 0f;
		Atom	newAtom;
		boolean newStrandPending = true;  // Will add the first strand with the first atom
		Residue newResidue = null;

		BufferedReader bufread = new BufferedReader(new InputStreamReader(is));

		curLine = bufread.readLine();
		if (curLine == null) { // stop if we've reached EOF
			bufread.close();
		}

	  	while(curLine != null) {
	  		// First pad line to 80 characters
	  		tmpInt = curLine.length();
	  		for (int i=0; i < (80-tmpInt); i++)
	  			curLine += " ";
	
	   		if ((curLine.regionMatches(true,0,"ATOM  ",0,6)) || (curLine.regionMatches(true,0,"HETATM",0,6))) {
					
	   			// Is an ATOM line
				tmpStg = curLine.substring(6,11);  // Snag atom serial number
				tmpStg = tmpStg.trim();
				modelAtomNumber = (new Integer(tmpStg)).intValue();
	   			atomName = curLine.substring(12,16);  // Snag atom name
	   			atomName = atomName.trim();
	   			residueName = curLine.substring(17,20);  // Snag short residue name
	   			residueName = residueName.trim();
	   			residueName = checkHis(residueName);
	   			fullResidueName = curLine.substring(17,26);  // Snag full residue atom name
	   			fullResidueName = fullResidueName.trim();
	
				if (!(fullResidueName.equals(lastResidueName))) {
					if (newResidue != null) {
						if (newStrandPending)
							m.addResidue(strandNumber-1, newResidue);
						else
							m.addResidue(strandNumber, newResidue);
					}
					newResidue = new Residue();
					newResidue.name = residueName;
					newResidue.fullName = fullResidueName;
					if (fullResidueName.length() >= 3)
						newResidue.name = checkHis(fullResidueName.substring(0,3));
					lastResidueName = fullResidueName;
					residueNumber++;
				}
	
				if (newStrandPending) {
					strandName = curLine.substring(21,22);  // Snag strand name
					m.addStrand(strandName);
				 	newStrandPending = false;
				}
				tmpStg = curLine.substring(30,38);  // Snag x coord
				x = (float) new Double(tmpStg).doubleValue();
				tmpStg = curLine.substring(38,46);  // Snag y coord
				y = (float) new Double(tmpStg).doubleValue();
				tmpStg = curLine.substring(46,54);  // Snag z coord
				z = (float) new Double(tmpStg).doubleValue();
	
				elementType = curLine.substring(76,78);  // Snag atom elementType
				elementType = elementType.trim();
				// If we can't get element type from substring(76,78) snag
				//  the first character of the atom name
				if (elementType.equalsIgnoreCase(""))
					elementType = getEleType(curLine.substring(12,15));
				newAtom = new Atom(atomName,x,y,z);
				newAtom.modelAtomNumber = modelAtomNumber;
				newAtom.strandNumber = strandNumber;
				newAtom.elementType = elementType;
				newResidue.addAtom(newAtom);
			} // end ATOM line
			else if (curLine.regionMatches(true,0,"TER   ",0,6)) {
				// Is the end of a strand
				lastResidueName = "ZZZ";
				residueNumber = 0;
				strandNumber++;
				newStrandPending = true;
			} // end TER line
			else   // is a line we skip
				;
			curLine = bufread.readLine();  // attempt to read next line
		}  // end while (curLine != null)
	
		if (newStrandPending)
			m.addResidue(strandNumber-1,newResidue);
		else
			m.addResidue(strandNumber,newResidue);
		bufread.close();  // close the buffer
		
		//Determine the bonds between the atoms in the molecule
		m.determineBonds();
		
		// Assign the molecule relative atom numbers
		m.updateMoleculeAtomNumbers();
	}
	
	// This function pulls the element type from
	//  the atom name
	private String getEleType(String str){
		
		int start=0, end=-1;
		int i=0;
		while( (str.charAt(i)==' ') || ((str.charAt(i)>='0') && (str.charAt(i)<='9')) ) {
			i++;
		}
		start = i;
		end = i++;
		if (i<str.length())
			if((str.charAt(i)>='a') && (str.charAt(i)<='z'))
				end = i;
		return(str.substring(start,end+1));	
	}
	
	//Check if the current residue is some version of HIS, and if so, rename it to HIP
	private String checkHis(String s){
		if ( (s.equalsIgnoreCase("HIS")) || (s.equalsIgnoreCase("HID")) || (s.equalsIgnoreCase("HIP")) || (s.equalsIgnoreCase("HIE")) )
				s = "HIP";
		return s;
	}

}
