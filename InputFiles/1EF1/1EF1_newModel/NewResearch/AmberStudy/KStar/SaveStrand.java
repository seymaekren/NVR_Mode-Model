////////////////////////////////////////////////////////////////////////////////////////////
// SaveStrand.java
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
 * SaveStrand created by Ryan Lilien (2001-2004),
 *   modified from SaveMolecule by Neill White
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

import java.io.PrintStream;
import java.io.BufferedInputStream;
import java.util.Hashtable;

class SaveStrand {

	PrintStream pw;

	// Saves strand number num of molecule m to the printstream pw.
	// Params include
	//  comment: (string) a one line comment for header
	//  printSegID: (boolean) print atom seg ids
	//  showConnect: (boolean) add connect terms after atom entries
	SaveStrand(Molecule m, PrintStream pw, int num, Hashtable params) throws Exception {
		Strand stnd = null;

		// Pull selected strand from molecule
		if (num < m.numberOfStrands) {
			stnd = m.strand[num];
			DoSaveStrand(m, pw, stnd, params);
		}
		
	}

	// Saves strand stnd of molecule m to the printstream pw.
	// Params include
	//  comment: (string) a one line comment for header
	//  printSegID: (boolean) print atom seg ids
	//  showConnect: (boolean) add connect terms after atom entries
	SaveStrand(Molecule m, PrintStream pw, Strand stnd, Hashtable params) throws Exception {
		DoSaveStrand(m, pw, stnd, params);
	}
	
	
	private void DoSaveStrand (Molecule m, PrintStream pw, Strand stnd,
			Hashtable params) throws Exception {
		this.pw = pw;
		int atomCounter = 1;
		int residueCounter = 0;
		
		if (stnd.numberOfAtoms == 0)
			return;
		// Pull out parameters
		String comment = (String) params.get("comment");
		boolean printSegID = ((Boolean) params.get("printSegID")).booleanValue();
		boolean showConnect = ((Boolean) params.get("showConnect")).booleanValue();

		String stndStg = stnd.name.toUpperCase();
		if (printSegID)
			stndStg = stnd.residue[0].atom[0].segID;
		int tmpLen = stndStg.length();
		if (tmpLen != 0)
			stndStg = new String(" "+stndStg.substring(0,1));
		else
			stndStg = new String("  ");
	
		if (m.connectivity12Valid == false)
			m.establishConnectivity(true);

		m.resolveCoordinates();

		pw.println("AUTHOR ");
		pw.println("REMARK   6 " + comment);

		for(int j=0; j<stnd.numberOfResidues; j++){
			Residue residue = stnd.residue[j];
			if (residue.numberOfAtoms == 0)
				continue;
			for(int k=0; k<residue.numberOfAtoms; k++){
				Atom atom = residue.atom[k];
				pw.print("ATOM");
				if (atomCounter < 10)
					pw.print("      ");
				else if (atomCounter < 100)
					pw.print("     ");
				else if (atomCounter < 1000)
					pw.print("    ");
				else if (atomCounter < 10000)
					pw.print("   ");
				else if (atomCounter < 100000)
					pw.print("  ");
				else if (atomCounter < 1000000)
					pw.print(" ");
				atom.modelAtomNumber = atomCounter;
				pw.print(atomCounter++ + " ");
				
				pw.print(getAtomField(atom));
				
				if (residue.name.length() == 1)
					pw.print("  " + residue.name.toUpperCase());
				else if (residue.name.length() == 2)
					pw.print(" " + residue.name.toUpperCase());
				else if (residue.name.length() == 3)
					pw.print(residue.name.toUpperCase());
				else if (residue.name.length() > 3)
					pw.print(residue.name.toUpperCase().substring(0, 4));
				pw.print(stndStg);
				residueCounter = residue.getResNumber();
				if (residueCounter < 10)
					pw.print("   ");
				else if (residueCounter < 100)
					pw.print("  ");
				else if (residueCounter < 1000)
					pw.print(" ");
				pw.print(residueCounter + "    ");
				
				pw.print(coordinate(atom.coord[0]));
				pw.print(coordinate(atom.coord[1]));
				pw.print(coordinate(atom.coord[2]));
				pw.print("  1.00  0.00      ");  // Print 0.00 for occupancy and temp factor
				if (printSegID)
					pw.println(atom.segID);
				else
					pw.println();
				
			} // end for atom loop
		} // end for residue loop

		if (showConnect){
			for(int i=0; i<m.numberOfAtoms; i++) {
				if (m.atom[i].strandNumber == stnd.number) {
					int numberOfConnections = m.connected[i][0];
					if (numberOfConnections > 0)
						pw.print("CONECT");
					else 
						continue;
					printAtom(m.atom[i].modelAtomNumber);	
					for(int j=1; j<=numberOfConnections; j++){
						int bondedToAtomNumber = m.connected[i][j];
						Atom bondedToAtom = m.atom[bondedToAtomNumber];
						printAtom(bondedToAtom.modelAtomNumber);
					}
					pw.println();
				}
			}
		}
		
		pw.println("END");
	}

	private void printAtom(int atomNumber){
		if (atomNumber < 10)
			pw.print("    " + atomNumber );
		else if (atomNumber < 100)
			pw.print("   " + atomNumber );
		else if (atomNumber < 1000)
			pw.print("  " + atomNumber);
		else if (atomNumber < 10000)
			pw.print(" " + atomNumber);
		else if (atomNumber < 100000)
			pw.print(atomNumber);
	}

	private String coordinate(float coord){
		if((coord<0.001) && (coord>-0.001))
			coord = 0.0f;
		String coordString = String.valueOf(coord);
		String intPart = " ";
		String floatPart = " ";
		String returnString = null;
		int radix = coordString.indexOf(".");
		if (radix != -1){
			intPart = coordString.substring(0, radix);
			if (radix != coordString.length())
				floatPart = coordString.substring(radix + 1);
		}
		else 
			return(coordString);
		if (intPart.length() == 1 )
			returnString = "   " + intPart;
		else if (intPart.length() == 2)
			returnString = "  " + intPart;
		else if (intPart.length() == 3)
			returnString = " " + intPart;
		else if (intPart.length() == 4)
			returnString = intPart;
		else if (intPart.length() > 4)
			returnString = intPart.substring(0, 3);
		returnString += ".";
		if (floatPart.length() == 1)
			returnString += floatPart + "00";
		else if (floatPart.length() == 2)
			returnString += floatPart + "0";
		else if (floatPart.length() == 3)
			returnString += floatPart;
		else if (floatPart.length() > 3)
			returnString += floatPart.substring(0, 3);
		return(returnString);
	}
	

	private String getAtomField(Atom at){
		if(at.elementType.length()==1) {
			if (at.name.length() == 1)
				return(" " + at.name + "   ");
			else if (at.name.length() == 2)
				return(" " + at.name + "  ");
			else if (at.name.length() == 3)
				return(" " + at.name + " ");
			else if ((at.name.length() >= 4) && (at.name.charAt(0)>='0')
				&& (at.name.charAt(0)<='9')) {
				if (at.name.length() == 4)
					return(at.name + " ");
				else if (at.name.length() > 4)
					return(at.name.substring(0, 4) + " ");
			}
			else if (at.name.length() > 3)
				return(" " + at.name.substring(0,3) + " ");
		}
		else
		{
			if (at.name.length() == 2)
				return(at.name + "   ");
			else if (at.name.length() == 3)
				return(at.name + "  ");
			else if ((at.name.length() >= 4) && (at.name.charAt(0)>='0')
				&& (at.name.charAt(0)<='9')) {
				if (at.name.length() == 4)
					return(at.name + " ");
				else if (at.name.length() > 4)
					return(at.name.substring(0, 4) + " ");
			}
			else if (at.name.length() > 3)
				return(at.name.substring(0,3) + "  ");
		}
		return("     ");
	}
}
