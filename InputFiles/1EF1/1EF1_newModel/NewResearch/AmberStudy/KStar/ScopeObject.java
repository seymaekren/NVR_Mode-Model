////////////////////////////////////////////////////////////////////////////////////////////
// ScopeObject.java
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
 * Written by Ryan Lilien (2001-2004)
 * 
 * This object is a wrapper for an object that goes into a scope
 *    hashtable
 *  i.e. it can wrap systems, molecules, atoms etc...
 * 
 *  Object types can be determined by looking at
 *    object.getClass().getName() which is a string and comparing
 *    to known types. ie. Molecule, Atom, java.lang.String,
 *	  java.lang.Integer
 *  Or by looking at the objType member:
 *  -1 = invalid (this object is not a valid object)
 *   0 = undefined
 *   1 = system     this is a molecule in the B class hierarchy
 *   2 = molecule   this is a strand in the B class hierarchy
 *					the secondary data here is the associated system
 *   3 = residue
 *   4 = atom
 *   5 = bond
 *  10 = string
 *  11 = integer
 *  12 = double
 *  13 = boolean
 *  20 = forcefield
 *					the secondary data here is the associated system
 *  21 = torsionList
 *  22 = torsion
 *  23 = grid
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

public class ScopeObject
{
	Object data;   // the object itself
	Object secondaryData;   // secondary data
	String objName;  // name of the object
	String objNotes;  // notes about the object
	int objType;   // the type of object
	
	ScopeObject(Object obj, String name, int numType) {
		data = obj;
		secondaryData = null;
		objName = name;
		objType = numType;
		objNotes = new String("");
	}
	ScopeObject(Object obj, String name) {
		data = obj;
		secondaryData = null;
		objName = name;
		objType = 0;
		objNotes = new String("");
	}
	ScopeObject() {
		data = null;
		secondaryData = null;
		objType = -1;
	}
	ScopeObject(Object obj, String name, String notes) {
		data = obj;
		secondaryData = null;
		objName = name;
		objNotes = notes;
		objType = 0;
	}
	ScopeObject(Object obj, String name, String notes, int numType) {
		data = obj;
		secondaryData = null;
		objName = name;
		objNotes = notes;
		objType = numType;
	}
}
