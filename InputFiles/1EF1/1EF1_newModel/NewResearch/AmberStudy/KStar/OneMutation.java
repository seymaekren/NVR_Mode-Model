////////////////////////////////////////////////////////////////////////////////////////////
// OneMutation.java
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
 * Written by Ryan Lilien (2002-2004)
 *
 * One mutation with its estimated energy
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

import java.math.BigDecimal;

public class OneMutation implements RyanComparable
{

	int mutNum = -1;
	BigDecimal score = new BigDecimal("0.0");
	float vol = 0.0f;
	String resTypes[] = null;
	int resMut[] = null;
	String flagMutType = null;
	
	OneMutation() {
	}
	
	/*public int compareTo(Object otherObject) {
		OneMutation mut = (OneMutation)otherObject;
		if (score > mut.score) return -1;
		if (score < mut.score) return 1;
		return 0;
	}*/
	
	public int compareTo(Object otherObject){
		OneMutation mut = (OneMutation)otherObject;
		String seq1 = "";
		String seq2 = "";
		if (resTypes!=null){
			for (int i=0; i<resTypes.length; i++){
				seq1 += resTypes[i];
				seq2 += mut.resTypes[i];
			}
		}
		else {
			for (int i=0; i<resMut.length; i++){
				seq1 += resMut[i];
				seq2 += mut.resMut[i];
			}
		}
		return (seq1.compareTo(seq2));
	}
	
	// Returns true if the passed mutation sequence and this
	//  mutation sequence are the same. Otherwise returns 0.
	public boolean isSame(String anotherMutation[]) {
		for(int i=0;i<anotherMutation.length;i++) {
			if (!anotherMutation[i].equalsIgnoreCase(resTypes[i]))
				return(false);
		}
		return(true);
	}	

}
