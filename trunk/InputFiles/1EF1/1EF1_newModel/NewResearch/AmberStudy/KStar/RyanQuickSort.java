////////////////////////////////////////////////////////////////////////////////////////////
// RyanQuickSort.java
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
 *
 * This software written by Ryan Lilien (2001-2004)
 * 
 * For some reason Microsoft doesn't implement any sorting functions
 *  so I have to do it myself. This doesn't make any sense, especially
 *  because the standard Java specifications include sorting.
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

public class RyanQuickSort
{

	Object dArray[] = null;

	RyanQuickSort() {}

	// Sorts the array a. Objects in array a must implement
	//  RyanComparable
	public void Sort(Object[] a) {
		
		RyanQuickSort r = new RyanQuickSort();
		r.dArray = a;
		if (r.isAlreadySorted()) return;
		r.quickSort(0, r.dArray.length-1);
		if (r.isAlreadySorted())
			return;
		else {
			System.out.println("SORT FAILED!");
		}
		return;	
	}

	private void quickSort(int q, int w) {
		if (q < w) {
			int p = partition(q,w);
			if (p == w)
				p--;
			quickSort(q,p);
			quickSort(p+1,w);
		}
	}

	private int partition(int b, int t) {
		Object pivot = dArray[b];
		while(true) {
			while ( (((RyanComparable)dArray[t]).compareTo(pivot) >= 0 ) && (b < t) )
				t--;
			while ( (((RyanComparable)dArray[b]).compareTo(pivot) < 0 ) && (b < t) )
				b++;
			if (b < t) {
				// exchange
				Object tmp = dArray[b];
				dArray[b] = dArray[t];
				dArray[t] = tmp;
			}
			else
			 return t;
		}
	}

	// Returns true if array is sorted
	private boolean isAlreadySorted() {
		for(int i=1;i<dArray.length;i++){
			if (((RyanComparable)dArray[i-1]).compareTo(dArray[i]) > 0)
				return false;
		}
		return true;
	}

}
