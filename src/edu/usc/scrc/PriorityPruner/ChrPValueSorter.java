/**
Copyright (c) 2014 Christopher K. Edlund, Malin Anker, Fredrick R. Schumacher, W. James Gauderman, David V. Conti,
University of Southern California,
Los Angeles, CA  90033, USA.

The MIT License (MIT)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 */

package edu.usc.scrc.PriorityPruner;

import java.util.Comparator;

/**
 * This sorter sorts SNPs by chromosome and base pair position. (It is being
 * used by SnpListFile when all chromosomes are included in the run.)
 */
public class ChrPValueSorter implements Comparator<SnpInfo> {

	@Override
	public int compare(SnpInfo x, SnpInfo y) {
		if (x.getChr().equals(y.getChr())) {
			if ((Double.compare(x.getPValue(), y.getPValue())) == 0) {
				// the name check was just added for comparison reasons
				return (x.getSnpName().compareTo(y.getSnpName()));
			} else {
				return Double.compare(x.getPValue(), y.getPValue());
			}
		} else {
			// first trying to parse the chromosomes as integers, otherwise they'll
			// get compared as Strings
			int xInt;
			int yInt;
			try {
				xInt = Integer.parseInt(x.getChr());
				yInt = Integer.parseInt(y.getChr());
				return Integer.valueOf(xInt).compareTo(Integer.valueOf(yInt));
			} catch (NumberFormatException e) {
				return x.getChr().compareTo(y.getChr());
			}
		}
	}
}