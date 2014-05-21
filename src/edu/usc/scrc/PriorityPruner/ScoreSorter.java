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
 * Sorts Result-objects in regards to if they are force included, have passed
 * design score, their number of bead types, their score and name. (Used by
 * Pruner.)
 */
public class ScoreSorter implements Comparator<Result> {
	@Override
	public int compare(Result x, Result y) {
		if (x.getPartnerSnp().getForceInclude() == y.getPartnerSnp().getForceInclude()) {
			if (x.getPartnerSnp().getScore() == y.getPartnerSnp()
					.getScore()) {
				return (y.getPartnerSnp().getSnpName().compareTo(x
						.getPartnerSnp().getSnpName()));
			} else {
				return (Double.compare(y.getPartnerSnp().getScore(), x
						.getPartnerSnp().getScore()));
			}

		} else {
			return (Boolean.valueOf(y.getPartnerSnp().getForceInclude())
					.compareTo(Boolean.valueOf(x.getPartnerSnp()
							.getForceInclude())));
		}
	}
}