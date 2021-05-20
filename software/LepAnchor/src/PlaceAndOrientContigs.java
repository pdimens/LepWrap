/**
    This file is part of Lep-Anchor.

    Lep-Anchor is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Lep-Anchor is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Lep-Anchor.  If not, see <http://www.gnu.org/licenses/>.

	Copyright (C) 2019 Pasi Rastas, pasi.rastas@gmail.com, University of Helsinki
	
*/

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.*;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

public class PlaceAndOrientContigs {
	
	private ArrayList<ContigInterval> bed = null; // store bed file 
	private HashMap<String, ArrayList<ContigInterval>> bedHash = new HashMap<String, ArrayList<ContigInterval>>(); // get intervals for each contig

	private HashMap<String, ArrayList<ContigInterval>> haplotypeHash = new HashMap<String, ArrayList<ContigInterval>>(); // get intervals for each haplotype contig

	private HashMap<String, Integer> chainLinkHash = new HashMap<String, Integer>(); // get score for each chain link...
	private HashMap<String, Long> chainLinkHashCut1 = new HashMap<String, Long>(); // get cut point for each link
	private HashMap<String, Long> chainLinkHashCut2 = new HashMap<String, Long>(); // get cut point for each link

	private HashMap<String, Long> chainLinkHashCut1a = new HashMap<String, Long>(); // get alternative cut point for each link
	private HashMap<String, Long> chainLinkHashCut2a = new HashMap<String, Long>(); // get alternative cut point for each link

	private HashMap<String, Integer> chainHaplotypeHash = new HashMap<String, Integer>(); // get score for each contig pair for second being haplotype of the first...
	private HashMap<String, long[]> chainHaplotypeCut = new HashMap<String, long[]>(); // get haplotype alignment end points 

	private HashMap<String, Integer> scaffoldingLink_tmp = new HashMap<String, Integer>(); // get scaffolding score linking contigs... 
	private HashMap<String, Integer> scaffoldingLink = new HashMap<String, Integer>(); // get scaffolding score linking contigs... 

	private HashMap<String, String> scaffoldingLinkInfo = new HashMap<String, String>(); // get scaffolding score linking contigs... 
	
	private ArrayList<String> contigs = new ArrayList<String>(); // store contigs
	
	//private HashMap<String, ArrayList<Long>> pafCoverage = new HashMap<String, ArrayList<Long>>();

	private int numErrorsPerContig = 3;
	private int numErrors = 40;

	private int numMaps = 0;
	private ArrayList<ContigInterval> intervalsInChr = null;
	
	private double scaleScore = 0.00001;
	private double orientationPenalty = 0.5;
	private double cutPenalty = 0.001;
	
	private double liftoverScoreDiff = 0.5;
	
	private boolean useChainAndPaf = false;
	
	private String printAnimation = "";
	
	private int maxBridge = 50000;
	private int maxIntersect = 2000;
	
	private int numRuns = 5;
	
	private int numThreads = 8;
	
	private int minHaplotypeAlignmentScore = -10;
	private int minLinkAlignmentScore = -10;
	private boolean keepEmptyIntervals;
	
	private boolean commentOutput; 
	private boolean compressMap;

	private Proximity prox = null;

	
	private void setCommentOutput(boolean commentOutput) {
		this.commentOutput = commentOutput;
	}
	public void setCompressMap(boolean value) {
		compressMap = value;
	}
	public void setPrintAnimation(String fn) {
		printAnimation = fn;
	}
	
	public void setUseChainAndPaf(boolean value) {
		useChainAndPaf = value;
	}

	public void setNumRuns(int runs) {
		numRuns = runs;
	}

	public void setNumErrors(int errors) {
		numErrors = errors;
	}

	public void setNumErrorsPerContig(int errors) {
		numErrorsPerContig = errors;
	}
	
	public void setKeepEmptyIntervals(boolean value) {
		keepEmptyIntervals = value;
	}

	public ArrayList<ArrayList<Marker>> getMarkers(ContigInterval ci_result){
		ArrayList<ArrayList<Marker>> ret =  new ArrayList<ArrayList<Marker>>();
		for (int map = 0; map < numMaps; ++map)
			ret.add(new ArrayList<Marker>());
		
		ArrayList<ContigInterval> cis = bedHash.get(ci_result.getContig());
		if (cis != null) {
			for (ContigInterval c2 : cis) {
				if (Misc.intersectIntervals(c2.getStart(), c2.getEnd(), ci_result.getStart(), ci_result.getEnd()))
					for (int map = 0; map < numMaps; ++map)
						ret.get(map).addAll(c2.getMarkers(map, ci_result.getOrientation()));
			}
		}
		return ret;
	}
	
	// forward calculation, increasing map positions
	private void forward1(ArrayList<Marker> markers, int S[][], int P[][], int startBin, int endBin) {
		forward1(markers, S, P, startBin, endBin, 0, markers.size());
	}
	
	private void forward1(ArrayList<Marker> markers, int S[][], int P[][], int startBin, int endBin, int start, int end)
	{
		for (int mi = start; mi < end; ++mi) {
			Marker m = markers.get(mi);

			//substract startBin
			int maxI  = 0;
			int Si[] = S[mi];
			int Si1[] = S[mi + 1];
			int Pi1[] = P[mi + 1];
			for (int b = 0; b <= endBin - startBin; ++b) { 
                if (Si[b] > Si[maxI])
                    maxI = b;
                Si1[b] = Si[maxI];// + m.inside(b + startBin); 
                Pi1[b] = maxI + startBin; 
            }
			int b = startBin;
			while (b <= endBin) {
				Si1[b - startBin] += m.inside(b);
				b = m.nextInside(b);
			}
		}		
	}
	
	private void forwardFast1(ArrayList<Marker> markers, int S[][], int startBin, int endBin) {
		forwardFast1(markers, S, startBin, endBin, 0, markers.size());
	}
	
	private void forwardFast1(ArrayList<Marker> markers, int S[][], int startBin, int endBin, int start, int end)
	{
		for (int mi = start; mi < end; ++mi) {
			Marker m = markers.get(mi);

			//substract startBin
			int maxS  = 0;
			int Si[] = S[mi];
			int Si1[] = S[mi + 1];
			for (int b = 0; b <= endBin - startBin; ++b) {
				int s = Si[b]; 
                if (s > maxS)
                    maxS = s;
                Si1[b] = maxS;// + m.inside(b + startBin); 
            }
			int b = startBin;
			while (b <= endBin) {
				Si1[b - startBin] += m.inside(b);
				b = m.nextInside(b);
			}

		}		
	}

	private void forward2(ArrayList<Marker> markers, int S[][], int P[][], int startBin, int endBin) {
		forward2(markers, S, P, startBin, endBin, 0, markers.size());
	}
	// forward calculation, decreasing map positions
	private void forward2(ArrayList<Marker> markers, int S2[][], int P2[][], int startBin, int endBin, int start, int end)
	{
		for (int mi = start; mi < end; ++mi) {
			Marker m = markers.get(mi);
			int maxI2  = endBin - startBin;
			int S2i[] = S2[mi];
			int S2i1[] = S2[mi + 1]; 
			int P2i1[] = P2[mi + 1];
			
			for (int b = endBin - startBin; b >= 0; --b) { 
                  
                if (S2i[b] > S2i[maxI2])
                    maxI2 = b;
              
                S2i1[b] = S2i[maxI2];// + m.inside(b + startBin); // add startBin as we subtract it from b 
                P2i1[b] = maxI2 + startBin; 
            }
			int b = startBin;
			while (b <= endBin) {
				S2i1[b - startBin] += m.inside(b);
				b = m.nextInside(b);
			}

		}		
	}

	private void forwardFast2(ArrayList<Marker> markers, int S[][], int startBin, int endBin) {
		forwardFast2(markers, S, startBin, endBin, 0, markers.size());
	}
	// forward calculation, decreasing map positions
	private void forwardFast2(ArrayList<Marker> markers, int S2[][], int startBin, int endBin, int start, int end)
	{
		for (int mi = start; mi < end; ++mi) {
			Marker m = markers.get(mi);
			int maxS = 0;
			int S2i[] = S2[mi];
			int S2i1[] = S2[mi + 1]; 
			
			for (int b = endBin - startBin; b >= 0; --b) { 
                int s = S2i[b]; 
                if (s > maxS)
                    maxS = s;
                S2i1[b] = maxS;// + m.inside(b + startBin); // add startBin as we subtract it from b 
            }
			int b = startBin;
			while (b <= endBin) {
				S2i1[b - startBin] += m.inside(b);
				b = m.nextInside(b);
			}

		}		
	}

	//positions increasing
	private int solveForward(ArrayList<Marker> markers) {
		return solveForward(markers, getMinBin(markers), getMaxBin(markers));
	}
	private int solveForward(ArrayList<Marker> markers, int startBin, int endBin) {
		//System.err.println(startBin);
		//System.err.println(endBin);
		
		int numMarkers = markers.size();
		int S[][] = new int[numMarkers + 1][endBin - startBin + 1];
		int P[][] = new int[numMarkers + 1][endBin - startBin + 1];

		forward1(markers, S, P, startBin, endBin);

		int maxI = startBin;
		for (int i = startBin; i <= endBin; ++i) { 
            if (S[numMarkers][i - startBin] > S[numMarkers][maxI - startBin])
                maxI = i;
		}
		
		int ret = S[numMarkers][maxI - startBin];
		
		 // + orientation
		for (int mi = numMarkers; mi > 0; --mi) {
			markers.get(mi - 1).pPlus = maxI;
			//System.err.println(markers.get(mi - 1) + "\t" + maxI);
			maxI = P[mi][maxI - startBin];
		}
		return ret;
	}
	private int solveForwardFast(ArrayList<Marker> markers) {
		return solveForwardFast(markers, getMinBin(markers), getMaxBin(markers));
	}
	private int solveForwardFast(ArrayList<Marker> markers, int startBin, int endBin) {
		//System.err.println(startBin);
		//System.err.println(endBin);
		
		int numMarkers = markers.size();
		int S[][] = new int[numMarkers + 1][endBin - startBin + 1];

		forwardFast1(markers, S, startBin, endBin);

		int ret = 0;
		int SnM[] = S[numMarkers];
		for (int b = startBin; b <= endBin; ++b) {
			int s = SnM[b - startBin];
            if (s > ret)
                ret = s;
		}
		
		return ret;
	}

	
	//positions decreasing
	private int solveBackward(ArrayList<Marker> markers) {
		return solveBackward(markers, getMinBin(markers), getMaxBin(markers));
	}
	private int solveBackward(ArrayList<Marker> markers, int startBin, int endBin) {
		//System.err.println(startBin);
		//System.err.println(endBin);
		
		int numMarkers = markers.size();

		int S2[][] = new int[numMarkers + 1][endBin - startBin + 1];
		int P2[][] = new int[numMarkers + 1][endBin - startBin + 1];
		
		forward2(markers, S2, P2, startBin, endBin);
		
		int maxI2 = endBin;
		for (int i = endBin; i >= startBin; --i) { 
            if (S2[numMarkers][i - startBin] > S2[numMarkers][maxI2 - startBin])
                maxI2 = i;
		}
		
		int ret = S2[numMarkers][maxI2 - startBin];
		
		 // - orientation
		for (int mi = numMarkers; mi > 0; --mi) {
			markers.get(mi - 1).pMinus = maxI2;
			//System.err.println(markers.get(mi - 1) + "\t" + maxI2);
			maxI2 = P2[mi][maxI2 - startBin];
		}
		return ret;
	}
	
	public static int[] markerScore(ArrayList<Marker> markers)
	{
		PlaceAndOrientContigs poc = new PlaceAndOrientContigs();
		return poc.solve(markers);
	}
	
	private int[] solve(ArrayList<Marker> markers)
	{
		int scorePlus = solveForward(markers);
		int scoreMinus = solveBackward(markers);
		int ret[] = new int[]{scorePlus, scoreMinus};
		return ret;
	}
	
	public void loadBed(String fn, int chromosome)
	{
		bed = InputData.loadBed(fn, chromosome);
		for (ContigInterval ci : bed) {
			if (!bedHash.containsKey(ci.getContig())) {
				bedHash.put(ci.getContig(), new ArrayList<ContigInterval>());
				contigs.add(ci.getContig());
			}
			bedHash.get(ci.getContig()).add(ci);
		}
		for (String c: contigs) {
			ArrayList<ContigInterval> aci = bedHash.get(c);
			Collections.sort(aci); // Put positions into physical order
			long prev = 0;
			for (ContigInterval ci: aci) {
				if (ci.getStart() > 0 && ci.getStart() > ci.getEnd() || ci.getStart() <= prev) { // check the consistency of the bed...
					System.err.println(c + "\t" + ci.getStart() + "\t" + ci.getEnd());
					System.err.println("Error: erroneous or overlapping intervals in the bed file");
					System.exit(-1);
				}
				prev = ci.getEnd();
			}
		}
		System.err.println(contigs.size() + " contigs loaded from the bed file");
		//TODO: 
	}

	//reverse map
	private void flip(ArrayList<Marker> markers) {
		int max = getMaxBin(markers);
		for (Marker m : markers)
			m.flip(max);
	}

	private long[] chain2OneBase(String length, String orientation, String start, String stop)
	{
		if ("+".equals(orientation)) {
			return new long[]{Long.parseLong(start) + 1, Long.parseLong(stop)};
		}
		else {
			long l = Long.parseLong(length);
			return new long[]{l - Long.parseLong(stop) + 1, l - Long.parseLong(start)};
		}
	}
	
	// first coordinate system to second
	private long mapPosition12(ArrayList<long[]> alignment, long pos1, boolean sameOrientation){
		int low = 0;
		int high = alignment.size() - 1;
		
		while (low <= high) {
			int mid = (low + high) / 2;
			long a[] = alignment.get(mid);
			long mpos = a[0];
			long len = a[2];
			
			if (mpos > pos1)
				high = mid - 1;
			else if (mpos <= pos1 && mpos + len > pos1)
				return a[1] + (sameOrientation ? (pos1 - mpos) : (mpos - pos1));
			else if (mpos < pos1)
				low = mid + 1;
			else
				return -1;
		}
		return -1; 
	}

	// second coordinate system to first
	private long mapPosition21(ArrayList<long[]> alignment, long pos2, boolean sameOrientation){
		int low = 0;
		int high = alignment.size() - 1;
		
		while (low <= high) {
			int mid = (low + high) / 2;
			long a[] = alignment.get(mid);
			long mpos = a[1];
			long len = a[2];
			if (sameOrientation) {
				if (mpos > pos2)
					high = mid - 1;
				else if (mpos <= pos2 && mpos + len > pos2)
					return a[0] + (pos2 - mpos);
				else if (mpos < pos2)
					low = mid + 1;
				else
					return -1;
			} else {
				if (mpos < pos2)
					high = mid - 1;
				else if (mpos >= pos2 && mpos - len < pos2)
					return a[0] + (pos2 - mpos);
				else if (mpos > pos2)
					low = mid + 1;
				else
					return -1;
			}
		}
		return -1; 
	}
	
	//c2 is full haplotype of ofHaplotype?
	private int calculateLiftOverHaplotype(ContigInterval ofHaplotype, int ofHaplotypeOrientation1, ContigInterval c2, ArrayList<long[]> alignment, boolean sameOrientation)
	{
		int score = 0;
		for (int map = 0; map < numMaps; ++map) {
			int numM2 = c2.getMarkers(map).size();

			if (numM2 > 0) {
				ArrayList<Marker> markers = new ArrayList<Marker>();
				markers.addAll(ofHaplotype.getMarkers(map, ofHaplotypeOrientation1));

				int score1 = solveForward(markers);

				ArrayList<Marker> markers2 = c2.getMarkers(map);
				for (Marker m: markers2) {
					long p21 = mapPosition21(alignment, m.getPosition(), sameOrientation);
					if (p21 > 0)
						markers.add(new Marker(m, p21));
				}
				Collections.sort(markers);
				if (ofHaplotypeOrientation1 != 0)
					Collections.reverse(markers);
				
				score += solveForward(markers) - score1;
			}
		}
		return score;
	}	

	private int[] calculateLiftOver(ContigInterval c1, int orientation1, ContigInterval c2, int orientation2, ArrayList<long[]> alignment, boolean sameOrientation, long interval1[], long interval2[])
	{
		int score1 = 0;
		int score2 = 0;
		for (int map = 0; map < numMaps; ++map) {
			int numM1 = c1.getMarkers(map).size();
			int numM2 = c2.getMarkers(map).size();
			if (numM2 > 0 || numM1 > 0){ // markers in one or both contigs...
				int start1 = 0;
				int end1 = numM1 - 1;
				if (orientation1 != 0) {
					end1 = 0;
					start1 = numM1 - 1;
				}
				
				int start2 = 0;
				int end2 = numM2 - 1;
				if (orientation2 != 0) {
					end2 = 0;
					start2 = numM2 - 1;
				}
				
				// c1|----------------|
				//              ||||| (interval 1&2 markers put to am1 and am2)
				//          c2|-------------|
				
				ArrayList<Marker> am1 = new ArrayList<Marker>(); 
				ArrayList<Marker> am2 = new ArrayList<Marker>(); 
				
				while (orientation1 == 0 && end1 >= start1 && c1.getMarker(map, end1).position >= interval1[0]) {
					if (c1.getMarker(map, end1).position <= interval1[1])
						am1.add(c1.getMarker(map, end1));
					--end1;
				}

				while (orientation1 != 0 && end1 <= start1 && c1.getMarker(map, end1).position <= interval1[1]) {
					if (c1.getMarker(map, end1).position >= interval1[0])
						am1.add(c1.getMarker(map, end1));
					++end1;
				}

				while (orientation2 == 0 && end2 >= start2 && c2.getMarker(map, start2).position <= interval2[1]) {
					if (c2.getMarker(map, start2).position >= interval2[0])
						am2.add(c2.getMarker(map, start2));
					++start2;
				}

				while (orientation2 != 0 && end2 <= start2 && c2.getMarker(map, start2).position >= interval2[0]) {
					if (c2.getMarker(map, start2).position <= interval2[1])
						am2.add(c2.getMarker(map, start2));
					--start2;
				}
				
				ArrayList<Marker> test1 = new ArrayList<Marker>();
				if (orientation1 == 0) 
					for (int i = start1; i <= end1; ++i)
						test1.add(c1.getMarker(map, i));
				else
					for (int i = start1; i >= end1; --i)
						test1.add(c1.getMarker(map, i));
				
				//DO MAGIC HERE...
				Collections.reverse(am1);

				ArrayList<Marker> mid1 = new ArrayList<Marker>();
				ArrayList<Marker> mid2 = new ArrayList<Marker>();
				ArrayList<Marker> merge1 = new ArrayList<Marker>();
				ArrayList<Marker> merge2 = new ArrayList<Marker>();
				
				//merge1
				for (Marker m : am1) {
					long p2 = mapPosition12(alignment, m.getPosition(), sameOrientation);
					if (p2 > 0) { // liftover possible
						Marker m2 = new Marker(m, p2);
						mid2.add(m2);
					}
				}
				merge1.addAll(mid2);
				merge1.addAll(am2);
				Collections.sort(merge1);
				if (orientation2 != 0)
					Collections.reverse(merge1);
				
				//merge2
				for (Marker m : am2) {
					long p1 = mapPosition21(alignment, m.getPosition(), sameOrientation);
					if (p1 > 0) { // liftover possible
						Marker m2 = new Marker(m, p1);
						mid1.add(m2);
					}
				}

				merge2.addAll(mid1);
				merge2.addAll(am1);
				Collections.sort(merge2);

				if (orientation1 != 0)
					Collections.reverse(merge2);
				
				//System.err.println(numM2 + "\t" + start2 + "\t" + end2 + "\t" +  numM1 + "\t" + start1 + "\t" + end1 );
				
				ArrayList<Marker> test3 = new ArrayList<Marker>();
				if (orientation2 == 0) 
					for (int i = start2; i <= end2; ++i)
						test3.add(c2.getMarker(map, i));
				else
					for (int i = start2; i >= end2; --i)
						test3.add(c2.getMarker(map, i));

				ArrayList<Marker> test = new ArrayList<Marker>();

				// c1 + c2 
				test.addAll(c1.getMarkers(map, orientation1));
				test.addAll(c2.getMarkers(map, orientation2));
				int s0 = solveForward(test);
			
				score1 -= s0;
				score2 -= s0;

				//System.err.println(solve(test)[0] + "\t(" + test.size() + ")\t");
				
				test.clear();
				test.addAll(test1);
				test.addAll(merge1);
				test.addAll(test3);
						
				// c1 + merge1 + c2
				int s1 = solveForward(test);
				score1 += s1;
				
				//System.err.print(c1.getContig() + "_" + orientation1 + "\t" + c2.getContig() + "_" + orientation2);
				
				//System.err.print("\tscores\t" + solve(test)[0] + "\t(" + test.size() + ")\t");

				test.clear();
				test.addAll(test1);
				test.addAll(merge2);
				test.addAll(test3);
				//System.err.print(solve(test)[0] + "\t(" + test.size() + ")\t");

				// c1 + merge2 + c2
				int s2 = solveForward(test);
				score2 += s2;

				//System.out.println(s0 + "\t" + s1 + "\t" + s2);
				
			}
		}
		return new int[]{score1, score2};
	}
	
	// change alignment from contig1=>contig2 to contig2=>contig1
	private ArrayList<long[]> reverseAlignment(ArrayList<long[]> a, boolean sameOrientation){
		ArrayList<long[]> ret = new ArrayList<long[]>();
		if (sameOrientation) {
			for (long[] aa : a)
				ret.add(new long[]{aa[1],aa[0], aa[2]});
		} else {
			int n = a.size();
			for (int i = n - 1; i>=0; --i) { // reverse and add length for each block
				long[] aa = a.get(i);
				ret.add(new long[]{aa[1] - aa[2] + 1, aa[0] + aa[2] - 1, aa[2]});
			}
		}
		return ret;
	}
	
	public void loadProximity(String fn, int bin, int maxDistance, double scale){
		prox = new Proximity(bin, maxDistance, scale);
		boolean nempty = prox.loadData(fn, bedHash);
		if (!nempty)
			prox = null;
	}
	
	public void loadPaf(String fn, int maxScore, double scalePaf){
 		System.err.println("loading paf...");
		try {
			BufferedReader br = null;
			if (fn.equals("-"))
				br = new BufferedReader(new InputStreamReader(System.in));
			else
				br = new BufferedReader(new FileReader(fn));
			
			ArrayList<ArrayList<String>> rows = new ArrayList<ArrayList<String>>(); 
			
			do {
				ArrayList<String> row = Input.loadTableRow(br, "\t ");
				if (row == null)
					break;

				if (row.size() < 9) {
					System.err.println("Warning: skipping row " + row);
					continue;
				}
				
				if (!bedHash.containsKey(row.get(0)) && bedHash.containsKey(row.get(5))) { // find "extra" contigs or single reads for scaffolding...
					if (rows.size() == 0 || rows.get(0).get(0).equals(row.get(0)))
						rows.add(row);
					else {
						processPAF(rows);
						rows.clear();
						rows.add(row);
					}
				}
				
			} while (true);
			processPAF(rows);
			System.err.println("Scaffolding links:"); 
			for (String key: scaffoldingLink.keySet()) {
				int value = (int)(0.5 + scalePaf * scaffoldingLink.get(key));
				if (value > maxScore) {
					value = maxScore;
				}
				scaffoldingLink.put(key, value);
				System.err.println(key + "\t" +  value);
			}
			
			br.close();
	    } catch (Exception e) {
	     	e.printStackTrace();
			System.err.println("Error in file " + fn);
      	}
	}
	
	//int tmps = 0;
	
	private void helperPAF1(ContigInterval ci1, ContigInterval ci2, ArrayList<String> row1, ArrayList<String> row2) {

		//check that alignments are in right order in the read...
		
		long rstart1 = Long.parseLong(row1.get(2));
		long rend1 = Long.parseLong(row1.get(3));
		
		long rstart2 = Long.parseLong(row2.get(2));
		long rend2 = Long.parseLong(row2.get(3));

		// require non-intersecting intervals...		
		if (Misc.intersectIntervalsLength(rstart1, rend1, rstart2, rend2) > maxIntersect)
			return;
			
		boolean plus = row1.get(4).equals("+");
		boolean r1first = rstart1 < rstart2;
		
		long start1 = Long.parseLong(row1.get(7));
		long end1 = Long.parseLong(row1.get(8));

		long start2 = Long.parseLong(row2.get(7));
		long end2 = Long.parseLong(row2.get(8));
		
		//++ orientation
		long gap1 = ci1.calculateGapLength(0, start1 + 1, end1);
		long gap2 = ci2.calculateGapLength(1, start2 + 1, end2);
		
		long length1 = end1 - start1;
		long length2 = end2 - start2;
		
		if (length1 - gap1 > 0 && length2 - gap2 > 0 ) {
			if (plus ^ !r1first) {
				//System.err.println(ci1.getContig() + "->" + ci2.getContig());
				String key = ci1.toString() + '+' + ci2.toString() + '+';
				scaffoldingLink_tmp.put(key, 1);
				scaffoldingLinkInfo.put(key, row1.get(0));
			} else { // read is not consistent
				//System.err.println("*");
				//System.err.println(row1);
				//System.err.println(ci1);
				//System.err.println(row2);
				//System.err.println(ci2);			
				//System.err.println("*");
			}
		}
		//-- orientation
		gap1 = ci1.calculateGapLength(1, start1 + 1, end1);
		gap2 = ci2.calculateGapLength(0, start2 + 1, end2);

		if (length1 - gap1 > 0 && length2 - gap2 > 0) {
			if  (plus ^ r1first) {
				//System.err.println(ci1.getContig() + "->" + ci2.getContig());
				String key = ci1.toString() + '-' + ci2.toString() + '-';
				scaffoldingLink_tmp.put(key, 1);
				scaffoldingLinkInfo.put(key, row1.get(0));
				
			} else { // read is not consistent
				//System.err.println("*");
				//System.err.println(row1);
				//System.err.println(ci1);
				//System.err.println(row2);
				//qSystem.err.println(ci2);			
				//System.err.println("*");
			}
		}
	}

	private void helperPAF2(ContigInterval ci1, ContigInterval ci2, ArrayList<String> row1, ArrayList<String> row2) {
		//check that alignments are in right order in the read...
		long rstart1 = Long.parseLong(row1.get(2));
		long rend1 = Long.parseLong(row1.get(3));
		long rstart2 = Long.parseLong(row2.get(2));
		long rend2 = Long.parseLong(row2.get(3));

		// require non-intersecting intervals...		
		if (Misc.intersectIntervalsLength(rstart1, rend1, rstart2, rend2) > maxIntersect)
			return;
			
		boolean plus = row1.get(4).equals("+");
		boolean r1first = rstart1 < rstart2;
		
		long start1 = Long.parseLong(row1.get(7));
		long end1 = Long.parseLong(row1.get(8));

		long start2 = Long.parseLong(row2.get(7));
		long end2 = Long.parseLong(row2.get(8));
		
		//+- orientation
		long gap1 = ci1.calculateGapLength(0, start1 + 1, end1);
		long gap2 = ci2.calculateGapLength(0, start2 + 1, end2);
		
		long length1 = end1 - start1;
		long length2 = end2 - start2;
		
		if (length1 - gap1 > 0 && length2 - gap2 > 0) {
		 	if (plus ^ !r1first) {
				//System.err.println(ci1.getContig() + "->" + ci2.getContig());
				String key = ci1.toString() + '+' + ci2.toString() + '-';
				scaffoldingLink_tmp.put(key, 1);
				scaffoldingLinkInfo.put(key, row1.get(0));
			} else { // read is not consistent
				//System.err.println("*");
				//System.err.println(row1);
				//System.err.println(ci1);
				//System.err.println(row2);
				//System.err.println(ci2);
				//System.err.println("*");
			}
		}
		//-+ orientation
		gap1 = ci1.calculateGapLength(1, start1 + 1, end1);
		gap2 = ci2.calculateGapLength(1, start2 + 1, end2);

		if (length1 - gap1 > 0 && length2 - gap2 > 0) {
		 	if (plus ^ r1first) {
				//System.err.println(ci1.getContig() + "->" + ci2.getContig());
				String key = ci1.toString() + '-' + ci2.toString() + '+';
				scaffoldingLink_tmp.put(key, 1);
				scaffoldingLinkInfo.put(key, row1.get(0));
			
			} else { // read is not consistent
				//System.err.println("*");
				//System.err.println(row1);
				//System.err.println(ci1);
				//System.err.println(row2);
				//System.err.println(ci2);
				//System.err.println("*");
			}
		}
	}	

	private void processPAF(ArrayList<ArrayList<String>> rows)
	{
		HashMap<String, ArrayList<ArrayList<String>>> alignmentHash = new HashMap<String, ArrayList<ArrayList<String>>>();
		for (ArrayList<String> row : rows) {
			ArrayList<ArrayList<String>> a = alignmentHash.get(row.get(5));
			if (a == null) {
				a = new ArrayList<ArrayList<String>>();
				alignmentHash.put(row.get(5), a);
			}
			a.add(row);
		}
		ArrayList<String> contigs = new ArrayList<String>();
		
		contigs.addAll(alignmentHash.keySet());
		
		int numContigs = contigs.size();
		
		for (int c1 = 0; c1 < numContigs; ++c1) {
			String contig1 = contigs.get(c1); 
			//for (int c2 = c1 + 1; c2 < numContigs; ++c2) {
			for (int c2 = c1; c2 < numContigs; ++c2) { //allow links between same contig...

				if (c1 == c2 && bedHash.get(contig1).size() <= 1) // no need to go further if only one possible ContigInterval and contig...
					continue;
				
				String contig2 = contigs.get(c2);
				
				
				ArrayList<ArrayList<String>> rows1 = alignmentHash.get(contig1); 
				ArrayList<ArrayList<String>> rows2 = alignmentHash.get(contig2);

				
				for (ArrayList<String> row1 : rows1)
					for (ArrayList<String> row2 : rows2) {
						long rstart1 = Long.parseLong(row1.get(2));
						long rend1   = Long.parseLong(row1.get(3));
						long rstart2 = Long.parseLong(row2.get(2));
						long rend2   = Long.parseLong(row2.get(3));
						
						// require non-intersecting intervals...		
						long il = Misc.intersectIntervalsLength(rstart1, rend1, rstart2, rend2);

						if (il > maxIntersect || il > 0.5 * Math.min(rend1 - rstart1, rend2 - rstart2))
							continue;

						//TODO: Take into account useChainAndPaf in the calculation of bridgeLength and other...
						//TODO: Have to take into account the contig overlap...

						//read  [rs1  re1] [rs2   re2], bridgeLength = rs2-re1
						//read  [rs2  re2] [rs1   re1], bridgeLength = rs1-re2
						long bridgeLength = Math.max(rstart2 - rend1, rstart1 - rend2);
						if (bridgeLength > maxBridge)
							continue;
					
						if (row1.get(4).equals(row2.get(4))) { // same orientation
							for (ContigInterval ci1 : bedHash.get(contig1)) {
								for (ContigInterval ci2 : bedHash.get(contig2)) {
									if (ci1 != ci2) {
										helperPAF1(ci1, ci2, row1, row2);
										helperPAF1(ci2, ci1, row2, row1);
									}
								}
							}
								
						} else { // different orientation
							for (ContigInterval ci1 : bedHash.get(contig1)) {
								for (ContigInterval ci2 : bedHash.get(contig2)) {
									if (ci1 != ci2) {
										helperPAF2(ci1, ci2, row1, row2);
										helperPAF2(ci2, ci1, row2, row1);
									}
								}
							}
						}
					}
			}
		}
		
		for (String key : scaffoldingLink_tmp.keySet()) {
			Integer os = scaffoldingLink.get(key);
			if (os == null)
				scaffoldingLink.put(key, scaffoldingLink_tmp.get(key));
			else
				scaffoldingLink.put(key, os + scaffoldingLink_tmp.get(key));
		}
		scaffoldingLink_tmp.clear();
		
		//if (numContigs > 1)
		//	System.err.println(numContigs + "\t" + (++tmps));
		
		//calculate coverage for each contig (reads spanning each position...)
	}
	
	private void addHaplotype(ContigInterval ofHaplotype, int ofHaplotypeOrientation, ContigInterval haplotype, long ha[], double ascore, ArrayList<long []> alignment, boolean sameOrientation)
	{
		
		int haplotypeScore = (int) (ascore - haplotype.calculateGapLengthHaplotype(ha[0], ha[1]) * cutPenalty + 0.5); //+0.5 = rounding
		if (haplotypeScore >= minHaplotypeAlignmentScore) {
			haplotypeScore += calculateLiftOverHaplotype(ofHaplotype, ofHaplotypeOrientation, haplotype, alignment, sameOrientation);
			if (haplotypeScore > 0) {
				String key = ofHaplotype.toString() + ((ofHaplotypeOrientation == 0) ? '+':'-') + haplotype.toString();
				Integer oscore = chainHaplotypeHash.get(key);
				if (oscore == null || oscore < haplotypeScore) {
					chainHaplotypeHash.put(key, (int) haplotypeScore);
					chainHaplotypeCut.put(key, ha);
				}
	
				//System.err.println(key + "\t" + (ascore * scaleScore - calculateGapLengthHaplotype(ci2, p2[0], p2[1])  * cutPenalty) + "\t" + haplotypeScore1 + "\thaplotype");
			}
		}
	}

	// add score between contigs ci1 and ci2 in these orientation based on alignment chain...
	private void addChainLink(ContigInterval ci1, int orientation1, ContigInterval ci2, int orientation2, long ai1[], long ai2[], double ascore, ArrayList<long []> alignment, boolean sameOrientation)
	{
		long l1 = ci1.calculateGapLength(orientation1, ai1[0], ai1[1]);
		long l2 = ci2.calculateGapLength(1 - orientation2, ai2[0], ai2[1]);
		
		double score = ascore;
		if ((orientation1 == orientation2) != sameOrientation)
			score *= orientationPenalty;
		
		score -= (l1 + l2) * cutPenalty;
		
		if (score >= minLinkAlignmentScore) { // is this needed?
			int los[] = calculateLiftOver(ci1, orientation1, ci2, orientation2, alignment, sameOrientation, ai1, ai2);
			int max1 = 0;
			if (los[1] > los[0])
				max1 = 1;

			score += los[max1];
			
			if (score >= 0.5) {
				String key = ci1.toString() + ((orientation1 == 0) ? '+':'-') + ci2.toString() + ((orientation2 == 0) ? '+':'-');
				//System.err.println(key + "\t" + score1);
				Integer oldScore = chainLinkHash.get(key);
				score = (int)(score + 0.5); // round up...
				if (oldScore == null || oldScore < score) {
					chainLinkHash.put(key, (int) score);
					
					//if ((orientation1 == orientation2) != sameOrientation) { // orientation do not match... make a palindrome...
					//	long start1 = (orientation1 == 0) ? ai1[1] : ai1[0];
					//	chainLinkHashCut1.put(key, start1);
					//	chainLinkHashCut1a.put(key, start1);
					//	long end2 = (orientation2 == 0) ? ai2[0] : ai2[1];
					//	chainLinkHashCut2.put(key, end2);
					//	chainLinkHashCut2a.put(key, end2);
					//} else 
					{
						if (orientation1 == 0) { // + orientation
							if (max1 == 1) {
								chainLinkHashCut1.put(key, ai1[1]); // take whole alignment
								chainLinkHashCut1a.put(key, ai1[0] - 1);
							}
							else {
								chainLinkHashCut1.put(key, ai1[0] - 1); // cut just before alignment
								chainLinkHashCut1a.put(key, ai1[1]);
							}
						}
						else { // - orientation
							if (max1 == 1) {
								chainLinkHashCut1.put(key, ai1[0]); // whole
								chainLinkHashCut1a.put(key, ai1[1] + 1);
							}
							else {
								chainLinkHashCut1.put(key, ai1[1] + 1); // after
								chainLinkHashCut1a.put(key, ai1[0]);
							}
						}
	 					if (orientation2 == 0) { // + orientation
							if (max1 == 1) {
								chainLinkHashCut2.put(key, ai2[1] + 1); // after
								chainLinkHashCut2a.put(key, ai2[0]);
							}
							else {
								chainLinkHashCut2.put(key, ai2[0]); // whole
								chainLinkHashCut2a.put(key, ai2[1] + 1);
							}
						}
						else { // - orientation
							if (max1 == 1) {
								chainLinkHashCut2.put(key, ai2[0] - 1); // before
								chainLinkHashCut2a.put(key, ai2[1]);
							}
							else {
								chainLinkHashCut2.put(key, ai2[1]); // whole
								chainLinkHashCut2a.put(key, ai2[0] - 1);
							}
						}
					}
				}
			}
		}
	}
	// Find largest index i, where array[i] <= value, array sorted ascending order
	private int binarySearchR(ArrayList<Integer> array, int value){
		int right = array.size() - 1;
		if (right < 0)
			return -1;
		int left = 0;
		
		while (right > left + 1) {
			int mid = (left + right) / 2;
			int am = array.get(mid); 
			if (am > value)
				right = mid - 1;
			else
				left = mid + ((am == value) ? 0 : 1);
		}
		if (array.get(right) <= value)
			return right;
		if (array.get(left) <= value)
			return left;
		return left - 1;
	}
	
	// Find smallest index i, where array[i] >= value, array sorted ascending order
	private int binarySearchL(ArrayList<Integer> array, int value){
		int right = array.size() - 1;
		if (right < 0)
			return 0;
		int left = 0;
		
		while (right > left + 1) {
			int mid = (left + right) / 2;
			int am = array.get(mid); 
			if (am < value)
				left = mid + 1;
			else
				right = mid - ((am == value) ? 0 : 1);
		}
		if (array.get(left) >= value)
			return left;
		if (array.get(right) >= value)
			return right;
		return right + 1;
	}
	
	//calculate how many markers in array between start and stop
	// array sorted
	private int calculateNumMarkers(ArrayList<Integer> array, int start, int end)
	{
		return  (binarySearchR(array, end) - binarySearchL(array, start) + 1);
	}
	
	private long[] getStartEnd1(ArrayList<long[]> alignment)
	{
		long first[] = alignment.get(0);
		long last[] = alignment.get(alignment.size() - 1);
		return new long[]{first[0], last[0] + last[2] - 1};
	}

	private long[] getStartEnd2(ArrayList<long[]> alignment, boolean sameOrientation)
	{
		long first[] = alignment.get(0);
		long last[] = alignment.get(alignment.size() - 1);
		if (sameOrientation)
			return new long[]{first[1], last[1] + last[2] - 1};
		else
			return new long[]{last[1] - last[2] + 1, first[1]}; 
	}
	
	//
	private double trimChain(ArrayList<long[]> alignment, boolean sameOrientation, ContigInterval ci1, ContigInterval ci2, ArrayList<long[]> ret) {
		if (alignment.size() == 0)
			return 0.0;

		long ap1[] = getStartEnd1(alignment); 
		long ap2[] = getStartEnd2(alignment, sameOrientation);

		// if contig intervals and alignment do not intersect, we are done here...
		if (!Misc.intersectIntervals(ap1[0], ap1[1], ci1.getMinStart(), ci1.getMaxEnd()) || !Misc.intersectIntervals(ap2[0], ap2[1], ci2.getMinStart(), ci2.getMaxEnd()))
			return 0.0;
		long trimLength1[] = new long[2];
		
		if (ap1[0] < ci1.getMinStart()) // alignment spans too far
			trimLength1[0] = ci1.getMinStart() - ap1[0];
		if (ap1[1] > ci1.getMaxEnd()) // alignment spans too far
			trimLength1[1] = ap1[1] - ci1.getMaxEnd();

		long trimLength2[] = new long[2];
		
		
		if (ap2[0] < ci2.getMinStart()) // alignment spans too far 
			trimLength2[0] = ci2.getMinStart() - ap2[0];  
		if (ap2[1] > ci2.getMaxEnd()) // alignment spans too far
			trimLength2[1] = ap2[1] - ci2.getMaxEnd();

		if (!sameOrientation) {
			long tmp =  trimLength2[0];
			trimLength2[0] = trimLength2[1];
			trimLength2[1] = tmp; 
		}
		
		long tl1 = trimLength1[0] + trimLength2[0];
		long tl2 = trimLength1[1] + trimLength2[1];
		if (tl1 + tl2 > 0) {
			//System.err.println("Alignment of " + ci1.getContig() + " and " + ci2.getContig() + " needs trimming of " + (tl1 + tl2));

			long first[] = alignment.get(0);
			long last[] = alignment.get(alignment.size() - 1);
			//System.err.println(first[0] + "\t" + first[1] + "\t" + last[0] + "\t" + last[1] + "\t" + last[2] + "\t" + ci1 + "\t" + ci2);
			
			long alignL = 0;
			for (long[] a : alignment) 
				alignL += a[2];
			
			int t1 = 0;
			long wt1 = 0;
			long totalW1 = 0;
			if (tl1 > 0) {
				for (long[] a : alignment) {
					long w = a[2];
					long dt1 = a[0] - ap1[0];
					long dt2 = ((sameOrientation) ? a[1] - ap2[0] : ap2[1] - a[1]);
					if (dt1 + w - 1 >= trimLength1[0] && dt2 + w - 1 >= trimLength2[0]) {
						if (dt1 >= trimLength1[0] && dt2 >= trimLength2[0]) {
							wt1 = 0;
							//System.err.print("*full ");
						} else {
							wt1 = Math.max(trimLength1[0] - dt1, trimLength2[0] - dt2);
							//System.err.print("*"+ wt1 + " of " +  w + " ");
						}
						totalW1 += wt1;
						//System.err.println("trimming " + t1 + " blocks (of " + alignment.size() +  ") with length " + dt1 + "\t" + dt2 + "\t" + ((alignL - totalW1) / (double) alignL));
						
						break;
					} else
						totalW1 += w;
					++t1;
				}
			}

			int t2 = 0;
			long wt2 = 0;
			long totalW2 = 0;
			if (tl2 > 0) {
				for (int ai = alignment.size() - 1; ai >= 0; --ai) {
					long[] a = alignment.get(ai);
					long w = a[2];

					long p1 = a[0] + w - 1;
					long p2 = ((sameOrientation) ? a[1] + w - 1 : a[1] - w + 1);
					
					long dt1 = ap1[1] - p1;
					long dt2 = ((sameOrientation) ? ap2[1] - p2 : p2 - ap2[0]);

					if (dt1 + w - 1 >= trimLength1[1] && dt2 + w - 1 >= trimLength2[1]) {
						if (dt1 >= trimLength1[1] && dt2 >= trimLength2[1]) {
							wt2 = 0;
							//System.err.print("full ");
						} else {
							wt2 = Math.max(trimLength1[1] - dt1, trimLength2[1] - dt2);
							//System.err.print(wt2 + " of " +  w + " ");
						}
						totalW2 += wt2;
						//System.err.println("trimming " + t2 + " blocks (of " + alignment.size() +  ") with length " + dt1 + "\t" + dt2 + "\t" + ((alignL - totalW2) / (double) alignL));
						
						break;
					} else
						totalW2 += w;
					++t2;
				}
			}
			if (totalW1 + totalW2 >= alignL) {
				//System.err.println("all trimmed of a chain");
				return 0.0; // all trimmed off
			} else {
				ret.clear();
				int n = alignment.size();
				for (int i = t1; i < n - t2; ++i) {
					long a[] = alignment.get(i);
					if (i == t1 && wt1 > 0) {
						ret.add(new long[]{a[0] + wt1, a[1] + ((sameOrientation) ? + wt1 : -wt1), a[2] - wt1});
					} else if (i == n - t2 - 1 && wt2 > 0) {
						ret.add(new long[]{a[0], a[1], a[2] - wt2});
					} else
						ret.add(a);
				}
				return (alignL - totalW1 - totalW2) / (double) alignL;
			}
			
		} else {
			//System.err.println("Alignment of " + ci1.getContig() + " and " + ci2.getContig() + " does not need trimming");
			return 1.0;
		}
		//total length of aligning bases...
	}
	
	public void loadChain(String fn){
		int numAlignments = 0;
		/*
		//calculate data structure to get the number of marker in each interval in O(log(n)) time
		
		HashMap<String, ArrayList<Integer>> markerHash = new HashMap<String, ArrayList<Integer>>();
		for (ContigInterval ci : intervalsInChr) {
			ArrayList<Integer> markers = new ArrayList<Integer>();
			markerHash.put(ci.getContig(), markers);
			for (int map = 0; map < numMaps; ++map)
				for (Marker m : ci.getMarkers(map))
					markers.add(m.getPosition());
			Collections.sort(markers);
		}*/

 		System.err.println("loading alignment chain...");
		try {
			BufferedReader br = null;
			if (fn.equals("-"))
				br = new BufferedReader(new InputStreamReader(System.in));
			else
				br = new BufferedReader(new FileReader(fn));
			do {
				ArrayList<String> row = Input.loadTableRow(br, "\t ");
				if (row == null)
					break;
				boolean skip = true;
				if (row.size() >= 12 && "chain".equals(row.get(0))) {
					String contig1 = row.get(2);
					String contig2 = row.get(7);
					
					if (bedHash.containsKey(contig1) && bedHash.containsKey(contig2)) {
						skip = false;
						long p1_[] = chain2OneBase(row.get(3), row.get(4), row.get(5), row.get(6));
						long p2_[] = chain2OneBase(row.get(8), row.get(9), row.get(10), row.get(11));
						if ("-".equals(row.get(4))) {
							System.err.println("Error: only ++, and +- orientation allowed in the chain file");
							continue;
						}
						double ascore_ = scaleScore * Double.parseDouble(row.get(1)); // alignment score
	
						//int numMarkers1 = calculateNumMarkers(markerHash.get(contig1), p1[0], p1[1]);
						//int numMarkers2 = calculateNumMarkers(markerHash.get(contig2), p2[0], p2[1]);
						//System.err.println(contig1 + "\t" + p1[0] + " " + p1[1] + " " + numMarkers1);
						//System.err.println(contig2 + "\t" + p2[0] + " " + p2[1] + " " + numMarkers2);
						
						boolean sameOrientation = row.get(4).equals(row.get(9));

						// load alignment for liftover...
						ArrayList<long[]> alignment_ = new ArrayList<long[]>();
						ArrayList<String> row2 = null;
						long v1 = p1_[0];
						long v2 = ((sameOrientation) ? p2_[0] : p2_[1]);
						int v3 = 0;
						do {
							row2 = Input.loadTableRow(br, "\t ");
							
							v3 = Integer.parseInt(row2.get(0));
							alignment_.add(new long[]{v1, v2, v3});
							v1 += v3;
							v2 += ((sameOrientation) ? v3 : -v3);
							if (row2.size() >= 3) {
								v1 += Integer.parseInt(row2.get(1));
								v2 += ((sameOrientation) ? Integer.parseInt(row2.get(2)) : -Integer.parseInt(row2.get(2)));
							}
						} while (row2 != null && row2.size() >= 3);
						
						ArrayList<long[]> revAlignment_ = reverseAlignment(alignment_, sameOrientation);

						//System.err.println(row);
						//for (int[] a : alignment)
						//	System.err.println(a[0] + "\t" +  a[1] + "\t" + a[2]);

						//for (int i = 1; i<=100; i+=1)
						//	System.err.println(mapPosition12(alignment, i, sameOrientation));
						//System.err.println("*");
						//for (int i = 1; i<=100; i+=1)
						//	System.err.println(mapPosition21(alignment, i, sameOrientation));
							
						//if (true)
						//	System.exit(-1);
					
						//check here
						
						for (ContigInterval ci1 : bedHash.get(contig1))
							for (ContigInterval ci2 : bedHash.get(contig2)) {
								
								ArrayList<long[]> alignment = new ArrayList<long[]>();
								
								double t = trimChain(alignment_, sameOrientation, ci1, ci2, alignment);
								if (t == 0.0)
									continue;
								
								ArrayList<long[]> revAlignment = new ArrayList<long[]>();
								long p1[]; 
								long p2[];
								double ascore = t * ascore_;
								
								if (t < 1.0) {
									double t2 = trimChain(revAlignment_, sameOrientation, ci2, ci1, revAlignment);
									if (t != t2) {
										System.err.println("Error in parsing alignment chain:" + t +"!="+  t2);
										System.exit(-1);
									}
									
									p1 = getStartEnd1(alignment);
									p2 = getStartEnd2(alignment, sameOrientation);
									if (p2[0] > p2[1]) {
										long tmp = p2[0];
										p2[0] = p2[1];
										p2[1] = tmp; 
									}
									
								} else {
									alignment = alignment_;
									revAlignment = revAlignment_;
									p1 = p1_;
									p2 = p2_;
								}
								for (int orientation1 = 0; orientation1 < 2; ++orientation1) { 
									//calculate haplotype scores...
									//ci2 is the haplotype of ci1?
									addHaplotype(ci1, orientation1, ci2, p2, ascore, alignment, sameOrientation);
									//ci1 is the haplotype of ci2?
									//use orientation1 for ci2 as well... (1 - orientation is not really needed)
									addHaplotype(ci2, 1 - orientation1, ci1, p1, ascore, revAlignment, sameOrientation);
									
									for (int orientation2 = 0; orientation2 < 2; ++orientation2) {  // four combinations, ++/+-/-+/--
										//add partial haplotype link scores
										addChainLink(ci1, orientation1, ci2, orientation2, p1, p2, ascore, alignment, sameOrientation);
										//add partial haplotype link scores (1 - orientation is not really needed)
										addChainLink(ci2, 1 - orientation2, ci1, 1 - orientation1, p2, p1, ascore, revAlignment, sameOrientation);
											//System.err.println("hiphei");										
									}
								}
							}
							
					++numAlignments;
					}
				} else {
					if (skip && (row.size() == 3 || row.size() == 1))
						;
					else
						System.err.println("Warning: skipping " + row);
				}
			} while (true);
			br.close();
	    } catch (Exception e) {
	     	e.printStackTrace();
			System.err.println("Error in file " + fn);
      	}
		System.err.println("loading " +  numAlignments + " alignments");
		
		System.err.println("Alignment links:");
		for (String key:chainLinkHash.keySet()) 
			System.err.println(key + "\t" + chainLinkHash.get(key));		
	}

	public void addMap(String fn, boolean flip, boolean nochromosome, boolean nointervals, int chromosome)
	{
		
		ArrayList<Marker> markers = InputData.loadMap(fn, nochromosome, nointervals, compressMap);
		
		if (flip) {
			flip(markers);
		}

		for (Marker m : markers) 
			if (chromosome <= 0 || m.getChromosome() == chromosome) { 
			boolean inside = false;
			if (bedHash.containsKey(m.contig))
				for (ContigInterval ci : bedHash.get(m.getContig())) {
					if (ci.inside(m.position)) {
						inside = true;
						ci.addMarker(m, numMaps);
						break;
					}
				}
			if (!inside) {
//				System.err.println(m);
//				System.err.println("Warning: marker position not listed in the bed file");
				//System.exit(-1);
			}
		}		
		System.err.println("#map = " + fn);
		System.err.println("contig\tstart_pos\tend_pos\tscore+orientation\tscore-orientation");
		for (ContigInterval ci : bed) {
			if (ci.markers.size() > numMaps) {
				Collections.sort(ci.markers.get(numMaps)); // Put positions into physical order

				int scores[] = solve(ci.markers.get(numMaps));	// Solve marker positions by dynamic programming...
				//ci.scores.set(numMaps, scores);
				
				System.err.println(ci + "\t" + scores[0] + "\t" + scores[1]);
			}
		}
		++numMaps;

		int maxChromosome = 0;
		for (ContigInterval ci : bed)
			maxChromosome = Math.max(maxChromosome, ci.getChromosome());

		int minChromosome = maxChromosome; 
		for (ContigInterval ci : bed)
			if (ci.getChromosome() > 0)
				minChromosome = Math.min(minChromosome, ci.getChromosome());
		
		if (minChromosome != maxChromosome) {
			System.err.println("Error: Only one chromosome allowed in the input maps!");
			System.exit(-1);
		}
		
		intervalsInChr = new ArrayList<ContigInterval>();

		for (ContigInterval ci : bed) {
			if (ci.getChromosome() == maxChromosome || ci.getChromosome() == 0)
				intervalsInChr.add(ci);
		}
		
	}
	
	public void combineMaps(boolean findOrientation, boolean randomOrder, boolean finalImprove)
	{
		//(remove) check intervals without any markers...
		ArrayList<ContigInterval> contigsWithMarkers = new ArrayList<ContigInterval>(); 

		ArrayList<ContigInterval> contigsWithoutMarkers = new ArrayList<ContigInterval>(); 

		for (ContigInterval ci : intervalsInChr) {
			if (ci.markers.size() > 0)
				contigsWithMarkers.add(ci);
			else {
				contigsWithoutMarkers.add(ci);
				System.err.println("Contig " + ci + " does not have any markers");
			}
		}
		intervalsInChr = contigsWithMarkers;
		if (keepEmptyIntervals)
			for (ContigInterval ci : contigsWithoutMarkers) {
				intervalsInChr.add(ci);
			}
		
		
		
		if (randomOrder) { // destroy initial solution... 
			Collections.shuffle(intervalsInChr);
			for (ContigInterval c : intervalsInChr)
				if (Math.random() < 0.5)
					c.flipOrientation();
		} else { // find contig order by ordering contigs on the most abundant marker

			ArrayList<Marker> allMarkers = getMarkers(intervalsInChr, null, 0);
			
			int maxBin = getMaxBin(allMarkers);
			int minBin = getMinBin(allMarkers);
			
			//int n = intervalsInChr.size();
			ArrayList<Double> sortValues = new ArrayList<Double>();
			
			for (ContigInterval ci : intervalsInChr) {
				ArrayList<Marker> markersInFirstMap = ((ci.markers.size() > 0) ? ci.markers.get(0) : new ArrayList<Marker>());
				int s[] = solve(markersInFirstMap);
				if (s[1] > s[0] || (s[1] == s[0] && Math.random() < 0.5))  
					ci.flipOrientation();
				
				HashMap<Integer, Integer> hm = new HashMap<Integer, Integer>();
				
				for (Marker m : markersInFirstMap) {
					int numI = m.intervals.length;
					for (int i = 0; i < numI; i+=2)
						for (int p = m.intervals[i]; p <= m.intervals[i + 1]; ++p) {
							Integer oldValue = hm.get(p);
							if (oldValue == null)
								hm.put(p, 1);
							else
								hm.put(p, oldValue + 1);
						}
				}
					
				int maxCount = -1;
				int bin = 0;
				int numMaxCount = 0;
				for (int p : hm.keySet()) {
					int count = hm.get(p); 
					if (count > maxCount) {
						numMaxCount = 1;
						bin = p;
						maxCount = count;
					} else if (count == maxCount) {
						++numMaxCount;
						if (Math.random() * numMaxCount < 1.0) // randomly pick one of multiple bins with maximum count
							bin = p;
					}
				}
				if (maxCount > 0) // there was at least one marker
					sortValues.add(bin + Math.random()); //bin + rand[0,1]
				else // no markers
					sortValues.add(Math.random() * (maxBin - minBin + 1) + minBin);
			}
		    Misc.ArrayIndexComparator<Double> comparator = new Misc.ArrayIndexComparator<Double>(sortValues);
		    Integer[] indexes = comparator.createIndexArray();
		    Arrays.sort(indexes, comparator);
		    
			ArrayList<ContigInterval> cis_new = new ArrayList<ContigInterval>(); 
		    for (int i : indexes)
		    	cis_new.add(intervalsInChr.get(i));
		    intervalsInChr = cis_new;
		}
		if (findOrientation && numMaps > 1) {
			String orientation = "+";
			String scores = "";

			//disable chains
			HashMap<String, Integer> oldChainLinkHash = chainLinkHash;
			chainLinkHash = new HashMap<String, Integer>();
			
			//use only first map
			int oldNumMaps = numMaps;
			numMaps = 1;
			
			ArrayList<ContigInterval> cis_new = new ArrayList<ContigInterval>(); 

			for (ContigInterval ci : intervalsInChr) {
				if (ci.markers.size() > 0 && ci.markers.get(0).size() > 0)
					cis_new.add(ci);
			}
			
			//Collections.shuffle(cis_new); // random order as init (could be based on most common marker in each contig...)
			
			for (int map = 1; map < oldNumMaps; ++map) {
				int oldNumRuns = numRuns;
				if (map == 1)
					numRuns = Math.min(2, oldNumRuns);
				else
					numRuns = 1; // this is probably enough...
				improveAnchoring(cis_new, false, true);
				numRuns = oldNumRuns;
						
				numMaps = map + 1;				

				int score1 = calculateScore(cis_new, numMaps, false);

				flip(getMarkers(intervalsInChr, null, map));
				
				int score2 = calculateScore(cis_new, numMaps, false);
				
				if (score1 >= score2) {
					flip(getMarkers(intervalsInChr, null, map));
					orientation += " +";
				} else
					orientation += " -";

				scores = scores + " " + (score1 - score2);


				//add new markers to cis_new to random places...

				ArrayList<ContigInterval> cis_new2 = new ArrayList<ContigInterval>(); 
				out: for (ContigInterval ci : intervalsInChr) {
					if (ci.markers.size() > map && ci.markers.get(map).size() > 0) {
						for (int m = 0; m < map; ++m)
							if (ci.markers.get(m).size() > 0)
								continue out;
						cis_new2.add(ci);
					}
				}
				//Collections.shuffle(cis_new2); no need to shuffle...
				cis_new = Misc.randomMerge(cis_new, cis_new2);
			}
			//System.out.println(intervalsInChr.size() + "\t" + cis_new.size() );
			
			//intervalsInChr.clear();
			//intervalsInChr.addAll(cis_new);
			
			//finally add contigs without markers 
			for (ContigInterval ci : intervalsInChr)
				if (ci.markers.size() == 0)
					cis_new.add(ci);
			
			assert(intervalsInChr.size() == cis_new.size());

			
			intervalsInChr = cis_new;
			
			//enable chains
			chainLinkHash = oldChainLinkHash;

			System.out.println("#found orientation=" + orientation + " ( support " + scores + " )");
		}
		if (finalImprove)
			improveAnchoring();
	}
	
	private ArrayList<Marker> getMarkers(ArrayList<ContigInterval> cis, ContigInterval excluded, int map)
	{
		ArrayList<Marker> ret = new ArrayList<Marker>();
		for (ContigInterval ci : cis) 
			if (ci != excluded) 
				ret.addAll(ci.getMarkers(map));
		return ret;
	}
	
	private String calculateKey(ContigInterval prev, boolean orientationPrev, ContigInterval ci, boolean orientationCi){
		return prev.toString() + ((orientationPrev) ? '+':'-') + ci.toString() + ((orientationCi) ? '+':'-');
	}
	
	private String calculateKey(ContigInterval prev, ContigInterval ci){
		return calculateKey(prev, (prev.getOrientation() >= 0), ci, (ci.getOrientation() >= 0));
	}

	//chainScore with possible flips for prev and ci
	private int calculateChainScore(ContigInterval prev, boolean flipPrev, ContigInterval ci, boolean flipCi) {
		String key = calculateKey(prev, (prev.getOrientation() >= 0) ^ flipPrev, ci, (ci.getOrientation() >= 0) ^ flipCi);
		return calculateChainScore(key);
	}
	private int calculateChainScore(ContigInterval prev, ContigInterval ci) {
		String key = calculateKey(prev, ci);
		return calculateChainScore(key);
	}

	private int calculateChainScore(String key) {
		if (chainLinkHash.containsKey(key)) {
			int ret = chainLinkHash.get(key);
			if (useChainAndPaf && scaffoldingLink.containsKey(key))
				ret += scaffoldingLink.get(key);
			return ret; 
		}
		if (scaffoldingLink.containsKey(key))
			return scaffoldingLink.get(key);
		return 0;
	}

	// haplotype is haplotype of ofHaplotype?
	private int calculateChainScoreHaplotype(ContigInterval ofHaplotype, ContigInterval haplotype) {
		String key = calculateKey(ofHaplotype, haplotype);
		key = key.substring(0, key.length() - 1); // last orientation 
		
		if (chainHaplotypeHash.containsKey(key))
			return chainHaplotypeHash.get(key);
		return 0;
	}

	private long[] calculateChainScoreHaplotypeCut(ContigInterval ofHaplotype, ContigInterval haplotype) {
		String key = calculateKey(ofHaplotype, haplotype);
		key = key.substring(0, key.length() - 1); // last orientation 
		
		return chainHaplotypeCut.get(key);
	}
	
	private int calculateChainScore(ArrayList<ContigInterval> cis)
	{
		int score = 0;
		ContigInterval prev = null;
		for (ContigInterval ci : cis) {
			if (prev != null)
				score += calculateChainScore(prev, ci);
			prev = ci;
		}
		return score;
	}

	private int calculateProximityScore(ArrayList<ContigInterval> cis)
	{
		if (prox != null)
			return prox.score(cis);
		return 0;
	}
	
	private int evaluateScore(ArrayList<ContigInterval> eval, boolean improve)
	{
		ArrayList<ContigInterval> eval2 = new ArrayList<ContigInterval>(); 
		for (ContigInterval e : eval) {
			boolean found = false;
			for (ContigInterval ci : intervalsInChr) {
				if (ci.equals(e)) {
					eval2.add(ci);
					ci.setOrientation(e.getOrientation());
					found = true;
					break;
				}
			}
			if (!found) {
				System.err.println("Contig " + e +  " not found");
			}
		}

		if (improve) {
			int score = calculateScore(eval2, false) + calculateNonMarkerScore(eval2);
			//System.out.print("#initial score " +  score + " ");
			improveAnchoring(eval2, true, true);
		} else {
			ArrayList<Double> orderSupport = calculateSupport(eval2);
			int score = calculateScore(eval2, false) + calculateNonMarkerScore(eval2);
			System.out.println("#final score " +  score);
			printAnchoring(eval2, orderSupport);
			findLikelyAssemblyErrors(eval2);
		}
		//change bed to the evaluated order...
		intervalsInChr.clear();
		intervalsInChr.addAll(eval2);
		
		return calculateScore(eval2, false) + calculateNonMarkerScore(eval2);
/*		
		int score1 = calculateScore(eval2, true);
		for (int map = 0; map < numMaps; ++map)
			flip(getMarkers(intervalsInChr, null, map));
		int score2 = calculateScore(eval2, true);
		
		return Math.max(score1, score2);*/ 
	}


	private int calculateScoreFast(ArrayList<ContigInterval> cis)
	{
		int score = 0;
		for (int map = 0; map < numMaps; ++map) {				
			int s = solveForwardFast(getMarkers(cis, null, map));
			score += s;
		}
		return score;
	}


	private int calculateScore(ArrayList<ContigInterval> cis, boolean verbose)
	{
		return calculateScore(cis, numMaps, verbose);
	}
	
	private int calculateScore(ArrayList<ContigInterval> cis, int numMaps, boolean verbose)
	{
		int score = 0;
		for (int map = 0; map < numMaps; ++map) {				
			int s = solveForward(getMarkers(cis, null, map));
			if (verbose)
				System.err.println("map " + map + " score " + s);
			score += s;
		}
		return score;
	}

	private int getMaxBin(ArrayList<Marker> markers) {
		int maxBin = Integer.MIN_VALUE;
		for (Marker m : markers)
			maxBin = Math.max(maxBin, m.maxBin());
		if (maxBin != Integer.MIN_VALUE)
			return maxBin;
		else 
			return 0;
	}

	private int getMinBin(ArrayList<Marker> markers) {
		int minBin = Integer.MAX_VALUE;
		for (Marker m : markers)
			minBin = Math.min(minBin, m.minBin());
		if (minBin != Integer.MAX_VALUE)
			return minBin;
		else 
			return 0;
	}

	//TODO: Remove calculateChainScore(cis) and calculateChainScore(c1, c2) add calculateNonMarkerScore(...) 
	private ArrayList<Double> calculateSupport(ArrayList<ContigInterval> cis) {
		
		//for (int iteration = 0; iteration < 1; ++iteration) { // add here "iteration < 2" if possible orientation errors should be corrected...
			
		
		
		//ArrayList<Integer> score1 = new ArrayList<Integer>(); 
		ArrayList<Integer> score2 = new ArrayList<Integer>(); 

		ArrayList<Integer> order1 = new ArrayList<Integer>(); 
		ArrayList<Integer> order2 = new ArrayList<Integer>(); 
		
		int rank = 0;
		for (ContigInterval ci : cis) {
			ci.setRank(rank++);
			score2.add(0);
			order1.add(0);
			order2.add(0);
		}
		
		for (int map = 0; map < numMaps; ++map) {
			ArrayList<Marker> markers = getMarkers(cis, null, map);

			int endBin = getMaxBin(markers);
			int startBin = getMinBin(markers);

			int numMarkers = markers.size();
			int S[][] = new int[numMarkers + 1][endBin - startBin + 1];
			int P[][] = new int[numMarkers + 1][endBin - startBin + 1];
		
			int Sb[][] = new int[numMarkers + 1][endBin - startBin + 1];
			//int Pb[][] = new int[numMarkers + 1][endBin - startBin + 1];

			int Sc[][] = new int[numMarkers + 1][endBin - startBin + 1];
			int Pc[][] = new int[numMarkers + 1][endBin - startBin + 1];

			//simple way to calculate (about correct) backward table...
			Collections.reverse(markers);
			forwardFast2(markers, Sb, startBin, endBin);
			Collections.reverse(markers);
			
			forward1(markers, S, P, startBin, endBin);
			int maxScore = -1;
			{
				int maxI = startBin;
				for (int i = startBin; i <= endBin; ++i) 
		            if (S[numMarkers][i - startBin] > maxScore) {
		            	maxScore = S[numMarkers][i - startBin];
		                maxI = i;
		            }
				
				for (int mi = numMarkers; mi > 0; --mi) {
					markers.get(mi - 1).pPlus = maxI;
					//System.err.println(markers.get(mi - 1) + "\t" + maxI);
					maxI = P[mi][maxI - startBin];
				}
			}
			
			ArrayList<ArrayList<Integer>> markerIndex = new ArrayList<ArrayList<Integer>>();
			for (ContigInterval ci : cis)
				markerIndex.add(new ArrayList<Integer>());

			for (int mi = 0; mi < numMarkers; ++mi) {
				Marker m = markers.get(mi);
				rank = m.ci.getRank();
				markerIndex.get(rank).add(mi);
			}
			for (ContigInterval ci : cis) {
				rank = ci.getRank();
				
				ArrayList<Marker> cMarkers = new ArrayList<Marker>();
		
				for (int i : markerIndex.get(rank))
					cMarkers.add(markers.get(i));
				
				int maxB = Integer.MIN_VALUE;
				int minB = Integer.MAX_VALUE;
				
				int numCMarkers = cMarkers.size(); 
				
				if (numCMarkers > 0) {
					for (Marker m : cMarkers) {
						maxB = Math.max(maxB, m.pPlus);
						minB = Math.min(minB, m.pPlus);
					}
					if (maxB > minB)
						order1.set(rank, order1.get(rank) + 1);
					
					Collections.reverse(cMarkers);
					int start = markerIndex.get(rank).get(0);
					int end = markerIndex.get(rank).get(numCMarkers - 1) + 1;
					
					for (int b = startBin; b <= endBin; ++b) {
						Sc[0][b - startBin] = S[start][b - startBin];
						Pc[0][b - startBin] = P[start][b - startBin];
					}
					
					forward1(cMarkers, Sc, Pc, startBin, endBin);

			
					int Scc[] = Sc[cMarkers.size()];
					for (int b = startBin + 1; b <= endBin; ++b) { // make table non-decreasing
						if (Scc[b - startBin] < Scc[b - startBin - 1])
							Scc[b - startBin] = Scc[b - startBin - 1];
					}
					int maxI = startBin;
					int max = -1; 
					for (int b = startBin; b <= endBin; ++b) {
						int s = Scc[b - startBin] + Sb[numMarkers - end][b - startBin]; 
						if (s > max) {
							max = s;
							maxI = b;
						}
					}
					score2.set(rank, score2.get(rank) + max);
					{
						for (int mi = numCMarkers; mi > 0; --mi) {
							cMarkers.get(mi - 1).pPlus = maxI;
							//System.err.println(markers.get(mi - 1) + "\t" + maxI);
							maxI = Pc[mi][maxI - startBin];
						}
					}
					maxB = Integer.MIN_VALUE;
					minB = Integer.MAX_VALUE;
					for (Marker m : cMarkers) {
						maxB = Math.max(maxB, m.pPlus);
						minB = Math.min(minB, m.pPlus);
					}
					if (maxB > minB)
						order2.set(rank, order2.get(rank) + 1);
				} else {
					score2.set(rank, score2.get(rank) + maxScore);
				}
			}
		}
		int oldNonMarkerScore = calculateNonMarkerScore(cis);
		int oldScore = calculateScore(cis, false) + oldNonMarkerScore;
		for (ContigInterval ci : cis) {
			rank = ci.getRank();
			int newNonMarkerScore = oldNonMarkerScore + nonMarkerScoreChange(cis, rank, rank, 0, true, 0);
			score2.set(rank, score2.get(rank) + newNonMarkerScore);
		}

		ArrayList<Double> ret = new ArrayList<Double>();
		//boolean changes = false; 
		for (ContigInterval c : cis) {
			rank = c.getRank();
			int support = oldScore - score2.get(rank);
			int order = 0;
			if (support == 0 && order1.get(rank) > order2.get(rank)) {
				if (c.getOrientation() >= 0)
					order = 1;
				else
					order = -1;
			}
			else if (support == 0 && order1.get(rank) < order2.get(rank)) {
				if (c.getOrientation() >= 0)
					order = -1;
				else
					order = 1;
			}
			//if (support < 0) {
			//	c.flipOrientation();
			//	changes = true;
			//}
			//System.err.println("order_support\t" + c +"\t" + ((support == 0) ? 0.5 * order: support));
			ret.add(((support == 0) ? (0.5 * order): support));
		}
		//if (!changes)
		//	break;
		//}
		return ret;
	}
	
	//flip(i..cr)
	//move(i...cr) to new_pos w/wo flipping
	//move cr to new position
	
	private int calculateNonMarkerScore(ArrayList<ContigInterval> cis)
	{
		return calculateChainScore(cis) + calculateProximityScore(cis);
	}
	
	
	//start...end => start+moveDirection...end+moveDirection
	public static void changeOrder(ArrayList<ContigInterval> cis, int start, int end, int moveDirection, boolean flip)
	{
		int n = cis.size();
		if (end + moveDirection >= n || start + moveDirection < 0) {
			System.err.println("change too large");
			return;
		}
		if (moveDirection == 0 && !flip) // nothing to do
			return;
		
		ArrayList<ContigInterval> cis_new = new ArrayList<ContigInterval>();
		if (moveDirection <= 0) {
			for (int i = 0; i < start + moveDirection; ++i)
				cis_new.add(cis.get(i));
			
			if (!flip)
				for (int i = start; i <= end; ++i)
					cis_new.add(cis.get(i));
			else
				for (int i = end; i >= start; --i) { //reverse and flip
					ContigInterval ci = cis.get(i);
					ci.flipOrientation();
					cis_new.add(ci);
				}

			for (int i = start + moveDirection; i < start; ++i)
				cis_new.add(cis.get(i));
			for (int i = end + 1; i < n; ++i)
				cis_new.add(cis.get(i));
		} else {
			for (int i = 0; i < start; ++i)
				cis_new.add(cis.get(i));

			for (int i = 1; i <= moveDirection; ++i)
				cis_new.add(cis.get(i + end));
			
			if (!flip)
				for (int i = start; i <= end; ++i)
					cis_new.add(cis.get(i));
			else
				for (int i = end; i >= start; --i) { //reverse and flip
					ContigInterval ci = cis.get(i);
					ci.flipOrientation();
					cis_new.add(ci);
				}

			for (int i = end + moveDirection + 1; i < n; ++i)
				cis_new.add(cis.get(i));
		}
		cis.clear();
		cis.addAll(cis_new);
	}
	
	//reverse of changeOrder
	public static void changeOrderReverse(ArrayList<ContigInterval> cis, int start, int end, int moveDirection, boolean flip)
	{
		changeOrder(cis, start + moveDirection, end + moveDirection, -moveDirection, flip);
	}
	public int nonMarkerScoreChange(ArrayList<ContigInterval> cis, int start, int end, int moveDirection, boolean flip, int flipScore)
	{
		int ret = ((flip) ? flipScore : 0);
		if (moveDirection == 0) { // flip 
			if (flip) {
				if (start > 0) {
					ret -= calculateChainScore(cis.get(start - 1), cis.get(start));
					ret += calculateChainScore(cis.get(start - 1), false, cis.get(end), true);
				}
				if (end + 1 < cis.size()) {
					ret -= calculateChainScore(cis.get(end), cis.get(end + 1));
					ret += calculateChainScore(cis.get(start), true, cis.get(end + 1), false);
				}
			} // else nothing to do...
		} else { // move one or multiple contigs...
			ContigInterval sc = cis.get(start);
			ContigInterval ec = cis.get(end);

			if (start > 0) {
				ret -= calculateChainScore(cis.get(start - 1), sc);
				if (end + 1 < cis.size())
					ret += calculateChainScore(cis.get(start - 1), cis.get(end + 1));
			}
			if (end + 1 < cis.size())
				ret -= calculateChainScore(ec, cis.get(end + 1));
			
			if (moveDirection <= 0) {
				if (flip) {
					if (start + moveDirection > 0) {
						ret -= calculateChainScore(cis.get(start + moveDirection - 1), cis.get(start + moveDirection));
						ret += calculateChainScore(cis.get(start + moveDirection - 1), false, ec, true);
					}
					ret += calculateChainScore(sc, true, cis.get(start + moveDirection), false);
				} else {
					if (start + moveDirection > 0) {
						ret -= calculateChainScore(cis.get(start + moveDirection - 1), cis.get(start + moveDirection));
						ret += calculateChainScore(cis.get(start + moveDirection - 1), sc);
					}
					ret += calculateChainScore(ec, cis.get(start + moveDirection));
				}
			} else { // moveDirection > 0
				if (flip) {
					if (end + moveDirection + 1 < cis.size()) {
						ret -= calculateChainScore(cis.get(end + moveDirection), cis.get(end + moveDirection + 1));
						ret += calculateChainScore(sc, true, cis.get(end + moveDirection + 1), false);
					}
					ret += calculateChainScore(cis.get(end + moveDirection), false, ec, true);
				} else {
					if (end + moveDirection + 1 < cis.size()) {
						ret -= calculateChainScore(cis.get(end + moveDirection), cis.get(end + moveDirection + 1));
						ret += calculateChainScore(ec, cis.get(end + moveDirection + 1));
					}
					ret += calculateChainScore(cis.get(end + moveDirection), sc);
				}
			}
		}
		if (prox != null)
			ret += prox.scoreChange(cis, start, end, moveDirection, flip);
		return ret;
	}

	private boolean linkedIntervals(ArrayList<ContigInterval> cis, int c1, int c2)
	{
		return (calculateChainScore(cis.get(c1), cis.get(c2)) > 0 || (prox != null && prox.linkedIntervals(cis, c1, c2)));
	}
	
	//join contig to linked contigs next to it and try to improve...
	//not thread safe... fix?: separate scaffoldingLink initialisation
	private int calculateBestLinkedContig(ArrayList<ContigInterval> cis, int cr, boolean verbose){

		int start = cr;
		if (start > 0 && linkedIntervals(cis, start - 1, start))
			return 0; // only start from the first contig of each linked chain of contigs 
	
		int numCis = cis.size();

		int end = cr;
		while (end + 1 < numCis && linkedIntervals(cis, end, end + 1))
			++end;
		
		boolean mapMarkers = false;
		for (int i = start; i <= end; ++i)
			if (cis.get(i).markers.size() > 0) {
				mapMarkers = true;
				break;
			}
		if (!mapMarkers) // no need to go further without markers...
			return 0;

		int numLinked = end - start + 1;
		
		if (start < cr || start == end || numLinked == numCis) // starts at c, at least two ContigIntervals and not connected to all
			return 0;
		
		ContigInterval ci_new = new ContigInterval("\t", 1, 1); // dummy name, this cannot be loaded from a file...
		
		int flipScore = 0; // calculate flipsScore to take into account possible decreasing score when flipping start..end
		for (int i = start; i < end; ++i)
			flipScore += calculateChainScore(cis.get(i + 1), true, cis.get(i), true) - calculateChainScore(cis.get(i), cis.get(i + 1));  
		
		//add scaffolding links, use scaffoldingLink_tmp not to mess other links... 
		scaffoldingLink_tmp.clear();
		for (int i = -1; i < numCis; ++i)
			if (i < start || i > end) {
				int moveDirection = ((i < start) ? i - start + 1 : i - end);
				int nmss = nonMarkerScoreChange(cis, start, end, moveDirection, false, flipScore);
				int nmss2 = nonMarkerScoreChange(cis, start, end, moveDirection, true, flipScore);
				if (i >= 0) {
					ContigInterval ci = cis.get(i);
					if (moveDirection != 0) // with links for moveDirection==0, we end up updating start-1...end as well
						scaffoldingLink_tmp.put(calculateKey(ci, ci_new), nmss);
					scaffoldingLink_tmp.put(calculateKey(ci, ci.getOrientation() >= 0, ci_new, false), nmss2);
				} else { // move to the very first position...
					int firstMarker = ((start > 0) ? 0 : end + 1);
					ContigInterval ci = cis.get(firstMarker);
					if (moveDirection != 0)  // with links for moveDirection==0, we end up updating start...end+1 as well
						scaffoldingLink_tmp.put(calculateKey(ci_new, ci), nmss);
					scaffoldingLink_tmp.put(calculateKey(ci_new, false, ci, ci.getOrientation() >= 0), nmss2);
				}
			}
		
		//and markers
		for (int map = 0; map < numMaps; ++map) {
			for (int i = start; i <= end; ++i) {
				for (Marker m : cis.get(i).getMarkers(map))
					ci_new.addMarker(new Marker(m), map);
			}
		}
		
		ArrayList<ContigInterval> cis_new = new ArrayList<ContigInterval>();
		for (int i = 0; i < start; ++i)
			cis_new.add(cis.get(i));
		cis_new.add(ci_new);
		for (int i = end + 1; i < numCis; ++i)
			cis_new.add(cis.get(i));

		//int os = this.calculateScore(cis, false) + calculateNonMarkerScore(cis);
		//int ns2 = this.calculateScore(cis_new, false) + calculateNonMarkerScore(cis_new);
		//int ns = this.calculateScore(cis_new, false) + calculateNonMarkerScore(cis_new);
		//System.err.println("linked " + ret + " " + start + "," + end + " " +ns2 + " " + os + " " + ns);

		Proximity tmpProx = prox;
		prox = null; // do not use proximity
		HashMap<String, Integer> tmpScaffoldingLink = scaffoldingLink;
		scaffoldingLink = scaffoldingLink_tmp;
		HashMap<String, Integer> tmpChainLink = chainLinkHash; 
		chainLinkHash = new HashMap<String, Integer>(); //nor chainLinks

		int ret = calculateBest(cis_new, ci_new, true);

		prox = tmpProx;
		scaffoldingLink = tmpScaffoldingLink;
		chainLinkHash = tmpChainLink;
		
		if (ret > 0) { // update cis
			//System.err.println("linked improvement " + ret);
			int crn = 0;
			for (ContigInterval ci : cis_new)
				if (ci == ci_new)
					break;
				else 
					++crn;
			
			cis_new.clear();
			//int os = calculateScore(cis, false) + calculateNonMarkerScore(cis);
			
			for (int i = 0; i < crn; ++i)
				if (i >= start)
					cis_new.add(cis.get(i + numLinked));
				else
					cis_new.add(cis.get(i));
			int orientation = ci_new.getOrientation();
			if (orientation >= 0)
				for (int i = start; i <= end; ++i)
					cis_new.add(cis.get(i));
			else
				for (int i = end; i >= start; --i) {
					cis.get(i).flipOrientation();
					cis_new.add(cis.get(i));
				}
			for (int i = crn; i < numCis - numLinked; ++i)
				if (i >= start)
					cis_new.add(cis.get(i + numLinked));
				else
					cis_new.add(cis.get(i));
			int ns = calculateScore(cis_new, false) + calculateNonMarkerScore(cis_new);
			//if (os >= ns) {
			//	if (orientation < 0)
			//		for (int i = start; i <= end; ++i)
			//			cis.get(i).flipOrientation();
			//	System.err.println("Error: Score not improving! " + start + "-" + end + "\t" + crn + "\t" + orientation);
			//	return 0;
			//}
			System.err.println("LINKED\t" + cis.get(start) + "\t" + (end - start + 1) + "\t" + ns);
			cis.clear();
			cis.addAll(cis_new);
		}
		return ret;
	}
	//join contig to linked contigs next to it and try to improve...
	//thread safe version...
	private int calculateBestLinkedContig_ts(ArrayList<ContigInterval> cis, int cr, boolean verbose){

		int start = cr;
		if (start > 0 && linkedIntervals(cis, start - 1, start))
			return 0; // only start from the first contig of each linked chain of contigs 
	
		int numCis = cis.size();

		int end = cr;
		while (end + 1 < numCis && linkedIntervals(cis, end, end + 1))
			++end;

		int numLinked = end - start + 1;

		if (start < cr || start == end || numLinked == numCis) // starts at c, at least two ContigIntervals and not connected to all
			return 0;
		
		boolean mapMarkers = false;
		for (int i = start; i <= end; ++i)
			if (cis.get(i).markers.size() > 0) {
				mapMarkers = true;
				break;
			}

		if (!mapMarkers) // no need to go further without markers...
			return 0;

		//update rank
		int rank = 0;
		for (ContigInterval ci : cis)
			ci.setRank(rank++);
		
		int bestRes[] = new int[numCis + 1];
		int bestRes2[] = new int[numCis + 1];
	
 		for (int map = 0; map < numMaps; ++map) {
 			
			ArrayList<Marker> cMarkers = new ArrayList<Marker>();
			for (int i = start; i <= end; ++i)
				cMarkers.addAll(cis.get(i).getMarkers(map));

			ArrayList<Marker> markers = new ArrayList<Marker>();
			for (int i = 0; i < start; ++i)
				markers.addAll(cis.get(i).getMarkers(map));
			for (int i = end + 1; i < numCis; ++i)
				markers.addAll(cis.get(i).getMarkers(map));

			int endBin = getMaxBin(markers);
			int startBin = getMinBin(markers);

			for (Marker m : cMarkers) {
				endBin = Math.max(endBin, m.maxBin());
				startBin = Math.min(startBin, m.minBin());
			}
			//System.err.println(startBin);
			//System.err.println(endBin);
			
			int numMarkers = markers.size();
			int S[][] = new int[numMarkers + 1][endBin - startBin + 1];
		
			int Sb[][] = new int[numMarkers + 1][endBin - startBin + 1];

			int Sc[][] = new int[1 + cMarkers.size()][endBin - startBin + 1];
		
			forwardFast1(markers, S, startBin, endBin);
			//simple way to calculate (about correct) backward table...
			Collections.reverse(markers);
			forwardFast2(markers, Sb, startBin, endBin);
			Collections.reverse(markers);
			
			int old = -1;
			for (int mi = 0; mi <= numMarkers; ++mi) {
				
				if (mi < numMarkers) {  
					Marker m = markers.get(mi); // this is the next marker
					numCis = m.ci.getRank();
				} else
					numCis = cis.size();
				
				if (numCis != old) {
					for (int b = 0; b <= endBin - startBin; ++b)
						Sc[0][b] = S[mi][b];
					
					forwardFast1(cMarkers, Sc, startBin, endBin);
					
					int Scc[] = Sc[cMarkers.size()];
					for (int b = 1; b <= endBin - startBin; ++b) { // make table non-decreasing
						if (Scc[b] < Scc[b - 1])
							Scc[b] = Scc[b - 1];
					}
					int Sbb[] = Sb[numMarkers - mi];
					int max = -1; 
					for (int b = 0; b <= endBin - startBin; ++b) {
						int s = Scc[b] + Sbb[b]; 
						if (s > max)
							max = s;
					}
					//System.out.println(max);
					for (int i = old + 1; i <= numCis; ++i)  
						bestRes[i] += max;
					
					Collections.reverse(cMarkers);
					forwardFast1(cMarkers, Sc, startBin, endBin);
					Collections.reverse(cMarkers);
					
					for (int b = 1; b <= endBin - startBin; ++b) { // make table non-decreasing
						if (Scc[b] < Scc[b - 1])
							Scc[b] = Scc[b - 1];
					}
					max = -1; 
					for (int b = 0; b <= endBin - startBin; ++b) {
						int s = Scc[b] + Sbb[b]; 
						if (s > max)
							max = s;
					}
					//System.out.println(max);
					for (int i = old + 1; i <= numCis; ++i)  
						bestRes2[i] += max;
				}
				old = numCis;
			}
		}

 		numCis = cis.size();
		int flipScore = 0; // calculate flipsScore to take into account possible decreasing score when flipping start..end
		for (int i = start; i < end; ++i)
			flipScore += calculateChainScore(cis.get(i + 1), true, cis.get(i), true) - calculateChainScore(cis.get(i), cis.get(i + 1));  

	 	int oldNonMarkerScore = calculateNonMarkerScore(cis);
	 	int oldScore = calculateScore(cis, false) + oldNonMarkerScore;
		for (int ci = 0; ci <= numCis; ++ci) {
			int mD = ci - start;
			if (ci > start)
				if (ci > end)
					mD = ci - end - 1;
				else
					mD = 0;
			bestRes[ci] += oldNonMarkerScore + nonMarkerScoreChange(cis, start, end, mD, false, flipScore); 
			bestRes2[ci]+= oldNonMarkerScore + nonMarkerScoreChange(cis, start, end, mD, true, flipScore); 
		}
		//for (int ci = 0; ci <= cis.size(); ++ci) {
		//	System.err.print(bestRes[ci] + " ");
		//}
		//System.err.println();
	
		int max = 0;
		for (int ci = 0; ci <= cis.size(); ++ci)
			if (bestRes[ci] > bestRes[max]) 
				max = ci;
		
		int max2 = 0;
		for (int ci = 0; ci <= cis.size(); ++ci)
			if (bestRes2[ci] > bestRes2[max2]) 
				max2 = ci;
		
		if (bestRes[max] > oldScore || bestRes2[max2] > oldScore) {
			int pos = max2;
			boolean flip = true;
			if (bestRes[max] >= bestRes2[max2]) {
				flip = false;
				pos = max;
				if (verbose)
					System.err.println("LINKED_MOVE\t" + cis.get(start) + "\t" +  numLinked + "\t" + bestRes[max]);// + "\t" + oldScore);

			} else {
				if (verbose)
					System.err.println(((max2 == start || max2 == start + 1) ? "LINKED_FLIP\t":"LINKED_MOVE+FLIP\t") + cis.get(start) + "\t" + numLinked + "\t" + bestRes2[max2]);// + "\t" + oldScore);
			}
			
			//table index is the new rank for the contig start... without flipping...
			int mD = pos - start;
			if (pos > start)
				if (pos > end)
					mD = pos - end - 1;
				else
					mD = 0;
			changeOrder(cis, start, end, mD, flip);
		}
		//int ns = (calculateNonMarkerScore(cis) + calculateScore(cis, false));
		//int ns2 = Math.max(bestRes[max], bestRes2[max2]);
		//if (ns != ns2)
		//	System.err.println("LINKED_DIFF " + oldScore + " " + ns + " " + ns2);
		return Math.max(bestRes[max], bestRes2[max2]) - oldScore;
	}

	
	private int calculateBest(ArrayList<ContigInterval> cis, ContigInterval c, boolean verbose) {
		
		//update rank
		int numCis = 0;
		for (ContigInterval ci : cis)
			ci.setRank(numCis++);

		boolean mapMarkers = (c.markers.size() > 0);

		int oldNonMarkerScore = calculateNonMarkerScore(cis);
		int oldScore = oldNonMarkerScore + ((mapMarkers) ? calculateScore(cis, false) : 0);

		int cr = c.getRank();

		//try flipping and/or moving multiple Contigs [i...cr]
		if (cr > 0) {
			int minFlipStart = cr - 1;
			while (minFlipStart >= 0 && linkedIntervals(cis, minFlipStart, minFlipStart + 1)) {
					--minFlipStart;
			}
			
			int flipS = 0; //calculate this to fix the nonsymmetrical pair-wise scores...
			for (int flipStart = cr - 1; flipStart > minFlipStart; --flipStart) {

				if (!mapMarkers && cis.get(flipStart).markers.size() > 0) {
					mapMarkers = true;
					oldScore += calculateScore(cis, false);
				}
				
				flipS -= calculateChainScore(cis.get(flipStart), cis.get(flipStart + 1)); //handles the nonsymmetrical pair-wise scores...
				flipS += calculateChainScore(cis.get(flipStart + 1), true, cis.get(flipStart), true); //handles the nonsymmetrical pair-wise scores...
				
				int flipScore = nonMarkerScoreChange(cis, flipStart, cr, 0, true, flipS); 

				if (flipScore >= 0) { //chainScore must not get worse, now calculate the marker score as well		
					
					changeOrder(cis, flipStart, cr, 0, true);
					int nms = calculateNonMarkerScore(cis);
					if (nms - oldNonMarkerScore != flipScore) {
						System.err.println("Error:scores differ " + (nms - oldNonMarkerScore) + " " + flipScore);
					}
					int newScore = ((mapMarkers) ? calculateScoreFast(cis) + nms : nms);
					if (newScore > oldScore) {// and the total score improves
						if (verbose)
							System.err.println("FLIP\t" + c + "\t" + (cr - flipStart + 1) + "\t" + (newScore));
						return newScore - oldScore;
					} else
						changeOrderReverse(cis, flipStart, cr, 0, true);
				}

				//try moving multiple Contigs [i...cr] with or without flipping...
	out:		for (int direction = -1; direction < 2;direction+=2){
					
					ArrayList<Integer> scores = new ArrayList<Integer>(); 
					for (int movePosition = direction; cr + movePosition < numCis && flipStart + movePosition >= 0; movePosition+=direction) {
						int moveScore1 = nonMarkerScoreChange(cis, flipStart, cr, movePosition, false, 0);
						int moveScore2 = nonMarkerScoreChange(cis, flipStart, cr, movePosition, true, flipS);
						scores.add(-Math.max(moveScore1, moveScore2));
					}
					//sort possible positions based on the chain (nonMarker) score... 
				    Misc.ArrayIndexComparator<Integer> comparator = new Misc.ArrayIndexComparator<Integer>(scores);
				    Integer[] indexes = comparator.createIndexArray();
				    Arrays.sort(indexes, comparator);  

					int maxEvaluations = 4;
					
					for (int i : indexes) {
						int movePosition = (i + 1) * direction; 
						int moveScore1 = nonMarkerScoreChange(cis, flipStart, cr, movePosition, false, 0);
						if (moveScore1 >= 0) {
							changeOrder(cis, flipStart, cr, movePosition, false);
							int score = ((mapMarkers) ? calculateScoreFast(cis) + calculateNonMarkerScore(cis) : calculateNonMarkerScore(cis));
							if (score > oldScore) {
								if (verbose)
									System.err.println("MOVE\t" + c + "\t" + (cr - flipStart + 1) + "\t" + score);
								return score - oldScore;
							} else
								changeOrderReverse(cis, flipStart, cr, movePosition, false);
						}
						int moveScore2 = nonMarkerScoreChange(cis, flipStart, cr, movePosition, true, flipS);
						if (moveScore2 >= 0) {
							changeOrder(cis, flipStart, cr, movePosition, true);
							int score = ((mapMarkers) ? calculateScoreFast(cis) + calculateNonMarkerScore(cis) : calculateNonMarkerScore(cis));
							if (score > oldScore) {
								if (verbose)
									System.err.println("MOVE+FLIP\t" + c + "\t" + (cr - flipStart + 1) + "\t" + score);
								return score - oldScore;
							} else
								changeOrderReverse(cis, flipStart, cr, movePosition, true);
						}
						//if (--maxEvaluations <= 0 || (moveScore1 < 0 && moveScore2 < 0)) // maximum number of evaluations reached...
						if (--maxEvaluations <= 0) // maximum number of evaluations reached...
							continue out;
					}
				}
			}
		}
		//Move only contig c into new position
		int bestRes[] = new int[cis.size() + 1];
		int bestRes2[] = new int[cis.size() + 1];
		
		//System.err.println(oldScore + "\t" + oldChainScore);
		if (c.markers.size() > 0) { // no need without markers in c...
	 		for (int map = 0; map < numMaps; ++map) {
				ArrayList<Marker> markers = getMarkers(cis, c, map);
				ArrayList<Marker> cMarkers = c.getMarkers(map);
	
				int endBin = getMaxBin(markers);
				int startBin = getMinBin(markers);
	
				for (Marker m : cMarkers) {
					endBin = Math.max(endBin, m.maxBin());
					startBin = Math.min(startBin, m.minBin());
				}
				//System.err.println(startBin);
				//System.err.println(endBin);
				
				int numMarkers = markers.size();
				int S[][] = new int[numMarkers + 1][endBin - startBin + 1];
			
				int Sb[][] = new int[numMarkers + 1][endBin - startBin + 1];
	
				int Sc[][] = new int[1 + cMarkers.size()][endBin - startBin + 1];
			
				forwardFast1(markers, S, startBin, endBin);
				//simple way to calculate (about correct) backward table...
				Collections.reverse(markers);
				forwardFast2(markers, Sb, startBin, endBin);
				Collections.reverse(markers);
				
				int old = -1;
				for (int mi = 0; mi <= numMarkers; ++mi) {
					
					if (mi < numMarkers) {  
						Marker m = markers.get(mi); // this is the next marker
						numCis = m.ci.getRank();
					} else
						numCis = cis.size();
					
					if (numCis != old) {
						for (int b = 0; b <= endBin - startBin; ++b)
							Sc[0][b] = S[mi][b];
						
						forwardFast1(cMarkers, Sc, startBin, endBin);
						
						int Scc[] = Sc[cMarkers.size()];
						for (int b = 1; b <= endBin - startBin; ++b) { // make table non-decreasing
							if (Scc[b] < Scc[b - 1])
								Scc[b] = Scc[b - 1];
						}
						int Sbb[] = Sb[numMarkers - mi];
						int max = -1; 
						for (int b = 0; b <= endBin - startBin; ++b) {
							int s = Scc[b] + Sbb[b]; 
							if (s > max)
								max = s;
						}
						//System.out.println(max);
						for (int i = old + 1; i <= numCis; ++i)  
							bestRes[i] += max;
						
						Collections.reverse(cMarkers);
						forwardFast1(cMarkers, Sc, startBin, endBin);
						Collections.reverse(cMarkers);
						
						for (int b = 1; b <= endBin - startBin; ++b) { // make table non-decreasing
							if (Scc[b] < Scc[b - 1])
								Scc[b] = Scc[b - 1];
						}
						max = -1; 
						for (int b = 0; b <= endBin - startBin; ++b) {
							int s = Scc[b] + Sbb[b]; 
							if (s > max)
								max = s;
						}
						//System.out.println(max);
						for (int i = old + 1; i <= numCis; ++i)  
							bestRes2[i] += max;
					}
					old = numCis;
				}
			}
		} else {
			if (mapMarkers) // if there are no markers in c, then we can only consider nonMarkerScore after this...
				oldScore = oldNonMarkerScore;
			mapMarkers = false;
		}
		for (int ci = 0; ci <= cis.size(); ++ci) {
			int mD = ci - cr;
			if (ci > cr)
				--mD;
			bestRes[ci] += oldNonMarkerScore + nonMarkerScoreChange(cis, cr, cr, mD, false, 0); 
			bestRes2[ci]+= oldNonMarkerScore + nonMarkerScoreChange(cis, cr, cr, mD, true, 0); 
		}
		//for (int ci = 0; ci <= cis.size(); ++ci) {
		//	System.err.print(bestRes[ci] + " ");
		//}
		//System.err.println();

		int max = 0;
		for (int ci = 0; ci <= cis.size(); ++ci)
			if (bestRes[ci] > bestRes[max]) 
				max = ci;
		
		int max2 = 0;
		for (int ci = 0; ci <= cis.size(); ++ci)
			if (bestRes2[ci] > bestRes2[max2]) 
				max2 = ci;
		
		if (bestRes[max] > oldScore || bestRes2[max2] > oldScore) {
			int pos = max2;
			if (bestRes[max] >= bestRes2[max2]) {
				pos = max;
				if (verbose)
					System.err.println("MOVE\t" + c + "\t1\t" + bestRes[max]);
			} else {
				if (verbose)
					System.err.println(((max2 == cr || max2 == cr + 1) ? "FLIP\t":"MOVE+FLIP\t") + c + "\t1\t" + bestRes2[max2]);
				c.flipOrientation();
			}
			cis.remove(cr);
			if (pos > cr)
				cis.add(pos - 1, c);
			else
				cis.add(pos, c);
		}
		return Math.max(bestRes[max], bestRes2[max2]) - oldScore;
 	}
	
	private ArrayList<ContigInterval> cloneContigs(ArrayList<ContigInterval> cis)
	{
		ArrayList<ContigInterval> ret = new ArrayList<ContigInterval>();
		for (ContigInterval c : cis)
			ret.add(new ContigInterval(c));
		return ret;
	}
	
	private void improveAnchoring()
	{
		improveAnchoring(intervalsInChr, true, true);
	}
	
	//parallel implementation of calculateBest
	private class CalculateBest implements Runnable{
		private AtomicInteger index;
		private AtomicBoolean stop;
		private ArrayList<ContigInterval> cis;
		private ArrayList<ContigInterval> cisOrig;
		private ArrayList<Integer> perm;
		private boolean verbose;
		private boolean linked;
		
		private CalculateBest next;
		private Thread thread; 
		
		private int score = 0;
		private int scoreIndex = 0;
		CalculateBest(ArrayList<ContigInterval> cis, ArrayList<Integer> perm, AtomicInteger index, AtomicBoolean stop, boolean verbose, boolean linked) {
			this.index = index;
			this.stop = stop;
			//clone cis
			this.cis = cloneContigs(cis);
			this.cisOrig = new ArrayList<ContigInterval>();
			cisOrig.addAll(this.cis);
			this.perm = perm;
			this.verbose = verbose;
			this.linked = linked;
		}
		public void setNext(CalculateBest next)
		{
			this.next = next;
		}
		public void setScore(int score)
		{
			this.score = score;
		}
		@Override
		public void run() {
			long time = System.currentTimeMillis();
			score = 0;
			thread = null;
			while (!stop.get()) { 
				if (thread == null && next != null && System.currentTimeMillis() - time >= 1000) { // start a new thread every 1 secs
					next.setCis(this.cis);
					if (stop.get()) //it is worth breaking here
						break; 
					thread = new Thread(next);
					if (stop.get()) //and maybe here as well even when we don't use thread...
						break;
					thread.start();
				}
				int i = index.getAndIncrement();
				if (i >= perm.size())
					break;
				if (linked)
					score = calculateBestLinkedContig_ts(cis, perm.get(i), verbose); // TODO: fix perm to cisOrig order
				else
					score = calculateBest(cis, cisOrig.get(perm.get(i)), verbose);
				if (score > 0) {
					scoreIndex = i;
					stop.set(true);
				}
			}
			if (thread != null && thread.isAlive())
				try {
					thread.join();
				} catch (Exception e) {
					e.printStackTrace();
					System.exit(-1);
				}
		}
		public void directRun(int i) {
			if (linked)
				score = calculateBestLinkedContig_ts(cis, perm.get(i), verbose);
			else
				score = calculateBest(cis, cisOrig.get(perm.get(i)), verbose);
		}
		public void setCis(ArrayList<ContigInterval> cis)
		{
			// avoid clone
			//this.cis = cloneContigs(cis);
			makeCisIdentical(this.cis, cis);
		}
		public int getScore()
		{
			return score;
		}
		public int getScoreIndex()
		{
			return scoreIndex;
		}
	}
	// make two lists of ContigIntervals indentical by placement and orientation, source=>target
	// lists must contain exactly the same contigs...
	private void makeCisIdentical(ArrayList<ContigInterval> target, ArrayList<ContigInterval> source)
	{
		HashMap<ContigInterval, ContigInterval> map1to1 = new HashMap<ContigInterval, ContigInterval>();
		for (ContigInterval c : target) {
			if (map1to1.containsKey(c)) {
				System.err.println("Error: same contig multiple times...");

				for (ContigInterval c2 : target) {
					System.err.println(c2);
					
				}
				System.exit(-1);
				
			}
			map1to1.put(c, c);
		}
		target.clear();
		for (ContigInterval c : source) {
			ContigInterval c2 = map1to1.get(c);
			if (c2 == null) {
				System.err.println("Error: some contig(s) not found...");
				System.exit(-1);
			}
			target.add(c2);
			c2.setOrientation(c.getOrientation());
		}
	}
	
	private void improveAnchoring(ArrayList<ContigInterval> cis, boolean finalRound, boolean verbose)
	{
		if (cis.size() == 0)
			return;

		if (finalRound)
			System.out.print("#initial score " +  (calculateScore(cis, false) + calculateNonMarkerScore(cis)));
		
		
		int numCis = cis.size(); 

		ArrayList<Integer> perm = new ArrayList<Integer>();
		for (int i = 0; i < numCis; ++i)
			perm.add(i);

		boolean foundBetter = false;

		
		int iteration = 0;
		
		int bestScore = -1;
		ArrayList<ContigInterval> cis_best = new ArrayList<ContigInterval>();
		ArrayList<Integer> cis_bestO = new ArrayList<Integer>();
		
		for (int run = 0; run < numRuns; ++run) {
			if (run > 0) {
				int n = cis.size();
				for (ContigInterval ci : cis)
					if (Math.random() < 0.05)
						ci.flipOrientation();
				if (n > 1)
					for (int flips = 0; flips < n; ++flips) { 
						int i = (int) (Math.random() * (n - 1));
						ContigInterval tmp = cis.get(i);
						cis.set(i, cis.get(i + 1));
						cis.set(i + 1, tmp);
					}
			}
			
			do {
				if (finalRound) {
					System.err.println("iteration " + (++iteration) + " score=" +  (calculateScore(cis, false) + calculateNonMarkerScore(cis)));
				}
				foundBetter = false;
				Collections.shuffle(perm); // random order
	
				if (numThreads <= 1) {
					ArrayList<ContigInterval> cisOrig = new ArrayList<ContigInterval>();
					cisOrig.addAll(cis);
					for (int i = 0; i < numCis; ++i)
						if (calculateBest(cis, cisOrig.get(perm.get(i)), verbose) > 0)
							foundBetter = true;
	
					//try linked version...
					if (!foundBetter)
						for (int i = 0; i < numCis; ++i) {
							//ArrayList<ContigInterval> cis_ = new ArrayList<ContigInterval>();
							//cis_.addAll(cis);
							//calculateBestLinkedContig_ts(cis_, perm.get(i), verbose);

							if (calculateBestLinkedContig_ts(cis, perm.get(i), verbose) > 0)
								foundBetter = true;
						}
				} else {
					for (int linked = 0; linked < 2; ++linked) {
						AtomicInteger index = new AtomicInteger();
						AtomicBoolean stop = new AtomicBoolean();
						
						CalculateBest[] cbs = new CalculateBest[numThreads];
						for (int t = 0; t < numThreads; ++t) {
							cbs[t] = new CalculateBest(cis, perm, index, stop, verbose, linked == 1);
							if (t > 0)
								cbs[t - 1].setNext(cbs[t]);
						}
						while (true) {
							cbs[0].run();
							int best = 0;
							for (int t = 1; t < numThreads; ++t) {
								if (cbs[t].getScore() > cbs[best].getScore())
									best = t;
							}
							if (cbs[best].getScore() > 0) {
								if (best > 0) {
									CalculateBest tmp = cbs[0]; //swap cbs[0] and cbs[best]
									cbs[0] = cbs[best];
									cbs[best] = tmp;
									for (int t = 1; t < numThreads; ++t)
										cbs[t - 1].setNext(cbs[t]);
									cbs[numThreads - 1].setNext(null);
								}
								for (int t = 1; t < numThreads; ++t) { //run other contigs that improved the score...
									if (cbs[t].getScore() > 0) {
										cbs[0].directRun(cbs[t].getScoreIndex());
									}
								}
								stop.set(false);
								foundBetter = true;
								for (int t = 0; t < numThreads; ++t)
									cbs[t].setScore(0);
							} else
								break;
						}
						if (foundBetter) {
							makeCisIdentical(cis, cbs[0].cis);
							break; // linked only if !foundBetter
						}
					}
					//if (!foundBetter) //calculateBestLinkedContig is not thread safe...
					//	for (int i = 0; i < numCis; ++i)
					//		if (calculateBestLinkedContig(cis, perm.get(i), verbose) > 0)
					//			foundBetter = true;

				}
			} while (foundBetter);
			int score = (calculateScore(cis, false) + calculateNonMarkerScore(cis));
			if (score > bestScore) {
				bestScore = score;
				cis_best.clear();
				cis_bestO.clear();
				for (ContigInterval ci : cis) {
					cis_best.add(ci);
					cis_bestO.add(ci.getOrientation());
				}
			} else {
				cis.clear();
				int i = 0;
				for (ContigInterval ci : cis_best) {
					cis.add(ci);
					ci.setOrientation(cis_bestO.get(i++));
				}
			}
		}
		
		if (finalRound) {
			ArrayList<Double> orderSupport = calculateSupport(cis);
			System.out.println(" final score " +  (calculateScore(cis, false) + calculateNonMarkerScore(cis)));
			printAnchoring(cis, orderSupport);
			findLikelyAssemblyErrors(cis);
		}
	}

	private void printAnchoring(ArrayList<ContigInterval> cis, ArrayList<Double> orderSupport) {
		printAnchoring(System.out, cis, orderSupport);
	}
	private void printAnchoring(PrintStream stream, ArrayList<ContigInterval> cis, ArrayList<Double> orderSupport) {
		stream.print("#contig\tstart\tend\torientation\tchr\torder_support\torig_start/alt_start\torig_end/alt_end\tlinked1\tvia_1\tlink_score1\tlinked2\tvia_2\tlink_score2\tpos_support");
		for (int map = 0; map < numMaps; ++map)
			stream.print("\tmap" + (map + 1));
		stream.println();
		
		boolean knownOrientation[] = new boolean[cis.size()]; // calculate whether global orientation is known...

		
		//calculate map info and follow linked2 links (fill forwards in knownOrientation)
		int prevC = 0;
		for (int cii = 0; cii < cis.size(); ++cii) {
			ContigInterval ci = cis.get(cii);
			boolean linked2 = false; // from chain file...
			if (cii + 1 < cis.size() && linkedIntervals(cis, cii, cii + 1))
				linked2 = true;
			if (!linked2) {
				for (int  map = 0; map < numMaps; ++map) {
					int max = Integer.MIN_VALUE;
					int min = Integer.MAX_VALUE;
					for (int c = prevC; c <= cii; ++c)
						for (Marker m : cis.get(c).getMarkers(map)) {
							max = Math.max(max, m.pPlus);
							min = Math.min(min, m.pPlus);
						}
					if (max > min)
						for (int c = prevC; c <= cii; ++c)
							knownOrientation[c] = true;
				}
				prevC = cii + 1;
			}
		}

		long  starts[] = new long[cis.size()]; // calculate start of each contig
		long  ends[] = new long[cis.size()]; // calculate end of each contig

		long  starts2[] = new long[cis.size()]; // calculate start of each contig
		long  ends2[] = new long[cis.size()]; // calculate end of each contig

		boolean prevSwapped = false; // fix for contigs cut out more than their length (a haplotype of two other contigs)
		for (int cii = 0; cii < cis.size(); ++cii) {
			ContigInterval ci = cis.get(cii);
			
			boolean linked1 = false; // from chain file...
			boolean linked2 = false; // from chain file...
			if (cii > 0 && linkedIntervals(cis, cii - 1, cii))
				linked1 = true;
			if (cii < cis.size() - 1 && linkedIntervals(cis, cii, cii + 1))
				linked2 = true;

			ContigInterval prev = ((linked1) ? cis.get(cii - 1): null);
			ContigInterval next = ((linked2) ? cis.get(cii + 1): null);
			
			long start = ci.getStart();
			long end = ci.getEnd();

			long start2 = ci.getStart();
			long end2 = ci.getEnd();

			if (prev != null) {
				String key = this.calculateKey(prev, ci);
				if (chainLinkHash.containsKey(key)) {
					if (ci.getOrientation() >= 0) {
						start = chainLinkHashCut2.get(key);
						start2 = chainLinkHashCut2a.get(key);
						if (prevSwapped) {
							long tmp = start;
							start = start2;
							start2 = tmp;
						}
					} else {
						end = chainLinkHashCut2.get(key);
						end2 = chainLinkHashCut2a.get(key);
						if (prevSwapped) {
							long tmp = end;
							end = end2;
							end2 = tmp;
						}
					}
				}
			}
			if (next != null) {
				String key = this.calculateKey(ci, next);
				if (chainLinkHash.containsKey(key)) {
					if (ci.getOrientation() >= 0) {
						end = chainLinkHashCut1.get(key);
						end2 = chainLinkHashCut1a.get(key);
						if (end < start - 1 && end2 > end) { // more than full length aligned...
							long tmp = end;
							end = end2;
							end2 = tmp;
							prevSwapped = true;
						} else
							prevSwapped = false;
					} else {
						start = chainLinkHashCut1.get(key);
						start2 = chainLinkHashCut1a.get(key);
						if (end < start - 1 && start2 < start) { // more than full length aligned...
							long tmp = start;
							start = start2;
							start2 = tmp;
							prevSwapped = true;
						} else
							prevSwapped = false;
					}
				} else
					prevSwapped = false;
			}

			starts[cii] = start;
			starts2[cii] = start2;

			ends[cii] = end;
			ends2[cii] = end2;
		}
		
		//finally try skipping negative length contigs and try to find new starts and ends... 
		for (int cii = 2; cii < cis.size(); ++cii) {
			if (ends[cii] >= starts[cii] - 1) { // non negative length
				int prevNN = cii - 1;
				while (prevNN >= 0 && ends[prevNN] < starts[prevNN] - 1) // negative length
					--prevNN;
				if (prevNN < cii - 1 && prevNN >= 0 && linkedIntervals(cis, prevNN, cii)) {
					ContigInterval ci = cis.get(cii);
					ContigInterval ci2 = cis.get(prevNN);
					String key = calculateKey(ci2, ci);
					if (chainLinkHash.containsKey(key)) {
						System.err.println("skipping contig " + cis.get(cii - 1));
						if (ci.getOrientation() >= 0) {
							starts[cii] = chainLinkHashCut2.get(key);
							starts2[cii] = chainLinkHashCut2a.get(key);
						} else {
							ends[cii] = chainLinkHashCut2.get(key);
							ends2[cii] = chainLinkHashCut2a.get(key);
						}
						if (ci2.getOrientation() >= 0) {
							ends[prevNN] = chainLinkHashCut1.get(key);
							ends2[prevNN] = chainLinkHashCut1a.get(key);
						} else {
							starts[prevNN] = chainLinkHashCut1.get(key);
							starts2[prevNN] = chainLinkHashCut1a.get(key);
						}
					}
				}
			}
		}
		
		boolean toggleScaffolding = true;
		for (int cii = 0; cii < cis.size(); ++cii) {
			ContigInterval ci = cis.get(cii);
			
			boolean linked1 = false; // from chain file...
			boolean linked2 = false; // from chain file...
			if (cii > 0 && linkedIntervals(cis, cii - 1, cii))
				linked1 = true;
			if (cii < cis.size() - 1 && linkedIntervals(cis, cii, cii + 1))
				linked2 = true;

			ContigInterval prev = ((linked1) ? cis.get(cii - 1): null);
			ContigInterval next = ((linked2) ? cis.get(cii + 1): null);
			
			String map_info = "";
			String orientation = "?";
			
			int posSupport = 0;
			
			for (int  map = 0; map < numMaps; ++map) {
				int max = Integer.MIN_VALUE;
				int min = Integer.MAX_VALUE;
				for (Marker m : ci.getMarkers(map)) {
					max = Math.max(max, m.pPlus);
					min = Math.min(min, m.pPlus);
					posSupport += m.inside(m.pPlus); 
				}
				if (max >= min)
					map_info = map_info + "\t" + min + "-" + max;
				else
					map_info = map_info + "\t-";
			}
			if (knownOrientation[cii]) 
				orientation = ((ci.getOrientation() >= 0) ?  "+" : "-");
			else if (linked1 || linked2) {
				if (!linked1)
					toggleScaffolding = !toggleScaffolding;
				if (!toggleScaffolding)
					orientation = ((ci.getOrientation() >= 0) ?  "++" : "--");
				else
					orientation = ((ci.getOrientation() >= 0) ?  "+++" : "---");
			}
			
			int scorePrev = 0;
			int scoreNext = 0;

			String prevScaffoldInfo = "null";		
			String nextScaffoldInfo = "null";		
			
			if (prev != null) {
				String key = this.calculateKey(prev, ci);
				if (chainLinkHash.containsKey(key)) {
					scorePrev = chainLinkHash.get(key);
					if (useChainAndPaf && scaffoldingLink.containsKey(key)) {
						scorePrev += scaffoldingLink.get(key);
						prevScaffoldInfo = "chain+" + scaffoldingLinkInfo.get(key);
					} else
						prevScaffoldInfo = "chain";
				} else if (scaffoldingLink.containsKey(key)) {
					scorePrev = scaffoldingLink.get(key);
					prevScaffoldInfo = scaffoldingLinkInfo.get(key);
				}
				int sp = ((prox == null) ? 0 : prox.linkScore(cis, cii - 1, cii));
				if (sp > 0) {
					scorePrev +=sp ;
					prevScaffoldInfo = prevScaffoldInfo + "+proximity";
				}
			}
			if (next != null) {
				String key = this.calculateKey(ci, next);
				if (chainLinkHash.containsKey(key)) {
					scoreNext = chainLinkHash.get(key);
					if (useChainAndPaf && scaffoldingLink.containsKey(key)) {
						scoreNext += scaffoldingLink.get(key);
						nextScaffoldInfo = "chain+" + scaffoldingLinkInfo.get(key);
					}
					else
						nextScaffoldInfo = "chain";
				} else if (scaffoldingLink.containsKey(key)) {
					scoreNext = scaffoldingLink.get(key);
					nextScaffoldInfo = scaffoldingLinkInfo.get(key);
				}
				int sp = ((prox == null) ? 0 : prox.linkScore(cis, cii, cii + 1));
				if (sp > 0) {
					scoreNext += prox.linkScore(cis, cii, cii + 1);
					nextScaffoldInfo = nextScaffoldInfo + "+proximity";
				}
			} // next != null

			//just get the start and end from the previously computed arrays...
			long start = starts[cii];
			long end = ends[cii];
			long start2 = starts2[cii];
			long end2 = ends2[cii];

			if (commentOutput)
				stream.print("#");
			stream.println(ci.getContig() + "\t" + start + "\t" + end + "\t" + orientation + "\t"
			        + ci.getChromosome() + "\t" + orderSupport.get(cii) + "\t" + ci.getStart() + (start2==ci.getStart() ? "":"/" + start2) +  "\t" + ci.getEnd() + (end2==ci.getEnd() ? "":"/" + end2)
					+ "\t" + (prev == null ? "null": prev.getContig()) + "\t" + prevScaffoldInfo + "\t" + scorePrev 
					+ "\t" + (next == null ? "null": next.getContig()) + "\t" + nextScaffoldInfo + "\t" + scoreNext + "\t" + posSupport + map_info);
		}
	}
	
	
	ArrayList<Double> logTable = new ArrayList<Double>();

	private double myLog(int x) {
		if (logTable.size() <= x) {
			for (int i = logTable.size(); i <= x; ++i)
				logTable.add(Math.log(i));
		}
		return logTable.get(x);
	}

	
	
	//markers is sorted within contig so no need to check orientation
	private String getOneInterval(ArrayList<Marker> markers, int m1, int m2)
	{
		Marker mm1 = markers.get(m1);
		Marker mm2 = markers.get(m2);
		ContigInterval ci = mm1.getContigInterval(); 
		long start = 0;
		long end = 0;
		
		if (m1 == 0 || !ci.equals(markers.get(m1 - 1).getContigInterval()))
				start = ci.getStart();
		else
				start = (mm1.getPosition() + markers.get(m1 - 1).getPosition()) / 2;

		if (m2 + 1 == markers.size() || !ci.equals(markers.get(m2 + 1).getContigInterval()))
				end = ci.getEnd();
		else
				end = (mm2.getPosition() + markers.get(m1 + 1).getPosition()) / 2;
		
		return ci.getContig() + "\t" + start + "\t" + end;  
	}

	//markers is sorted
	private String getIntervals(ArrayList<Marker> markers, int m1, int m2)
	{
		String ret = "";
		ContigInterval prev = null;
		int start = -1;
		int end = -1;
		for (int i = m1; i <= m2; ++i) {
			Marker m = markers.get(i);
			ContigInterval c = m.getContigInterval();
			if (c.equals(prev))
				end = i;
			else {
				if (start >= 0)
					ret += '\t' + getOneInterval(markers, start, end);
				start = i;
				end = i;
			}
			prev = c;
		}
		return ret + '\t' + getOneInterval(markers, start, end);
	}
	
	private class PossibleError implements Comparable<PossibleError>{
		private double ll;
		private String info;
		private int ci;
		private long pos1[], pos2[];
		public  PossibleError(double ll, int ci, String info) {
			this.ll = ll;
			this.info = info;
			this.ci = ci;
		}
		public  PossibleError(double ll, int ci, String info, long pos1[], long pos2[]) {
			this.ll = ll;
			this.info = info;
			this.ci = ci;
			this.pos1 = pos1;
			this.pos2 = pos2;
		}
		public String getInfo(){return info;}
		@Override
		public int compareTo(PossibleError other) {
			if (ll < other.ll)
				return 1;
			if (ll > other.ll)
				return -1;
			return 0;
		}
		public int getCi() {
			return ci;
		}
		public long[] getPos1() {
			return pos1;
		}
		public long[] getPos2() {
			return pos2;
		}
	}
	
	// log( (a/(a+b))^a * (b/(a+b))^b)
	private double logLike(int a, int b)
	{
		double ret = 0;
		if (a > 0)
			ret += a * (myLog(a) - myLog(a + b));
		if (b > 0)
			ret += b * (myLog(b) - myLog(a + b));
		return ret;
	}

	// log( (a/(a+b))^c * (b/(a+b))^d)
	private double logLike(int a, int b, int c, int d)
	{
		double ret = 0;
		if (c > 0)
			ret += c * (myLog(a) - myLog(a + b));
		if (d > 0)
			ret += d * (myLog(b) - myLog(a + b));
		return ret;
	}

	// find max maxNumErrors regions with more errors than expected, error rate is e0/(e0 + e1)
	private void calculateErrors2(ArrayList<Marker> markers_, ArrayList<PossibleError> errors, int e0, int e1, int maxNumErrors)
	{
		boolean inside[] = new boolean[markers_.size()];
		int cumulative[] = new int[markers_.size() + 1];
		
		int sum = 0;
		for (int mi = 0; mi < markers_.size(); ++mi) {
			Marker m = markers_.get(mi);
			if (m.inside(m.pPlus) > 0)
				inside[mi] = true;
			else
				++sum;
			cumulative[mi + 1] = sum;
		}
		if (sum == 0) // all inside or no markers
			return;
		
		ArrayList<Integer> markers = new ArrayList<Integer>();

		for (int mi = 0; mi < markers_.size(); ++mi) { // take only first and last marker of multiple adjacent markers !inside
			if ((mi == 0 || inside[mi - 1]) && !inside[mi])
				markers.add(mi);
			else
				if ((mi + 1 == markers_.size() || inside[mi + 1]) && !inside[mi])
					markers.add(mi);
		}
		//System.err.println(markers.size() + "\t" + markers_.size());

		int numMarkers = markers.size();
		
		ArrayList<int[]> found = new ArrayList<int[]>(); 
		for (int iteration = 0; iteration < maxNumErrors; ++iteration) {
			double max = 0;
			int max1 = -1;
			int max2 = -1;
			int maxo0 = 0;
			int maxo1 = 0;
out2:		for (int m1i = 0; m1i < numMarkers; ++m1i) {
				int m1 = markers.get(m1i); 
				for (int[] f : found)
					if (f[1] >= m1 && f[0] <= m1) // found intervals contain m1
						continue out2;
				for (int m2i = m1i; m2i < numMarkers; ++m2i) {
					int m2 = markers.get(m2i);
					for (int[] f : found)
						if (f[1] >= m2 && f[0] <= m2) // found intervals intersect m2 ([m1...m2])
							continue out2;
					int o0 = cumulative[m2 + 1] - cumulative[m1]; // number of markers outside interval
					int o1 = (m2 - m1 + 1) - o0; // number of markers within interval 
					//o0 * (e0 + e1) > e0 * (o0 + o1) (positive e and o) <=> o0 / (o0 + o1) > e0 / (e0 + e1)
					if ( o0 * (long)(e0 + e1) > e0 * (long)(o0 + o1)) { //more errors than expected
						double l = logLike(o0, o1) - logLike(e0, e1, o0, o1); //likelihood ratio
						if (l > max) {
							max = l;
							max1 = m1; 
							max2 = m2;
							maxo0 = o0;
							maxo1 = o1;
						}
					} else
						continue out2; // no need to continue if we fall short
				}
			}
			if (max > 0) {
				found.add(new int[]{max1, max2});
				ContigInterval ci = markers_.get(max1).getContigInterval();
				long pos1 = markers_.get(max1).getPosition();
				long pos1n = ((max1 > 0) ? markers_.get(max1 - 1).getPosition() + 1 : ci.getStart());
				
				long pos2 = markers_.get(max2).getPosition();
				long pos2n = ((max2 + 1 < markers_.size()) ? markers_.get(max2 + 1).getPosition() - 1 : ci.getEnd());
				
				//if (pos1 < pos1n || pos2 > pos2n) {
				//	System.exit(-1);
				//}

				String info = max + "\t" + maxo0 + "\t" + maxo1 + "\t" + ci.getContig() + "\t" + pos1n + "-" + pos1 + "\t" + pos2 + "-" + pos2n;
				errors.add(new PossibleError(max, ci.getRank(), info , new long[]{pos1n, pos1}, new long[]{pos2, pos2n}));
			}
		}
	}
	
	
	
	private ArrayList<Integer> calculateBinSum(ArrayList<Marker> markers, int endBin, boolean plusOrientation){
		//calculate map position distribution
		ArrayList<Integer> binSum = new ArrayList<Integer>();
		for (int i = 0; i <= endBin; ++i)
			binSum.add(0);

		//int plus = 0;
		//int minus = 0;
		for (Marker m : markers) {
			if (plusOrientation)
				binSum.set(m.pPlus, binSum.get(m.pPlus) + 1);
			else
				binSum.set(m.pMinus, binSum.get(m.pMinus) + 1);
		}
		return binSum;
	}

	private class LongComparator0 implements Comparator<long[]>
	{
		@Override
		public int compare(long[] l1, long[] l2){
			return Long.compare(l1[0], l2[0]);
		}
	}

	private class LongComparator1 implements Comparator<long[]>
	{
		@Override
		public int compare(long[] l1, long[] l2){
			return Long.compare(l1[1], l2[1]);
		}
	}
	
	//calculate how much score could be increasing if we split a contig...
	//make a gap in map positions L..R, and see how many points we could add there
	private ArrayList<Integer> increaseScore(ArrayList<Marker> markers, ArrayList<Integer> binSum){
		
		ArrayList<Integer> ret = new ArrayList<Integer>(); 

		int endBin = getMaxBin(markers);
		if (endBin < binSum.size() - 1)
			endBin = binSum.size() - 1;
		int startBin = getMinBin(markers);
		if (startBin > 0)
			startBin = 0;
		
		int numMarkers = markers.size();
		
		int S[][] = new int[numMarkers + 1][endBin - startBin + 1];

		int Sp[] = new int[endBin - startBin + 1]; 
	
		int Sb[][] = new int[numMarkers + 1][endBin - startBin + 1];

		//simple way to calculate (about correct) backward table...
		Collections.reverse(markers);
		forwardFast2(markers, Sb, startBin, endBin);
		Collections.reverse(markers);
		forwardFast1(markers, S, startBin, endBin);
		//Should we make S and Sb monotonic (+/- 1) ?

		for (int i = 0; i <= numMarkers; ++i) {
			int max = 0;
			//int maxL = 0;
			//int maxR = 0;

			//reuse P as PSp
			for (int b = 0; b < endBin - startBin; ++b) { // add binSum to S
				int bs = 0;
				if (b + startBin < binSum.size())
					bs = binSum.get(b + startBin);

				//either we take score S[i][b]+bs or skip b and get bs + score Sp[b-1]
				if (b == 0) {
					Sp[b] = S[i][b] + bs;
					//P[i][b] = b;
				} else {
					int s1 = S[i][b] + bs;
					int s2 = Sp[b - 1] + bs;
					if (s1 >= s2) {
						Sp[b] = s1;
						//P[i][b] = b;
					} else {
						Sp[b] = s2;
						//P[i][b] = P[i][b - 1];
					}
				}
			}
			

			for (int b = 0; b <= endBin - startBin; ++b) {
				int s = Sp[b] + Sb[numMarkers - i][b];
				if (s > max) {
					max = s;
					//maxR = P[i][b] + startBin;
					//maxL = b + startBin;
				}
			}
			ret.add(max);
		}
		return ret;
	}
	
	
	private ArrayList<PossibleError> findNonHaplotypeAssemblyErrors(ArrayList<ContigInterval> cis, boolean verbose)
	{
		ArrayList<PossibleError> ret = new ArrayList<PossibleError>(); 

		ArrayList<PossibleError> errors = new ArrayList<PossibleError>(); 
		int e0 = 0;
		int e1 = 0;
		for (int iteration = 0; iteration < 2; ++iteration) { //first iteration: calculate e0 and e1, second iteration find errors...
			for (int cii = 0; cii < cis.size(); ++cii) {
				ContigInterval ci = cis.get(cii);
				
				//Haplotypes are typical sources of errors, 
  				//by adjusting start and end, one could reduce errors caused by them...
				//best way would be to remove markers from one pair (with less markers) of haplotype and re-evaluate the maps...
/*				
				long start = ci.getStart();
				long end = ci.getEnd();
                int orientation = ci.getOrientation();
				if (cii > 0) {
					ContigInterval prev = cis.get(cii - 1);
					String key = calculateKey(prev, ci);
					Long cut = chainLinkHashCut2.get(key); //Cut2a
					if (cut != null)
						if (orientation >= 0)
							start = cut;
						else
							end = cut;
				}
				if (cii + 1  < cis.size()) {
					ContigInterval next = cis.get(cii + 1);
					String key = calculateKeyInverse(ci, next);
					Long cut = chainLinkHashCut2.get(key); //Cut2a
					if (cut != null)
						if (orientation >= 0)
							end = cut;
						else
							start = cut;
				}*/
	
				if (iteration == 0) {
					for (int  map = 0; map < numMaps; ++map)
						for (Marker m : ci.getMarkers(map))
							//if (m.getPosition() >= start && m.getPosition() <= end)
								if (m.inside(m.pPlus) > 0) // TODO: allow any scores, now 0 vs >0
									++e1;
								else
									++e0;
							
				} else {
					//System.err.println(e1 + "\t" + e0);
					ArrayList<Marker> markersCi = new ArrayList<Marker>();			
					for (int  map = 0; map < numMaps; ++map)
						for (Marker m : ci.getMarkers(map))
							markersCi.add(m);
					Collections.sort(markersCi);
					//calculateErrors(markersCi, errors, e0, e1);
					calculateErrors2(markersCi, errors, e0, e1, numErrorsPerContig);
				}
			}
		}

		if (verbose) {
			System.err.println("#*** possible errors ***");
			System.err.println("#total outside/inside:\t" + e0 + "\t" + e1);
			System.err.println("#X2\t0\t>0\tcontig\tstart\tend\tcontig2\tpos");
		}

		Collections.sort(errors);
		if (errors != null)
			for (int i = 0; i < numErrors && i < errors.size(); ++i) {
				ret.add(errors.get(i));
				PossibleError ff = findFix(cis, errors.get(i), false); 
				ret.add(ff);
				if (verbose) {
					System.err.println(errors.get(i).getInfo() + ((ff == null) ? "" : "\t" + ff.getInfo()) + "\terror");
				}
			}
		return ret;
	}
	
	//TODO: Check that calculateChainScore is ok
	//note: assumes that score has been evaluated for the markers in cis and in the corresponding order and orientation
	private void findLikelyAssemblyErrors(ArrayList<ContigInterval> cis)
	{
		System.err.println("*** possible haplotypes ***");
		System.err.println("score\tcontig\tstart\tend\tof_contig\tstart\tend\tcontig_aligment_start\tcontig_aligment_end");
		
		ArrayList<PossibleError> haplotypes = new ArrayList<PossibleError>(); 
		
		for (int cii = 0; cii < cis.size(); ++cii) {
			ContigInterval ci = cis.get(cii);
			int score = 0;
			for (int  map = 0; map < numMaps; ++map)
				for (Marker m : ci.getMarkers(map))
					score += m.inside(m.pPlus);
			if (cii > 0)
				score += calculateChainScore(cis.get(cii - 1), ci);
			if (cii + 1  < cis.size())
				score += calculateChainScore(ci, cis.get(cii + 1));

			//int scoreH = 0;
			//ContigInterval haplotypeOf = null;
			for (ContigInterval ci2 : cis) {
				Integer s = calculateChainScoreHaplotype(ci2, ci);
				if (s != null && s > score) {
					long ha[] = calculateChainScoreHaplotypeCut(ci2, ci);
					haplotypes.add(new PossibleError((s - score), 0, (s - score) + "\t" + ci + "\t" + ci2 + "\t" + ha[0] + "\t" + ha[1]));
//					scoreH = s;
//					haplotypeOf = ci2;
				}
			}
//			if (scoreH > score) {
//				haplotypes.add(new PossibleError((scoreH - score), (scoreH - score) + "\t" + ci + "\t" + haplotypeOf));
//				//System.err.println((scoreH - score) + "\t" + ci + "\t" + haplotypeOf); 
//			}

		}
		Collections.sort(haplotypes);
		for (PossibleError pe: haplotypes)
			System.err.println(pe.getInfo() + "\thaplotype");
		
		findNonHaplotypeAssemblyErrors(cis, true);
		
	}
	void setMaxIntersect(int parseInt) {
		maxIntersect = parseInt;
	}

	void setMinLinkAlignmentScore(int parseInt) {
		minLinkAlignmentScore = parseInt;
		
	}

	private void setMinHaplotypeAlignmentScore(int parseInt) {
		minHaplotypeAlignmentScore = parseInt;		
	}

	//autogenerated stubs...
	void setMaxBridge(int parseInt) {
		maxBridge = parseInt;		
	}

	void setCutPenalty(double parseDouble) {
		cutPenalty = parseDouble;
	}

	void setOrientationPenalty(double parseDouble) {
		orientationPenalty = parseDouble;
	}

	void setScaleScore(double parseDouble) {
		scaleScore = parseDouble;
	}
	void setNumThreads(int nt) {
		numThreads = nt;
	}

	private ArrayList<ArrayList<String>> myLiftover(ArrayList<ArrayList<String>> markers, ArrayList<long[]> alignment, boolean sameOrientation, String newContig)
	{
		ArrayList<ArrayList<String>> ret = new ArrayList<ArrayList<String>>();
		for (ArrayList<String> row: markers) {
			//System.err.println(row);
			long position = InputData.myParseLong(row.get(1));
			long position_new = mapPosition12(alignment, position, sameOrientation);
			//System.err.println(position + "->" + position_new);
			if (position_new > 0) {
				//System.err.println("HIPHEI");
				ArrayList<String> nrow = new ArrayList<String>(); 
				for (int i = 0; i < row.size(); ++i)
					if (i == 0)
						nrow.add(newContig);
					else if (i == 1)
						nrow.add("" + position_new);
					else
						nrow.add(row.get(i));
				ret.add(nrow);
			}
		}
		return ret;
	}
	
	private void liftover__(ContigInterval hap, ContigInterval c2, ArrayList<long[]> alignment, boolean sameOrientation, long ascore, HashMap<ContigInterval, Long> bestScore, HashMap<ContigInterval, ContigInterval> bestScoreContig, HashMap<ContigInterval, ArrayList<ArrayList<String>>> map, HashMap<ContigInterval, ArrayList<ArrayList<String>>> liftoverMap) {
		Long prevScore = bestScore.get(hap);
		if (prevScore != null) {
			ContigInterval c3 = bestScoreContig.get(hap);
			if (c3.getContig().equals(c2.getContig()) && Misc.intersectIntervals(c2.getStart(), c2.getEnd(), c3.getStart(), c3.getEnd())) { // same alignment...
				if (ascore > prevScore) {
					liftoverMap.put(hap, myLiftover(map.get(hap), alignment, sameOrientation, c2.getContig())); // do liftover
					bestScoreContig.put(hap, c2);
					bestScore.put(hap, ascore);					
				}
			} else {
				if (ascore > prevScore) {
					bestScoreContig.put(hap, c2);
					bestScore.put(hap, ascore);					
				}
				if (ascore > 2 * prevScore) {
					liftoverMap.put(hap, myLiftover(map.get(hap), alignment, sameOrientation, c2.getContig())); // do liftover
				} 
				else if (ascore > liftoverScoreDiff * prevScore) {// too little score difference...
					liftoverMap.get(hap).clear(); // remove markers...
					//System.err.println("remove");
				}
			}
			
		} else {
			//System.err.println("liftover");
			liftoverMap.put(hap, myLiftover(map.get(hap), alignment, sameOrientation, c2.getContig())); // do liftover
			bestScore.put(hap, ascore);
			bestScoreContig.put(hap, c2);
		}
	}
	
	private void liftover_(String fn, HashMap<ContigInterval, ArrayList<ArrayList<String>>> map) {
		System.err.println("loading alignment chain...");
		
		HashMap<ContigInterval, Long> bestScore = new HashMap<ContigInterval, Long>(); // score of used chain   
		HashMap<ContigInterval, ContigInterval> bestScoreContig = new HashMap<ContigInterval, ContigInterval>(); // contigInterval of used chain   
		
		HashMap<ContigInterval, ArrayList<ArrayList<String>>> liftoverMap = new HashMap<ContigInterval, ArrayList<ArrayList<String>>>(); // store liftover markers...
		
		for (ContigInterval ci : map.keySet()) {
			liftoverMap.put(ci, new ArrayList<ArrayList<String>>());
			//System.err.println(ci);
		}
		
		try {
			BufferedReader br = null;
			if (fn.equals("-"))
				br = new BufferedReader(new InputStreamReader(System.in));
			else
				br = new BufferedReader(new FileReader(fn));
			do {
				ArrayList<String> row = Input.loadTableRow(br, "\t ");
				if (row == null)
					break;
				boolean skip = true;
				if (row.size() >= 12 && "chain".equals(row.get(0))) {
					String contig1 = row.get(2);
					String contig2 = row.get(7);
					
					boolean someMarkers1 = false;
					if (haplotypeHash.containsKey(contig1))
						for (ContigInterval hap : haplotypeHash.get(contig1))
							if (map.containsKey(hap)) {
								someMarkers1 = true;
								break;
							}
					boolean someMarkers2 = false;
					if (haplotypeHash.containsKey(contig2))
						for (ContigInterval hap : haplotypeHash.get(contig2)) {
							if (map.containsKey(hap)) {
								someMarkers2 = true;
								break;
							}
						}
					if (!someMarkers1 && !someMarkers2) // no markers to do liftover...
						continue;

					skip = false;
					long p1_[] = chain2OneBase(row.get(3), row.get(4), row.get(5), row.get(6));
					long p2_[] = chain2OneBase(row.get(8), row.get(9), row.get(10), row.get(11));
					if ("-".equals(row.get(4))) {
						System.err.println("Error: only ++, and +- orientation allowed in the chain file");
						continue;
					}
					
					//System.err.println("running");
					
					boolean sameOrientation = row.get(4).equals(row.get(9));

					// load alignment for liftover...
					ArrayList<long[]> alignment_ = new ArrayList<long[]>();
					ArrayList<String> row2 = null;
					long v1 = p1_[0];
					long v2 = ((sameOrientation) ? p2_[0] : p2_[1]);
					int v3 = 0;
					do {
						row2 = Input.loadTableRow(br, "\t ");
						
						v3 = Integer.parseInt(row2.get(0));
						alignment_.add(new long[]{v1, v2, v3});
						v1 += v3;
						v2 += ((sameOrientation) ? v3 : -v3);
						if (row2.size() >= 3) {
							v1 += Integer.parseInt(row2.get(1));
							v2 += ((sameOrientation) ? Integer.parseInt(row2.get(2)) : -Integer.parseInt(row2.get(2)));
						}
					} while (row2 != null && row2.size() >= 3);
					
					long ascore_ = Long.parseLong(row.get(1));

					if (someMarkers1) //  && haplotypeHash.containsKey(contig1)
						for (ContigInterval hap : haplotypeHash.get(contig1)) 
							if (map.containsKey(hap)) {
								ContigInterval c2 = new ContigInterval(contig2, p2_[0], p2_[1]);
								
								ArrayList<long[]> alignment = new ArrayList<long[]>();
								double t = trimChain(alignment_, sameOrientation, hap, c2, alignment);
								if (t == 0.0) // all trimmed
									continue;

								if (t < 1.0) { // update c2 start and end
									long p[] = getStartEnd2(alignment, sameOrientation);
									c2 = new ContigInterval(contig2, p[0], p[1]);
								} else {
									alignment = alignment_;
								}

								long ascore = (long)(t * ascore_);
								
								boolean doLift = true;
								if (haplotypeHash.containsKey(contig2))
									for (ContigInterval hap2 : haplotypeHash.get(contig2)) // no liftover between haplotypes...
										if (Misc.intersectIntervals(hap2.getStart(), hap2.getEnd(), c2.getStart(), c2.getEnd())) {
											doLift = false;
											break;
										}
								if (doLift)
									liftover__(hap, c2, alignment, sameOrientation, ascore, bestScore, bestScoreContig, map, liftoverMap);
							}
					//otherway round...
					if (someMarkers2) { //  && haplotypeHash.containsKey(contig2)
						ArrayList<long[]> revAlignment_ = reverseAlignment(alignment_, sameOrientation);
						for (ContigInterval hap : haplotypeHash.get(contig2))
							if (map.containsKey(hap)) {
								ContigInterval c2 = new ContigInterval(contig1, p1_[0], p1_[1]);
								
								ArrayList<long[]> alignment = new ArrayList<long[]>();
								double t = trimChain(revAlignment_, sameOrientation, hap, c2, alignment);
								if (t == 0.0) // all trimmed
									continue;
								if (t < 1.0) { // update c2 start and end
									long p[] = getStartEnd2(alignment, sameOrientation);			 						
									c2 = new ContigInterval(contig1, p[0], p[1]);
								} else
									alignment = revAlignment_; 
								long ascore = (long)(t * ascore_);

								boolean doLift = true;
								if (haplotypeHash.containsKey(contig1))
									for (ContigInterval hap2 : haplotypeHash.get(contig1)) // no liftover between haplotypes...
										if (Misc.intersectIntervals(hap2.getStart(), hap2.getEnd(), c2.getStart(), c2.getEnd())) {
											doLift = false;
											break;
										}
								if (doLift)
									liftover__(hap, c2, alignment, sameOrientation, ascore, bestScore, bestScoreContig, map, liftoverMap);
							}
					}
				} else {
					if (skip && (row.size() == 3 || row.size() == 1))
						;
					else
						System.err.println("Warning: skipping " + row);
				}
			} while (true);
			br.close();
			int numLifoverMarkers = 0;
			for (ContigInterval ci : liftoverMap.keySet()) {
				ArrayList<ArrayList<String>> lom = liftoverMap.get(ci);
				numLifoverMarkers += lom.size();
				StringBuilder sb = new StringBuilder(); 
				for (ArrayList<String> row : lom) {
					sb.append(row.get(0));
					for (int i = 1; i < row.size(); ++i) {
						sb.append('\t');
						sb.append(row.get(i));
					}
					sb.append('\n');
				}
				System.out.print(sb);
			}
			System.err.println("Lifting over " + numLifoverMarkers + " markers");
		} catch (Exception e) {
	     	e.printStackTrace();
			System.err.println("Error in file " + fn);
      	}
	}

	public void liftover(String haplotypeFile, String chainFile, String mapFile) {
		ArrayList<ContigInterval> haplotypes = InputData.loadHaplotypes(haplotypeFile);
		
		for (ContigInterval hapci: haplotypes) {
			String key = hapci.getContig(); 
			ArrayList<ContigInterval> list = haplotypeHash.get(key);
			if (list == null) {
				list = new ArrayList<ContigInterval>(); 
				haplotypeHash.put(key, list);
			}
			list.add(hapci);
		}

		HashMap<ContigInterval, ArrayList<ArrayList<String>>> map = InputData.loadRaw(mapFile, haplotypeHash);
		//System.out.println(map);
		
 		liftover_(chainFile, map);
		
	}

	private void loadCutSites(String fn) {
		try {
			HashMap<String, ArrayList<ContigInterval>> bedHash2 = new HashMap<String, ArrayList<ContigInterval>>();
			for (ContigInterval ci : bed) {
				if (ci.getStart() != ci.getMinStart() || ci.getEnd() != ci.getMaxEnd()) { // ContigInterval with uncertainty in its start or end
					String contig = ci.getContig();
					if (!bedHash2.containsKey(contig))
						bedHash2.put(contig, new ArrayList<ContigInterval>());
					bedHash2.get(contig).add(ci);
				}
			}
			int numCutSites = 0;
			BufferedReader	br = new BufferedReader(new FileReader(fn));
			do {
				ArrayList<String> row = Input.loadTableRow(br, "\t ");
				if (row == null)
					break;
				if (row.size() >= 2) {
					String contig = row.get(0);
					long position = 0;
					long position2 = 0;
					if (bedHash2.containsKey(contig))
						for (ContigInterval ci : bedHash2.get(contig)) {
							if (position == 0) {
								position = Long.parseLong(row.get(1));
								if (row.size() >= 3) {
									position2 = Long.parseLong(row.get(2));
									--position; //N gap 
									++position2;
								} else
									position2 = position + 1;
							}
							
							
							if (position2 >= ci.getMinStart() && position2 <= ci.getMaxStart()) {
								++numCutSites;
								ci.setStart(position2);
								break;
							}
							if (position >= ci.getMinEnd() && position <= ci.getMaxEnd()) {
								++numCutSites;
								ci.setEnd(position);
								break;
							}
						}
				}
			} while (true);
			br.close();
			System.err.println("Trimming " + numCutSites + " contig ends based on cutSites file");
		} catch (Exception e) {
	     	e.printStackTrace();
			System.err.println("Error in file " + fn);
      	}
	}
	
	//find how to fix likely erroneous genomic interval (e)
	//assumes solve has been called on the cis 
	private PossibleError findFix(ArrayList<ContigInterval> cis, PossibleError e, boolean verbose){
		//update rank and numContigIntervals
		int numContigIntervals = 0; // update rank (might be corrupted)
		for (ContigInterval c : cis)
			c.setRank(numContigIntervals++);

		//calculate full anchoring score
		ArrayList<Integer> oldPos = new ArrayList<Integer>(); 
		int oldFullScore = 0;
		for (int map = 0; map < numMaps; ++map)
			for (Marker m : getMarkers(cis, null, map)) {
				oldFullScore += m.inside(m.pPlus);
				oldPos.add(m.pPlus); // store mPlus, as this is changed...
			}

		//get error region (region inside one contigInterval)
		ContigInterval errorContig = cis.get(e.getCi());
		long pos1[] = e.getPos1();
		long pos2[] = e.getPos2();
	
		int s0 = 0; // old score for region e
		int s1 = 0; // score for e in + orient
		int s2 = 0; // score for e in - orient
		
		ArrayList<ArrayList<Integer>> iss = new ArrayList<ArrayList<Integer>>(); //store increaseScores... 

		for (int map = 0; map < numMaps; ++map) {
			//get markers spanning error region for each map separately 

			ArrayList<Marker> errorMarkers = new ArrayList<Marker>();
			ArrayList<Marker> nonErrorMarkers = new ArrayList<Marker>();
			
			for (ContigInterval c : cis)
				for (Marker m : c.getMarkers(map)) {
					if ((c == errorContig) && m.getPosition() >= pos1[0] && m.getPosition() <= pos2[1]) 
						errorMarkers.add(m);
					else
						nonErrorMarkers.add(m);
				}

			int numMarkers = errorMarkers.size();
			//calculate old score
			//int oldScore[] = new int[numMarkers];
			for (int i = 0; i < numMarkers; ++i) {
				Marker m = errorMarkers.get(i);
				int s = m.inside(m.pPlus);
				s0 += s;
			}
			s1 += solveForward(errorMarkers);
			s2 += solveBackward(errorMarkers);

			ArrayList<Marker> mm2 = new ArrayList<Marker>();
			for (int i = 0; i < numMarkers; ++i) {
				Marker m = errorMarkers.get(i);
				if (m.inside(m.pPlus) > 0) // marker is now contributing to the score
					mm2.add(m);
			}
			int maxBin = getMaxBin(errorContig.getMarkers(map));
			
			iss.add(increaseScore(nonErrorMarkers, calculateBinSum(mm2, maxBin, true))); // mm2 in + orientation
			
			mm2.clear();
			for (int i = 0; i < numMarkers; ++i) {
				Marker m = errorMarkers.get(i);
				if (m.inside(m.pMinus) > 0) // marker is now contributing to the score
					mm2.add(m);
			}
			
			iss.add(increaseScore(nonErrorMarkers, calculateBinSum(mm2, maxBin, false))); // mm2 in - orientation
			
			//System.err.println("size=" + nonErrorMarkers.size());

		}

		int pi = 0;
		for (int map = 0; map < numMaps; ++map)
			for (Marker m : getMarkers(cis, null, map))
				m.pPlus = oldPos.get(pi++); // restore mPlus, this is changed by increaseScore
		
		if ((s0 >= s1 && s0 >= s2))
			return null; // no fix found...
				
		//System.err.println("max score delta = " + (s1 - s0) + " or " + (s2 - s0) + "\t" + s0 + "\t" + s1 + "\t" + s2 + "\t" + e.getInfo());
		//merge results from each increaseScore

		long maxScore = oldFullScore;
		int maxC = 0;
		long maxP1 = 0;
		long maxP2 = 0;
		
		for (int orient = 0; orient < 2; ++orient) {
			int startMarker[] = new int[numMaps];
			for (ContigInterval c : cis) {
				//store cut intervals
				//if interval [a,b] is in cov, then prefix can be [1..a], [1..a+1],...,[1..b]
				ArrayList<long []> cov = new ArrayList<long []>();
			
				for (int map = 0; map < numMaps; ++map) {
					ArrayList<Integer> is = iss.get(map + map + orient);
					
					//System.err.println("size=" + is.size());
					
					ArrayList<Marker> markers = new ArrayList<Marker>(); 
					for (Marker m : c.getMarkers(map))
						if ((errorContig == c) && m.getPosition() >= pos1[0] && m.getPosition() <= pos2[1]) // errorMarker
							;
						else
							markers.add(m);

					//System.err.println("size=" + markers.size());
					
					int numMarkers = markers.size();
					if (numMarkers > 0) {
						long isi = is.get(startMarker[map]);
						if (c.getOrientation() >= 0) { // + orientation
							cov.add(new long[]{c.getStart() - 1, markers.get(0).getPosition() - 1, isi});
							for (int mi = 1; mi < numMarkers; ++mi) {
								//Marker m = markers.get(mi);
								isi = is.get(startMarker[map] + mi);
								cov.add(new long[]{markers.get(mi - 1).getPosition(), markers.get(mi).getPosition() - 1, isi});
							}
							isi = is.get(startMarker[map] + numMarkers);
							cov.add(new long[]{markers.get(numMarkers - 1).getPosition(), c.getEnd(), isi});
						} else {
							cov.add(new long[]{markers.get(0).getPosition(), c.getEnd(), isi});
							for (int mi = 1; mi < numMarkers; ++mi) {
								//Marker m = markers.get(mi);
								isi = is.get(startMarker[map] + mi);
								cov.add(new long[]{markers.get(mi).getPosition(), markers.get(mi - 1).getPosition() - 1, isi});
							}
							isi = is.get(startMarker[map] + numMarkers);
							cov.add(new long[]{c.getStart() - 1, markers.get(numMarkers - 1).getPosition() - 1, isi});
						}
					} else { // 0 markers
						int isi = is.get(startMarker[map]);
						cov.add(new long[]{c.getStart() - 1, c.getEnd(), isi});
					}
					startMarker[map] += numMarkers;
				}
				ArrayList<Long> cov_ret = Misc.cov(cov);
				//System.err.println(cov_ret);
				
				for (int i = 0; i < cov_ret.size(); i+=2) {
					long pos = cov_ret.get(i);
					long cr = cov_ret.get(i + 1);
					//System.err.println(pos + "\t" + cr);
					if (cr > maxScore) {
						maxScore = cr;
						maxC = c.getRank();
						maxP1 = pos;
						maxP2 = ((i + 2 < cov_ret.size()) ? cov_ret.get(i + 2) - 1 : maxP1); 
					}
				}
			}
		}

		if (verbose)
			System.err.println(cis.get(maxC).getContig() + "\t" + maxP1 + "-" + maxP2 + "\t" + errorContig.getContig() + "\t" + pos1[0] + "-" + pos1[1] + "\t" + pos2[0] + "-" + pos2[1] + "\t" + (maxScore - oldFullScore) + "\tfix");
		
		String info = cis.get(maxC).getContig() + "\t" + maxP1 + "-" + maxP2;

		if (maxScore <= oldFullScore) {
			return null;
		}
		return new PossibleError((maxScore - oldFullScore), maxC, info , new long[]{maxP1, maxP2}, null);		
	}

	// find a region inside an error with even higher error rate 
	private PossibleError splitError(ContigInterval c, PossibleError e)
	{
		//PossibleError ret = new PossibleError(max, ci.getRank(), info , new long[]{pos1n, pos1}, new long[]{pos2, pos2n}
		long pos1[] = e.getPos1();
		long pos2[] = e.getPos2();
		ArrayList<Marker> markersC = new ArrayList<Marker>();
		int e0 = 0;
		int e1 = 0;
		for (int  map = 0; map < numMaps; ++map)
			for (Marker m : c.getMarkers(map))
				if (m.getPosition() >= pos1[0] && m.getPosition() <= pos2[1]) {
					if (m.inside(m.pPlus) > 0)
						++e1;
					else
						++e0;
					markersC.add(m);
				}
		Collections.sort(markersC);
		//System.err.println("e0=" + e0 + "\te1=" + e1);
		
		//calculateErrors(markersCi, errors, e0, e1);
		ArrayList<PossibleError> pes = new ArrayList<PossibleError>(); 
		calculateErrors2(markersC, pes, e0, e1, 1);
		if (pes.size() > 0) {
			PossibleError ret = pes.get(0);
			long p1[] = ret.getPos1(); // adjust region so that it is within the original error...
			if (p1[0] < pos1[0])
				p1[0] = pos1[0];
			long p2[] = ret.getPos2();
			if (p2[1] > pos2[1])
				p2[1] = pos2[1];
			ret.ci = e.ci; // rank might be corrupted so calculateErrors2 could fail to get this right
			return ret;
		}

		return null;
	}
	//iterative version of findLikelyAssemblyErrors
	private void findContigErrors(int minImprovement) {
		System.err.println("Finding contig errors...");
		ArrayList<ContigInterval> newcis = new ArrayList<ContigInterval>();
		
		int animation = 0;
		boolean foundFix = true;
		
		HashMap<String, Double> missedErrors = new HashMap<String, Double>();
		
		while (foundFix) {
			foundFix = false;
			int scoreOld = calculateScore(intervalsInChr, false) + calculateNonMarkerScore(intervalsInChr);

			int orientation[] = new int[intervalsInChr.size()]; // store orientaions
			for (int ci = 0; ci < orientation.length; ++ci)
				orientation[ci] = intervalsInChr.get(ci).getOrientation();
			
			//ArrayList<PossibleError> pe = findNonHaplotypeAssemblyErrors(intervalsInChr, false);
			ArrayList<PossibleError> pe = findNonHaplotypeAssemblyErrors(intervalsInChr, false);
			
			System.err.println("score = " + scoreOld);

		out:for (int e = 0; e + e < pe.size(); ++e) {
				PossibleError pe1 = pe.get(e + e);
				ContigInterval c1 = intervalsInChr.get(pe1.getCi());
				
				PossibleError pe2 = pe.get(e + e + 1);
				String errorId = pe1.getInfo() + ((pe2 == null) ? "" : "\t" + pe2.getInfo());
				
				// see if we have already tried to fix this error, skip with decreasing prob
				String errorId2 = errorId.substring(errorId.indexOf('\t'));
				if (missedErrors.containsKey(errorId2)) {
					double newValue = missedErrors.get(errorId2) * 0.5;
					if (Math.random() >= newValue)
						continue;
					missedErrors.put(errorId2, newValue);
				} else
					missedErrors.put(errorId2, 2.0); // try twice, then with pr 1/2, 1/4, 1/8, ...
				
				System.err.println(errorId);
				for (int fixIteration = 0; fixIteration < 2; ++fixIteration) { // try to split each region in fixIteration==1
					if (fixIteration == 1) {
						calculateScore(intervalsInChr, false); // update pPlus...
						pe1 = splitError(c1, pe1);
						if (pe1 == null) {
							//System.err.println("no split");
							break;
						}
						pe2 = findFix(intervalsInChr, pe1, false);
						errorId = pe1.getInfo() + ((pe2 == null) ? "" : "\t" + pe2.getInfo());
						System.err.println(errorId);
					}
					ContigInterval sp1[] = null;
					ContigInterval sp2[] = null;
					ContigInterval c2 = null;
					
					long p11[] = pe1.getPos1();
					p11[0] -= 1; //-1
					p11[1] -= 1; //-1
					if (pe2 != null) {
						c2 = intervalsInChr.get(pe2.getCi());
						if (c1 == c2) {
							sp1 = c1.splitContigInterval(p11, pe1.getPos2(), pe2.getPos1());
						} else {
							sp1 = c1.splitContigInterval(p11, pe1.getPos2());
							sp2 = c2.splitContigInterval(pe2.getPos1());
						}
					} else {
						sp1 = c1.splitContigInterval(p11, pe1.getPos2());
					}
		
					newcis.clear();
					for (ContigInterval ci : intervalsInChr) {
						if (ci != c1 && ci != c2)
							newcis.add(ci);
						else if (ci == c1)
							for (ContigInterval c : sp1)
								newcis.add(c);
						else
							for (ContigInterval c : sp2)
								newcis.add(c);
					}
					//System.err.println(newcis);
					int splitSize1 = sp1.length;
					int splitSize2 = ((sp2 == null) ? 1 : sp2.length);
					if (splitSize1 == 1 && splitSize2 == 1) // no split...
						continue;
					
	
					//add scaffolding links
					for (int ci = 1; ci < sp1.length; ++ci) {
						String key1 = calculateKey(sp1[ci - 1], true, sp1[ci],     true);
						String key2 = calculateKey(sp1[ci],    false, sp1[ci - 1], false);
						scaffoldingLink.put(key1, 1);
						scaffoldingLink.put(key2, 1);
					}
					if (sp2 != null)
						for (int ci = 1; ci < sp2.length; ++ci) {
							String key1 = calculateKey(sp2[ci - 1], true, sp2[ci],      true);
							String key2 = calculateKey(sp2[ci],    false, sp2[ci - 1], false);
							scaffoldingLink.put(key1, 1);
							scaffoldingLink.put(key2, 1);
						}
					//keep old scaffoldingLinks
					for (ContigInterval c : intervalsInChr) {
						String k1 = calculateKey(c1, true, c, true);
						if (scaffoldingLink.containsKey(k1)) {
							scaffoldingLink.put(calculateKey(sp1[sp1.length - 1], true, c,true), 1);
							scaffoldingLink.put(calculateKey(c, false, sp1[sp1.length - 1], false), 1);
						}
						String k2 = calculateKey(c, true, c1, true);
						if (scaffoldingLink.containsKey(k2)) {
							scaffoldingLink.put(calculateKey(c, true, sp1[0], true), 1);
							scaffoldingLink.put(calculateKey(sp1[0], false, c, false), 1);
						}
						if (sp2 != null) {
							k1 = calculateKey(c2, true, c, true);
							if (scaffoldingLink.containsKey(k1)) {
								scaffoldingLink.put(calculateKey(sp2[sp2.length - 1], true, c,true), 1);
								scaffoldingLink.put(calculateKey(c, false, sp2[sp2.length - 1], false), 1);
							}
							k2 = calculateKey(c, true, c2, true);
							if (scaffoldingLink.containsKey(k2)) {
								scaffoldingLink.put(calculateKey(c, true, sp2[0], true), 1);
								scaffoldingLink.put(calculateKey(sp2[0], false, c, false), 1);
							}
						}
					}
					//keep old chainLinks
					for (ContigInterval c : intervalsInChr) {
						for (int o = 0; o < 2; ++o)
							for (int o1 = 0; o1 < 2; ++o1) {
								String k1 = calculateKey(c, o == 0, c1, o1 == 0);
								if (chainLinkHash.containsKey(k1)) {
									String knew = calculateKey(c, o == 0, (o1==0) ? sp1[0] : sp1[sp1.length - 1], o1 == 0);
									chainLinkHash.put(knew, chainLinkHash.get(k1));
								}
								String k2 = calculateKey(c1, o1 == 0, c, o == 0);
								if (chainLinkHash.containsKey(k2)) {
									String knew = calculateKey((o1==0) ? sp1[sp1.length - 1]: sp1[0], o1 == 0, c, o == 0);
									chainLinkHash.put(knew, chainLinkHash.get(k2));
								}
								if (sp2 != null) {
									String k3 = calculateKey(c, o == 0, c2, o1 == 0);
									if (chainLinkHash.containsKey(k3)) {
										String knew = calculateKey(c, o == 0, (o1==0) ? sp2[0] : sp2[sp2.length - 1], o1 == 0);
										chainLinkHash.put(knew, chainLinkHash.get(k3));
									}
									String k4 = calculateKey(c2, o1 == 0, c, o == 0);
									if (chainLinkHash.containsKey(k4)) {
										String knew = calculateKey((o1==0) ? sp2[sp2.length - 1]: sp2[0], o1 == 0, c, o == 0);
										chainLinkHash.put(knew, chainLinkHash.get(k4));
									}
									
								}

							}
					}
					if (sp2 != null) // keep chain links between c1 to c2
						for (int o1 = 0; o1 < 2; ++o1)
							for (int o2 = 0; o2 < 2; ++o2) {
								String k1 = calculateKey(c1, o1 == 0, c2, o2 == 0); 
								if (chainLinkHash.containsKey(k1)) {
									String knew = calculateKey((o1==0) ? sp1[sp1.length - 1]: sp1[0], o1 == 0, (o2==0) ? sp2[0]: sp2[sp2.length - 1], o2 == 0);
									chainLinkHash.put(knew, chainLinkHash.get(k1));
								}
								String k2 = calculateKey(c2, o2 == 0, c1, o1 == 0); 
								if (chainLinkHash.containsKey(k2)) {
									String knew = calculateKey((o2==0) ? sp2[sp2.length - 1]: sp2[0], o2 == 0, (o1==0) ? sp1[0]: sp1[sp1.length - 1], o1 == 0);
									chainLinkHash.put(knew, chainLinkHash.get(k2));
								}
							}
					
					
							
					//improveAnchoring(newcis, false, false);
	
					//faster version of improveAnchoring 
					{
						int bestScore = -1;
						ArrayList<ContigInterval> newcontigs = new ArrayList<ContigInterval>();
						for (ContigInterval c : sp1) 
							newcontigs.add(c);
						if (sp2 != null)
							for (ContigInterval c : sp2)
								newcontigs.add(c);
						ArrayList<ContigInterval> newcis2 = new ArrayList<ContigInterval>();
						newcis2.addAll(newcis);
						
						ArrayList<Integer> bestOrientation = new ArrayList<Integer>(); 
	
						for (int run = 0; run < numRuns; ++run) {
							Collections.shuffle(newcontigs);
							for (ContigInterval c : newcontigs) {
								c.setOrientation((int)(2.0 * Math.random()) - 1);
							}
							boolean foundBetter = true;
							
							while (foundBetter) {
								foundBetter = false;
								for (ContigInterval c : newcontigs) {
									if (calculateBest(newcis2, c, false) > 0)
										foundBetter = true;
								}
							}
							int score = calculateScore(newcis2, false) + calculateNonMarkerScore(newcis2);
							if (score > bestScore) {
								bestScore = score;
	
								newcis.clear();
								newcis.addAll(newcis2);
								
								bestOrientation.clear();
								for (ContigInterval c : newcis)
									bestOrientation.add(c.getOrientation());
							}
						}
						int oi = 0;
						for (ContigInterval c : newcis)
							c.setOrientation(bestOrientation.get(oi++));
					}
	
					int score = calculateScore(newcis, false) + calculateNonMarkerScore(newcis);
					//System.err.println(score);
					if (score >= minImprovement + scoreOld + 2 * (splitSize1 + splitSize2 - 2)) { // at least +minImprovement in score and +2 for each cut...
						System.err.println("split into " + splitSize1 + (splitSize2 <= 1 ? "" : "+" + splitSize2));
						System.err.print("number of markers");
						for (ContigInterval c : sp1)
							System.err.print("\t" + c.getNumMarkers());
						if (sp2 != null) {
							System.err.print("\t|");
							for (ContigInterval c : sp2)
								System.err.print("\t" + c.getNumMarkers());
						}
						System.err.println();
						
						if (!printAnimation.equals("")) {
							if (animation == 0) { 
								String fn = printAnimation + animation + ".la";
								
								ArrayList<Double> orderSupport = calculateSupport(intervalsInChr);
								try {
									PrintStream ps = new PrintStream(fn);
									printAnchoring(ps, intervalsInChr, orderSupport);
									ps.close();
								}
								catch (Exception ex){
									ex.printStackTrace();
								}
							}
							++animation;							
							String fn = printAnimation + animation + ".la";
							ArrayList<Double> orderSupport = calculateSupport(newcis);
							try {
								PrintStream ps = new PrintStream(fn);
								printAnchoring(ps, newcis, orderSupport);
								ps.close();
							}
							catch (Exception ex){
								ex.printStackTrace();
							}
						}
	
						System.err.println("Score improvement = " + (score - scoreOld) + " score = " + score);
						//System.err.println("score = " + score);
						intervalsInChr.clear();
						intervalsInChr.addAll(newcis);
						foundFix = true;
						
						break out; // found fix...
					}
					// set orientations back to stored state
					for (int ci = 0; ci < orientation.length; ++ci)
						intervalsInChr.get(ci).setOrientation((orientation[ci]));
				} // fixIteration
			} // for (int inr e=0; ...
		} // while foundFix
		int score = calculateScore(intervalsInChr, false);
		System.out.println("#final score with corrected contigs " +  score + " " + calculateNonMarkerScore(intervalsInChr));
		ArrayList<Double> orderSupport = calculateSupport(intervalsInChr);

		printAnchoring(intervalsInChr, orderSupport);

		//print bed for another run of Lep-Anchor
		ContigInterval start = null;
		for (int ci = 0; ci < intervalsInChr.size(); ++ci) {
			ContigInterval c = intervalsInChr.get(ci);
			ContigInterval next = ((ci + 1 < intervalsInChr.size()) ? intervalsInChr.get(ci + 1) : null);
			boolean linked = false;
			if (next != null) { 						 
				String key = calculateKey(c, next);
				if (scaffoldingLink.containsKey(key) && scaffoldingLink.get(key) > 0)
					linked = true;
			}
			if (start == null)
				start = c;
			
			if (!linked) {
				long pos1[] = c.getStartI();
				long pos2[] = c.getEndI();
				if (c != start) // more than one linked contigsIntervals...
					if (c.getOrientation() >= 0) {
						pos1 = start.getStartI();
					} else {
						pos2 = start.getEndI();
					}
				System.err.println(c.getContig() + "\t" + pos1[0] + "-" + pos1[1] + "\t" + pos2[0]+ "-" + pos2[1] + ((pos2.length==2) ? "" : "*") + "\t?\t" + c.getChromosome() + "\tbed");

				start = null;
			}
 		}
	}
	
	private static void usageInfo()
	{
        System.err.println("usage: java PlaceAndOrientContigs bed=bed.file map=map1.txt [map2 [map3 ...]] options");
        System.err.println("       bed=file                    a file containing (contig start stop) intervals in 1-based coordinates");
        System.err.println("       map=file1 [file2 [...]]     linkage map file(s)");
        System.err.println("       columns contig, pos, chromosome, map_bin_start [map_bin_stop [map_bin_start2 map_bin_stop2] [...]]");
        System.err.println("       orientation=+/- [+/- [...]] manual orientation for each map file (found automatically by default)");
        System.err.println("       chain=file                  chain file ");
        System.err.println("       noChromosome=1              input file does not have chromosome column ");
        System.err.println("       noIntervals=1               input file does not have intervals but map positions");
        System.err.println("       numRuns=NUM                 run this many runs to find better anchoring [5]");
        System.err.println("       chromosome=NUM              take only chromosome NUM from the bed (chr=column 5) [not set]");
        System.err.println("                                   and from the map(s)");
        System.err.println("       numThreads=NUM              number of threads [1]");

        System.err.println("       compressMap=0	           Do not compress map positions [1]");

        System.err.println("       randomOrder=1               Start with a random anchoring [not set]");
        
        System.err.println("       keepEmptyIntervals=1        Keep (contig)intervals without any markers [not set]");
        
        System.err.println("       numErrors=NUM               List at most this many potential errors [40]");
        System.err.println("       numErrorsPerContig=NUM      List at most this many potential errors from one contig [3]\n");
        
        System.err.println("       paf=file                    load alignment file in paf (minimap2) format");
        System.err.println("       maxBridge=NUM               maximum scaffolding bridge length (for paf input) [50000]");
        System.err.println("       scalePaf=NUM                multiply scaffolding links by NUM [1]");
        System.err.println("       maxIntersect=NUM            maximum alignment intersection from paf [2000]");
        System.err.println("       maxPafScore=NUM             maximum link score from the paf [not set]\n");

        System.err.println("       scaleScore=NUM              scale aligment scores to markers [0.00001] (0.00001 = 1kb 100% identity = 1)");
        System.err.println("       orientationPenalty=NUM      if an aligment is in wrong orientation, multiply score by this [0.5]");
        System.err.println("       cutPenalty=NUM              alignment cut penalty [0.001] (0.001 = 1kb gap = -1");
        
        System.err.println("       useChainAndPaf=0            do not use both chain and paf score between contigs when available [not set]");
        
        System.err.println("       proximity=file NUM1 NUM2 NUM3  load proximity data, NUM1=bin size [10000]");
        System.err.println("                                      NUM2=max distance in bins[25], NUM3=scale score [1.0]");


        System.err.println("       minHaplotypeAlignmentScore=NUM   min alignment score required to consider haplotype [-10]");
        System.err.println("       minLinkAlignmentScore=NUM   min alignment score required to consider contig link [-10]\n");

        System.err.println("       evaluateAnchoring=FILE      load initial anchring from a FILE (experimental)");
        System.err.println("       improveAnchoring=1          improve loaded initial anchring (experimental)");

        System.err.println("       cutSites=FILE               list possible contig cut sites for contigs");
        System.err.println("                                   the first cut site within each cut region is taken");
        
        System.err.println("       findContigErrors=1          Iteratively find possible contig errors based on the map only");
        System.err.println("                                   paf, proximity and keepEmptyIntervals not allowed");

        System.err.println("       minImprovement=NUM          minimum improvement for findContigErrors [1]");
        System.err.println("                                   improvement of (NUM + 2*number_of_cuts) required");

        System.err.println("       printAnimation=file	       print iterative solutions of findContigErrors to file0.la,...,fileN.la");

        System.err.println("       alternativeJoins=1	       prints alternative scaffolding joins, does not consider map information (yet)");
        System.err.println("       linksWithin=FILE	           only keep links between contig pairs listed in FILE");

	}
	
	public static void main(String[] args)
	{
		
        if (args.length == 0) {
        	usageInfo();
            System.exit(0);
        }
	    String extraParameters = "";
	    for (int i = 0; i < args.length; ++i) {
	            extraParameters +=  " " + args[i];
	    }
		ParameterParser pp = new ParameterParser();
		if (!pp.init(extraParameters)) {
        	usageInfo();
            System.exit(0);
		}
		pp.warning(new String[]{"findContigErrors","maxPafScore", "map", "bed", "orientation", "chain", "noIntervals", "noChromosome", "randomOrder", "paf", "evaluateAnchoring", "numRuns", "numErrors", "numErrorsPerContig", "scaleScore", "orientationPenalty", "cutPenalty", "maxBridge", "scalePaf", "minLinkAlignmentScore","minHaplotypeAlignmentScore", "maxIntersect", "chromosome", "keepEmptyIntervals", "cutSites", "useChainAndPaf", "proximity", "printAnimation", "compressMap", "minImprovement", "improveAnchoring", "alternativeJoins", "numThreads", "linksWithin"});
		
//		ArrayList<Marker> m = InputData.loadMap(pp.getValueAsString("map", null));
		//PlaceAndOrientContigs.solve(m);
		
		PlaceAndOrientContigs poc = new PlaceAndOrientContigs();
		
		System.out.println("#java PlaceAndOrientContigs" + extraParameters);
		
		boolean findContigErrors = pp.getValueAsString("findContigErrors", "0").equals("1");
		
		String chain = pp.getValueAsString("chain", null);
		String paf = pp.getValueAsString("paf", null);
		String prox = pp.getValueAsString("proximity", 0, null);

		if (findContigErrors && (paf != null || prox != null || pp.getValueAsString("keepEmptyIntervals", "0").equals("1"))) {
			System.err.println("parameters paf, proximity and keepEmptyIntervals not allowed with findContigErrors!");
			System.exit(-1);
		}

		poc.setNumThreads(Integer.parseInt(pp.getValueAsString("numThreads", "1")));
		poc.setScaleScore(Double.parseDouble(pp.getValueAsString("scaleScore", "0.00001")));

		poc.setMaxIntersect(Integer.parseInt(pp.getValueAsString("maxIntersect", "2000")));

		poc.setOrientationPenalty(Double.parseDouble(pp.getValueAsString("orientationPenalty", "0.5")));
		poc.setCutPenalty(Double.parseDouble(pp.getValueAsString("cutPenalty", "0.001")));
		poc.setMaxBridge(Integer.parseInt(pp.getValueAsString("maxBridge", "50000")));

		int numMaps = pp.getNumberOfValues("map"); 

		int chromosome = Integer.parseInt(pp.getValueAsString("chromosome", "-1"));
		String bed = pp.getValueAsString("bed", null);
		if (bed != null && numMaps > 0) 
			poc.loadBed(bed, chromosome);
		else {
			System.err.println("You have to provide bed and one or more map file(s)!");
			System.exit(-1);
		}

		String cutSites = pp.getValueAsString("cutSites", null);
		if (cutSites != null) {
			poc.loadCutSites(cutSites);
		}
		
		boolean findOrientation = false;
		for (int i = 0; i < pp.getNumberOfValues("map"); ++i) {
			String o = pp.getValueAsString("orientation", i, "?");
			if (!o.equals("+") && !o.equals("-"))
				findOrientation = true;
		}
		if (pp.getNumberOfValues("map") != pp.getNumberOfValues("orientation") && pp.getNumberOfValues("orientation") > 0) {
			System.err.println("You have to provide orientation for all maps or none!");
			System.exit(-1);
		}

		poc.setCompressMap(pp.getValueAsString("compressMap", "1").equals("1"));
		
		for (int i = 0; i < pp.getNumberOfValues("map"); ++i) {
			String o = pp.getValueAsString("orientation", i, "?");
			poc.addMap(pp.getValueAsString("map", i, null), !findOrientation && o.equals("-"), pp.getValueAsString("noChromosome", "0").equals("1"), pp.getValueAsString("noIntervals", "0").equals("1"), chromosome);
		}

		if (chain != null) {
			poc.loadChain(chain);
		}
		
		int maxPafScore = Integer.parseInt(pp.getValueAsString("maxPafScore", "" + Integer.MAX_VALUE));
		
		if (paf != null)
			poc.loadPaf(paf, maxPafScore, Double.parseDouble(pp.getValueAsString("scalePaf", "1")));

		if (prox != null) {
			int bin = Integer.parseInt(pp.getValueAsString("proximity", 1, "10000"));
			int maxD = Integer.parseInt(pp.getValueAsString("proximity", 2, "25"));
			double scale = Double.parseDouble(pp.getValueAsString("proximity", 3, "1.0"));
			poc.loadProximity(prox, bin, maxD, scale);
		}
		
		poc.setNumRuns(Integer.parseInt(pp.getValueAsString("numRuns", "5")));

		poc.setNumErrors(Integer.parseInt(pp.getValueAsString("numErrors", "40")));
		poc.setNumErrorsPerContig(Integer.parseInt(pp.getValueAsString("numErrorsPerContig", "3")));

		poc.setMinHaplotypeAlignmentScore(Integer.parseInt(pp.getValueAsString("minHaplotypeAlignmentScore", "-10")));
		poc.setMinLinkAlignmentScore(Integer.parseInt(pp.getValueAsString("minLinkAlignmentScore", "-10")));
		
		poc.setKeepEmptyIntervals(pp.getValueAsString("keepEmptyIntervals", "0").equals("1"));
		poc.setUseChainAndPaf(pp.getValueAsString("useChainAndPaf", "1").equals("1"));
		poc.setPrintAnimation(pp.getValueAsString("printAnimation", ""));
		
		String lwf = pp.getValueAsString("linksWithin", null);
		if (lwf != null)
			poc.linksWithin(lwf);
		
		poc.setCommentOutput(findContigErrors);
		
		String eval = pp.getValueAsString("evaluateAnchoring", null);
		if (eval != null) {
  			 poc.combineMaps(findOrientation, false, false);
			 ArrayList<ContigInterval> eval_result = InputData.loadLa(eval);
			 poc.evaluateScore(eval_result, pp.getValueAsString("improveAnchoring", "0").equals("1"));
		} else {
			poc.combineMaps(findOrientation, pp.getValueAsString("randomOrder", "0").equals("1"), true);
		}

		if (findContigErrors) {
			poc.setCommentOutput(false);
			poc.findContigErrors(Integer.parseInt(pp.getValueAsString("minImprovement", "1")));
		}
		
		if (pp.getValueAsString("alternativeJoins", "0").equals("1")) {
			poc.printAlternativeJoins();
		}
		
		
		//System.out.println(m);
	}
	
	private String[] contigs(String key)
	{
		String split[] = key.split("\t");
		int ip = split[2].indexOf('+');
		int im = split[2].indexOf('-');
		int p = Math.max(ip, im);
		if (ip > 0 && im > 0)
			p = Math.min(ip, im);
		return new String[]{split[0], split[3].substring(p + 1)};
		
	}
	
	private void linksWithin(String lwf) {
		// TODO Clear proximity as well...
		HashMap<String, Integer> pairs = new HashMap<String, Integer>();
		ArrayList<ArrayList<String>> table = Input.loadTable(lwf, "\t ");
		for (ArrayList<String> row : table) {
			if (row.size() >= 2)
				pairs.put(row.get(0) + "\t" + row.get(1), 1);
		}
		
		for (String key : chainLinkHash.keySet()) {
			String cs[] = contigs(key);
			String p1 = cs[0] + "\t" + cs[1];
			if (!pairs.containsKey(p1)) {
				String p2 = cs[1] + "\t" + cs[0];
				if (!pairs.containsKey(p2))
					chainLinkHash.put(key, 0);
			}
		}
		for (String key : scaffoldingLink.keySet()) {
			String cs[] = contigs(key);
			String p1 = cs[0] + "\t" + cs[1];
			if (!pairs.containsKey(p1)) {
				String p2 = cs[1] + "\t" + cs[0];
				if (!pairs.containsKey(p2))
					scaffoldingLink.put(key, 0);
			}
		}
	}
	private int findLastIndex(ArrayList<ContigInterval> cis, int start, int maxDistance, boolean right) {
		if (prox == null)
			return start;
		int n = cis.size();
		int end = start;
		int d = prox.binLength(cis.get(start));
		if (right)
			while (end + 1 < n && d < maxDistance) {
				++end;
				d += prox.binLength(cis.get(end));
			}
		else
			while (end > 0 && d < maxDistance) {
				--end;
				d += prox.binLength(cis.get(end));
			}
		
		return end;
	}
	
	//score between c and c + 1, c in 0,1,...,cis.size() - 2
	private int linkScore(ArrayList<ContigInterval> cis, int c) {
		int ret = calculateChainScore(cis.get(c), cis.get(c + 1));
		if (prox != null)
			ret += prox.linkScore(cis, c, c + 1);
		return ret;
	}
	
	//score of [start,end] => [start2,end2] in all orientations (++, -+, +-, --)
	private int[] scores(ArrayList<ContigInterval> cis, int start, int end, int start2, int end2){
		ArrayList<ContigInterval> tmp1 = new ArrayList<ContigInterval>();
		for (int i = start; i <= end; ++i)
			tmp1.add(cis.get(i));

		ArrayList<ContigInterval> tmp2 = new ArrayList<ContigInterval>();
		for (int i = start2; i <= end2; ++i)
			tmp2.add(cis.get(i));
		
		int s[] = new int[4];
		for (int o2 = 0; o2 < 2; ++o2) { 
			for (int o1 = 0; o1 < 2; ++o1) { 
				ArrayList<ContigInterval> tmp = new ArrayList<ContigInterval>();
				tmp.addAll(tmp1);
				tmp.addAll(tmp2);
				s[o1 + 2 * o2] = linkScore(tmp, end - start);

				Collections.reverse(tmp1);
				for (ContigInterval c : tmp1)
					c.flipOrientation();
			}
			Collections.reverse(tmp2);
			for (ContigInterval c : tmp2)
				c.flipOrientation();
		}
		return s;
	}
	
	private void printAlternativeJoins() {
		int maxD = 0;
		if (prox != null)
				maxD = prox.getMaxDistance();
		int n = intervalsInChr.size();
		if (n <= 1)
			return;
		
		int jTable[][] = new int[2 * n][2 * n];
		int start = 0;
		while (start < n) {
			int end_ = findLastIndex(intervalsInChr, start, maxD, true); // this could be faster...
			for (int end = ((start == 0) ? 0 : end_); end <= end_; ++end)
				for (int start2 = end + 1; start2 < n; ++start2) {
					int end2 = findLastIndex(intervalsInChr, start2, maxD, true); // and this
					int s[] = scores(intervalsInChr, start, end, start2, end2);

					//end2 is always maximal
					boolean endMaximal = (end == end_);
					boolean startMaximal = (findLastIndex(intervalsInChr, end, maxD, false) == start);
					boolean start2Maximal = (findLastIndex(intervalsInChr, end2, maxD, false) == start2);
					//remove non-maximal scores, otherwise you count each score multiple times for same endpoint with other end not maximal...

					if (startMaximal)
						jTable[2 * end]      [2 * start2]   += s[0]; // += is not needed, = is enough
					if (endMaximal)
						jTable[2 * start + 1][2 * start2]   += s[1];
					if (startMaximal && start2Maximal)
						jTable[2 * end]      [2 * end2 + 1] += s[2];
					if (endMaximal && start2Maximal)
						jTable[2 * start + 1][2 * end2 + 1] += s[3];
							
					int s2[] = scores(intervalsInChr, start2, end2, start, end);

					if (start2Maximal && endMaximal)
						jTable[2 * end2]      [2 * start]   += s2[0];
					if (endMaximal)
						jTable[2 * start2 + 1][2 * start]   += s2[1];
					if (start2Maximal && startMaximal)
						jTable[2 * end2]      [2 * end + 1] += s2[2];
					if (startMaximal)
						jTable[2 * start2 + 1][2 * end + 1] += s2[3];
				}
			++start;
		}
/*		for (int i = 0; i < n + n; ++i) {
			StringBuilder sb = new StringBuilder();
			for (int j = 0; j < n + n; ++j) {
				if (j != 0) 
					sb.append('\t');
				if (i == j)
					sb.append('X');
				else
				sb.append(jTable[i][j]);
			}
			System.err.println(sb);
		}*/
		System.err.println("*** alternativeJoins ***");
		for (int ci = 0; ci < n; ++ci) {
			//ArrayList<Integer> list = new ArrayList<Integer>();
			for (int i = 0; i < n + n; ++i) {
				int ip2 = i >> 1;
				if (jTable[2 * ci][i] > 0) {
					System.err.println(intervalsInChr.get(ci) + "\t+\t" + intervalsInChr.get(ip2) + "\t" + ("+-".substring(i & 1, (i & 1) + 1)) + "\t" + jTable[2 * ci][i] + "\t" + ci + "\t" + ip2);
				}
				if (jTable[2 * ci + 1][i] > 0) {
					System.err.println(intervalsInChr.get(ci) + "\t-\t" + intervalsInChr.get(ip2) + "\t" + ("+-".substring(i & 1, (i & 1) + 1)) + "\t" + jTable[2 * ci + 1][i] + "\t" + ci + "\t" + ip2);
				}
			}
			
		}
	}
}
