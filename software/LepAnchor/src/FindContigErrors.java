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

import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.*;

//Calculate score over positions of each contig (paf + ?)
//...and to different position...

public class FindContigErrors {
	//calculate coverage for read mapping intervals...
	//could be used to find (about) exact cut positions
	
	//prox, use borderScore
	//map,

	private void addPositions(String contig, HashMap<String, ArrayList<long[]>> positionHash, long start, long end, long weight)
	{
		ArrayList<long[]> list = positionHash.get(contig);
		if (list == null) {
			list = new ArrayList<long[]>();
			positionHash.put(contig, list);
		}
		list.add(new long[]{start, end, weight});
	}
	
	ArrayList<String> contigs = new ArrayList<String>();
	ArrayList<Long> contigLengths = new ArrayList<Long>();
	
	HashMap<String, Integer> contigHash = new HashMap<String, Integer>();
	
	private void addContig(String name, long length){
		if (!contigHash.containsKey(name)) {
			contigHash.put(name, contigs.size());
			contigs.add(name);
			contigLengths.add(length);
		}
	}

	//TODO: is this needed? Could we use cov instead?
	private void mergeCovs(String contig, ArrayList<Long> cov, ArrayList<Long> cov2) {
		System.out.println(contig);
		int covi = 0;
		int cov2i = 0;
		//merge in cov and cov2

		long pos = 0;
		long cp = 0;
		long cp2 = 0;
		while (covi < cov.size() || cov2i < cov2.size()) {
			if (covi >= cov.size() ) {
				pos = cov2.get(cov2i);
				cp2 = cov2.get(cov2i + 1);
				cov2i+=2;
			}
			else if (cov2i >= cov2.size() ) {
				pos = cov.get(covi);
				cp = cov.get(covi + 1);
				covi+=2;
			}
			else {
				if (cov.get(covi) <= cov2.get(cov2i)) {
						pos = cov.get(covi);
						cp = cov.get(covi + 1);
						covi+=2;
				}
				else {
					pos = cov2.get(cov2i);
					cp2 = cov2.get(cov2i + 1);
					cov2i+=2;
				}
			}
			System.out.println(pos + "\t" + cp + "\t" + cp2);
		}
	}
	
	HashMap<String, ArrayList<long[]>> pafPositions = new HashMap<String, ArrayList<long[]>>();
	HashMap<String, ArrayList<long[]>> pafPositions2 = new HashMap<String, ArrayList<long[]>>();
	public void processErrorsFromPaf(String fn)
	{
 		System.err.println("loading paf...");
 		ArrayList<ArrayList<String>> rows = new ArrayList<ArrayList<String>>(); 
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
				if (rows.size() == 0 || rows.get(0).get(0).equals(row.get(0)))
					rows.add(row);
				else {
					processPAF(rows);
					rows.clear();
					rows.add(row);
				}

			} while (true);
			processPAF(rows);

			for (String contig : contigs) {
				ArrayList<long[]> list = pafPositions.get(contig);
				ArrayList<Long> cov = null;
				if (list != null)
					cov = Misc.cov(list);
				else {
					cov = new ArrayList<Long>();
					cov.add(0L);
					cov.add(0L);
				}

				ArrayList<long[]> list2 = pafPositions2.get(contig);
				ArrayList<Long> cov2 = null;
				if (list2 != null)
					cov2 = Misc.cov(list2);
				else {
					cov2 = new ArrayList<Long>();
					cov2.add(0L);
					cov2.add(0L);
				}

				mergeCovs(contig, cov, cov2);
				//for (int i = 0; i < cov.size(); i+=2) {
				//	System.out.println(cov.get(i) + "\t" + cov.get(i + 1));
				//}
			}

		}
		catch (Exception e) {
			e.printStackTrace();
			System.err.println("Error in file " + fn);
  		}
	}
	
	private void processPAF(ArrayList<ArrayList<String>> rows) {
		for (ArrayList<String> row : rows) {
			String contig = row.get(5); 
			long length = Long.parseLong(row.get(6));
			addContig(contig, length);
			long start = Long.parseLong(row.get(7)) + 1;
			long end = Long.parseLong(row.get(8));
			addPositions(contig, pafPositions, start, end, 1); // calculate coverage...
		}
		//take two longest alignments to two different contigs...  
		long maxL = 0;
		ArrayList<String> maxRow = null;
		for (ArrayList<String> row : rows) {
			long aLen = Long.parseLong(row.get(8)) - Long.parseLong(row.get(7));
			if (aLen > maxL) {
				maxL = aLen;
				maxRow = row;
			}
		}
		long maxL2 = 0;
		ArrayList<String> maxRow2 = null;
		for (ArrayList<String> row : rows)
			if (!row.get(5).equals(maxRow.get(5))) {
				long aLen = Long.parseLong(row.get(8)) - Long.parseLong(row.get(7));
				if (aLen > maxL2) {
					maxL2 = aLen;
					maxRow2 = row;
				}
			}
		if (maxRow2 != null) { // two aligments found to two contigs
			long rstart1 =  Long.parseLong(maxRow.get(2));
			long rstart2 =  Long.parseLong(maxRow2.get(2));
			if (rstart1 > rstart2) { // make sure maxRow starts first...
				long tmp = rstart2;
				rstart2 = rstart1;
				rstart1 = tmp;
				
				tmp = maxL;
				maxL = maxL2;
				maxL2 = tmp;
				
				ArrayList<String> tmp2 = maxRow;
				maxRow = maxRow2;
				maxRow2 = tmp2;
			}
			long rend1 =  Long.parseLong(maxRow.get(3));
			//long rend2 =  Long.parseLong(maxRow2.get(3));
			
			long distance = rstart2 - rend1 + 2; //
			
			long minL = Math.min(maxL, maxL2);
			
			if (distance < -minL || distance > 2 * minL) // no overlap distance 2 x the length of the shorter alignment
				return;
			if (distance < 0)
				distance = 0;
//			TODO: figure out what happens if distance < 0, maybe distance can be -2 minL
			
			String o1 = maxRow.get(4);
			if ("+".equals(o1)) {
				long end1 = Long.parseLong(maxRow.get(8));
				long length1 = Long.parseLong(maxRow.get(6));
				addPositions(maxRow.get(5), pafPositions2, end1 + 1, Math.min(length1 + 1, end1 + 1 + distance), 1); // calculate coverage...
			} else {
				long start1 = Long.parseLong(maxRow.get(7));
				addPositions(maxRow.get(5), pafPositions2, Math.max(0, start1 - distance), start1, 1); // calculate coverage...
			}
			String o2 = maxRow2.get(4);
			if ("+".equals(o2)) {
				long start2 = Long.parseLong(maxRow2.get(7));
				addPositions(maxRow2.get(5), pafPositions2, Math.max(0, start2 - distance), start2, 1); // calculate coverage...
			} else {
				long end2 = Long.parseLong(maxRow2.get(8));
				long length2 = Long.parseLong(maxRow2.get(6));
				addPositions(maxRow2.get(5), pafPositions2, end2 + 1, Math.min(length2 + 1, end2 + 1 + distance), 1); // calculate coverage...
			}
			
		}
	}

	HashMap<String, ArrayList<long[]>> chainPositions = new HashMap<String, ArrayList<long[]>>(); 
	public void processErrorsFromChain(String fn)
	{
 		System.err.println("loading chain...");
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
				if (row.get(4).equals("+")) {
					String contig1 = row.get(2);
					long start1 = Long.parseLong(row.get(5)) + 1;
					long end1 = Long.parseLong(row.get(6));
					long aScore = Long.parseLong(row.get(1));

					long length1 = Long.parseLong(row.get(3));

					addContig(contig1, length1);
					addPositions(contig1, chainPositions, start1, end1, aScore / length1);

					String contig2 = row.get(7);
					long length2 = Long.parseLong(row.get(8));
					addContig(contig2, length2);
					long start2 = Long.parseLong(row.get(10));
					long end2 = Long.parseLong(row.get(11));
					if (row.get(9).equals("+"))
						addPositions(contig2, chainPositions, start2 + 1, end2, aScore / length2);
					else {
						addPositions(contig2, chainPositions, length2 - end2 + 1, length2 - start2, aScore / length2);
					}
				}
			} while (true);

		}
		catch (Exception e) {
			e.printStackTrace();
			System.err.println("Error in file " + fn);
  		}
	}

//	1 (n starts)
//	2,3 (n/2)
//	4,5,6,7 (n/4)
//	8,9,10,11,12,13,14,15 (n/8)
// 16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31 (n/16)
// ...

	private int getPosition(Marker m, boolean plusOrientation) {
		if (plusOrientation)
			return m.pPlus;
		else
			return m.pMinus;
	}
	
	HashMap<String, ArrayList<long[]>> proximityPositions = new HashMap<String, ArrayList<long[]>>(); 
	HashMap<String, ArrayList<long[]>> proximityPositions2 = new HashMap<String, ArrayList<long[]>>(); 
	public void processErrorsFromProximity(String fn, int bin, int maxD, double scale)
	{
 		System.err.println("loading proximity data...");
 		long maxDistance = bin * maxD;
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
				if (row.size() < 5) {
					System.err.println("Warning: skipping " + row);
					continue;
				}
				
				String contig1 = row.get(0);
				String contig2 = row.get(2);
				long pos1 = Long.parseLong(row.get(1));
				long pos2 = Long.parseLong(row.get(3));
				long weight = (long) Double.parseDouble(row.get(4));
				if (contig1.equals(contig2)) {
					if (pos1 < pos2 && pos2 - pos1 < maxDistance) {
						addPositions(contig1, proximityPositions, pos1, pos2, weight);
					}
				} else {
					addPositions(contig1, proximityPositions2, pos1, pos1 + bin, weight);
					addPositions(contig2, proximityPositions2, pos2, pos2 + bin, weight);
				}
			} while (true);

			for (String key : proximityPositions.keySet())
				if (!contigHash.containsKey(key))
					addContig(key, 0);

			for (String key : proximityPositions2.keySet())
				if (!contigHash.containsKey(key))
					addContig(key, 0);
				
			
			for (String contig : contigs) {
				ArrayList<long[]> list = proximityPositions.get(contig);
				ArrayList<Long> cov = null;
				if (list != null)
					cov = Misc.cov(list);
				else {
					cov = new ArrayList<Long>();
					cov.add(0L);
					cov.add(0L);
				}

				ArrayList<long[]> list2 = proximityPositions2.get(contig);
				ArrayList<Long> cov2 = null;
				if (list2 != null)
					cov2 = Misc.cov(list2);
				else {
					cov2 = new ArrayList<Long>();
					cov2.add(0L);
					cov2.add(0L);
				}

				mergeCovs(contig, cov, cov2);
				//for (int i = 0; i < cov.size(); i+=2) {
				//	System.out.println(cov.get(i) + "\t" + cov.get(i + 1));
				//}
			}

		}
		catch (Exception e) {
			e.printStackTrace();
			System.err.println("Error in file " + fn);
  		}
		
		
	}
		
	public void findErrors() {
		int numContigs = contigs.size();
		for (int ci = 0; ci < numContigs; ++ci) {
			String contig = contigs.get(ci);
			long length = contigLengths.get(ci);
			ArrayList<long[]> list = pafPositions.get(contig);
		}
		
		/*ArrayList<Long> t = new ArrayList<Long>();
		t.add(-2l);
		t.add(-2l);
		t.add(-1l);
		t.add(-1l);
		t.add(-1l);
		t.add(10l);
		t.add(10l);
		t.add(10l);
		t.add(11l);
		t.add(12l);
		cov(t);*/
		//paf
		//for (String contig: contigs) {
		//	ArrayList<Long> list = pafCoverage.get(contig);
		//	if (list != null) {
		//		Collections.sort(list);
				//System.err.println(list);
				//cov(list);
				//ArrayList<Long> new_list = new ArrayList<Long>();
		//	}
		//}
		//chain
		
		//prox
		
		//map
		
	}
	/*
	if (processErrors) {
		for (int c1 = 0; c1 < numContigs; ++c1) {
			String contig1 = contigs.get(c1); 
			ArrayList<Long> list = pafCoverage.get(contig1);
			if (list == null) {
				list = new ArrayList<Long>();
				pafCoverage.put(contig1, list);
			}
			ArrayList<ArrayList<String>> rows1 = alignmentHash.get(contig1);
			
			for (ArrayList<String> row1 : rows1) {
				long start1 = Long.parseLong(row1.get(7)) + 1;
				long end1 = Long.parseLong(row1.get(8));
				list.add(-start1); // add -start
				list.add(end1);    // and end position
			}

			//TODO: Try joining adjacent aligment blocks if the gap <= min{maxBridge, length1, length2} 
			/*ArrayList<Long> comp = new ArrayList<Long>(); // sort row by start position... 
			/for (ArrayList<String> row1 : rows1)
				comp.add(Long.parseLong(row1.get(7))); // start
		    Misc.ArrayIndexComparator<Long> comparator = new Misc.ArrayIndexComparator<Long>(comp);
		    Integer[] indexes = comparator.createIndexArray();
		    Arrays.sort(indexes, comparator);
		    
			ArrayList<Long> list = pafCoverage.get(contig1);
			if (list == null) {
				list = new ArrayList<Long>();
				pafCoverage.put(contig1, list);
			}
		    
			for (int orientation = 0; orientation < 2; ++orientation) {
		    	long minStart = Long.MAX_VALUE;
		    	long maxEnd = 0;

		    	long minRStart = Long.MAX_VALUE;
		    	long maxREnd = 0;
			    for (int i : indexes) {
			    	ArrayList<String> row1 = rows1.get(i);
			    	
			    	if ((orientation != 0) ^ ("+".equals(row1.get(4)))) { // take only + or - orientation...
				    	long start1 = Long.parseLong(row1.get(7)) + 1;
						long end1 = Long.parseLong(row1.get(8));
						
						long oldStart = minStart;
						long oldEnd = maxEnd;
						
						minStart = Math.min(minStart, start1);
						maxEnd = Math.max(maxEnd, end1);

						long rstart1 = Long.parseLong(row1.get(2));
						long rend1 = Long.parseLong(row1.get(3));

						minRStart = Math.min(minRStart, rstart1);
						maxREnd = Math.max(maxREnd, rend1);
						
						long RL = maxREnd - minRStart; // read length
						long L = maxEnd - minStart;    // contig length
						
						if (oldEnd == 0 || RL <= 1.5 * L && L <= 1.5 * RL) {
							 
						} else {
							list.add(-oldStart); // add -start
							list.add(oldEnd);    // and end position
						}
						
			    	}
				}
			}*/

	
	private static void usageInfo()
	{
        System.err.println("usage: java FindContigErrors [options]");        
        System.err.println("       chain=file                  chain file ");
        
        System.err.println("       paf=file                    load alignment file in paf (minimap2) format");

        System.err.println("       proximity=file NUM1 NUM2 NUM3  load proximity data, NUM1=bin size [10000]");
        System.err.println("                                      NUM2=max distance in bins[25], NUM3=scale score [1.0]");
	}
	
	private static void test()
	{
		ArrayList<long[]> list = new ArrayList<long[]>();
		list.add(new long[]{0, 9, 1});
		list.add(new long[]{0, 10, 1});
		list.add(new long[]{0, 11, 1});
		list.add(new long[]{5, 10, 1});
		FindContigErrors fce = new FindContigErrors(); 
		System.err.println(Misc.cov(list));
	}
	
	public static void main(String[] args)
	{
		//FindContigErrors.test();
		
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
		pp.warning(new String[]{"paf", "chain", "proximity"});
				
		FindContigErrors fe = new FindContigErrors();
		System.out.println("#java FindContigErrors" + extraParameters);
		
		//for (int i = 0; i < pp.getNumberOfValues("map"); ++i) {
		//	fe.processErrorsFromMap(pp.getValueAsString("map", i, null), pp.getValueAsString("noIntervals", "0").equals("1"));
		//}

		String chain = pp.getValueAsString("chain", null);
		if (chain != null) {
			fe.processErrorsFromChain(chain);
		}

		String paf = pp.getValueAsString("paf", null);
		if (paf != null)
			fe.processErrorsFromPaf(paf);

		String prox = pp.getValueAsString("proximity", 0, null);
		if (prox != null) {
			int bin = Integer.parseInt(pp.getValueAsString("proximity", 1, "10000"));
			int maxD = Integer.parseInt(pp.getValueAsString("proximity", 2, "25"));
			double scale = Double.parseDouble(pp.getValueAsString("proximity", 3, "1.0"));
			fe.processErrorsFromProximity(prox, bin, maxD, scale);
		}
		fe.findErrors();
	}

}
