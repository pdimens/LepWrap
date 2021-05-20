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
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;



public class InputData {
	
	public static long myParseLong(String s)
	{
		if (s.charAt(s.length() - 1) == '*')
			return Long.parseLong(s.substring(0, s.length() - 1));
		else
			return Long.parseLong(s);
			 
	}
	private static long[] myParseLongInterval(String s)
	{
		long ret[] = null;
		int pos = s.indexOf('-');
		if (pos < 0) {
			long p = myParseLong(s);
			return new long[]{p, p};		
		}
		if (s.charAt(s.length() - 1) == '*') {
			ret = new long[3];
			ret[1] =  Long.parseLong(s.substring(pos + 1, s.length() - 1));			
		}
		else {
			ret = new long[2];
			ret[1] =  Long.parseLong(s.substring(pos + 1));
		}
		ret[0] =  Long.parseLong(s.substring(0, pos));
		return ret;
	}
	
	public static HashMap<ContigInterval, ArrayList<ArrayList<String>>> loadRaw(String fn, HashMap<String, ArrayList<ContigInterval>> haplotypeHash)
	{
		HashMap<ContigInterval, ArrayList<ArrayList<String>>> ret = new HashMap<ContigInterval, ArrayList<ArrayList<String>>>();
		try {
			
			int numMarkers = 0;

			BufferedReader br = null;
			if (fn.equals("-"))
				br = new BufferedReader(new InputStreamReader(System.in));
			else
				br = new BufferedReader(new FileReader(fn));
			
			do {
				String row = Input.loadRow(br);
				if (row == null)
					break;
				ArrayList<String> srow = Input.splitRow(row, "[\t ]");
				if (srow.size() >= 3) {
					String contig = srow.get(0);
					boolean added = false;
					if (haplotypeHash != null && haplotypeHash.containsKey(contig)) {  
						long pos = myParseLong(srow.get(1));
						for (ContigInterval ci : haplotypeHash.get(contig)) {
							if (ci.inside(pos)) {
								ArrayList<ArrayList<String>> list = ret.get(ci);
								if (list == null) {
									list = new ArrayList<ArrayList<String>>();
									ret.put(ci, list);
								}
								list.add(srow);
								++numMarkers;
								added = true;
								break;
							}
						}
					}
					if (!added)
						System.out.println(row);
				} else
					System.err.println("Warning: skipping " + row);
			} while (true);
			
			System.err.println("Liftover for " + numMarkers + " markers over " + ret.size() + " regions");
			return ret;
		} catch (Exception e) {
	     	e.printStackTrace();
			System.err.println("Error in file " + fn);
			return null;
      	}
	}
	
	public static ArrayList<Marker> loadMap(String fn, boolean nochromosome, boolean nointervals, boolean compressMap)
	{
		ArrayList<Double> pos = new ArrayList<Double>(); 
		ArrayList<Integer> chr = new ArrayList<Integer>(); 

		ArrayList<Marker> ret = new ArrayList<Marker>(); 
		try {
			BufferedReader br = null;
			if (fn.equals("-"))
				br = new BufferedReader(new InputStreamReader(System.in));
			else
				br = new BufferedReader(new FileReader(fn));
			
			do {
				ArrayList<String> row = Input.loadTableRow(br, "[\t ]");
				if (row == null)
					break;
				if (row.size() > 0) {
					if (nochromosome)
						row.add(2, "1");
					
					if (row.size() >= 4 && (nointervals || (row.size() & 1) == 1) || row.size() == 4) {
						int intervals[] = new int[Math.max(2, row.size() - 3)];
						if (row.size() == 4 || nointervals) {
							if (nointervals) {
								intervals = new int[2];
								double p = Double.parseDouble(row.get(3));
								if (row.size() >= 5)
									p += Double.parseDouble(row.get(4));
								pos.add(p);
								chr.add(Integer.parseInt(row.get(2)));
							}
							else
								try {
									intervals[0] = intervals[1] =  Integer.parseInt(row.get(3));
								} catch (NumberFormatException e) {
									System.err.println("Error: non-integer value for interval (check parameter noIntervals=1)" + row);
									System.exit(-1);
								}
						}
						else
							for (int i = 3; i < row.size(); i+=2) {
								try {
									intervals[i-3] = Integer.parseInt(row.get(i));
									intervals[i-2] = Integer.parseInt(row.get(i + 1));
								} catch (NumberFormatException e) {
									System.err.println("Error: non-integer value for interval (check parameter noIntervals=1) " + row);								
									System.exit(-1);
								}
							}
						ret.add(new Marker(row.get(0), myParseLong(row.get(1)), Integer.parseInt(row.get(2)), intervals));
					} else
						System.err.println("Warning: skipping " + row);
				}
			} while (true);
			br.close();
			if (nointervals) {
				//put markers in each chr to separate list... and sort...
				HashMap<Integer, ArrayList<Double>> chrMarkers = new HashMap<Integer, ArrayList<Double>>(); // store chromosome info
				HashMap<Integer, ArrayList<Integer>> indexMarkers = new HashMap<Integer, ArrayList<Integer>>(); // store original index in ret...

				for (int mi = 0; mi < chr.size(); ++mi) {
					int c = chr.get(mi);
					ArrayList<Double> list = chrMarkers.get(c); 
					ArrayList<Integer> list2 = indexMarkers.get(c); 
					if (list == null) {
						list = new ArrayList<Double>();
						chrMarkers.put(c, list);
						list2 = new ArrayList<Integer>();
						indexMarkers.put(c, list2);
					}
					list.add(pos.get(mi));
					list2.add(mi);
				}
				for (int c : chrMarkers.keySet()) {
					ArrayList<Double> list = chrMarkers.get(c); 
					ArrayList<Integer> list2 = indexMarkers.get(c); 
				    Misc.ArrayIndexComparator<Double> comparator = new Misc.ArrayIndexComparator<Double>(list);
				    Integer[] indexes = comparator.createIndexArray();
				    java.util.Arrays.sort(indexes, comparator);
				    
				    double prev = Double.NEGATIVE_INFINITY;
				    int ip = -1;
				    for (int i : indexes) {
				    	double p = list.get(i); 
				    	if (p > prev)
				    		++ip;
				    	Marker m = ret.get(list2.get(i));
				    	m.intervals[0] = ip;
				    	m.intervals[1] = ip;
				    	prev = p;
				    }
					
				}
			} else if (compressMap){ // no noIntervals && compress map positions

				// compressing map speeds up computation... 
				HashMap<Integer, ArrayList<Marker>> chrMarkers = new HashMap<Integer, ArrayList<Marker>>();

				for (Marker m : ret) {
					int c = m.getChromosome();
					ArrayList<Marker> list = chrMarkers.get(c); 
					if (list == null) {
						list = new ArrayList<Marker>();
						chrMarkers.put(c, list);
					}
					list.add(m);
				}
				for (ArrayList<Marker> markers: chrMarkers.values()) {
					ArrayList<Integer> positions = new ArrayList<Integer>(); 
					for (Marker m : markers) {
						int[] interval = m.getIntervals();
						for (int i : interval)  
							positions.add(i);
					}
					Collections.sort(positions);
					HashMap<Integer, Integer> compressHash = new HashMap<Integer, Integer>();
					int numPositions = 0;
					for (int i : positions) {
						if (!compressHash.containsKey(i))
							compressHash.put(i, numPositions++);
					}
					for (Marker m : markers) {
						int[] interval = m.getIntervals();
						for (int ii = 0; ii < interval.length; ++ii)  
							interval[ii] = compressHash.get(interval[ii]); 
					}
				}
			}
			return ret;
	    } catch (Exception e) {
	     	e.printStackTrace();
			System.err.println("Error in file " + fn);
			return null;
      	}
  	}
	public static ArrayList<ContigInterval> loadBedAndComments(String fn, StringBuilder comments)
	{
		Input.setKeepComments(true);
		ArrayList<ContigInterval> ret = loadBed(fn, -1);
		Input.setKeepComments(false);
		comments.append(Input.getComments());
		return ret;
	}
	
	public static ArrayList<ContigInterval> loadBed(String fn)
	{
		return loadBed(fn, -1);
	}
	public static ArrayList<ContigInterval> loadBed(String fn, int chr)
	{
		ArrayList<ContigInterval> ret = new ArrayList<ContigInterval>(); 
		try {
			BufferedReader br = null;
			if (fn.equals("-"))
				br = new BufferedReader(new InputStreamReader(System.in));
			else
				br = new BufferedReader(new FileReader(fn));
			
			do {
				ArrayList<String> row = Input.loadTableRow(br, "[\t ]");
				if (row == null)
					break;
					if (row.size() >= 3) {
						String start = row.get(1);
						String end = row.get(2);
						ContigInterval ci = null;
						if (start.indexOf('-') >= 0 || end.indexOf('-') >= 0)
							ci = new ContigInterval(row.get(0), myParseLongInterval(start), myParseLongInterval(end));
						else
							ci = new ContigInterval(row.get(0), myParseLong(start), myParseLong(end));
						if (row.size() >= 4 && (row.get(3).equals("-")  || row.get(3).equals("--")  || row.get(3).equals("---")))
							ci.flipOrientation();
						
						if (row.size() >= 5) {
							int bChr = Integer.parseInt(row.get(4));
							ci.setChromosome(bChr);
						}
						if (chr < 0 || ci.getChromosome() == chr)
							ret.add(ci);
					}
					else
						System.err.println("Warning: skipping " + row);
			} while (true);
			
			int prevC = 0;
			for (ContigInterval ci : ret) {
				if (prevC != 0 && ci.getChromosome() != prevC) {
					System.err.println("Error: multiple cromosomes in the bed");
					System.err.println("Try providing parameter chromosome or split bed into chromosomes");
					System.exit(-1);
				}
			}
			br.close();
			return ret;
	    } catch (Exception e) {
	     	e.printStackTrace();
			System.err.println("Error in file " + fn);
			return null;
      	}
	}
	
	public static ArrayList<ContigInterval> loadLa(String fn)
	{
		ArrayList<ContigInterval> ret = new ArrayList<ContigInterval>(); 
		try {
			BufferedReader br = null;
			if (fn.equals("-"))
				br = new BufferedReader(new InputStreamReader(System.in));
			else
				br = new BufferedReader(new FileReader(fn));
			
			do {
				ArrayList<String> row = Input.loadTableRow(br, "[\t ]");
				if (row == null)
					break;
					if (row.size() >= 3) {
						String start = row.get(1);
						String end = row.get(2);
						if (row.size() >= 8) {
							start = row.get(6);
							int sis = start.indexOf('/');
							if (sis >= 0)
								start = start.substring(0, sis);
							end = row.get(7);
							int sie = end.indexOf('/');
							if (sie >= 0)
								end = end.substring(0, sie);
						}
						ContigInterval ci = null;
						if (start.indexOf('-') >= 0 || end.indexOf('-') >= 0)
							ci = new ContigInterval(row.get(0), myParseLongInterval(start), myParseLongInterval(end));
						else
							ci = new ContigInterval(row.get(0), myParseLong(start), myParseLong(end));
						if (row.size() >= 4 && (row.get(3).equals("-")  || row.get(3).equals("--")  || row.get(3).equals("---")))
							ci.flipOrientation();
						
						if (row.size() >= 5) {
							int bChr = Integer.parseInt(row.get(4));
							ci.setChromosome(bChr);
						}
						ret.add(ci);
					}
					else
						System.err.println("Warning: skipping " + row);
			} while (true);
			
			int prevC = 0;
			for (ContigInterval ci : ret) {
				if (prevC != 0 && ci.getChromosome() != prevC) {
					System.err.println("Error: multiple cromosomes in the bed");
					System.err.println("Try providing parameter chromosome or split bed into chromosomes");
					System.exit(-1);
				}
			}
			br.close();
			return ret;
	    } catch (Exception e) {
	     	e.printStackTrace();
			System.err.println("Error in file " + fn);
			return null;
      	}
	}

	
	public static ArrayList<ContigInterval> loadHaplotypes(String fn)
	{
		ArrayList<ContigInterval> ret = new ArrayList<ContigInterval>(); 
		try {
			BufferedReader br = null;
			if (fn.equals("-"))
				br = new BufferedReader(new InputStreamReader(System.in));
			else
				br = new BufferedReader(new FileReader(fn));
			
			do {
				ArrayList<String> row = Input.loadTableRow(br, "[\t ]");
				if (row == null)
					break;
					if (row.size() >= 4) {
						ContigInterval ci = new ContigInterval(row.get(1), myParseLong(row.get(2)), myParseLong(row.get(3)));
						ret.add(ci);
					}
					else
						System.err.println("Warning: skipping haplotype " + row);
			} while (true);
			br.close();
			return ret;
	    } catch (Exception e) {
	     	e.printStackTrace();
			System.err.println("Error in file " + fn);
			return null;
      	}
	}
	

	public ArrayList<Marker> loadLGMap(String fn)
	{
		ArrayList<Marker> ret = new ArrayList<Marker>(); 
		try {
			BufferedReader br = null;
			if (fn.equals("-"))
				br = new BufferedReader(new InputStreamReader(System.in));
			else
				br = new BufferedReader(new FileReader(fn));
			String line = null;
			do {
				line = br.readLine();
				if (line == null)
					break;
		
				int ic = line.indexOf('#');
				if (ic >= 0)
					line = line.substring(0, ic); 
				
				if (line.length() > 0) {
					String[] row = line.split("\t");
					if (row.length > 2 || row.length == 3 || (row.length & 1) == 0) {
						int intervals[] = new int[Math.max(2, row.length - 2)];
						if (row.length == 3)
							intervals[0] = intervals[1] =  Integer.parseInt(row[2]);
						else
							for (int i = 2; i < row.length; i+=2) {
									intervals[i-2] = Integer.parseInt(row[i]);
									intervals[i-1] = Integer.parseInt(row[i + 1]);
							}
						//ret.add(new Marker(row[0], myParseInt(row[1]), intervals));
					} else
						System.err.println("Warning: skipping " + line);
				}
			} while (true);
			br.close();
			return ret;
	    } catch (Exception e) {
	     	e.printStackTrace();
			System.err.println("Error in file " + fn);
			return null;
      	}
  	}
	

}
