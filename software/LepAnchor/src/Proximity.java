import java.io.BufferedReader;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;


public class Proximity {
	private int binSize;
	private int maxDistance;
	private double scale;
	//eg. binSize=10000 and maxDistance=25 =>10kb bins and max of 250kb
	
	//TODO: find possible assembly errors...

	public Proximity(int binSize, int maxDistance, double scale){
		this.binSize = binSize;
		this.maxDistance = maxDistance;
		this.scale = scale;
	}
	HashMap<String, int[]> distanceHash = new HashMap<String, int[]>();
	
	public int getMaxDistance(){
		return maxDistance;
	}
	
	public boolean loadData(String fn, HashMap<String, ArrayList<ContigInterval>> bedHash)
	{
		boolean nonEmpty = false;
		System.err.println("loading proximity data...");
		HashMap<String, String> symmetryHash = new HashMap<String, String>();
		try {
			HashMap<String, double[]> distanceHash_tmp = new HashMap<String, double[]>();
			
			BufferedReader br = null;
			if (fn.equals("-"))
				br = new BufferedReader(new InputStreamReader(System.in));
			else
				br = new BufferedReader(new FileReader(fn));
			do {
				ArrayList<String> row = Input.loadTableRow(br, "\t ");
				if (row == null)
					break;
				String contig1 = row.get(0);
				long pos1 = Long.parseLong(row.get(1));
			if (bedHash.containsKey(contig1))
				for (ContigInterval ci1 : bedHash.get(contig1))
					if (ci1.inside(pos1)) {
						
						int p1 = bin(pos1);
						
						int de1 = bin(ci1.getEnd()) - p1;
						int ds1 = p1 - bin(ci1.getStart());
						
						if (de1 >= maxDistance && ds1 >= maxDistance)
							continue;

						String contig2 = row.get(2);
						long pos2 = Long.parseLong(row.get(3));
					if (bedHash.containsKey(contig2))
						for (ContigInterval ci2 : bedHash.get(contig2))
							if (ci1 != ci2 && ci2.inside(pos2)) {
								//DO MAGIC HERE
								
								int p2 = bin(pos2);
								int de2 = bin(ci2.getEnd()) - p2;
								int ds2 = p2 - bin(ci2.getStart()); 
							
								if (de2 >= maxDistance && ds2 >= maxDistance || Math.min(de1, ds1) + Math.min(de2, ds2) >= maxDistance)
									continue;
								
								nonEmpty = true;
								
								double value = Double.parseDouble(row.get(4)) * scale;
								
								if (de1 + ds2 < maxDistance) { // ++
									String key = calculateKey(ci1, true, ci2, true);
									symmetryHash.put(key, calculateKey(ci2, false, ci1, false));
									
									double list[] = distanceHash_tmp.get(key);
									if (list == null) {
										list = new double[maxDistance];
										distanceHash_tmp.put(key, list);
									}
									list[de1 + ds2] += value;  
								}
								if (de1 + de2 < maxDistance) { // +-
									String key = calculateKey(ci1, true, ci2, false);
									symmetryHash.put(key, calculateKey(ci2, true, ci1, false));
									
									double list[] = distanceHash_tmp.get(key);
									if (list == null) {
										list = new double[maxDistance];
										distanceHash_tmp.put(key, list);
									}
									list[de1 + de2] += value;  
								}
								if (ds1 + ds2 < maxDistance) { // -+
									String key = calculateKey(ci1, false, ci2, true);
									symmetryHash.put(key, calculateKey(ci2, false, ci1, true));

									double list[] = distanceHash_tmp.get(key);
									if (list == null) {
										list = new double[maxDistance];
										distanceHash_tmp.put(key, list);
									}
									list[ds1 + ds2] += value;  
								}
								if (ds1 + de2 < maxDistance) { // --
									String key = calculateKey(ci1, false, ci2, false);
									symmetryHash.put(key, calculateKey(ci2, true, ci1, true));
									
									double list[] = distanceHash_tmp.get(key);
									if (list == null) {
										list = new double[maxDistance];
										distanceHash_tmp.put(key, list);
									}
									list[ds1 + de2] += value;  
								}
							}
					}
			} while (true);
			br.close();
			
			//transfer and truncate data from distanceHash_tmp (double) to distanceHash (int)
			//make scores symmetric and do rounding
			ArrayList<String> missingKeys = new ArrayList<String>(); 
			for (String key : distanceHash_tmp.keySet()) {
				String sk = symmetryHash.get(key);
				if (!distanceHash_tmp.containsKey(sk))
					missingKeys.add(key);
			}
			for (String key : missingKeys) { // put missing symmetric keys to the hashes
				String sk = symmetryHash.get(key);
				distanceHash_tmp.put(sk, distanceHash_tmp.get(key));
				symmetryHash.put(sk, key);
			}
			
			for (String key : distanceHash_tmp.keySet()) {
				double list[] = distanceHash_tmp.get(key);
				double list2[] = distanceHash_tmp.get(symmetryHash.get(key)); // symmetric key

				int new_list[] = new int[maxDistance];
				
				double sum = 0.0;
				for (int i = 0; i < maxDistance; ++i) {
					sum += list[i] + list2[i];
					new_list[maxDistance - i - 1] = (int)(0.5 * sum + 0.5);  
				}
				int maxLength = maxDistance;
				while (maxLength > 0 && new_list[maxLength - 1] == 0)
					--maxLength;
				if (maxLength > 0)
					distanceHash.put(key, Arrays.copyOf(new_list, maxLength));
				
/*				int maxIndex = maxDistance - 1;
				double sum = 0.0;
				while (maxIndex >= 0 && sum + list[maxIndex] + list2[maxIndex] < 1.0) {
					sum += list[maxIndex] + list2[maxIndex];
					--maxIndex;
				}
				if (maxIndex >= 0) {
					int new_list[] = new int[maxIndex + 1];
					while (maxIndex >= 0) {
						sum += list[maxIndex] + list2[maxIndex];
						new_list[maxIndex] = (int)(0.5 * sum + 0.5); 
						--maxIndex;
					}
					distanceHash.put(key, new_list);
				}*/
			}

		} catch (Exception e) {
	     	e.printStackTrace();
			System.err.println("Error in file " + fn);
      	}
		
		System.err.println("Proximity links:"); 
		for (String key: distanceHash.keySet())
			System.err.println(key + "\t" +  distanceHash.get(key)[0]);
		return nonEmpty;
	}
	
	
	
	//"bin" a contig position
	private int bin(long position)
	{
		return (int) ((position - 1) / binSize); 
	}
	
	//how many bins a contig spans...
	public int binLength(ContigInterval ci)
	{
		int s = bin(ci.getStart()); 
		int e = bin(ci.getEnd());
		return (e - s + 1); 
	}

	private String calculateKey(ContigInterval c1, ContigInterval c2)
	{
		return calculateKey(c1, c1.getOrientation() >= 0, c2, c2.getOrientation() >= 0);
	}
	
	private String calculateKey(ContigInterval c1, boolean orientation1, ContigInterval c2, boolean orientation2)
	{
		return c1.toString() + (orientation1 ? '+':'-') + c2.toString() + (orientation2 ? '+':'-');
	}

	// calculate a list of contigs that should be evaluate against cis.get(ci)
	private ArrayList<Integer> next(ArrayList<ContigInterval> cis, int ci, int last) {
		ArrayList<Integer> ret = new ArrayList<Integer>();
		//int n = cis.size();
		int i = ci + 1;
		int d = 0;
		while (i <= last && d < maxDistance) {
			ContigInterval c = cis.get(i);
			ret.add(i);
			d += binLength(c);
			++i;
		}
		return ret;
	}

	// calculate a list of contigs that should be evaluate against cis.get(ci)
	private ArrayList<Integer> prev(ArrayList<ContigInterval> cis, int ci, int first) {
		ArrayList<Integer> ret = new ArrayList<Integer>();
		int i = ci - 1;
		int d = 0;
		while (i >= first && d < maxDistance) {
			ContigInterval c = cis.get(i);
			ret.add(i);
			d += binLength(c);
			--i;
		}
		return ret;
	}

	private int borderScore(ArrayList<ContigInterval> cis, ArrayList<Integer> prev, ArrayList<Integer> next){
		int ret = 0;
		int dp = 0;
		for (int pci : prev) {
			ContigInterval pc = ((pci < 0) ? cis.get(-pci - 1) : cis.get(pci)); 
			int dnp = dp;
			dp += binLength(pc);
			for (int nci : next) {
				if (dnp >= maxDistance)
					break;
				ContigInterval nc = ((nci < 0) ? cis.get(-nci - 1) : cis.get(nci)); 
				
				ret += score(pc, pci < 0, nc, nci < 0, dnp);
				dnp += binLength(nc); 
			}
		}
		return ret;
	}
	// difference in score if cis[start...end] are moved to new_position and possibly flipped
	public int scoreChange_slow(ArrayList<ContigInterval> cis, int start, int end, int moveDirection, boolean flip)
	{
		int scoreOld = score(cis);
		PlaceAndOrientContigs.changeOrder(cis, start, end, moveDirection, flip);
		int scoreNew = score(cis); 
		PlaceAndOrientContigs.changeOrderReverse(cis, start, end, moveDirection, flip);
		//System.err.println(scoreNew - scoreOld);
		return scoreNew - scoreOld;
	}
	
	// difference in score if cis[start...end] are moved to new_position and possibly flipped
	public int scoreChange(ArrayList<ContigInterval> cis, int start, int end, int moveDirection, boolean flip)
	{
		int ret = 0;
		int n = cis.size();
		if (moveDirection == 0) { // flip 
			if (flip) {
				//first part
				ArrayList<Integer> next = next(cis, start - 1, end); // start, start+1, ..., end
				ArrayList<Integer> prev = prev(cis, start, 0); // start-1, start-2, ..., 0
				
				ArrayList<Integer> new_next = new ArrayList<Integer>();
				int d = 0;
				for (int i = end; d < maxDistance && i >= start; --i) {
					ContigInterval nc = cis.get(i); 
					new_next.add(-i - 1);
					d += binLength(nc);
				}
				//subtract old score for start-1 -> start boundary, 
				ret -= borderScore(cis, prev, next);
				//...and add new score for this boundary
				ret += borderScore(cis, prev, new_next);
				
				//second part
				
				next = next(cis, end, n - 1); // end+1, end+2, ..., n-1
				prev = prev(cis, end + 1, start); // end, end-1, ..., start
				
				ArrayList<Integer> new_prev = new ArrayList<Integer>();
				d = 0;
				for (int i = start; d < maxDistance && i <= end; ++i) {
					ContigInterval nc = cis.get(i); 
					new_prev.add(-i - 1);
					d += binLength(nc);
				}
				//subtract old score for end -> end+1 boundary,
				ret -= borderScore(cis, prev, next);
				//..and add new score for this boundary
				ret += borderScore(cis, new_prev, next);

			} // else nothing to do...
		} else { // move one or multiple contigs...
			if (moveDirection >= 0) { // actually moveDirection > 0
				//int first = start 
				//int second = end; 
				//int third = end + moveDirection;
				//first part
				ArrayList<Integer> next = next(cis, start - 1, n - 1); // start, start+1, ..., n-1
				ArrayList<Integer> prev = prev(cis, start, 0); // start-1, start-2, ..., 0
				ret -= borderScore(cis, prev, next);
				
				ArrayList<Integer> new_next = new ArrayList<Integer>();
				int d = 0;

				for (int i = end + 1; d < maxDistance && i <= end + moveDirection; ++i) {
					ContigInterval nc = cis.get(i); 
					new_next.add(i);
					d += binLength(nc);
				}
				
				if (flip)
					for (int i = end; d < maxDistance && i >= start; --i) {
						ContigInterval nc = cis.get(i); 
						new_next.add(-i - 1);
						d += binLength(nc);
					}
				else
					for (int i = start; d < maxDistance && i <= end; ++i) {
						ContigInterval nc = cis.get(i); 
						new_next.add(i);
						d += binLength(nc);
					}

				for (int i = end + moveDirection + 1; d < maxDistance && i < n; ++i) {
					ContigInterval nc = cis.get(i); 
					new_next.add(i);
					d += binLength(nc);
				}
				//flip:  new_next = end + 1, ..., end + moveDirection,end...start, end + moveDirection + 1 ... n - 1
				//!flip: new_next = end + 1, ..., end + moveDirection,start...end, end + moveDirection + 1 ... n - 1
				//prev = start - 1, start - 2, ..., 0

				ret += borderScore(cis, prev, new_next);
				
				//second part
				next = next(cis, end, n - 1); // second, second+1, ..., n-1
				prev = prev(cis, end + 1, start); // second-1, second-2, ..., first

				ret -= borderScore(cis, prev, next);
				
				new_next = new ArrayList<Integer>();
				d = 0;
				if (flip)
					for (int i = end; d < maxDistance && i >= start; --i) {
						ContigInterval nc = cis.get(i); 
						new_next.add(-i - 1);
						d += binLength(nc);
					}
				else
					for (int i = start; d < maxDistance && i <= end; ++i) {
						ContigInterval nc = cis.get(i); 
						new_next.add(i);
						d += binLength(nc);
					}
				for (int i = end + moveDirection + 1; d < maxDistance && i < n; ++i) {
					ContigInterval nc = cis.get(i); 
					new_next.add(i);
					d += binLength(nc);
				}
				//flip:  new_next = end...start, end + moveDirection + 1 ... n - 1
				//!flip: new_next = start..end, ..., end + moveDirection + 1 ... n - 1

				ArrayList<Integer> new_prev = new ArrayList<Integer>();
				d = 0;

				for (int i = end + moveDirection; d < maxDistance && i > end; --i) { //TODO: check if (>= end or) > is proper
					ContigInterval nc = cis.get(i); 
					new_prev.add(i);
					d += binLength(nc);
				}
				//new_prev = end + moveDirection...end
				ret += borderScore(cis, new_prev, new_next);

				//third part
				next = next(cis, end + moveDirection, n - 1); // third, third+1, ..., n-1
				prev = prev(cis, end + moveDirection + 1, end); // third-1, third-2, ..., second

				ret -= borderScore(cis, prev, next);
				new_prev = new ArrayList<Integer>();
				
				//TODO: calculate new_prev
				d = 0;
				if (flip)
					for (int i = start; d < maxDistance && i <= end; ++i) {
						ContigInterval nc = cis.get(i); 
						new_prev.add(-i - 1);
						d += binLength(nc);
					}
				else
					for (int i = end; d < maxDistance && i >= start; --i) {
						ContigInterval nc = cis.get(i); 
						new_prev.add(i);
						d += binLength(nc);
					}
				ret += borderScore(cis, new_prev, next);
				//System.err.println(borderScore(cis, prev, next));
				//System.err.println(prev);
				//System.err.println(next);
				
			} else { //moveDirection < 0
				//int first = start + moveDirection; 
				//int second = start; 
				//int third = end;
				//first part
				ArrayList<Integer> next = next(cis, start + moveDirection - 1, n - 1); // first, first+1, ..., n-1
				ArrayList<Integer> prev = prev(cis, start + moveDirection, 0); // first-1, first-2, ..., 0
				ret -= borderScore(cis, prev, next);
				
				ArrayList<Integer> new_next = new ArrayList<Integer>();
				int d = 0;
				if (flip)
					for (int i = end; d < maxDistance && i >= start; --i) {
						ContigInterval nc = cis.get(i); 
						new_next.add(-i - 1);
						d += binLength(nc);
					}
				else
					for (int i = start; d < maxDistance && i <= end; ++i) {
						ContigInterval nc = cis.get(i); 
						new_next.add(i);
						d += binLength(nc);
					}
				
				for (int i = start + moveDirection; d < maxDistance && i < start; ++i) {
					ContigInterval nc = cis.get(i); 
					new_next.add(i);
					d += binLength(nc);
				}

				for (int i = end + 1; d < maxDistance && i < n; ++i) {
					ContigInterval nc = cis.get(i); 
					new_next.add(i);
					d += binLength(nc);
				}

				ret += borderScore(cis, prev, new_next);
				//System.err.println("");
				//System.err.println(prev);
				//System.err.println(next);
				//System.err.println("->");
				//System.err.println(prev);
				//System.err.println(new_next);
				
				//second part
				next = next(cis, start - 1, n - 1); // second, second+1, ..., n-1
				prev = prev(cis, start, start + moveDirection); // second-1, second-2, ..., first

				ret -= borderScore(cis, prev, next);
				
				new_next = new ArrayList<Integer>();

				d = 0;
				for (int i = start + moveDirection; d < maxDistance && i < start; ++i) {
					ContigInterval nc = cis.get(i); 
					new_next.add(i);
					d += binLength(nc);
				}

				for (int i = end + 1; d < maxDistance && i < n; ++i) {
					ContigInterval nc = cis.get(i); 
					new_next.add(i);
					d += binLength(nc);
				}

				ArrayList<Integer> new_prev = new ArrayList<Integer>();

				d = 0;
				if (flip)
					for (int i = start; d < maxDistance && i <= end; ++i) {
						ContigInterval nc = cis.get(i); 
						new_prev.add(-i - 1);
						d += binLength(nc);
					}
				else
					for (int i = end; d < maxDistance && i >= start; --i) {
						ContigInterval nc = cis.get(i); 
						new_prev.add(i);
						d += binLength(nc);
					}
				
				ret += borderScore(cis, new_prev, new_next);
				//System.err.println("");
				//System.err.println(prev);
				//System.err.println(next);
				//System.err.println("->");
				//System.err.println(new_prev);
				//System.err.println(new_next);

				//third part
				next = next(cis, end, n - 1); // third, third+1, ..., n-1
				prev = prev(cis, end + 1, start); // third-1, third-2, ..., second

				ret -= borderScore(cis, prev, next);
				new_prev = new ArrayList<Integer>();
				//TODO: calculate new_prev
				d = 0;
				for (int i = start - 1; d < maxDistance && i >= start + moveDirection; --i) {
					ContigInterval nc = cis.get(i); 
					new_prev.add(i);
					d += binLength(nc);
				}
				ret += borderScore(cis, new_prev, next);
				//System.err.println("");
				//System.err.println(prev);
				//System.err.println(next);
				//System.err.println("->");
				//System.err.println(new_prev);
				//System.err.println(next);
			}

		}
		//if (true) {
		//	int value = scoreChange_slow(cis, start, end, moveDirection, flip);
		//	if (value != ret)
		//		System.err.println(ret + "!=" + value + " " + start + " " + end + " " + moveDirection + " " + flip);
		//	return value;
		//}

		return ret;
	}
	
	public int score(ArrayList<ContigInterval> cis)
	{
		int ret = 0;
		for (int i = 1; i < cis.size(); ++i) {
			int j = i - 1;
			int d = 0;
			while (j >= 0 && d < maxDistance) {
				ret += score(cis.get(j), cis.get(i), d);
				d += binLength(cis.get(j));
				--j;
			}
		}
		return ret;
	}
	
	private int score(ContigInterval c1, ContigInterval c2, int distance)
	{
		int d[] = distanceHash.get(calculateKey(c1, c2));
		if (d != null && d.length > distance) 
			return d[distance];
		return 0;
	}

	public int score(String key, int distance)
	{
		int d[] = distanceHash.get(key);
		if (d != null && d.length > distance) 
			return d[distance];
		return 0;
	}

	private int score(ContigInterval c1, boolean flip1, ContigInterval c2, boolean flip2, int distance)
	{
		int d[] = distanceHash.get(calculateKey(c1, (c1.getOrientation() >= 0) ^ flip1, c2, (c2.getOrientation() >= 0) ^ flip2));
		if (d != null && d.length > distance) 
			return d[distance];
		return 0;
	}
	
	public static void main(String args[])
	{
		//Test Proximity class
		ArrayList<ContigInterval> bed = InputData.loadBed("proximity.bed");
		
		HashMap<String, ArrayList<ContigInterval>> bedHash2 = new HashMap<String, ArrayList<ContigInterval>>();
		for (ContigInterval ci : bed) {
			String contig = ci.getContig();
			if (!bedHash2.containsKey(contig))
				bedHash2.put(contig, new ArrayList<ContigInterval>());
			bedHash2.get(contig).add(ci);
		}
		Proximity p = new Proximity(10000, 2000, 0.01);
		p.loadData("ld.txt", bedHash2);
		ContigInterval c1 = bedHash2.get("000604F|quiver|pilon").get(0);
		ContigInterval c2 = bedHash2.get("000000F|quiver|pilon").get(0);
		ContigInterval c3 = bedHash2.get("000310F|quiver|pilon").get(0);
		ArrayList<ContigInterval> cis = new ArrayList<ContigInterval>();
		cis.add(c1);
		cis.add(c2);
		cis.add(c3);
		Collections.reverse(cis);
		for (int orient = 0; orient < 8; ++orient) {			
			c1.setOrientation((orient & 4) - 1);
			c2.setOrientation((orient & 2) - 1);
			c3.setOrientation((orient & 1) - 1);
			System.err.println(Integer.toBinaryString(orient + 8).substring(1) + ":\t" + p.score(cis));
		}
		System.err.println();
		Collections.reverse(cis);
		for (int orient = 0; orient < 8; ++orient) {			
			c1.setOrientation((orient & 4) - 1);
			c2.setOrientation((orient & 2) - 1);
			c3.setOrientation((orient & 1) - 1);
			System.err.println(Integer.toBinaryString(orient + 8).substring(1) + ":\t" + p.score(cis));
		}

		c1.setOrientation(1);
		c2.setOrientation(-1);
		c3.setOrientation(1);
		System.err.println(p.score(cis));
		System.err.println();
		for (int i = 0; i < 3; ++i) {
			System.err.println(p.scoreChange(cis, i, i, 0, true));
			System.err.println(p.scoreChange_slow(cis, i, i, 0, true));
			System.err.println();
		}

		for (int i = 0; i < 2; ++i) {
			System.err.println(p.scoreChange(cis, i, i, 1, true));
			System.err.println(p.scoreChange_slow(cis, i, i, 1, true));
			System.err.println();
		}
		System.err.println(p.scoreChange(cis, 0, 0, 2, true));
		System.err.println(p.scoreChange_slow(cis, 0, 0, 2, true));
		System.err.println();

		System.err.println(p.scoreChange(cis, 1, 2, -1, true));
		System.err.println(p.scoreChange_slow(cis, 1, 2, -1, true));
		System.err.println();

		
		System.err.println(p.scoreChange(cis, 0, 1, 1, true));
		System.err.println(p.scoreChange_slow(cis, 0, 1, 1, true));
		System.err.println();
		
	}
	// linkScore between ... ,prevC-1, prevC and nextC, nextC + 1, ...
	public int linkScore(ArrayList<ContigInterval> cis, int prevC, int nextC) {
		ArrayList<Integer> next = next(cis, nextC - 1, cis.size() - 1);
		ArrayList<Integer> prev = prev(cis, prevC + 1, 0);
		return borderScore(cis, prev, next);
	}

	public boolean linkedIntervals(ArrayList<ContigInterval> cis, int c1, int c2) {
		return (linkScore(cis, c1, c2) > 0);
	}
}
