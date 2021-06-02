/**
    This file is part of Lep-MAP3.

    Lep-MAP is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Lep-MAP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Lep-MAP.  If not, see <http://www.gnu.org/licenses/>.

	Copyright (C) 2017 Pasi Rastas, pasi.rastas@gmail.com, University of Helsinki
	
*/


import java.util.*;
import java.io.*;
//TODO: Add support for missing phased haplotypes (-) 
public class LMPlot{

	private int numNodes;
	private ArrayList<String> nodeNames = new ArrayList<String>();
	private ArrayList<Integer> nodeWeights = new ArrayList<Integer>();
	private ArrayList<long[][][]> prints = new ArrayList<long[][][]>();
	
	UnionFind uf = null;
	int cname[] = null;
	int size[] = null;
	int numComponents = 0;
	int familySize[];
	
	boolean selfingPhase = false;
	boolean outputAll = false;
	
	HashMap<String, Integer> hm = new HashMap<String, Integer>();

	public LMPlot() {
	}

	private int addNode(String name, long segPrint[][][], int weight)
	{
		if (!hm.containsKey(name)) {
			hm.put(name, numNodes);
			nodeNames.add(name);
			nodeWeights.add(weight);
			prints.add(segPrint);
			return numNodes++;
		} else
			return hm.get(name);
	}

	
	public void  connectedComponenets() {
		cname = new int[numNodes];
		size = new int[numNodes];
        numComponents = 0;
        for (int n = 0; n < numNodes; ++n) {
                int c = uf.find(n);
                if (cname[c] == 0)
                	cname[c] = ++numComponents;
                ++size[cname[c]];
        }
        System.out.println("#Connected components = " + numComponents);
	}
	
	private long[][] string2BitPrints(String s)
	{
		int n = s.length();
		long ret[][] = new long[2][(n+63) / 64];
		
		int index = 0;
		long bitIndex = 1; 
		
		for (int i = 0; i < n; ++i) {
			char c = s.charAt(i);
			if (c == '0' || c == '1') {
				ret[1][index] |= bitIndex;
				if (c == '1')
					ret[0][index] |= bitIndex;
			}
				
			bitIndex *= 2;
			if (bitIndex == 0) {
				bitIndex = 1;
				++index;
			}
		}
		return ret;
	}
	
	//mask out if both p1 and p2 are heterozygote
	private int[] hammingDistance3(long p1[][][], long p2[][][])
	{
		int ret[] = new int[3];
		// store different individual...
		int difFam = -1;
		for (int f = 0; f < p1.length; ++f) {
			int n = p1[f][0].length;
			int d1 = 0;
			int d2 = 0;
			for (int i = 0; i < n; ++i) {
	            long mask1 = (p1[f][0][i] ^ p1[f ^ 1][0][i]); // 1 if heterozygote
	            long mask2 = (p2[f][0][i] ^ p2[f ^ 1][0][i]); // 1 if heterozygote
	            long mask12 = ~(mask1 & mask2);  // 0 if both heterozygote

	            long mask = p1[f][1][i] & p2[f][1][i] & mask12;
	            long b1 = (p1[f][0][i] ^ p2[f][0][i]) & mask;
	            d1 += Long.bitCount(b1);
	            d2 += Long.bitCount(b1 ^ mask);
			}
			if (d1 == 1 || d2 == 1)
				difFam = f;
			if (d1 >= d2) {
				ret[0] += d2;
				ret[1] += d1;
			} else {
				ret[0] += d1;
				ret[1] += d2;
			}
		}
		
		if (ret[0] == 1) {
			int n = p1[difFam][0].length;
			for (int i = 0; i < n; ++i) {
	            long mask = p1[difFam][1][i] & p2[difFam][1][i];
	            long b1 = (p1[difFam][0][i] ^ p2[difFam][0][i]) & mask;
	            int d1 = Long.bitCount(b1);
	            int d2 = Long.bitCount(b1 ^ mask);
	            if (d1 == 1 || d2 == 1) {
		            if (d2 == 1)
		            	b1 = b1 ^ mask;
		            
		            ret[2] = familySize[difFam/2] + i * 64 + Long.numberOfTrailingZeros(b1) + 1;
		            return ret;
	            }
			}
			
		}
		return ret;
	}

	private int[] hammingDistance(long p1[][][], long p2[][][])
	{
		if (selfingPhase)
			return hammingDistance3(p1, p2);
		int ret[] = new int[3];
		// store different individual...
		int difFam = -1;
		for (int f = 0; f < p1.length; ++f) {
			int n = p1[f][0].length;
			int d1 = 0;
			int d2 = 0;
			for (int i = 0; i < n; ++i) {
	            long mask = p1[f][1][i] & p2[f][1][i];
	            long b1 = (p1[f][0][i] ^ p2[f][0][i]) & mask;
	            d1 += Long.bitCount(b1);
	            d2 += Long.bitCount(b1 ^ mask);
			}
			if (d1 == 1 || d2 == 1)
				difFam = f;
			if (d1 >= d2) {
				ret[0] += d2;
				ret[1] += d1;
			} else {
				ret[0] += d1;
				ret[1] += d2;
			}
		}
		
		if (ret[0] == 1) {
			int n = p1[difFam][0].length;
			for (int i = 0; i < n; ++i) {
	            long mask = p1[difFam][1][i] & p2[difFam][1][i];
	            long b1 = (p1[difFam][0][i] ^ p2[difFam][0][i]) & mask;
	            int d1 = Long.bitCount(b1);
	            int d2 = Long.bitCount(b1 ^ mask);
	            if (d1 == 1 || d2 == 1) {
		            if (d2 == 1)
		            	b1 = b1 ^ mask;
		            
		            ret[2] = familySize[difFam/2] + i * 64 + Long.numberOfTrailingZeros(b1) + 1;
		            return ret;
	            }
			}
			
		}
			
		return ret;
	}
	
	public String normalize(String s){
				return "";
	}
	
	public boolean loadPrints(String fn, int limit2)
	{
		HashMap<String, Integer> printMap = new HashMap<String, Integer>();
		int numPrints = 0;
		
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
			if (line.charAt(0) != '#') {
				String[] row = line.split("\t");
				StringBuilder sb = new StringBuilder();
				for (int i = 4; i < row.length; ++i)
					sb.append(row[i].substring(0, row[i].indexOf(" ")));
				String str = sb.toString();

				if (!printMap.containsKey(str)) {
					printMap.put(str, numPrints++);
					long ps[][][] = new long[(row.length - 4) * 2][][];
					if (familySize == null) {
						familySize = new int[row.length - 4 + 1];
						for (int i = 4; i < row.length; ++i)
							familySize[i - 3] = familySize[i - 4] + row[i].substring(0, row[i].indexOf(" ")).length() / 2; 
					}
					for (int i = 4; i < row.length; ++i) {
						String str12 = row[i].substring(0, row[i].indexOf(" "));
						int mid = str12.length() / 2;
						ps[(i-4) * 2] =  string2BitPrints(str12.substring(0, mid));
						ps[(i-4) * 2 + 1] =  string2BitPrints(str12.substring(mid));
					}
					addNode("" + numPrints, ps, 1);					
				} else {
					int n = printMap.get(str);
					nodeWeights.set(n, nodeWeights.get(n) + 1);
				}
			}
		} while (true);
		System.err.println("#nodes=" + numNodes);

        br.close();

        if (outputAll)
	        for (int n1 = 0; n1 < numNodes; ++n1)
	            for (int n2 = n1 + 1; n2 < numNodes; ++n2) {
	            	int d[] = hammingDistance(prints.get(n1), prints.get(n2));
	            	System.err.println((n1 + 1) + "\t" + (n2 + 1) + "\t" + d[0] + "\t" + d[1]);
	            }

        
        uf = new UnionFind(numNodes);
        
        System.err.println("Computing pair-wise distances");
        System.out.println("graph g {");
        System.out.println("node [fontsize=40,penwidth=4]");
        System.out.println("edge [penwidth=3]");
        
        for (int n1 = 0; n1 < numNodes; ++n1) {
        	double area = Math.sqrt(nodeWeights.get(n1)) / 30.0;
			System.out.println(nodeNames.get(n1) + "[width=" + area + ",height=" + area + ",fixedsize=true]");
        }
        boolean dir = Math.random() < 0.5;
    	
        int minDist[] = new int[numNodes];
        int union[] = new int[numNodes];
        while (true) {
        	boolean allDone = true;
        	Arrays.fill(minDist, Integer.MAX_VALUE);
        	//calculate minDist
	        for (int n1 = 0; n1 < numNodes; ++n1) {
	        	int u1 = uf.find(n1);
	            for (int n2 = n1 + 1; n2 < numNodes; ++n2) {
	            	int u2 = uf.find(n2);
            		if (u1 != u2) { 
            			int d[] = hammingDistance(prints.get(n1), prints.get(n2));
	            		if (d[1] >= limit2) {
	            			allDone = false;
	            			minDist[u1] = Math.min(minDist[u1], d[0]);
	            			minDist[u2] = Math.min(minDist[u2], d[0]);
	            		}
            		}
	            }
	        }
	        if (allDone)
	        	break;

	        for (int n1 = 0; n1 < numNodes; ++n1) {
	        	int u = uf.find(n1);
	        	union[n1] = u;
	        }

	        for (int n1 = 0; n1 < numNodes; ++n1) {
	        	int u1 = union[n1];
	        	int md1 = minDist[u1];
	            for (int n2 = n1 + 1; n2 < numNodes; ++n2) {
	            	int u2 = union[n2];
            		if (u1 != u2) { 
            			int d[] = hammingDistance(prints.get(n1), prints.get(n2));
	            		if (d[1] >= limit2 && (d[0] <= md1 || d[0] <= minDist[u2])) {
	            			uf.union(u1, u2);
	            			
	            			if (Math.random() < 0.1)
	            				dir = !dir;
	            			
	            			String color = "";
	            			try {
	            				int name1 = Integer.parseInt(nodeNames.get(n1));
	            				int name2 = Integer.parseInt(nodeNames.get(n2));
	            				if (Math.abs(name1 - name2) > 1)
	            					color = "color=red";
	            			} catch (Exception e)
	            			{
	            			}
            				int limit = d[0];
            				if (limit == 0) // add edge for distance 0 as well...
            					++limit;
	            			if (dir)
		            			for (int i = 0; i < limit; ++i)
		            				System.out.println(nodeNames.get(n1) + "--" + nodeNames.get(n2) + "[" + ((d[0] == 1) ? "label=" + d[2] + "": "") + color + "]");
	            			else
		            			for (int i = 0; i < limit; ++i)
		            				System.out.println(nodeNames.get(n2) + "--" + nodeNames.get(n1) + "[" + ((d[0] == 1) ? "label=" + d[2] + "": "") + color + "]");
	            		}
            		}
	            }
	        }
        }
        System.out.println("}");
        return true;
      } catch (Exception e) {
    	  	e.printStackTrace();
			System.err.println("Error");
			return false;
      }
	}
	

	public static void main(String args[])
	{
		if (args.length == 0)
			Error.error(1201);
		
		String extraParameters = "";
		System.out.print("#java LMPlot " + args[0]);
		for (int i = 1; i < args.length; ++i) {
			extraParameters += args[i] + " ";
			System.out.print(" " + args[i]);
		}
		System.out.println();
		
		ParameterParser pp = new ParameterParser();
		if (!pp.init(extraParameters))
			Error.error(801);
		pp.warning(new String[]{"limit2", "selfingPhase", "outputAllDistances"});

		
		LMPlot sp = new LMPlot();
		sp.setOutputAll(pp.getValueAsString("outputAllDistances", "0").equals("1"));
		sp.setSelfingPhase(pp.getValueAsString("selfingPhase", "0").equals("1"));
		sp.loadPrints(args[0], Integer.parseInt(pp.getValueAsString("limit2", "1")));
		//sp.printDot();
	}

	private void setOutputAll(boolean value) {
		outputAll = value;
	}

	private void setSelfingPhase(boolean value) {
		selfingPhase = value;
	}
}
