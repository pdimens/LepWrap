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

	Copyright (C) 2013-2016 Pasi Rastas, pasi.rastas@gmail.com, University of Helsinki, University of Cambridge
	
*/
import java.util.*;
import java.io.*;

public class ShortPath{

	private ArrayList<int[]> edges = new ArrayList<int[]>();
	private int numNodes;
	private ArrayList<String> nodeNames = new ArrayList<String>();

	private ArrayList<Integer> edgeWeights = new ArrayList<Integer>();
	private ArrayList<Integer> nodeWeights = new ArrayList<Integer>();
	private ArrayList<long[][]> prints = new ArrayList<long[][]>();
	
	private int alpha = 0;
	private int beta = 1;
	
	UnionFind uf = null;
	int cname[] = null;
	int size[] = null;
	int numComponents = 0;

	HashMap<String, Integer> hm = new HashMap<String, Integer>();

	public ShortPath(int alpha, int beta) {
		this.alpha = alpha;
		this.beta = beta;
	}

	private int addNode(String name, long segPrint[][], int weight)
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
	
	
	private void allPairsShortestPaths(String name1, String name2, ArrayList<Integer> nodeIndex, ArrayList<Integer> edgeIndex, int component)
	{
		int numNodes = nodeIndex.size();
		
		int sp[][] = new int[numNodes][numNodes];
		int sp2[][] = new int[numNodes][numNodes];

		int next[][] = new int[numNodes][numNodes];
		//int count[][] = new int[numNodes][numNodes];

		for (int i = 0; i < numNodes; ++i)
			for (int j = 0; j < numNodes; ++j) {
				sp[i][j] = (i == j) ? 0 : Integer.MAX_VALUE;
				sp2[i][j] = nodeWeights.get(nodeIndex.get(i));
			}

		for (int ei : edgeIndex) {
			int n1 = edges.get(ei)[0];
			int n2 = edges.get(ei)[1];
			
			int w = edgeWeights.get(ei);
			int nw =  nodeWeights.get(n1) + nodeWeights.get(n2);

			for (int i = 0; i < numNodes; ++i) {
				int ni = nodeIndex.get(i);
				if (ni == n1) {
					n1 = i;
					break;
				}
			}
			for (int i = 0; i < numNodes; ++i) {
				int ni = nodeIndex.get(i);
				if (ni == n2) {
					n2 = i;
					break;
				}
			}

			
			if (w == 0) {
				sp[n1][n2] = (n1 < n2) ? 0 : 1;
				sp[n2][n1] = (n1 < n2) ? 1 : 0;
			} else {
				sp[n1][n2] = w;
				sp[n2][n1] = w;
			}
			sp2[n1][n2] = nw;
			sp2[n2][n1] = nw;
			next[n1][n2] = n2;
			next[n2][n1] = n1;
		}
		
		for (int k = 0; k < numNodes; ++k)
			for (int i = 0; i < numNodes; ++i)
				if (sp[i][k] < Integer.MAX_VALUE)
					for (int j = 0; j < numNodes; ++j)
						if (sp[k][j] < Integer.MAX_VALUE) {
							int newValue = sp[i][k] + sp[k][j];
							int newValue2 = sp2[i][k] + sp2[k][j] - nodeWeights.get(nodeIndex.get(k));
							if (newValue < sp[i][j] || (newValue == sp[i][j] && newValue2 > sp2[i][j]) ) {									
								sp[i][j] = newValue;
								sp2[i][j] = newValue2;
								next[i][j] = next[i][k];
								//count[i][j] = 1;
							}/* else if (newValue == sp[i][j]) {
								++count[i][j];
								if (Math.random() < 1.0 / count[i][j])
									next[i][j] = next[i][k];
							}*/
						}

					
		int longestSP = 0;
		int maxI = 0;
		int maxJ = 0;
		int countLongest = 0;
		for (int i = 0; i < numNodes; ++i)
			for (int j = 0; j < numNodes; ++j)
				if (sp2[i][j] >= longestSP) {
					if (sp2[i][j] > longestSP)
						countLongest = 1;
					else
						++countLongest;
					if (Math.random() < 1.0 / countLongest) {
						maxI = i;
						maxJ = j;
					}
					longestSP = sp2[i][j];
				}
		if (name1 != null) {
			for (int i = 0; i < numNodes; ++i)
				if (nodeNames.get(nodeIndex.get(i)).equals(name1))
					maxI = i;
		}
		if (name2 != null) {
			for (int i = 0; i < numNodes; ++i)
				if (nodeNames.get(nodeIndex.get(i)).equals(name2))
					maxJ = i;
		}
		

		System.out.println("#Path length = " + sp[maxI][maxJ]);
		System.out.println("#Path weight = " + sp2[maxI][maxJ]);
		
		System.err.println("Heviest Shortest Path = " + longestSP);
		System.err.println("Number of solutions = " + countLongest);
		while (maxI != maxJ) {
			//System.out.println("^" + nodeNames.get(next[maxI][maxJ]) + "\t" + nodeNames.get(maxI));
			//System.out.println("^" + nodeNames.get(maxI) + "\t" + nodeNames.get(next[maxI][maxJ]));
			System.out.println(nodeNames.get(nodeIndex.get(maxI)) + "\t" + component + "\t" + nodeWeights.get(nodeIndex.get(maxI)));
			int oldMaxI = maxI;
			maxI = next[maxI][maxJ];
			if (oldMaxI != maxI) {
				int d[] = hammingDistance(prints.get(nodeIndex.get(maxI)), prints.get(nodeIndex.get(oldMaxI)));
				if (d[0] != 1)
					System.out.println("#distance = " + d[0]);
			}
		}
		System.out.println(nodeNames.get(nodeIndex.get(maxI)) + "\t" + component + "\t" + nodeWeights.get(nodeIndex.get(maxI)));
	}
	
	public void allPairsShortestPaths(String name1, String name2)
	{
		System.err.println("#Nodes = " + numNodes);
		System.err.println("#Edges = " + edges.size());
		
		ArrayList<ArrayList<Integer>> nodeIndex = new ArrayList<ArrayList<Integer>>(); 
		ArrayList<ArrayList<Integer>> edgeIndex = new ArrayList<ArrayList<Integer>>(); 
		
		for (int component = 0; component <= numComponents; ++component) {
			nodeIndex.add(new ArrayList<Integer>());
			edgeIndex.add(new ArrayList<Integer>());
		}
		
		for (int i = 0; i < numNodes; ++i) {
			int c = uf.find(i);
			nodeIndex.get(cname[c]).add(i);
		}

		for (int ei = 0; ei < edges.size(); ++ei) {
			int c = uf.find(edges.get(ei)[0]);
			edgeIndex.get(cname[c]).add(ei);
		}
			
			
		for (int component = 1; component <= numComponents; ++component) {
			System.err.println("#Component = " + component);
			System.err.println("#Nodes = " + nodeIndex.get(component).size());
			allPairsShortestPaths(name1, name2, nodeIndex.get(component), edgeIndex.get(component), component);
			
			
		}
	}
	
	private long[][] string2BitPrints(String s)
	{
		int n = s.length();
		long ret[][] = new long[2][(n+63) / 64];
		
		int index = 0;
		long bitIndex = 1; 
		
		for (int i = 0; i < n; ++i) {
			char c = s.charAt(i);
			if (c == '0' || c == '2' || c == '1') {
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
	
	private int[] hammingDistance(long p1[][], long p2[][])
	{
		int n = p1[0].length;
		int d1 = 0;
		int d2 = 0;
		for (int i = 0; i < n; ++i) {
            long mask = p1[1][i] & p2[1][i];
            long b1 = (p1[0][i] ^ p2[0][i]) & mask;
            d1 += Long.bitCount(b1);
            d2 += Long.bitCount(b1 ^ mask);
		}
		if (d1 >= d2)
			return new int[]{d2, d1};
		else
			return new int[]{d1, d2};
	}
	
	
	public boolean loadPrints(String fn, int type, int sizeLimit, int limit1, int limit2)
	{
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
				if (row.length >= 5 && Integer.parseInt(row[1]) == type && Integer.parseInt(row[2]) > sizeLimit)
					addNode(row[0], string2BitPrints(row[(type == 2) ? 4 : 3]), alpha + beta * Integer.parseInt(row[2]));
			}
		} while (true);

        br.close();
        
        uf = new UnionFind(numNodes);
        
        System.err.println("Computing pair-wise distances");
        for (int n1 = 0; n1 < numNodes; ++n1)
            for (int n2 = n1 + 1; n2 < numNodes; ++n2) {
            		int d[] = hammingDistance(prints.get(n1), prints.get(n2));
            		if (d[0] <= limit1 && d[1] >= limit2) {
            			edges.add(new int[]{n1, n2});
            			edgeWeights.add(Math.max(0, 4 * d[0] - 3));
            			//edgeWeights.add(d[0]);
            			uf.union(n1, n2);
            		}
            }

        connectedComponenets();
        
        
        return true;

      } catch (Exception e) {
			System.err.println("Error");
			return false;
      }
	}
	

	public static void main(String args[])
	{
		if (args.length == 0)
			Error.error(801);
		
		String extraParameters = "";
		System.out.print("#java ShortPath " + args[0]);
		for (int i = 1; i < args.length; ++i) {
			extraParameters += args[i] + " ";
			System.out.print(" " + args[i]);
		}
		System.out.println();
		
		ParameterParser pp = new ParameterParser();
		if (!pp.init(extraParameters))
			Error.error(801);
		pp.warning(new String[]{"informativeMask", "sizeLimit", "limit1", "limit2", "alpha", "beta", "begin", "end"});

		
		ShortPath sp = new ShortPath(Integer.parseInt(pp.getValueAsString("alpha", "0")), Integer.parseInt(pp.getValueAsString("beta", "1")));
		sp.loadPrints(args[0], Integer.parseInt(pp.getValueAsString("informativeMask", "1")), Integer.parseInt(pp.getValueAsString("sizeLimit", "1")), Integer.parseInt(pp.getValueAsString("limit1", "1")), Integer.parseInt(pp.getValueAsString("limit2", "1")));

		sp.allPairsShortestPaths(pp.getValueAsString("begin", null), pp.getValueAsString("end", null));

	}
}
