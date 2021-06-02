/**
    This file is part of Lep-MAP.

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

	Copyright (C) 2013 Pasi Rastas, pasi.rastas@helsinki.fi, University of Helsinki
*/
// Union-Find data structure

public class UnionFind {
	
	private volatile int rank[];
	private volatile int uf[];
	
	public UnionFind(int size)
	{
		rank = new int[size];
		uf = new int[size];
		for (int i = 0; i < size; ++i)
			uf[i] = i;
	}
	public int find(int element) {
		int j = element;
		while (j != uf[j])
			j = uf[j];

		int i = element;
		while (i != uf[i]) {
			int t = uf[i];
			uf[i] = j;
			i = t;
		}

		return j;
	}
	public int union(int element1, int element2) {
		int i = find(element1);
		int j = find(element2);
		if (i != j) {
			if (rank[j] > rank[i])
				uf[i] = j;
			else if (rank[j] < rank[i])
				uf[j] = i;
			else {
				uf[j] = i;
				++rank[i];
			}
		}
		return uf[i];
	}
	
	public static void main(String args[])
	{
		final int size = 6000000;
		final UnionFind uf = new UnionFind(size);
		
		for (int i = 0; i < size / 2; ++i)
			uf.union(size / 2 + i, i/32);
		for (int i = 0; i < 64; ++i)
			System.err.println(uf.find(i) + "\t " + uf.rank[i]);
	}
}
