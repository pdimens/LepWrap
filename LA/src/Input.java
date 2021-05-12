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

	Copyright (C) 2013 Pasi Rastas, pasi.rastas@gmail.com, University of Helsinki
*/

//Reads files to ArrayLists of Strings 
import java.io.*;
import java.util.*;

public class Input {
  private Input() { };

  private static boolean keepComments = false;
  private static StringBuilder comments = new StringBuilder();
  public static String getComments()
  {
	  return comments.toString();
  }
  public static void setKeepComments(boolean value)
  {
	  keepComments = value;
  }
  
  public static String loadRow(BufferedReader br) throws Exception
  {
	  String s = "";
	  do { 
		  s = br.readLine();
		  //System.err.println(s);
		  if (s != null) {
			  int index = s.indexOf('#');
			  if (index >= 0) {
				  s = s.substring(0, index);
				  if (keepComments) {
					  comments.append(s.substring(index + 1));
					  comments.append('\n');
				  }
			  }
		  }
	  }	while ( s != null && (s.length() == 0) );
	  return s;
  }
  
  public static ArrayList<String> splitRow(String s, String delim)
  {
	  StringTokenizer st = new StringTokenizer(s, delim, false);
	  ArrayList<String> row = new ArrayList<String>(); 
	  while (st.hasMoreTokens()) {
		  String nt = st.nextToken();
		  row.add(nt);
	  }
	  return row;
  }
  
  
  public static ArrayList<String> loadTableRow(BufferedReader br, String delim) throws Exception
  {
	  String s = loadRow(br);
	  if (s == null)
		  return null;
	  
	  /*String sa[] = s.split(delim);
	  ArrayList<String> row = new ArrayList<String>();
	  for (String tmp : sa) 
		  row.add(tmp);*/
	  return splitRow(s, delim);
  }

  public static ArrayList<ArrayList<String>> loadTable(BufferedReader br, String delim) throws Exception {
	  ArrayList<ArrayList<String>> ret = new ArrayList<ArrayList<String>>();
	  //ArrayList<String> comments = new ArrayList<String>(); 
	  while (true) {
		  ArrayList<String> row =  loadTableRow(br, delim);
		  if (row == null)
			  break;
		  else
			  ret.add(row);
	  }
	  return ret;
  }  
  
  public static ArrayList<ArrayList<String>> loadTable(BufferedReader br, String delim, String returnDelim) throws Exception {
	  ArrayList<ArrayList<String>> ret = new ArrayList<ArrayList<String>>();
	  
	  //ArrayList<String> comments = new ArrayList<String>(); 
	  while (true) {
		  StringTokenizer st = null;
		  String s = null;
		  do { 
			  s = loadRow(br);
			  if (s != null)
				  st = new StringTokenizer(s, delim, true);
		  }	while ( s != null && !st.hasMoreTokens());
		  if (s == null)
			  break;
		  ArrayList<String> row = new ArrayList<String>(); 
		  while (st.hasMoreTokens()) {
			  String nt = st.nextToken();
			  if (delim.indexOf(nt) >= 0) {
				  if (returnDelim.indexOf(nt) >= 0)
					  row.add(nt);
			  } else
				  row.add(nt);
		  }
		  ret.add(row);
	  }
	  return ret;
  }
 
  public static ArrayList<ArrayList<String>> loadTable(String filename, String delim, String returnDelim) {
	    try {
	      BufferedReader br = new BufferedReader(new FileReader(filename));
	      ArrayList<ArrayList<String>> t = loadTable(br, delim, returnDelim);
	      br.close();
	      return t;

	    } catch (Exception e) {
	      System.err.println(e);
	      return null;
	    }
	  }
  
  
  public static ArrayList<ArrayList<String>> loadTable(String filename, String delim) {
	    try {
	      BufferedReader br = new BufferedReader(new FileReader(filename));
	      ArrayList<ArrayList<String>> t = loadTable(br, delim);
	      br.close();
	      return t;

	    } catch (Exception e) {
	      System.err.println(e);
	      return null;
	    }
	  }
  
  public static int[][] loadIntTable(String filename, String delim) {
	  	
	  ArrayList<int[]> ret = new ArrayList<int[]>(); 
	    try {
	      BufferedReader br = new BufferedReader(new FileReader(filename));
	      while (true) {
	    	  ArrayList<String> t = loadTableRow(br, delim);
	    	  if (t == null)
	    		  break;
	    	  else {
	    		  int it[] = new int[t.size()];
	    		  int i = 0;
	    		  for (String s : t)
	    			  it[i++] = Integer.parseInt(s);
	    		  ret.add(it);
	    	  }
	    		  
	      }
	      br.close();
	      return ret.toArray(new int[ret.size()][]);

	    } catch (Exception e) {
	      System.err.print(e);
	      return null;
	    }
	  }

  public static void main(String args[])
  {
	  ArrayList<ArrayList<String>> test = Input.loadTable("pedigree.txt", " \t");
	  for (ArrayList<String> row : test) {
		  for (String item : row)
			  System.err.print(item + " ");
		  System.err.println();
	  }
		  
	  
  }
 
}
