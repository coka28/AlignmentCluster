package core;

import java.io.BufferedReader;  
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

class SegmentRunner implements Runnable {
	
	int maxIndex,seqNr;
	double[] scoring;
	boolean[][] pointers;
	String[] strings;
	int[] strLengths;
	int[] factorials;
	int[] seqParts;
	double mt,mm,pt;
	int minLength;
	int[] starts;
	String thrdID;
	String[] seqPaths;
	private Thread thrd;
	
	int posIndex, threadID;
	
	SegmentRunner(int positionalIndex, int threadNr) {
		posIndex = positionalIndex;
		thrdID = String.valueOf(threadNr);
	}
	
	
	public void start() throws IOException {
		File build = new File("project.build");
		BufferedReader br = new BufferedReader(new FileReader(build));
		mt = Double.parseDouble(br.readLine());
		mm = -Double.parseDouble(br.readLine());
		pt = Double.parseDouble(br.readLine());
		minLength = Integer.parseInt(br.readLine());
		seqNr = 0;
		while (br.readLine() != null) seqNr++;
		br.close();
		build = new File("project.build");
		br = new BufferedReader(new FileReader(build));
		for (int i=0; i<4; i++) br.readLine();
		seqPaths = new String[seqNr];
		starts = new int[seqNr];
		strings = new String[seqNr];
		strLengths = new int[seqNr];
		factorials = new int[seqNr];
		seqParts = new int[seqNr];
		
		for (int i=0; i<seqNr; i++) seqPaths[i] = br.readLine();
		br.close();
		for (int i=0; i<seqNr; i++) {
			strLengths[i] = Integer.parseInt(seqPaths[i].substring(seqPaths[i].lastIndexOf("*")+1));
			seqParts[i] = Integer.parseInt(seqPaths[i].substring(seqPaths[i].lastIndexOf(";")+1,seqPaths[i].indexOf("*")));
			seqPaths[i] = seqPaths[i].substring(0,seqPaths[i].indexOf(";"));
		}

		File[] seqFiles = new File[seqNr];  
		BufferedReader[] brs = new BufferedReader[seqNr]; 
		for (int i=0; i<seqNr; i++) {
			seqFiles[i] = new File(seqPaths[i]);
			brs[i] = new BufferedReader(new FileReader(seqFiles[i]));
			strings[i] = "";
			String tmp = "";
			while ((tmp = brs[i].readLine()) != null) strings[i] += tmp;
			brs[i].close();
		}
		
		factorials[0] = 1;
		for (int i=1; i<strings.length; i++) 
			factorials[i] = factorials[i-1]*seqParts[i-1];
		
		for (int i=0; i<seqNr; i++) 
			starts[i] = ((posIndex/factorials[i])%seqParts[i])*(strLengths[i]-minLength+1);
	
		for (int i=0; i<seqNr; i++)
			strings[i] = ","+strings[i].substring(starts[i],starts[i]+strLengths[i]);
		
		for (int i=0; i<strings.length; i++)
			strLengths[i] = strings[i].length();
		
		maxIndex = strLengths[0];
		for (int i=1; i<seqNr; i++) {
			factorials[i] = factorials[i-1]*strLengths[i-1];
			maxIndex *= strLengths[i];
		}
		
		scoring = new double[maxIndex];
		pointers = new boolean[maxIndex][seqNr];
		
		thrd = new Thread(this,thrdID);
		thrd.start();
	}
	
	@Override
	public void run() {
		double min=0;
		for (int i=0; i<maxIndex; i++) {
			int[] thisPos = position(i);
			boolean[][] p = legitPointers(thisPos);
			boolean[] maxPointer = new boolean[thisPos.length];
			double mtc=match(thisPos), max=min;
			int maxindex = 0;
			double[] scores = new double[p.length];
			for (int k=0; k<p.length; k++) {
				int zeros = 0;
				for (int n=0; n<thisPos.length; n++)
					if (!p[k][n]) zeros++;
				double plty = zeros*pt;
				scores[k] = scoring[index(step(thisPos,p[k]))]+mtc-plty;
				if (scores[k]>=scores[maxindex]) {
					max = scores[k];
					maxPointer = p[k];
					maxindex = k;
				} else if (scores[k]<min) min = scores[k];
			}
			if (max<=0) {
				scoring[i] = 0;
				pointers[i] = null;
			} else {
				scoring[i] = max;
				pointers[i] = maxPointer;
			}
		}
		try {
			(new File("alignments."+String.valueOf(posIndex))).createNewFile();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		try {
			while (getOne(seqNr*3,5)); // this is a problem
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	boolean getOne(double minScore, int minLen) throws ClassNotFoundException, IOException {
		int maxindex = 0;
		for (int i=0; i<maxIndex; i++)
			if (scoring[maxindex]<=scoring[i])
				maxindex = i;
		
		if (scoring[maxindex] < minScore) return false;
		
		boolean[] currentpointer = pointers[maxindex];
		ArrayList<boolean[]> alignmentPointers = new ArrayList<boolean[]>();
		ArrayList<int[]> alignmentPositions = new ArrayList<int[]>();
		
		double newscore = 0;
		int currentindex = maxindex;
		int[] endpos = position(maxindex), startpos = new int[seqNr];
		for (int i=0; i<endpos.length; i++) endpos[i] += starts[i];
		boolean borderRegion = false;
		
		do {
			startpos = position(currentindex);
			scoring[currentindex] = 0;
			pointers[currentindex] = null;
			int zeros = 0;
			for (boolean i: currentpointer) if (!i) zeros++;
			alignmentPointers.add(currentpointer);
			alignmentPositions.add(position(currentindex));
			int[] gpos = position(currentindex);
			newscore += match(gpos)-zeros*pt; 
			for (int i=0; currentpointer!=null && i<currentpointer.length; i++) if (currentpointer[i]) {
				gpos[i]--;
				currentindex -= factorials[i];
			}
			currentpointer = pointers[currentindex];
		} while(currentpointer!=null);
		
		if (newscore < minScore || alignmentPositions.size() < minLen) return true;
		
		
		for (int i=0; i<alignmentPositions.size(); i++)
			for (int k=0; k<alignmentPositions.get(i).length; k++)
				alignmentPositions.get(i)[k] += starts[k];
		print(alignmentPositions);
		return true;
		
		/*
		for (int i=0; i<startpos.length; i++) {
			startpos[i] += starts[i];
			endpos[i]++;
		}
		
		String[] moddedstrings = new String[strings.length];
		for (int i=0; i<moddedstrings.length; i++) moddedstrings[i]="";
		
		int[] currentpos;
		for (int i=alignmentPointers.size()-1; i>=0; i--) {
			currentpointer = alignmentPointers.get(i);
			currentpos = alignmentPositions.get(i);
			alignmentPointers.remove(i);
			alignmentPositions.remove(i);
			for (int k=0; k<currentpos.length; k++) {
				if (currentpointer[k]) {
					moddedstrings[k] += strings[k].substring(currentpos[k], currentpos[k]+1);
				} else moddedstrings[k] += "-";
			}
		}
		
		print(moddedstrings, startpos, endpos, newscore);
		return true;
		**/
	}
	
	void print(ArrayList<int[]> pos) throws IOException, ClassNotFoundException {
		BufferedWriter fw = new BufferedWriter(new FileWriter("alignments."+String.valueOf(posIndex), true));
		
		int[] delta = new int[seqNr];
		String newlines = "";
		for (int i=0; i<pos.get(pos.size()-1).length && pos.size()>1; i++) {
			newlines += String.valueOf(pos.get(pos.size()-1)[i]);
			if (i==pos.get(pos.size()-1).length-1) newlines += "\n";
			else newlines += ";";
		}
		
		if (pos.size()-1 > 0)
			for (int i=0; i<pos.get(pos.size()-2).length; i++)
				delta[i] = pos.get(pos.size()-2)[i] - pos.get(pos.size()-1)[i];
		
		for (int i=pos.size()-2; i>=1; i--) {
			boolean changeOfDirection = false;
			for (int k=0; k<pos.get(i).length; k++) 
				if (delta[k] != pos.get(i-1)[k]-pos.get(i)[k])
					changeOfDirection = true;
				
			if (changeOfDirection)
				for (int k=0; k<pos.get(i).length; k++) {
					delta[k] = pos.get(i-1)[k]-pos.get(i)[k];
					newlines += String.valueOf(pos.get(i)[k]);
					if (k==pos.get(i).length-1) newlines += "\n";
					else newlines += ";";
				}
		}
		
		for (int i=0; i<pos.get(0).length; i++) {
			newlines += String.valueOf(pos.get(0)[i]);
			if (i==pos.get(0).length-1) newlines += "\n";
			else newlines += ";";
		}
		
		newlines += "\n";
		fw.write(newlines);
		fw.close();
	}
	
	void print(String[] alignedStrings, int[] startpos, int[] endpos, double maxscore) throws IOException, ClassNotFoundException {
		BufferedWriter fw = new BufferedWriter(new FileWriter("alignments."+String.valueOf(posIndex), true));
		
		for (int i=0; i<alignedStrings.length-1; i++)
			for (int k=i+1; k<alignedStrings.length; k++) {
				String[] tempstr = alignedStrings.clone();
				String str1 = "";
				String map = "";
				String str2 = "";
				for (int n=0; n<tempstr[i].length(); n++) if (tempstr[i].charAt(n)==tempstr[k].charAt(n)) {
					if (tempstr[i].charAt(n)=='-') {
						tempstr[i] = tempstr[i].substring(0, n) + tempstr[i].substring(n+1);
						tempstr[k] = tempstr[k].substring(0, n) + tempstr[k].substring(n+1);
						n--;
					} else {
						str1 += tempstr[i].substring(n,n+1);
						map += "|";
						str2 += tempstr[k].substring(n,n+1);
					}
				} else {
					str1 += tempstr[i].substring(n,n+1);
					map += "*";
					str2 += tempstr[k].substring(n,n+1);
				}
				String newlines = "seq"+String.valueOf(i+1)+" "+String.valueOf(startpos[i])+".."+String.valueOf(endpos[i])+", "+
						"seq"+String.valueOf(k+1)+" "+String.valueOf(startpos[k])+".."+String.valueOf(endpos[k])+", score : "+
						String.valueOf(maxscore+"\n");
				int strlength = str1.length();
				int n=0;
				while (str1.length()>70*n) {
					String str1part = str1.substring(70*n, Math.min(strlength, 70*(n+1)));
					String str2part = str2.substring(70*n, Math.min(strlength, 70*(n+1)));
					String mappart = map.substring(70*n, Math.min(strlength, 70*(n+1)));
					newlines += str1part+"\n"+mappart+"\n"+str2part+"\n";
					n++;
				}
				fw.write(newlines);
			}
		fw.write("\n");
			
		fw.close();
	}

	double match(int[] position) {
		double res = 0;
		for (int i=0; i<strings.length; i++)
			for (int k=i+1; k<strings.length; k++) {
				if (strings[i].charAt(position[i])==strings[k].charAt(position[k]))
					res += mt; 
				else res += mm; }
		return res;
	}
	
	int[] step(int[] a, boolean[] b) {
		int[] res = new int[a.length];
		for (int i=0; i<a.length; i++)
			if (b[i]) res[i] = a[i]-1; 
			else res[i] = a[i];
		return res;
	}
	
	int index(int[] position) {
		int res = 0;
		for (int i=0; i<position.length; i++) {
			int tmp = 1;
			for (int k=0; k<i; tmp*=strLengths[k], k++);
			res += tmp*position[i];
		}
		return res;
	}
	
	int[] position(int index) {
		int[] res = new int[strings.length];
		for (int i=0; i<res.length; i++) {
			int j = 1;
			for (int k=0; k<i; j*=strLengths[k], k++);
			int tmp = index/j;
			res[i] = tmp%strLengths[i];
		}
		return res;
	}
	
	boolean[][] legitPointers(int[] position) {
		int[] maxPos = new int[position.length];
		boolean[][] res;
		int n=0, m=1;
		for (int i=0; i<maxPos.length; i++) if (position[i]>0) {
			maxPos[n] = i;
			n++;
			m *= 2;}
		m--;
		res = new boolean[m][maxPos.length];
		int j = 1;
		for (int i=0; i<m; i++,j++) for (int k=0; k<n; k++)
			res[i][maxPos[k]] = ((int)(j/Math.pow(2, k))%2==1);
		return res;
	}
	
}

public class Exec {

	public static void main(String[] args) throws NumberFormatException, IOException {
		SegmentRunner[] thrds = new SegmentRunner[args.length];
		for (int i=0; i<args.length; i++) {
			thrds[i] = new SegmentRunner(Integer.parseInt(args[i]),i+1);
			thrds[i].start();
		}
	}

}
