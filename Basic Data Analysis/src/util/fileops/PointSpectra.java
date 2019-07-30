package util.fileops;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.JFileChooser;

import util.ArrayOps;
import util.FieldOps;

public class PointSpectra {
	public double[][] data; //data is [nspec][nlayers]
	
	public double[] x;

	public double[] y;
	public int nspec;
	public double[] v;
//	public String[] names;
	public int nlayers;
	String extra;
	public double[] average;

	public PointSpectra(double[][] data, double[] v, double[] x, double[] y) {
		this.data = data;
		this.v = v;
		this.x = x;
		this.y = y;
		
		this.nspec = x.length;
		nlayers = v.length;
		average = new double [nlayers];
		for (int j = 0; j < nlayers; j++){
			for (int i = 0; i < nspec; i++)
				average[j] += data[i][j];
			average[j]/=nspec;
		}
	}
	public void resetAverage()
	{
		average = new double [nlayers];
		for (int j = 0; j < nlayers; j++){
			for (int i = 0; i < nspec; i++)
				average[j] += data[i][j];
			average[j]/=nspec;
		}
		
	}
	public String[] toLines(){
		String[] lines = new String[nlayers+2];
		lines[0] = "X";
		lines[1] = "Y";
		for (int j = 0; j < nspec; j++)
		{
			lines[0]+= "\t" + x[j];
			lines[1]+= "\t" + y[j];
		}
		for (int i = 0; i < nlayers; i++)
		{
			lines[i+2] = "" + v[i];
			for (int j = 0; j < nspec; j++)
			{
				lines[i+2] += "\t" + data[j][i];
			}
		}
		return lines;
	}
	public double[] getSpectrum(int i)
	{
		double[] s = new double [nlayers];
		for (int k = 0; k < nlayers; k++)
			s[k] = data[i][k];
		return s;
	}public int getByteLength()
	{
		int l = 0;
		l += 8;
		l += 8*(nspec + nspec + nlayers);
		l += 8*(nspec*nlayers);
		return l;
	}
	public PointSpectra copy()
	{
		return new PointSpectra(FieldOps.copy(data), FieldOps.copy(v), FieldOps.copy(x), FieldOps.copy(y));
	}
	//This is especially intended for ivsK spectra where v is k and y is energyy
	public Layer toLayer()
	{
		return new Layer(FieldOps.transpose(data), v, y, 1, 1);
	}
	public double evaluateAt(double voltage, int index)
	{
		double dv = (v[v.length-1] - v[0])/(v.length-1);
		double vi = (voltage - v[0])/dv;
		double vres = vi - (int)vi;
		return (1-vres)*data[index][(int)vi] + vres*data[index][(int)vi+1];
	}
	//the center is used because that is what appears in the tooltip in RHK
	//write nx, ny, nlateyers.
	//then x, y, v.
	//then data.
	//then names
	public static void writeBIN(PointSpectra t, String filepath)
	{
		if (!filepath.endsWith(".bin")) filepath += ".bin";
		File file = new File(filepath);
		
		FileOutputStream outf = null;
		BufferedOutputStream outbuff = null;
		DataOutputStream outd = null;
		
		try {
			outf = new FileOutputStream(file);
			outbuff = new BufferedOutputStream(outf, t.getByteLength());
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		outd = new DataOutputStream(outbuff);
		try {
			outd.writeInt(t.nspec);
			outd.writeInt(t.nlayers);
			for (int i = 0; i < t.nspec; i++)
			{
				outd.writeDouble(t.x[i]);
				outd.writeDouble(t.y[i]);
			}
			for (int i = 0; i < t.nlayers; i++)
				outd.writeDouble(t.v[i]);
			for (int i = 0; i < t.nlayers; i++){ //System.out.print(" " + i);
				for (int j = 0; j < t.nspec; j++)
					outd.writeDouble(t.data[j][i]);
			}
//			System.out.println();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try {
//			outbuff.flush();
//			outf.flush();
			outd.close();
			outbuff.close();
			outf.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	public static void writeContents(DataOutputStream outd, PointSpectra t)
	{
		try {
			outd.writeInt(t.nspec);
			outd.writeInt(t.nlayers);
			for (int i = 0; i < t.nspec; i++)
			{
				outd.writeDouble(t.x[i]);
				outd.writeDouble(t.y[i]);
			}
			for (int i = 0; i < t.nlayers; i++)
				outd.writeDouble(t.v[i]);
			for (int i = 0; i < t.nlayers; i++){ //System.out.print(" " + i);
				for (int j = 0; j < t.nspec; j++)
					outd.writeDouble(t.data[j][i]);
			}
//			System.out.println();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	public static void writeBIN(PointSpectra t)
	{
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
		{
			String s = fc.getSelectedFile().toString();
			if (!s.endsWith(".bin")) s += ".bin";
			PointSpectra.writeBIN(t, s);
		}
	}
	public static void writeBIN(PointSpectra t, JFileChooser fc)
	{
		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
		{
			String s = fc.getSelectedFile().toString();
			if (!s.endsWith(".bin")) s += ".bin";
			PointSpectra.writeBIN(t, s);
		}
	}
	
	public static void writeASCII(PointSpectra t, JFileChooser fc)
	{
		double[][] table = new double [t.nspec+1][t.nlayers+2];
		for (int i = 0; i < t.nlayers; i++)
			table[0][i+2] = t.v[i];
		for (int i = 0; i < t.nspec; i++)
			table[i+1][0] = t.x[i];
		for (int i = 0; i < t.nspec; i++)
			table[i+1][1] = t.y[i];
		for (int i = 0; i < t.nspec; i++)
			for (int j = 0; j < t.nlayers; j++)
				table[i+1][j+2] = t.data[i][j];
		
		FileOps.writeTableASCII(fc, table);

	}
	public static PointSpectra readBIN(String filepath)
	{
		int nspec, nlayers;
		double[] x = null, y = null, v = null;
		double[][] data = null;
		double[] origin = new double[2];
		File file = new File(filepath);
		FileInputStream inf = null;
		BufferedInputStream inbuff = null;
		DataInputStream ind = null;
		PointSpectra t = null;
		try {
			inf = new FileInputStream(file);
			inbuff = new BufferedInputStream(inf, (int)file.length()/4);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
		ind = new DataInputStream(inbuff);
		int temp;
		try {
			nspec = ind.readInt();
			nlayers = ind.readInt();
			x = new double[nspec]; y = new double [nspec]; v = new double [nlayers];
			data = new double [nspec][nlayers];
			
			for (int i = 0; i < nspec; i++){
				x[i] = ind.readDouble();
				y[i] = ind.readDouble();
			}
			for (int i = 0; i < nlayers; i++)
				v[i] = ind.readDouble();
			for (int i = 0; i < nlayers; i++){//System.out.print(" " + i);
				for (int k = 0; k < nspec; k++)
						data[k][i] = ind.readDouble();
			}
//			System.out.println();
			t = new PointSpectra(data, v, x, y);
		
			ind.close();
			inbuff.close();
			inf.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return t;
	}
	//the user opens the topomap from a file of his choosing.
	public static PointSpectra open()
	{
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		return open(fc);
	}
	public static PointSpectra open(JFileChooser fc)
	{
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
			return open(fc.getSelectedFile().toString());
		return null;
	}
	public static PointSpectra open(String path)
	{
		if (path.endsWith(".bin"))
			return PointSpectra.readBIN(path);
		else if (path.endsWith(".txt"))
			return RHKFileOps.getPointSpectra(new File(path));
		else if (path.endsWith(".asc"))
			return NanonisFileOps.getFromASC(path);
		else if (path.endsWith(".dat"))
			return NanonisFileOps.getFromDat(path);
		else;
		return null;
	}
	//This replaces the data block of a topomap while keeping the voltages and horizontal positions the same
	public static PointSpectra newPointSpectra(PointSpectra t, double[][] data)
	{
		return new PointSpectra(data, t.v, t.x, t.y);
	}
	public static PointSpectra truncateTo(PointSpectra t, int imin, int imax)
	{
		int dx = imax - imin;
		double[] v = new double [dx];
		double[][] data = new double [t.nspec][dx];
		for (int i = 0; i < dx; i++)
		{
			for (int j = 0; j < t.nspec; j++)
				data[j][i] = t.data[j][i+imin];
			v[i] = t.v[i+imin];
		}
		return new PointSpectra(data, v, t.x, t.y);
	}
	
	public PointSpectra[] splitByPosition()
	{
		int npos = 0;
		int[] nspecs;
		int[][] specindices;
		ArrayList<Double> xs = new ArrayList<Double>();
		ArrayList<Double> ys = new ArrayList<Double>();
		boolean posInList = false;
		for (int i = 0; i < nspec; i++)
		{
			posInList = false;
			for (int j = 0; j < xs.size() && !posInList; j++)
			{
				posInList = posInList || (x[i] == xs.get(j) && y[i] == ys.get(j)); 
			}
			if (!posInList)
			{
				xs.add(x[i]);
				ys.add(y[i]);
			}
		}
		
		npos = xs.size();
		nspecs = new int [npos];
		PointSpectra[] answer = new PointSpectra[npos];
		for (int j = 0; j < xs.size(); j++)
		{
			for (int i = 0; i < nspec; i++)
				if (x[i] == xs.get(j) && y[i] == ys.get(j))
					nspecs[j]++;
		}
		specindices = new int [npos][];
		int n = 0;
		for (int j = 0; j < npos; j++)
		{
			n = 0;
			specindices[j] = new int [nspecs[j]];
			for (int i = 0; i < nspec; i++)
				if (x[i] == xs.get(j) && y[i] == ys.get(j))
				{
						specindices[j][n++] = i;
				}
		}
		
		double[][] newx = new double [npos][];
		double[][] newy = new double [npos][];
		double[][][] newdata = new double [npos][][];
		for (int i = 0; i < npos; i++)
		{
			newx[i] = new double [nspecs[i]];
			newy[i] = new double [nspecs[i]];
			newdata[i] = new double [nspecs[i]][nlayers];
			for (int j = 0; j < nspecs[i]; j++)
			{
				newx[i][j] = x[specindices[i][j]];
				newy[i][j] = y[specindices[i][j]];
				for (int k = 0; k < nlayers; k++)
					newdata[i][j][k] = data[specindices[i][j]][k];
			}
			answer[i] = new PointSpectra(newdata[i], v, newx[i], newy[i]);
		}
		return answer;
	}
	
	public double[] getMetricCoords(int i) {
		return new double[] {x[i], y[i]};
	}
	
	public String[] toTextTable()
	{
		String firstLine = "Spec:";
		String secondLine = "X:";
		String thirdLine = "Y:";
		String[] lines = new String[nlayers+3];
		for (int i = 0; i < nspec; i++)
		{
			firstLine += "\t" + i;
			secondLine += "\t" + x[i];
			thirdLine += "\t" + y[i];
		}
		for (int i = 0; i < nlayers; i++){lines[i+3] = "" + v[i];
			for (int j = 0 ; j < nspec; j++)
			{
				lines[i+3] += "\t" + data[j][i]; 
			}
		}
		return lines;

	}
	public void flipV() {
		ArrayOps.flip(v);
		ArrayOps.flipY(data);
	}

}
