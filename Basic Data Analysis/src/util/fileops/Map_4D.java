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

import javax.swing.JFileChooser;

import schrodinger.MovieMaker;
import util.FieldOps;
import util.fourier.FFTOps;
import util.geom.MapRectCoordSystem;

/**
 * This is intended to be exactly like Topomap except that multiple pieces of data can be associated with each 
 * (v, x, y) point.
 * The data array is arranged [npieces][nv][nx][ny] so that topomaps may be generated from data[i].
 * 
 * @author madhavanlab2011
 *
 */
public class Map_4D extends MapRectCoordSystem{
	
	public static String currentDir = "";
//	public static String stddir = "C:\\";
	

	public double[][][][] data;
	
	public double[][] mean;
	
	public double[] v;
	public String[] names;
	public String[] pieceNames;
	public int nlayers;
	
	public int npieces;
	String extra;

	public Map_4D(double[][][][] data, double[] v, double[] x, double[] y, String[] names, String[] pieceNames) {
		super(x, y);
		this.data = data;
		this.v = v;
		
		this.npieces = pieceNames.length;
		this.pieceNames = pieceNames;
		nlayers = v.length;
		if (names != null)
			this.names = names;
		else{
			this.names = new String [nlayers];
			for (int i = 0; i < nlayers; i++)
				this.names[i] = MovieMaker.fromInt(i);
		}
		mean = new double [npieces][nlayers];
		for (int i = 0; i < nlayers; i++)
			for (int j = 0; j < npieces; j++)
			mean[j][i] = FieldOps.mean(data[j][i]);
	}
	
	public Topomap getMap(int i)
	{
		return new Topomap(data[i], v, x, y, names);
	}
	public double[] getKValuesPixel(int length)
	{
		double[] k = new double [length];
		double dk = 2*Math.PI/x.length; //assume isotropix
		for (int i = 0; i < length; i++)
			k[i] = i*dk;
		return k;
	}
	public double[] getKValuesReal(int length)
	{
		double[] k = new double [length];
		double dx = x[1]-x[0];
		double L = dx*x.length;
		double dk = 2*Math.PI/L; //assume isotropix
		for (int i = 0; i < length; i++)
			k[i] = i*dk;
		return k;
	}
	
	public double[] getSpectrum(int i, int j, int pieceIndex)
	{
		double[] s = new double [nlayers];
		for (int k = 0; k < nlayers; k++)
			s[k] = data[pieceIndex][k][i][j];
		return s;
	}
	public void putSpectrum(int i, int j, double[] s, int pieceIndex)
	{
		for (int k = 0; k < nlayers; k++)
			s[k] = data[pieceIndex][k][i][j];
	}
	public void putSpectra(int i, int j, double[][] s)
	{
		for (int h = 0; h < npieces; h++)
			for (int k = 0; k < nlayers; k++)
				s[h][k] = data[h][k][i][j];
	}
	public int getByteLength()
	{
		int l = 0;
		l += 16; //the 4 size parameters nx ny nv and npieces
		l += 8*(nx + ny + nlayers); //the doubles which are the arrays x, y, and v.
		l += 8*(nx*ny*nlayers*npieces); //the dataset
		l += 4*nlayers; //the int for the length of the layername
		for (int i = 0; i < nlayers; i++)
			l += 2*names[i].length(); //the layername chars
		l += 4*npieces;
		for (int i = 0; i < npieces; i++)
			l += 2*names[i].length(); //the pieceName chars
		
		
		l += 16; //for the origin
		return l;
	}
	public double getValueAt(double[] point, int index, int pieceIndex)
	{
		return FieldOps.getValueAt(data[pieceIndex][index], point[0], point[1], mean[pieceIndex][index]);
	}
	public Layer getLayer(int i, int pieceIndex)
	{
		return new Layer(data[pieceIndex][i], x, y, v[i], 0);
	}
	public Layer getLayerPixels(int i, int pieceIndex)
	{
		double[] x = new double [nx], y = new double [ny];
		for (int j = 0; j < nx; j++)
			x[j] = j;
		for (int j = 0; j < ny; j++)
			y[j] = j;
		return new Layer(data[pieceIndex][i], x, y, v[i], 0);
		
	}
	//the center is used because that is what appears in the tooltip in RHK
	//write nx, ny, nlateyers.
	//then x, y, v.
	//then data.
	//then names
	
	public double[][][] getFFT(boolean log, int pieceIndex)
	{
		double[][][] fftmag = new double [nlayers][nx][ny];
		for (int i = 0; i < nlayers; i++)
		{
			FFTOps.obtainFFTmagCent(data[pieceIndex][i], fftmag[i]);
			if (log)FieldOps.log(fftmag[i]);
		}
		return fftmag;

	}
	public static void writeBIN(Map_4D t, String filepath)
	{
		File file = new File(filepath);
		
		FileOutputStream outf = null;
		BufferedOutputStream outbuff = null;
		DataOutputStream outd = null;
		
		try {
			outf = new FileOutputStream(file);
			outbuff = new BufferedOutputStream(outf, t.getByteLength()/4);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		outd = new DataOutputStream(outbuff);
		try {
			outd.writeInt(t.nx);
			outd.writeInt(t.ny);
			outd.writeInt(t.nlayers);
			outd.writeInt(t.npieces);
			for (int i = 0; i < t.nx; i++)
				outd.writeDouble(t.x[i]);
			for (int i = 0; i < t.ny; i++)
				outd.writeDouble(t.y[i]);
			for (int i = 0; i < t.nlayers; i++)
				outd.writeDouble(t.v[i]);
			for (int h = 0; h < t.npieces; h++){System.out.print(" " + h);
				for (int i = 0; i < t.nlayers; i++)
					for (int j = 0; j < t.nx; j++)
						for (int k = 0; k < t.ny; k++)
							outd.writeDouble(t.data[h][i][j][k]);
			}
			System.out.println();
			for (int i = 0; i < t.nlayers; i++){
				outd.writeInt(t.names[i].length());
				outd.writeChars(t.names[i]);
			}
			for (int i = 0; i < t.npieces; i++)
			{
				outd.writeInt(t.pieceNames[i].length());
				outd.writeChars(t.pieceNames[i]);
			}
			
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
	public static void writeBIN(Map_4D t)
	{
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
		{
			String s = fc.getSelectedFile().toString();
			if (!s.endsWith(".bin")) s += ".bin";
			Map_4D.writeBIN(t, s);
		}
	}
	public static void writeBIN(Map_4D t, JFileChooser fc)
	{
		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
		{
			String s = fc.getSelectedFile().toString();
			if (!s.endsWith(".bin")) s += ".bin";
			Map_4D.writeBIN(t, s);
		}
	}
	public static Map_4D readBIN(String filepath)
	{
		int nx, ny, nlayers, npieces;
		double[] x = null, y = null, v = null;
		double[][][][] data = null;
		String[] names = null;
		String[] pieceNames = null;
		double[] origin = new double[2];
		File file = new File(filepath);
		FileInputStream inf = null;
		BufferedInputStream inbuff = null;
		DataInputStream ind = null;
		Map_4D t = null;
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
			nx = ind.readInt();
			ny = ind.readInt();
			nlayers = ind.readInt();
			npieces = ind.readInt();
			x = new double[nx]; y = new double [ny]; v = new double [nlayers];
			data = new double [npieces][nlayers][nx][ny];
			names = new String[nlayers];
			pieceNames = new String[npieces];
			
			for (int i = 0; i < nx; i++)
				x[i] = ind.readDouble();
			for (int i = 0; i < ny; i++)
				y[i] = ind.readDouble();
			for (int i = 0; i < nlayers; i++)
				v[i] = ind.readDouble();
			for (int h = 0; h < npieces; h++){System.out.print(" " + h);
				for (int i = 0; i < nlayers; i++)
					for (int j = 0; j < nx; j++)
						for (int k = 0; k < ny; k++)
							data[h][i][j][k] = ind.readDouble();
			}
			System.out.println();
			for (int i = 0; i < nlayers; i++){
				temp = ind.readInt();
				names[i] = "";
				for (int j = 0; j < temp; j++)
					names[i] += ind.readChar();
			}
			for (int i = 0; i < npieces; i++){
				temp = ind.readInt();
				pieceNames[i] = "";
				for (int j = 0; j < temp; j++)
					pieceNames[i] += ind.readChar();
			}
			t = new Map_4D(data, v, x, y, names, pieceNames);
		
			ind.close();
			inbuff.close();
			inf.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return t;
	}
	public static void writeBINsImages(Map_4D t, String dir, boolean oneScale)
	{
		if (!oneScale)
			for (int h = 0; h < t.npieces; h++)
				for (int i = 0; i < t.nlayers; i++)
				{
					RHKFileOps.write3Files(dir+"\\", dir + "\\fft\\", t.data[h][i], t.pieceNames[h] + "_" + t.names[i]);
				}
	}
	public static void writeBINsImages(Map_4D t, boolean oneScale)
	{
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION);
		String dir = fc.getCurrentDirectory().toString();
		if (!oneScale)
			for (int h = 0; h < t.npieces; h++)
				for (int i = 0; i < t.nlayers; i++)
				{
					RHKFileOps.write3Files(dir+"\\", dir + "\\fft\\", t.data[h][i], t.pieceNames[h] + "_" + t.names[i]);
				}
	}

	//the user opens the topomap from a file of his choosing.
	public static Map_4D open()
	{
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION){
			currentDir = fc.getCurrentDirectory().toString() + "\\";
			return Map_4D.readBIN(fc.getSelectedFile().toString());
		}
		return null;
	}
	public static Map_4D open(JFileChooser fc)
	{
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
			return Map_4D.readBIN(fc.getSelectedFile().toString());
		return null;
	}
	//This replaces the data block of a topomap while keeping the voltages and horizontal positions the same
	public static Map_4D newTopomap(Map_4D t, double[][][][] data)
	{
		return new Map_4D(data, t.v, t.x, t.y, t.names, t.pieceNames);
	}
	public static void setStdDir()
	{
		String s = ColumnIO.readString("directory.txt");
		if (s != null) Topomap.stddir = s.trim();
	}
	
}
