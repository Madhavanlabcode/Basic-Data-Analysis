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
import javax.swing.JOptionPane;

import schrodinger.MovieMaker;
import util.ArrayOps;
import util.FieldOps;
import util.fourier.FFTOps;
import util.geom.MapRectCoordSystem;

public class Topomap extends MapRectCoordSystem{
	
	public static String stddir = "C:\\data\\analysis\\SrIrO 327\\241_2\\";
	public static String currentDir = "";
//	public static String stddir = "C:\\";
	public boolean fileIsComplete = true;
	

	public double[][][] data;
	
	public double[] mean;
	
	public double[] v;
	public String[] names;
	public int nlayers;
	String extra;

	public Topomap(double[][][] data, double[] v, double[] x, double[] y, String[] names) {
		super(x, y);
		this.data = data;
		this.v = v;
		
		nlayers = v.length;
		if (names != null)
			this.names = names;
		else{
			this.names = new String [nlayers];
			for (int i = 0; i < nlayers; i++)
				this.names[i] = MovieMaker.fromInt(i);
		}
		mean = new double [nlayers];
		for (int i = 0; i < nlayers; i++)
			mean[i] = FieldOps.mean(data[i]);
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
		double dk = Math.abs(2*Math.PI/xLength); //assume isotropix
		for (int i = 0; i < length; i++)
			k[i] = i*dk;
		return k;
	}
	
	public double[] getSpectrum(int i, int j)
	{
		double[] s = new double [nlayers];
		for (int k = 0; k < nlayers; k++)
			s[k] = data[k][i][j];
		return s;
	}
	public double[] getAverageSpectrum()
	{
		double[] s = new double [nlayers];
		for (int k = 0; k < nlayers; k++)
			s[k] = FieldOps.mean(data[k]);
		return s;
	}
	public void putSpectrum(int i, int j, double[] s)
	{
		for (int k = 0; k < nlayers; k++)
			s[k] = data[k][i][j];
	}
	public int getByteLength()
	{
		int l = 0;
		l += 12;
		l += 8*(nx + ny + nlayers);
		l += 8*(nx*ny*nlayers);
		l += 4*nlayers;
		for (int i = 0; i < nlayers; i++)
			l += 2*names[i].length();
		l += 16; //for the origin
		return l;
	}
	public void flipY(){
		for (int i = 0; i < data.length; i++)
			ArrayOps.flipY(data[i]);
		y = ArrayOps.getFlip(y);
		resetCoords();
	}
	public double getValueAt(double[] point, int index)
	{
		return FieldOps.getValueAt(data[index], point[0], point[1], mean[index]);
	}
	public double getValueAt(double x, double y, int index)
	{
		return FieldOps.getValueAt(data[index], x, y, mean[index]);
	}
	public Layer getLayer(int i)
	{
		Layer f = new Layer(data[i], x, y, v[i], 0);
		f.fileIsComplete = fileIsComplete;
		return f;
	}
	public Layer getLayerPixels(int i)
	{
		double[] x = new double [nx], y = new double [ny];
		for (int j = 0; j < nx; j++)
			x[j] = j;
		for (int j = 0; j < ny; j++)
			y[j] = j;
		return new Layer(data[i], x, y, v[i], 0);
		
	}
	//the center is used because that is what appears in the tooltip in RHK
	//write nx, ny, nlateyers.
	//then x, y, v.
	//then data.
	//then names
	
	public double[][][] getFFT(boolean log)
	{
		double[][][] fftmag = new double [nlayers][nx][ny];
		for (int i = 0; i < nlayers; i++)
		{
			FFTOps.obtainFFTmagCent(data[i], fftmag[i]);
			if (log)FieldOps.log(fftmag[i]);
		}
		return fftmag;

	}
	public static void writeBIN(Topomap t, String filepath)
	{
		File file = new File(filepath.endsWith(".bin") ? filepath : filepath+".bin");
		putStdDir(file.getParent() + "\\");
		
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
			for (int i = 0; i < t.nx; i++)
				outd.writeDouble(t.x[i]);
			for (int i = 0; i < t.ny; i++)
				outd.writeDouble(t.y[i]);
			for (int i = 0; i < t.nlayers; i++)
				outd.writeDouble(t.v[i]);
			for (int i = 0; i < t.nlayers; i++){ System.out.print(" " + i);
				for (int j = 0; j < t.nx; j++)
					for (int k = 0; k < t.ny; k++)
						outd.writeDouble(t.data[i][j][k]);
			}
			System.out.println();
			for (int i = 0; i < t.nlayers; i++){
				outd.writeInt(t.names[i].length());
				outd.writeChars(t.names[i]);
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
	
	public static Topomap getCropped(Topomap t, int xi, int xf, int yi, int yf)
	{
		int nx = xf-xi;
		int ny = yf-yi;
		double[] x = new double[nx];
		double[] y = new double[ny];
		
		double[][][] data = new double [t.nlayers][nx][ny];
		
		for (int k = 0; k < t.nlayers; k++)
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
					data[k][i][j] = t.data[k][xi+i][yi+j];
		for (int i = 0; i < nx; i++)
			x[i] = t.x[i+xi];
		for (int i = 0; i < ny; i++)
			y[i] = t.y[i+yi];
		
		return new Topomap(data, t.v, x, y, null);
	}
	public static void writeBIN(Topomap t)
	{
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
		{
			String s = fc.getSelectedFile().toString();
			if (!s.endsWith(".bin")) s += ".bin";
			Topomap.writeBIN(t, s);
		}
	}
	public static void writeBIN(Topomap t, JFileChooser fc)
	{
		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
		{
			String s = fc.getSelectedFile().toString();
			if (!s.endsWith(".bin")) s += ".bin";
			Topomap.writeBIN(t, s);
		}
	}
	public static Topomap readBIN(String filepath)
	{
		int nx, ny, nlayers;
		double[] x = null, y = null, v = null;
		double[][][] data = null;
		String[] names = null;
		File file = new File(filepath);
		putStdDir(file.getParent() + "\\");
		FileInputStream inf = null;
		BufferedInputStream inbuff = null;
		DataInputStream ind = null;
		Topomap t = null;
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
			//if nLayers reasonable, proceed, if not, assume Layer format
			nlayers = ind.readInt();
			if(nlayers>256 || nlayers<1){
				System.out.println("nLayers>256, trying as layer");
				ind.close();
				inbuff.close();
				inf.close();
				Layer[] lays = {Layer.readBIN(filepath)};
				return Topomap.newTopomap(lays);
			}
			
			x = new double[nx]; y = new double [ny]; v = new double [nlayers];
			data = new double [nlayers][nx][ny];
			names = new String[nlayers];
			
			for (int i = 0; i < nx; i++)
				x[i] = ind.readDouble();
			for (int i = 0; i < ny; i++)
				y[i] = ind.readDouble();
			for (int i = 0; i < nlayers; i++)
				v[i] = ind.readDouble();
			for (int i = 0; i < nlayers; i++){System.out.print(" " + i);
				for (int j = 0; j < nx; j++)
					for (int k = 0; k < ny; k++)
						data[i][j][k] = ind.readDouble();
			}
			System.out.println();
			for (int i = 0; i < nlayers; i++){
				temp = ind.readInt();
				names[i] = "";
				for (int j = 0; j < temp; j++)
					names[i] += ind.readChar();
			}
			t = new Topomap(data, v, x, y, names);
		
			ind.close();
			inbuff.close();
			inf.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return t;
	}
	public static void writeBINsImages(Topomap t, String dir, boolean oneScale)
	{
		if (!oneScale)
			for (int i = 0; i < t.nlayers; i++)
			{
				RHKFileOps.write3Files(dir+"\\", dir + "\\fft\\", t.data[i], t.names[i]);
			}
	}
	public static void writeBINsImages(Topomap t, boolean oneScale)
	{
		JFileChooser fc = new JFileChooser(stddir);
		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION);
		String dir = fc.getCurrentDirectory().toString();
		if (!oneScale)
			for (int i = 0; i < t.nlayers; i++)
			{
				RHKFileOps.write3Files(dir+"\\", dir + "\\fft\\", t.data[i], t.names[i]);
			}
	}

	//the user opens the topomap from a file of his choosing.
	public static Topomap open()
	{
		JFileChooser fc = new JFileChooser(stddir);
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION){
			currentDir = fc.getCurrentDirectory().toString() + "\\";
			if (fc.getSelectedFile().toString().endsWith(".bin"))
				return Topomap.readBIN(fc.getSelectedFile().toString());
			else if (fc.getSelectedFile().toString().endsWith(".txt"))
				return RHKFileOps.getTopoDIDV(fc.getSelectedFile());
			else;
		}
		return null;
	}
	public static Topomap open(JFileChooser fc)
	{
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
			if (fc.getSelectedFile().toString().endsWith(".bin"))
				return Topomap.readBIN(fc.getSelectedFile().toString());
			else if (fc.getSelectedFile().toString().endsWith(".txt"))
				return RHKFileOps.getTopoDIDV(fc.getSelectedFile());
			else;
		return null;
	}
	public static Topomap openFromLayers(JFileChooser fc)
	{
		int nlayers = Integer.parseInt(JOptionPane.showInputDialog("Enter the nummber of layers."));
		Layer[] l = new Layer[nlayers];
		for (int i = 0; i < nlayers; i++)
			l[i] = Layer.open(fc);
		return Topomap.newTopomap(l);
	}
	public static Topomap openFromLayersAuto(JFileChooser fc)
	{
		int nlayers = Integer.parseInt(JOptionPane.showInputDialog("Enter the nummber of layers."));
		Layer[] l = new Layer[nlayers];
		l[0] = Layer.open(fc);
		String s = fc.getSelectedFile().toString();
		File base = new File(s.toString().substring(0, s.length() - 5)); //This assumes that the file ends with "0.bin"
		for (int i = 1; i < nlayers; i++)
			l[i] = Layer.readBIN(base.toString() + i + ".bin");
		return Topomap.newTopomap(l);
	}
	//This replaces the data block of a topomap while keeping the voltages and horizontal positions the same
	public static Topomap newTopomap(Topomap t, double[][][] data)
	{
		if (t != null)
			return new Topomap(data, t.v, t.x, t.y, null);
		else return new Topomap(data, ArrayOps.generateArrayInclBoth(0, data.length-1, data.length),
				ArrayOps.generateArrayInclBoth(0, data[0].length-1, data[0].length),
				ArrayOps.generateArrayInclBoth(0, data[0][0].length-1, data[0][0].length),
				null);
		
	}
	public static Topomap newTopomap(Layer[] l)
	{
		double[][][] data = new double [l.length][l[0].nx][l[0].ny];
		for (int i = 0; i < data.length; i++)
			data[i] = FieldOps.copy(l[i].data);
		
		double[] v = new double[l.length];
		for (int i = 0; i < l.length; i++)
			v[i] = l[i].v;
		
		double[] x = new double [l[0].nx];
		for (int i = 0; i < l[0].nx; i++)
			x[i] = l[0].x[i];
		double[] y = new double [l[0].ny];
		for (int i = 0; i < l[0].ny; i++)
			y[i] = l[0].y[i];
		
		
		return new Topomap(data, v, x, y, null);
	}
	public static void setStdDir()
	{
		String s = ColumnIO.readString("directory.txt");
		if (s != null) stddir = s.trim();
		else stddir = "C:\\";
	}
	public static void putStdDir(String dir)
	{
		stddir = dir;
		ColumnIO.writeString(dir, "directory.txt");
	}
	
}
