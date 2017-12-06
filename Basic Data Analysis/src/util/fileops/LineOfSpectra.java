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

public class LineOfSpectra extends PointSpectra{ //THis is really a list of spectra, not necessarily a line.
	

	public double[] z;

	public LineOfSpectra(double[][] data, double[] v, double[] x, double[] y, double[] z) {
		super(data, v, x, y);
		this.z = z;
		// TODO Auto-generated constructor stub
	}

	public static void writeBIN(LineOfSpectra t, String filepath)
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
		
		writeContents(outd, t);
		try {
			for (int i = 0; i < t.nspec; i++)
			outd.writeDouble(t.z[i]);
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

	public static LineOfSpectra readBIN(String filepath)
	{
		int nspec, nlayers;
		double[] x = null, y = null, v = null, z = null;
		double[][] data = null;
		double[] origin = new double[2];
		File file = new File(filepath);
		FileInputStream inf = null;
		BufferedInputStream inbuff = null;
		DataInputStream ind = null;
		LineOfSpectra t = null;
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
			x = new double[nspec]; y = new double [nspec]; v = new double [nlayers]; z = new double [nspec];
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
			
			for (int i = 0; i < nspec; i++)
				z[i] = ind.readDouble();
				//			System.out.println();
			t = new LineOfSpectra(data, v, x, y, z);
		
			ind.close();
			inbuff.close();
			inf.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return t;
	}
}

