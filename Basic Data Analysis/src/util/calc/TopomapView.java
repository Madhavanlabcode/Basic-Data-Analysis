package util.calc;

import javax.swing.JFileChooser;

import drawing.FieldCalcDrawerP;
import drawing.GraphDrawerCart;
import util.ArrayOps;
import util.FieldOps;
import util.fileops.ColumnIO;
import util.fileops.RHKFileOps;
import util.fileops.Topomap;
import util.fourier.FFTOps;

public class TopomapView implements Calculation2D_1P{
	String info;
	
	int N;
	public int nlayers;
	int index = 0;
	double[][][] topomap;
	double[][][] fftmag;
	
	GraphDrawerCart spec;
	
	Topomap t;
	
	int current = 0;
	final int DISP_FIELDS = 2;
	
	public TopomapView(double[][][] topomap)
	{
		this.topomap = topomap;
		N = topomap[0].length;
		nlayers = topomap.length;
		fftmag = new double[nlayers][N][N];
		for (int i = 0; i < nlayers; i++)
		{
			FFTOps.obtainFFTmagCent(topomap[i], fftmag[i]);
			FieldOps.log(fftmag[i]);
		}
	}
	public TopomapView(Topomap topomap)
	{
		this.t = topomap;
		this.topomap = topomap.data;
		N = this.topomap[0].length;
		nlayers = this.topomap.length;
		fftmag = new double[nlayers][N][N];
		for (int i = 0; i < nlayers; i++)
		{
			FFTOps.obtainFFTmagCent(topomap.data[i], fftmag[i]);
			FieldOps.log(fftmag[i]);
		}
	}

	
	@Override
	public double[][] realField() {
		switch(current)
		{
		case 0:
			return topomap[index];
		case 1:
			return fftmag[index];
		}
		return null;
	}

	@Override
	public double[][][] complexField() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getInfo(int x, int y) {
		info = "r = " + "[" + x + ", " + y + "] + \r\n";
		info += "value = " + realField()[x][y];
		return info;
	}

	@Override
	public String getInfo(double x, double y) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double[] getMinMax(boolean real) {
		return new double [] {ArrayOps.min(realField()), ArrayOps.max(realField())};
	}

	@Override
	public void redoCalc(double p) {
		int oldindex = index;
		index = (int)p;
		if (oldindex != index)
			setSpecialCurrent();
	}

	@Override
	public void resize(double xmin, double xmax, double ymin, double ymax) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public double[] getFieldBounds() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void switchDisplayField() {
		current++;
		current %= DISP_FIELDS;
		setSpecialCurrent();
	}

	@Override
	public void save(String dir) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void run(String dir, String name) {
		// TODO Auto-generated method stub
	}
	private void setSpecialCurrent()
	{
		if (current == 1){
		}
	}
	
	public static TopomapView getNew(String filetable, String dir, String suffix)
	{
		String[] lines = ColumnIO.readString(filetable).split("\r\n");
		int n = lines.length - 1;
		double[][][] topomap = new double [n][][];
		for (int i = 0; i < n; i++)
		{
			topomap[i] = ColumnIO.readSquareTable(dir + lines[i+1] + suffix + ".dat");
			System.out.println(lines[i+1]);
		}
		return new TopomapView(topomap);
	}
	public static void main(String[] args)
	{
		String dir = "C:\\data\\analysis\\SrIrO 327\\241\\";
//		dir = "C:\\data\\lawlerhex\\flat8022010\\subloc\\";
//		TopomapView t = TopomapView.getNew(dir + "topomap.txt", dir + "", "");
		Topomap tm = null;
		JFileChooser fc = new JFileChooser(dir);
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
			if (fc.getSelectedFile().toString().endsWith(".bin"))
				tm = Topomap.readBIN(fc.getSelectedFile().toString());
			else if (fc.getSelectedFile().toString().endsWith(".txt"))
				tm = RHKFileOps.getTopoDIDV(fc.getSelectedFile());
		TopomapView t = new TopomapView(tm);
		new FieldCalcDrawerP(t, dir, new double[] {0, t.nlayers-1}, true, t.nlayers-1);
	}
	@Override
	public String getInfo(double p) {
		return String.format("%.1f", t.v[(int)p]*1000) + " mV";
	}

}
