package util.movie;

import java.io.File;

import util.color.ColorScale;
import util.robot.Robo;
import main.SRAW;

public abstract class BasicMovie {

	public static final String DIR = "C:\\Users\\madhavanlab2011\\Documents\\BMP to AVI\\";
	public static final int waitTime = 60000; //1 minutes to make the movie before the files are deleted.
	
	//to save memory:
	public double[][] getRField(double p){return null;}
	public double[][][] getCField(double p){return null;}
	public void makeComplexMovie(double pmin, double pmax, int npts, String picturename, boolean minzero, boolean log)
	{
		double[][][] x;
		double p;
		for (int i = 0; i < npts; i++)
		{
			p = pmin + i*(pmax-pmin)/(npts-1);
			x = getCField(p);
			SRAW.writeImage(DIR + picturename + fromInt(i), x, minzero, log);
		}
		typeBMPtoAVICommand(picturename, 0, npts-1, 20);
		Robo.wait(waitTime);
		
		File f;
		for (int i = 0; i < npts; i++)
		{
			f = new File(DIR + picturename + fromInt(i) + ".bmp");
			if (f.exists()) f.delete();
		}
	}
	public void makeRealMovie(double pmin, double pmax, int npts, String picturename, boolean log)
	{
		double[][] x;
		double p;
		for (int i = 0; i < npts; i++)
		{
			p = pmin + i*(pmax-pmin)/(npts-1);
			x = getRField(p);
			SRAW.writeImage(DIR + picturename + fromInt(i), x, log);
		}
		typeBMPtoAVICommand(picturename, 0, npts-1, 20);
		Robo.wait(waitTime);
		
		File f;
		for (int i = 0; i < npts; i++)
		{
			f = new File(DIR + picturename + fromInt(i) + ".bmp");
			if (f.exists()) f.delete();
		}
	}
	public static String fromInt(int i)	{
		if (i < 10) return "00000" + i;
		if (i < 100) return "0000" + i;
		if (i < 1000) return "000" + i;
		if (i < 10000) return "00" + i;
		if (i < 100000) return "0" + i;
		else return "" + i;
	}
	public static void typeBMPtoAVICommand(String title, int start, int stop, int framerate)
	{
		Robo r = new Robo();
		Robo.wait(1000);
		r.typeString("EasyBMPtoAVI -filebase " + title + " -start " + start + " -end " + stop + " -framerate " + framerate + " -output " + title + ".avi");
		r.pressEnter();
	}
	
	public static void writeFrames(double[][][] data, ColorScale scale, String picturename, int framerate)
	{
		for (int i = 0; i < data.length; i++)
			SRAW.writeImage(DIR + picturename + fromInt(i), data[i], scale, true);
		typeBMPtoAVICommand(picturename, 0, data.length-1, framerate);
		Robo.wait(waitTime);
		
		File f;
		for (int i = 0; i < data.length; i++)
		{
			f = new File(DIR + picturename + fromInt(i) + ".bmp");
			if (f.exists()) f.delete();
		}
		
	}
}
