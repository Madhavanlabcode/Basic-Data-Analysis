package schrodinger;

import java.io.File;

import main.SRAW;
import util.ArrayOps;
import util.FieldOps;
import util.color.ColorScales;
import util.fileops.Topomap;
import util.robot.Robo;

public class MovieMaker {
	static double hbar = 1, m = 1;
	
	public static String avidir = "D:\\Java\\movie maker\\";
	public static void main(String[] args)
	{
//		int N = 512;
//		double[][][] psi0 = new double[N][N][2];
//		
//		psi0[0][0][0] = 1;
//		int nframes = 2700;
//		double dt = 4/(Math.PI*nframes*2); //dt*nframes = (4/pi) corresponds to complete periodicity.
//		makeCompleteMovies(psi0, dt, nframes, "C:\\Program Files\\BMP to AVI\\");
		makeMovie("imp", 0, 395, 20, 100000);
	}

	
	static void makeMovie(double[][][] psi0, double dt, int nsteps, String dir)
	{
		int N = psi0.length;
		FreePEvolver f = new FreePEvolver(psi0, dt, hbar, m);
		System.out.println("time step " + f.dt + "; time constant " + f.tau + "; total time " + (f.dt*nsteps));
		double[][] mag = new double [N][N], phase = new double [N][N];
		double[][] magHat = new double [N][N];
		double[][] phaseHat = new double[N][N];
		
		double[][][] psiHatCentered = new double[N][N][2];
		double[][][] psiCentered = new double[N][N][2];
		f.getPsiCentered(psiCentered);
		FieldOps.magnitude(psiCentered, mag);
		FieldOps.phase(psiCentered, phase);
		f.getPsiHatCentered(psiHatCentered);
		FieldOps.magnitude(psiHatCentered, magHat);
		FieldOps.phase(psiHatCentered, phaseHat);
		
		ColorScales.MYC2d scale = new ColorScales.MYC2d(ArrayOps.max(mag), ArrayOps.min(mag), 2*Math.PI);
		ColorScales.LinearBRYW magscale = new ColorScales.LinearBRYW(ArrayOps.max(mag), ArrayOps.min(mag));
		ColorScales.MYC2d scalehat = new ColorScales.MYC2d(ArrayOps.max(magHat), 0, 2*Math.PI);
		
		double magmax, magmin;
		for (int i = 0; i < nsteps; i++)
		{
			SRAW.writeImage(dir + "psi" + fromInt(i) + ".bmp", mag, phase, scale);
			SRAW.writeImage(dir + "mag" + fromInt(i) + ".bmp", mag, magscale);
			SRAW.writeImage(dir + "psihat" + fromInt(i) + ".bmp", magHat, phaseHat, scalehat);
			
			
			f.incrementPsiHat();
			f.obtainPsi();
			f.getPsiCentered(psiCentered);
			FieldOps.magnitude(psiCentered, mag);
			FieldOps.phase(psiCentered, phase);
			f.getPsiHatCentered(psiHatCentered);
			FieldOps.phase(psiHatCentered, phaseHat);  //since it is a free particle, all 
			//magnitudes of the fourier transform are constant.
			magmax = ArrayOps.max(mag); magmin = ArrayOps.min(mag);
			scale.renormalize(magmax, magmin);
			magscale.renormalize(magmax, magmin);
		}
		SRAW.writeImage(dir + "psi" + fromInt(nsteps) + ".bmp", mag, phase, scale);
		SRAW.writeImage(dir + "mag" + fromInt(nsteps) + ".bmp", mag, magscale);
		SRAW.writeImage(dir + "psihat" + fromInt(nsteps) + ".bmp", magHat, phaseHat, scalehat);
	}
	
	static void makeCompleteMovies(double[][][] psi0, double dt, int nsteps, String dir)
	{
		Robo r = new Robo();
		makeMovie(psi0, dt, nsteps, dir);
		//Enter the following sentence into Command Prompt. We assume it is selected.

		r.typeString(BMPtoAVICommand("psi", 0, nsteps, 20));
		r.pressEnter();
		Robo.wait(900000); //wait 15 mins

		r.typeString(BMPtoAVICommand("mag", 0, nsteps, 20));
		r.pressEnter();
		Robo.wait(900000); //wait 15 mins

		r.typeString(BMPtoAVICommand("psihat", 0, nsteps, 20));
		r.pressEnter();
		Robo.wait(900000); //wait 15 mins
		
		//Delete the original bmps.
		File f = null;
		for (int i = 0; i < nsteps; i++)
		{
			f = new File(dir + "psihat" + fromInt(i) + ".bmp");
			if (f.exists())
				f.delete();
			f = new File(dir + "psi" + fromInt(i) + ".bmp");
			if (f.exists())
				f.delete();
			f = new File(dir + "mag" + fromInt(i) + ".bmp");
			if (f.exists())
				f.delete();
		}
	}
	
	public static String BMPtoAVICommand(String title, int start, int stop, int framerate)
	{
		return "EasyBMPtoAVI -filebase " + title + " -start " + start + " -end " + stop + " -framerate " + framerate + " -output " + title + ".avi";
	}
	public static void typeBMPtoAVICommand(String title, int start, int stop, int framerate)
	{
		Robo r = new Robo();
		Robo.wait(1000);
		r.typeString("EasyBMPtoAVI -filebase " + title + " -start " + start + " -end " + stop + " -framerate " + framerate + " -output " + title + ".avi");
		r.pressEnter();
	}
	
	public static String fromInt( int i)
	{
		if (i < 10) return "00000" + i;
		if (i < 100) return "0000" + i;
		if (i < 1000) return "000" + i;
		if (i < 10000) return "00" + i;
		if (i < 100000) return "0" + i;
		else return "" + i;
	}
	
	public static void makeMovie(Topomap t, boolean singleScale)
	{
		String title = "Topomap";
		if (singleScale)
		{
			ColorScales.LinearBRYW scale = new ColorScales.LinearBRYW(ArrayOps.max(t.data), ArrayOps.min(t.data));
			for (int i = 0; i < t.nlayers; i++)
			{
				SRAW.writeImage(avidir + title + fromInt(i), t.data[i], scale);
			}
			typeBMPtoAVICommand(title, 0, t.nlayers-1, 10);
		}
		else{
			for (int i = 0; i < t.nlayers; i++)
			{
				SRAW.writeImage(avidir + title + fromInt(i), t.data[i]);
			}
			typeBMPtoAVICommand(title, 0, t.nlayers-1, 10);
		}
		
		Robo.wait(300000);
		File f = null;
		for (int i = 0; i < t.nlayers; i++)
		{
			f = new File(avidir + title + fromInt(i) + ".bmp");
			if (f.exists())
				f.delete();
		}

	}
	public static void makeMovie(String title, int first, int last, int fps, int waitTimeMillis)
	{
		Robo.wait(3000);
		typeBMPtoAVICommand(title, first, last, fps);
		Robo.wait(waitTimeMillis);
		File f = null;
		for (int i = first; i <= last; i++)
		{
			f = new File(avidir + title + fromInt(i) + ".bmp");
			if (f.exists())
				f.delete();
		}

	}
}
