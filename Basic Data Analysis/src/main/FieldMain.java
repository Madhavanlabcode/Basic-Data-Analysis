package main;

import java.awt.image.BufferedImage;

import util.FieldOps;

import drawing.*;

public class FieldMain {

	public double[][] data;
	public FieldDrawer drawer;
	FieldPanel panel;
	DataZoomPanel dataPanel;
	
	public FieldMain()
	{
		String name = "";
		String txt = "C:\\MinGW\\data\\veccomps\\";
		double[][] real = load(txt + name + "X3.txt");
		double[][] imag = load(txt + name + "Y3.txt");
		
		double[][] real2 = load(txt + name + "X1.txt");
		double[][] im2 = load(txt + name + "Y1.txt");
		
		
		int N = 512;
////		double[][] imag = new double[N][N];
//		
		
		FieldOps.plusEquals(real, -512);
		FieldOps.plusEquals(imag, -512);

		double[][] data2 = new double[256][256];
//		data = DataManip.magnitude(real, imag);
		data = DataManip.phase(real, imag);
		FieldOps.plusEquals(real2, -512);
		FieldOps.plusEquals(im2, -512);
		data2 = DataManip.phase(real2, im2);
		for (int i = 0; i < 256; i++)
			for (int j = 0; j < 256; j++)
			{
				data[i][j] -= data2[i][j];
				data[i][j] = -data[i][j];
			}
//		double[][] phase = DataManip.phase(real, imag);
//		
//		data = load(txt + name + ".txt");
//		for (int i = 0; i < N; i++)
//			for(int j = 0; j < N; j++)
//				data[i][j] = Math.log(data[i][j] + .01);

		dataPanel = new DataZoomPanel(this);
// 		drawer = new FieldDrawer(dataPanel.field, phase);
		drawer = new FieldDrawer(dataPanel.field);
		drawer.panel = dataPanel;
		panel = new FieldPanel(drawer);
		saveBMP(txt + name + "th13.bmp", drawer.getImage());

	}
	
	public FieldDrawer drawer()
	{
		return drawer;
	}
	public static void main(String[] args){
		new FieldMain();
	}
	
	public static double[][] load(String filepath)
	{
		return SRAW.getData(filepath);
	}
	static void saveBMP(String path, double[][] data)
	{
		SRAW.writeImage(path, data);
	}
	static void saveBMP(String path, BufferedImage image)
	{
		SRAW.writeImage(path, image);
	}
}
