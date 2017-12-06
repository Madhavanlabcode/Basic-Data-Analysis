package util;

import image.GifSequenceWriter;
import image.ImageEditing;

import java.awt.Graphics;
import java.awt.image.BufferedImage;

import util.color.ColorScale;
import util.color.ColorScales;
import util.fileops.Layer;
import util.fileops.Topomap;
import util.fourier.FFTOps;

//This is intended to make nice images of things for export
public class ExportUtil {
	
	
	//The main thing is that the aspect ratio will be 4:3. So the larger image will be 3:3, the smaller
	//will be 1:1, and the other 2:1 will be reserved for energy labels.
	//This 4:3 aspect ratio assumes the topomap is square
	public static void exportGIFForPPTSlide(Topomap t, Layer topography, ColorScale fft, ColorScale main, boolean showColorScale, boolean fftIsBigger, boolean fftLog, int delay, String path, int scaleInt)
	{
		BufferedImage[] output = new BufferedImage[t.nlayers];
		double[][] fftData;
		
		BufferedImage imfft, immain;
		BufferedImage bigger, smaller;
		bigger = new BufferedImage(3*t.nx, 3*t.ny, BufferedImage.TYPE_INT_RGB);
		smaller = new BufferedImage(t.nx, t.ny, BufferedImage.TYPE_INT_RGB);
		BufferedImage topo;
		Graphics g;
		boolean mainWasNull = false;
		
		for (int i = 0; i < t.nlayers; i++)
		{
//			output[i] = new BufferedImage(4*t.nx, 3*t.ny, BufferedImage.TYPE_INT_RGB);
				
			
			if (main == null || mainWasNull)
			{
				mainWasNull = true;
				main = ColorScales.getNew(FieldOps.min(t.data[i]), FieldOps.max(t.data[i]), scaleInt);
			}
			immain = ImageEditing.getBufferedImage(t.data[i], main);
			fftData = FFTOps.obtainFFTmagCent(t.data[i]);
			if (fftLog) FieldOps.log(fftData);
			if (fft == null) fft = ColorScales.getNew(fftData, scaleInt);
			imfft = ImageEditing.getBufferedImage(fftData, fft);
			if (fftIsBigger){
				ImageEditing.enlargeBasic(imfft, bigger, 3);
				smaller = immain;
			}
			else
			{
				ImageEditing.enlargeBasic(immain, bigger, 3);
				smaller = imfft;
			}
			output[i] = ImageEditing.weldHorizontal(bigger, smaller);
			g = output[i].getGraphics();
			g.setColor(java.awt.Color.WHITE);
			g.fillRect(3*t.nx, t.ny, t.nx, 3*t.ny);
			g.setColor(java.awt.Color.BLACK);
			g.drawString(String.format("%.1f", t.v[i]*1000) + " mV", 3*t.nx + 10, t.ny + 50);
			
			//Map characteristics:
			String[] info = new String [3];
			info[0] = ("" + String.format("%.1f", t.xLength*1e9) + " nm x " + String.format("%.1f", t.yLength*1e9) + " nm");
			info[1] = "" + t.nx + " x " + t.ny + " pixels";
			info[2] = "" + t.nlayers + " Layers; " + String.format("%.1f", t.v[0]*1000) + " mV to " + String.format("%.1f", t.v[t.nlayers-1]*1000) + " mV";
			for (int j = 0; j < info.length; j++)
				g.drawString(info[j], 3*t.nx + 10, t.ny + 90 + 20*j);
			
		
			if (topography != null)
			{
				g.drawString("Topography", 3*t.nx + 50, 2*t.ny - 20);
				topo = ImageEditing.getBufferedImage(topography.data, scaleInt);
				g.drawImage(topo, 3*t.nx, 2*t.ny, null);
			}
			
		}
//		for (int i = 0; i < 20; i++)
			GifSequenceWriter.writeGifSequence(output, delay/10, path.endsWith(".gif") ? path : path + ".gif");
	}
}
