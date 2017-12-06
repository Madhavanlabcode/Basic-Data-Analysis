package util.color;

import java.awt.Color;
import java.awt.image.BufferedImage;

//Linear 1D color scale
public abstract class ColorScale1D implements ColorScale{

	Color[] colors;
	int nc;
	double max, min;
	double delta;
	
	public ColorScale1D(Color[] colors, double min, double max)
	{
		this.colors = colors;
		this.nc = colors.length;
		this.min = min;
		this.max = max;
		delta = max - min;
	}
	public Color of(double x) {
		int c = (int)((x - min)*nc/delta);
		c = Math.min(nc-1, c);
		c = Math.max(0, c);
		return colors[c];
	}
	public void renormalize(double max, double min)
	{
		this.max = max;
		this.min = min;
		delta = max - min;
	}

	public double getMax() {
		return max;
	}
	public void setMax(double max) {
		this.max = max;
	}
	public double getMin() {
		return min;
	}
	public void setMin(double min) {
		this.min = min;
	}
	public String toString()
	{
		return "Min: " + min + "; Max: " + max;
	}
	
	public java.awt.image.BufferedImage getScaleImage(int ratio)
	{
		java.awt.image.BufferedImage image = new java.awt.image.BufferedImage(nc/ratio, nc, BufferedImage.TYPE_INT_RGB);
		for (int i = 0; i < nc; i++)
			for (int j = 0; j < nc/ratio; j++)
				image.setRGB(j, i, colors[i].getRGB());
		return image;
	}
	
}
