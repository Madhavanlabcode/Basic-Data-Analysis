package util.color;

import java.awt.Color;

public interface ColorScale {

	public Color of(double x);
	public void renormalize(double max, double min);
	public double getMax();
	public double getMin();
}