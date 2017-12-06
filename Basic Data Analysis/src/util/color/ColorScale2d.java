package util.color;

import java.awt.Color;

public interface ColorScale2d {

	public Color of(double a, double b);

	public void renormalize(double scalemin, double scalemax);
}
