package util.color;

public class ColorScaleHolder {

	private ColorScale1D cs;
	int currentCS = 1;
	
	public ColorScaleHolder(double min, double max)
	{
		resetColorScaleChange(min, max, false);
	}
	public void resetColorScaleChange(double min, double max, boolean oldBounds)
	{
		if (oldBounds)
			cs = ColorScales.getNew(cs.getMin(), cs.getMax(), currentCS);
		else cs = ColorScales.getNew(min, max, currentCS);
	}
	public void reBoundColorScale(double min, double max)
	{
		cs.renormalize(max, min);
	}
	
	public ColorScale1D getScale()
	{
		return cs;
	}
	public void setCurrentCS(int currentCS)
	{
		this.currentCS = currentCS;
		resetColorScaleChange(cs.getMin(), cs.getMax(), true);
	}
	
	public java.awt.Color getUnusedColor()
	{
		return ColorScales.getUnusedColor(currentCS);
	}
	public java.awt.Color getUnusedColor2()
	{
		return ColorScales.getUnusedColor2(currentCS);
	}
	public void incrementScaleIndex()
	{
		currentCS++;
		currentCS += ColorScales.NSCALES;
		currentCS %= ColorScales.NSCALES;
		resetColorScaleChange(0, 1, true);
	}
	public void decrementScaleIndex()
	{
		currentCS--;
		currentCS += ColorScales.NSCALES;
		currentCS %= ColorScales.NSCALES;
		resetColorScaleChange(0, 1, true);
	}
	
	public int getCurrentCS()
	{
		return currentCS;
	}

}
