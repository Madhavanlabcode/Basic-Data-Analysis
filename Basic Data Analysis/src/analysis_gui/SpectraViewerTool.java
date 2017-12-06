package analysis_gui;

import drawing.SpectraDrawer;

public abstract class SpectraViewerTool {
	SpectraDrawer sd;
	
	public SpectraViewerTool(SpectraDrawer sd)
	{
		this.sd = sd;
	}
	
	public abstract void processKeyStroke(char c);
	public abstract void processSpectra();
}
