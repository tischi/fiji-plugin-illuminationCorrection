package de.embl.cba.illuminationcorrection.basic;

import ij.IJ;

import de.embl.cba.illuminationcorrection.Log;

import net.imagej.DatasetService;
import net.imagej.ops.OpService;

import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import java.io.File;


@Plugin(type = Command.class, menuPath = "Plugins>Restoration>Illumination>BaSiC" )
public class BaSiCCommand implements Command
{
	@Parameter
	public UIService uiService;

	@Parameter
	public DatasetService datasetService;

	@Parameter
	public LogService logService;

	@Parameter
	public OpService opService;

	@Parameter
	public StatusService statusService;

	@Parameter( label = "Input file" )
	public File inputFile;

	BaSiCSettings settings = new BaSiCSettings();

	@Parameter( choices = { BaSiCSettings.FLAT_FIELD_ONLY, BaSiCSettings.FLAT_FIELD_AND_DARK_FIELD } )
	public String myShadingModelChoice = settings.myShadingModelChoice;

	public void run()
	{
		setSettingsFromUI();

		final BaSiC baSiC = new BaSiC( settings );
		baSiC.run();

		Log.info( "Done!" );

	}

	public boolean acceptFile( String fileNameEndsWith, String file )
	{
		final String[] fileNameEndsWithList = fileNameEndsWith.split( "," );

		for ( String endsWith : fileNameEndsWithList )
		{
			if ( file.endsWith( endsWith.trim() ) )
			{
				return true;
			}
		}

		return false;
	}


	public void setSettingsFromUI()
	{
		settings.imp = IJ.openImage( inputFile.getAbsolutePath() );
		settings.imp.show();
		settings.myShadingModelChoice = myShadingModelChoice;
	}


}
