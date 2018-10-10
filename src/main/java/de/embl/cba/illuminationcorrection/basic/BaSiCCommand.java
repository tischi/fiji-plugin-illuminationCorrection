package de.embl.cba.illuminationcorrection.basic;

import ij.IJ;

import de.embl.cba.illuminationcorrection.Log;

import ij.ImagePlus;
import ij.io.FileSaver;
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
	public static final String CORRECTED_SUFFIX = "-corrected.tif";
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

	@Parameter( label = "Output directory", style = "directory")
	public File outputDirectory;

	BaSiCSettings settings = new BaSiCSettings();

	@Parameter( choices = { BaSiCSettings.FLAT_FIELD_ONLY, BaSiCSettings.FLAT_FIELD_AND_DARK_FIELD } )
	public String myShadingModelChoice = settings.myShadingModelChoice;

	public void run()
	{
		setSettingsFromUI();

		final BaSiC baSiC = new BaSiC( settings );
		baSiC.run();
		final ImagePlus correctedImage = baSiC.getCorrectedImage();
		saveCorrectedImage( correctedImage );
	}

	private void saveCorrectedImage( ImagePlus correctedImage )
	{
		final FileSaver fileSaver = new FileSaver( correctedImage );
		final String outputPath = outputDirectory + File.separator + inputFile.getName() + CORRECTED_SUFFIX;
		fileSaver.saveAsTiff( outputPath );
		Log.info( "Saved file to: " + outputPath );
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
		settings.myShadingModelChoice = myShadingModelChoice;
	}


}
