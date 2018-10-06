package de.embl.cba.illuminationcorrection.basic;

import ij.ImagePlus;

public class BaSiCSettings
{

	public static final String FLAT_FIELD_ONLY = "Estimate flat-field only (ignore dark-field)";
	public static final String FLAT_FIELD_AND_DARK_FIELD = "Estimate both flat-field and dark-field";

	public static final String[] shadingEstimationOptions = {"Skip estimation and use predefined shading profiles","Estimate shading profiles"};
	public static final String [] shadingModelOptions = { FLAT_FIELD_ONLY, FLAT_FIELD_AND_DARK_FIELD };
	public static final String [] parameterSettingOptions = {"Automatic","Manual"};
	public static final String [] driftOptions = {"Ignore","Replace with zero","Replace with temporal mean"};
	public static final String [] correctionOptions = {"Compute shading and correct images","Compute shading only"};
	public static final String none = "None";

	public ImagePlus imp;
	public ImagePlus imp_flat;
	public ImagePlus imp_dark;
	public String myShadingEstimationChoice = shadingEstimationOptions[ 1 ];
	public String myShadingModelChoice = shadingModelOptions[ 0 ];
	public String myParameterChoice = parameterSettingOptions[ 0 ];
	public double lambda_flat = 0.5;
	public double lambda_dark = 0.5;
	public String myDriftChoice = driftOptions[ 0 ];
	public String myCorrectionChoice = correctionOptions[ 0 ] ;

}
