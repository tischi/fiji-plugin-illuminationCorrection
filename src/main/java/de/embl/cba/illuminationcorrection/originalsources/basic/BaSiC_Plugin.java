package de.embl.cba.illuminationcorrection.originalsources.basic;

import java.util.*;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NewImage;
import ij.macro.Interpreter;
import ij.measure.ResultsTable;
import ij.plugin.ContrastEnhancer;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.plugin.frame.Recorder;
import ij.plugin.filter.Analyzer;
//import ij.plugin.filter.PlugInFilter;
import ij.process.*;
import ij.gui.*;

import java.awt.AWTEvent;
import java.awt.Choice;
import java.awt.Font;
import java.awt.TextField;
import java.awt.Checkbox;
import java.awt.Button;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.image.*;

import cern.colt.list.tdouble.DoubleArrayList;
import cern.colt.list.tint.IntArrayList;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleSingularValueDecomposition;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import cern.jet.math.tdouble.DoubleFunctions;
import cern.jet.math.tdouble.DoublePlusMultSecond;

//public class TestPlugin_ implements PlugInFilter{
public class BaSiC_Plugin implements PlugIn, DialogListener{
	
	private ImageStack stack;
	private int noOfSlices;
	
	@SuppressWarnings("unused")
	private int noOfChannels;
	
	
	private class Shading{
		public ImagePlus flatfield;
		public ImagePlus darkfield;
		//public double[] basefluor;
		//public double[] ratioflat;
		}
	private class Baseline{
		public double[] basefluor;
		//public double[] ratioflat;
	}
	
	private class DecomposedMatrix{
		public DoubleMatrix2D LowrankComponent;
		public DoubleMatrix2D SparseComponent;
		public DoubleMatrix1D Offset;
		public DoubleMatrix1D Coeff;
	}

	private static final class Parameters{
		private static final double epslon = 0.1;
		//private static final double myOptions.lambda = 0.5;
		private static final int reweightingIteration = 5;
		private static final int processingWidth = 128;
		private static final int processingHeight = 128;
		private static final double tolerance = 1e-6;
		private static final double maxIter = 500;
	}
	private class Options{
		// setting default value
		boolean lambda_auto = true;
		boolean lambdadark_auto = true;
		double lambda = 0.5;
		double lambda_dark = 0.5;
		boolean shadingEst = true;
		boolean darkfieldEst = false;
		boolean imageCorr = false;
		int driftOpt = 0;
	}
	
	private static String[] shadingEstimationOptions = {"Skip estimation and use predefined shading profiles","Estimate shading profiles"};
	private static String [] shadingModelOptions = {"Estimate flat-field only (ignore dark-field)", "Estimate both flat-field and dark-field"};
	private static String [] parameterSettingOptions = {"Automatic","Manual"};
	private static String [] driftOptions = {"Ignore","Replace with zero","Replace with temporal mean"};
	private static String [] correctionOptions = {"Compute shading and correct images","Compute shading only"};
	private static String none = "None";
//	private static int maxChannels = 3;
//	private boolean autoFillDisabled;
//	private String firstChannelName;
//	private static String[] colors = {"Input_sequence","Precomputed_flatfield","Precomputed_darkfield"};

	public void run(String arg) {
		showDialog();
	}
	
	void showDialog(){
		int[] wList = WindowManager.getIDList();
        if (wList==null) {
            IJ.noImage();
            return;
        }
        String[] titles = new String[wList.length+1];
        //titles[0] = ""; //default no measured darkfield;
        for (int i=0; i<wList.length; i++) {
            ImagePlus imp_temp = WindowManager.getImage(wList[i]);
            if (imp_temp!=null)
                titles[i] = imp_temp.getTitle();
            else
                titles[i] = "";    
        titles[wList.length] = none;     
        }
        //String[] names = getInitialNames(titles);
       // boolean macro = IJ.macroRunning();
		//GUI
		GenericDialog gd = new GenericDialog("BaSiC",IJ.getInstance());
		gd.addMessage("BaSiC: " + "A Tool for Background and Shading Correction of Optical Microscopy Images",new Font( Font.SANS_SERIF, Font.BOLD, 13 ));
		gd.addMessage("");
		gd.addChoice("Processing_stack", titles, titles[0]);
		gd.addMessage("");
		gd.addRadioButtonGroup("Shading_estimation", shadingEstimationOptions, 1, 2, shadingEstimationOptions[1]);
		gd.addChoice("Flat-field", titles, none);
		gd.addChoice("Dark-field", titles, none);
		gd.addMessage("");
    	gd.addRadioButtonGroup("Shading_model:", shadingModelOptions, 1, 2, shadingModelOptions[0]);
    	gd.addMessage(""); 	
    	gd.addRadioButtonGroup("Setting_regularisationparametes:", parameterSettingOptions, 1, 2, parameterSettingOptions[0]);
    	gd.addSlider("lambda_flat",0.0,5.0,0.5);
    	gd.addSlider("lambda_dark",0.0,5.0,0.5);
    	gd.addMessage("");   	
    	gd.addRadioButtonGroup("Temporal_drift of baseline (ignore or remove & replace):", driftOptions,1,3,driftOptions[0]);
    	gd.addMessage("");
    	gd.addRadioButtonGroup("Correction_options:", correctionOptions, 1, 2, correctionOptions[0]);
    	//gd.addCheckbox("Image correction (untick this option will only compute shading profiles, but not correct images)", true);
		gd.addMessage("");
		gd.addMessage("Copyright \u00a9 2016 Tingying Peng, Helmholtz Zentrum München and TUM, Germany. All rights reserved.");
		final TextField flat_lambda = (TextField) gd.getNumericFields().get(0);
        final TextField dark_lambda = (TextField) gd.getNumericFields().get(1);
        dark_lambda.setEnabled(false);
        flat_lambda.setEnabled(false);
		// add dialog listener
		gd.addDialogListener(this);
		gd.showDialog();
		if (gd.wasCanceled()) return;
		processDialog(gd,wList);
	}
		public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
		{
			if (gd.wasCanceled()) return false;
			int stackindex = gd.getNextChoiceIndex();
			int flatfieldindex = gd.getNextChoiceIndex();
			int darkfieldindex = gd.getNextChoiceIndex();
			String myShadingEstimationChoice = gd.getNextRadioButton();
			String myShadingModelChoice = gd.getNextRadioButton();
			String myParameterChoice = gd.getNextRadioButton();
			String myDriftChoice = gd.getNextRadioButton();
			String myCorrectionChoice = gd.getNextRadioButton();
			double lambda_flat = gd.getNextNumber();
			double lambda_dark = gd.getNextNumber();
	        final TextField flat_lambda = (TextField) gd.getNumericFields().get(0);
	        final TextField dark_lambda = (TextField) gd.getNumericFields().get(1);


	        if ((myShadingEstimationChoice.equals(shadingEstimationOptions[0]))||(myShadingModelChoice.equals(shadingModelOptions[0]))||(myParameterChoice.equals(parameterSettingOptions[0])))
				dark_lambda.setEnabled(false);
			else
				dark_lambda.setEnabled(true);
			
		   if ((myShadingEstimationChoice.equals(shadingEstimationOptions[0]))||(myParameterChoice.equals(parameterSettingOptions[0])))
			   flat_lambda.setEnabled(false);
		   else
			   flat_lambda.setEnabled(true);
		   
		   /*if(Recorder.record){
				 // if (myOptions.shadingEst)
				  {
					  Recorder.record("BaSiC","ShadingEstimationOptions=["+myShadingEstimationChoice+"] ShadingModelOptions=["+myShadingModelChoice+"] parameterSettingOptions=["+myParameterChoice+"] lambda_flat="+lambda_flat+" lambda_dark="+lambda_dark+" BaselineOptions=["+myDriftChoice+"] CorrectionOptions=["+myCorrectionChoice+"]");
				  }
				  Recorder.setCommand(null);
				  
			  }*/
		   return true;
		}
		
        
	   void processDialog(GenericDialog gd, int[] wList)
	   {
	    // get settings from gd
	    int StackIndex = gd.getNextChoiceIndex();
	    int FlatIndex = gd.getNextChoiceIndex();
	    int DarkIndex = gd.getNextChoiceIndex();
	    ImagePlus imp = null;
	    ImagePlus imp_flat = null;
	    ImagePlus imp_dark = null;
	    
	    if (StackIndex==wList.length){
			IJ.error("Please select an image sequence as the input");
			return;
	    }
	    else
	    {	    	
	    	imp = WindowManager.getImage(wList[StackIndex]);
			if (imp.getBitDepth() == 24)
			{
				IJ.error("Please decompose RGB images into single channels.");
				return;
			} 	
			noOfSlices = imp.getNSlices();
			if (noOfSlices == 1){
				IJ.error("Input must be an image stack, not a single image");
				return;
			}
			if (FlatIndex!=wList.length) 
			{
		        imp_flat = WindowManager.getImage(wList[FlatIndex]);
		        if (imp_flat.getNSlices() !=1)
		        {
		        	IJ.error("Flat-field input must be a single image.");
					return;
		        }
			}
			if (DarkIndex!=wList.length) 
			{
		        imp_dark = WindowManager.getImage(wList[DarkIndex]);
		        if (imp_dark.getNSlices() !=1)
		        {
		        	IJ.error("Dark-field input must be a single image.");
					return;
		        }
			}
			
			
			String myShadingEstimationChoice = gd.getNextRadioButton();
		    String myShadingModelChoice = gd.getNextRadioButton();
		    String myParameterChoice = gd.getNextRadioButton();
		    String myDriftChoice = gd.getNextRadioButton();
		    String myCorrectionChoice = gd.getNextRadioButton();
		    double lambda_flat = gd.getNextNumber();
	    	double lambda_dark = gd.getNextNumber();
        
			//* execute the command
			exec(imp,imp_flat,imp_dark, myShadingEstimationChoice, myShadingModelChoice,myParameterChoice,lambda_flat,lambda_dark, myDriftChoice,myCorrectionChoice);		
		   }
	   }

/*public void run(String params, ImagePlus imp)
{
	if (imp==null){
		IJ.error("Please select an image sequence as the input");
		return;
	}
	ImagePlus imp_flat = null;
	ImagePlus imp_dark = null;
	
	exec(imp,imp_flat,imp_dark, myShadingEstimationChoice, myShadingModelChoice,myParameterChoice,lambda_flat,lambda_dark, myDriftChoice,myCorrectionChoice);		
}*/
	
public void exec( ImagePlus imp,
				  ImagePlus imp_flat,
				  ImagePlus imp_dark,
				  String myShadingEstimationChoice,
				  String myShadingModelChoice,
				  String myParameterChoice,
				  double lambda_flat,
				  double lambda_dark,
				  String myDriftChoice,
				  String myCorrectionChoice)
{
	    stack = imp.getStack();
		int outputWidth = stack.getWidth(); 
		int outputHeight = stack.getHeight();
		String title = imp.getTitle();
		Options myOptions = new Options();
		if (myShadingEstimationChoice.equals(shadingEstimationOptions[1]))
			myOptions.shadingEst = true;
		else if (myShadingEstimationChoice.equals(shadingEstimationOptions[0]))
			myOptions.shadingEst = false;
		else{
			IJ.showMessage("Your inputs for Shading_estimation does not fit any of the available options; use default setting instead");
			myOptions.shadingEst = true;
			}
		
		if (myShadingModelChoice.equals(shadingModelOptions[0]))
			myOptions.darkfieldEst = false;
		else if (myShadingModelChoice.equals(shadingModelOptions[1]))
			myOptions.darkfieldEst = true;
		else{
			myOptions.darkfieldEst = false;
			IJ.showMessage("Your inputs for Shading_model does not fit any of the available options; use default setting instead");
			}
			
		
		
		if(myParameterChoice.equals(parameterSettingOptions[0])){
			myOptions.lambdadark_auto = true;
			myOptions.lambda_auto = true;
			}
		else if (myParameterChoice.equals(parameterSettingOptions[1])){
			myOptions.lambdadark_auto = false;
			myOptions.lambda_auto = false;
		    myOptions.lambda_dark = lambda_dark;
		    myOptions.lambda = lambda_flat;
			}
		else{
			myOptions.lambdadark_auto = true;
			myOptions.lambda_auto = true;
			IJ.showMessage("Your inputs for Setting_regularisationparameters does not fit any of the available options; use default setting instead");
			}
			
		
		
		if (myDriftChoice.equals(driftOptions[0]))
			myOptions.driftOpt=0;
		else if(myDriftChoice.equals(driftOptions[1]))
			myOptions.driftOpt=1;
		else if(myDriftChoice.equals(driftOptions[2]))
			myOptions.driftOpt = 2;
		else{
			myOptions.driftOpt = 0;
			IJ.showMessage("Your inputs for Temporal_drift does not fit any of the available options; use default setting instead");
		}
		
		
		if (myCorrectionChoice.equals(correctionOptions[0]))
			myOptions.imageCorr = true;
		else if (myCorrectionChoice.equals(correctionOptions[1]))
			myOptions.imageCorr = false;
		else{
			myOptions.imageCorr = true;
			IJ.showMessage("Your inputs for Correction_options does not fit any of the available options; use default setting instead");
		}
		
		
		 //*  Iterate through each of the slices of the stack to get the respective pixels array from each slice
		 	
	    // Create a new stack object to store resized stack and for processing on each iteration (D)
	    ImageStack processingStack = new ImageStack(Parameters.processingWidth, Parameters.processingHeight);
	    //IJ.log("processingStackDimension" +processingStackDimension);	
	    for(int j=1; j <= noOfSlices; j++)
	    {
			ImageProcessor imageResized = stack.getProcessor(j).convertToFloat();
//			if (!myOptions.darkfieldEst){
//				float[] imageResizedPixels = (float[]) imageResized.getPixels();
//				for (int i = 0; i<imageResized.getWidth()*imageResized.getHeight();i++){
//				    imageResizedPixels[i] = imageResizedPixels[i]-darkfield_preset.getf(i);
//				}
//			}
			imageResized.setInterpolate(true);
			imageResized.setInterpolationMethod(ImageProcessor.BILINEAR);
			ImageProcessor temp = imageResized.resize(Parameters.processingWidth, Parameters.processingWidth,true);
			processingStack.addSlice(temp);
		}
	    
	    Shading processingShading = new Shading();
	    ImageProcessor ip_flatfield = null;
	    ImageProcessor ip_darkfield = null;
	    if (myOptions.shadingEst==false)
	    {
		    if (imp_flat!=null)
		    {
		    	ImageProcessor imageResized = imp_flat.getProcessor().convertToFloat();
		    	imageResized.setInterpolate(true);
				imageResized.setInterpolationMethod(ImageProcessor.BILINEAR);
				DoubleMatrix2D flatfieldMatrix = imageToMatrix( imageResized.resize(Parameters.processingWidth, Parameters.processingWidth,true));
				flatfieldMatrix.normalize();
				flatfieldMatrix.assign(DoubleFunctions.mult((double) Parameters.processingWidth*Parameters.processingHeight));
				processingShading.flatfield = new ImagePlus("flatfield",matrixToImage(flatfieldMatrix));
	    		ip_flatfield = processingShading.flatfield.getProcessor();
		    }
		    else
		    {
		    	 processingShading.flatfield = NewImage.createFloatImage("flatfield", Parameters.processingWidth, Parameters.processingWidth,1, NewImage.GRAY32);
		    	 ip_flatfield = processingShading.flatfield.getProcessor();
		 	     ip_flatfield.set(1.0);
		    }
		    if (imp_dark!=null)
		    {
		    	ImageProcessor imageResized = imp_dark.getProcessor().convertToFloat();
		    	imageResized.setInterpolate(true);
				imageResized.setInterpolationMethod(ImageProcessor.BILINEAR);
				processingShading.darkfield = new ImagePlus("darkfield", imageResized.resize(Parameters.processingWidth, Parameters.processingWidth,true));	 
		    	ip_darkfield = processingShading.darkfield.getProcessor();
		    }
		    else
		    {
		    	processingShading.darkfield = NewImage.createFloatImage("darkfield", Parameters.processingWidth, Parameters.processingWidth,1, NewImage.GRAY32);
		    	ip_darkfield = processingShading.darkfield.getProcessor();
		 	    ip_darkfield.set(0.0);
		    }	    	
	    }
	    	
	    else
	    {   
	        // correction starts: firstly shading and then baseline  
		    processingShading = ShadingCorrection(processingStack,myOptions);
		    ip_flatfield = processingShading.flatfield.getProcessor();
		    ip_darkfield = processingShading.darkfield.getProcessor();
	    }
	    
	    ip_flatfield.setInterpolate(true);
	    ip_flatfield.setInterpolationMethod(ImageProcessor.BICUBIC);
	    ImageProcessor ip_outputflatfield = ip_flatfield.resize(outputWidth,outputHeight);
	    ip_darkfield.setInterpolate(true);
	    ip_darkfield.setInterpolationMethod(ImageProcessor.BICUBIC);
	    ImageProcessor ip_outputdarkfield = ip_darkfield.resize(outputWidth,outputHeight);
	    Baseline processingBaseline = BaselineCorrection(processingStack,myOptions,ip_flatfield,ip_darkfield);
	    float [] basefluor = new float[noOfSlices];
	    float meanbasefluor = 0;
	    for(int i=0; i<noOfSlices; i++)
	    {
	    	basefluor[i]=(float) processingBaseline.basefluor[i];
	    	meanbasefluor = meanbasefluor+basefluor[i];	
	    }
	    meanbasefluor = meanbasefluor/noOfSlices;
	    
	    Shading outputShading = new Shading();
	    
	    outputShading.flatfield = new ImagePlus("Flat-field:".concat(title),ip_outputflatfield); 
	    ContrastEnhancer ce = new ContrastEnhancer();
	    ce.stretchHistogram(ip_outputflatfield,0.35);
	    IJ.run("Enhance Contrast","saturated = 0.35");
	    outputShading.flatfield.show();
	    if (myOptions.darkfieldEst)
	    {
	    	outputShading.darkfield = new ImagePlus("Dark-field:".concat(title),ip_outputdarkfield);
	    	ce.stretchHistogram(ip_outputdarkfield, 0.01);
	    	IJ.run("Enhance Contrast","saturated = 0.01");
	    	outputShading.darkfield.show();
	    }
	    if(myOptions.driftOpt!=0)
	    {
	    	IJ.log("Temporal components:");
	    	float [] x = new float[noOfSlices];	
	    	ResultsTable temporal_rt = new ResultsTable();
	    	Analyzer.setResultsTable(temporal_rt);
	  
	    	for (int j = 1; j <= noOfSlices; j++) {
	    		x[j-1] = j;
	    	    temporal_rt.incrementCounter();
	    	    temporal_rt.addValue("FrameNo", j);
	    	    temporal_rt.addValue("Basefluor", basefluor[j-1]);
//	    	    if (myOptions.varyingflatfieldOpt)
//	    	    	temporal_rt.addValue("Ratioflat",ratioflat[j-1]);
	    	}
	    	temporal_rt.show("Temporal components");
	    	Plot basefluor_plot = new Plot("Basefluor","FrameNo","Basefluor",x,basefluor);
	    	basefluor_plot.show();
//	    	if (myOptions.varyingflatfieldOpt)
//	    	{
//		    	Plot ratioflat_plot = new Plot("Ratioflat","FrameNo","Ratioflat",x,ratioflat);
//		    	ratioflat_plot.show();
//	    	}
	    }
	    // image correction
	    	    
	    if (myOptions.imageCorr)
	    {
		    ImageStack corrected_stack = new ImageStack(imp.getWidth(), imp.getHeight());
		    //IJ.log("processingStackDimension" +processingStackDimension);
		      
		    for(int j=1; j <= noOfSlices; j++)
		    {
		    	ImageProcessor imageCorrected = imp.getStack().getProcessor(j).duplicate().convertToFloat();
				
				float[] imageCorrectedPixels = (float[]) imageCorrected.getPixels();
				if(myOptions.darkfieldEst)
				{for (int i = 0; i<imageCorrected.getWidth()*imageCorrected.getHeight();i++){
					imageCorrectedPixels[i] = (float) ((imageCorrectedPixels[i]-outputShading.darkfield.getProcessor().getf(i))/(outputShading.flatfield.getProcessor().getf(i)));	
					}
				}
				else{
				for (int i = 0; i<imageCorrected.getWidth()*imageCorrected.getHeight();i++){
					imageCorrectedPixels[i] = (float) (imageCorrectedPixels[i]/(outputShading.flatfield.getProcessor().getf(i)));
					}
				}
				if(myOptions.driftOpt!=0)
				{
					//IJ.log("Slice "+j+":basefluor"+basefluor[j-1]+":ratioflat"+ratioflat[j-1]);
					
					
					if(myOptions.driftOpt==1)
					{
					for (int i = 0; i<imageCorrected.getWidth()*imageCorrected.getHeight();i++){
						imageCorrectedPixels[i] = (float) (imageCorrectedPixels[i]-basefluor[j-1]);
						}
					}
					else if(myOptions.driftOpt==2)
					{
						for (int i = 0; i<imageCorrected.getWidth()*imageCorrected.getHeight();i++){
							imageCorrectedPixels[i] = (float) (imageCorrectedPixels[i]-basefluor[j-1]+meanbasefluor);
							}
					}
				}
				// covert to the same bits as the input image
		    	if (imp.getBitDepth()==8)
		    		imageCorrected = imageCorrected.convertToByte(false);
		    	else if (imp.getBitDepth()==16)
		    		imageCorrected = imageCorrected.convertToShort(false);
//		    	else if (imp.getBitDepth()==32)
//		    		imageCorrected = imageCorrected.convertToFloat();
				corrected_stack.addSlice(imp.getStack().getShortSliceLabel(j), imageCorrected);
			}
		    String str1 = "Corrected:";
		    ImagePlus corrected_imp = new ImagePlus(str1.concat(title),corrected_stack);
	    	ce.stretchHistogram(corrected_imp.getProcessor(),0.01);
		    IJ.run("Enhance Contrast","saturated = 0.01");
	    	corrected_imp.show();
	    }
	    
	  
		    
	}
	
	
	private final Shading ShadingCorrection(ImageStack stack,Options myOptions)
	{
		final int stackWidth = stack.getWidth();
		final int stackHeight = stack.getHeight();
		final int stackDim = stackWidth*stackHeight;
		final int stackZSlice = stack.getSize();
//		ImagePlus DStackImage = new ImagePlus("D",stack);
//		DStackImage.show();
		DoubleMatrix2D D =  stackToMatrix(stack);
//		IJ.log("mean of D:"+matrixMean(D.viewDice()));
//		ImageStack DStack = matrixToStack(D,stackWidth,stackHeight);
//		IJ.log("Dimension of DStack:"+DStack.getSize());
//		ImagePlus DStackImage = new ImagePlus("D",DStack);
//		DStackImage.show();
		//setting lambda and lambda dark automatically if automatic parameter setting is used
		if (myOptions.lambda_auto||myOptions.lambdadark_auto)
		{
			DoubleMatrix1D meanD = matrixMean(D);
			meanD.assign(DoubleFunctions.div(meanD.zSum()/meanD.size()));
			DenseDoubleMatrix2D W_meanD = new DenseDoubleMatrix2D(Parameters.processingWidth,Parameters.processingHeight);
			W_meanD.assign(meanD.reshape(Parameters.processingWidth, Parameters.processingHeight));
			W_meanD.dct2(true);
			Double abs_W = W_meanD.assign(DoubleFunctions.abs).zSum();
			if (myOptions.lambda_auto)
				myOptions.lambda = abs_W/(400.0)*0.5;
			if (myOptions.lambdadark_auto)
				myOptions.lambda_dark = abs_W/(400.0)*0.2;
			IJ.log("Smooth regularisation parameters:\n");
			if (myOptions.darkfieldEst)
				IJ.log("\u03BB_flat = " + myOptions.lambda + "; \u03BB_dark = " + myOptions.lambda_dark);
			else
				IJ.log("\u03BB_flat = " + myOptions.lambda);
			
		}
		
		IJ.log("Compute spatial shading profiles...");
		DoubleMatrix2D Dsorted;
	    Dsorted = new DenseDoubleMatrix2D(stackDim,stackZSlice);
		for (int i=0; i<D.rows(); i++)
		{
			Dsorted.viewRow(i).assign(D.viewRow(i).viewSorted()); 
		}
       
       //initiate weight matrix and assign all elements to be one
       DoubleMatrix2D weight = new DenseDoubleMatrix2D(Dsorted.rows(),Dsorted.columns());
       weight.assign(1.0);
		
		DecomposedMatrix D_decompose = new DecomposedMatrix();
	    // estimate flatfield and darkfield by iterative reweighting l1 minimisation
	    for (int iter=1; iter<=Parameters.reweightingIteration; iter++)
	    {
	    	IJ.log("Reweighting Iteration:"+iter);
	    	//L1 minimisation based low rand and sparse decomposition
	    	D_decompose = inexactAlmRspcaL1(Dsorted,weight,myOptions);
//	    	if (iter==2)
//	    		int break_point = 1;
	    	//update the weight
	     	DoubleMatrix1D XARowMean = matrixMean(D_decompose.LowrankComponent.viewDice());
	     	//DoubleMatrix2D XENormalized = new DenseDoubleMatrix2D(stackDim,stackZSlice);
	     	weight.assign(D_decompose.SparseComponent);
	     	for (int u=0; u<weight.columns();u++)
	     	{  
	     		weight.viewColumn(u).assign(DoubleFunctions.div(XARowMean.getQuick(u)+1e-6)); //to avoid zero dividing
	     	}
	     	weight.assign(DoubleFunctions.chain(DoubleFunctions.inv, DoubleFunctions.chain(DoubleFunctions.plus(Parameters.epslon), DoubleFunctions.abs)));
	     	//normalize weight to make sure mean of weights to be one
	     	weight.normalize();
	     	weight.assign(DoubleFunctions.mult((double)weight.size()));
//	     	IJ.log("Mean of the weight matrix:"+weight.zSum()/(double)weight.size());
	    	
	    }
		DoubleMatrix1D temp_XA = matrixMean(D_decompose.LowrankComponent);
		temp_XA.assign(D_decompose.Offset, DoubleFunctions.minus);
		DoubleMatrix2D flatfieldMatrix = new DenseDoubleMatrix2D(Parameters.processingWidth,Parameters.processingHeight);
		flatfieldMatrix.assign(temp_XA.reshape(Parameters.processingWidth,Parameters.processingHeight));
		flatfieldMatrix.normalize();
		flatfieldMatrix.assign(DoubleFunctions.mult((double) Parameters.processingWidth*Parameters.processingHeight));
		//IJ.log("Mean of the flatfield matrix:"+flatfieldMatrix.zSum()/(double)flatfieldMatrix.size());
		DoubleMatrix2D darkfieldMatrix = new DenseDoubleMatrix2D(Parameters.processingWidth,Parameters.processingHeight);
		darkfieldMatrix.assign(D_decompose.Offset.reshape(Parameters.processingWidth, Parameters.processingHeight));
		Shading shadingExp = new Shading();
		shadingExp.darkfield = new ImagePlus("darkfield",matrixToImage(darkfieldMatrix));
		shadingExp.flatfield = new ImagePlus("flatfield",matrixToImage(flatfieldMatrix));
		return shadingExp;
		
		
	}
	
	private final Baseline BaselineCorrection(ImageStack stack,Options myOptions,ImageProcessor flatfield, ImageProcessor darkfield)
	{
		final int stackWidth = stack.getWidth();
		final int stackHeight = stack.getHeight();
		final int stackDim = stackWidth*stackHeight;
		final int stackZSlice = stack.getSize();
//		ImagePlus DStackImage = new ImagePlus("D",stack);
//		DStackImage.show();
		DoubleMatrix2D D =  stackToMatrix(stack);
		DoubleMatrix2D flatfieldMatrix = imageToMatrix(flatfield);
		flatfieldMatrix.normalize();
		flatfieldMatrix.assign(DoubleFunctions.mult((double) Parameters.processingWidth*Parameters.processingHeight));
		DoubleMatrix2D darkfieldMatrix = imageToMatrix(darkfield); 		
		DoubleMatrix1D baseFluor = new DenseDoubleMatrix1D(D.rows());
		DoubleMatrix2D weight = new DenseDoubleMatrix2D(D.rows(),D.columns());
	    weight.assign(1.0);
	    DecomposedMatrix D_decompose = new DecomposedMatrix();
	//	DoubleMatrix1D ratioFlat = new DenseDoubleMatrix1D(D.rows());
		//ratioFlat.assign(1.0);
	    //estimate basefluor and ratioflatamp taking darkfield and flatfield as input
		if(myOptions.driftOpt!=0)
		{
			IJ.log("Compute temporal components...");
			for (int iter=1; iter<=2; iter++)
		    {
		    	IJ.log("Reweighting Iteration:"+iter);
		    	//L1 minimisation based low rand and sparse decomposition
		    	D_decompose = baseFluorEst(D,weight,flatfieldMatrix,darkfieldMatrix,myOptions);
		     	baseFluor = matrixMean(D_decompose.LowrankComponent.viewDice());
		     	baseFluor.assign(D_decompose.Coeff);
		     	weight.assign(D_decompose.SparseComponent);
		     	for (int u=0; u<weight.columns();u++)
		     	{  
		     		weight.viewColumn(u).assign(DoubleFunctions.div(baseFluor.getQuick(u)+1e-8)); //to avoid zero dividing
		     	}
		     	weight.assign(DoubleFunctions.chain(DoubleFunctions.inv, DoubleFunctions.chain(DoubleFunctions.plus(Parameters.epslon), DoubleFunctions.abs)));
		     	//normalize weight to make sure mean of weights to be one
		     	weight.normalize();
		     	weight.assign(DoubleFunctions.mult((double)weight.size()));
//		     	IJ.log("Mean of the weight matrix:"+weight.zSum()/(double)weight.size());
		    	
		    }
		}
				
		Baseline baselineExp = new Baseline();
		baselineExp.basefluor = baseFluor.toArray();
		return baselineExp;		
	}
	
	private final DecomposedMatrix inexactAlmRspcaL1(DoubleMatrix2D D, DoubleMatrix2D weight,Options myOptions)
	{
		//define all internal parameters
		final int ent1 = 1;//for flat field and E
		final int ent2 = 10; //for darkfield
		DenseDoubleSingularValueDecomposition svd = new DenseDoubleSingularValueDecomposition(D,false,false);
		final double norm_two = svd.norm2();
		DenseDoubleAlgebra Algebra = new DenseDoubleAlgebra();
//		ImagePlus DStackImage = new ImagePlus("Dsorted",matrixToStack(D,Parameters.processingWidth,Parameters.processingHeight));
//		DStackImage.show();
		
		//final double norm_two = Algebra.norm2(D);
		double mu = 12.5/norm_two;
		double mu_bar = mu*1e7;
		final double rho = 1.5;
		final double normF_D= Algebra.normF(D);
		
			
		// initiate the low rank and sparse components
		DoubleMatrix2D A_hat = new DenseDoubleMatrix2D(D.rows(),D.columns());
		DoubleMatrix2D E_hat = new DenseDoubleMatrix2D(D.rows(),D.columns());
		DoubleMatrix1D A_offset = new DenseDoubleMatrix1D(D.rows());
		A_offset.assign(0.0);
		DoubleMatrix1D A_coeff = new DenseDoubleMatrix1D(D.columns());
		A_coeff.assign(1.0);
		
		//temporary variables
		double stopCriterion;
		DoubleMatrix2D Y1 = new DenseDoubleMatrix2D(D.rows(),D.columns());
		//DoubleMatrix2D R1 = new DenseDoubleMatrix2D(D.rows(),D.columns());
		DoubleMatrix2D Z1 = new DenseDoubleMatrix2D(D.rows(),D.columns());
		//DoubleMatrix2D temp_Z1 = new DenseDoubleMatrix2D(D.rows(),D.columns());
		//DoubleMatrix2D temp_W = new DenseDoubleMatrix2D(D.rows(),D.columns());
		DenseDoubleMatrix2D temp_Wmean = new DenseDoubleMatrix2D(Parameters.processingHeight,Parameters.processingWidth);
		DoubleMatrix2D W_hat = new DenseDoubleMatrix2D(Parameters.processingWidth,Parameters.processingHeight);
		DenseDoubleMatrix2D W_idct_hat = new DenseDoubleMatrix2D(Parameters.processingWidth,Parameters.processingHeight);
		double[] D_min = D.getMinLocation();
		double B1_uplimit = D_min[0];
		double B1_offset = 0.0;
			
		//variables needed for darkfield estimation
		
//    	IntMatrix2D Mask = new DenseIntMatrix2D(Parameters.processingHeight,Parameters.processingWidth);
//    	Mask.assign(-1);
//		for (int i = (int)Math.round(Mask.rows()/6.0)-1; i<(int)Math.round(Mask.rows()*5.0/6.0); i++)
//			for (int j = (int)Math.round(Mask.columns()/6.0)-1; j<(int)Math.round(Mask.columns()*5.0/6.0); j++)
//				Mask.setQuick(i, j, 1);	
//	    
//		IntArrayList InMaskList = new IntArrayList(); 
//		IntArrayList OutMaskList = new IntArrayList();
//		IntArrayList ValueList = new IntArrayList();
//		Mask.vectorize().getPositiveValues(InMaskList, ValueList);
//		Mask.vectorize().getNegativeValues(OutMaskList, ValueList);
				
		
		
		// main iteration loop starts, the detailed update refer to inexact_alm_rspca_l1.m
		int iter = 0;
		boolean converged = false;
		while (!converged)
		{
			iter = iter + 1;
			W_idct_hat.assign(W_hat);
			W_idct_hat.idct2(true);
			Algebra.multOuter(W_idct_hat.vectorize(), A_coeff, A_hat);
			for (int v=0; v<A_hat.rows();v++)
	     		A_hat.viewRow(v).assign(DoubleFunctions.plus(A_offset.getQuick(v)));
			
		  //  temp_W.assign(computeResidual(D,A_hat,E_hat,Y1,mu,ent1));
		    temp_Wmean.assign(matrixMean(computeResidual(D,A_hat,E_hat,Y1,mu,ent1)).reshape(Parameters.processingWidth, Parameters.processingHeight));
		    temp_Wmean.dct2(true);
		    W_hat.assign(temp_Wmean, DoubleFunctions.plus);
		    W_hat.assign(shrinkageOperator(W_hat,myOptions.lambda/(ent1*mu)));
		    W_idct_hat.assign(W_hat);
		    W_idct_hat.idct2(true);
		    A_hat = Algebra.multOuter(W_idct_hat.vectorize(), A_coeff, A_hat);
		    for (int v=0; v<A_hat.rows();v++)
	     		A_hat.viewRow(v).assign(DoubleFunctions.plus(A_offset.getQuick(v)));
		    E_hat.assign(computeResidual(D,A_hat,E_hat,Y1,mu,ent1), DoubleFunctions.plus);
		   // E_hat.assign(shrinkageOperator(E_hat,1.0/(ent2*mu)));
		    E_hat.assign(shrinkageOperator(E_hat, weight, 1.0/(ent1*mu)));
		    //Z1 is R1 now
		    Z1.assign(D).assign(E_hat,DoubleFunctions.minus);
		    A_coeff.assign(matrixMean(Z1.viewDice()).assign(DoubleFunctions.div((Z1.zSum()/Z1.size()))));
		    A_coeff.assign(DoubleFunctions.max(0.0));
		    
		    // update the offset
		    if(myOptions.darkfieldEst)
		    {
		    	//DoubleMatrix1D B_coeff = new DenseDoubleMatrix1D(D.columns());
		    	DoubleMatrix1D B_offset = new DenseDoubleMatrix1D(D.rows());
		    	DoubleMatrix1D A1_offset = new DenseDoubleMatrix1D(D.rows());
		    	DoubleMatrix1D A_coeffminus1 = new DenseDoubleMatrix1D(D.columns());
		    	IntArrayList validAcoeffList = new IntArrayList();
		    	DoubleArrayList valueList = new DoubleArrayList();
		    	A_coeffminus1.assign(A_coeff).assign(DoubleFunctions.minus(1.0));
		    	A_coeffminus1.getNegativeValues(validAcoeffList,valueList); //A_coeff<1
		    	double mean_W_idct = W_idct_hat.zSum()/W_idct_hat.size();
		    	IntArrayList InMaskList = new IntArrayList(); 
				IntArrayList OutMaskList = new IntArrayList();
				DoubleMatrix1D compare_W_idct_hat = new DenseDoubleMatrix1D(D.rows());
				compare_W_idct_hat.assign(W_idct_hat.vectorize()).assign(DoubleFunctions.minus(mean_W_idct-1e-6));
				compare_W_idct_hat.getPositiveValues(InMaskList, valueList);
				compare_W_idct_hat.assign(W_idct_hat.vectorize()).assign(DoubleFunctions.minus(mean_W_idct+1e-6));
				compare_W_idct_hat.getNegativeValues(OutMaskList, valueList);
				int[] validAcoeff = validElements(validAcoeffList);
				int[] InMask = validElements(InMaskList);
				int[] OutMask = validElements(OutMaskList);
				DoubleMatrix1D B_coeff = new DenseDoubleMatrix1D(validAcoeff.length);  	
		    	B_coeff.assign(matrixMean(Z1.viewSelection(InMask, validAcoeff).viewDice()));
		    	B_coeff.assign(matrixMean(Z1.viewSelection(OutMask, validAcoeff).viewDice()),DoubleFunctions.minus);
		    	B_coeff.assign(DoubleFunctions.div(Z1.zSum()/(double)Z1.size()));
		    	double temp1 = A_coeff.viewSelection(validAcoeff).aggregate(DoubleFunctions.plus,DoubleFunctions.square);
		    	double temp2 = A_coeff.viewSelection(validAcoeff).zSum();
		    	double temp3 = B_coeff.zSum();
		    	double temp4 = A_coeff.viewSelection(validAcoeff).aggregate(B_coeff, DoubleFunctions.plus, DoubleFunctions.mult);
		    	double temp6 = temp2*temp3-validAcoeff.length*temp4;
		    	if (temp6!=0.0)
		    			B1_offset = (temp1*temp3-temp2*temp4)/temp6;
		    	else
		    		B1_offset = 0.0;	    		
		    	// restrict B1_offset value between [0 min(D(:)]
		    	B1_offset = Math.max(B1_offset, 0.0);
		    	B1_offset = Math.min(B1_offset, B1_uplimit/Math.max(mean_W_idct,1e-6));
		    	B_offset.assign(W_idct_hat.vectorize()).assign(DoubleFunctions.mult((-1.0)*B1_offset)).assign(DoubleFunctions.plus(B1_offset*mean_W_idct));
		    	
		    	//A1_offset in matlab
		    	double temp5 = A_coeff.viewSelection(validAcoeff).zSum()/validAcoeffList.size();
		    	A1_offset.assign(W_idct_hat.vectorize()).assign(DoubleFunctions.mult(-1.0*temp5));
		    	A1_offset.assign(matrixMean(Z1.viewSelection(null, validAcoeff)),DoubleFunctions.plus);
//		    	ImagePlus A1_offsetImage = new ImagePlus("A_offset step 1",matrixToImage(A1_offset.reshape(Parameters.processingHeight,Parameters.processingWidth)));
//				A1_offsetImage.show();
		    	// now it is A_offset
		    	A_offset.assign(A1_offset).assign(DoubleFunctions.minus(A1_offset.zSum()/A1_offset.size())).assign(B_offset,DoubleFunctions.minus);
//		    	ImagePlus A2_offsetImage = new ImagePlus("A_offset step 2",matrixToImage(A_offset.reshape(Parameters.processingHeight,Parameters.processingWidth)));
//				A2_offsetImage.show();
		    	// smooth Aoffset
		    	DenseDoubleMatrix2D W_offset = new DenseDoubleMatrix2D(Parameters.processingHeight,Parameters.processingWidth);
		    	W_offset.assign(A_offset.reshape(Parameters.processingHeight,Parameters.processingWidth));
		    	W_offset.dct2(true);
		    	W_offset.assign(shrinkageOperator(W_offset,myOptions.lambda_dark/(ent2*mu)));
		    	W_offset.idct2(true);
//		    	ImagePlus A3_offsetImage = new ImagePlus("A_offset step 3, after smooth",matrixToImage(W_offset));
//				A3_offsetImage.show();
		    	A_offset.assign(W_offset.vectorize());
		    	// encourage sparse A_offset
		    	A_offset.assign(shrinkageOperator(A_offset,myOptions.lambda_dark/(ent2*mu)));
		    	A_offset.assign(B_offset,DoubleFunctions.plus);	
//		    	ImagePlus A4_offsetImage = new ImagePlus("A_offset step 4",matrixToImage(A_offset.reshape(Parameters.processingWidth, Parameters.processingHeight)));
//				A4_offsetImage.show();
		    	//A_offset.assign(B_offset);
		    }
		    // now Z1 is Z1
		    Z1.assign(D).assign(A_hat,DoubleFunctions.minus).assign(E_hat, DoubleFunctions.minus);
		    double normF_Z1 = Algebra.normF(Z1);
		    
		    Y1.assign(Z1,DoublePlusMultSecond.plusMult(mu));
		   
		   // Y1.assign(Z1).assign(DoubleFunctions.mult(mu)), DoubleFunctions.plus);
		    
		    mu = ((mu*rho) < mu_bar) ? mu*rho : mu_bar;  // min(mu*rho, mu_bar)
		    
		    //Stopping criteria
		   // DenseDoubleAlgebra stop_Algebra = new DenseDoubleAlgebra();	    
		    stopCriterion = normF_Z1/normF_D; 
		   // IJ.log("normF of Z1"+normF_Z1+"normF of D"+normF_D);
		    //System.out.println("Stop Criterion" +iteration +"  " +stopCriterion);
			//IJ.log("Stop Criterion" +iter +"  " +stopCriterion); // + "B1_offset" + B1_offset);
			
			if(stopCriterion < Parameters.tolerance)
				converged = true;
			
			if(!converged && iter >= Parameters.maxIter)
			{
				IJ.log("Maximum iteration reached");
				converged = true;
			}
			
		}
							
		DecomposedMatrix D_decomposed = new DecomposedMatrix();
		D_decomposed.LowrankComponent = A_hat;
		D_decomposed.SparseComponent = E_hat;
		D_decomposed.Offset = A_offset;
		D_decomposed.Offset.assign(W_idct_hat.vectorize(),DoublePlusMultSecond.plusMult(B1_offset));
		D_decomposed.Coeff = A_coeff;
		return D_decomposed;
	}
	
	private final DecomposedMatrix baseFluorEst(DoubleMatrix2D D, DoubleMatrix2D weight,DoubleMatrix2D flatfield, DoubleMatrix2D darkfield,Options myOptions)
	{
		//define all internal parameters
		final int ent1 = 1;//for flat field and E
		DenseDoubleSingularValueDecomposition svd = new DenseDoubleSingularValueDecomposition(D,false,false);
		final double norm_two = svd.norm2();
		DenseDoubleAlgebra Algebra = new DenseDoubleAlgebra();
//		ImagePlus DStackImage = new ImagePlus("Dsorted",matrixToStack(D,Parameters.processingWidth,Parameters.processingHeight));
//		DStackImage.show();
		
		//final double norm_two = Algebra.norm2(D);
		double mu = 12.5/norm_two;
		double mu_bar = mu*1e7;
		final double rho = 1.5;
		final double normF_D= Algebra.normF(D);
		
			
		// initiate the low rank and sparse components
		DoubleMatrix2D A_hat = new DenseDoubleMatrix2D(D.rows(),D.columns());
		DoubleMatrix2D E_hat = new DenseDoubleMatrix2D(D.rows(),D.columns());
		DoubleMatrix1D A_offset = darkfield.vectorize();
		DoubleMatrix1D A_coeff = matrixMean(D.viewDice());
		DenseDoubleMatrix2D W_idct_hat = new DenseDoubleMatrix2D(Parameters.processingWidth,Parameters.processingHeight);
		W_idct_hat.assign(flatfield);
		
		//temporary variables
		double stopCriterion;
		DoubleMatrix2D Y1 = new DenseDoubleMatrix2D(D.rows(),D.columns());
		//DoubleMatrix2D R1 = new DenseDoubleMatrix2D(D.rows(),D.columns());
		DoubleMatrix2D Z1 = new DenseDoubleMatrix2D(D.rows(),D.columns());
						
		// main iteration loop starts, the detailed update refer to inexact_alm_rspca_l1.m
		int iter = 0;
		boolean converged = false;
		while (!converged)
		{
			iter = iter + 1;
//			W_idct_hat.assign(W_hat);
//			W_idct_hat.idct2(true);
			Algebra.multOuter(W_idct_hat.vectorize(), A_coeff, A_hat);
			for (int v=0; v<A_hat.rows();v++)
	     		A_hat.viewRow(v).assign(DoubleFunctions.plus(A_offset.getQuick(v)));
		    E_hat.assign(computeResidual(D,A_hat,E_hat,Y1,mu,ent1), DoubleFunctions.plus);
		   // E_hat.assign(shrinkageOperator(E_hat,1.0/(ent2*mu)));
		    E_hat.assign(shrinkageOperator(E_hat, weight, 1.0/(ent1*mu)));
		    //Z1 is R1 now
		    Z1.assign(D).assign(E_hat,DoubleFunctions.minus);
		    A_coeff.assign(matrixMean(Z1.viewDice()));
		    A_coeff.assign(DoubleFunctions.minus(A_offset.zSum()/A_offset.size()));
		    A_coeff.assign(DoubleFunctions.max(0.0));
		    
		    // now Z1 is Z1
		    Z1.assign(D).assign(A_hat,DoubleFunctions.minus).assign(E_hat, DoubleFunctions.minus);
		    double normF_Z1 = Algebra.normF(Z1);
		    
		    Y1.assign(Z1,DoublePlusMultSecond.plusMult(mu));
		   
		   // Y1.assign(Z1).assign(DoubleFunctions.mult(mu)), DoubleFunctions.plus);
		    
		    mu = ((mu*rho) < mu_bar) ? mu*rho : mu_bar;  // min(mu*rho, mu_bar)
		    
		    //Stopping criteria
		   // DenseDoubleAlgebra stop_Algebra = new DenseDoubleAlgebra();	    
		    stopCriterion = normF_Z1/normF_D; 
		   // IJ.log("normF of Z1"+normF_Z1+"normF of D"+normF_D);
		    //System.out.println("Stop Criterion" +iteration +"  " +stopCriterion);
			//IJ.log("Stop Criterion" +iter +"  " +stopCriterion); // + "B1_offset" + B1_offset);
			
			if(stopCriterion < Parameters.tolerance)
				converged = true;
			
			if(!converged && iter >= Parameters.maxIter)
			{
				IJ.log("Maximum iteration reached");
				converged = true;
			}
			
		}							
		DecomposedMatrix D_decomposed = new DecomposedMatrix();
		D_decomposed.LowrankComponent = A_hat;
		D_decomposed.SparseComponent = E_hat;
		D_decomposed.Offset = A_offset;
		D_decomposed.Coeff = A_coeff;
		return D_decomposed;
	}
	
	/*private final DecomposedMatrix ratioFlatEst(DoubleMatrix2D D, DoubleMatrix2D weight,DoubleMatrix2D flatfield, DoubleMatrix2D darkfield,Options myOptions)
	{
		//define all internal parameters
		final int ent1 = 1;//for flat field and E
		DenseDoubleSingularValueDecomposition svd = new DenseDoubleSingularValueDecomposition(D,false,false);
		final double norm_two = svd.norm2();
		DenseDoubleAlgebra Algebra = new DenseDoubleAlgebra();
//		ImagePlus DStackImage = new ImagePlus("Dsorted",matrixToStack(D,Parameters.processingWidth,Parameters.processingHeight));
//		DStackImage.show();
		
		//final double norm_two = Algebra.norm2(D);
		double mu = 12.5/norm_two;
		double mu_bar = mu*1e7;
		final double rho = 1.5;
		final double normF_D= Algebra.normF(D);
		
			
		// initiate the low rank and sparse components
		DoubleMatrix2D A_hat = new DenseDoubleMatrix2D(D.rows(),D.columns());
		DoubleMatrix2D E_hat = new DenseDoubleMatrix2D(D.rows(),D.columns());
		DoubleMatrix1D A_offset = darkfield.vectorize();
		DoubleMatrix1D A_coeff = new DenseDoubleMatrix1D(D.columns());
		A_coeff.assign(0.0);
		DenseDoubleMatrix2D W_idct_hat = new DenseDoubleMatrix2D(Parameters.processingWidth,Parameters.processingHeight);
		W_idct_hat.assign(flatfield).assign(DoubleFunctions.minus(1.0));
		
		//temporary variables
		double stopCriterion;
		DoubleMatrix2D Y1 = new DenseDoubleMatrix2D(D.rows(),D.columns());
		//DoubleMatrix2D R1 = new DenseDoubleMatrix2D(D.rows(),D.columns());
		DoubleMatrix2D Z1 = new DenseDoubleMatrix2D(D.rows(),D.columns());
						
		// main iteration loop starts, the detailed update refer to inexact_alm_rspca_l1.m
		int iter = 0;
		boolean converged = false;
		//int [] select_idx = new int[Parameters.processingHeight*Parameters.processingWidth/10];
//		for (int i = 0; i<select_idx.length; i++)
//			select_idx[i] = i*10;
		while (!converged)
		{
			iter = iter + 1;
//			W_idct_hat.assign(W_hat);
//			W_idct_hat.idct2(true);
			Algebra.multOuter(W_idct_hat.vectorize(), A_coeff, A_hat);
//			for (int v=0; v<A_hat.rows();v++)
//	     		A_hat.viewRow(v).assign(DoubleFunctions.plus(A_offset.getQuick(v)));
		    E_hat.assign(computeResidual(D,A_hat,E_hat,Y1,mu,ent1), DoubleFunctions.plus);
		   // E_hat.assign(shrinkageOperator(E_hat,1.0/(ent2*mu)));
		    E_hat.assign(shrinkageOperator(E_hat, weight, 1.0/(ent1*mu)));
		    //Z1 is R1 now
		    Z1.assign(D).assign(E_hat,DoubleFunctions.minus);
		    for (int v=0; v<Z1.columns();v++)
	     		Z1.viewColumn(v).assign(W_idct_hat.vectorize(),DoubleFunctions.div);     
		    
		    //A_coeff.assign(matrixMedian(Z1.viewSelection(select_idx,null).viewDice()));
		    A_coeff.assign(matrixMedian(Z1.viewDice()));
		    //A_coeff.assign(DoubleFunctions.max(0.0));
		    
		    // now Z1 is Z1
		    Z1.assign(D).assign(A_hat,DoubleFunctions.minus).assign(E_hat, DoubleFunctions.minus);
		    double normF_Z1 = Algebra.normF(Z1);
		    
		    Y1.assign(Z1,DoublePlusMultSecond.plusMult(mu));
		   
		   // Y1.assign(Z1).assign(DoubleFunctions.mult(mu)), DoubleFunctions.plus);
		    
		    mu = ((mu*rho) < mu_bar) ? mu*rho : mu_bar;  // min(mu*rho, mu_bar)
		    
		    //Stopping criteria
		   // DenseDoubleAlgebra stop_Algebra = new DenseDoubleAlgebra();	    
		    stopCriterion = normF_Z1/normF_D; 
		   // IJ.log("normF of Z1"+normF_Z1+"normF of D"+normF_D);
		    //System.out.println("Stop Criterion" +iteration +"  " +stopCriterion);
			IJ.log("Stop Criterion" +iter +"  " +stopCriterion); // + "B1_offset" + B1_offset);
			
			if(stopCriterion < Parameters.tolerance)
				converged = true;
			
			if(!converged && iter >= Parameters.maxIter)
			{
				IJ.log("Maximum iteration reached");
				converged = true;
			}
			
		}							
		DecomposedMatrix D_decomposed = new DecomposedMatrix();
		D_decomposed.LowrankComponent = A_hat;
		D_decomposed.SparseComponent = E_hat;
		D_decomposed.Offset = A_offset;
		D_decomposed.Coeff = A_coeff;
		return D_decomposed;
	}*/
	
	public final int[] validElements(IntArrayList array) {
			    int[] myElements = new int[array.size()];
		        for (int i=0; i<array.size();i++ ) 
		        	myElements[i]=array.getQuick(i);
		        
			     return myElements;
		}
	
	public final DoubleMatrix1D matrixMean(DoubleMatrix2D A)
	// compute the row means of a matrix, use dice view to compute the column means
	{
		DoubleMatrix1D average = new DenseDoubleMatrix1D(A.rows());
		for(int i = 0; i<A.rows(); i++)
			average.setQuick(i,A.viewRow(i).zSum()/A.columns());
		return average;
	}
	
	public final DoubleMatrix1D matrixMedian(DoubleMatrix2D A)
	// compute the row medians of a matrix, use dice view to compute the column means
	{
		DoubleMatrix1D median = new DenseDoubleMatrix1D(A.rows());
		//DoubleMatrix1D temp = new DenseDoubleMatrix1D(A.rows());
		for(int i = 0; i<A.rows(); i++)
			//temp.assign(A.viewRow(i));
			median.setQuick(i,A.viewRow(i).viewSorted().getQuick(A.columns()/2));
		return median;
	}
	
	private final DoubleMatrix2D computeResidual(DoubleMatrix2D D, DoubleMatrix2D A_hat, DoubleMatrix2D E_hat, DoubleMatrix2D Y1, double mu, double ent)
	//Residual = (D-A_hat-E_hat+(1/mu)*Y1)/ent
	{
		DoubleMatrix2D Residual = new DenseDoubleMatrix2D(D.rows(),D.columns());
		Residual.assign(Y1).assign(DoubleFunctions.div(mu)).assign(D,DoubleFunctions.plus).assign(A_hat,DoubleFunctions.minus).assign(E_hat,DoubleFunctions.minus).assign(DoubleFunctions.div(ent));
		return Residual;
	}
	
	private final DoubleMatrix2D shrinkageOperator(DoubleMatrix2D a_offset, double epsilon)
	//max(A-epsilon,0)+min(A+epsilon,0);
	{
		DoubleMatrix2D A_shrink = new DenseDoubleMatrix2D(a_offset.rows(),a_offset.columns());
		//DoubleMatrix2D A_shrinkpositive = new DenseDoubleMatrix2D(A.rows(),A.columns());
		DoubleMatrix2D A_shrinknegative = new DenseDoubleMatrix2D(a_offset.rows(),a_offset.columns());
		A_shrink.assign(a_offset).assign(DoubleFunctions.minus(epsilon)).assign(DoubleFunctions.max(0));
		A_shrinknegative.assign(a_offset).assign(DoubleFunctions.plus(epsilon)).assign(DoubleFunctions.min(0));
		A_shrink.assign(A_shrinknegative,DoubleFunctions.plus);
		return A_shrink;
	}
	
	private final DoubleMatrix1D shrinkageOperator(DoubleMatrix1D a_offset, double epsilon)
	//overload method max(A-epsilon,0)+min(A+epsilon,0);
	{
		DoubleMatrix1D A_shrink = new DenseDoubleMatrix1D((int) a_offset.size());
		//DoubleMatrix2D A_shrinkpositive = new DenseDoubleMatrix2D(A.rows(),A.columns());
		DoubleMatrix1D A_shrinknegative = new DenseDoubleMatrix1D((int) a_offset.size());
		A_shrink.assign(a_offset).assign(DoubleFunctions.minus(epsilon)).assign(DoubleFunctions.max(0));
		A_shrinknegative.assign(a_offset).assign(DoubleFunctions.plus(epsilon)).assign(DoubleFunctions.min(0));
		A_shrink.assign(A_shrinknegative,DoubleFunctions.plus);
		return A_shrink;
	}
	
	private final DoubleMatrix2D shrinkageOperator(DoubleMatrix2D A, DoubleMatrix2D weight, double epsilon)
	//overload method
	//max(A-weight.*epsilon,0)+min(A+weight.*epsilon,0);
	{
		DoubleMatrix2D A_shrink = new DenseDoubleMatrix2D(A.rows(),A.columns());
		//DoubleMatrix2D A_shrinkpositive = new DenseDoubleMatrix2D(A.rows(),A.columns());
		DoubleMatrix2D A_shrinknegative = new DenseDoubleMatrix2D(A.rows(),A.columns());
		A_shrink.assign(weight).assign(DoubleFunctions.mult((-1.0)*epsilon)).assign(A,DoubleFunctions.plus).assign(DoubleFunctions.max(0));
		A_shrinknegative.assign(weight).assign(DoubleFunctions.mult(epsilon)).assign(A,DoubleFunctions.plus).assign(DoubleFunctions.min(0));
		A_shrink.assign(A_shrinknegative,DoubleFunctions.plus);
		return A_shrink;
	}
	
	
	public final DoubleMatrix2D stackToMatrix(ImageStack stack)
	{	
		final int rows = stack.getHeight();
		final int columns = stack.getWidth();
		final int nSlices = stack.getSize(); 
		DoubleMatrix2D stackMatrix = new DenseDoubleMatrix2D(rows*columns,nSlices);
		
		for(int k = 1; k<=nSlices; k++)
		{
			for (int i = 0; i<stack.getWidth(); i++)
			{
				for (int j = 0; j<stack.getHeight(); j++)
				{
					double temp1 = stack.getVoxel(i,j,k-1);
					stackMatrix.setQuick(i*rows+j, k-1, temp1);
				}
			}
				
		}
		
		return stackMatrix;
	}
	
	
	
	public final ImageStack matrixToStack(DoubleMatrix2D matrix,int stackWidth, int stackHeight)
	{
		int stackSize = stackWidth*stackHeight;
		if (matrix.rows()!=stackSize){
			throw new EmptyStackException();
		}
		int stackZSlices = matrix.columns();	
		
		
		ImageStack stack = new ImageStack(stackWidth, stackHeight, stackZSlices);
		ImagePlus imgPlus = NewImage.createFloatImage(null, stackWidth, stackHeight, stackZSlices, NewImage.GRAY32);
		stack = imgPlus.getStack();
		for(int k = 1; k<=stackZSlices; k++)
		{
			ImageProcessor ip = matrixToImage(matrix.viewColumn(k-1).reshape(stackHeight, stackWidth));
			stack.setProcessor(ip,k);		
				
		}
		
		return stack;
	}
	
	public final ImageProcessor matrixToImage(DoubleMatrix2D matrix)
	{
		
		int imageWidth = matrix.columns();
		int imageHeight = matrix.rows();
		//ImageProcessor ip= new ImageProcessor(imageWidth,imageHeight);
		ImagePlus imgPlus = NewImage.createFloatImage(null, imageWidth, imageHeight,1, NewImage.GRAY32);
		ImageProcessor ip = imgPlus.getProcessor();
		
		for (int i = 0; i<imageWidth; i++)
		{
			for (int j = 0; j<imageHeight;j++)
			{
				double temp1 = matrix.getQuick(i,j);
				ip.putPixelValue(j, i, temp1);
			}
		
		}
		return ip;
	}
	
	public final DoubleMatrix2D imageToMatrix(ImageProcessor ip)
	{
		
		int rows = ip.getHeight();
		int columns = ip.getWidth();
		DoubleMatrix2D matrix = new DenseDoubleMatrix2D(rows,columns);
		
		for (int i = 0; i<ip.getWidth(); i++)
		{
			for (int j = 0; j<ip.getHeight();j++)
			{
                		
				//double temp1 = ip.get(i,j);
				matrix.setQuick(j, i, (double)ip.getf(i,j));
			}
		
		}
		//IJ.log("Imagepixel+" ip.getf);
		return matrix;
	}
	
}
			
			
			
