package de.embl.cba.illuminationcorrection;

import ij.IJ;

public abstract class Log
{

	public static void info( String message )
	{
		IJ.log( message );
	}


}
