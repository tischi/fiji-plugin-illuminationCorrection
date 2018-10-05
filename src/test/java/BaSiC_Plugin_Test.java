import de.embl.cba.illuminationcorrection.originalsources.basic.BaSiC_Plugin;
import ij.IJ;
import net.imagej.ImageJ;

public class BaSiC_Plugin_Test
{
	public static void main( String[] args )
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		IJ.open( "/Users/tischer/Documents/fiji-plugin-illuminationCorrection/src/test/resources/transmission-3d.zip" );

		final BaSiC_Plugin baSiC_plugin = new BaSiC_Plugin();
		baSiC_plugin.run( "" );
	}
}
