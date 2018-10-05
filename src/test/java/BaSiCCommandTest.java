import de.embl.cba.illuminationcorrection.basic.BaSiCCommand;
import net.imagej.ImageJ;

public class BaSiCCommandTest
{

	public static void main(final String... args) throws Exception
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		// invoke the plugin
		ij.command().run( BaSiCCommand.class, true );
	}

}
