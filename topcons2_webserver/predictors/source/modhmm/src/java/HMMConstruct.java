/* the idea is that if there is only one alphabet, the hmm is saved in the
 * original one-alphabet format, if there are 2-4 alphabets, the hmm
 * is saved in multialphabet format */

class modhmmc
{
    public static void main(String[] args)
    {
	if(args.length == 0) {
	    DataReaderMultiAlpha datareader = new DataReaderMultiAlpha();
	}
	else {
	    System.out.println("Wrong number of arguments");
	}
    }
}
