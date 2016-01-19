import java.io.*;
import java.util.*;
import java.text.*;

class HMMSaver
{

    HMM theHMM;
    BufferedWriter writer = null;
    
    public HMMSaver(HMM hmm)
    {
	theHMM = hmm;

    }

    public int saveHMM(String outdir)
    {
	String filename = theHMM.getFileName();
	/* if file exists, return EXISTS */
	
	/* start writing */
	try {
	    if(outdir != null && outdir.length() > 0) {
		writer = new BufferedWriter(new FileWriter(outdir + filename));
		System.out.println("Saving " + outdir + filename);
	    }
	    else {
		writer = new BufferedWriter(new FileWriter(filename));
	    }
	    /* write header */
	    writeln("***********************Header*****************************");
	    writeln("NAME: " + theHMM.getName());
	    writeln("TIME OF CREATION: " + getTime());
	    int nrOfAlphabets = theHMM.getNrOfAlphabets();
	    if(nrOfAlphabets == 1) {
		writeln("ALPHABET: "+ theHMM.getAlphabet());
		writeln("ALPHABET LENGTH: " + theHMM.getAlphabetSize());
	    }
	    else {
		writeln("NR OF ALPHABETS: "+ nrOfAlphabets);
		writeln("ALPHABET 1: "+ theHMM.getAlphabet());
		writeln("ALPHABET LENGTH 1: " + theHMM.getAlphabetSize());
		writeln("ALPHABET 2: "+ theHMM.getAlphabet2());
		writeln("ALPHABET LENGTH 2: " + theHMM.getAlphabetSize2());
		if(nrOfAlphabets > 2) {
		    writeln("ALPHABET 3: "+ theHMM.getAlphabet3());
		    writeln("ALPHABET LENGTH 3: " + theHMM.getAlphabetSize3());
		}
		if(nrOfAlphabets > 3) {
		    writeln("ALPHABET 4: "+ theHMM.getAlphabet4());
		    writeln("ALPHABET LENGTH 4: " + theHMM.getAlphabetSize4());
		}
	    }
	    writeln("NR OF MODULES: " + theHMM.getNrModules());
	    writeln("NR OF VERTICES: " + theHMM.getNrVertices());
	    writeln("NR OF TRANSITIONS: " + theHMM.getNrOfTransitions());
	    writeln("NR OF DISTRIBUTION GROUPS: " + theHMM.getNrOfDistributionGroups());
	    //+ "/" + theHMM.getNrOfDistributionGroupVertices());
	    writeln("NR OF TRANSITION TIE GROUPS: " + theHMM.getNrOfTransTieGroups());
	    //+ "/" + theHMM.getNrOfTransTieGroupVertices());
	    if(nrOfAlphabets == 1) {
		writeln("NR OF EMISSION PRIORFILES: " + theHMM.getNrOfPriorfiles());
		write("EMISSION PRIORFILES: ");
		for(Iterator i = theHMM.getPriorfiles(); i.hasNext();) {
		    String s = (String)i.next();
		    write(s + " ");
		}
		writeln("");
	    }
	    else {
		writeln("NR OF EMISSION PRIORFILES 1: " + theHMM.getNrOfPriorfiles());
		write("EMISSION PRIORFILES_1: ");
		for(Iterator i = theHMM.getPriorfiles(); i.hasNext();) {
		    String s = (String)i.next();
		    write(s + " ");
		}
		writeln("");
		writeln("NR OF EMISSION PRIORFILES 2: " + theHMM.getNrOfPriorfiles(2));
		write("EMISSION PRIORFILES_2: ");
		for(Iterator i = theHMM.getPriorfiles(2); i.hasNext();) {
		    String s = (String)i.next();
		    write(s + " ");
		}
		writeln("");
		if(nrOfAlphabets > 2) {
		    writeln("NR OF EMISSION PRIORFILES 3: " + theHMM.getNrOfPriorfiles(3));
		    write("EMISSION PRIORFILES_3: ");
		    for(Iterator i = theHMM.getPriorfiles(3); i.hasNext();) {
			String s = (String)i.next();
			write(s + " ");
		    }
		    writeln("");
		}
		if(nrOfAlphabets > 3) {
		    writeln("NR OF EMISSION PRIORFILES 4: " + theHMM.getNrOfPriorfiles(4));
		    write("EMISSION PRIORFILES_4: ");
		    for(Iterator i = theHMM.getPriorfiles(4); i.hasNext();) {
			String s = (String)i.next();
			write(s + " ");
		    }
		    writeln("");
		}
		
	    }
	    writeln("NR OF TRANSITION PRIORFILES: " + theHMM.getNrOfTransPriorfiles());
	    write("TRANSITION PRIORFILES: ");
	    for(Iterator i = theHMM.getTransPriorfiles(); i.hasNext();) {
		String s = (String)i.next();
		write(s + " ");
	    }
	    writeln("");
	    writeln("");
	    writeln("");
	    

	    /* write the modules giving them the writer for better efficiency */
	    writeln("***********************Modules****************************");

	    for(ListIterator i = theHMM.getModules(); i.hasNext();) {
		 String name =  ((Module)i.next()).getName();
		 if(!theHMM.writeModule(name, writer)) {
		     P.MESSAGE("IO error: could not write HMM to file");
		     System.exit(0);
		 }
	    }

	    /* write the distribution groups */
	    writeln("");
	    writeln("***************Emission distribution groups******************");
	    int q = 1;
	    for(ListIterator i = theHMM.getDistributionGroups();i.hasNext();) {
		write("Group " + q + ": ");
		LinkedList group = (LinkedList)i.next();
		for(ListIterator j = (ListIterator)(group.iterator());j.hasNext();) {
		    Module module = theHMM.getModule((String)(j.next()));
		    LinkedList vertices = module.getVertices();
		    for(ListIterator k = (ListIterator)(vertices.iterator());k.hasNext();) {
			int nr = ((Vertex)(k.next())).getNumber();
			write(nr + " ");
		    }
		}
		writeln("");
		q++;
	    }

	    /* write the transition tie groups */
	    writeln("");
	    writeln("***************Transition tie groups******************");
	    int q2 = 1;
	    for(ListIterator i = theHMM.getTransTieGroups();i.hasNext();) {
		LinkedList group = (LinkedList)i.next();
		LinkedList[] moduleVertices = new LinkedList[group.size()];
		int k = 0;
		for(ListIterator j = (ListIterator)(group.iterator());j.hasNext();) {
		    moduleVertices[k] = theHMM.getModule((String)(j.next())).getVertices();
		    k++;
		}
		int lLimit = getLimit(moduleVertices);


		for(int l = 0; l < lLimit; l++) {
		    LinkedList endTransVertices = new LinkedList();
		    LinkedList noEndTransVertices = new LinkedList();
		    for(k = 0; k < moduleVertices.length; k++) {
			if(((Vertex)moduleVertices[k].get(l)).getNrOfEndTransitions() == 0) {
			    noEndTransVertices.add(((Vertex)moduleVertices[k].get(l)));
			}
			else if(((Vertex)moduleVertices[k].get(l)).getNrOfEndTransitions() == 1) {
			    endTransVertices.add(((Vertex)moduleVertices[k].get(l)));
			}
			else {
			    P.INTERNAL_ERROR("HMMSaver: More than one end transition");
			}
		    }
		    
		    int mLimit = getLimit(noEndTransVertices);
		    for(int m = 0; m < mLimit; m++) {
			//if(noEndTransVertices.size() > 1) {
			if(!differentOutDegree(noEndTransVertices)) {
			    write("Tie " + q2 + ": ");
			    for(k = 0; k < noEndTransVertices.size(); k++) {
				int from = 0;
				from = ((Vertex)noEndTransVertices.get(k)).getNumber();
				Transition toT = ((Transition)((Vertex)noEndTransVertices.get(k)).transitionList.get(m));
				int to = toT.getToVertex();
				write(from + "->" + to + " ");
				writer.flush();
			    }
			    writeln("");
			    q2++;
			}
		    }
		    
		    mLimit = getLimit(endTransVertices);
		    for(int m = 0; m < mLimit; m++) {
			//if(noEndTransVertices.size() > 1) {
			if(!differentOutDegree(endTransVertices)) {
			    write("Tie " + q2 + ": ");
			    for(k = 0; k < endTransVertices.size(); k++) {
				int from = 0;
				from = ((Vertex)endTransVertices.get(k)).getNumber();
				Transition toT = ((Transition)((Vertex)endTransVertices.get(k)).transitionList.get(m));
				int to = toT.getToVertex();
			    write(from + "->" + to + " ");
			    writer.flush();
			    }
			    writeln("");
			    q2++;
			}
		    }
		}
			
	    }
	    //((Vertex)moduleVertices[0].get(l)).transitionList.size()



	    writeln("");
	    
	    writer.flush();
	    
	}
	catch(IOException e) {
	    /* check what went wrong and return different values */
	    System.out.println("Could not save hmm");
	}

	return 0;
    }

    private void write(String s)
    {
	try {
	    writer.write(s);
	}
	catch(IOException e) {
	    P.INTERNAL_ERROR("HMMSaver.write: Could not write to buffer");
	    }
    }

    private void writeln(String s)
    {
	try {
	    writer.write(s);
	    writer.newLine();
	}
	catch(IOException e) {
	    P.INTERNAL_ERROR("HMMSaver.writeln: Could not write to buffer");
	    }	
    }
    
    private String getTime()
    {
	DateFormat df = DateFormat.getDateTimeInstance(DateFormat.MEDIUM, DateFormat.MEDIUM);
	df.setTimeZone(TimeZone.getTimeZone("ECT"));
	return df.format(new Date());
    }
 

    
    private int getLimit(LinkedList[] moduleVertices)
    {
	int limit = 0;
	for(int i = 0; i < moduleVertices.length; i++) {
	    if(moduleVertices[i].size() > limit) {
		limit = moduleVertices[i].size();
	    }
	}
	return limit;
    }

    
    private int getLimit(LinkedList vertices)
    {
	int limit = 0;
	for(int i = 0; i < vertices.size(); i++) {
	    if(((Vertex)vertices.get(i)).transitionList.size() > limit) {
		limit = ((Vertex)vertices.get(i)).transitionList.size();
	    }
	}
	return limit;
    }

    private boolean differentOutDegree(LinkedList vertices)
    {
	int outDegree = 0;
	boolean res = false;
	if(vertices.size() > 0) {
	    outDegree = ((Vertex)vertices.get(0)).getNrOfTransitions();
	}
	for(int i = 1; i < vertices.size(); i++) {
	    if(((Vertex)vertices.get(i)).getNrOfTransitions() != outDegree) {
		res = true;
	    }
	}
	return res;
    }
}
