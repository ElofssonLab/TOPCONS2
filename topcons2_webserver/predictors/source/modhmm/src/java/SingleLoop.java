import java.util.*;
import java.io.*;

class SingleLoop extends Module
{
    String name;
    Vertex theVertex;
    int vertexType;
    int distribType;
    int distribType_2;
    int distribType_3;
    int distribType_4;
    int length;
    String label;
    int size;
    String priorfile;
    String priorfile_2;
    String priorfile_3;
    String priorfile_4;
    String transPriorfile;

    public SingleLoop(String name, int distribType, int vertexType, int length, String label)
    {
	this.name = name;
	this.vertexType = vertexType;
	this.distribType = distribType;
	this.distribType_2 = distribType;
	this.distribType_3 = distribType;
	this.distribType_4 = distribType;
	this.length = length;
	this.label = label;
	theVertex = new Vertex(name, vertexType, distribType, label);
	inVertices = new int[1];
	inVertices[0] = theVertex.getNumber();
	outVertices = new int[1];
	outVertices[0] = theVertex.getNumber();
	priorfile = null;
	priorfile_2 = null;
	priorfile_3 = null;
	priorfile_4 = null;
	transPriorfile = null;
	theVertex.addTransition(theVertex.getNumber(),((double)length)/((double)(length+1)),true);
	size = 1;
    }

    public SingleLoop(String name, double[] initDistrib, int vertexType, int length, String label)
    {
	this.name = name;
	this.vertexType = vertexType;
	this.distribType = HMM.MANUAL;
	this.distribType_2 = HMM.MANUAL;
	this.distribType_3 = HMM.MANUAL;
	this.distribType_4 = HMM.MANUAL;
	this.length = length;
	this.label = label;
	theVertex = new Vertex(name, vertexType, initDistrib, label);
	inVertices = new int[1];
	inVertices[0] = theVertex.getNumber();
	outVertices = new int[1];
	outVertices[0] = theVertex.getNumber();
	priorfile = null;
	priorfile_2 = null;
	priorfile_3 = null;
	priorfile_4 = null;
	transPriorfile = null;
	theVertex.addTransition(theVertex.getNumber(),((double)length)/((double)(length+1)),true);
	size = 1;
    }

    public int getSize()
    {
	return size;
    }

    public String getName()
    {
	return name;
    }

    public String getLabel()
    {
	return label;
    }
    
    public String getPriorfile()
    {
	return priorfile;
    }
    
    public LinkedList getVertices()
    {
	LinkedList l = new LinkedList();
	l.add(theVertex);
	return l;
    }

    public int[] getInVertices()
    {
	return inVertices;
    }

    public int[] getOutVertices()
    {
	return outVertices;
    }

    public int getVertexType()
    {
	return vertexType;
    }
    
    public int getDistribType()
    {
	return distribType;
    }
    
    public Vertex getVertex(int nr)
    {
	Vertex v = theVertex;
	if(v.getNumber() == nr) {
	    return v;
	}

	return null;
    }

    public double[] getEmissionProbs()
    {
	return theVertex.getEmissionProbs();
    }
    
    public int getNrOfTransitions() {
	return theVertex.getNrOfTransitions() + theVertex.getNrOfEndTransitions();
    }
    
    public int getNrOfRegularTransitions() {
	return theVertex.getNrOfTransitions();
    }
    
    public int getNrOfEndTransitions() {
	return theVertex.getNrOfEndTransitions();
    }

    public void addVerticesToVertexHash(Hashtable theVerticesHash)
    {
	theVerticesHash.put(new Integer(theVertex.getNumber()), theVertex);
    }
    
    public boolean addTransition(int fromVertex, int toVertex)
    {
	if(vertexType == HMM.END) {
	    return false;
	}
	else {
	    return theVertex.addTransition(toVertex);
	}
    }

    public void setTransPriorScaler(double d)
    {
	theVertex.setTransPriorScaler(d);
    }

    public void setEmissPriorScaler(double d)
    {
	theVertex.setEmissPriorScaler(d);
    }
    
    public void setEmissPriorScaler(int nr, double d)
    {
	theVertex.setEmissPriorScaler(nr,d);
    }
    
    public void lockVertexEmissions()
    {
	theVertex.lock();
    }

    public boolean addTransition(int toVertex)
    {
	if(vertexType == HMM.END) {
	    return false;
	}
	else {
	    return theVertex.addTransition(toVertex);
	}  
    }

    public void setDistribType(int distribType, double[] distribution)
    {
	this.distribType = distribType;
	theVertex.setInitialEmissionProbs(distribution);
    }

     public void setDistribType(int nr, int distribType, double[] distribution)
    {
	switch(nr) {
	case 1:
	    this.distribType = distribType;
	    
	    break;
	case 2:
	    this.distribType_2 = distribType;
	    break;
	case 3:
	    this.distribType_3 = distribType;
	    break;
	case 4:
	    this.distribType_4 = distribType;
	    break;
	}
	theVertex.setInitialEmissionProbs(nr, distribution);
    }

    public void setDistribType(int distribType)
    {
	this.distribType = distribType;
	theVertex.setInitialEmissionProbs(distribType);
    }

    public void setDistribType(int nr, int distribType)
    {
	switch(nr) {
	case 1:
	    this.distribType = distribType;
	    break;
	case 2:
	    this.distribType_2 = distribType;
	    break;
	case 3:
	    this.distribType_3 = distribType;
	    break;
	case 4:
	    this.distribType_4 = distribType;
	    break;
	}
	theVertex.setInitialEmissionProbs(nr,distribType);
    }

    public void setPriorfile(String pf)
    {
	priorfile = pf;
    }

    public void setPriorfile(int nr, String priorfile)
    {
	switch(nr) {
	case 1:
	    this.priorfile = priorfile;
	    break;
	case 2:
	    this.priorfile_2 = priorfile;
	    break;
	case 3:
	    this.priorfile_3 = priorfile;
	    break;
	case 4:
	    this.priorfile_4 = priorfile;
	}
    }

    public void setTransPriorfile(String priorfile)
    {
	this.transPriorfile = priorfile;
    }

    public boolean addEndTransition(int toVertex)
    {
	return theVertex.addEndTransition(toVertex);
    }

    public void initializeTransitionProbabilities()
    {
	theVertex.initializeTransitionProbabilities();
    }

    public void setInternalInitDistrib(InternalInitDistrib iid)
    {
	/* not implemented for this module */
    }


    public boolean write(BufferedWriter writer)
    {
	try {
	    String s = "Module: " + name + "\n";
	    s = s + "Type: Singleloop\n";
	    s = s + "NrVertices: " + "1\n";
	    s = s + "Emission prior file: " + priorfile + "\n";
	    s = s + "Transition prior file: " + transPriorfile + "\n\n";
	    writer.write(s);
	    s = "Vertex " + theVertex.getNumber() + ":\n";
	    s = s + "Vertex type: " + theVertex.typeAsString() + "\n";
	    s = s + "Vertex label: " + label + "\n";
	    s = s + "Transition prior scaler: " + theVertex.getTransPriorScaler() + "\n";
	    s = s + "Emission prior scaler: " + theVertex.getEmissPriorScaler() + "\n";
	    s = s + "Nr transitions = " + theVertex.getNrOfTransitions() + "\n";
	    writer.write(s);
	    s = "Nr end transitions = " + theVertex.getNrOfEndTransitions() + "\n";
	    s = s + "Nr emissions = " + HMM.alphabet.length + "\n";
	    s = s + "Transition probabilities\n";
	    writer.write(s);
	    for(ListIterator i = theVertex.getTransitions(); i.hasNext();) {
		Transition t = (Transition)i.next();
		s = "\tVertex " + t.toVertex + ": " + t.probability + "\n";
		writer.write(s);
		
	    }
	    
	    s = "End transition probabilities\n";
	    writer.write(s);
	    for(ListIterator i = theVertex.getEndTransitions(); i.hasNext();) {
		Transition t = (Transition)i.next();
		s = "\tVertex " + t.toVertex + ": " + t.probability + "\n";
		writer.write(s);
	    }
	    
	    s = "Emission probabilities\n";
	    writer.write(s);
	    
	    for(int i = 0; i < HMM.alphabet.length; i++) {
		s = "\t" + HMM.alphabet[i] + ": " + theVertex.getEmissionProb(i) + "\n";
		writer.write(s);
	    }
	    s = "\n";
	    s = s + "-------------------------------------------------------\n";
	    writer.write(s);
	    return true;
	}
	catch (IOException e) {
	    return false;
	}
    }

    public boolean write(int nrOfAlphabets, BufferedWriter writer)
    {
	try {
	    String s = "Module: " + name + "\n";
	    s = s + "Type: Singleloop\n";
	    s = s + "NrVertices: " + "1\n";
	    if(nrOfAlphabets == 1) {
		s = s + "Emission prior file: " + priorfile + "\n";
	    }
	    else {
		s = s + "Emission prior file 1: " + priorfile + "\n";
		s = s + "Emission prior file 2: " + priorfile_2 + "\n";
		if(nrOfAlphabets > 2) {
		  s = s + "Emission prior file 3: " + priorfile_3 + "\n";  
		}
		if(nrOfAlphabets > 3) {
		  s = s + "Emission prior file 4: " + priorfile_4 + "\n";  
		}
	    }
	    s = s + "Transition prior file: " + transPriorfile + "\n\n";
	    writer.write(s);
	    s = "Vertex " + theVertex.getNumber() + ":\n";
	    s = s + "Vertex type: " + theVertex.typeAsString() + "\n";
	    s = s + "Vertex label: " + label + "\n";
	    s = s + "Transition prior scaler: " + theVertex.getTransPriorScaler() + "\n";
	    if(nrOfAlphabets == 1) {
		    s = s + "Emission prior scaler: " + theVertex.getEmissPriorScaler() + "\n";
		}
	    else {
		s = s + "Emission prior scaler 1: " + theVertex.getEmissPriorScaler(1) + "\n";
		s = s + "Emission prior scaler 2: " + theVertex.getEmissPriorScaler(2) + "\n";
		if(nrOfAlphabets > 2) {
		    s = s + "Emission prior scaler 3: " + theVertex.getEmissPriorScaler(3) + "\n";	
		}
		if(nrOfAlphabets > 3) {
		    s = s + "Emission prior scaler 4: " + theVertex.getEmissPriorScaler(4) + "\n";	
		}
	    }
	    
	    s = s + "Nr transitions = " + theVertex.getNrOfTransitions() + "\n";
	    writer.write(s);
	    s = "Nr end transitions = " + theVertex.getNrOfEndTransitions() + "\n";
	    if(nrOfAlphabets == 1) {
		s = s + "Nr emissions = " + HMM.alphabet.length + "\n";
	    }
	    else {
		s = s + "Nr emissions 1 = " + HMM.alphabet.length + "\n";
		s = s + "Nr emissions 2 = " + HMM.alphabet_2.length + "\n";
		if(nrOfAlphabets > 2) {
		    s = s + "Nr emissions 3 = " + HMM.alphabet_3.length + "\n";
		}
		if(nrOfAlphabets > 3) {
		    s = s + "Nr emissions 4 = " + HMM.alphabet_4.length + "\n";
		}
	    }
	    s = s + "Transition probabilities\n";
	    writer.write(s);
	    for(ListIterator i = theVertex.getTransitions(); i.hasNext();) {
		Transition t = (Transition)i.next();
		s = "\tVertex " + t.toVertex + ": " + t.probability + "\n";
		writer.write(s);
		
	    }
	    
	    s = "End transition probabilities\n";
	    writer.write(s);
	    for(ListIterator i = theVertex.getEndTransitions(); i.hasNext();) {
		Transition t = (Transition)i.next();
		s = "\tVertex " + t.toVertex + ": " + t.probability + "\n";
		writer.write(s);
	    }
	    
	    if(nrOfAlphabets == 1) {
		    s = "Emission probabilities\n";
		    writer.write(s);
		    for(int i = 0; i < HMM.alphabet.length; i++) {
			s = "\t" + HMM.alphabet[i] + ": " + theVertex.getEmissionProb(i) + "\n";
			writer.write(s);
		    }
		    s = "\n";
		    writer.write(s);
	    }
	    else {
		s = "Emission probabilities 1\n";
		writer.write(s);
		for(int i = 0; i < HMM.alphabet.length; i++) {
		    s = "\t" + HMM.alphabet[i] + ": " + theVertex.getEmissionProb(1,i) + "\n";
		    writer.write(s);
		}
		s = "\n";
		writer.write(s);
		
		s = "Emission probabilities 2\n";
		writer.write(s);
		for(int i = 0; i < HMM.alphabet_2.length; i++) {
		    s = "\t" + HMM.alphabet_2[i] + ": " + theVertex.getEmissionProb(2,i) + "\n";
		    writer.write(s);
		}
		s = "\n";
		writer.write(s);
		
		if(nrOfAlphabets > 2) {
		    s = "Emission probabilities 3\n";
		    writer.write(s);
		    for(int i = 0; i < HMM.alphabet_3.length; i++) {
			s = "\t" + HMM.alphabet_3[i] + ": " + theVertex.getEmissionProb(3,i) + "\n";
			writer.write(s);
		    }
		    s = "\n";
		    writer.write(s);
		}
		
		if(nrOfAlphabets > 3) {
		    s = "Emission probabilities 4\n";
		    writer.write(s);
		    for(int i = 0; i < HMM.alphabet_4.length; i++) {
			s = "\t" + HMM.alphabet_4[i] + ": " + theVertex.getEmissionProb(4,i) + "\n";
			writer.write(s);
		    }
		    s = "\n";
		    writer.write(s);
		}
	    }
	    s = "-------------------------------------------------------\n";
	    writer.write(s);
	    return true;
	}
	catch (IOException e) {
	    return false;
	}
    }
}
