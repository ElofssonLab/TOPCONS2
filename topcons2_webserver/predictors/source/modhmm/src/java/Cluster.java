import java.util.*;
import java.io.*;

class Cluster extends Module
{
    String name;
    LinkedList theVertices;
    int nrVertices;
    int vertexType;
    int distribType;
    int distribType_2;
    int distribType_3;
    int distribType_4;
    String label;
    int size;
    String priorfile;
    String priorfile_2;
    String priorfile_3;
    String priorfile_4;
    String transPriorfile;


    /* Default constructor making all vertices in/out vertices */
    public Cluster(String name, int distribType, int vertexType, int size, String label)
    {
	this.name = name;
	this.vertexType = vertexType;
	this.distribType = distribType;
	this.distribType_2 = distribType;
	this.distribType_3 = distribType;
	this.distribType_4 = distribType;
	this.label = label;
	theVertices = new LinkedList();
	nrVertices = size;
	inVertices = new int[nrVertices];
	outVertices = new int[nrVertices];
	priorfile = null;
	priorfile_2 = null;
	priorfile_3 = null;
	priorfile_4 = null;
	transPriorfile = null;
	
	/* create the vertices */
	Vertex v = new Vertex(name, vertexType, distribType, label);
	theVertices.add(v);
	inVertices[0] = v.getNumber();
	outVertices[0] = v.getNumber();
	for(int i = 1; i < nrVertices; i++) {
	    Vertex w = new Vertex(name, vertexType, v.getEmissionProbs(), label);
	    theVertices.add(w);
	    inVertices[i] = w.getNumber();
	    outVertices[i] = w.getNumber();
	}

	/* add transitions between all vertices */
	for(ListIterator i = theVertices.listIterator();i.hasNext();) {
	    int nr = ((Vertex)i.next()).getNumber();
	    if(!addTransition(nr)) {
		P.INTERNAL_ERROR("Cluster.construct: Could not add intraconnections");
	    }
	}
	size = theVertices.size();
    }

    public Cluster(String name, double[] initDistrib, int vertexType, int size, String label)
    {
	this.name = name;
	this.vertexType = vertexType;
	this.distribType = HMM.MANUAL;
	this.distribType_2 = HMM.MANUAL;
	this.distribType_3 = HMM.MANUAL;
	this.distribType_4 = HMM.MANUAL;
	this.label = label;
	theVertices = new LinkedList();
	nrVertices = size;
	inVertices = new int[nrVertices];
	outVertices = new int[nrVertices];
	priorfile = null;
	priorfile_2 = null;
	priorfile_3 = null;
	priorfile_4 = null;
	transPriorfile = null;
	
	/* create the vertices */
	for(int i = 0; i < nrVertices; i++) {
	    Vertex v = new Vertex(name, vertexType, initDistrib, label);
	    theVertices.add(v);
	    inVertices[i] = v.getNumber();
	    outVertices[i] = v.getNumber();
	    priorfile = null;
	}

	/* add transitions between all vertices */
	for(ListIterator i = theVertices.listIterator();i.hasNext();) {
	    int nr = ((Vertex)i.next()).getNumber();
	    if(!addTransition(nr)) {
		P.INTERNAL_ERROR("Cluster.construct: Could not add intraconnections");
	    }
	}
	size = theVertices.size();
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

    public int getVertexType()
    {
	return vertexType;
    }
    
    public int getDistribType()
    {
	return distribType;
    }

    public String getPriorfile()
    {
	return priorfile;
    }

    public double[] getEmissionProbs()
    {
	Vertex v = (Vertex)theVertices.get(0);
	return v.getEmissionProbs();
    }
    public int[] getInVertices()
    {
	return inVertices;
    }

    public LinkedList getVertices()
    {
	return theVertices;
    }

    public Vertex getVertex(int nr)
    {
	for(int i = 0; i < theVertices.size(); i++) {
	    Vertex v = ((Vertex)theVertices.get(i));
	    if(v.getNumber() == nr) {
		return v;
	    }
	}

	return null;
    }

    public int getNrOfTransitions()
    {
	int nrTransitions = 0;
	for(ListIterator i = (ListIterator)theVertices.iterator(); i.hasNext();) {
	    Vertex v = (Vertex)i.next();
	    nrTransitions = nrTransitions + v.getNrOfTransitions();
	    nrTransitions = nrTransitions + v.getNrOfEndTransitions();
	}
	return nrTransitions;
    }

    public int getNrOfRegularTransitions()
    {
	int nrTransitions = 0;
	for(ListIterator i = (ListIterator)theVertices.iterator(); i.hasNext();) {
	    Vertex v = (Vertex)i.next();
	    nrTransitions = nrTransitions + v.getNrOfTransitions();
	}
	return nrTransitions;
    }

    public int getNrOfEndTransitions()
    {
	int nrTransitions = 0;
	for(ListIterator i = (ListIterator)theVertices.iterator(); i.hasNext();) {
	    Vertex v = (Vertex)i.next();
	    nrTransitions = nrTransitions + v.getNrOfEndTransitions();
	}
	return nrTransitions;
    }

    public int[] getOutVertices()
    {
	return outVertices;
    }

    public void setTransPriorScaler(double d)
    {
	for(ListIterator i = (ListIterator)theVertices.iterator(); i.hasNext();) {
	    Vertex v = (Vertex)i.next();
	    v.setTransPriorScaler(d);
	}
    }

    public void setEmissPriorScaler(double d)
    {
	for(ListIterator i = (ListIterator)theVertices.iterator(); i.hasNext();) {
	    Vertex v = (Vertex)i.next();
	    v.setEmissPriorScaler(d);
	}
    }

    public void setEmissPriorScaler(int nr, double d)
    {
	for(ListIterator i = (ListIterator)theVertices.iterator(); i.hasNext();) {
	    Vertex v = (Vertex)i.next();
	    v.setEmissPriorScaler(nr, d);
	}
    }

    public void lockVertexEmissions()
    {
	for(ListIterator i = (ListIterator)theVertices.iterator(); i.hasNext();) {
	    Vertex v = (Vertex)i.next();
	    v.lock();
	}
    }

    public void addVerticesToVertexHash(Hashtable theVerticesHash)
    {
	for(int i = 0; i < nrVertices; i++) {
	    Vertex v = ((Vertex)theVertices.get(i));
	    theVerticesHash.put(new Integer(v.getNumber()), v);
	}
    }

    /* adds transition from specified vertex to specified vertex */
    public boolean addTransition(int fromVertex, int toVertex)
    {
	for(ListIterator i = (ListIterator)theVertices.iterator(); i.hasNext();) {
	    Vertex v = (Vertex)i.next();
	    if(v.getNumber() == fromVertex) {
		return v.addTransition(toVertex);
	    }
	}
	
	/* could not find the vertex */
	P.INTERNAL_ERROR("Cluster.addTansition: vertex not in module");
	return false;
    }
    
    /* adds transitions from all out vertices to specified vertex */
    public boolean addTransition(int toVertex)
    {
	boolean res = true;
	for(ListIterator i = (ListIterator)theVertices.iterator(); i.hasNext();) {
	    if(!((Vertex)i.next()).addTransition(toVertex)) {
		res = false;
	    }
	}
	return res;
    }

    /* adds end transitions from all out vertices to specified vertex */
    public boolean addEndTransition(int toVertex)
    {
	boolean res = true;
	for(ListIterator i = (ListIterator)theVertices.iterator(); i.hasNext();) {
	    if(!((Vertex)i.next()).addEndTransition(toVertex)) {
		res = false;
	    }
	}
	return res;
    }

    public void initializeTransitionProbabilities()
    {
	for(ListIterator i = (ListIterator)theVertices.iterator(); i.hasNext();) {
	    Vertex v = (Vertex)i.next();
	    v.initializeTransitionProbabilities();
	}
    }

    public void setDistribType(int distribType, double[] distribution)
    {
	this.distribType = distribType;
	for(ListIterator i = (ListIterator)theVertices.iterator();i.hasNext();) {
	    Vertex v = (Vertex)i.next();
	    v.setInitialEmissionProbs(distribution);
	}
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
	
	for(ListIterator i = (ListIterator)theVertices.iterator();i.hasNext();) {
	    Vertex v = (Vertex)i.next();
	    v.setInitialEmissionProbs(nr,distribution);
	}
    }
    
    public void setDistribType(int distribType)
    {
	this.distribType = distribType;
	for(ListIterator i = (ListIterator)theVertices.iterator();i.hasNext();) {
	    Vertex v = (Vertex)i.next();
	    v.setInitialEmissionProbs(distribType);
	}
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
	for(ListIterator i = (ListIterator)theVertices.iterator();i.hasNext();) {
	    Vertex v = (Vertex)i.next();
	    v.setInitialEmissionProbs(nr,distribType);
	}
    }

    public void setPriorfile(String priorfile)
    {
	this.priorfile = priorfile;
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

    public void setInternalInitDistrib(InternalInitDistrib iid)
    {
	/* not implemented for this module */
    }
    
   

    public boolean write(int nrOfAlphabets, BufferedWriter writer)
    {
	try {
	    String s = "Module: " + name + "\n";
	    s = s + "Type: Cluster\n";
	    s = s + "NrVertices: " + nrVertices + "\n";
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
	    
	    for(ListIterator li = (ListIterator)theVertices.iterator(); li.hasNext();) {
		Vertex v = (Vertex)li.next();
		s = "Vertex " + v.getNumber() + ":\n";
		s = s + "Vertex type: " + v.typeAsString() + "\n";
		s = s + "Vertex label: " + v.getLabel() + "\n";
		s = s + "Transition prior scaler: " + v.getTransPriorScaler() + "\n";
		if(nrOfAlphabets == 1) {
		    s = s + "Emission prior scaler: " + v.getEmissPriorScaler() + "\n";
		}
		else {
		    s = s + "Emission prior scaler 1: " + v.getEmissPriorScaler(1) + "\n";
		    s = s + "Emission prior scaler 2: " + v.getEmissPriorScaler(2) + "\n";
		    if(nrOfAlphabets > 2) {
			s = s + "Emission prior scaler 3: " + v.getEmissPriorScaler(3) + "\n";	
		    }
		    if(nrOfAlphabets > 3) {
			s = s + "Emission prior scaler 4: " + v.getEmissPriorScaler(4) + "\n";	
		    }
		}
		s = s + "Nr transitions = " + v.getNrOfTransitions() + "\n";
		s = s + "Nr end transitions = " + v.getNrOfEndTransitions() + "\n";
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
		
		for(ListIterator i = v.getTransitions(); i.hasNext();) {
		    Transition t = (Transition)i.next();
		    s = "\tVertex " + t.toVertex + ": " + t.probability + "\n";
		    writer.write(s);
		}
		s = "End transition probabilities\n";
		writer.write(s);
		for(ListIterator i = v.getEndTransitions(); i.hasNext();) {
		    Transition t = (Transition)i.next();
		    s = "\tVertex " + t.toVertex + ": " + t.probability + "\n";
		    writer.write(s);
		}
		
		if(nrOfAlphabets == 1) {
		    s = "Emission probabilities\n";
		    writer.write(s);
		    for(int i = 0; i < HMM.alphabet.length; i++) {
			s = "\t" + HMM.alphabet[i] + ": " + v.getEmissionProb(i) + "\n";
			writer.write(s);
		    }
		    s = "\n";
		    writer.write(s);
		}
		else {
		    s = "Emission probabilities 1\n";
		    writer.write(s);
		    for(int i = 0; i < HMM.alphabet.length; i++) {
			s = "\t" + HMM.alphabet[i] + ": " + v.getEmissionProb(1,i) + "\n";
			writer.write(s);
		    }
		    s = "\n";
		    writer.write(s);

		    s = "Emission probabilities 2\n";
		    writer.write(s);
		    for(int i = 0; i < HMM.alphabet_2.length; i++) {
			s = "\t" + HMM.alphabet_2[i] + ": " + v.getEmissionProb(2,i) + "\n";
			writer.write(s);
		    }
		    s = "\n";
		    writer.write(s);

		    if(nrOfAlphabets > 2) {
			s = "Emission probabilities 3\n";
			writer.write(s);
			for(int i = 0; i < HMM.alphabet_3.length; i++) {
			    s = "\t" + HMM.alphabet_3[i] + ": " + v.getEmissionProb(3,i) + "\n";
			    writer.write(s);
			}
			s = "\n";
			writer.write(s);
		    }

		    if(nrOfAlphabets > 3) {
			s = "Emission probabilities 4\n";
			writer.write(s);
			for(int i = 0; i < HMM.alphabet_4.length; i++) {
			    s = "\t" + HMM.alphabet_4[i] + ": " + v.getEmissionProb(4,i) + "\n";
			    writer.write(s);
			}
			s = "\n";
			writer.write(s);
		    }
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














