import java.util.*;
import java.io.*;

class Profile7 extends Module
{

    final double INIT_DD_PROB = 0.8;
    final double INIT_DM_PROB = 0.2;
    final double INIT_MD_PROB = 0.025;
    final double INIT_MM_PROB = 0.925;
    final double INIT_MI_PROB = 0.025;
    final double INIT_EXIT_PROB = 0.025;
    final double INIT_IM_PROB = 0.2;
    final double INIT_II_PROB = 0.8;
    final double INIT_NA_CLOSE_PROB = 0.01;
    final double INIT_NA_LOOP_PROB = 0.99;
    final double INIT_NA_OPEN_PROB = 0.01;
    final double INIT_A_OPEN_PROB = 0.99;



    String name;
    LinkedList theVertices;
    LinkedList theMatchVertices;
    LinkedList theDeleteVertices;
    LinkedList theInsertVertices;
    Vertex[] theOutVertices;
    int nrVertices;
    int vertexType;
    int distribType;
    int distribType_2;
    int distribType_3;
    int distribType_4;
    String label;
    String priorfile;
    String priorfile_2;
    String priorfile_3;
    String priorfile_4;
    String transPriorfile;
    int size;
    


    /* Default constructor making all vertices in/out vertices */
    public Profile7(String name, int distribType, int vertexType, int size, String label, boolean global)
    {
	createProfile7Module(name, vertexType, size, label, distribType, null, true, global);
    }

    public Profile7(String name, double[] initDistrib, int vertexType, int size, String label, boolean global)
    {
	createProfile7Module(name, vertexType, size, label, 0, initDistrib, false, global);
    }

    private void createProfile7Module(String name, int vertexType, int size, String label, int distribType,
				      double[] initDistrib, boolean autoDistrib, boolean global)
    {
	this.name = name;
	this.vertexType = vertexType;
	if(autoDistrib) {
	    this.distribType = distribType;
	    this.distribType_2 = distribType;
	    this.distribType_3 = distribType;
	    this.distribType_4 = distribType;
	}
	this.label = label;
	theVertices = new LinkedList();
	theDeleteVertices = new LinkedList();
	theMatchVertices = new LinkedList();
	theInsertVertices = new LinkedList();
	nrVertices = size * 3 + 4;
	inVertices = new int[1];
	outVertices = new int[1];
	theOutVertices = new Vertex[1];
	priorfile = null;
	priorfile_2 = null;
	priorfile_3 = null;
	priorfile_4 = null;
	transPriorfile = null;
	
	

	/* create the vertices */
	Vertex d0 = new Vertex(name, HMM.SILENT, HMM.ZERO, "D");
	theVertices.add(d0);
	theDeleteVertices.add(d0);
	inVertices[0] = d0.getNumber();
	Vertex N = null;
	if(autoDistrib) {
	    N = new Vertex(name, HMM.STANDARD, distribType, "I");
	}
	else {
	    N = new Vertex(name, HMM.STANDARD, initDistrib, "I");
	}
	theVertices.add(N);
	theInsertVertices.add(N);
	
	Vertex d1 = new Vertex(name, HMM.SILENT, HMM.ZERO, "D");
	theVertices.add(d1);
	theDeleteVertices.add(d1);
	Vertex m1 = null;
	if(autoDistrib) {
	    m1 = new Vertex(name, HMM.STANDARD, distribType, "M");
	}
	else {
	    m1 = new Vertex(name, HMM.STANDARD, initDistrib, "M");
	}
	theVertices.add(m1);
	theMatchVertices.add(m1);
	Vertex i1 = null;
	if(autoDistrib) {
	    i1 = new Vertex(name, HMM.STANDARD, distribType, "I");
	}
	else {
	    i1 = new Vertex(name, HMM.STANDARD, initDistrib, "I");
	}
	theVertices.add(i1);
	theInsertVertices.add(i1);
	for(int j = 1; j < (nrVertices - 3)/3; j++) {
	    Vertex d = new Vertex(name, HMM.SILENT, HMM.ZERO, "D");
	    theVertices.add(d);
	    theDeleteVertices.add(d);
	    Vertex m = new Vertex(name, HMM.STANDARD, m1.getEmissionProbs(), "M");
	    theVertices.add(m);
	    theMatchVertices.add(m);
	    if(j < (nrVertices - 6)/3) {
		Vertex i = new Vertex(name, HMM.STANDARD, i1.getEmissionProbs(), "I");
		theVertices.add(i);
		theInsertVertices.add(i);
	    }
	}
	Vertex d_last =  new Vertex(name, HMM.SILENT, HMM.ZERO, "D");
	theVertices.add(d_last);
	theDeleteVertices.add(d_last);
	
	Vertex i = new Vertex(name, HMM.STANDARD, i1.getEmissionProbs(), "I");
	theVertices.add(i);
	theInsertVertices.add(i);
	
	Vertex d = new Vertex(name, HMM.SILENT, HMM.ZERO, "D");
	theVertices.add(d);
	theDeleteVertices.add(d);
	outVertices[0] = d.getNumber();
	theOutVertices[0] = d;
	int d_last_nr = d_last.getNumber();


	/* add transitions between vertices according to plan7 architecture, that is 
	 * no delete-insert or insert-delete transitions */
	ListIterator ds = theDeleteVertices.listIterator();
	ListIterator ms = theMatchVertices.listIterator();
	ListIterator is = theInsertVertices.listIterator();
	d = null;
	Vertex d_first = null;
	Vertex m = null;
	i = null;
	int d_first_nr = 0;
	int m_nr = 0;
	int i_nr = 0;
	int d_nr = 0;


	if(ds.hasNext() && ms.hasNext() && is.hasNext()) {
	    d_first = (Vertex)ds.next();
	    if(ds.hasNext()) {
		d = (Vertex)ds.next();
	    }
	    m = (Vertex)ms.next();
	    i = (Vertex)is.next();
	    d_first_nr = d_first.getNumber();
	    m_nr = m.getNumber();
	    i_nr = i.getNumber();
	    d_nr = d.getNumber();
	    i.addTransition(d_first_nr, INIT_NA_CLOSE_PROB, true);
	    i.addTransition(i_nr, INIT_NA_LOOP_PROB, true);
	    d_first.addTransition(i_nr, INIT_NA_OPEN_PROB, true);
	    d_first.addTransition(m_nr, INIT_A_OPEN_PROB / ((double)((size + 1))), true);
	    d_first.addTransition(d_nr, INIT_A_OPEN_PROB / ((double)((size + 1))), true);
	}
	
	while(ds.hasNext() && ms.hasNext() && is.hasNext()) {
	    Vertex d_next = (Vertex)ds.next();
	    Vertex m_next = (Vertex)ms.next();
	    Vertex i_next = (Vertex)is.next();
	    int d_next_nr = d_next.getNumber();
	    int m_next_nr = m_next.getNumber();
	    int i_next_nr = i_next.getNumber();
	    d.addTransition(d_next_nr, INIT_DD_PROB, true);
	    d.addTransition(m_next_nr, INIT_DM_PROB, true);
	    m.addTransition(d_next_nr, INIT_MD_PROB, true);
	    m.addTransition(m_next_nr, INIT_MM_PROB, true);
	    m.addTransition(i_next_nr, INIT_MI_PROB, true);
	    if(!global) {
		m.addTransition(d_last_nr, INIT_EXIT_PROB, true);
	    }
	    i_next.addTransition(m_next_nr, INIT_IM_PROB, true);
	    i_next.addTransition(i_next_nr, INIT_II_PROB, true);
	    if(!global) {
		d_first.addTransition(m_next_nr, INIT_A_OPEN_PROB / ((double)(size + 1)), true);
	    }
	    d = d_next;
	    m = m_next;
	    i = i_next;
	    d_nr = d_next_nr;
	    m_nr = m_next_nr;
	    i_nr = i_next_nr;
	}
	
	if(ds.hasNext()) {
	    d_last = (Vertex)ds.next();
	    d_last_nr = d_last.getNumber();
	    m.addTransition(d_last_nr, 1.0, true);
	    d.addTransition(d_last_nr, 1.0, true);
	}

	if(is.hasNext() && ds.hasNext()) {
	    i = (Vertex)is.next();
	    i_nr = i.getNumber();
	    d = (Vertex)ds.next();
	    d_nr = d.getNumber();
	    i.addTransition(d_last_nr, INIT_NA_CLOSE_PROB, true);
	    i.addTransition(i_nr, INIT_NA_LOOP_PROB, true);
	    d_last.addTransition(i_nr, INIT_NA_OPEN_PROB, true);
	    d_last.addTransition(d_nr, INIT_EXIT_PROB, true);
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
    
    public void setEmissPriorScalerInsert(int nr, double d)
    {
	for(ListIterator i = (ListIterator)theVertices.iterator(); i.hasNext();) {
	    Vertex v = (Vertex)i.next();
	    if(v.getLabel().equals("I")) {
		v.setEmissPriorScaler(nr, d);
	    }
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
	for(int i = 0; i < outVertices.length; i++) {
	    if(!theOutVertices[i].addTransition(toVertex)) {
		res = false;
	    }
	}
	return res;
    }

    /* adds end transitions from all out vertices to specified vertex */
    public boolean addEndTransition(int toVertex)
    {
	boolean res = true;
	for(int i = 0; i < outVertices.length; i++) {
	    if(!theOutVertices[i].addEndTransition(toVertex)) {
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
	for(ListIterator i = (ListIterator)theMatchVertices.iterator();i.hasNext();) {
	    Vertex v = (Vertex)i.next();
	    v.setInitialEmissionProbs(distribution);
	}
	for(ListIterator i = (ListIterator)theInsertVertices.iterator();i.hasNext();) {
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
	for(ListIterator i = (ListIterator)theMatchVertices.iterator();i.hasNext();) {
	    Vertex v = (Vertex)i.next();
	    v.setInitialEmissionProbs(nr,distribution);
	}
	for(ListIterator i = (ListIterator)theInsertVertices.iterator();i.hasNext();) {
	    Vertex v = (Vertex)i.next();
	    v.setInitialEmissionProbs(nr,distribution);
	}
    }
    
    public void setDistribType(int distribType)
    {
	this.distribType = distribType;
	for(ListIterator i = (ListIterator)theMatchVertices.iterator();i.hasNext();) {
	    Vertex v = (Vertex)i.next();
	    v.setInitialEmissionProbs(distribType);
	}
	for(ListIterator i = (ListIterator)theInsertVertices.iterator();i.hasNext();) {
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
	for(ListIterator i = (ListIterator)theMatchVertices.iterator();i.hasNext();) {
	    Vertex v = (Vertex)i.next();
	    v.setInitialEmissionProbs(nr,distribType);
	}
	for(ListIterator i = (ListIterator)theInsertVertices.iterator();i.hasNext();) {
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

    public boolean write(BufferedWriter writer)
    {
	try {
	    String s = "Module: " + name + "\n";
	    s = s + "Type: Profile7\n";
	    s = s + "NrVertices: " + nrVertices + "\n";
	    s = s + "Emission prior file: " + priorfile + "\n";
	    s = s + "Transition prior file: " + transPriorfile + "\n\n";
	    writer.write(s);
	    
	    int v_nr = 0;
	    for(ListIterator li = (ListIterator)theVertices.iterator(); li.hasNext(); v_nr++) {
		Vertex v = (Vertex)li.next();
		s = "Vertex " + v.getNumber() + ":\n";
		s = s + "Vertex type: " + v.typeAsString() + "\n";
		s = s + "Vertex label: " + v.getLabel() + "\n";
		s = s + "Transition prior scaler: " + v.getTransPriorScaler() + "\n";
		s = s + "Emission prior scaler: " + v.getEmissPriorScaler() + "\n";
		s = s + "Nr transitions = " + v.getNrOfTransitions() + "\n";
		s = s + "Nr end transitions = " + v.getNrOfEndTransitions() + "\n";
		s = s + "Nr emissions = " + HMM.alphabet.length + "\n";
		
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
		
		s = "Emission probabilities\n";
		writer.write(s);
		for(int i = 0; i < HMM.alphabet.length; i++) {
		    s = "\t" + HMM.alphabet[i] + ": " + v.getEmissionProb(i) + "\n";
		    writer.write(s);
		}
		s = "\n";
		writer.write(s);
	    }
	    
	    s = "-------------------------------------------------------\n";
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
	    s = s + "Type: Profile7\n";
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
	    
	    int v_nr = 0;
	    for(ListIterator li = (ListIterator)theVertices.iterator(); li.hasNext(); v_nr++) {
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
