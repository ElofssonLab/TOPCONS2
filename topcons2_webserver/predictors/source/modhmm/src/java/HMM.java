import java.util.*;
import java.io.*;


class HMM
{

    /*alpabets*/
    public static final int AMINO = 1;
    public static final int DNA = 2;
    public static final int AMINO_GAP = 3;
    public static final int AMINO_DOUBLE_GAP = 4;
    public static final int OWN = 100;

    public static final String[] A = {"A","C","D","E","F","G","H","I","K","L",
				      "M","N","P","Q","R","S","T","V","W","Y"};
    public static final String[] D = {"A","C","G","T"};

    public static final String[] A_gap = {"A","C","D","E","F","G","H","I","K","L",
					  "M","N","P","Q","R","S","T","V","W","Y","g"};
    public static final String[] A_double_gap = {"A","C","D","E","F","G","H","I","K","L",
						 "M","N","P","Q","R","S","T","V","W","Y","g","b"};

    /* module types */
    public static final int NOTYPE = -1;
    public static final int SINGLENODE = 0;
    public static final int CLUSTER = 100;
    public static final int HIGHWAY = 107;
    public static final int PROFILE7 = 101;
    public static final int PROFILE9 = 105;
    public static final int SINGLELOOP = 102;
    public static final int FORWARD_STD = 103;
    public static final int FORWARD_ALT = 104;
    public static final int U_TURN = 106;
    public static final int STARTNODE = 1;
    public static final int ENDNODE = 2;

    
    /* init distribution types (for emission probs) */
    public static final int EVEN = 0;
    public static final int ZERO = 1;
    public static final int RANDOM = 2;
    public static final int MANUAL = 3;
    public static final int LOCKED_EVEN = 4;
    public static final int LOCKED_MANUAL = 5;

    /* vertex types */
    public static final int STANDARD = 1;
    public static final int SILENT = 2;
    public static final int LOCKED = 5;
    public static final int INSERT = 3;
    public static final int DELETE = 4;
    public static final int START = 0;
    public static final int END = 99;

   

    /* other constants */
    private final int OKNAME = 0;
    private final int DOUBLENAME = 1;
    private final int INCORRECT_SIZE = 2;

    /* Counter for numbering the vertices */
    private static int vertexNumber = -1;
    public static int GETNUMBER()
    {
        vertexNumber++;
        return vertexNumber;
    }

    

    public static String[] alphabet;
    public static String[] alphabet_2;
    public static String[] alphabet_3;
    public static String[] alphabet_4;


    private Hashtable theModules; /* for quick retrieval of the modules */
    private Hashtable theVertices;
    private LinkedList theModuleList; /* for easier traversals of the modules */
    private LinkedList theDistributionGroups;
    private LinkedList theTransTieGroups;
    private TreeSet thePriorfiles;
    private TreeSet thePriorfiles_2;
    private TreeSet thePriorfiles_3;
    private TreeSet thePriorfiles_4;
    private TreeSet theTransPriorfiles;
    private String name;
    private String filename;
    int nrOfDistributionGroups;
    int nrOfDistributionGroupVertices;
    int nrOfTransTieGroups;
    int nrOfTransTieGroupVertices;
    int nrOfAlphabets;

    public HMM(String name)
    {
        theModules = new Hashtable();
	theVertices = new Hashtable();
        theModuleList = new LinkedList();
        theDistributionGroups = new LinkedList();
	theTransTieGroups = new LinkedList();
	thePriorfiles = new TreeSet();
	thePriorfiles_2 = new TreeSet();
	thePriorfiles_3 = new TreeSet();
	thePriorfiles_4 = new TreeSet();
	theTransPriorfiles = new TreeSet();
        nrOfDistributionGroups = 0;
	nrOfDistributionGroupVertices = 0;
	nrOfTransTieGroups = 0;
	nrOfTransTieGroupVertices = 0;
        this.name = name;
        filename = name + ".hmg";
        P.DEBUG("HMM.C: created HMM '" + name + "'");
	nrOfAlphabets = 1;
	HMM.alphabet = new String[0];
	HMM.alphabet_2 = new String[0];
	HMM.alphabet_3 = new String[0];
	HMM.alphabet_4 = new String[0];
    }

    /* *****************STATIC methods ************************************/
    public static void resetStatics()
    {
	vertexNumber = -1;
    }


    /* *****************GET methods****************************************/ 
    public boolean transitionExists(int fromNr, int toNr)
    {
	Vertex v = getVertex(fromNr);
	if(v == null) {
	    P.MESSAGE("Error: no such vertex exists");
	    System.exit(0);
	}
	else {
	    return v.hasTransitionTo(toNr);
	}
	return false;
    }

    public int getNrModules()
    {
        return theModules.size();
    }

    public int getNrOfDistributionGroups()
    {
        return nrOfDistributionGroups;
    }
    
    public int getNrOfDistributionGroupVertices()
    {
	return nrOfDistributionGroupVertices;
    }

    public int getNrOfTransTieGroups()
    {
	return nrOfTransTieGroups;
    }

    public int getNrOfTransTieGroupVertices()
    {
	return nrOfTransTieGroupVertices;
    }
    
    public int getNrVertices()
    {
        return HMM.vertexNumber + 1;
    }

    public int getNrOfTransitions()
    {
        int nrTransitions = 0;
        for(ListIterator i = (ListIterator)theModuleList.iterator();i.hasNext();) {
            Module m = (Module)i.next();
            nrTransitions = nrTransitions + m.getNrOfTransitions();
        }
        return nrTransitions;
    }
    
    public String getName()
    {
        return name;
    }

    public String getFileName()
    {
        return filename;
    }

    public int getNrOfAlphabets()
    {
	return nrOfAlphabets;
    }

    public String getAlphabet()
    {
        String a = "";
        for(int i = 0; i < alphabet.length;i++) {
            a = a + alphabet[i] + ";";
        }
        return a;
    }
    
    public int getAlphabetSize()
    {
        return alphabet.length;
    }

    public String getAlphabet2()
    {
        String a = "";
        for(int i = 0; i < alphabet_2.length;i++) {
            a = a + alphabet_2[i] + ";";
        }
        return a;
    }
    
    public int getAlphabetSize2()
    {
        return alphabet_2.length;
    }
    
    public String getAlphabet3()
    {
        String a = "";
        for(int i = 0; i < alphabet_3.length;i++) {
            a = a + alphabet_3[i] + ";";
        }
        return a;
    }
    
    public int getAlphabetSize3()
    {
        return alphabet_3.length;
    }
    
    public String getAlphabet4()
    {
        String a = "";
        for(int i = 0; i < alphabet_4.length;i++) {
            a = a + alphabet_4[i] + ";";
        }
        return a;
    }
    
    public int getAlphabetSize4()
    {
        return alphabet_4.length;
    }

    public Enumeration getModuleKeys()
    {
        return theModules.keys();
    }

    public Module getModule(String name)
    {
        return (Module)(theModules.get(name));
    }

    public Vertex getVertex(int nr)
    {
	
	Vertex v = ((Vertex)theVertices.get(new Integer(nr)));
	return v;
	
	/*
	for(ListIterator i = getModules(); i.hasNext();) {
	    Vertex v =  ((Module)i.next()).getVertex(nr);
	    if(v != null) {
		return v;
	    }
	}
	P.MESSAGE("Could not find vertex " + nr);
	System.exit(0);
	return null;
	*/
    }

    public ListIterator getModules()
    {
        return (ListIterator)theModuleList.iterator();
    }

    public ListIterator getDistributionGroups()
    {
        return (ListIterator)theDistributionGroups.iterator();
    }

    public ListIterator getTransTieGroups()
    {
        return (ListIterator)theTransTieGroups.iterator();
    }
    
    public int getNrOfPriorfiles()
    {
	return thePriorfiles.size();
    }

    public int getNrOfPriorfiles(int nr)
    {
	switch(nr) {
	case 1:
	    return thePriorfiles.size();
	case 2:
	    return thePriorfiles_2.size();
	case 3:
	    return thePriorfiles_3.size();
	case 4:
	    return thePriorfiles_4.size();
	}
	P.INTERNAL_ERROR("HMM.getNrOfPriorfiles: incorrect alphabet nr");
	return 0;
    }

    public int getNrOfTransPriorfiles()
    {
	return theTransPriorfiles.size();
    }
    
    public Iterator getPriorfiles()
    {
	return (Iterator)thePriorfiles.iterator();
    }

    public Iterator getPriorfiles(int nr)
    {
	switch(nr) {
	case 1:
	    return (Iterator)thePriorfiles.iterator();
	case 2:	
	    return (Iterator)thePriorfiles_2.iterator();
	case 3:
	    return (Iterator)thePriorfiles_3.iterator();
	case 4:
	    return (Iterator)thePriorfiles_4.iterator();
	}
	P.INTERNAL_ERROR("HMM.getPriorfiles: incorrect alphabet nr");
	return null;
    }

    public Iterator getTransPriorfiles()
    {
	return (Iterator)theTransPriorfiles.iterator();
    }


    /* ************************SET methods********************************/ 
    public void setFileName(String filename)
    {
        this.filename = filename;
    }
    
    public void setNrOfAlphabets(int nr)
    {
	nrOfAlphabets = nr;
    }
    
    public void setAlphabet(int alphabet)
    {
        switch(alphabet){
        case AMINO:
            HMM.alphabet = A;
            break;
        case DNA:
            HMM.alphabet = D;
            break;
	case AMINO_GAP:
	    HMM.alphabet = A_gap;
	    break;
	case AMINO_DOUBLE_GAP:
	    HMM.alphabet = A_double_gap;
	    break;
        default:
            P.INTERNAL_ERROR("HMM.setAlphabet: Unrecognized alphabet type");
        }
    }             
    
    public void setAlphabet(int nr, int alphabet)
    {
        switch(alphabet){
        case AMINO:
	    switch(nr) {
	    case 1:
		HMM.alphabet = A;
		break;
	    case 2:
		HMM.alphabet_2 = A;
		break;
	    case 3:
		HMM.alphabet_3 = A;
		break;
	    case 4:
		HMM.alphabet_4 = A;
		break;	
	    }
	    break;
        case DNA:
	    switch(nr) {
	    case 1:
		HMM.alphabet = D;
		break;
	    case 2:
		HMM.alphabet_2 = D;
		break;
	    case 3:
		HMM.alphabet_3 = D;
		break;
	    case 4:
		HMM.alphabet_4 = D;
		break;	
	    }
	    break;
	case AMINO_GAP:
	    switch(nr) {
	    case 1:
		HMM.alphabet = A_gap;
		break;
	    case 2:
		HMM.alphabet_2 = A_gap;
		break;
	    case 3:
		HMM.alphabet_3 = A_gap;
		break;
	    case 4:
		HMM.alphabet_4 = A_gap;
		break;	
	    }
	    break;
	case AMINO_DOUBLE_GAP:
	    switch(nr) {
	    case 1:
		HMM.alphabet = A_double_gap;
		break;
	    case 2:
		HMM.alphabet_2 = A_double_gap;
		break;
	    case 3:
		HMM.alphabet_3 = A_double_gap;
		break;
	    case 4:
		HMM.alphabet_4 = A_double_gap;
		break;	
	    }
	    break;
        default:
	    System.out.println("nr = " + nr + ": alpha = " + alphabet);
            P.INTERNAL_ERROR("HMM.setAlphabet multi: Unrecognized alphabet type");
        }
    }             
    

    public void setAlphabet(String[] alphabet)
    {
        HMM.alphabet=alphabet;
    }

    public void setAlphabet(int nr, String[] alphabet)
    {
	switch(nr) {
	case 1: 
	    HMM.alphabet=alphabet;
	    break;
	case 2:
	    HMM.alphabet_2=alphabet;
	    break;
	case 3: 
	    HMM.alphabet_3=alphabet;
	    break;
	case 4:
	    HMM.alphabet_4=alphabet;
	    break;
	}
    }
    

    public void setTransition(String fromModule, String toModule)
    {
        /* Get out-vertices from the from-module and in-vertices
           from the to-module and add transition */
        
        Module from = (Module)theModules.get(fromModule);
        Module to = (Module)theModules.get(toModule);
        if(to.getVertexType() == HMM.START) {
            P.MESSAGE("Transition from " + fromModule + " to " + toModule + 
                      " already exists or is illegal");
            return;
        }

        int[] inVertices = to.getInVertices();
        if(to.getVertexType() == HMM.END  && from.getVertexType() != HMM.END) {
            for(int i = 0; i < inVertices.length; i++) {
                if(from.addEndTransition(inVertices[i])) {
                    P.MESSAGE("Added end transition from " + fromModule + " to " + toModule);
                    
                }
                else {
                    P.MESSAGE("End transition from " + fromModule + " to " + toModule + 
                              " already exists or is illegal");
                }
            }
        }
        else {
            for(int i = 0; i < inVertices.length; i++) {
                if(from.addTransition(inVertices[i])) {
                    P.MESSAGE("Added transition from " + fromModule + " to " + toModule);
                }
                else {
                    P.MESSAGE("Transition from " + fromModule + " to " + toModule + 
                              " already exists or is illegal");
                }
            }
        }
    }

    public void setInternalInitDistrib(String name, InternalInitDistrib iid)
    {
	Module m = ((Module)theModules.get(name));
	m.setInternalInitDistrib(iid);
    }

    /********************CREATE and ADD methods*****************************************/
    public int createModule(String name, int moduleType, int distribType, int size, String label, boolean global)
    {
        /* Check for doubles */
        if(theModules.containsKey(name)) {
            return DOUBLENAME;
        }
	if(size <= 0) {
	    return INCORRECT_SIZE;
	}
        Module m;
        switch (moduleType) {
        case HMM.SINGLENODE:
            m = new SingleNode(name, distribType, STANDARD, label);
            theModuleList.add(m);
            theModules.put(name, m);
	    m.addVerticesToVertexHash(theVertices);
            P.MESSAGE("Created singlenode module with name '"+name+"'");
            break;
        case HMM.STARTNODE:
            m = new SingleNode(name, distribType, START, label);
            theModuleList.add(m);
            theModules.put(name, m);
	    m.addVerticesToVertexHash(theVertices);
            P.MESSAGE("Created startnode module with name '"+name+"'");
            break;
        case HMM.ENDNODE:
            m = new SingleNode(name, distribType, END, label);
            theModuleList.add(m);
            theModules.put(name, m);
	    m.addVerticesToVertexHash(theVertices);
            P.MESSAGE("Created endnode module with name '"+name+"'");
            break;
        case HMM.CLUSTER:
            m = new Cluster(name, distribType, STANDARD, size, label);
            theModuleList.add(m);
            theModules.put(name, m);
	    m.addVerticesToVertexHash(theVertices);
            P.MESSAGE("Created cluster module with name '"+name+"'");
            break;
	case HMM.HIGHWAY:
            m = new Highway(name, distribType, STANDARD, size, label);
            theModuleList.add(m);
            theModules.put(name, m);
	    m.addVerticesToVertexHash(theVertices);
            P.MESSAGE("Created highway module with name '"+name+"'");
            break;
	case HMM.SINGLELOOP:
	    m = new SingleLoop(name, distribType, STANDARD, size, label);
	    theModuleList.add(m);
            theModules.put(name, m);
	    m.addVerticesToVertexHash(theVertices);
	    P.MESSAGE("Created singleloop module with name '"+name+"'");
            break;
	case HMM.PROFILE7:
	    if(size < 2) {
		return INCORRECT_SIZE;
	    }
	    m = new Profile7(name, distribType, STANDARD, size, label, global);
	    theModuleList.add(m);
            theModules.put(name, m);
	    m.addVerticesToVertexHash(theVertices);
	    P.MESSAGE("Created profile module with name '"+name+"'");
	    break;
	case HMM.PROFILE9:
	    if(size < 2) {
		return INCORRECT_SIZE;
	    }
	    m = new Profile9(name, distribType, STANDARD, size, label, global);
	    theModuleList.add(m);
            theModules.put(name, m);
	    m.addVerticesToVertexHash(theVertices);
	    P.MESSAGE("Created profile module with name '"+name+"'");
	    break;
	case HMM.U_TURN:
	    if(size < 1) {
		return INCORRECT_SIZE;
	    }
	    m = new U_Turn(name, distribType, STANDARD, size, label);
	    theModuleList.add(m);
            theModules.put(name, m);
	    m.addVerticesToVertexHash(theVertices);
	    P.MESSAGE("Created u-turn module with name '"+name+"'");
	    break;

	    
	    /* add more modules here */
	    
        default:
            P.INTERNAL_ERROR("HMM.createModule: Unrecognized module type");
            
        }
        return OKNAME;
    }
    
    public int createModule(String name, int moduleType, int distribType,
			    int intervalStart, int intervalEnd, String label)
    {
        /* Check for doubles */
        if(theModules.containsKey(name)) {
            return DOUBLENAME;
        }
	if(intervalStart >= intervalEnd) {
	    return INCORRECT_SIZE;
	}
        Module m;
        switch (moduleType) {
	case HMM.FORWARD_STD:
            m = new Forward_std(name, distribType, STANDARD, intervalStart, intervalEnd, label);
            theModuleList.add(m);
            theModules.put(name, m);
            P.MESSAGE("Created forward_std  module with name '"+name+"'");
            break;
	case HMM.FORWARD_ALT:
	    m = new Forward_alt(name, distribType, STANDARD, intervalStart, intervalEnd, label);
            theModuleList.add(m);
            theModules.put(name, m);
            P.MESSAGE("Created forward_std  module with name '"+name+"'");
            break;
	default:
            P.INTERNAL_ERROR("HMM.createModule: Unrecognized module type");
        }
        return OKNAME;
    }
    
    public void addDistributionGroup(LinkedList group)
    {
	theDistributionGroups.add(group);
	/* also check the emission probabilities for the nodes in the same
	 * distribution group, they should be the same, and if they're not,
	 * make them the same by forcing them all to have the same probs as
	 * the first node in the list */
	Module cur = null;
	Module last = null;
	ListIterator i = (ListIterator)group.iterator();
	if(i.hasNext()) {
	    cur = (Module)theModules.get(i.next());
	    nrOfDistributionGroupVertices = nrOfDistributionGroupVertices + cur.getSize();
	}
	for(;i.hasNext();) {
	    last = cur;
	    cur = (Module)theModules.get(i.next());
	    nrOfDistributionGroupVertices = nrOfDistributionGroupVertices + cur.getSize();
	    if(cur.getDistribType() != last.getDistribType()) {
		P.MESSAGE("Warning: Different initial probability types for modules in same " +
			  "distribution group detected - autocorrecting");
		double[] distribution = last.getEmissionProbs();
		int distribType = last.getDistribType();
		cur.setDistribType(distribType, distribution);
	    }
	    else if(cur.getDistribType() == HMM.RANDOM && last.getDistribType() == HMM.RANDOM) {
		double[] distribution = last.getEmissionProbs();
		int distribType = last.getDistribType();
		cur.setDistribType(distribType, distribution);
	    }
	    
	}
	nrOfDistributionGroups++;
    }
    
    public void addPriorfile(String s)
    {
	if(s != null) {
	    thePriorfiles.add(s);
	}
    }

    public void addPriorfile(int nr, String s)
    {
	if(s != null) {
	    switch(nr){
	    case 1: 
		thePriorfiles.add(s);
		break;
	    case 2:
		thePriorfiles_2.add(s);
		break;
	    case 3: 
		thePriorfiles_3.add(s);
		break;
	    case 4:
		thePriorfiles_4.add(s);
		break;	
	    }
	}
    }

    public void addTransPriorfile(String s)
    {
	if(s != null) {
	    theTransPriorfiles.add(s);
	}
    }

    public void addTransTieGroup(LinkedList group)
    {
	theTransTieGroups.add(group);
	ListIterator li = group.listIterator();
	LinkedList nrTransList = new LinkedList();
	while(li.hasNext()) {
	    Module m = getModule(((String)li.next()));
	    if(nrTransList.contains(new Integer(m.getNrOfRegularTransitions()))) {
		nrOfTransTieGroupVertices = nrOfTransTieGroupVertices + m.getNrOfRegularTransitions();
	    }
	    else {
		nrOfTransTieGroups = nrOfTransTieGroups + m.getNrOfRegularTransitions();
		nrTransList.add(new Integer(m.getNrOfRegularTransitions()));
		nrOfTransTieGroupVertices = nrOfTransTieGroupVertices + m.getNrOfRegularTransitions();
	    }
	   
	}
    }



    
    /*************************MISC methods*************************************/
    public void initializeTransitionProbabilities(String name)
    {
	Module module = (Module)theModules.get(name);
	
	if(module == null) /* just double checking */{
	    P.INTERNAL_ERROR("HMM.initializeTransitionProbabilities: module doesn't exist");
	}

	else {
	    module.initializeTransitionProbabilities();
	}
    }


    public boolean writeModule(String moduleName, BufferedWriter writer)
    {
	Module module = (Module)theModules.get(moduleName);
	if(module == null) /* just checking */{
	    P.INTERNAL_ERROR("HMM.moduleToString: module doesn't exist");
	    return false;
	}
	
	else {
	    return module.write(nrOfAlphabets, writer);
	}
    }
    
    public boolean identicalModules(String moduleA, String moduleB)
    {
	/* modules are considered identical if they are of the same type
	 * and the same size */
	Module A = getModule(moduleA);
	Module B = getModule(moduleB);
	if(A instanceof SingleNode && B instanceof SingleNode) {
	    return true;
	}
	else if(A instanceof SingleLoop && B instanceof SingleLoop) {
	    return true;
	}
	else if(A instanceof Forward_std && B instanceof Forward_std) {
	    if(A.getSize() == B.getSize()) {
		return true;
	    }
	}
	else if(A instanceof Forward_alt && B instanceof Forward_alt) {
	    if(A.getSize() == B.getSize()) {
		return true;
	    }
	}
	else if(A instanceof Cluster && B instanceof Cluster) {
	    if(A.getSize() == B.getSize()) {
		return true;
	    }
	}
	else if(A instanceof Profile7 && B instanceof Profile7) {
	    if(A.getSize() == B.getSize()) {
		return true;
	    }
	}
	else if(A instanceof Profile9 && B instanceof Profile9) {
	    if(A.getSize() == B.getSize()) {
		return true;
	    }
	}
	return false;
    }
    
}
