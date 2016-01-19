import java.util.*;
import java.lang.*;

class Vertex
{

    Hashtable transitions, endTransitions;
    LinkedList transitionList, endTransitionList;
    double[] emissions;
    double[] emissions_2;
    double[] emissions_3;
    double[] emissions_4;
    boolean hasEmissions;
    int number;
    int nrOfTransitions;
    int nrOfEndTransitions;
    int vertexType;
    double transPriorScaler;
    double emissPriorScaler;
    double emissPriorScaler_2;
    double emissPriorScaler_3;
    double emissPriorScaler_4;
    String label;

    public Vertex(String moduleName, int vertexType, int distribType, String label)
    {
	number = HMM.GETNUMBER();
	this.vertexType = vertexType;
	transitions = new Hashtable();
	endTransitions = new Hashtable();
	transitionList = new LinkedList();
	endTransitionList = new LinkedList();
	nrOfTransitions = 0;
	nrOfEndTransitions = 0;
	this.hasEmissions = true;
	emissions = new double[HMM.alphabet.length];
	emissions_2 = new double[HMM.alphabet_2.length];
	emissions_3 = new double[HMM.alphabet_3.length];
	emissions_4 = new double[HMM.alphabet_4.length];
	setPreliminaryInitialEmissionProbs(distribType);
	transPriorScaler = 1.0;
	emissPriorScaler = 1.0;
	emissPriorScaler_2 = 1.0;
	emissPriorScaler_3 = 1.0;
	emissPriorScaler_4 = 1.0;
	this.label = label;
    }

    public Vertex(String moduleName, int vertexType, double[] initDistrib, String label)
    {
	number = HMM.GETNUMBER();
	this.vertexType = vertexType;
	transitions = new Hashtable();
	endTransitions = new Hashtable();
	transitionList = new LinkedList();
	endTransitionList = new LinkedList();
	nrOfTransitions = 0;
	nrOfEndTransitions = 0;
	this.hasEmissions = true;
	emissions = new double[HMM.alphabet.length];
	emissions_2 = new double[HMM.alphabet_2.length];
	emissions_3 = new double[HMM.alphabet_3.length];
	emissions_4 = new double[HMM.alphabet_4.length];
	setInitialEmissionProbs(initDistrib);
	transPriorScaler = 1.0;
	emissPriorScaler = 1.0;
	emissPriorScaler_2 = 1.0;
	emissPriorScaler_3 = 1.0;
	emissPriorScaler_4 = 1.0;
	this.label = label;
    }

    public boolean hasTransitionTo(int vNr)
    {
	if(transitions.get(new Integer(vNr)) != null) {
	    return true;
	}
	else {
	    return false;
	}
    }

    public int getNumber()
    {
	return number;
    }
    
    public String getLabel()
    {
	return label;
    }

    public int getVertexType()
    {
	return vertexType;
    }

    public String typeAsString()
    {
	switch(vertexType) {
	case HMM.STANDARD:
	    return "standard";
	case HMM.SILENT:
	    return "silent";
	case HMM.START:
	    return "start";
	case HMM.END:
	    return "end";
	case HMM.LOCKED:
	    return "locked";
	default:
	    P.INTERNAL_ERROR("Vertex.typeAsString: No such type exists");
	    return null;
	    
	}
    }

    public double getTransPriorScaler()
    {
	return transPriorScaler;
    }
    
    public double getEmissPriorScaler()
    {
	return emissPriorScaler;
    }
    
    public double getEmissPriorScaler(int nr)
    {
	switch(nr) {
	case 1:
	    return emissPriorScaler;
	case 2:
	    return emissPriorScaler_2;
	case 3:
	    return emissPriorScaler_3;
	case 4:
	    return emissPriorScaler_4;
	}
	P.INTERNAL_ERROR("Vertex.getEmissPriorScaler: wrong alphabet nr");
	return 0;
    }


    public ListIterator getTransitions()
    {
	return (ListIterator)transitionList.iterator();
    }
    
    public ListIterator getEndTransitions()
    {
	return (ListIterator)endTransitionList.iterator();
    }
    
    public int getNrOfTransitions()
    {
	return nrOfTransitions;
    }
    
    public int getNrOfEndTransitions()
    {
	return nrOfEndTransitions;
    }
    
    public double getEmissionProb(int pos)
    {
	return emissions[pos];
    }

    public double getEmissionProb(int nr, int pos)
    {
	switch(nr) {
	case 1:
	    return emissions[pos];
	case 2:
	    return emissions_2[pos];
	case 3:
	    return emissions_3[pos];
	case 4:
	    return emissions_4[pos];
	}
	P.INTERNAL_ERROR("Vertex.getEmissionProb: incorrect alphabet nr");
	return 0;
    }

    public double[] getEmissionProbs()
    {
	return emissions;
    }


    /*
    public void setEmissionProbability(char letter, double probability)
    {
	boolean foundLetter = false;
	for(int i = 0; i < HMM.alphabet.length; i++) {
	    if(HMM.alphabet[i] == letter) {
		emissions[i] = probability;
		foundLetter = true;
		break;
	    }
	}
	if(!foundLetter) {
	    P.INTERNAL_ERROR("Vertex.setEmissionProbability: No letter '" + letter +
			     "' in alphabet");
	}
	
    }
    */

    public void setLabel(String label)
    {
	this.label = label;
    }

    public void setTransPriorScaler(double d)
    {
	transPriorScaler = d;
    }

    public void setEmissPriorScaler(double d)
    {
	emissPriorScaler = d;
    }

    public void setEmissPriorScaler(int nr, double d)
    {
	switch(nr){
	case 1:
	    emissPriorScaler = d;
	    break;
	case 2:
	    emissPriorScaler_2 = d;
	    break;
	case 3:
	    emissPriorScaler_3 = d;
	    break;
	case 4:
	    emissPriorScaler_3 = d;
	    break;
	}
    }

    public void setTransitionProbability(int toVertex, double probability)
    {
	Transition t = (Transition)transitions.get(new Integer(toVertex));
	if(t != null) {
	    t.probability = probability;
	}
	else {
	    P.INTERNAL_ERROR("Vertex.setTransitionProbability: Transistion from " + number +
			     " to " + toVertex + " does not exist");
	}
    }

    public void setAndLockTransitionProbability(int toVertex, double probability)
    {
	Transition t = (Transition)transitions.get(new Integer(toVertex));
	if(t != null) {
	    t.probability = probability;
	    t.lock();
	}
	else {
	    P.INTERNAL_ERROR("Vertex.setTransitionProbability: Transistion does not exist");
	}
    }

    public boolean addTransition(int toVertex, double probability, boolean isLocked)
    {
	if(!transitions.containsKey(new Integer(toVertex))) {
	    Transition t = new Transition(toVertex, probability, isLocked);
	    transitions.put(new Integer(toVertex), t);
	    transitionList.add(t);
	    nrOfTransitions++;
	    return true;
	}
	else {
	    return false;
	}
    }

    public boolean addTransition(int toVertex)
    {
	if(!transitions.containsKey(new Integer(toVertex))) {
	    Transition t = new Transition(toVertex, 0.0);
	    transitions.put(new Integer(toVertex), t);
	    transitionList.add(t);
	    nrOfTransitions++;
	    return true;
	}
	else {
	    return false;
	}
	
    }

    public boolean addEndTransition(int toVertex)
    {
	if(!endTransitions.containsKey(new Integer(toVertex))) {
	    Transition t = new Transition(toVertex, 1.0);
	    endTransitions.put(new Integer(toVertex),t);
	    endTransitionList.add(t);
	    nrOfEndTransitions++;
	    return true;
	}
	else {
	    return false;
	}
	
    }

    public void initializeTransitionProbabilities()
    {
	/* Initialize transition probabilities evenly for all possible transitions
	 * that are not already locked */
	int nrTransitions = transitions.size();
	if(nrTransitions == 0) {
	    return;
	}
	double initProb = 1.0;
	int nrUnlocked = 0;
	for(ListIterator i = (ListIterator)transitionList.iterator(); i.hasNext();) {
	    Transition t = (Transition)i.next();
	    if(t.isLocked) {
		initProb = initProb - t.probability;
	    }
	    else {
		nrUnlocked++;
	    }
	}
	if(nrUnlocked == 0 || initProb < 0.0 - 0.0000001) {
	    /* if all probs are locked and the sum is not equal to 1.0 or if
	     * not all probs are locked and their sum is larger than 1.0, autocorrect + warn */
	    if(initProb < 0.0 - 0.0000001) {
		P.MESSAGE("Warning: transition probabilities in vertex " + number +
			  " sum to more than one. The error will be autocorrected, but you " +
			  "should check this yourself too");
		P.DEBUG("initprob = " + initProb);
	    }
	    for(ListIterator i = (ListIterator)transitionList.iterator(); i.hasNext();) {
		Transition t = (Transition)i.next();
		if(t.isLocked) {
		    t.probability = t.probability / (1.0 - initProb);
		}
	    }
	    initProb = 0.0;
	    
	}
	initProb = initProb / nrUnlocked;
	
	for(ListIterator i = (ListIterator)transitionList.iterator(); i.hasNext();) {
	    Transition t = (Transition)i.next();
	    if(!t.isLocked) {
		t.probability = initProb;
	    }
	}
    }

    public void setInitialEmissionProbs(double[] initDistrib)
    {
	vertexType = HMM.STANDARD;
	for(int i = 0; i < emissions.length;i++) {
	    emissions[i] = initDistrib[i];
	}
    }

    public void setInitialEmissionProbs(int nr, double[] initDistrib)
    {
	vertexType = HMM.STANDARD;
	switch(nr) {
	case 1:
	    for(int i = 0; i < emissions.length;i++) {
		emissions[i] = initDistrib[i];
	    }
	    break;
	case 2:
	    for(int i = 0; i < emissions_2.length;i++) {
		emissions_2[i] = initDistrib[i];
	    }
	    break;
	case 3:
	    for(int i = 0; i < emissions_3.length;i++) {
		emissions_3[i] = initDistrib[i];
	    }
	    break;
	case 4:
	    for(int i = 0; i < emissions_4.length;i++) {
		emissions_4[i] = initDistrib[i];
	    }
	    break;
	}
    }
    
    public void lock()
    {
	if(vertexType == HMM.STANDARD) {
	    vertexType = HMM.LOCKED;
	}
    }

    public void normalizeEmissProbs(int alphabet)
    {
	double[] emiss = null;
	
	switch(alphabet) {
	case 1:
	    emiss = emissions;
	    break;
	case 2:
	    emiss = emissions_2;
	    break;
	case 3:
	    emiss = emissions_3;
	    break;
	case 4:
	    emiss = emissions_4;
	    break;  
	}
	double sum = 0.0;
	for(int i = 0; i < emiss.length; i++) {
	    sum += emiss[i];
	}
	if(sum != 0.0) {
	    for(int i = 0; i < emiss.length; i++) {
		emiss[i] = emiss[i] / sum;
	    }
	}
    }

    public void normalizeTransProbs()
    {
	double sum = 0.0;
	for(int i = 0; i < transitionList.size(); i++) {
	    sum += ((Transition)transitionList.get(i)).probability;
	}
	if(sum != 0.0) {
	    for(int i = 0; i < transitionList.size(); i++) {
		((Transition)transitionList.get(i)).probability =  ((Transition)transitionList.get(i)).probability / sum;
	    }
	}
    }
    
    public void setInitialEmissionProbs(int distribType)
    {
	/* Initialize emission probabilities to zero for all possible emissions */
	if(distribType == HMM.ZERO) {
	    if(vertexType != HMM.START && vertexType != HMM.END) {
		vertexType = HMM.SILENT;
	    }
	    for(int i = 0; i < emissions.length; i++) {
		emissions[i] = 0.0;
	    }
	}
	/* Initialize emission probabilities evenly for all possible emissions */
	else if (distribType == HMM.EVEN) {
	    vertexType = HMM.STANDARD;
	    int nrEmissions = emissions.length;
	    double initProb = (double)1/nrEmissions;
	    for(int i = 0; i < emissions.length; i++) {
		emissions[i] = initProb;
	    }
	}

	/* Initialize emission probabilities randomly */
	else if (distribType == HMM.RANDOM) {
	    vertexType = HMM.STANDARD;
	    int nrEmissions = emissions.length;
	    double sum = 0;
	    for(int i = 0; i < emissions.length; i++) {
		double initProb = Math.random();
		emissions[i] = initProb;
		sum = sum + initProb;
	    }
	    for(int i = 0; i < emissions.length; i++) {
		emissions[i] = emissions[i] / sum;
	    }
	    
	}
    }

    public void setInitialEmissionProbs(int nr, int distribType)
    {
	switch(nr) {
	case 1:
	    /* Initialize emission probabilities to zero for all possible emissions */
	    if(distribType == HMM.ZERO) {
		if(vertexType != HMM.START && vertexType != HMM.END) {
		    vertexType = HMM.SILENT;
		}
		for(int i = 0; i < emissions.length; i++) {
		    emissions[i] = 0.0;
		}
	    }
	    /* Initialize emission probabilities evenly for all possible emissions */
	    else if (distribType == HMM.EVEN) {
		vertexType = HMM.STANDARD;
		int nrEmissions = emissions.length;
		double initProb = (double)1/nrEmissions;
		for(int i = 0; i < emissions.length; i++) {
		    emissions[i] = initProb;
		}
	    }
	    
	    /* Initialize emission probabilities randomly */
	    else if (distribType == HMM.RANDOM) {
		vertexType = HMM.STANDARD;
		int nrEmissions = emissions.length;
		double sum = 0;
		for(int i = 0; i < emissions.length; i++) {
		    double initProb = Math.random();
		    emissions[i] = initProb;
		    sum = sum + initProb;
		}
		for(int i = 0; i < emissions.length; i++) {
		    emissions[i] = emissions[i] / sum;
		}
		
	    }
	    break;
	case 2:
	    /* Initialize emission probabilities to zero for all possible emissions */
	    if(distribType == HMM.ZERO) {
		if(vertexType != HMM.START && vertexType != HMM.END) {
		    vertexType = HMM.SILENT;
		}
		for(int i = 0; i < emissions_2.length; i++) {
		    emissions_2[i] = 0.0;
		}
	    }
	    /* Initialize emission probabilities evenly for all possible emissions */
	    else if (distribType == HMM.EVEN) {
		vertexType = HMM.STANDARD;
		int nrEmissions = emissions_2.length;
		double initProb = (double)1/nrEmissions;
		for(int i = 0; i < emissions_2.length; i++) {
		    emissions_2[i] = initProb;
		}
	    }
	    
	    /* Initialize emission probabilities randomly */
	    else if (distribType == HMM.RANDOM) {
		vertexType = HMM.STANDARD;
		int nrEmissions = emissions_2.length;
		double sum = 0;
		for(int i = 0; i < emissions_2.length; i++) {
		    double initProb = Math.random();
		    emissions_2[i] = initProb;
		    sum = sum + initProb;
		}
		for(int i = 0; i < emissions_2.length; i++) {
		    emissions_2[i] = emissions_2[i] / sum;
		}
		
	    }
	    break;
	case 3:
	    /* Initialize emission probabilities to zero for all possible emissions */
	    if(distribType == HMM.ZERO) {
		if(vertexType != HMM.START && vertexType != HMM.END) {
		    vertexType = HMM.SILENT;
		}
		for(int i = 0; i < emissions_3.length; i++) {
		    emissions_3[i] = 0.0;
		}
	    }
	    /* Initialize emission probabilities evenly for all possible emissions */
	    else if (distribType == HMM.EVEN) {
		vertexType = HMM.STANDARD;
		int nrEmissions = emissions_3.length;
		double initProb = (double)1/nrEmissions;
		for(int i = 0; i < emissions_3.length; i++) {
		    emissions_3[i] = initProb;
		}
	    }
	    
	    /* Initialize emission probabilities randomly */
	    else if (distribType == HMM.RANDOM) {
		vertexType = HMM.STANDARD;
		int nrEmissions = emissions_3.length;
		double sum = 0;
		for(int i = 0; i < emissions_3.length; i++) {
		    double initProb = Math.random();
		    emissions_3[i] = initProb;
		    sum = sum + initProb;
		}
		for(int i = 0; i < emissions_3.length; i++) {
		    emissions_3[i] = emissions_3[i] / sum;
		}
		
	    }
	    break;
	case 4:
	    /* Initialize emission probabilities to zero for all possible emissions */
	    if(distribType == HMM.ZERO) {
		if(vertexType != HMM.START && vertexType != HMM.END) {
		    vertexType = HMM.SILENT;
		}
		for(int i = 0; i < emissions_4.length; i++) {
		    emissions_4[i] = 0.0;
		}
	    }
	    /* Initialize emission probabilities evenly for all possible emissions */
	    else if (distribType == HMM.EVEN) {
		vertexType = HMM.STANDARD;
		int nrEmissions = emissions_4.length;
		double initProb = (double)1/nrEmissions;
		for(int i = 0; i < emissions_4.length; i++) {
		    emissions_4[i] = initProb;
		}
	    }
	    
	    /* Initialize emission probabilities randomly */
	    else if (distribType == HMM.RANDOM) {
		vertexType = HMM.STANDARD;
		int nrEmissions = emissions_4.length;
		double sum = 0;
		for(int i = 0; i < emissions_4.length; i++) {
		    double initProb = Math.random();
		    emissions_4[i] = initProb;
		    sum = sum + initProb;
		}
		for(int i = 0; i < emissions_4.length; i++) {
		    emissions_4[i] = emissions_4[i] / sum;
		}
		
	    }
	    break;
	}
    }

    public void setPreliminaryInitialEmissionProbs(int distribType)
    {
	/* Initialize emission probabilities to zero for all possible emissions */
	if(distribType == HMM.ZERO) {
	    if(vertexType != HMM.START && vertexType != HMM.END) {
		vertexType = HMM.SILENT;
	    }
	    for(int i = 0; i < emissions.length; i++) {
		emissions[i] = 0.0;
	    }
	    for(int i = 0; i < emissions_2.length; i++) {
		emissions_2[i] = 0.0;
	    }
	    for(int i = 0; i < emissions_3.length; i++) {
		emissions_3[i] = 0.0;
	    }
	    for(int i = 0; i < emissions_4.length; i++) {
		emissions_4[i] = 0.0;
	    }
	}
	/* Initialize emission probabilities evenly for all possible emissions */
	else if (distribType == HMM.EVEN) {
	    vertexType = HMM.STANDARD;
	    int nrEmissions = emissions.length;
	    double initProb = (double)1/nrEmissions;
	    for(int i = 0; i < emissions.length; i++) {
		emissions[i] = initProb;
	    }
	    nrEmissions = emissions_2.length;
	    if(nrEmissions != 0) {
		initProb = (double)1/nrEmissions;
	    }
	    else {
		initProb = 0.0;
	    }
	    for(int i = 0; i < emissions_2.length; i++) {
		emissions_2[i] = initProb;
	    }
	    nrEmissions = emissions_3.length;
	    if(nrEmissions != 0) {
		initProb = (double)1/nrEmissions;
	    }
	    else {
		initProb = 0.0;
	    }
	    for(int i = 0; i < emissions_3.length; i++) {
		emissions_3[i] = initProb;
	    }
	    nrEmissions = emissions_4.length;
	    if(nrEmissions != 0) {
		initProb = (double)1/nrEmissions;
	    }
	    else {
		initProb = 0.0;
	    }
	    for(int i = 0; i < emissions_4.length; i++) {
		emissions_4[i] = initProb;
	    }
	    
	}

	/* Initialize emission probabilities randomly */
	else if (distribType == HMM.RANDOM) {
	    vertexType = HMM.STANDARD;
	    int nrEmissions = emissions.length;
	    double sum = 0;
	    for(int i = 0; i < emissions.length; i++) {
		double initProb = Math.random();
		emissions[i] = initProb;
		sum = sum + initProb;
	    }
	    for(int i = 0; i < emissions.length; i++) {
		emissions[i] = emissions[i] / sum;
	    }
	    nrEmissions = emissions_2.length;
	    sum = 0;
	    for(int i = 0; i < emissions_2.length; i++) {
		double initProb = Math.random();
		emissions_2[i] = initProb;
		sum = sum + initProb;
	    }
	    for(int i = 0; i < emissions_2.length; i++) {
		emissions_2[i] = emissions_2[i] / sum;
	    }
	    nrEmissions = emissions_3.length;
	    sum = 0;
	    for(int i = 0; i < emissions_3.length; i++) {
		double initProb = Math.random();
		emissions_3[i] = initProb;
		sum = sum + initProb;
	    }
	    for(int i = 0; i < emissions_3.length; i++) {
		emissions_3[i] = emissions_3[i] / sum;
	    }
	    nrEmissions = emissions_4.length;
	    sum = 0;
	    for(int i = 0; i < emissions_4.length; i++) {
		double initProb = Math.random();
		emissions_4[i] = initProb;
		sum = sum + initProb;
	    }
	    for(int i = 0; i < emissions_4.length; i++) {
		emissions_4[i] = emissions_4[i] / sum;
	    }
	    
	}
    }
}
