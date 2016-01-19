import java.io.*;
import java.util.*;

class DataReaderMultiAlpha
{
    

    BufferedReader reader;
    ModelMaker modelmaker;
    
    public DataReaderMultiAlpha()
    {
	reader = new BufferedReader(new InputStreamReader(System.in));
	modelmaker = new ModelMaker();
	
	println("Welcome to modhmmc!");
	createHMM();
	int nrAlpha = specifyNrOfAlphabets();
	for(int i = 0; i < nrAlpha; i++) {
	    specifyAlphabet(i+1);
	}
	specifyModules();
	specifyInterConnectivity();
	specifyDistributionGroups();
	specifyTransitionTies();
	for(int i = 0; i < nrAlpha; i++) {
	    specifyInitvalues(i+1);
	}
	specifyTransitionInitvalues();
	cleanUp();
	println("Creating untrained HMM...");
	saveHMM("");
	println("HMM '"+ modelmaker.getName() + "' saved to '" +
		modelmaker.getFileName() + "'.");
    }


    /****************input/output handling**********************************/
    private void println(String s)
    {
	System.out.println(s);
    }

    private void newln(){
	System.out.println("");
    }
    private void printerr(String s)
    {
	System.out.println("Error: " + s);
    }
    private void print(String s)
    {
	System.out.print(s);
	System.out.flush();
    }

    private String readln()
    {
	try {
	String s = reader.readLine();
	s = s.trim();
	return s;
	}
	catch(IOException e) {
	    P.INTERNAL_ERROR("DataReader.readLine: IOException");
	    return null;
	}
    }


    /*******************************HMM methods****************************/
    private void createHMM()
    {
	print("Name of HMM? ");
	String name = readln();
	name = name.trim();
	modelmaker.createHMM(name);
	newln();
    }

    /********************************************************************************************
     ********************** methods for getting alphabet specification ****************************
     ********************************************************************************************/
    private int specifyNrOfAlphabets()
    {
	print("Nr of alphabets (1-4)? ");
	String s = readln();
	s = s.trim();
	try {
	    int nr = Integer.parseInt(s);
	    if(nr > 0 && nr <= 4) {
		modelmaker.setNrOfAlphabets(nr);
		return nr;
	    }
	    else {
		printerr("Unallowed nr of alphabets");
		return specifyNrOfAlphabets();
	    }
	}
	catch (NumberFormatException e) {
	    printerr("Not an integer");
	    return specifyNrOfAlphabets();
	}
    }
    
    private void specifyAlphabet(int nr)
    {

	print("Alphabet " + nr + " (A:N:G:S)? ");
	String s = readln();
	if(s.equals("A")){
	    modelmaker.setAlphabet(nr,HMM.AMINO);
	    println("Alphabet set to amino");
	    newln();
	}
	else if(s.equals("N")){
	    modelmaker.setAlphabet(nr,HMM.DNA);
	    println("Alphabet set to nucleotides");
	    newln();
	}
	else if(s.equals("G")){
	    modelmaker.setAlphabet(nr,HMM.AMINO_GAP);
	    println("Alphabet set to amino + gap");
	    newln();
	}
	else if(s.equals("S")){
	    modelmaker.setAlphabet(nr,getSpecAlphabet());
	    println("Alphabet set to custom");
	    newln();
	}
	
	else {
	    printerr("Invalid option");
	    specifyAlphabet(nr);
	}	
    }

    private String[] getSpecAlphabet()
    {
	print("Specify alphabet (separate with ';'): ");
	String s = readln();
	StringTokenizer st = new StringTokenizer(s,";");
	String[] alphabet = new String[st.countTokens()];
	int i=0; /* counter for the alphabet array */
	while (st.hasMoreTokens()){
	    String l = st.nextToken();
	    l = l.trim();
	    if(l.length() <= 4) /* letter OK*/ {
		alphabet[i] = l;
		i++;
	    }
	    else /* letter not OK */ {
		printerr("Only symbols of 4 characters or less are allowed in alphabet");
		return getSpecAlphabet();
		
	    }
	}
	
	/* Check for doubles */
	for(int k = 0 ; k < alphabet.length - 1; k++) {
	    for (int j = k+1; j < alphabet.length; j++) {
		if(alphabet[k].equals( alphabet[j])) {
		    printerr("Alphabet contains doubles");
		    return getSpecAlphabet();
		}
	    }
	}
	return alphabet;
    }

    /********************************************************************************************
     ********************** methods for getting module specification ****************************
     ********************************************************************************************/
    private void specifyModules()
    {
	println("Specify modules");

	/* start node */
	print("Start node: ");
	String str = readln();
	if(str.length() == 0) {
	    int res = modelmaker.createModule("s", HMM.STARTNODE, HMM.ZERO, 1, "0");
	}
	else {
	    int res = modelmaker.createModule(str, HMM.STARTNODE, HMM.ZERO, 1, "0");
	}
	
	/* intermediate nodes */
	int i = 1;
	while(true) /* loop until user is done */ {
	    print("Module " + i + " (type,name,label): ");
	    String s = readln();
	    if(s.length() == 0) {
		break;
	    }
	    else if(s.equals("d")){
		break;
	    }
	    else if(s.equals("done")){
		break;
	    }
	    StringTokenizer st = new StringTokenizer(s,",");
	    if(st.countTokens() != 2 && st.countTokens() != 3) {
		printerr("illegal choice");
		continue;
	    }
	    
	    String t = st.nextToken();
	    t = t.trim();
	    String name = st.nextToken();
	    name = name.trim();
	    String label = "1";
	    if(st.hasMoreTokens()) {
		label = st.nextToken();
		label = label.substring(0,1);
	    }
	    int type = parseType(t);
	    int size = 0;
	    int intervalStart = 0; /* only for forward and forward_std modules */
	    int intervalEnd = 0; /*   ------------------ " -------------------------*/
	    InternalInitDistrib  internalInitDistrib = new InternalInitDistrib(); /* different meaning for different */
	    switch(type) {
	    case HMM.SINGLENODE: size = 1; break;
	    case HMM.CLUSTER:
	    case HMM.HIGHWAY:	
	    case HMM.U_TURN:
	    case HMM.PROFILE7:
	    case HMM.PROFILE9: size = getSize(); break;
	    case HMM.SINGLELOOP: size = getLength(); break;
	    case HMM.FORWARD_STD:
	    case HMM.FORWARD_ALT:
		while(true) {
		    String intervals = getIntervals();
		    intervalStart = getIntervalStart(intervals);
		    intervalEnd = getIntervalEnd(intervals);
		    if(intervalStart >= 0 && intervalEnd >= 0) {
			break;
		    }
		    else {
			printerr("illegal choice");
		    }
		}
		getInternalInitDistrib(HMM.FORWARD_STD, internalInitDistrib);
		break;
	    }
	    if(type != HMM.NOTYPE && type != HMM.FORWARD_STD && type != HMM.FORWARD_ALT) {
		/* create most module types */
		int res = modelmaker.createModule(name, type, HMM.EVEN, size, label);
		if(res == 1) /* name exists already */{
		    printerr("module names may not be identical");
		    continue;
		}
		else if(res == 1) {
		    printerr("impossible module size");
		    continue;
		}
		else if(res == 2) {
		    printerr("Unallowed module size");
		    continue;
		}
		/* more errors may be added */
		else /* Everything OK*/{
		    i++;
		}
	    }
	    else if(type != HMM.NOTYPE) {
		/* create forward_std and forward_alt modules */
		int res = modelmaker.createModule(name, type, HMM.EVEN, intervalStart, intervalEnd, label);
		if(res == 1) /* name exists already */{
		    printerr("module names may not be identical");
		    continue;
		}
		else if(res == 2) {
		    printerr("Unallowed module size");
		    continue;
		}
		/* more errors may be added */
		else /* Everything OK*/{
		    i++;
		}
		modelmaker.setInternalInitDistrib(name, internalInitDistrib);
	    }
	    else /* notype module */ {
		printerr("Nonexisting module- or distribution type");
		continue;
	    }
	}
				
	/* end node */
	print("End node: ");
	String e = readln();
	if(e.length() == 0) {
	    int res = modelmaker.createModule("e", HMM.ENDNODE, HMM.ZERO, 1, "0");
	}
	else {
	    int res = modelmaker.createModule(e, HMM.ENDNODE, HMM.ZERO, 1, "0");
	}
	newln();
 
    }

    private int parseType(String t) 
    {
	t = t.toLowerCase();
	if(t.equals("s")){
	    return HMM.SINGLENODE;
	}
	else if(t.equals("singlenode")) {
	    return HMM.SINGLENODE;
	}
	else if(t.equals("c")) {
	    return HMM.CLUSTER;
	}
	else if(t.equals("cluster")) {
	    return HMM.CLUSTER;
	}
	else if(t.equals("h")) {
	    return HMM.HIGHWAY;
	}
	else if(t.equals("highway")) {
	    return HMM.HIGHWAY;
	}
	else if(t.equals("u")) {
	    return HMM.U_TURN;
	}
	else if(t.equals("uturn")) {
	    return HMM.U_TURN;
	}
	else if(t.equals("l")) {
	    return HMM.SINGLELOOP;
	}
	else if(t.equals("loop")) {
	    return HMM.SINGLELOOP;
	}
	else if(t.equals("forward")) {
	    return HMM.FORWARD_STD;
	}
	else if(t.equals("f")) {
	    return HMM.FORWARD_STD;
	}
	else if(t.equals("fa")) {
	    return HMM.FORWARD_ALT;
	}
	else if(t.equals("forward_alternative")) {
	    return HMM.FORWARD_ALT;
	}
	else if(t.equals("p9")) {
	    return HMM.PROFILE9;
	}
	else if(t.equals("profile9")) {
	    return HMM.PROFILE9;
	}
		else if(t.equals("p7")) {
	    return HMM.PROFILE7;
	}
	else if(t.equals("profile7")) {
	    return HMM.PROFILE7;
	}
	/* Free to define more valid types here */

	else {
	    /* type not defined */
	    return HMM.NOTYPE;
	}
    }

    private int getSize() {
	print("Specify module size: ");
	String s = readln();
	try {
	    int size = Integer.parseInt(s);
	    if(size > 0) {
		return size;
	    }
	    else {
		printerr("Unallowed module size");
		return getSize();
	    }
	}
	catch (NumberFormatException e) {
	    printerr("Not an integer");
	    return getSize();
	}
    }
    
    private int getLength()
    {
	print("Specify expected loop length: ");
	String s = readln();
	try {
	    int size = Integer.parseInt(s);
	    if(size > 0) {
		return size;
	    }
	    else {
		printerr("Negative loop length");
		return getLength();
	    }
	}
	catch (NumberFormatException e) {
	    printerr("Not an integer");
	    return getLength();
	}
    }

    private String getIntervals()
    {
	print("Specify min and max length (min,max): ");
	String s = readln();
	return s;
    }
    
    private int getIntervalStart(String s)
    {
	StringTokenizer st = new StringTokenizer(s, ",");
	String start = "";
	if(st.hasMoreTokens()) {
	    start = st.nextToken();
	}
	else {
	    return -1;
	}
	try {
	    int i = Integer.parseInt(start);
	    if(i > 0) {
		return i;
	    }
	    else {
		printerr("Negative min value");
		return -1;
	    }
	}
	catch (NumberFormatException e) {
	    printerr("Not an integer");
	    return -1;
	}
    }
    
    private int getIntervalEnd(String s)
    {
	StringTokenizer st = new StringTokenizer(s, ",");
	String end = "";
	if(st.countTokens() == 2) {
	    st.nextToken();
	    end = st.nextToken();
	}
	else {
	    return -1;
	}
	try {
	    int i = Integer.parseInt(end);
	    if(i > 0) {
		return i;
	    }
	    else {
		printerr("Negative min value");
		return -1;
	    }
	}
	catch (NumberFormatException e) {
	    printerr("Not an integer");
	    return -1;
	}
    }

    private void getInternalInitDistrib(int moduleType, InternalInitDistrib iid)
    {
	if(moduleType == HMM.FORWARD_STD || moduleType == HMM.FORWARD_ALT) {
	    while(true) {
		print("Specify initial distribution (U,B,P): ");
		String s = readln();
		if(s.equals("U") || s.equals("Uniform")) {
		    iid.setDistrib("U");
		    break;
		}
		else if(s.equals("B") || s.equals("Binomial")) {
		    iid.setDistrib("B");
		    while(true) {
			print("Specify value of parameter p (n depends on the module length): ");
			s = readln();
			try {
			    double p = Double.parseDouble(s);
			    if(p < 0.0) {
				continue;
			    }
			    else {
				iid.setPar1(p);
				break;
			    }
			}
			catch(NumberFormatException e) {
			    continue;
			}
		    }
		    break;
		}
		else if(s.equals("P") || s.equals("Poisson")) {
		    iid.setDistrib("P");
		    while(true) {
			print("Specify value of parameter lambda and possible reversing (separate by ','): ");
			s = readln();
			StringTokenizer st = new StringTokenizer(s, ",");
			if(st.countTokens() == 0) {
			    continue;
			}
			else {
			    try {
				double p = Double.parseDouble(st.nextToken());
				if(p < 0.0) {
				    continue;
				}
				else {
				    iid.setPar1(p);
				}
			    }
			    catch(NumberFormatException e) {
				continue;
			    }
			    if(st.hasMoreTokens()) {
				String rev = st.nextToken();
				rev = rev.trim();
				if(rev.equals("Y") || rev.equals("r")) {
				    iid.setPar2(-1);
				    
				}
				else {
				    iid.setPar2(1);
				}
			    }
			    else {
				iid.setPar2(1);
			    }
			    
			    break;
			}
		    }
		    break;
		}
		else {
		    
		}
	    }
	}
	else {
	    /* not implemented yet */
	}
    }

    /********************************************************************************************
     ********************** methods for getting module interconnectivity ****************************
     ********************************************************************************************/
    private void specifyInterConnectivity()
    {
	int nrModules = modelmaker.getNrModules();
	String[] theNames = new String[nrModules];
	println("Specify interconnectivity");
	print("(Possible modules are:");
	int j = 0;
	for(ListIterator i = modelmaker.getModules();i.hasNext();) {
	    String name = ((Module)i.next()).getName();
	    print(" " + name);
	    theNames[j] = name;
	    j++;
	}
	print(")");
	newln();

	int i = 0;
	boolean nextModule;
	while(i < nrModules) {
	    nextModule = true;
	    print("Connection from " + theNames[i] +
		  " to (separate by ';'): ");
	
	    String s = readln();
	    s = s.trim();
	    StringTokenizer st = new StringTokenizer(s,";");
	    while(st.hasMoreTokens()) {
		String toModule = st.nextToken();
		toModule = toModule.trim();
		if(isModuleName(toModule, theNames)) {
		    modelmaker.setTransition(theNames[i], toModule);
		}
		else {
		    printerr(toModule + " is not a specified module");
		    nextModule = false;
		}
	    }
	    if(nextModule) {
		modelmaker.initializeTransitionProbabilities(theNames[i]);
		i++;
	    }
	}
	
    }
    
    private boolean isModuleName(String s, String[] sa)
    {
	for(int i = 0; i < sa.length; i++) {
	    if(s.equals(sa[i])) {
		return true;
	    } 
	}
	return false;
    }

    /********************************************************************************************
     ********************** methods for getting distribution groups ****************************
     ********************************************************************************************/
    private void specifyDistributionGroups()
    {
	int nrModules = modelmaker.getNrModules();
	String[] theNames = new String[nrModules];
	Hashtable addedModules = new Hashtable();
	println("Specify emisssion distribution groups");
	print("(Possible modules are:");
	int j = 0;
	for(ListIterator i = modelmaker.getModules();i.hasNext();) {
	    String name = ((Module)i.next()).getName();
	    print(" " + name);
	    theNames[j] = name;
	    j++;
	}
	print(")");
	newln();

	boolean done = false;
	boolean createDistributionGroup;
	j = 1;
	while(!done) {
	    LinkedList distribGroup = new LinkedList();
	    print("Emission distribution group " + j + " (separate by ';'):");
	    String s = readln();
	    if(s.equals("") || s.equals("d") || s.equals("done")) {
		done = true;
		continue;
	    }
	    s = s.trim();
	    StringTokenizer st = new StringTokenizer(s,";");
	    createDistributionGroup = true;
	    while(st.hasMoreTokens()) {
		String module = st.nextToken();
		module = module.trim();
		if(isModuleName(module, theNames)) {
		    if(!addedModules.containsKey(module)) {
			distribGroup.add(module);
			addedModules.put(module, module);
		    }
		    else {
			printerr("module " + module + " is already in a distribution group");
			createDistributionGroup = false;
			break;
		    }
		}
		else {
		    printerr(module + " is not a specified module");
		    createDistributionGroup = false;
		    break;
		}
	    }
	    if(createDistributionGroup) {
		modelmaker.addDistributionGroup(distribGroup);
		j++;
	    }
	    else {
		/* do not add group, let user try again */
		for(ListIterator i = (ListIterator)distribGroup.iterator();i.hasNext();) {
		    String n = (String)(i.next());
		    addedModules.remove(n);
		}
	    }
	}
    }

    
    /********************************************************************************************
     ********************** methods for getting transition ties ****************************
     ********************************************************************************************/
    private void specifyTransitionTies()
    {
	int nrModules = modelmaker.getNrModules();
	String[] theNames = new String[nrModules];
	Hashtable addedModules = new Hashtable();
	println("Specify transition-tie groups");
	print("(Possible modules are:");
	int j = 0;
	for(ListIterator i = modelmaker.getModules();i.hasNext();) {
	    String name = ((Module)i.next()).getName();
	    print(" " + name);
	    theNames[j] = name;
	    j++;
	}
	print(")");
	newln();

	boolean done = false;
	boolean createTransTieGroup;
	j = 1;
	while(!done) {
	    LinkedList transTieGroup = new LinkedList();
	    print("Transition-tie group " + j + " (separate by ';'):");
	    String s = readln();
	    if(s.equals("") || s.equals("d") || s.equals("done")) {
		done = true;
		continue;
	    }
	    s = s.trim();
	    StringTokenizer st = new StringTokenizer(s,";");
	    createTransTieGroup = true;
	    while(st.hasMoreTokens()) {
		String module = st.nextToken();
		module = module.trim();
		if(isModuleName(module, theNames)) {
		    if(!addedModules.containsKey(module)) {
			transTieGroup.add(module);
			addedModules.put(module, module);
			if(!modelmaker.identicalModules(((String)transTieGroup.get(0)), module)) {
			    printerr("module " + module + " is of incorrect type");
			    createTransTieGroup = false;
			    break;
			}
		    }
		    else {
			printerr("module " + module + " is already in a distribution group");
			createTransTieGroup = false;
			break;
		    }
		}
		else {
		    printerr(module + " is not a specified module");
		    createTransTieGroup = false;
		    break;
		}
	    }
	    if(createTransTieGroup) {
		modelmaker.addTransTieGroup(transTieGroup);
		j++;
	    }
	    else {
		/* do not add group, let user try again */
		for(ListIterator i = (ListIterator)transTieGroup.iterator();i.hasNext();) {
		    String n = (String)(i.next());
		    addedModules.remove(n);
		}
	    }
	}
    }
    
    



    /********************************************************************************************
     ********************** methods for getting initial values ****************************
     ********************************************************************************************/
    private void specifyInitvalues()
    {
	LinkedList restModules = new LinkedList();
	for(ListIterator i = modelmaker.getModules();i.hasNext();) {
	    restModules.add(i.next());
	}
	println("Specify initial emission probabilities and prior distributions for each " + 
		"module/distribution group");
	for(ListIterator i = modelmaker.getDistributionGroups();i.hasNext();) {
	    print("[ ");
	    LinkedList distribGroup = (LinkedList)i.next();
	    for(ListIterator j = (ListIterator)distribGroup.iterator();
		j.hasNext();) {
		Module m = modelmaker.getModule(((String)j.next()));
		print(m.getName() + " ");
		restModules.remove(m);
	    }
	    print("] ");
	    boolean done = false;
	    while(!done) {
		print("(initprobtype, priorfile, emisspriorscaler):");
		String s = readln();
		StringTokenizer st = new StringTokenizer(s,",");
		int initDistribType = HMM.EVEN;
		boolean locked = false;
		double[] initDistrib = null;
		if(st.countTokens() > 3) {
		    printerr("illegal choice");
		    continue;
		}
		else if(st.countTokens() == 0) {
		    s = null;
		    
		}
		else {
		    /* read init distribution instruction */
		    s = st.nextToken();
		    s = s.trim();
		    initDistribType = parseDistribType(s);
		    initDistrib = null;
		    if(initDistribType == HMM.NOTYPE) {
			continue;
		    }
		    else if(initDistribType == HMM.MANUAL || initDistribType == HMM.LOCKED_MANUAL) {
			initDistrib = getInitDistribution();
			if(initDistrib == null) {
			    continue;
			}
			if(initDistribType == HMM.LOCKED_MANUAL) {
			    locked = true;
			}
		    }
		    else if(initDistribType == HMM.RANDOM) {
			initDistrib = getRandomDistrib();
		    }
		    else if(initDistribType == HMM.LOCKED_EVEN) {
			locked = true;
		    }
		    else {
			
		    }
		}
		double priScaleValue = 1.0;
		if(st.countTokens() == 1) {
		    /* read prior file */
		    s = st.nextToken();
		    s = s.trim();
		}
		else if(st.countTokens() == 2) {
		    /* read prior file */
		    s = st.nextToken();
		    s = s.trim();

		    /* read prior scaler */
		    String priScale = st.nextToken();
		    priScale = priScale.trim();
		    priScaleValue = 1.0;
		    try {
			priScaleValue = Double.parseDouble(priScale);
			if(priScaleValue < 0.0) {
			    printerr("illegal value, setting priorscaler to 1.0");
			    priScaleValue = 1.0;
			}
			
		    }
		    catch(NumberFormatException e){
			printerr("illegal value, setting priorscaler to 1.0");
		    }
		}
		else if(st.countTokens() == 0) {
		    s = null;
		}
		/* give info to modules in distribution group */
		for(ListIterator k = (ListIterator)distribGroup.iterator();
		    k.hasNext();) {
		    Module m = modelmaker.getModule((String) k.next());
		    m.setPriorfile(s);
		    m.setEmissPriorScaler(priScaleValue);
		    modelmaker.addPriorfile(s);
		    if(initDistribType == HMM.MANUAL || initDistribType == HMM.RANDOM) {
			m.setDistribType(initDistribType, initDistrib);
		
		    }
		    else {
			m.setDistribType(initDistribType);
		    }
		    if(locked) {
			m.lockVertexEmissions();
		    }
		}
		done = true;
	    }
	}
	for(ListIterator i = (ListIterator)restModules.iterator();i.hasNext();) {
	    Module m = (Module)i.next();
	    if(m.getVertexType() == HMM.START || m.getVertexType() == HMM.END) {
		continue;
	    }
	    print(m.getName() + " ");
	    boolean done = false;
	    while(!done) {
		print("(initprobtype, priorfile, emisspriorscaler):");
		String s = readln();
		StringTokenizer st = new StringTokenizer(s,",");
		int initDistribType = HMM.EVEN;
		boolean locked = false;
		double[] initDistrib = null;
		if(st.countTokens() > 3 ) {
		    printerr("illegal choice");
		    continue;
		}
		else if(st.countTokens() == 0) {
		    s = null;
		    
		}
		else {
		    /* read init distribution instruction */
		    s = st.nextToken();
		    s = s.trim();
		    initDistribType = parseDistribType(s);
		    initDistrib = null;
		    if(initDistribType == HMM.NOTYPE) {
			continue;
		    }
		    else if(initDistribType == HMM.MANUAL || initDistribType == HMM.LOCKED_MANUAL) {
			initDistrib = getInitDistribution();
			if(initDistrib == null) {
			    continue;
			}
			if(initDistribType == HMM.LOCKED_MANUAL) {
			    locked = true;
			}
		    }
		    else if(initDistribType == HMM.RANDOM) {
			initDistrib = getRandomDistrib();
		    }
		    else if(initDistribType == HMM.LOCKED_EVEN) {
			locked = true;
		    }
		    else {
			
		    }
		}   
		double priScaleValue = 1.0;
		if(st.countTokens() == 1) {
		    /* read prior file */
		    s = st.nextToken();
		    s = s.trim();
		}
		else if(st.countTokens() == 2) {
		    /* read prior file */
		    s = st.nextToken();
		    s = s.trim();

		    /* read prior scaler */
		    String priScale = st.nextToken();
		    priScale = priScale.trim();
		    priScaleValue = 1.0;
		    try {
			priScaleValue = Double.parseDouble(priScale);
			if(priScaleValue < 0.0) {
			    printerr("illegal value, setting priorscaler to 1.0");
			    priScaleValue = 1.0;
			}
			
		    }
		    catch(NumberFormatException e){
			printerr("illegal value, setting priorscaler to 1.0");
		    }
		}
		else if(st.countTokens() == 0) {
		    s = null;
		}
		/* give info to module */
		m.setPriorfile(s);
		m.setEmissPriorScaler(priScaleValue);
		modelmaker.addPriorfile(s);
		if(initDistribType == HMM.MANUAL) {
		    m.setDistribType(initDistribType, initDistrib);
		}
		else {
		    m.setDistribType(initDistribType);
		}
		if(locked) {
		    m.lockVertexEmissions();
		}
		done = true;
		
	    }
	}
    }

    private void specifyInitvalues(int nr)
    {
	LinkedList restModules = new LinkedList();
	for(ListIterator i = modelmaker.getModules();i.hasNext();) {
	    restModules.add(i.next());
	}
	println("Specify initial emission probabilities and prior distributions for each " + 
		"module/distribution group");
	for(ListIterator i = modelmaker.getDistributionGroups();i.hasNext();) {
	    print("[ ");
	    LinkedList distribGroup = (LinkedList)i.next();
	    for(ListIterator j = (ListIterator)distribGroup.iterator();
		j.hasNext();) {
		Module m = modelmaker.getModule(((String)j.next()));
		print(m.getName() + " ");
		restModules.remove(m);
	    }
	    print("] ");
	    boolean done = false;
	    while(!done) {
		print(" (alphabet nr " + nr + ") (initprobtype, priorfile, emisspriorscaler):");
		String s = readln();
		StringTokenizer st = new StringTokenizer(s,",");
		int initDistribType = HMM.EVEN;
		boolean locked = false;
		double[] initDistrib = null;
		if(st.countTokens() > 3) {
		    printerr("illegal choice");
		    continue;
		}
		else if(st.countTokens() == 0) {
		    s = null;
		    
		}
		else {
		    /* read init distribution instruction */
		    s = st.nextToken();
		    s = s.trim();
		    initDistribType = parseDistribType(s);
		    initDistrib = null;
		    if(initDistribType == HMM.NOTYPE) {
			continue;
		    }
		    else if(initDistribType == HMM.MANUAL || initDistribType == HMM.LOCKED_MANUAL) {
			initDistrib = getInitDistribution();
			if(initDistrib == null) {
			    continue;
			}
			if(initDistribType == HMM.LOCKED_MANUAL) {
			    locked = true;
			}
		    }
		    else if(initDistribType == HMM.RANDOM) {
			initDistrib = getRandomDistrib(nr);
		    }
		    else if(initDistribType == HMM.LOCKED_EVEN) {
			locked = true;
		    }
		    else {
			
		    }
		}
		double priScaleValue = 1.0;
		if(st.countTokens() == 1) {
		    /* read prior file */
		    s = st.nextToken();
		    s = s.trim();
		}
		else if(st.countTokens() == 2) {
		    /* read prior file */
		    s = st.nextToken();
		    s = s.trim();

		    /* read prior scaler */
		    String priScale = st.nextToken();
		    priScale = priScale.trim();
		    priScaleValue = 1.0;
		    try {
			priScaleValue = Double.parseDouble(priScale);
			if(priScaleValue < 0.0) {
			    printerr("illegal value, setting priorscaler to 1.0");
			    priScaleValue = 1.0;
			}
			
		    }
		    catch(NumberFormatException e){
			printerr("illegal value, setting priorscaler to 1.0");
		    }
		}
		else if(st.countTokens() == 0) {
		    s = null;
		}
		/* give info to modules in distribution group */
		for(ListIterator k = (ListIterator)distribGroup.iterator();
		    k.hasNext();) {
		    Module m = modelmaker.getModule((String) k.next());
		    m.setPriorfile(nr,s);
		    m.setEmissPriorScaler(nr,priScaleValue);
		    modelmaker.addPriorfile(nr,s);
		    if(initDistribType == HMM.MANUAL || initDistribType == HMM.RANDOM) {
			m.setDistribType(nr, initDistribType, initDistrib);
		    }
		    else {
			m.setDistribType(nr, initDistribType);
		    }
		    if(locked) {
			m.lockVertexEmissions();
		    }
		}
		done = true;
	    }
	}
	for(ListIterator i = (ListIterator)restModules.iterator();i.hasNext();) {
	    Module m = (Module)i.next();
	    if(m.getVertexType() == HMM.START || m.getVertexType() == HMM.END) {
		continue;
	    }
	    print(m.getName() + " ");
	    boolean done = false;
	    while(!done) {
		print(" (alphabet nr " + nr + ") (initprobtype, priorfile, emisspriorscaler):");
		String s = readln();
		StringTokenizer st = new StringTokenizer(s,",");
		int initDistribType = HMM.EVEN;
		boolean locked = false;
		double[] initDistrib = null;
		if(st.countTokens() > 3 ) {
		    printerr("illegal choice");
		    continue;
		}
		else if(st.countTokens() == 0) {
		    s = null;
		    
		}
		else {
		    /* read init distribution instruction */
		    s = st.nextToken();
		    s = s.trim();
		    initDistribType = parseDistribType(s);
		    initDistrib = null;
		    if(initDistribType == HMM.NOTYPE) {
			continue;
		    }
		    else if(initDistribType == HMM.MANUAL || initDistribType == HMM.LOCKED_MANUAL) {
			initDistrib = getInitDistribution(nr);
			if(initDistrib == null) {
			    continue;
			}
			if(initDistribType == HMM.LOCKED_MANUAL) {
			    locked = true;
			}
		    }
		    else if(initDistribType == HMM.RANDOM) {
			initDistrib = getRandomDistrib(nr);
		    }
		    else if(initDistribType == HMM.LOCKED_EVEN) {
			locked = true;
		    }
		    else {
			
		    }
		}   
		double priScaleValue = 1.0;
		if(st.countTokens() == 1) {
		    /* read prior file */
		    s = st.nextToken();
		    s = s.trim();
		}
		else if(st.countTokens() == 2) {
		    /* read prior file */
		    s = st.nextToken();
		    s = s.trim();

		    /* read prior scaler */
		    String priScale = st.nextToken();
		    priScale = priScale.trim();
		    priScaleValue = 1.0;
		    try {
			priScaleValue = Double.parseDouble(priScale);
			if(priScaleValue < 0.0) {
			    printerr("illegal value, setting priorscaler to 1.0");
			    priScaleValue = 1.0;
			}
			
		    }
		    catch(NumberFormatException e){
			printerr("illegal value, setting priorscaler to 1.0");
		    }
		}
		else if(st.countTokens() == 0) {
		    s = null;
		}
		/* give info to module */
		m.setPriorfile(nr,s);
		m.setEmissPriorScaler(nr,priScaleValue);
		modelmaker.addPriorfile(nr,s);
		if(initDistribType == HMM.MANUAL) {
		    m.setDistribType(nr,initDistribType, initDistrib);
		}
		else {
		    m.setDistribType(nr,initDistribType);
		}
		if(locked) {
		    m.lockVertexEmissions();
		}
		done = true;
		
	    }
	}
    }

    
    /********************************************************************************************
     ********************** methods for getting initial values ****************************
     ********************************************************************************************/
    private void specifyTransitionInitvalues()
    {
	LinkedList restModules = new LinkedList();
	for(ListIterator i = modelmaker.getModules();i.hasNext();) {
	    restModules.add(i.next());
	}
	println("Specify initial transition probabilities and prior distributions for each " + 
		"module/transition distribution group");
	for(ListIterator i = modelmaker.getTransTieGroups();i.hasNext();) {
	    print("[ ");
	    LinkedList distribGroup = (LinkedList)i.next();
	    for(ListIterator j = (ListIterator)distribGroup.iterator();
		j.hasNext();) {
		Module m = modelmaker.getModule(((String)j.next()));
		print(m.getName() + " ");
		restModules.remove(m);
	    }
	    print("] ");
	    boolean done = false;
	    while(!done) {
		print("(priorfile, transpriorscaler):");
		String s = readln();
		StringTokenizer st = new StringTokenizer(s,",");
		int initDistribType = HMM.EVEN;
		double[] initDistrib = null;
		double priScaleValue = 1.0;
		if(st.countTokens() > 2) {
		    printerr("illegal choice, try again");
		    continue;
		}
		else if(st.countTokens() == 0) {
		    s = null;
		    
		}
		else {
		    
		    if(st.countTokens() == 1) {
			/* read prior file */
			s = st.nextToken();
			s = s.trim();
		    }
		    else if(st.countTokens() == 2) {
			/* read prior file */
			s = st.nextToken();
			s = s.trim();
			
			/* read prior scaler */
			String priScale = st.nextToken();
			priScale = priScale.trim();
			priScaleValue = 1.0;
			try {
			    priScaleValue = Double.parseDouble(priScale);
			    if(priScaleValue < 0.0) {
				printerr("illegal value, setting priorscaler to 1.0");
				priScaleValue = 1.0;
			    }
			    
			}
			catch(NumberFormatException e){
			    printerr("illegal value, setting priorscaler to 1.0");
			}
		    }
		}
		/* give info to modules in distribution group */
		for(ListIterator k = (ListIterator)distribGroup.iterator(); k.hasNext();) {
		    Module m = modelmaker.getModule((String) k.next());
		    m.setTransPriorfile(s);
		    m.setTransPriorScaler(priScaleValue);
		    modelmaker.addTransPriorfile(s);
		}
		done = true;
	    }
	}
	for(ListIterator i = (ListIterator)restModules.iterator();i.hasNext();) {
	    Module m = (Module)i.next();
	    if(m.getVertexType() == HMM.START || m.getVertexType() == HMM.END) {
		continue;
	    }
	    print(m.getName() + " ");
	    boolean done = false;
	    while(!done) {
		print("(priorfile, transpriorscaler):");
		String s = readln();
		StringTokenizer st = new StringTokenizer(s,",");
		int initDistribType = HMM.EVEN;
		double[] initDistrib = null;
		double priScaleValue = 1.0;
		if(st.countTokens() > 3 ) {
		    printerr("illegal choice");
		    continue;
		}
		else if(st.countTokens() == 0) {
		    s = null;
		    
		}
		else {
		    if(st.countTokens() == 1) {
			/* read prior file */
			s = st.nextToken();
			s = s.trim();
		    }
		    else if(st.countTokens() == 2) {
			/* read prior file */
			s = st.nextToken();
			s = s.trim();
			
			/* read prior scaler */
			String priScale = st.nextToken();
			priScale = priScale.trim();
			priScaleValue = 1.0;
			try {
			    priScaleValue = Double.parseDouble(priScale);
			    if(priScaleValue < 0.0) {
				printerr("illegal value, setting priorscaler to 1.0");
				priScaleValue = 1.0;
			    }
			    
			}
			catch(NumberFormatException e){
			    printerr("illegal value, setting priorscaler to 1.0");
			}
		    }
		}
		/* give info to module */
		m.setTransPriorfile(s);
		m.setTransPriorScaler(priScaleValue);
		modelmaker.addTransPriorfile(s);
		done = true;
		
	    }
	}
	
	/* mojliggor specificering av transitionssannolikheter inom modul,
	   och fran modul, specifikt for varje modul */
    }
    
    private int parseDistribType(String t) 
    {
	t = t.toLowerCase();
	if(t.equals("u")){
	    return HMM.EVEN;
	}
	else if(t.equals("uniform")) {
	    return HMM.EVEN;
	}
	else if(t.equals("z")) {
	    return HMM.ZERO;
	}
	else if(t.equals("zero")) {
	    return HMM.ZERO;
	}
	else if(t.equals("r")) {
	    return HMM.RANDOM;
	}
	else if(t.equals("random")) {
	    return HMM.RANDOM;
	}
	else if(t.equals("m")) {
	    return HMM.MANUAL;
	}
	else if(t.equals("man")) {
	    return HMM.MANUAL;
	}
	else if(t.equals("lm")) {
	    return HMM.LOCKED_MANUAL;
	}
	else if(t.equals("locked_man")) {
	    return HMM.LOCKED_MANUAL;
	}
	else if(t.equals("lu")) {
	    return HMM.LOCKED_EVEN;
	}
	else if(t.equals("locked_uniform")) {
	    return HMM.LOCKED_EVEN;
	}
	/* Free to define more valid types here */

	else {
	    /* type not defined */
	    return HMM.NOTYPE;
	}
    }

    private double[] getInitDistribution()
    {
	print("Specify initial emission probabilities file ('k' to write on keyboard): ");
	String s = readln();
	if(s.equals("k")) {
	    print("Enter probabilities (separate by blank): ");
	    s = readln();
	    return makeProbArray(s);
	}
	else {
	    /* read probs from file */
	    try {
		BufferedReader probReader = new BufferedReader(new FileReader(s));
		String probs = probReader.readLine();
		probs = probs.trim();
		return makeProbArray(probs);
	    }
	    catch(IOException e) {
		P.MESSAGE("I/O error: Couldn't read from file '" + s + "'");
		return null;
	    }
	}
    }

    private double[] getInitDistribution(int nr)
    {
	print("Specify initial emission probabilities file ('k' to write on keyboard): ");
	String s = readln();
	if(s.equals("k")) {
	    print("Enter probabilities (separate by blank): ");
	    s = readln();
	    return makeProbArray(nr, s);
	}
	else {
	    /* read probs from file */
	    try {
		BufferedReader probReader = new BufferedReader(new FileReader(s));
		String probs = probReader.readLine();
		probs = probs.trim();
		return makeProbArray(nr, probs);
	    }
	    catch(IOException e) {
		P.MESSAGE("I/O error: Couldn't read from file '" + s + "'");
		return null;
	    }
	}
    }

    private double[] getRandomDistrib()
    {
	double[] probs = new double[modelmaker.getAlphabetSize()];
	int nrEmissions = modelmaker.getAlphabetSize();
	double sum = 0;
	for(int i = 0; i < nrEmissions; i++) {
	    double initProb = Math.random();
	    probs[i] = initProb;
	    sum = sum + initProb;
	}
	for(int i = 0; i < nrEmissions; i++) {
	    probs[i] = probs[i] / sum;
	}
	return probs;
    }

    private double[] getRandomDistrib(int nr)
    {
	double[] probs = new double[modelmaker.getAlphabetSize(nr)];
	int nrEmissions = modelmaker.getAlphabetSize(nr);
	double sum = 0;
	for(int i = 0; i < nrEmissions; i++) {
	    double initProb = Math.random();
	    probs[i] = initProb;
	    sum = sum + initProb;
	}
	for(int i = 0; i < nrEmissions; i++) {
	    probs[i] = probs[i] / sum;
	}
	return probs;
    }
    
    private double[] makeProbArray(String s)
    {
	StringTokenizer st = new StringTokenizer(s, " ");
	double[] probs = new double[modelmaker.getAlphabetSize()];
	if(st.countTokens() != modelmaker.getAlphabetSize()) {
	    P.MESSAGE("Distribution has incorrect format for this alphabet");
	    return null;
	}
	else {
	    try {
		double sum = 0;
		for(int i = 0; i < modelmaker.getAlphabetSize();i++) {
		    String p = st.nextToken();
		    double prob = Double.parseDouble(p);
		    sum = sum + prob;
		    probs[i] = prob;
		}
		if(sum != 1.0) {
		    P.MESSAGE("Warning: sum of probabilities not equal" +
			      " to 1.0: autocorrecting");
		    for(int i = 0; i < modelmaker.getAlphabetSize();i++) {
			probs[i] = probs[i] / sum;
		    }
		}
	    }
	    catch(NumberFormatException e) {
		P.MESSAGE("Some distribution object is not a number");
		return null;
	    }
	}
	return probs;
    }

    private double[] makeProbArray(int nr, String s)
    {
	StringTokenizer st = new StringTokenizer(s, " ");
	double[] probs = new double[modelmaker.getAlphabetSize(nr)];
	if(st.countTokens() != modelmaker.getAlphabetSize(nr)) {
	    P.MESSAGE("Distribution has incorrect format for this alphabet");
	    return null;
	}
	else {
	    try {
		double sum = 0;
		for(int i = 0; i < modelmaker.getAlphabetSize(nr);i++) {
		    String p = st.nextToken();
		    double prob = Double.parseDouble(p);
		    sum = sum + prob;
		    probs[i] = prob;
		}
		if(sum != 1.0) {
		    P.MESSAGE("Warning: sum of probabilities not equal" +
			      " to 1.0: autocorrecting");
		    for(int i = 0; i < modelmaker.getAlphabetSize(nr);i++) {
			probs[i] = probs[i] / sum;
		    }
		}
	    }
	    catch(NumberFormatException e) {
		P.MESSAGE("Some distribution object is not a number");
		return null;
	    }
	}
	return probs;
    }










    
    private void cleanUp()
    {
	
    }


    /********************************************************************************************
     ********************** method for saving HMM ****************************
     ********************************************************************************************/

    private void saveHMM(String outdir)
    {
	int res = modelmaker.saveHMM(outdir);
    }

}
