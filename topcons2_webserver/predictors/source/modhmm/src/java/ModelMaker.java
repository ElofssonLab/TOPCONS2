/* This class is an interface between the input data and the actual HMM-model.
   The purpose of this is to simplify the addition of a grafical user interface. */

import java.util.*;
 
class ModelMaker
{
    
    HMM theHMM;
    HMMSaver hmmsaver;
    
    public ModelMaker()
    {
       
    }
    /*******************get methods***************************/
    public int getAlphabetSize()
    {
	return theHMM.getAlphabetSize();
    }
    
    public int getAlphabetSize(int nr)
    {
	switch(nr) {
	case 1:
	    return theHMM.getAlphabetSize();
	case 2:
	    return theHMM.getAlphabetSize2();
	case 3:
	    return theHMM.getAlphabetSize3();
	case 4:
	    return theHMM.getAlphabetSize4();
	}
	return 0;
    }

    public int getNrVertices()
    {
        return theHMM.getNrVertices();
    }

    public int getNrModules()
    {
	return theHMM.getNrModules();
    }
    
    public String getName()
    {
	return theHMM.getName();
    }

    public String getFileName()
    {
	return theHMM.getFileName();
    }

    public Vertex getVertex(int nr)
    {
	return theHMM.getVertex(nr);
    }
    
    public Enumeration getModuleKeys()
    {
	return theHMM.getModuleKeys();
    }
    
    public ListIterator getModules()
    {
	return theHMM.getModules();
    }

    public Module getModule(String name)
    {
        return theHMM.getModule(name);
    }
    
    public ListIterator getDistributionGroups()
    {
        return theHMM.getDistributionGroups();
    }

    public ListIterator getTransTieGroups()
    {
        return theHMM.getTransTieGroups();
    }

    public boolean transitionExists(int fromNr, int toNr)
    {
	return theHMM.transitionExists(fromNr, toNr);
    }

    /************create-methods*****************************/
    public void createHMM(String name)
    {
	theHMM = new HMM(name);
	hmmsaver = new HMMSaver(theHMM);
	HMM.resetStatics();
    }

    public int createModule(String name, int type, int distribType, int size, String label)
    {
	return theHMM.createModule(name, type, distribType, size, label, false);
    }
     
    public int createModule(String name, int type, int distribType, int intervalStart, int intervalEnd, String label)
    {
	return theHMM.createModule(name, type, distribType, intervalStart, intervalEnd, label);
    }

    public int createModule(String name, int type, int distribType, int size, String label, boolean global)
    {
	return theHMM.createModule(name, type, distribType, size, label, global);
    }
    
    /* not needed for now
    public int createModule(String name, int type, double[] initDistrib, int size)
    {
	return theHMM.createModule(name, type, initDistrib, size);
    }
    */


    /*******************set methods ********************************/
    public void setNrOfAlphabets(int nr)
    {
	theHMM.setNrOfAlphabets(nr);
    }
    
    public void setAlphabet(int alphabet)
    {
	theHMM.setAlphabet(alphabet);
    }

    public void setAlphabet(String[] alphabet)
    {
	theHMM.setAlphabet(alphabet);
    }

    public void setAlphabet(int nr, int alphabet)
    {
	theHMM.setAlphabet(nr, alphabet);
    }

    public void setAlphabet(int nr, String[] alphabet)
    {
	theHMM.setAlphabet(nr, alphabet);
    }
    public void setFileName(String s)
    {
	theHMM.setFileName(s);
    }

    public void setTransition(String fromModule, String toModule)
    {
	theHMM.setTransition(fromModule, toModule);
    }
    
    public void setInternalInitDistrib(String name, InternalInitDistrib iid)
    {
	theHMM.setInternalInitDistrib(name, iid);
    }

    /*******************add methods ****************************/
    public void addPriorfile(String s)
    {
	theHMM.addPriorfile(s);
    }

    public void addPriorfile(int nr, String s)
    {
	theHMM.addPriorfile(nr, s);
    }

    public void addTransPriorfile(String s)
    {
	theHMM.addTransPriorfile(s);
    }

    public void addDistributionGroup(LinkedList group)
    {
	theHMM.addDistributionGroup(group);
    }

    public void addTransTieGroup(LinkedList transTieGroup)
    {
	theHMM.addTransTieGroup(transTieGroup);
    }

    /*******************misc methods ***************************/   
    public boolean identicalModules(String moduleA, String moduleB)
    {
	return theHMM.identicalModules(moduleA, moduleB);
    }

    public void initializeTransitionProbabilities(String moduleName)
    {
	theHMM.initializeTransitionProbabilities(moduleName);
    }

    public int saveHMM(String outdir)
    {
	return hmmsaver.saveHMM(outdir);
    }
    
}
