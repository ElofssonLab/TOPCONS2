import java.io.*;
import java.util.*;

abstract class Module
{
    int[] inVertices;
    int[] outVertices;
    String priorfile;

    abstract String getName();
    abstract boolean addTransition(int fromVertex, int toVertex);
    abstract boolean addTransition(int toVertex);
    abstract int[] getInVertices();
    abstract int[] getOutVertices();
    abstract int getNrOfTransitions();
    abstract int getNrOfRegularTransitions();
    abstract int getNrOfEndTransitions();
    abstract void initializeTransitionProbabilities();
    abstract boolean write(int nrOfAlphabets, BufferedWriter writer);
    abstract int getVertexType();
    abstract boolean addEndTransition(int toVertex);
    abstract LinkedList getVertices();
    abstract int getDistribType();
    abstract double[] getEmissionProbs();
    abstract void setDistribType(int distribType, double[] distribution);
    abstract void setDistribType(int distribType);
    abstract void setDistribType(int nr, int distribType, double[] distribution);
    abstract void setDistribType(int nr, int distribType);
    //abstract void setTransDistribType(int transDistribType, double[] distribution);
    //abstract void setTransDistribType(int transDistribType);
    abstract void setPriorfile(String s);
    abstract void setPriorfile(int nr, String s);
    abstract void setTransPriorfile(String s);
    abstract String getPriorfile();
    abstract int getSize();
    abstract void setEmissPriorScaler(double d);
    abstract void setEmissPriorScaler(int nr, double d);
    abstract void setTransPriorScaler(double d);
    abstract void setInternalInitDistrib(InternalInitDistrib iid);
    abstract void lockVertexEmissions();
    abstract Vertex getVertex(int nr);
    abstract void addVerticesToVertexHash(Hashtable h);
}
