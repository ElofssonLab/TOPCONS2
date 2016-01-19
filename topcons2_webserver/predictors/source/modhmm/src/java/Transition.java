class Transition
{
    public double probability;
    public int toVertex;
    public boolean isLocked;
    public Transition(int toVertex, double probability)
    {
	this.toVertex = toVertex;
	this.probability = probability;
	isLocked = false;
    }

    public Transition(int toVertex, double probability, boolean isLocked)
    {
	this.toVertex = toVertex;
	this.probability = probability;
	this.isLocked = isLocked;
    }

    public double getProbability()
    {
	return probability;
    }

    public int getToVertex()
    {
	return toVertex;
    }

    public void setProbability(double newProb)
    {
	probability = newProb;
    }
    
    public boolean isLocked()
    {
	return isLocked;
    }

    public void lock()
    {
	isLocked = true;
    }

    public void unLock()
    {
	isLocked = false;
    }
}
