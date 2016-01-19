class Emission
{


    public char letter;
    public double probability;

    public Emission(double probability, char letter)
    {
	this.letter = letter;
	this.probability = probability;
    }

    public double getProbability()
    {
	return probability;
    }

    public char getLetter()
    {
	return letter;
    }

    public void setProbability(double newProbability)
    {
	probability = newProbability;
    }
}
