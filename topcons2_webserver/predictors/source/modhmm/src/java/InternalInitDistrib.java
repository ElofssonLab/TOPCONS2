class InternalInitDistrib
{
    /* possible distributions :
       U = uniform, uses 0 parameters
       P = poisson distributed, uses par1 for lambda and par2 for reverseness (negative value means reverse)
       B = binomial distributed, uses par1 for p

    */

    double par1;
    double par2;
    double par3;
    String distrib;
    
    public InternalInitDistrib()
    {
	distrib = "U";
	par1 = 1.0;
	par2 = 1.0;
	par3 = 1.0;
    }

    public void setDistrib(String d)
    {
	distrib = d;
    }
    
    public void setPar1(double p)
    {
	par1 = p;
    }
    
    public void setPar2(double p)
    {
	par2 = p;
    }
    
    public void setPar3(double p)
    {
	par2 = p;
    }

    public String getDistrib()
    {
	return distrib;
    }

    public double getPar1()
    {
	return par1;
    }
    
    public double getPar2()
    {
	return par1;
    }
    
    public double getPar3()
    {
	return par1;
    }
    

    
}
