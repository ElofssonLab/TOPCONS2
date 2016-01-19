class P
{
    public static void DEBUG(String s)
    {
	System.out.println(s);
    }

    public static void INTERNAL_ERROR(String s)
    {
	System.out.println(s);
	System.exit(0);
    }
    
    public static void MESSAGE(String s)
    {
	System.out.println(s);
    }
}
