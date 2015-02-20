#!/usr/bin/perl -w

if($#ARGV != 1) {
    printf "Usage: mod_upper_caser.pl <protnamefile> <modfiledir>\n";
    exit;
}

my $protnamefile = $ARGV[0];
my $modfiledir = $ARGV[1]."/";


open(PROTNAMEFILE, "$protnamefile")
    or die "could not open $protnamefile\n";

while(<PROTNAMEFILE>) {
    chomp;
    my $protname = "$_";
    my $modfile = "$modfiledir"."$protname".".mod";
    
    open(MODFILE, "$modfile")
	or die "could not open $modfile\n";
    
    my $mod = "";
    while(<MODFILE>) {
	if($_ =~ '/') {
	}
	else {
	    $_ =~ tr/a-z/A-Z/;
	    $mod .= $_;
	    
	}
    }
  
    close MODFILE;
    
    open(MODFILE, ">"."$modfile")
	or die "could not open $modfile\n";
    print MODFILE $mod;
    close MODFILE;

    print "Done upper_casing ${modfile}\n";
}

close(PROTNAMEFILE);
