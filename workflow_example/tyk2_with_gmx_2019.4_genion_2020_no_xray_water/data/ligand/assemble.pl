#!/usr/bin/perl

use strict;
use warnings;
use Cwd;

my $case = "tyk2";

my $base = getcwd;
my $toppath = "/home/ykhalak/Projects/schroedinger_set/absolute/$case/top_gaff2_sigmahole";
my $ligpath = "/home/ykhalak/Projects/schroedinger_set/absolute/$case/ligands";


open(FOO,"$ligpath/ligands.txt");
my @cont = <FOO>;
close(FOO);

foreach my $lig(@cont)
{
    chomp($lig);
    print("$lig\n");

    #my $folder = "lig_$lig";
    
    system("mkdir -p $base/$lig > /dev/null 2>&1");

    system("cp $toppath/lig_$lig/mol_gmx.pdb $base/$lig/ligand.pdb");
    system("cat $toppath/lig_$lig/ffMOL.itp > $base/$lig/lig.itp");
    system("cat $toppath/lig_$lig/MOL.itp >> $base/$lig/lig.itp");

}
print "Done.\n";
