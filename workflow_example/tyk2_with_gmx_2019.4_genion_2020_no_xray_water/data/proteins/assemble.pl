#!/usr/bin/perl

use strict;
use warnings;
use Cwd;

my $case = "tyk2";

my $ff="amber";

my $base = getcwd;
my $toppath = "/home/ykhalak/Projects/schroedinger_set/absolute/$case/top_gaff2_sigmahole";
my $ligpath = "/home/ykhalak/Projects/schroedinger_set/absolute/$case/ligands";
my $protpath = "/home/ykhalak/Projects/schroedinger_set/absolute/$case/protein/$ff";


open(FOO,"$ligpath/ligands.txt");
my @cont = <FOO>;
close(FOO);


system("mkdir -p $base/$case > /dev/null 2>&1");
system("grep ATOM $protpath/prot.pdb > $base/$case/prot.pdb");
system("cp $protpath/topol.itp $base/$case/prot.itp");


foreach my $lig(@cont)
{
    chomp($lig);
    print("$lig\n");

    my $folder = "lig_$lig";
    
    system("grep ATOM $toppath/lig_$lig/mol_gmx.pdb > $base/$case/temp_$lig.pdb");
    system("cat $base/$case/prot.pdb $base/$case/temp_$lig.pdb > $base/$case/prot_$lig.pdb");
    system("rm -r $base/$case/temp_$lig.pdb");
    
}
print "Done.\n";
