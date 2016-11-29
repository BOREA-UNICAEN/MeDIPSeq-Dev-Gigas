#!/usr/bin/perl

$INPUT = "AllmRNAs.txt";
$OUTPUT = "GOmappings.txt";

open IN, $INPUT;
open OUT, ">$OUTPUT";
while (<IN>)
   {
   if (/(CGI_\d+).+\"(.+)\"/)
      {
      $go=$2;
      print OUT "$1\t";
      @GOS = ($go =~ m/GO:(\d+)/g);
      print OUT "GO:", join(", GO:", @GOS), "\n";
      }
   else { print OUT $_; }
   }
close IN;
close OUT;
