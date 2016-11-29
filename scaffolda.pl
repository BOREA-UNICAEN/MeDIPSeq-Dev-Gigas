#!/usr/bin/perl

unless ($ARGV[1]) {
   print "Uso:  $0  \"/percorso/file(s)/input\"  /file/di/filtro\nVirgolettare il primo percorso se contiene caratteri speciali.\n";
   exit
   }

$dir=$ARGV[0];
@inputs = split /\n/, `ls -1 $dir`;
if ($#inputs==-1)  {
   print "$dir non contiene nessun file o non è un percorso valido.\nVirgolettare il primo percorso se contiene caratteri speciali.\n";
   exit
   }

unless (-e $ARGV[1])  {
   print "$ARGV[1] non esiste o non è un file valido.\n";
   exit
   }

@patterns=();
open FILTRO, $ARGV[1];
while (<FILTRO>)
   {
   /^(.+)\s+.+$/;
   push @patterns, $1;
   }
close FILTRO;

foreach $in (@inputs)
   {
   print "$in -> X".$in."\n";
   open IN, $in;
   open OUT, "> X".$in;
   while (<IN>)
      {
      foreach $pat (@patterns)
         {
         print OUT if /$pat/;
         }
      }
   close OUT;
   close IN;
   }
