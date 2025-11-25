#!/usr/bin/perl
use strict;
my $maxrand  = 50;
my $maxframe = 10000;
my $maxatom  = 250;
srand();
for(my $frame=0; $frame<$maxframe; $frame++)
{
    print ">$frame\n";
    for(my $atom=0; $atom<$maxatom; $atom++)
    {
        my $x = rand($maxrand);
        my $y = rand($maxrand);
        my $z = rand($maxrand);
        
        printf("%.3f %.3f %.3f\n", $x, $y, $z);
    }
}
