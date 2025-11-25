#!/usr/bin/perl
use strict;
srand();

for(my $frame=0; $frame<10000; $frame++)
{
    print ">$frame\n";
    for(my $atom=0; $atom<250; $atom++)
    {
        my $x = rand(100);
        my $y = rand(100);
        my $z = rand(100);
        printf("%.3f %.3f %.3f\n", $x, $y, $z);
    }
}
