#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Math::Complex;

my $add   = "P";
my $core; 

my $set  = 1;
my $msh  = "4";
my $type = "60";
my $damp = "16000";
my $pos  = 0.3;
my $str  = "shell_young_modulus";
my $young = 1e9;	
my $pres = 0;
my $thick = 1e-2 / $type;

my $sp    = 0.1;
my $bulk  = 4 * $young / ($type**2 * sqrt(12*(1-$pos**2))*0.10);
my $vol   = 0.3;
my $ir    = 0.0085;
my $stop  = 0;


if ($stop == 1) 
{ 
	$sp = 1;
    $add = "R";
    $core = "true";
} else {
	$core = "false";
}

for (my $a = 0; $a < 1; $a = $a + 1) 
{
	my $ort = 1;
	$young *= $ort;	
	my $shear = (3846 - (1 - $ort)*2e3)*1e5;

	if ($set == 1) {$pos *= $ort;} 
	else { $str .= "2"; }

	my $name;
	if ($stop == 1)
	{
		$name = "data$set\_$msh\_T$type\_D$damp\_O$ort\_V$vol\_$add$ir";
	} else {
		$name = "data$set\_$msh\_T$type\_D$damp\_O$ort\_V$vol\_$add$sp";
	}

	system("OMP_NUM_THREADS=3 OMP_PROC_BIND=true ./shell-opt -c config_initial.cfg -randomness 0.01 -$str $young ".
			"-shell_poisson_ratio $pos -shell_shear_modulus $shear -bulk_modulus $bulk -thickness $thick ".
			"-Rx $ir -Ry $ir -Rz $ir -damping_alpha $damp -pressure -$pres -volume_growth $vol -stop_pressure ".
			"$sp -shell_core_contact $core -m icosasphere$msh.off -d $name");
}
