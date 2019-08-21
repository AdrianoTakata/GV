#!/usr/bin/perl -w

use strict;

my $filein  = 'par.for';
my $fileout = 'constants';
my $fileres = 'pre_processing/constants.f90';

my %map;

open(FD, "$filein");
while(<FD>) {
	if(m/^c/) {
		next;
	}
	if(m/my_form = (\d+)/) {
		$map{"myform"} = $1;
	}
	if(m/U_1 = ([-0-9.de]+)/) {
		$map{"u1"} = $1;
	}
	if(m/L_1 = ([-0-9.de]+)/) {
		$map{"l1"} = $1;
	}
	if(m/N_1 = ([-0-9.de\/]+)/) {
		$map{"n1"} = $1;
	}
	if(m/Re = ([-0-9.den*_LUN\/]+)/) {
		$map{"re"} = $1;
	}
	if(m/fac_y = ([-0-9.dRe]+)/) {
		$map{"fac_y"} = $1;
	}
	if(m/imax = (\d+)/) {
		$map{"imax"} = $1;
	}
	if(m/dx = ([-0-9.de]+)/) {
		$map{"dx"} = $1;
	}
	if(m/jmax = (\d+)/) {
		$map{"jmax"} = $1;
	}
	if(m/dy = ([-0-9.de*dsqrt(fac_y)\/]+)/) {
		$map{"dy"} = $1;
	}
	if(m/stf = ([-0-9.de]+)/) {
		$map{"stf"} = $1;
	}
	if(m/x0 = ([-0-9.de]+)/) {
		$map{"x0"} = $1;
	}
	if(m/Pr = ([-0-9.de]+)/) {
		$map{"pr"} = $1;
	}
	if(m/msh = (\d+)/) {
		$map{"lvl"} = $1;
	}
}

close FD;

open(FD, "$fileout");
open(FDRES, ">$fileres");
while(<FD>) {
	s/__MYFORM__/$map{"myform"}/g;
	s/__U1__/$map{"u1"}/g;
	s/__L1__/$map{"l1"}/g;
	s/__N1__/$map{"n1"}/g;
	s/__UN__/$map{"un"}/g;
	s/__RE__/$map{"re"}/g;
	s/__FAC_Y__/$map{"fac_y"}/g;
	s/__IMAX__/$map{"imax"}/g;
	s/__DX__/$map{"dx"}/g;
	s/__JMAX__/$map{"jmax"}/g;
	s/__DY__/$map{"dy"}/g;
	s/__STF__/$map{"stf"}/g;
	s/__X0__/$map{"x0"}/g;
	s/__PR__/$map{"pr"}/g;
	s/__LVL__/$map{"lvl"}/g;
	print FDRES;
}
close FD;
close FDRES;
