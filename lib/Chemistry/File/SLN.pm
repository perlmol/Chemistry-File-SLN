package Chemistry::File::SLN;

$VERSION = "0.10";
# $Id$

use 5.006;
use strict;
use warnings;
use base "Chemistry::File";
use Chemistry::Mol;
use Chemistry::File::SLN::Parser;

=head1 NAME

Chemistry::File::SLN - SLN linear notation parser/writer

=head1 SYNOPSYS

    #!/usr/bin/perl
    use Chemistry::File::SLN;

    # parse a SLN string
    my $s = 'C1CC1(=O)[O-]';
    my $mol = Chemistry::Mol->parse($s, format => 'sln');

    # print a SLN string
    print $mol->print(format => 'sln');

    # print a unique (canonical) SLN string
    print $mol->print(format => 'sln', unique => 1);

    # parse a SLN file
    my @mols = Chemistry::Mol->read("file.smi", format => 'sln');

    # write a multiline SLN file
    Chemistry::Mol->write("file.smi", mols => [@mols]);


=head1 DESCRIPTION

This module parses a SLN (Sybyl Line Notation) string. This is a File I/O
driver for the PerlMol project.  L<http://www.perlmol.org/>. It registers the
'sln' format with Chemistry::Mol.

=cut

# INITIALIZATION
Chemistry::Mol->register_format('sln');
my $Parser = Chemistry::File::SLN::Parser->new;

sub name_is {
    my ($self, $name) = @_;
    $name =~ /\.sln/;
}

sub parse_string {
    my ($self, $string, %opts) = @_;
    my $mol_class = $opts{mol_class} || "Chemistry::Mol";

    # these two are not used yet...
    my $atom_class = $opts{atom_class} || $mol_class->atom_class;
    my $bond_class = $opts{bond_class} || $mol_class->bond_class;

    my $mol = $mol_class->new;
    my $tree = $Parser->run($string);
    use Data::Dumper; print Dumper $tree;
    my @nodes = @{$tree->{chain}};
    my $node = shift @nodes;
    my %closures;
    my $last_atom = compile_atom($mol, $node, \%closures);
    my @stack;
    
    while ($node = shift @nodes) {
        if ($node eq '(') {
            push @stack, $last_atom;
        } elsif ($node eq ')') {
            $last_atom = pop @stack;
        } else { # bond
            my $next = shift @nodes;
            if ($next->{closure}) {
                my $atom = $closures{$next->{closure}};
                compile_bond($mol, $node, $last_atom, $atom);
            } else {
                my $atom = compile_atom($mol, $next, \%closures);
                compile_bond($mol, $node, $last_atom, $atom);
                $last_atom = $atom;
            }
        }
    }
    $mol;
}

sub compile_atom {
    my ($mol, $node, $closures) = @_;
    my $atom = $mol->new_atom(
        symbol => $node->{symbol},
        hydrogens => $node->{hcount},
    );
    $closures->{$node->{id}} = $atom if $node->{id};
    $atom;
}

my %TYPE_TO_ORDER = (
    '-' => 1,
    '=' => 2,
    '#' => 3,
    ':' => 1, 
    '.' => 0,
);

sub compile_bond {
    my ($mol, $node, $atom1, $atom2) = @_;
    my $order = $TYPE_TO_ORDER{$node->{type}};
    if ($order) {
        my $bond = $mol->new_bond(
            type => $node->{type}, 
            atoms=>[$atom1, $atom2],
            order => $order,
        );
        if ($node->{type} eq ':') { 
            $_->aromatic(1) for ($atom1, $atom2, $bond);
        }
    }
}

1;

=head1 CAVEATS

Reading branches that start before an atom, such as (OC)C, which should be
equivalent to C(OC) and COC, according to some variants of the SMILES
specification. Many other tools don't implement this rule either.

=head1 VERSION

0.10

=head1 SEE ALSO

L<Chemistry::Mol>, L<Chemistry::File>

The PerlMol website L<http://www.perlmol.org/>

=head1 AUTHOR

Ivan Tubert E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2004 Ivan Tubert. All rights reserved. This program is free
software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

