package Chemistry::File::SLN;

$VERSION = "0.10";
# $Id$

use 5.006;
use strict;
use warnings;
use base "Chemistry::File";
use Chemistry::Mol;
use Chemistry::File::SLN::Parser;
use Chemistry::Bond::Find 'assign_bond_orders';
use List::Util qw(sum);

=head1 NAME

Chemistry::File::SLN - SLN linear notation parser/writer

=head1 SYNOPSYS

    #!/usr/bin/perl
    use Chemistry::File::SLN;

    # parse a SLN string for benzene
    my $s = 'C[1]H:CH:CH:CH:CH:CH@1';
    my $mol = Chemistry::Mol->parse($s, format => 'sln');

    # print a SLN string
    print $mol->print(format => 'sln');

    # print a unique (canonical) SLN string
    print $mol->print(format => 'sln', unique => 1);

    # parse a multiline SLN file
    my @mols = Chemistry::Mol->read("file.sln", format => 'sln');

    # write a multiline SLN file
    Chemistry::Mol->write("file.sln", mols => [@mols]);


=head1 DESCRIPTION

This module parses a SLN (Sybyl Line Notation) string. This is a File I/O
driver for the PerlMol project.  L<http://www.perlmol.org/>. It registers the
'sln' format with Chemistry::Mol.

=head1 OPTIONS

=head1 CAVEATS

This version does not implement the full SLN specification. It supports
simple structures and some attributes, but it does not support any of the
following:

=over

=item Macro atoms

=item Pattern matching options

=item Markush structures

=back

If the parser doesn't understand a string, it only says "syntax error", which
may not be very helpful.

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

    # call the actual yapp-generated parser
    my $tree = $Parser->run($string) or return;
    #use Data::Dumper; print Dumper $tree;
    my @nodes = @{$tree->{chain}};
    #my $node = shift @nodes;
    my %closures;
    my $last_atom;
    my @stack;
    
    while (my $node = shift @nodes) {
        if ($node eq '(') {
            push @stack, $last_atom;
        } elsif ($node eq ')') {
            $last_atom = pop @stack;
        } elsif($last_atom) { # bond
            my $next = shift @nodes;
            if ($next->{closure}) {
                my $atom = $closures{$next->{closure}};
                compile_bond($mol, $node, $last_atom, $atom);
            } else {
                my $atom = compile_atom($mol, $next, \%closures);
                compile_bond($mol, $node, $last_atom, $atom);
                $last_atom = $atom;
            }
        } else {  # first atom
            $last_atom = compile_atom($mol, $node, \%closures);
        }
    }
    if ($opts{kekulize}) {
        assign_bond_orders($mol, method => "itub", use_coords => 0, 
            scratch => 0, charges => 0);
    }
    my @sln_attr;
    while (my ($attr, $value) = each %{$tree->{attr}}) {
        if ($attr eq 'name') {
            $mol->name($value);
        } elsif ($attr eq 'type') {
            $mol->type($value);
        } elsif ($attr eq 'coord3d') {
            read_coords($mol, $value);
        } else {
            push @sln_attr, $attr, $value;
        }
    }
    $mol->attr("sln/attr", {@sln_attr}) if @sln_attr;
    $mol;
}

sub compile_atom {
    my ($mol, $node, $closures) = @_;
    my $atom = $mol->new_atom(
        symbol          => $node->{symbol},
        hydrogens       => $node->{hcount},
        formal_charge   => $node->{attr}{charge},
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

sub read_coords {
    my ($mol, $coords_str) = @_;
    $coords_str =~ s/[()]//g;
    my (@coords) = split /,/, $coords_str;
    my $fh = $mol->formula_hash;
    my $n = sum(values %$fh);
    my $sprout = (@coords == 3*$n);
    for my $atom ($mol->atoms) {
        $atom->coords(splice @coords, 0, 3);
        if ($sprout) {
            for (1 .. $atom->implicit_hydrogens) {
                my $H = $mol->new_atom(symbol => 'H', 
                    coords => [splice @coords, 0, 3]);
                $mol->new_bond(atoms => [$atom, $H]);
            }
            $atom->implicit_hydrogens(0);
        } 
    }
}


########### WRITER #################


sub write_string {
    my ($self, $mol_ref, %opts) = @_;

    my $eol;
    my @mols;
    if ($opts{mols}) {
        @mols = @{$opts{mols}};
        $eol = "\n";
    } else {
        @mols = $mol_ref; 
        $eol = "";
    }

    my $sln;
    for my $mol (@mols) {
        my $oldmol = $mol;
        $mol = $mol->clone; 
        collapse_hydrogens($mol);
        my @atoms = $mol->atoms; 

        my @id_log;
        if (@atoms) {
            if ($opts{unique}) {
                unless ($atoms[0]->attr("canon/class")) {
                    require Chemistry::Canonicalize;
                    Chemistry::Canonicalize::canonicalize($mol);
                }
                $opts{aromatic} = 1; # all unique sln have to be aromatic
                @atoms = sort {
                    $a->attr("canon/class") <=> $b->attr("canon/class")
                } @atoms;
            }

            if ($opts{aromatic}) {
                require Chemistry::Ring;
                Chemistry::Ring::aromatize_mol($mol);
            }

            my $visited = {};
            my @s;
            for my $atom (@atoms) {
                next if $visited->{$atom};
                my $ring_atoms = {};

                # first pass to find and number the ring bonds
                find_ring_bonds($mol, \%opts, $atom, undef, {}, $ring_atoms);

                # second pass to actually generate the sln string
                push @s, branch($mol, \%opts, $atom, undef, $visited, 
                    $ring_atoms, \@id_log);
            }
            $sln .= join '.', @s;
        }

        if ($opts{name} or $opts{attr} or $opts{coords} or $opts{coord3d}) {
            my @attr;
            push @attr, 'name="' . $mol->name . '"' 
                if $opts{name} and length $mol->name;
            my @coords;
            if ($opts{coord3d} or $opts{coords}) {
                my @all_atoms = map { 
                    (
                        $oldmol->by_id($_), 
                        grep {$_->symbol eq 'H'}
                            $oldmol->by_id($_)->neighbors
                    )
                } @id_log;
                push @coords, sprintf("(%.3f,%.3f,%.3f)",$_->coords->array)
                    for @all_atoms;
                push @attr, 'coord3d=' . join(',',@coords);
            }
            if ($opts{attr}) {
                while (my ($key, $val) = each %{$mol->attr("sln/attr")||{}}) {
                    push @attr, "$key" . ($val eq 'TRUE' ? "" : "=$val");
                }
            }
            $sln .= '<' . join(';', @attr) . '>' if @attr;
        }
        $sln .= $eol;
    }
    $sln;
}

sub find_ring_bonds {
    my ($mol, $opts, $atom, $from_bond, $visited, $ring_atoms) = @_;

    $visited->{$atom}  = 1;
    for my $bn (sorted_bonds_neighbors($atom, $opts)) {
        my $nei  = $bn->{to};
        my $bond = $bn->{bond};
        next if $visited->{$bond};
        $visited->{$bond}  = 1;
        if ($visited->{$nei}) { # closed ring
            #print "closing ring\n";
            $ring_atoms->{$nei}++;
        } else {
            find_ring_bonds($mol, $opts, $nei, $bond, $visited, $ring_atoms);
        }
    }
}

sub branch {
    my ($mol, $opts, $atom, $from_bond, $visited, $digits, $id_log) = @_;

    my $prev_branch = "";
    my $sln;
    $sln .= bond_symbol($from_bond, $opts);
    my $digit;
    if ($digits->{$atom}) {  # opening a ring
        $digit = next_digit($digits);
        $digits->{$atom} = $digit;
    }
    $sln .= format_atom($atom, $opts, $digit);
    push @$id_log, $atom->id;

    $visited->{$atom}  = 1;
    my @bns = sorted_bonds_neighbors($atom, $opts);

    for my $bn (@bns) {
        my $nei  = $bn->{to};
        my $bond = $bn->{bond};
        next if $visited->{$bond};
        $visited->{$bond} = 1;
        if ($visited->{$nei}) { # closed a ring
            if ($prev_branch) {
                $sln .= "($prev_branch)";
            }
            $prev_branch = bond_symbol($bond, $opts) . '@' . $digits->{$nei};
            $visited->{$bond} = 1;
        } else {
            my $branch = branch($mol, $opts, $nei, $bond, $visited, 
                $digits, $id_log);
            if ($prev_branch) {
                $sln .= "($prev_branch)";
            }
            $prev_branch = $branch;
        }
    }
    $sln .= "$prev_branch";
    $sln;
}

sub next_digit {
    my ($digits) = @_;
    ++$digits->{used_digits};
}

sub collapse_hydrogens {
    my ($mol) = @_;

    for my $atom (grep {$_->symbol eq 'H'} $mol->atoms) {
        my ($neighbor) = $atom->neighbors or next;
        $atom->delete;
        my $h_count = $neighbor->hydrogens;
        $h_count++;
        $neighbor->hydrogens($h_count);
    }
}

sub sorted_bonds_neighbors {
    my ($atom, $opts) = @_;
    my @bn = $atom->bonds_neighbors;
    if ($opts->{unique}) {
        @bn = sort { 
            $a->{to}->attr("canon/class") <=> $b->{to}->attr("canon/class") 
        } @bn;
    }
    @bn;
}

my %ORDER_TO_TYPE = (
    2 => '=', 1 => '', 3 => '#',
);

sub bond_symbol {
    my ($bond, $opts) = @_;
    return '' unless $bond;
    return ':' if $bond->aromatic;
    return $ORDER_TO_TYPE{$bond->order};
}

sub format_atom {
    my ($atom, $opts, $digit) = @_;
    my $s;
    no warnings 'uninitialized';
    my $h_count = $atom->hydrogens;
    my $charge  = $atom->formal_charge;
    my $symbol  = $atom->symbol;

    $charge  = $charge ? sprintf("%+d", $charge): '';
    $h_count = $h_count ? ($h_count > 1 ? "H$h_count" : 'H') : '';

    $s = $symbol;
    if ($charge or $digit) {
        $s .= '['; 
        $s .= $digit;
        $s .= ':' if $charge and $digit;
        $s .= $charge;
        $s .= ']';
    }
    $s .= $h_count;
    $s;
}


1;

=head1 VERSION

0.10

=head1 SEE ALSO

L<Chemistry::Mol>, L<Chemistry::File>, L<Chemistry::File::SMILES>

The PerlMol website L<http://www.perlmol.org/>

Ash, S.; Cline, M. A.; Homer, R. W.; Hurst, T.; Smith, G. B., SYBYL Line
Notation (SLN): A Versatile Language for Chemical Structure Representation. J.
Chem. Inf. Comput. Sci; 1997; 37(1); 71-79.  DOI: 10.1021/ci960109j

=head1 AUTHOR

Ivan Tubert E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2004 Ivan Tubert. All rights reserved. This program is free
software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

