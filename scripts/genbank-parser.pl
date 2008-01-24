#!/usr/local/bin/perl

# vim: tw=78: sw=4: ts=4: et: 

# $Id: $

use strict;
use warnings;
use Bio::GenBankParser;
use English qw( -no_match_vars );
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Readonly;

Readonly my $VERSION => sprintf '%d.%02d', 
                        qq$Revision: 1.01$ =~ /(\d+)\.(\d+)/;

my ( $help, $man_page, $show_version );
GetOptions(
    'help'    => \$help,
    'man'     => \$man_page,
    'version' => \$show_version,
) or pod2usage(2);

if ( $help || $man_page ) {
    pod2usage({
        -exitval => 0,
        -verbose => $man_page ? 2 : 1
    });
}; 

if ( $show_version ) {
    my $prog = basename( $PROGRAM_NAME );
    print "$prog v$VERSION\n";
    exit 0;
}

my @files   = @ARGV or die 'No input files';
my $grammar = join( '', <DATA> );
my $parser  = Bio::GenBankParser->new;

for my $file ( @files ) {
    print "Processing '$file'\n";
    open my $fh, '<', $file or die "Can't read $file: $!\n";

    local $/ = "//\n";

    while ( my $rec = <$fh> ) {
        print "record = ", Dumper($parser->startrule( $rec )), "\n";
        last;
    }
}

__END__

# ----------------------------------------------------

=pod

=head1 NAME

genbank-parser.pl - a script

=head1 VERSION

This documentation refers to version $Revision: 1.01$

=head1 SYNOPSIS

  genbank-parser.pl 

Options:

  --help        Show brief help and exit
  --man         Show full documentation
  --version     Show version and exit

=head1 DESCRIPTION

Describe what the script does, what input it expects, what output it
creates, etc.

=head1 SEE ALSO

perl.

=head1 AUTHOR

Ken Youens-Clark E<lt>kclark@cshl.eduE<gt>.

=head1 COPYRIGHT

Copyright (c) 2008 Cold Spring Harbor Laboratory

This module is free software; you can redistribute it and/or
modify it under the terms of the GPL (either version 1, or at
your option, any later version) or the Artistic License 2.0.
Refer to LICENSE for the full license text and to DISCLAIMER for
additional warranty disclaimers.

=cut
