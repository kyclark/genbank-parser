package Bio::GenBankParser;

use warnings;
use strict;
use Carp qw( croak );
use File::Spec::Functions;
use Parse::RecDescent;
use Readonly;

Readonly my $GENBANK_RECORD_SEPARATOR => "//\n";

=pod

=head1 NAME

Bio::GenBankParser - Parse::RecDescent parser for a GenBank record

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

This module aims to improve on the BioPerl GenBank parser by using
a grammar-based approach with Parse::RecDescent.

    use Bio::GenBankParser;

    my $parser = Bio::GenBankParser->new();

    local $/ = "//\n";
    while ( my $rec = <$input> ) {
        my $gb_rec = $parser->parse( $rec );
    }

Or:

    my $parser = Bio::GenBankParser->new( file => $file );
    while ( my $seq = $parser->next_seq ) {
        ...
    }

=head1 METHODS

=cut

# ----------------------------------------------------------------
sub new {

=pod

=head2 new

=cut

    my $class = shift;
    my %args  = ref $_[0] eq 'HASH' ? %{ $_[0] } : @_;
    my $self  = bless \%args, $class;

    if ( $args{'file'} ) {
        $self->file( $args{'file'} );
    }

    return $self;
}

# ----------------------------------------------------------------
sub DESTROY {
    my $self = shift;

    if ( my $fh = $self->{'fh'} ) {
        close $fh;
    }
}

# ----------------------------------------------------------------
sub file {

=pod

=head2 file

    $parser->file('/path/to/file');

Informs the parser to read sequentially from a file.

=cut

    my $self = shift;

    if ( my $file = shift ) {
        $file = canonpath( $file );

        if ( -e $file && -s _ && -r _ ) { 
            open my $fh, '<', $file or croak("Can't read file '$file'; $!\n");

            $self->{'file'} = $file;
            $self->{'fh'}   = $fh;
        }
        else {
            croak("Non-existent, empty or unreadable file: '$file'");
        }
    }

    return 1;
}

# ----------------------------------------------------------------
sub current_record {

=pod

=head2 current_record

    my $genbank_record = $parser->current_record;

Returns the current unparsed GenBank record.

=cut

    my $self = shift;

    return $self->{'current_record'};
}

# ----------------------------------------------------------------
sub next_seq {

=pod

=head2 next_seq

    my $seq = $parser->next_seq;

Returns the next sequence from the C<file>.

=cut

    my $self = shift;

    if ( my $fh = $self->{'fh'} ) {
        local $/ = $GENBANK_RECORD_SEPARATOR; 

        if ( my $rec = <$fh> ) {
            return $self->parse( $rec );
        }
        else {
            return undef;
        }
    }
    else {
        croak("Can't call 'next_seq' without a 'file' argument");
    }
}

# ----------------------------------------------------------------
sub parse {

=pod

=head2 parse

    my $rec = $parser->parse( $text );
    print $rec->{'ACCESSION'};

Parses a (single) GenBank record into a hash(ref).

=cut

    my $self   = shift;
    my $text   = shift() or croak('No input to parse');
    my $parser = $self->parser or croak('No parser');

    $self->{'current_record'} = $text;

    return $parser->startrule( $text );
}

# ----------------------------------------------------------------
sub parser {

=pod

=head2 parser

Returns the Parse::RecDescent object.

=cut

    my $self = shift;

    if ( !defined $self->{'parser'} ) {
        my $grammar = $self->grammar or croak('No grammar');
        $self->{'parser'} = Parse::RecDescent->new( $grammar );
    }

    return $self->{'parser'};
}

# ----------------------------------------------------------------
sub grammar {

=pod

=head2 grammar

Returns the Parse::RecDescent grammar for a GenBank record.

=cut

    my $self = shift;
    return <<'END_OF_GRAMMAR';
{
    my $ref_num  = 1;
    my %record   = ();
    my %ATTRIBUTE_PROMOTE = map { $_, 1 } qw[ 
        mol_type 
        cultivar 
        variety 
        strain 
    ];

    $::RD_ERRORS; # report fatal errors
#    $::RD_TRACE  = 0;
#    $::RD_WARN   = 0; # Enable warnings. This will warn on unused rules &c.
#    $::RD_HINT   = 0; # Give out hints to help fix problems.
}

startrule: section(s) eofile 
    { 
        if ( !$record{'ACCESSION'} ) {
            $record{'ACCESSION'} = $record{'LOCUS'}->{'genbank_accession'};
        }

        if ( ref $record{'SEQUENCE'} eq 'ARRAY' ) {
            $record{'SEQUENCE'} = join('', @{ $record{'SEQUENCE'} });
        }

        $return = { %record };
        %record = ();
    }
    | <error>

section: header
    | locus
    | definition
    | accession_line
    | version_line
    | keywords
    | source_line
    | organism
    | reference
    | features
    | origin
    | comment
    | record_delimiter
    | <error>

header: /.+(?=\nLOCUS)/xms

locus: /LOCUS/xms genbank_accession sequence_length molecule_type 
    genbank_division modification_date
    {
        $record{'LOCUS'} = {
            genbank_accession => $item{'genbank_accession'},
            sequence_length   => $item{'sequence_length'},
            molecule_type     => $item{'molecule_type'},
            genbank_division  => $item{'genbank_division'},
            modification_date => $item{'modification_date'},
        }
    }

space: /\s+/

genbank_accession: /[A-Z]{1}\d{5}/ | /[A-Z]{2}\d{6}/ 

genbank_version: /[A-Z]{1}\d{5}\.\d+/ | /[A-Z]{2}\d{6}\.\d+/

genbank_gi: /GI:\d+/

sequence_length: /\d+/ 'bp' { $return = "$item[1] $item[2]" }

molecule_type: /\w+/ (/\w+/)(?) 
    { $return = join(' ', map { $_ || () } $item[1], $item[2][0] ) }

genbank_division: 
    /(PRI|ROD|MAM|VRT|INV|PLN|BCT|VRL|PHG|SYN|UNA|EST|PAT|STS|GSS|HTG|HTC|ENV)/

modification_date: /\d+-[A-Z]{3}-\d{4}/

definition: /DEFINITION/ section_continuing_indented
    {
        ( $record{'DEFINITION'} = $item[2] ) =~ s/\n\s+/ /g;
    }

section_continuing_indented: /.*?(?=\n[A-Z]+\s+)/xms

accession_line: /ACCESSION/ genbank_accession(s)
    {
        my @accs = @{ $item[2] };

        if ( scalar @accs <= 1 ) {
            $record{'ACCESSION'} = shift @accs;
        }
        else {
            push @{ $record{'VERSION'} }, @accs;
        }
    }

version_line: /VERSION/ synonym(s)
    {
        push @{ $record{'VERSION'} }, @{ $item[2] };
    }

synonym: genbank_version | genbank_gi

keywords: /KEYWORDS/ keyword_value
    { 
        $record{'KEYWORDS'} = $item[2];
    }

keyword_value: /([^.]+)[.]?(?=\n)/ { $return = [ split(/,\s*/, $1) ] }
    | PERIOD { $return = [] }

source_line: /SOURCE/ source_value 
    { 
        $record{'SOURCE'} = $item[2];
    }

source_value: /(.+?)(?=\n\s{0,2}[A-Z]+)/xms { $return = $1 }

organism: organism_line classification_line
    { 
        $record{'ORGANISM'} = $item[1];
        $record{'CLASSIFICATION'} = $item[2];
    }

organism_line: /ORGANISM/ organism_value { $return = $item[2] }

organism_value: /([^\n]+)(?=\n)/xms { $return = $1 }

classification_line: /([^.]+)[.]/xms { $return = [ split(/;\s*/, $1) ] }

word: /\w+/

reference: /REFERENCE/ NUMBER(?) parenthetical_phrase(?) authors title journal pubmed(?)
    {
        my $num = $item[2][0] || $ref_num++;

        push @{ $record{'REFERENCES'} }, {
            number  => $num,
            authors => $item{'authors'},
            title   => $item{'title'},
            journal => $item{'journal'},
            pubmed  => $item[7][0],
            note    => $item[3][0],
        };
    }

parenthetical_phrase: /\( ([^)]+) \)/xms
    {
        $return = $1;
    }

authors: /AUTHORS/ author_value { $return = $item[2] }

author_value: /(.+?)(?=\n\s{0,2}[A-Z]+)/xms 
    { 
        $return = [ 
            grep  { !/and/ }      
            map   { s/,$//; $_ } 
            split /\s+/, $1
        ];
    }

title: /TITLE/ /.*?(?=\n\s{0,2}[A-Z]+)/xms
    { ( $return = $item[2] ) =~ s/\n\s+/ /; }

journal: /JOURNAL/ /.*?(?=\n\s+PUBMED)/xms 
    { 
        ( $return = $item[2] ) =~ s/\n\s+/ /g;
    }
    | /JOURNAL/ journal_value { $return = $item[2] }

journal_value: /(.+?)(?=\n[A-Z]+)/xms { $return = $1; $return =~ s/\n\s+/ /g; }

pubmed: /PUBMED/ NUMBER
    { $return = $item[2] }

features: /FEATURES/ section_continuing_indented
    { 
        my ( $location, $cur_feature_name, %cur_features, $cur_key );
        for my $fline ( map { s/^\s+|\s+$//g; $_ } split(/\n/, $item[2]) ) {
            next if $fline eq 'Location/Qualifiers';
            next if $fline !~ /\S+/;

            if ( $fline =~ /^\/ (\w+?) = (.+)$/xms ) {
                my ( $key, $value )   = ( $1, $2 );
                $value                =~ s/^"|"$//g;
                $cur_key              = $key;
                $cur_features{ $key } = $value;

                if ( $key eq 'db_xref' && $value =~ /^taxon:(\d+)$/ ) {
                    $record{'NCBI_TAXON_ID'} = $1;
                }

                if ( $ATTRIBUTE_PROMOTE{ $key } ) {
                    $record{ uc $key } = $value;
                }
            }
            elsif ( $fline =~ /^(\S+) \s+ (.+)$/xms ) {
                my ( $this_feature_name, $this_location ) = ( $1, $2 );
                $cur_key = '';

                if ( $cur_feature_name ) {
                    push @{ $record{'FEATURES'} }, {
                        name     => $cur_feature_name,
                        location => $location,
                        feature  => { %cur_features },
                    };

                    %cur_features = ();
                }

                ( $cur_feature_name, $location ) = 
                    ( $this_feature_name, $this_location );
            }
            elsif ( $fline =~ /^(\w+)["]?$/ ) {
                if ( $cur_key ) {
                    $cur_features{ $cur_key } .= $1;
                }
            }
        }

        push @{ $record{'FEATURES'} }, {
            name     => $cur_feature_name,
            location => $location,
            feature  => { %cur_features },
        };
    }

origin: /ORIGIN/ origin_value { $record{'ORIGIN'} = $item[2] }

origin_value: /(.*?)(?=\n\/\/)/xms
    {
        my $seq = $1;
        $record{'SEQUENCE'} = [];
        while ( $seq =~ /([actg]+)/gc ) {
            push @{ $record{'SEQUENCE'} }, $1;
        }

        $return = $seq;
    }

comment: /COMMENT/ comment_value

comment_value: /(.+?)(?=\n[A-Z]+)/xms
    { 
        $record{'COMMENT'} = $1;
    }

NUMBER: /\d+/

PERIOD: /\./

record_delimiter: /\/\/\s*/xms

eofile: /^\Z/

END_OF_GRAMMAR
}

# ----------------------------------------------------------------
=pod

=head1 AUTHOR

Ken Youens-Clark E<lt>kclark at cshl.eduE<gt>.

=head1 BUGS

Please report any bugs or feature requests to C<bug-bio-genbankparser
at rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-GenBankParser>.
I will be notified, and then you'll automatically be notified of
progress on your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::GenBankParser

=over 4

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-GenBankParser>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Bio-GenBankParser>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Bio-GenBankParser>

=item * Search CPAN

L<http://search.cpan.org/dist/Bio-GenBankParser>

=back

=head1 ACKNOWLEDGEMENTS

Lincoln Stein, Doreen Ware and everyone at Cold Spring Harbor Lab.

=head1 COPYRIGHT & LICENSE

Copyright 2008 Ken Youens-Clark, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of Bio::GenBankParser
