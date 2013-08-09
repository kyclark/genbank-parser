# NAME

Bio::GenBankParser - Parse::RecDescent parser for a GenBank record

# VERSION

Version 0.05

# SYNOPSIS

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

# METHODS

## new

    use Bio::GenBankParser;
    my $parser = Bio::GenBankParser->new;

## file

    $parser->file('/path/to/file');

Informs the parser to read sequentially from a file.

## current\_record

    my $genbank_record = $parser->current_record;

Returns the current unparsed GenBank record.

## next\_seq

    my $seq = $parser->next_seq;

Returns the next sequence from the `file`.

## parse

    my $rec = $parser->parse( $text );
    print $rec->{'ACCESSION'};

Parses a (single) GenBank record into a hash reference.

## parser

Returns the Parse::RecDescent object.

## grammar

Returns the Parse::RecDescent grammar for a GenBank record.

# AUTHOR

Ken Youens-Clark <kclark at cpan.org>.

# BUGS

Please report any bugs or feature requests to `bug-bio-genbankparser
at rt.cpan.org`, or through the web interface at
[http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-GenBankParser](http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-GenBankParser).
I will be notified, and then you'll automatically be notified of
progress on your bug as I make changes.

# SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::GenBankParser

- RT: CPAN's request tracker

    [http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-GenBankParser](http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-GenBankParser)

- AnnoCPAN: Annotated CPAN documentation

    [http://annocpan.org/dist/Bio-GenBankParser](http://annocpan.org/dist/Bio-GenBankParser)

- CPAN Ratings

    [http://cpanratings.perl.org/d/Bio-GenBankParser](http://cpanratings.perl.org/d/Bio-GenBankParser)

- Search CPAN

    [http://search.cpan.org/dist/Bio-GenBankParser](http://search.cpan.org/dist/Bio-GenBankParser)

# ACKNOWLEDGEMENTS

Lincoln Stein, Doreen Ware and everyone at Cold Spring Harbor Lab.

# COPYRIGHT & LICENSE

Copyright 2008 Ken Youens-Clark, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.
